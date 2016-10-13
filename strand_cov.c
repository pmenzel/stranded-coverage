#include "htslib/hfile.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <getopt.h>
#include <inttypes.h>
#include <string.h>

void usage();

typedef struct { // helper data structure
	samFile * fp; // file handle to one BAM file
	bam_hdr_t * hdr; // BAM header
	hts_itr_t * iter; // NULL if a region not specified
	int min_mapq;
	double norm_factor;
} mplp_data;

int filter_func(void *data, bam1_t *b) {
	int ret;
	mplp_data * d = (mplp_data *) data;
	// get next useful read
	while(1) {
		ret = d->iter ? sam_itr_next(d->fp, d->iter, b) : sam_read1(d->fp, d->hdr, b);
		if(ret < 0) break;
		if(b->core.tid < 0 || b->core.flag & BAM_FUNMAP) continue; // UNMAPPED
		if(b->core.flag & (BAM_FQCFAIL | BAM_FDUP) ) continue; // QCFAIL and DUPLICATE
		if(b->core.flag & BAM_FSECONDARY) continue; // SECONDARY, that would skip alternative alignments for multi-mapping reads
		if((int)b->core.qual < d->min_mapq) continue; // min mapping quality (option -m)
		//TODO why set the 0x2 flag (==read mapped in proper pair) for discordant pairs here ?
		if((b->core.flag & 0x9) == 0x1) b->core.flag |= 0x2; //Discordant pairs can cause double counts
		break;
	}
	return ret;
}

int processDepths(mplp_data ** data, int n_files, char * reg, int min_phred, FILE *of1, FILE *of2, double norm_factor) {
	bam_mplp_t mplp;
	int tid;
	int pos;
	int f,i;
	int * n_plp; // number of covering reads per position
	const bam_pileup1_t **plp = NULL;
	uint32_t plus_depth, minus_depth;
	int32_t ctid = -1;
	int status = 0;
	bam_hdr_t * h = NULL; // BAM header of the 1st input
	//int reg_tid = 0;
	int beg = 0, end = INT_MAX;  // set the default region


	h = data[0]->hdr; // easy access to the header of the 1st BAM
	if(reg) {
		beg = data[0]->iter->beg; // and to the parsed region coordinates
		end = data[0]->iter->end;
		//reg_tid = data[0]->iter->tid;
	}

	// Initialize mplpator over positions
	mplp = bam_mplp_init(n_files, filter_func, (void **)data);
	bam_mplp_init_overlaps(mplp); // set the quality scores of overlapping bases to 0 in one of the reads
	bam_mplp_set_maxcnt(mplp, 1e6); // TODO why 1e6 ?

	n_plp = calloc((size_t)n_files, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
	plp = calloc((size_t)n_files, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)

	// iterate over all covered positions
	int ret;
	while((ret = bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) {
		if (pos < beg || pos >= end) continue; // out of range; skip
		plus_depth = 0; minus_depth = 0;
		// check if entering another reference sequence (e.g. chromosome)
		// then output the new name
		if(tid != ctid) {
			ctid = tid;
			fprintf(of1, "variableStep chrom=%s span=1\n", h->target_name[tid]);
			fprintf(of2, "variableStep chrom=%s span=1\n", h->target_name[tid]);
		}
		for(f=0;f<n_files;f++) {
			for(i=0; i<n_plp[f]; i++) { // go through each read that covers position pos
				// do not count read towards coverage if it has a deletion or skips the reference (split reads)
				if(plp[f][i].is_del || plp[f][i].is_refskip) continue;
				// do not count read if below minimum quality
				if(bam_get_qual(plp[f][i].b)[plp[f][i].qpos] < min_phred) continue;

				// determine strand of read
				if(plp[f][i].b->core.flag & BAM_FREAD2) { //Read #2
					if(plp[f][i].b->core.flag & BAM_FREVERSE) { //- strand
						minus_depth++;
					} else { //+ strand
						plus_depth++;
					}
				} else {
					if(plp[f][i].b->core.flag & BAM_FREVERSE) { //+ strand
						plus_depth++;
					} else { //- strand
						minus_depth++;
					}
				}
			}
		}

		if(norm_factor == 0.0 ) {
			if(plus_depth>0)
				fprintf(of1, "%"PRId32"\t%"PRIu32"\n", pos+1, plus_depth); //wig is 1-based whereas BAM is 0-based, thus add 1
			if(minus_depth>0)
				fprintf(of2, "%"PRId32"\t%"PRIu32"\n", pos+1, minus_depth);
		}
		else {
			if(plus_depth>0)
				fprintf(of1, "%"PRId32"\t%.6f\n", pos+1, plus_depth * norm_factor);
			if(minus_depth>0)
				fprintf(of2, "%"PRId32"\t%.6f\n", pos+1, minus_depth * norm_factor);
		}

	} // end while
	if(ret < 0)
		status = 1;

	bam_mplp_destroy(mplp);
	free(plp);
	free(n_plp);
	return status;
}


// ========================= MAIN =====================================

int main(int argc, char *argv[]) {
	int strand = 2; // either 1 or 2, default is 2 (e.g. dUTP protocols), Ovattion protocol would be 1
	int min_mapq = 0;
	int min_phred = 1;
	FILE * of1 = NULL;
	FILE * of2 = NULL;
	char * reg = NULL; // specified region
	char * prefix = NULL; // prefix for output files
	int ret = 0; // return value for main()
	int n_files = 0;
	mplp_data ** data;
	int debug = 0;
	int verbose = 0;
	int norm_rpm = 0;
	uint64_t mapped_reads = 0;

	int c;
	while((c = getopt(argc, argv, "s:m:p:r:o:hdvn")) >= 0) {
		switch(c) {
			case 'o' :
				prefix = malloc(strlen(optarg)+1);
				if(prefix==0) {
					fprintf(stderr,"Error while copying prefix string :(\n");
					return 1;
				}
				strcpy(prefix,optarg);
				break;
			case 'r' :
				reg = malloc(strlen(optarg)+1);
				if(reg==NULL) {
					fprintf(stderr,"Error while copying region string :(\n");
					return 1;
				}
				strcpy(reg,optarg);
				break;
			case 'v':
				verbose = 1;
				break;
			case 'n':
				norm_rpm = 1;
				break;
			case 'd':
				debug = 1;
				break;
			case 'h' :
				usage(argv[0]);
				return 0;
			case 's' :
				strand = atoi(optarg);
				if(strand != 1 && strand != 2) {
					fprintf(stderr, "-s must be either 1 or 2. We'll use the default of 2.\n");
					fflush(stderr);
					strand = 2;
				}
				break;
			case 'm' :
				min_mapq = atoi(optarg);
				if(min_mapq < 0)
					min_mapq = 0;
				break;
			case 'p' :
				min_phred = atoi(optarg);
				if(min_phred < 1)
					min_phred = 1;
				break;
			default :
				fprintf(stderr, "Invalid option '%c'\n", c);
				usage(argv[0]);
				return 1;
		}
	}

	n_files = argc - optind; // the number of dditional arguments denote the input BAM files
	if(debug) fprintf(stderr,"Number of BAM files to read = %i\n",n_files);
	if(n_files < 1) { // at least one file argument is needed
		fprintf(stderr,"Error: Need one or more BAM file names.\n");
		usage(argv[0]);
		return 1;
	}

	if(prefix==NULL) {
		fprintf(stderr,"Error: Need to specify the prefix for output files.\n");
		usage(argv[0]);
		return 1;
	}

	if(norm_rpm && reg!=NULL) {
		fprintf(stderr,"Error: Options -r and -n cannot be used together.\n");
		usage(argv[0]);
		return 1;
	}

	// -------- open files and build data_mplp structures ----------------


	if((data = calloc((size_t)n_files, sizeof(mplp_data*)))==NULL) { // data[i] for the i-th input
		fprintf(stderr,"Error: Could not allocate memory\n");
		ret = 1;
		goto the_end;
	}

	int i;
	for(i = 0; i < n_files; i++) {
		data[i] = calloc(1, sizeof(mplp_data));
		data[i]->iter = NULL;
		data[i]->min_mapq = min_mapq;  // set the mapQ filter
		data[i]->norm_factor = 1.0;  // default is no normalisation

		// when using RPM normalization, need to count the number of mapped reads in the BAM file first
		if(norm_rpm) {
			bam1_t * aln;
			aln = bam_init1();
			samFile * in = sam_open(argv[optind+i], "r");
			bam_hdr_t *hdr = NULL;
			if((hdr = sam_hdr_read(in))==NULL) {    // read the BAM header
				fprintf(stderr, "Fail to read header from BAM file %s.\n", argv[optind+i]);
				ret = 1;
				goto the_end;

			}
			while(sam_read1(in, hdr, aln)>1) {
					if(aln->core.flag & BAM_FUNMAP || aln->core.tid < 0) continue; // UNMAPPED
					if(aln->core.flag & (BAM_FQCFAIL | BAM_FDUP) )continue; // QCFAIL and DUPLICATE
					if(aln->core.flag & BAM_FSECONDARY) continue; // SECONDARY
					if((int)aln->core.qual < min_mapq) continue;

					mapped_reads++;
			}
			sam_close(in);
		}

		if(verbose) fprintf(stderr,"Reading file %s\n",argv[optind+i]);
		if((data[i]->fp = sam_open(argv[optind+i], "r"))==NULL) { // open i'th BAM file
			fprintf(stderr, "Fail to open file %s for reading.\n", argv[optind+i]);
			ret = 1;
			goto the_end;
		}
		if((data[i]->hdr = sam_hdr_read(data[i]->fp))==NULL) {    // read the BAM header
			fprintf(stderr, "Fail to read header from BAM file %s.\n", argv[optind+i]);
			ret = 1;
			goto the_end;
		}
		if(reg) { // if a region is specified
			hts_idx_t *idx = sam_index_load(data[i]->fp, argv[optind+i]);  // load the index
			if(idx == NULL) {
				fprintf(stderr, "Fail to read index for BAM file %s.\n", argv[optind+i]);
				ret = 1;
				goto the_end;
			}
			data[i]->iter = sam_itr_querys(idx, data[i]->hdr, reg); // set the iterator
			hts_idx_destroy(idx); // the index is not needed any more; free the memory
			if(data[i]->iter == NULL) {
				fprintf(stderr, "Cannot parse region \"%s\".\n", reg);
				ret = 1;
				goto the_end;
			}
		}
	}

	//open output files using the prefix
	char *oname;
	if((oname = malloc(sizeof(char)*(strlen(prefix) + strlen(".minus.wig "))))==NULL) {
		fprintf(stderr, "Fail to allocate memory for output file names.");
		ret = 1;
		goto the_end;
	}

	sprintf(oname, "%s.plus.wig", prefix);
	if(verbose) fprintf(stderr,"Writing (+)strand coverage to file %s\n",oname);
	if(strand == 2)
		of1 = fopen(oname, "w");
	else
		of2 = fopen(oname, "w");

	sprintf(oname, "%s.minus.wig", prefix);
	if(verbose) fprintf(stderr,"Writing (-)strand coverage to file %s\n",oname);
	if(strand == 1)
		of1 = fopen(oname, "w");
	else
		of2 = fopen(oname, "w");

	free(oname);

	if(of1 == 0 || of2 == 0) {
		fprintf(stderr, "Fail to open output file for writing.\n");
		ret = 1;
		goto the_end;
	}

	if(verbose) fprintf(stderr,"Counted %lu mapped reads\n",mapped_reads);
	double norm_factor = (norm_rpm) ? 1.0e6 / (double)mapped_reads : 0.0;  // reads per million mapped reads
	ret = processDepths(data, n_files, reg, min_phred, of1, of2, norm_factor);

the_end:

	for (i = 0; i < n_files && data[i]; ++i) {
		bam_hdr_destroy(data[i]->hdr);
		if (data[i]->fp)
			sam_close(data[i]->fp);
		hts_itr_destroy(data[i]->iter);
		free(data[i]);
	}
	free(data);
	free(reg);
	if(of1)
		fclose(of1);
	if(of2)
		fclose(of2);

	return ret;
}


void usage(char *prog) {
	fprintf(stderr, "Usage: %s -o prefix [-s 1|2] [-r STRING] [-p INT] [-m INT] <in1.bam> [<in2.bam> ...]\n",prog);
	fprintf(stderr, "\n");
	fprintf(stderr, "This program reads sorted BAM files containing mapped reads from stranded RNASeq libraries\n");
	fprintf(stderr, "and produces two WIG files (prefix.plus.wig and prefix.minus.wig) with the\n");
	fprintf(stderr, "per-base depth of each strand.\n");
	fprintf(stderr, "\n");
	fprintf(stderr,	"Options:\n");
	fprintf(stderr,	"-o STRING  prefix for output files\n");
	fprintf(stderr,	"\n");
	fprintf(stderr,	"-s 1|2     Indicate which read in a pair denotes the strand.\n");
	fprintf(stderr,	"           The default is 2, meaning that read #2 denotes the strand of a pair.\n");
	fprintf(stderr,	"           This appropriate for dUTP-based libraries (for example Illumina TruSeq).\n");
	fprintf(stderr,	"           Use 1 for NuGEN stranded libraries.\n");
	fprintf(stderr,	"\n");
	fprintf(stderr,	"-m INT     Minimum required MAPQ value to count a read.\n");
	fprintf(stderr, "           The default is 0, meaning that all aligned reads are counted.\n");
	fprintf(stderr,	"\n");
	fprintf(stderr,	"-p INT     Minimum Phred score to include a base in an alignment. The default is 1.\n");
	fprintf(stderr,	"\n");
	fprintf(stderr, "-r STRING  Allows to specify a target region, e.g. 'chr3L:1000-2500' or '2L:1,000,000-2,000,000'\n");
	fprintf(stderr, "           This option requires a bam index file (.bai).\n");
	fprintf(stderr, "           See: samtools index \n");
	fprintf(stderr,	"\n");
	fprintf(stderr, "-n         normalization by million mapped reads (RPM) \n");
	fprintf(stderr,	"\n");
	fprintf(stderr, "-v         verbose output: samtools index \n");
	fprintf(stderr,	"\n");
	fprintf(stderr, "-d         debug output\n");
}

