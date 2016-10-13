strand_cov
==========

**strand_cov** reads sorted BAM files containing mapped paired-end RNASeq reads that originate
from a strand-specific library and calculates the coverage on each strand.

The output are two wig files, one for the plus strand and one for the minus strand.

The implementation is based on the file `bam2depth.c` from
the [samtools package](https://github.com/samtools/samtools)
and the strand splitting code from [here](https://github.com/dpryan79/Answers/tree/master/SEQanswers_48599).

### Installation

strand_cov depends on the [htslib](https://github.com/samtools/htslib), which should
be downloaded to `../htslib`. It is automatically build when building strand_cov.

Example:
```
git clone https://github.com/samtools/htslib.git
git clone https://github.com/pmenzel/strand_cov.git
cd strand_cov
make
```

### Usage
Example using the file `Aligned.out.bam`, which needs to be sorted by coordinates:
```
strand_cov -o coverage Aligned.out.bam
```
This will create the output files `coverage.plus.wig` and `coverage.minus.wig`.

### Conversion to bigWig format
The ouput files in wiggle format can be converted to bigWig format using the [wigToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig) program from UCSC, which
also requires the chromosome sizes in an extra file `chrom.sizes`.

For example:
```
wigToBigWig coverage.plus.wig chrom.sizes coverage.plus.bw
```


### Counting
Example of the counts for one mapped read-pair (the first read is split on a splice junction):
```
Genome       GGCCGGCATCGCATCAGCACATGCACACTGACACACACTGACTGGCTGCTGACTGACTGACTGCTGCTGCGCTATGCATGCCTGCTGAC
Read pair         GCATCGCA-----------ACACTGA>                  <ACTGACTGACTGCTGCTGCGCTATGCA
wig file          1111111100000000000111111100000000000000000000111111111111111111111111111
```


###License

See the file LICENSE.


