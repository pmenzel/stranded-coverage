stranded-coverage
=================

This program reads sorted and indexed BAM files with single-end or paired-end RNA-Seq reads that originate
from a strand-specific library and calculates the coverage on each strand.

The output are two wig files, one for the plus strand and one for the minus strand.

The implementation is based on the file `bam2depth.c` from
the [samtools package](https://github.com/samtools/samtools)
and the strand splitting code from [here](https://github.com/dpryan79/Answers/tree/master/SEQanswers_48599).

### Installation

The program depends on the [htslib](https://github.com/samtools/htslib), which should
be downloaded to `../htslib`.

Example:
```
git clone https://github.com/samtools/htslib.git
git clone https://github.com/pmenzel/stranded-coverage.git
cd stranded-coverage
make
```
This will produce the executable file `strand_cov`.

### Usage
Example using the file `Aligned.out.bam`, which needs to be sorted by coordinates:
```
strand_cov -o coverage Aligned.out.bam
```
This will create the output files `coverage.plus.wig` and `coverage.minus.wig`.

It is also possible to specify a region to be counted using the option `-r`:
```
strand_cov -r '3L:100000-200000' -o coverage Aligned.out.bam
```

### Counting
Example of the counts for one mapped read-pair (the first read is split on a splice junction):
```
Genome       GGCCGGCATCGCATCAGCACATGCACACTGACACACACTGACTGGCTGCTGACTGACTGACTGCTGCTGCGCTATGCATGCCTGCTGAC
Read pair         GCATCGCA-----------ACACTGA>                  <ACTGACTGACTGCTGCTGCGCTATGCA
wig file          1111111100000000000111111100000000000000000000111111111111111111111111111
```

The program uses the smart overlapping code from `htslib` to only count coverage once for overlapping regions of paired-end mates.
```
Genome       GGCCGGCATCGCATCAGCACATGCACACTGACACACACTGACTGGCTGCTGACTGACTGACTGCTGCTGCGCTATGCATGCCTGCTGAC
Read pair         GCATCGCATCAGCACATGCACACTGA>  
                                     <CACTGACACACACTGACTGGCTG
wig file          1111111111111111111111111111111111111111111
```

### RPM Normalization
The option `-n` enables normalization of the coverage using the reads per million mapped reads (RPM).

### Conversion to bigWig format
The ouput files in wiggle format can be converted to bigWig format using the [wigToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig) program from UCSC, which
also requires the chromosome sizes in an extra file `chrom.sizes`.

For example:
```
wigToBigWig coverage.plus.wig chrom.sizes coverage.plus.bw
```

###License

See the file LICENSE.


