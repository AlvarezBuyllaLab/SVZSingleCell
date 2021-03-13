
#Loading the original FastQ files was not possible due their huge size (R2 is 27gb, while the preprocessing is optimized for files under 2gb). We have to downsample our barcode fastqfiles to 100M reads:

gunzip TGACCAAT_S5_L001_R2_001.fastq.gz
gunzip TGACCAAT_S5_L001_R1_001.fastq.gz
seqtk sample -s100 TGACCAAT_S5_L001_R2_001.fastq 100000000 > TGACCAAT_S5_L001_R2_001_subsampled.fastq
seqtk sample -s100 TGACCAAT_S5_L001_R1_001.fastq 100000000 > TGACCAAT_S5_L001_R1_001_subsampled.fastq
gzip TGACCAAT_S5_L001_R2_001_subsampled.fastq
gzip TGACCAAT_S5_L001_R1_001_subsampled.fastq

gunzip ACAGTGAT_S6_L001_R1_001.fastq.gz
gunzip ACAGTGAT_S6_L001_R2_001.fastq.gz
seqtk sample -s100 ACAGTGAT_S6_L001_R1_001.fastq 100000000 > ACAGTGAT_S6_L001_R1_001_subsampled.fastq
seqtk sample -s100 ACAGTGAT_S6_L001_R2_001.fastq 100000000 > ACAGTGAT_S6_L001_R2_001_subsampled.fastq
gzip ACAGTGAT_S6_L001_R2_001_subsampled.fastq
gzip ACAGTGAT_S6_L001_R1_001_subsampled.fastq
