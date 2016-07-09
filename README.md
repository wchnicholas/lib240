##Amplicon sequencing 240 bp single nucleotide mutant library
* script/Mapper1.py and script/Mapper2.py: Have to run sequentially to generate count data from raw reads. Make sure the fastq files are not zipped.
* Fasta/Barcode: For demultiplexing the samples (first 3 bp of both forward read and reverse read
* Fasta/flu1amp.fa: Reference sequence for flu1 (PB2 segment)
* Fasta/flu1offset: position offset of each amplicon
* Fasta/fluTranslate: Translation info for point mutations
