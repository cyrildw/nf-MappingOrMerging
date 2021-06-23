# nf-MappingOrMerging
This workflow does the following : 
- In case of mapping : 
    1. Report the number of sequenced reads, trim reads, mapped reads and duplicates
    2. Trims the reads with trim_galore
    3. Mapps the reads with Bowtie or Subread (only one)
    4. Sort and index the mapping file (samtools)
    5. Creates a mapping file with nonduplicated reads (samtools)
    6. Creates bigwig files from bam files (deeptools bamCoverage)

- In case of bam file merging :
    1. Merges bam files (with gatk MergSamFiles)
    2. Sort and index the mapping file (samtools)
    3. Creates a mapping file with nonduplicated reads (samtools)
    4. Creates bigwig files from bam files (deeptools bamCoverage)


In case of mapping, this workflow only accepts paired-end sequencing files. This will hopefully be changed in further versions.

To run the workflow : 
nextflow run -c nf-MappingOrMerging.config --name 'BT_mapping' nf-MappingOrMerging.nf -bg --input_design "InputFiles.csv"

