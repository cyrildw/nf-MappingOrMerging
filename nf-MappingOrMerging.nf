/* TODO :
1] usage
   - option names & descriptions
   - InputFiles.csv : 
      * add more columns ?
      * add a description in usage.
2] reporting nb mapped reads etc.
   - multiqc configuration (software not yet in the container)
   - use samtools stats ?
3] Think about how to integrate the spike in normalization.
   - add a spike_in_genome option
   - first map on the spike in genome
   - map on the normal genome
   - calculate SpikeIn normalization ratio
4] Add more in the configuration file
5] Think, if necessary ? how to automatically configure for CNRS vs CEA cluster, for file localization.

Merging bam files : 
nohup nextflow run -c nf-MappingOrMerging.config --name 'L10_L14_Merge' nf-MappingOrMerging.nf -with-trace -resume -bg --merge_bam --input_design "L10_L14_toMerge.csv"

PERFORMS ONLY ONE MAPPING : BOWTIE as default
Mapping with bowtie : 
nohup nextflow run -c nf-MappingOrMerging.config --name 'Toto_BT_only' nf-MappingOrMerging.nf -with-trace -resume -bg --input_design "InputFiles.csv"

Mapping with subread : 
nohup nextflow run -c nf-MappingOrMerging.config --name 'Toto_SR_only' nf-MappingOrMerging.nf -with-trace -resume -bg --bowtie_mapping false --subread_mapping true --input_design "InputFiles.csv"
*/


if (!params.merge_bam){
   if(params.bowtie_mapping){
      params.subread_mapping = false
      params.mapper_id = params.bowtie_id
   }
   else if(params.subread_mapping) { 
      params.bowtie_mapping = false
      params.mapper_id = params.subread_id
   }
   /* Reading design.csv file
    getting 3 values: LibName, LibFastq1, LibFastq2
    adding 4th value : prefix = LibName.mapper_id.genome_prefix.pe
   */

   Channel
      .fromPath(params.input_design)
      .splitCsv(header:true, sep:';')
      .map { row -> [ row.LibName,  file("$params.input_dir/$row.LibFastq1", checkIfExists: true), file("$params.input_dir/$row.LibFastq2", checkIfExists: true), "$row.LibName.${params.mapper_id}.${params.genome_prefix}"+".pe" ] }
      .into { design_reads_csv; ch_Toreport_reads_nb }
   
   /* Reporting at multiple steps
   - 1 nb of sequenced reads
   - 2 nb of reads post trimming
   - 3 nb of mapped reads w/ filter
   - 4 nb of non pcr duplicate reads
   - Merge all and write in csv
   */
   process _report_Nbseqreads {
      tag "$LibName"
      input:
      tuple LibName, file(LibFastq1), file(LibFastq2), MappingPrefix from ch_Toreport_reads_nb
      output:
      tuple LibName, MappingPrefix, stdout into ch_Toreport_trim

      script:
      """
      nb_line=`gunzip -dc ${LibFastq1} | wc -l`
      let nb_reads=\$nb_line/2
      echo -n \$nb_reads
      """
   }

if(params.bowtie_mapping){
   /*
   * Step 1. Trim the reads 
   *   - using trim_galore
   */
   process trimming {
      tag "$LibName"
      label "multiCpu"
      publishDir "${params.outdir}/TrimedReads", mode: 'copy', //params.publish_dir_mode,
      saveAs: { filename ->
               if (filename.endsWith('.fq.gz')) "./$filename"
               else null
      }

      input:
      tuple LibName, file(LibFastq1), file(LibFastq2), MappingPrefix from design_reads_csv
      output:
      tuple LibName, file("${LibName}_val_1.fq.gz"), file("${LibName}_val_2.fq.gz"), MappingPrefix into design_mapping_ch

      script:
      """
      trim_galore ${params.trim_galore_options} \
      --cores ${task.cpus} \
      --basename ${LibName} \
      ${LibFastq1} ${LibFastq2}
      """

   }

   /*
   * Step 2. Builds the genome index required by the mapping process
   TODO   - check if genome is already indexed
   TODO   - checkIfExists file basedir(params.genome)/params.genome_prefix.1.bt2,2.bt2, 3.bt2, 4.bt2, rev.1.bt2, rev.2.bt2
   TODO    - see https://github.com/SciLifeLab/NGI-smRNAseq/blob/master/main.nf
   */
   process buildIndexBT {
      tag "$genome.baseName"
      label "multiCpu"
      input:
      path genome from params.genome
         
      output:
      path 'genome.index*' into index_ch
         
      """
      bowtie2-build --threads ${task.cpus} ${genome} genome.index
      """
   }
   process mapping_Bowtie2 {
      echo true
      tag "$LibName"
      label 'multiCpu'
      input:
      tuple LibName, file(LibFastq1), file(LibFastq2), MappingPrefix, SequencedReads from design_mapping_ch
      path genome from params.genome
      file index from index_ch

      output:
      tuple val(LibName), val(MappingPrefix), val(SequencedReads), file("${MappingPrefix}.bam") into mapping_ch
      val LibName into libName_ch

      """
      bowtie2 \
      --very-sensitive \
      --threads ${task.cpus} \
      -x genome.index \
      -1 ${LibFastq1} \
      -2 ${LibFastq2} | samtools view -bSh -F 4 -f 3 -q ${params.samtools_q_filter} - > ${MappingPrefix}.bam
      """
      //-S ${MappingPrefix}.raw.sam 
      /*
      f 3 includes mapped reads and properly paired
      F 4 excludes unmapped reads. [F 256 excludes secondary alignments]-> not used anymore
      */

   }
}
else if(params.subread_mapping){
   process buildIndexSR {
      tag "$genome.baseName"
      input:
      path genome from params.genome


      output:
      path 'genome.index*' into index_ch
      file("log")
      """
      subread-buildindex -o genome.index ${genome} 2>log 1>>log
      """
   }
   process mapping_Subread {
      echo true
      tag "$LibName"
      label 'multiCpu_short'
      maxForks 8
      input:
      tuple LibName, file(LibFastq1), file(LibFastq2), MappingPrefix, SequencedReads from design_mapping_ch
      path genome from params.genome
      file index from index_ch

      output:
      tuple val(LibName), val(MappingPrefix), val(SequencedReads), file("${MappingPrefix}.bam") into mapping_ch
      val LibName into libName_ch
      file("log")
      file("${MappingPrefix}.tmp.bam*")
      // --sv is raising some issues with L15_28_F_UNG1D_B150 for example. 
      // add && rm ${MappingPrefix}.tmp.bam when it's working
      """
      subread-align \
      -t 1 \
      -T ${task.cpus} \
      -i genome.index \
      -r ${LibFastq1} \
      -R ${LibFastq2} \
      --multiMapping -B 3 \
      -o ${MappingPrefix}.tmp.bam &>log && samtools view -bh -F 4 -f 3 -q ${params.samtools_q_filter} ${MappingPrefix}.tmp.bam  > ${MappingPrefix}.bam
      """

      //-S ${MappingPrefix}.raw.sam 
      /*
      f 3 includes mapped reads and properly paired
      F 4 excludes unmapped reads. [F 256 excludes secondary alignments]-> not used anymore
      */

   }
}   

}
else {
   Channel
      .fromPath(params.input_design)
      .splitCsv(header:true, sep:';')
      .map { row -> [ row.LibName, row.LibExp, file("$params.input_dir/$row.LibBam", checkIfExists: true), file("$params.input_dir/${row.LibBam}.bai", checkIfExists: true), "$row.LibName.${params.mapper_id}.${params.genome_prefix}"+".merged" ] }
      .into { design_bam_csv; test_design }
      test_design.view()

   design_bam_csv
      .map{ it -> [it[1], it[2] ]}
      .groupTuple(by: 0)
      .map { it-> [ it[0], it[1].flatten()]}
      .set {design_bam_merged}

   /* Merge bam from same experiment
   *      
   */

   process mergeBamFiles {
      echo true
      tag "$LibExp"
      label 'usePicard'
      input:
      tuple val(LibExp), path(bams) from design_bam_merged
      output:
      tuple val(LibExp), val(MappingPrefix), val("NA"), file("${MappingPrefix}.bam") into mapping_ch
      script:
      MappingPrefix="${LibExp}.${params.mapper_id}.${params.genome_prefix}.pe.merged"
      bam_files=bams
      //bam_files=bams.findAll { it.toString().endsWith('.bam')}.sort()
      """
      gatk MergeSamFiles ${'-I='+bam_files.join(' -I=')} -O=${MappingPrefix}.bam --TMP_DIR=tmp
      """
   }

}



/*
 * Step 2. Map paired end reads on the genome
 *   ONGOING - couple with a samtools view direclty to avoid high weight file.
 *    *    configure optionnal parameters for filtering
  */

/*
 * Step 3. Filters the mapping file with samtools
 TODO    add flagstat
 TODO    add idxstats
 */
process samtools {
   echo true
   tag "$LibName"
   

   publishDir "${params.outdir}/Mapping", mode: 'copy' //params.publish_dir_mode,

   input:
   tuple val(LibName), val(prefix), val(SequencedReads), file(RawMapping) from mapping_ch
   output:
   file("${prefix}.*.bam*") 
   tuple val(LibName), val(prefix), val(SequencedReads), file("${prefix}.sorted.bam*") into samtooled_ch
   tuple val(LibName), val(prefix), val(SequencedReads), file("${prefix}.sorted.rmdup.bam*") into samtooled_rmdup_ch

	//samtools view -bSh -F 4 -f 3 -q ${params.samtools_q_filter} ${prefix}.raw.sam > ${prefix}.bam
   script:
	"""
   samtools sort ${prefix}.bam ${prefix}.sorted
   samtools index ${prefix}.sorted.bam

   samtools rmdup ${prefix}.sorted.bam ${prefix}.sorted.rmdup.bam
   samtools index ${prefix}.sorted.rmdup.bam
   """
}   

/*
 * Step 4. Get the number of mapped reads.
 *    Parallelized with two different processes
 *       - for .bam           -for .rmdup.bam
 */
process nb_mapped_reads {
	tag "$LibName .bam"
	input:
	tuple val(LibName), val(prefix), val(SequencedReads), path(bamFiles) from samtooled_ch
	output:
	tuple val(LibName), val(prefix), val(SequencedReads), stdout, path(bamFiles) into forGenomeCov_ch

	script:
	"""
	mapped_reads=`samtools view -c ${bamFiles[0]}`
	echo -n \$mapped_reads
	"""
}

process nb_mapped_reads_rmdup {
	tag "$LibName rmdup.bam"
	input:
	tuple val(LibName), val(prefix), val(SequencedReads), path(bamFiles) from samtooled_rmdup_ch
	output:
	tuple val(LibName), val(prefix), val(SequencedReads), stdout, path(bamFiles) into forGenomeCov_rmdup_ch

	script:
	"""
	mapped_reads=`samtools view -c ${bamFiles[0]}`
	echo -n \$mapped_reads
	"""
}

//process nb_reads_concat {
// Concatenate all nb_reads stats into one file
// requires nb_reads output
// collect()
// output a single file.
//}

process genome_coverage_bam {
   tag "$LibName genome coverage .bam"
   label 'multiCpu'
   publishDir "${params.outdir}/GenomeCoverage", mode: 'copy', //params.publish_dir_mode,
      saveAs: { filename ->
               if (filename.endsWith('.bw')) "./$filename"
               else null
      }
   input:
  	tuple val(LibName), val(prefix), val(SequencedReads), val(MappedReads), path(bamFiles) from forGenomeCov_ch
   output:
	tuple val(LibName), val(prefix), val(SequencedReads), val(MappedReads), path(bamFiles), path("${prefix}.bin${params.bin_size}.RPM.bamCoverage.bw") into genCoved_ch

   """
   bamCoverage \
   -b ${bamFiles[0]} \
   -o ${prefix}.bin${params.bin_size}.RPM.bamCoverage.bw -of bigwig \
   --extendReads --centerReads --samFlagInclude 3 \
   --scaleFactor 1 --normalizeUsing CPM  --exactScaling \
   --binSize ${params.bin_size} -p ${task.cpus}
   """
}
//--effectiveGenomeSize 12157105 

process genome_coverage_rmdup {
   tag "$LibName genome coverage rmdup.bam"
   label 'multiCpu'
   publishDir "${params.outdir}/GenomeCoverage", mode: 'copy', //params.publish_dir_mode,
      saveAs: { filename ->
               if (filename.endsWith('.bw')) "./$filename"
               else null
      }
   input:
  	tuple val(LibName), val(prefix), val(SequencedReads), val(MappedReads), path(bamFiles) from forGenomeCov_rmdup_ch
   output:
	tuple val(LibName), val(prefix), val(SequencedReads), val(MappedReads), path(bamFiles), path("${prefix}.bin${params.bin_size}.RPM.rmdup.bamCoverage.bw") into genCoved_rmdup_ch

   """
   bamCoverage \
   -b ${bamFiles[0]} \
   -o ${prefix}.bin${params.bin_size}.RPM.rmdup.bamCoverage.bw -of bigwig \
   --extendReads --centerReads --samFlagInclude 3 \
   --scaleFactor 1 --normalizeUsing CPM  --exactScaling \
   --binSize ${params.bin_size} -p ${task.cpus}
   """
}

process report_stats {
   tag "$LibName .bam"
   publishDir "${params.outdir}/Stats", mode: 'copy'
   input:
   tuple val(LibName), val(prefix), val(SequencedReads), val(MappedReads), path(bamFiles), path(bwFile) from genCoved_ch
   output:
   path("${LibName}.mapped_reads.txt")
   script:
   """
   echo "${LibName};${SequencedReads};${MappedReads};" > ${LibName}.mapped_reads.txt
   """
}
process report_stats_rmdup {
   tag "$LibName .rmdup.bam"
   publishDir "${params.outdir}/Stats", mode: 'copy'
   input:
   tuple val(LibName), val(prefix), val(SequencedReads), val(MappedReads), path(bamFiles), path(bwFile) from genCoved_rmdup_ch
   output:
   path("${LibName}.rmdup.mapped_reads.txt")
   script:
   """
   echo "${LibName};${SequencedReads};${MappedReads};" > ${LibName}.rmdup.mapped_reads.txt
   """
}
/* 
TODO  - Channel.collectFile to output what's needed for nf-AnalysesOnCoordinates (BigwigDesign.csv)