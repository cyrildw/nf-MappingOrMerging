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
    adding the index to keep track of the input order. This value will only be kept in the reporting channels
    adding 5th value : prefix = LibName.mapper_id.genome_prefix.pe
   */
i=0;
   Channel
      .fromPath(params.input_design)
      .splitCsv(header:true, sep:';')
      .map { row -> [ row.LibName, i++, file("$params.input_dir/$row.LibFastq1", checkIfExists: true), file("$params.input_dir/$row.LibFastq2", checkIfExists: true), "$row.LibName.${params.mapper_id}.${params.genome_prefix}"+".pe" ] }
      .into { design_reads_csv; ch_Toreport_reads_nb }
   
   /* Reporting at multiple steps
   - 1 nb of sequenced reads (_report_Nbseqreads)
   - 2 nb of reads post trimming (_report_Nbtrimreads)
   - 3 nb of mapped reads w/ filter (_report_Nbmappedreads)
   - 4 nb of non pcr duplicate reads (_report_Nbuniqreads)
   - 5 insert size
   - Merge all and write in csv
   */
   process _report_Nbseqreads {
      tag "$LibName"
      input:
      tuple val(LibName), LibIdx, file(LibFastq1), file(LibFastq2), MappingPrefix from ch_Toreport_reads_nb
      output:
      tuple val(LibName), LibIdx, stdout into ( ch_Toreport_trim_nb, test1_ch )
      script:
      """
      nb_line1=`gunzip -dc ${LibFastq1} | wc -l`
      nb_line2=`gunzip -dc ${LibFastq2} | wc -l`
      let nb_reads1=\$nb_line1/4
      let nb_reads2=\$nb_line2/4
      let nb_reads=\$nb_reads1+\$nb_reads2
      echo -n \$nb_reads
      """
   }
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
      tuple val(LibName), val(LibIdx), file(LibFastq1), file(LibFastq2), MappingPrefix from design_reads_csv
      output:
      tuple val(LibName), file("${LibName}_val_1.fq.gz"), file("${LibName}_val_2.fq.gz"), MappingPrefix into design_mapping_ch 
      tuple val(LibName), file("${LibName}_val_1.fq.gz"), file("${LibName}_val_2.fq.gz") into trimed_reads_ch //don't need the LibIdx since it's already in the ch_Toreport_trim_nb
      script:
      """
      trim_galore ${params.trim_galore_options} \
      --cores ${task.cpus} \
      --basename ${LibName} \
      ${LibFastq1} ${LibFastq2}
      """
   }

ch_Toreport_trim_nb.join(trimed_reads_ch)
   .set{ch_report_trim_nb}

   process _report_Nbtrimreads {
      tag "$LibName"
      input:
      tuple val(LibName), val(LibIdx),  NbSeqReads, file(LibFastq1), file(LibFastq2) from ch_report_trim_nb
      output:
      tuple val(LibName), val(LibIdx),  NbSeqReads, stdout into (ch_Toreport_mapped_nb, ch_Toreport_uniq_nb)

      script:
      """
      nb_line1=`gunzip -dc ${LibFastq1} | wc -l`
      nb_line2=`gunzip -dc ${LibFastq2} | wc -l`
      let nb_reads1=\$nb_line1/4
      let nb_reads2=\$nb_line2/4
      let nb_reads=\$nb_reads1+\$nb_reads2
      echo -n \$nb_reads
      """
   }


if(params.bowtie_mapping){
   
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

   /*
   * Step 3. Mapping */
   process mapping_Bowtie2 {
      echo true
      tag "$LibName"
      label 'multiCpu'
      input:
      tuple val(LibName), file(LibFastq1), file(LibFastq2), MappingPrefix from design_mapping_ch
      path genome from params.genome
      file index from index_ch

      output:
      tuple val(LibName),  val(MappingPrefix), file("${MappingPrefix}.bam") into mapping_ch
      val LibName into libName_ch

      """
      bowtie2 \
      --very-sensitive \
      --threads ${task.cpus} \
      -x genome.index \
      -1 ${LibFastq1} \
      -2 ${LibFastq2} 2>/dev/null | samtools view -bSh -F 4 -f 3 -q ${params.samtools_q_filter} - > ${MappingPrefix}.bam 
      """
      //-S ${MappingPrefix}.raw.sam 
      /*
      f 3 includes mapped reads and properly paired
      F 4 excludes unmapped reads. [F 256 excludes secondary alignments]-> not used anymore
      */

   }
}
else if(params.subread_mapping){
/* Step 2 indexing genome */
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
      /*
   * Step 3. Mapping */
   process mapping_Subread {
      echo true
      tag "$LibName"
      label 'multiCpu_short'
      maxForks 8
      input:
      tuple val(LibName), file(LibFastq1), file(LibFastq2), MappingPrefix from design_mapping_ch
      path genome from params.genome
      file index from index_ch

      output:
      tuple val(LibName), val(MappingPrefix), file("${MappingPrefix}.bam") into mapping_ch
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
   i=0;
   Channel
      .fromPath(params.input_design)
      .splitCsv(header:true, sep:';')
      .map { row -> [ row.LibName, row.LibExp,i++, file("$params.input_dir/$row.LibBam", checkIfExists: true), file("$params.input_dir/${row.LibBam}.bai", checkIfExists: true), "$row.LibName.${params.mapper_id}.${params.genome_prefix}"+".merged" ] }
      .into { design_bam_csv; test_design }
      test_design.view()

   design_bam_csv
      .map{ it -> [it[1],it[2], it[3] ]}
      .groupTuple(by: 0)
      .map { it-> [ it[0], it[1], it[3].flatten()]}
      .set {design_bam_merged}

   /* Merge bam from same experiment
   *      
   */

   process mergeBamFiles {
      echo true
      tag "$LibExp"
      label 'usePicard'
      input:
      tuple val(LibExp), val(LibIdx), path(bams) from design_bam_merged
      output:
      tuple val(LibExp), val(MappingPrefix), file("${MappingPrefix}.bam") into mapping_ch
      tuple val(LibExp), val(LibIdx), val("NA") , val("NA") into (ch_Toreport_mapped_nb, ch_Toreport_uniq_nb) // creating the report channels with unavailable data for nbseq & nb_trim
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
 * Step 4. Filters the mapping file with samtools
 TODO    add flagstat
 TODO    add idxstats
 */
process samtools {
   echo true
   tag "$LibName"
   publishDir "${params.outdir}/Mapping", mode: 'copy' //params.publish_dir_mode,

   input:
   tuple val(LibName),  val(prefix),  file(RawMapping) from mapping_ch
   output:
   file("${prefix}.*.bam*") 
   tuple val(LibName), val(prefix), file("${prefix}.sorted.bam*") into samtooled_ch
   tuple val(LibName), file("${prefix}.sorted.bam*") into mapped_reads_ch
   tuple val(LibName), val(prefix),  file("${prefix}.sorted.rmdup.bam*") into samtooled_rmdup_ch
   tuple val(LibName), file("${prefix}.sorted.rmdup.bam*") into mapped_uniq_reads_ch

	//samtools view -bSh -F 4 -f 3 -q ${params.samtools_q_filter} ${prefix}.raw.sam > ${prefix}.bam
   script:
	"""
   samtools sort ${prefix}.bam ${prefix}.sorted
   samtools index ${prefix}.sorted.bam

   samtools rmdup ${prefix}.sorted.bam ${prefix}.sorted.rmdup.bam
   samtools index ${prefix}.sorted.rmdup.bam
   """
}   




process genome_coverage_bam {
   tag "$LibName genome coverage .bam"
   label 'multiCpu'
   publishDir "${params.outdir}/GenomeCoverage", mode: 'copy', //params.publish_dir_mode,
      saveAs: { filename ->
               if (filename.endsWith('.bw')) "./$filename"
               else null
      }
   input:
  	tuple val(LibName), val(prefix), path(bamFiles) from samtooled_ch
   output:
	tuple val(LibName), val(prefix), bamFiles, val("${prefix}.bin${params.bin_size}.RPM.bamCoverage.bw") into genCoved_ch

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
  	tuple val(LibName), val(prefix), path(bamFiles) from samtooled_rmdup_ch
   output:
	tuple val(LibName), val(prefix), bamFiles, val("${prefix}.bin${params.bin_size}.RPM.rmdup.bamCoverage.bw") into genCoved_uniq_ch

   """
   bamCoverage \
   -b ${bamFiles[0]} \
   -o ${prefix}.bin${params.bin_size}.RPM.rmdup.bamCoverage.bw -of bigwig \
   --extendReads --centerReads --samFlagInclude 3 \
   --scaleFactor 1 --normalizeUsing CPM  --exactScaling \
   --binSize ${params.bin_size} -p ${task.cpus}
   """
}

/*
 * Step 5. Get the number of mapped reads.
 *    Parallelized with two different processes
 *       - for .bam           -for .rmdup.bam
 */
ch_Toreport_mapped_nb
   .join(mapped_reads_ch)
   .set{ch_report_mapped_nb}
   

process _report_nb_mapped_reads {
	tag "$LibName "
	input:
	tuple val(LibName), val(LibIdx),  val(NbSeqReads), val(NbTrimReads), path(bamFiles) from ch_report_mapped_nb
	output:
	tuple val(LibName), val(LibIdx),  val(NbSeqReads), val(NbTrimReads), stdout, path(bamFiles) into ch_Toreport_insert_size

	script:
	"""
	mapped_reads=`samtools view -c ${bamFiles[0]}`
	echo -n \$mapped_reads
	"""
}

process _report_insert_size {
   tag "$LibName"
   input:
   tuple val(LibName), val(LibIdx),  val(NbSeqReads), val(NbTrimReads), val(NbMapReads), path(bamFiles) from ch_Toreport_insert_size
   output:
   tuple val(LibName), val(LibIdx),  val(NbSeqReads), val(NbTrimReads), val(NbMapReads), stdout into (ch_Toreport_all_stats, ch_ToAoC)
   file(table)
   script:
   """
   bamPEFragmentSize --bamfiles ${bamFiles[0]} --table table >/dev/null 2>&1
   ins_size=`tail -1 table | awk '{ print \$6}'`
   echo -n \$ins_size
   """
}

ch_Toreport_all_stats.view()
.toSortedList( { a, b -> a[1] <=> b[1] })//sorting by input order
.map{ it -> [it[0],it[2],it[3],it[4],it[5] ]}//removing the index number used for sorting
.map{it -> [it.join(";")]}.collect().set{ ch_report_all_stats} //Joining stats with ";" then use collect to have a single entry channel

process _report_mapping_stats_csv {
   publishDir "${params.outdir}/Stats", mode: 'copy'
   input:
   val x from ch_report_all_stats
   output:
   path("mapping_stats.txt")
   // echoing all the channel with join('\n') into the mapping_stats.txt file
   script:
   """
   echo "LibName;Nb_sequenced_read;Nb_trimmed_reads;Nb_mapped_reads;Median_insert_size" > mapping_stats.txt
   echo "${x.join('\n')}" >> mapping_stats.txt
   """
}



ch_Toreport_uniq_nb
   .join(mapped_uniq_reads_ch)
   .set{ch_report_uniq_nb}

process _report_nb_uniq_reads {
	tag "$LibName rmdup.bam"
	input:
	tuple val(LibName), val(LibIdx),  val(NbSeqReads), val(NbTrimReads), path(bamFiles) from ch_report_uniq_nb
	output:
	tuple val(LibName), val(LibIdx),  val(NbSeqReads), val(NbTrimReads), stdout, path(bamFiles) into ch_Toreport_uniq_insert_size
	script:
	"""
	mapped_reads=`samtools view -c ${bamFiles[0]}`
	echo -n \$mapped_reads
	"""
}
process _report_uniq_insert_size {
   tag "$LibName"
   input:
   tuple val(LibName), val(LibIdx),  val(NbSeqReads), val(NbTrimReads), val(NbMapReads), path(bamFiles) from ch_Toreport_uniq_insert_size
   output:
   tuple val(LibName), val(LibIdx),  val(NbSeqReads), val(NbTrimReads), val(NbMapReads), stdout into (ch_Toreport_uniq_stats, ch_ToAoC_uniq)
   file(table_uniq)
   script:
   """
   bamPEFragmentSize --bamfiles ${bamFiles[0]} --table table_uniq >/dev/null 2>&1
   ins_size_uniq=`tail -1 table_uniq | awk '{ print \$6}'`
   echo -n \$ins_size_uniq
   """
}

ch_Toreport_uniq_stats.view()
.toSortedList( { a, b -> a[1] <=> b[1]}) //sorting by input order
.map{ it -> [it[0],it[2],it[3],it[4],it[5] ]} //removing the index number used for sorting
.map{it -> [it.join(";")]}.collect().set{ ch_report_uniq_stats} //Joining stats with ";" then use collect to have a single entry channel

process _report_mapping_uniq_stats_csv {
   publishDir "${params.outdir}/Stats", mode: 'copy'
   input:
   val x from ch_report_all_stats
   output:
   path("mapping_uniq_stats.txt")
   // echoing all the channel with join('\n') into the mapping_stats.txt file
   script:
   """
   echo "LibName;Nb_sequenced_read;Nb_trimmed_reads;Nb_mapped_reads;Median_insert_size" > mapping_uniq_stats.txt
   echo "${x.join('\n')}" >> mapping_uniq_stats.txt
   """
}


genCoved_ch.join(ch_ToAoC, by:0).view()
.toSortedList( { a, b -> a[4] <=> b[4] }) //sorting by input order
.map{ it -> [it[0], it[2][0], it[3], it[5], it[7], 'NA',  1, '', '', '', '', '', '', '']}
.map{ it -> [it.join(";")]}
.collect()
.set {ch_report_Aoc}

genCoved_uniq_ch.join(ch_ToAoC_uniq, by:0).view()
.toSortedList( { a, b -> a[4] <=> b[4]}) //sorting by input order
.map{ it -> [it[0], it[2][0], it[3], it[5], it[7], it[7],  1, '', '', '', '', '', '', '']}
.map{ it -> [it.join(";")]}
.collect()
.set {ch_report_Aoc_uniq}

process _report_AoC_csv {
   publishDir "${params.outdir}", mode: 'copy'
   input:
   val x from ch_report_Aoc
   output:
   path("${params.name}.bigwigDesign.csv")
   // echoing all the channel with join('\n') into the "${params.name}.bigwigDesign.csv" file
   script:
   """
   echo "LibName;LibBam;LibBW;LibSequenced;LibMapped;LibUnique;LibInsertSize;LibQpcrNorm;LibType;LibProj;LibExp;LibCondition;LibOrder;LibIsControl;LibControl" > ${params.name}.bigwigDesign.csv
   echo "${x.join('\n')}" >> ${params.name}.bigwigDesign.csv
   """
}

process _report_AoC_uniq_csv {
   publishDir "${params.outdir}", mode: 'copy'
   input:
   val x from ch_report_Aoc_uniq
   output:
   path("${params.name}.rmdup.bigwigDesign.csv")
   // echoing all the channel with join('\n') into the "${params.name}.rmdup.bigwigDesign.csv" file
   script:
   """
   echo "LibName;LibBam;LibBW;LibSequenced;LibMapped;LibUnique;LibInsertSize;LibQpcrNorm;LibType;LibProj;LibExp;LibCondition;LibOrder;LibIsControl;LibControl" > ${params.name}.rmdup.bigwigDesign.csv
   echo "${x.join('\n')}" >> ${params.name}.rmdup.bigwigDesign.csv
   """
}
