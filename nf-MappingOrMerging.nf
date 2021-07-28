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
   - first map on the concatenation of ref_genome and spike_in_genome
   - extract and count reads
   - calculate SpikeIn normalization ratio for the bamCoverage
4] Allow single end mapping and or multiple fastq pairs (Fq1a Fq1b Fq2a Fq2b)
   - when reading the csv file, the fq colomn is no longer a file but a list of file (, separated)
   - adapt read counting, mapping, samtools filtering (in case of single end)
   - reject hybrid setups with 1 paired end and 1 single end sequencing for the same sample
   - for single end sequencing : estimate insert size for bamCoverage
TODO - samtools options (-F 4 -f 3) to include in the merge options
OK   - bamCoverage options (indluding bins)
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
    adding 5th value : prefix = LibName.mapper_id.ref_genome_prefix.pe
   */
i=0;
   Channel
      .fromPath(params.input_design)
      .splitCsv(header:true, sep:';')
      .map { row -> [ row.LibName, i++, file("$params.input_dir/$row.LibFastq1", checkIfExists: true), file("$params.input_dir/$row.LibFastq2", checkIfExists: true), "$row.LibName.${params.mapper_id}.${params.ref_genome_prefix}"+".pe" ] }
      .into { design_reads_csv; ch_Toreport_reads_nb }


/*
#############################################################
#
#     READ PROCESSING : 
# 1. Counting reads in the fastq files
# 2. Triming reads with Trim Galore
# 3. Counting trimed reads.
#
#
##############################################################
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


/*
#############################################################
#
#     READ MAPPING : 
# 1. if Spike in : 
#        + from individual reference genome : get the seq_ids
#        + Create the concatenated reference'
# 2. Indexing the reference
# 3. Read mapping with the choosen mapper
#
#
##############################################################
*/   

/*TODO : 
- If spike_in_norm = true : make genome concatenation => channel for genome ref and value for concatenated genome prefix
- Create a new channel for the reference mapping_ref_ch
   * if(!spike_in_norm){mapping_ref_ch = ${params.genome_ref} }
*/
if(params.spike_in_norm){
   process hybrid_genome {
      // Concatenate the two genome (ref & si) and extract names of the seq_ids in both files.
      label 'noContainer'
      input:
      path ref_genome from params.ref_genome
      path spike_in_genome from params.spike_in_genome
      output:
      path "${params.ref_genome_prefix}_${params.spike_in_genome_prefix}.fa" into ref_to_index_ch
      path "${params.ref_genome_prefix}_${params.spike_in_genome_prefix}.fa" into mapping_ref_ch
      path "${params.ref_genome_prefix}.seq_ids.txt" into ref_seq_id_File_ch 
      path "${params.spike_in_genome_prefix}.seq_ids.txt" into spike_in_seq_id_File_ch
      
      """
      cat ${ref_genome} ${spike_in_genome} > ${params.ref_genome_prefix}_${params.spike_in_genome_prefix}.fa
      grep ">" ${ref_genome} | sed 's/>\\([^[:blank:]]*\\)[[:blank:]]*.*/\\1/g' > ${params.ref_genome_prefix}.seq_ids.txt
      grep ">" ${spike_in_genome} | sed 's/>\\([^[:blank:]]*\\)[[:blank:]]*.*/\\1/g' > ${params.spike_in_genome_prefix}.seq_ids.txt
      """
   }
   process ref_seq_id_parsing {
      label 'noContainer'
      input:
      path myFile from ref_seq_id_File_ch
      output:
      stdout into (ref_genome_seq_id_ch, ref_genome_seq_id_4uniq_ch)
      
      """
      cat ${myFile}
      """
   }
   process spike_in_seq_id_parsing {
      label 'noContainer'
      input:
      path myFile from spike_in_seq_id_File_ch
      output:
      stdout into (spike_in_genome_seq_id_ch, spike_in_genome_seq_id_4uniq_ch)
      
      """
      cat ${myFile}
      """
   }

   
}
else if(!params.spike_in_norm){
   //In absence of spike in the mapping_ref_ch contains the ref_genome
   mapping_ref_ch = Channel.fromPath( params.ref_genome)
   ref_to_index_ch = Channel.fromPath( params.ref_genome)
}
if(params.bowtie_mapping){
   
   /*
   * Step 2. Builds the genome index required by the mapping process
   TODO   - check if genome is already indexed
   TODO   - checkIfExists file basedir(params.genome)/params.ref_genome_prefix.1.bt2,2.bt2, 3.bt2, 4.bt2, rev.1.bt2, rev.2.bt2
   TODO    - see https://github.com/SciLifeLab/NGI-smRNAseq/blob/master/main.nf
   */
   process buildIndexBT {
      tag "$genome.baseName"
      label "multiCpu"
      input:
      path genome from ref_to_index_ch
         
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
      path genome from mapping_ref_ch
      file index from index_ch

      output:
      tuple val(LibName),  val(MappingPrefix), file("${MappingPrefix}.bam") into mapping_ch
      val LibName into libName_ch

      """
      bowtie2 \
      ${params.bowtie_options} \
      --threads ${task.cpus} \
      -x genome.index \
      -1 ${LibFastq1} \
      -2 ${LibFastq2} 2>/dev/null | samtools view -bSh ${params.samtools_flag_filter} -q ${params.samtools_q_filter} - > ${MappingPrefix}.bam 
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
      path genome from ref_to_index_ch


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
      path genome from mapping_ref_ch
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
      -t 1 ${params.subread_options} \
      -T ${task.cpus} \
      -i genome.index \
      -r ${LibFastq1} \
      -R ${LibFastq2} \
      -o ${MappingPrefix}.tmp.bam &>log && samtools view -bh ${params.samtools_flag_filter} -q ${params.samtools_q_filter} ${MappingPrefix}.tmp.bam  > ${MappingPrefix}.bam
      """

      //-S ${MappingPrefix}.raw.sam 
      /*
      f 3 includes mapped reads and properly paired
      F 4 excludes unmapped reads. [F 256 excludes secondary alignments]-> not used anymore
      */

   }
}

}

/*
#############################################################
#
#     MAPPING FILE MERGING : 
# 1. Importing bam files
# 2. Use MergeSamFiles for merging
# TODO Filtering with samtools view ?
#
#
##############################################################
*/   

else {
   i=0; //Counter used to keep the  input files order
   Channel
      .fromPath(params.input_design)
      .splitCsv(header:true, sep:';')
      .map { row -> [ row.LibName, row.LibExp,i++, file("$params.input_dir/$row.LibBam", checkIfExists: true), file("$params.input_dir/${row.LibBam}.bai", checkIfExists: true), "$row.LibName.${params.mapper_id}.${params.ref_genome_prefix}"+".merged" ] }
      .into { design_bam_csv; test_design }
      test_design.view()

   design_bam_csv
      .map{ it -> [it[1],it[2], it[3] ]}
      .groupTuple(by: 0)
      .map { it-> [ it[0], it[1], it[3].flatten()]}
      .set {design_bam_merged}

   /* Merge bam from same experiment
   *      TODO : perform samtools view with filtering parameters
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
      MappingPrefix="${LibExp}.${params.mapper_id}.${params.ref_genome_prefix}.pe.merged"
      bam_files=bams
      //bam_files=bams.findAll { it.toString().endsWith('.bam')}.sort()
      """
      gatk MergeSamFiles ${'-I='+bam_files.join(' -I=')} -O=${MappingPrefix}.bam --TMP_DIR=tmp
      """
   }

}



/*
#############################################################
#
#     MAPPING file filtering : 
# 1. Samtools index
# 2. if Spike in : 
      + split the mapping file according to both genomes
      + rehead the bam file with correct header (only reference seq_ids)
      + count the mapped reads on each genomes
# 2. Indexing the reference
# 3. Read mapping with the choosen mapper
#
#
##############################################################
*/   
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

if(params.spike_in_norm){
   /* From samtools process : 
         1. Count reads on the concatenated genome
            1.1. if nb reads mapped on the spike_in_genome < spike_in_min_read_nb : set the channels to the output of samtools process.
         2. Split the alignments based on the genome
         3. Get the normalization factor.
   */
   process si_mapping_split{
      input:
      tuple val(LibName), val(prefix), path(bamFiles) from samtooled_ch
      val(ref_seq_ids) from ref_genome_seq_id_ch.splitText().map{ it.replaceAll("\n", "")}.collect()
      val(si_seq_ids) from spike_in_genome_seq_id_ch.splitText().map{ it.replaceAll("\n", "")}.collect()
      output:
      tuple val(LibName), val(prefix), file("${prefix}.split_ref.sorted.bam*"), stdout into to_bamCov_ch
      tuple val(LibName), file("${prefix}.sorted.bam*") into to_count_mapped_reads_ch
      path "${prefix}.split_spike_in.sorted.bam*"
      path "tmp.bam"
      path "header.txt"
      path "header_ref.txt"
      path "header_spike_in.txt"
            
      // Getting the header & counting the total number of mapped reads.
      // Extracting reads mapping on the ref_genome, creating a correct header then rehead and index the mapping file and count mapped reads
      // Extracting reads mapping on the spike_in_genome, creating a correct header then go through a sam file intermediate (otherwise reheading bugs) and index the mapping file and count mapped reads
      // Calculating the NormFactor for bamCoverage
      """
      samtools view -H ${bamFiles[0]} > header.txt
      NB_TOTAL_MAPPED=`samtools view -c ${bamFiles[0]}`

      samtools view -bh -o tmp.bam ${bamFiles[0]} ${ref_seq_ids.join(' ')}
      grep -vP "${'SN:'+si_seq_ids.join('\\s|SN:')}" header.txt > header_ref.txt
      samtools reheader header_ref.txt tmp.bam > ${prefix}.split_ref.sorted.bam && samtools index ${prefix}.split_ref.sorted.bam && rm tmp.bam
      NB_REF_MAPPED=`samtools view -c ${prefix}.split_ref.sorted.bam`
      
      samtools view -o tmp.bam ${bamFiles[0]} ${si_seq_ids.join(' ')}
      grep -vP "${'SN:'+ref_seq_ids.join('\\s|SN:')}" header.txt > header_spike_in.txt
      cat header_spike_in.txt tmp.bam | samtools view -bSh - >  ${prefix}.split_spike_in.sorted.bam && samtools index ${prefix}.split_spike_in.sorted.bam && rm tmp.bam
      NB_SPIKE_IN_MAPPED=`samtools view -c ${prefix}.split_spike_in.sorted.bam`

      NORM_FACTOR=`R --slave -q -e "cat(round((10e6/$NB_TOTAL_MAPPED)*(${params.spike_in_fraction}/($NB_SPIKE_IN_MAPPED/$NB_TOTAL_MAPPED)),6), '\n')"`
      echo \$NORM_FACTOR
      """
   }
   process si_mapping_uniq_split{
   input:
   tuple val(LibName), val(prefix), path(bamFiles) from samtooled_rmdup_ch
   val(ref_seq_ids) from ref_genome_seq_id_4uniq_ch.splitText().map{ it.replaceAll("\n", "")}.collect()
   val(si_seq_ids) from spike_in_genome_seq_id_4uniq_ch.splitText().map{ it.replaceAll("\n", "")}.collect()
   output:
   tuple val(LibName), val(prefix), file("${prefix}.split_ref.sorted.rmdup.bam*"), stdout into to_bamCov_rmdup_ch
   tuple val(LibName), file("${prefix}.sorted.rmdup.bam*") into to_count_uniq_mapped_reads_ch
   path "${prefix}.split_spike_in.sorted.rmdup.bam*"
   path "tmp.bam"
   path "header.txt"
   path "header_ref.txt"
   path "header_spike_in.txt"
         
   // Getting the header & counting the total number of mapped reads.
   // Extracting reads mapping on the ref_genome, creating a correct header then rehead and index the mapping file and count mapped reads
   // Extracting reads mapping on the spike_in_genome, creating a correct header then go through a sam file intermediate (otherwise reheading bugs) and index the mapping file and count mapped reads
   // Calculating the NormFactor for bamCoverage
   """
   samtools view -H ${bamFiles[0]} > header.txt
   NB_TOTAL_MAPPED=`samtools view -c ${bamFiles[0]}`

   samtools view -bh -o tmp.bam ${bamFiles[0]} ${ref_seq_ids.join(' ')}
   grep -vP "${'SN:'+si_seq_ids.join('\\s|SN:')}" header.txt > header_ref.txt
   samtools reheader header_ref.txt tmp.bam > ${prefix}.split_ref.sorted.rmdup.bam && samtools index ${prefix}.split_ref.sorted.rmdup.bam && rm tmp.bam
   NB_REF_MAPPED=`samtools view -c ${prefix}.split_ref.sorted.rmdup.bam`
   
   samtools view -o tmp.bam ${bamFiles[0]} ${si_seq_ids.join(' ')}
   grep -vP "${'SN:'+ref_seq_ids.join('\\s|SN:')}" header.txt > header_spike_in.txt
   cat header_spike_in.txt tmp.bam | samtools view -bSh - > ${prefix}.split_spike_in.sorted.rmdup.bam && samtools index ${prefix}.split_spike_in.sorted.rmdup.bam && rm tmp.bam
   NB_SPIKE_IN_MAPPED=`samtools view -c ${prefix}.split_spike_in.sorted.rmdup.bam`

   NORM_FACTOR=`R --slave -q -e "cat(round((10e6/$NB_TOTAL_MAPPED)*(${params.spike_in_fraction}/($NB_SPIKE_IN_MAPPED/$NB_TOTAL_MAPPED)),6), '\n')"`
   echo \$NORM_FACTOR
   """
}
}
else{
   //Set the channels as the ouput of samtools process, adding 1 as scale factor
   samtooled_ch.map{ it -> [it[0],it[1],it[2],0 ] }.set{to_bamCov_ch}
   mapped_reads_ch.set{to_count_mapped_reads_ch}
   samtooled_rmdup_ch.map{ it -> [it[0],it[1],it[2],0 ] }.set{to_bamCov_rmdup_ch}
   mapped_uniq_reads_ch.set{to_count_uniq_mapped_reads_ch}
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
  	tuple val(LibName), val(prefix), path(bamFiles), val(scaleF) from to_bamCov_ch
	output:
	tuple val(LibName), val(prefix), bamFiles, val("${prefix}.bin${params.bin_size}.RPM.bamCoverage.bw") into genCoved_ch
   file("${prefix}.bin${params.bin_size}.RPM.bamCoverage.bw")
  	 
	"""
	if(! ${params.spike_in_norm} ); then
		bamCoverage \
	  	 -b ${bamFiles[0]} \
	  	 -o ${prefix}.bin${params.bin_size}.RPM.bamCoverage.bw -of bigwig \
	  	 ${params.bamcoverage_options} --binSize ${params.bin_size} -p ${task.cpus}
	else
		bamCoverage \
		-b ${bamFiles[0]} \
		-o ${prefix}.bin${params.bin_size}.RPM.bamCoverage.bw -of bigwig \
		${params.bamcoverage_options} --scaleFactor ${scaleF} --binSize ${params.bin_size} -p ${task.cpus}
	fi	
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
  	tuple val(LibName), val(prefix), path(bamFiles), val(scaleF) from to_bamCov_rmdup_ch
   output:
	tuple val(LibName), val(prefix), bamFiles, val("${prefix}.bin${params.bin_size}.RPM.rmdup.bamCoverage.bw") into genCoved_uniq_ch
   file("${prefix}.bin${params.bin_size}.RPM.rmdup.bamCoverage.bw") 

	"""
	if(! ${params.spike_in_norm} ); then
	   bamCoverage \
	   -b ${bamFiles[0]} \
	   -o ${prefix}.bin${params.bin_size}.RPM.rmdup.bamCoverage.bw -of bigwig \
	   ${params.bamcoverage_options} --binSize ${params.bin_size} -p ${task.cpus}
	else
	    bamCoverage \
	   -b ${bamFiles[0]} \
	   -o ${prefix}.bin${params.bin_size}.RPM.rmdup.bamCoverage.bw -of bigwig \
	   ${params.bamcoverage_options} --scaleFactor ${scaleF}  --binSize ${params.bin_size} -p ${task.cpus}
	fi
   """
}

/*
 * Step 5. Get the number of mapped reads.
 *    Parallelized with two different processes
 *       - for .bam           -for .rmdup.bam
 */
ch_Toreport_mapped_nb
   .join(to_count_mapped_reads_ch)
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

ch_Toreport_all_stats
.map{ it -> [it[1],it[0],it[2],it[3],it[4],it[5]  ] }
.map{it -> [it.join("\t")]}.collect().set{ ch_report_all_stats} //Joining stats with "\t" then use collect to have a single entry channel

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
   echo "${x.join('\n')}" | sort -k1 | awk '{for(i=2;i<=NF;i++) printf \$i";"; print ""}' >> mapping_stats.txt
   """
}



ch_Toreport_uniq_nb
   .join(to_count_uniq_mapped_reads_ch)
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

ch_Toreport_uniq_stats
.map{ it -> [it[1],it[0],it[2],it[3],it[4],it[5]  ] }
.map{it -> [it.join("\t")]}.collect().set{ ch_report_uniq_stats} //Joining stats with ";" then use collect to have a single entry channel

process _report_mapping_uniq_stats_csv {
   publishDir "${params.outdir}/Stats", mode: 'copy'
   input:
   val x from ch_report_uniq_stats
   output:
   path("mapping_uniq_stats.txt")
   // echoing all the channel with join('\n') into the mapping_stats.txt file
   script:
   """
   echo "LibName;Nb_sequenced_read;Nb_trimmed_reads;Nb_mapped_reads;Median_insert_size" > mapping_uniq_stats.txt
   echo "${x.join('\n')}" | sort -k1 | awk '{for(i=2;i<=NF;i++) printf \$i";"; print ""}' >> mapping_uniq_stats.txt
   """
}

/*
   PRODUCING csv files that serves as input design for nf-AnalysesOnCoordinates.nf
      This should be updated in aggreement with the subsequent pipeline.
*/

genCoved_ch.join(ch_ToAoC)
.map{ it -> [it[4], it[0], it[2][0], it[3], it[5], it[7], 'NA', it[8],  1, '', '', '', '', '', '', '']}
.map{ it -> [it.join("\t")]}
.collect()
.set {ch_report_Aoc}
 
genCoved_uniq_ch.join(ch_ToAoC_uniq)
.map{ it -> [it[4], it[0], it[2][0], it[3], it[5], it[7], it[7], it[8], 1, '', '', '', '', '', '', '']}
.map{ it -> [it.join("\t")]}
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
   echo "${x.join('\n')}" | sort -k1 | awk '{for(i=2;i<=NF;i++) printf \$i";"; print ""}' >> ${params.name}.bigwigDesign.csv
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
  echo "${x.join('\n')}" | sort -k1 | awk '{for(i=2;i<=NF;i++) printf \$i";"; print ""}' >> ${params.name}.rmdup.bigwigDesign.csv
   """
}
