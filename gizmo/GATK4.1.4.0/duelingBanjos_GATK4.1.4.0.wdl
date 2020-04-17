version 1.0
## Consensus variant calling workflow for human hybrid capture-based paired-sample panel based DNA sequencing.
## Input requirements:
## - For both "sample" and "reference", pair-end sequencing data in unmapped BAM (uBAM) format that comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile (a Picard tool)
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
##
## Output Files:
## - An analysis-ready recalibrated bam and it's index for both the "sample" and the "reference"
## - GATK-Mutect2 vcf from GATK 4.1.4.0
## - Strelka vcf
## - Annovar annotated vcfs and tabular variant list for each variant caller
## - Basic QC stats from bedtools for mean coverage over regions in panel
## - QC stats from Picard HSMetrics tool
## 
workflow HybridCap_BWA_Mutect2_Strelka_AnnotatedVariants {
  input {
    # Batch information
    File batchFile
    # Reference Genome information
    String ref_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative". Leave blank in JSON for legacy
    # references such as b37 and hg19.
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    # Annovar package
    String annovarDIR
    String annovar_protocols
    String annovar_operation
    # Gizmo Modules
    String GATKModule
    String samtoolsModule
    String bedtoolsModule
    String perlModule
    String bwaModule
    String strelkaModule
    String bcftoolsModule
    String RModule
    # consensus script
    String githubRepoURL
    String githubTag
  }

Array[Object] batchInfo = read_objects(batchFile)

scatter (job in batchInfo){
  # variables
  String sampleName = job.sampleName
  String referenceName = job.referenceName
  String sampleID = job.molecular_id
  String referenceID = job.ref_molecular_id
  # s3 files to retrieve
  File refBam = job.refBamLocation
  File sampleBam = job.sampleBamLocation
  File bedFile = job.bedLocation


# Get the basename, i.e. strip the filepath and the extension
  String bam_basename = basename(sampleBam, ".unmapped.bam")
  String ref_basename = basename(refBam, ".unmapped.bam")
  String base_file_name = bam_basename + "." + ref_name
  String ref_file_name = ref_basename + "." + ref_name


  # Prepare bed file and check sorting
  call SortBed {
    input:
      unsorted_bed = bedFile,
      ref_dict = ref_dict,
      modules = GATKModule
  }

  # Convert unmapped bam to interleaved fastq
  call SamToFastq as sampleSamToFastq {
    input:
      input_bam = sampleBam,
      base_file_name = base_file_name,
      modules = GATKModule
  }

  # Convert unmapped bam to interleaved fastq
  call SamToFastq as refSamToFastq {
    input:
      input_bam = refBam,
      base_file_name = ref_file_name,
      modules = GATKModule
  }

  #  Map reads to reference
  call BwaMem as sampleBwaMem {
    input:
      input_fastq = sampleSamToFastq.output_fastq,
      base_file_name = base_file_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_alt = ref_alt,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      modules = samtoolsModule + " " + bwaModule
  }

  #  Map reads to reference
  call BwaMem as refBwaMem {
    input:
      input_fastq = refSamToFastq.output_fastq,
      base_file_name = ref_file_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_alt = ref_alt,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      modules = samtoolsModule + " " + bwaModule
  }

  # Merge original uBAM and BWA-aligned BAM
  call MergeBamAlignment as sampleMergeBamAlignment {
    input:
      unmapped_bam = sampleBam,
      aligned_bam = sampleBwaMem.output_bam,
      base_file_name = base_file_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      modules = GATKModule
  }

  # Merge original uBAM and BWA-aligned BAM
  call MergeBamAlignment as refMergeBamAlignment {
    input:
      unmapped_bam = refBam,
      aligned_bam = refBwaMem.output_bam,
      base_file_name = ref_file_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      modules = GATKModule
  }

  # Generate the recalibration model by interval and apply it
  call ApplyBaseRecalibrator  as sampleApplyBaseRecalibrator {
    input:
      input_bam = sampleMergeBamAlignment.output_bam,
      input_bam_index = sampleMergeBamAlignment.output_bai,
      base_file_name = base_file_name,
      intervals = SortBed.intervals,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      modules = samtoolsModule + " " + GATKModule
    }

  # Generate the recalibration model by interval and apply it
  call ApplyBaseRecalibrator  as refApplyBaseRecalibrator {
    input:
      input_bam = refMergeBamAlignment.output_bam,
      input_bam_index = refMergeBamAlignment.output_bai,
      base_file_name = ref_file_name,
      intervals = SortBed.intervals,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      modules = samtoolsModule + " " + GATKModule
    }

    call bedToolsQC as samplebedToolsQC {
      input: 
        input_bam = sampleApplyBaseRecalibrator.recalibrated_bam,
        genome_sort_order = sampleApplyBaseRecalibrator.sortOrder,
        bed_file = SortBed.sorted_bed,
        base_file_name = base_file_name,
        modules = bedtoolsModule
    }
    call bedToolsQC as refbedToolsQC {
      input: 
        input_bam = refApplyBaseRecalibrator.recalibrated_bam,
        genome_sort_order = refApplyBaseRecalibrator.sortOrder,
        bed_file = SortBed.sorted_bed,
        base_file_name = ref_file_name,
        modules = bedtoolsModule
    }

    call CollectHsMetrics as sampleHSMetrics {
      input: 
        input_bam = sampleApplyBaseRecalibrator.recalibrated_bam,
        base_file_name = base_file_name,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        intervals = SortBed.intervals,
        modules = GATKModule
    }

    call CollectHsMetrics as refHSMetrics {
      input: 
        input_bam = refApplyBaseRecalibrator.recalibrated_bam,
        base_file_name = ref_file_name,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        intervals = SortBed.intervals,
        modules = GATKModule
    }

  # Generate Mutect2 vcf
  call Mutect2 {
    input:
      input_bam = sampleApplyBaseRecalibrator.recalibrated_bam,
      base_file_name = base_file_name,
      sampleID = sampleID,
      intervals = SortBed.intervals,
      input_ref_bam = refApplyBaseRecalibrator.recalibrated_bam,
      ref_file_name = ref_file_name,
      referenceID = referenceID,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      modules = samtoolsModule + " " + GATKModule
    }
  # Use Strelka somatic variant caller
  call StrelkaSomatic {
    input:
        tumorBam = sampleApplyBaseRecalibrator.recalibrated_bam,
        tumorBamIndex = sampleApplyBaseRecalibrator.recalibrated_bai,
        base_file_name = base_file_name,
        normalBam = refApplyBaseRecalibrator.recalibrated_bam,
        normalBamIndex = refApplyBaseRecalibrator.recalibrated_bai,
        ref_file_name = ref_file_name,
        referenceFasta = ref_fasta,
        referenceFastaFai = ref_fasta_index,
        bed_file = SortBed.sorted_bed,
        modules = bcftoolsModule + " " + samtoolsModule + " " + strelkaModule
    }

    # Annotate variants
    call annovar as mutectAnnovar {
      input:
        input_vcf = Mutect2.output_vcf,
        ref_name = ref_name,
        annovarDIR = annovarDIR,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols,
        annovarDIR = annovarDIR,
        modules = perlModule
    }
    # Annotate Strelka variants
    call annovar as strelkaAnnovar {
      input:
        input_vcf = StrelkaSomatic.vcf,
        ref_name = ref_name,
        annovarDIR = annovarDIR,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols,
        annovarDIR = annovarDIR,
        modules = perlModule
    }

    call consensusFiltering {
      input:
        MutectVars = mutectAnnovar.output_annotated_table,
        StrVars = strelkaAnnovar.output_annotated_table,
        base_file_name = base_file_name + "." + ref_file_name,
        githubRepoURL = githubRepoURL,
        githubTag = githubTag,
        molecular_id = sampleID,
        ref_molecular_id = referenceID,
        modules = RModule
    }

  # End scatter 
  }
  # Outputs that will be retained when execution is complete
  output {
    Array[File] analysisReadyBam = sampleApplyBaseRecalibrator.recalibrated_bam 
    Array[File] analysisReadyIndex = sampleApplyBaseRecalibrator.recalibrated_bai
    Array[File] analysisReadyRefBam = refApplyBaseRecalibrator.recalibrated_bam
    Array[File] analysisReadyRefIndex = refApplyBaseRecalibrator.recalibrated_bai
    Array[File] mutectVcf = Mutect2.output_vcf
    Array[File] mutectVcfIndex = Mutect2.output_vcf_index
    Array[File] mutectAnnotatedVcf = mutectAnnovar.output_annotated_vcf
    Array[File] mutectAnnotatedTable = mutectAnnovar.output_annotated_table
    Array[File] strelkaVcf = StrelkaSomatic.vcf
    Array[File] strelkaVcfIndex = StrelkaSomatic.vcfIndex
    Array[File] strelkaAnnotatedVcf = strelkaAnnovar.output_annotated_vcf
    Array[File] strelkaAnnotatedTable = strelkaAnnovar.output_annotated_table
    Array[File] sampleBedQC = samplebedToolsQC.meanQC
    Array[File] refBedQC = refbedToolsQC.meanQC
    Array[File] samplePicardQC = sampleHSMetrics.picardMetrics
    Array[File] samplePicardQCpertarget = sampleHSMetrics.picardPerTarget
    Array[File] refPicardQC = refHSMetrics.picardMetrics
    Array[File] refPicardQCpertarget = refHSMetrics.picardPerTarget
    Array[File] consensus = consensusFiltering.consensusfile
  }
# End workflow
}
#### TASK DEFINITIONS

# Prepare bed file and check sorting
task SortBed {
  input {
    File unsorted_bed
    File ref_dict
    String modules
  }
  command {
    set -eo pipefail

    echo "Sort bed file"
    sort -k1,1V -k2,2n -k3,3n ${unsorted_bed} > sorted.bed

    echo "Transform bed file to intervals list with Picard----------------------------------------"
    gatk --java-options "-Xms4g" \
      BedToIntervalList \
      -I=sorted.bed \
      -O=sorted.interval_list \
      -SD=${ref_dict}
  }
  output {
    File intervals = "sorted.interval_list"
    File sorted_bed = "sorted.bed"
  }
  runtime {
    memory: "8GB"
    modules: "${modules}"
  }
}

# Read unmapped BAM, convert to FASTQ
task SamToFastq {
  input {
    File input_bam
    String base_file_name
    String modules
  }
  command {
    set -eo pipefail

    gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g" \
      SamToFastq \
      --INPUT=${input_bam} \
      --FASTQ=${base_file_name}.fastq \
      --INTERLEAVE=true \
      --INCLUDE_NON_PF_READS=true 
  }
  output {
    File output_fastq = "${base_file_name}.fastq"
  }
  runtime {
    memory: "16GB"
    modules: "${modules}"
  }
}

# align to genome
## Currently uses -M but GATK uses -Y and no -M
task BwaMem {
  input {
    File input_fastq
    String base_file_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    String modules
  }
  command {
    set -eo pipefail

    bwa mem \
      -p -v 3 -t 16 -M \
      ${ref_fasta} ${input_fastq} > ${base_file_name}.sam 
    samtools view -1bS -@ 15 -o ${base_file_name}.aligned.bam ${base_file_name}.sam
  }
  output {
    File output_bam = "${base_file_name}.aligned.bam"
  }
  runtime {
    memory: "48GB"
    cpu: 16
    partition: "largenode"
    modules: "${modules}"
  }
}


# Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
  input {
    File unmapped_bam
    File aligned_bam
    String base_file_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String modules
  }
  command {
    set -eo pipefail

    gatk --java-options "-Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Xmx8g" \
      MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ALIGNED_BAM ${aligned_bam} \
      --UNMAPPED_BAM ${unmapped_bam} \
      --OUTPUT ${base_file_name}.merged.bam \
      --REFERENCE_SEQUENCE ${ref_fasta} \
      --PAIRED_RUN true \
      --SORT_ORDER coordinate \
      --IS_BISULFITE_SEQUENCE false \
      --ALIGNED_READS_ONLY false \
      --CLIP_ADAPTERS false \
      --MAX_RECORDS_IN_RAM 200000 \
      --ADD_MATE_CIGAR true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
      --ALIGNER_PROPER_PAIR_FLAGS true \
      --CREATE_INDEX true
  }
  output {
    File output_bam = "${base_file_name}.merged.bam"
    File output_bai = "${base_file_name}.merged.bai"
  }
  runtime {
    memory: "16GB"
    modules: "${modules}"
  }
}

# Generate Base Quality Score Recalibration (BQSR) model and apply it
task ApplyBaseRecalibrator {
  input {
    File input_bam
    File intervals 
    File input_bam_index
    String base_file_name
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String modules
  }
  command {
  set -eo pipefail

  samtools index ${input_bam}

  gatk --java-options "-Xms4g" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${base_file_name}.recal_data.csv \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
      --intervals ${intervals} \
      --interval-padding 100 

  gatk --java-options "-Xms4g" \
      ApplyBQSR \
      -bqsr ${base_file_name}.recal_data.csv \
      -I ${input_bam} \
      -O ${base_file_name}.recal.bam \
      -R ${ref_fasta} \
      --intervals ${intervals} \
      --interval-padding 100 

  #finds the current sort order of this bam file
  samtools view -H ${base_file_name}.recal.bam | grep @SQ | sed 's/@SQ\tSN:\|LN://g' > ${base_file_name}.sortOrder.txt

  }
  output {
    File recalibrated_bam = "${base_file_name}.recal.bam"
    File recalibrated_bai = "${base_file_name}.recal.bai"
    File sortOrder = "${base_file_name}.sortOrder.txt"
  }
  runtime {
    memory: "33GB"
    cpu: 6
    partition: "largenode"
    modules: "${modules}"
  }
}

# use bedtools to find basic QC data
task bedToolsQC {
  input {
    File input_bam
    File bed_file
    File genome_sort_order
    String base_file_name
    String modules
  }
  command {
  set -eo pipefail

  bedtools sort -g ${genome_sort_order} -i ${bed_file} > correctly.sorted.bed
  bedtools coverage -mean -sorted -g ${genome_sort_order} -a correctly.sorted.bed \
      -b ${input_bam} > ${base_file_name}.bedtoolsQCMean.txt
  }
  output {
    File meanQC = "${base_file_name}.bedtoolsQCMean.txt"
  }
  runtime {
    modules: "${modules}"
  }
}


# Mutect 2 calling
task Mutect2 {
  input {
    File input_bam
    String base_file_name
    File input_ref_bam
    String ref_file_name
    File intervals
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String sampleID
    String referenceID
    String modules
  }
  command {
    set -eo pipefail

  # Just in case the bam header isn't the referenceID or sampleID
    gatk --java-options "-Xmx24g" \
      AddOrReplaceReadGroups \
        -I=${input_ref_bam} \
        -O=reference.bam \
        -SM=${referenceID} \
        -LB=unknown \
        -PL=illumina \
        -PU=unknown 


    gatk --java-options "-Xmx24g" \
      AddOrReplaceReadGroups \
        -I=${input_bam} \
        -O=tumor.bam \
        -SM=${sampleID} \
        -LB=unknown \
        -PL=illumina \
        -PU=unknown 
    

    samtools index tumor.bam
    samtools index reference.bam
  

  ## Mutect with FFPE artifact bias filter
    gatk --java-options "-Xmx24g" \
      Mutect2 \
        -R ${ref_fasta} \
        -I tumor.bam \
        -I reference.bam \
        -normal ${referenceID} \
        -tumor ${sampleID} \
        -O preliminary.vcf.gz \
        --intervals ${intervals} \
        --interval-padding 100  \
        --f1r2-tar-gz f1r2.tar.gz 

  # This creates an output with raw data used to learn the orientation bias model
    gatk --java-options "-Xmx24g" \
      LearnReadOrientationModel \
        -I f1r2.tar.gz \
        -O read-orientation-model.tar.gz

  # Finally, pass the learned read orientation model to FilterMutectCallswith the -ob-priors argument:
    gatk --java-options "-Xmx24g" \
      FilterMutectCalls \
        -V preliminary.vcf.gz \
        --ob-priors read-orientation-model.tar.gz \
        -O ${base_file_name}_${ref_file_name}.mutect2.vcf.gz \
        -R ${ref_fasta} 
    }
  output {
    File output_vcf = "${base_file_name}_${ref_file_name}.mutect2.vcf.gz"
    File output_vcf_index = "${base_file_name}_${ref_file_name}.mutect2.vcf.gz.tbi"
  }
  runtime {
    memory: "30GB"
    cpu: 4
    modules: "${modules}"
  }
}

# annotate with annovar
task annovar {
  input {
    File input_vcf
    String ref_name
    String annovar_protocols
    String annovar_operation
    String annovarDIR
    String modules
  }
  String base_vcf_name = basename(input_vcf, ".vcf")
  command {
  set -eo pipefail

  perl ${annovarDIR}/annovar/table_annovar.pl ${input_vcf} ${annovarDIR}/annovar/humandb/ \
    -buildver ${ref_name} \
    -outfile ${base_vcf_name} \
    -remove \
    -protocol ${annovar_protocols} \
    -operation ${annovar_operation} \
    -nastring . -vcfinput
  }
  output {
    File output_annotated_vcf = "${base_vcf_name}.${ref_name}_multianno.vcf"
    File output_annotated_table = "${base_vcf_name}.${ref_name}_multianno.txt"
  }
  runtime {
    modules: "${modules}"
  }
}


task StrelkaSomatic {
  input {
    File normalBam
    File normalBamIndex
    File tumorBam
    File tumorBamIndex
    File referenceFasta
    File referenceFastaFai
    File bed_file
    String base_file_name
    String ref_file_name
    String modules
  }
  command {
    set -eo pipefail
    mkdir tempDir

    # zip and index the bed file for Strelka
    bgzip ${bed_file}
    tabix -0 -p bed -s 1 ${bed_file}.gz
      

    # run strelka
    configureStrelkaSomaticWorkflow.py \
        --normalBam ${normalBam} \
        --tumorBam ${tumorBam} \
        --ref ${referenceFasta} \
        --runDir tempDir \
        --callRegions ${bed_file}.gz \
        --exome

    tempDir/runWorkflow.py \
        -m local \
        -j 16 \
        -g 24
        ## j is the number of cores available
        ## g is the amount of memory available


    # Concatenate the SNVs and indels
    bcftools concat -a \
      --output-type v \
      --output ${base_file_name}_${ref_file_name}.strelka.noGT.vcf \
      tempDir/results/variants/somatic.indels.vcf.gz tempDir/results/variants/somatic.snvs.vcf.gz

    # Add in a dummy GT tag b/c Strelka doesn't do that and Annovar needs it.
    first_format_num=$(grep -n -m 1 '##FORMAT' "${base_file_name}_${ref_file_name}.strelka.noGT.vcf" | cut -d : -f 1)
    sed "$first_format_num"'i##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' "${base_file_name}_${ref_file_name}.strelka.noGT.vcf" > "${base_file_name}_${ref_file_name}.strelka.vcf"
    sed -ri 's|(DP:)|GT:\1|g' "${base_file_name}_${ref_file_name}.strelka.vcf"
    sed -ri 's|(:BCN50\t)|\10/0:|g' "${base_file_name}_${ref_file_name}.strelka.vcf"
    sed -ri 's|(:BCN50\t[^\t]*\t)|\10/1:|g' "${base_file_name}_${ref_file_name}.strelka.vcf"

    # Zip it and index it
    bgzip ${base_file_name}_${ref_file_name}.strelka.vcf
    tabix -0 -p vcf -s 1 ${base_file_name}_${ref_file_name}.strelka.vcf.gz

    }
    output {
      File vcf = "${base_file_name}_${ref_file_name}.strelka.vcf.gz"
      File vcfIndex = "${base_file_name}_${ref_file_name}.strelka.vcf.gz.tbi"
    }
    runtime {
      memory: "48GB"
      cpu: 16
      partition: "largenode"
      modules: "${modules}"
    }
}

# get hybrid capture based QC metrics via Picard
task CollectHsMetrics {
  input {
    File input_bam
    String base_file_name
    File ref_fasta
    File ref_fasta_index
    File intervals
    String modules
  }
  command {
    set -eo pipefail

    gatk --java-options "-Xmx4g" \
      CollectHsMetrics \
      --INPUT=${input_bam} \
      --OUTPUT=${base_file_name}.picard.metrics.txt \
      --REFERENCE_SEQUENCE=${ref_fasta} \
      --ALLELE_FRACTION=0.01 \
      --BAIT_INTERVALS=${intervals} \
      --TARGET_INTERVALS=${intervals} \
      --PER_TARGET_COVERAGE=${base_file_name}.picard.pertarget.txt 
  }
  output {
    File picardMetrics = "${base_file_name}.picard.metrics.txt"
    File picardPerTarget = "${base_file_name}.picard.pertarget.txt"
  }
  runtime {
    memory: "4GB"
    modules: "${modules}"
  }
}

task consensusFiltering {
  input {
    File MutectVars
    File StrVars
    String base_file_name
    String githubRepoURL
    String githubTag
    String molecular_id
    String ref_molecular_id
    String modules
  }
  command {
    set -eo pipefail
    git clone --branch ${githubTag} ${githubRepoURL}
    Rscript ./tg-wdl-LILAC-workflow/gizmo/ConsensusFiltering/consensus-Mutect2-Strelka.R \
      ${MutectVars} ${StrVars} ${base_file_name} ${molecular_id} ${ref_molecular_id}
  }
  output {
    File consensusfile = "${base_file_name}.consensus.tsv"
  }
  runtime {
    memory: "8GB"
    modules: "${modules}"
  }
}
