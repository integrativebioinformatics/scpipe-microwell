process DGE_MATRIX {
    label 'dropseq_tools'
    publishDir "results/${sample_id}/", mode: 'copy'

    input:
        tuple path(read1), path(read2), path(gtf), path(genome_fasta), path(genome_dict), val(sample_id)
    output:
        path "${sample_id}.dge.txt.gz"
        path "${sample_id}.dge.summary.txt"
        path "metadata/outputs/${sample_id}_tagged_cell_summary.txt"
        path "metadata/outputs/${sample_id}_tagged_mol_summary.txt"
        path "metadata/outputs/${sample_id}_substitution_report.txt"
        path "metadata/outputs/${sample_id}_synthesis_stats.summary.txt"

    script:
    """
    mkdir -p metadata/intervals_dir

    # ---- 0. FASTQ TO BAM ----
    picard FastqToSam F1=${read1} F2=${read2} O=${sample_id}.bam SM=${sample_id}

    # ---- 1. CreateSequenceDictionary ----
    if [ ! -f metadata/genome.dict ]; then
        echo "[Dropseq] Creando metadata/genome.dict..."
        picard CreateSequenceDictionary REFERENCE=${genome_fasta} OUTPUT=metadata/genome.dict
    fi

    # ---- 2. ConvertToRefFlat ----
    if [ ! -f metadata/genes.refFlat ]; then
        echo "[Dropseq] Creando metadata/genes.refFlat..."
        java -jar /opt/dropseq-3.0.2/lib/dropseq.jar ConvertToRefFlat \
            ANNOTATIONS_FILE=${gtf} \
            SEQUENCE_DICTIONARY=metadata/genome.dict \
            OUTPUT=metadata/genes.refFlat
    fi

    # ---- 3. ReduceGTF ----
    if [ ! -f metadata/reduced.gtf ]; then
        echo "[Dropseq] Creando metadata/reduced.gtf..."
        java -jar /opt/dropseq-3.0.2/lib/dropseq.jar ReduceGtf \
            SEQUENCE_DICTIONARY=metadata/genome.dict \
            GTF=${gtf} \
            OUTPUT=metadata/reduced.gtf
    fi

    # ---- 4. CreateIntervalsFiles ----
    if [ ! -d metadata/intervals_dir ] || [ ! -f metadata/intervals_dir/genome.genes.intervals ]; then
        echo "[Dropseq] Creando archivos de intervalos..."
        java -jar /opt/dropseq-3.0.2/lib/dropseq.jar CreateIntervalsFiles \
            SEQUENCE_DICTIONARY=metadata/genome.dict \
            REDUCED_GTF=metadata/reduced.gtf \
            PREFIX=genome \
            OUTPUT=metadata/intervals_dir \
            MT_SEQUENCE=${params.mt_sequence}
    fi

    # ---- 5. STAR genomeGenerate ----
    echo "Usando CPUs: ${task.cpus}"
    if [ ! -d star_index ] || [ ! -f star_index/Genome ]; then
        echo "[Dropseq] Creando índice STAR..."
        mkdir -p star_index
        STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles ${genome_fasta} --sjdbGTFfile ${gtf}
    fi

    mkdir -p metadata/outputs

    # ---- 6. DROP-SEQ WORKFLOW ----

    java -jar /opt/dropseq-3.0.2/lib/dropseq.jar TagBamWithReadSequenceExtended INPUT=${sample_id}.bam OUTPUT=tagged_cell.bam \
        SUMMARY=metadata/outputs/${sample_id}_tagged_cell_summary.txt \
        BASE_RANGE=${params.tag1_base_range} \
        BARCODED_READ=${params.tag1_barcode_read} \
        TAG_NAME=${params.tag1_tag_name} \
        BASE_QUALITY=${params.tag1_base_quality} \
        NUM_BASES_BELOW_QUALITY=${params.tag1_num_bases_below_quality}

    java -jar /opt/dropseq-3.0.2/lib/dropseq.jar TagBamWithReadSequenceExtended INPUT=tagged_cell.bam OUTPUT=tagged_mol.bam \
        SUMMARY=metadata/outputs/${sample_id}_tagged_mol_summary.txt \
        BASE_RANGE=${params.tag2_base_range} \
        BARCODED_READ=${params.tag2_barcode_read} \
        TAG_NAME=${params.tag2_tag_name} \
        BASE_QUALITY=${params.tag2_base_quality} \
        NUM_BASES_BELOW_QUALITY=${params.tag2_num_bases_below_quality} \
        DISCARD_READ=True

    java -jar /opt/dropseq-3.0.2/lib/dropseq.jar FilterBam INPUT=tagged_mol.bam OUTPUT=filtered.bam TAG_REJECT=XQ

    java -jar /opt/dropseq-3.0.2/lib/dropseq.jar TrimStartingSequence INPUT=filtered.bam OUTPUT=trimmed_smart.bam \
        SEQUENCE=${params.smart_adapter_sequence} \
        NUM_BASES=${params.smart_adapter_num_bases} \
        MISMATCHES=${params.smart_adapter_mismatches}

    java -jar /opt/dropseq-3.0.2/lib/dropseq.jar PolyATrimmer INPUT=trimmed_smart.bam OUTPUT=polyA_trimmed.bam \
        NUM_BASES=${params.polyA_num_bases} \
        MISMATCHES=${params.polyA_mismatches} \
        USE_NEW_TRIMMER=${params.polyA_use_new_trimmer}

    picard SamToFastq INPUT=polyA_trimmed.bam FASTQ=reads.fastq

    STAR --runThreadN ${task.cpus} --genomeDir star_index --readFilesIn reads.fastq --outFileNamePrefix star_

    picard SortSam I=star_Aligned.out.sam O=aligned.sorted.bam SO=queryname

    picard MergeBamAlignment REFERENCE_SEQUENCE=${genome_fasta} UNMAPPED_BAM=polyA_trimmed.bam \
        ALIGNED_BAM=aligned.sorted.bam OUTPUT=merged.bam INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false

    java -jar /opt/dropseq-3.0.2/lib/dropseq.jar TagReadWithGeneFunction I=merged.bam O=gene_tagged.bam ANNOTATIONS_FILE=metadata/genes.refFlat

    java -jar /opt/dropseq-3.0.2/lib/dropseq.jar DetectBeadSubstitutionErrors I=gene_tagged.bam O=subst_corrected.bam \
        OUTPUT_REPORT=metadata/outputs/${sample_id}_substitution_report.txt

    java -jar /opt/dropseq-3.0.2/lib/dropseq.jar DetectBeadSynthesisErrors I=subst_corrected.bam O=synth_corrected.bam \
        REPORT=indel_report.txt \
        OUTPUT_STATS=synthesis_stats.txt \
        SUMMARY=metadata/outputs/${sample_id}_synthesis_stats.summary.txt

    java -jar /opt/dropseq-3.0.2/lib/dropseq.jar DigitalExpression I=synth_corrected.bam \
        O=${sample_id}.dge.txt.gz \
        SUMMARY=${sample_id}.dge.summary.txt \
        NUM_CORE_BARCODES=${params.num_core_barcodes} \
        READ_MQ=${params.read_mq} \
        EDIT_DISTANCE=${params.edit_distance} \
        LOCUS_FUNCTION_LIST=${params.locus_function_list} \
        STRAND_STRATEGY=${params.strand_strategy}
    """
}

