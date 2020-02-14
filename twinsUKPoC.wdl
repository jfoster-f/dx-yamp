version 1.0

import "yamp.wdl" as yamp

workflow TwinsUKPoC{
    input{
        File reads1_file
        File? reads2_file
        File adapters_file
        File artifacts_file
        File phix174ill_file
        File ref_foreign_genome_qz
        File mpa_pkl_file
        File bowtie2db_gz_file
        File chocophlan_gz_file
        File uniref_file
        String bowtiedb_files_name
        String bt2_options = "very-sensitive"
        String library_layout = "paired"
        String prefix
        Int qin = 33
        Int kcontaminants = 23
        Int mink = 11
        Int hdist = 1
        Int phred = 10
        Int min_length = 60
        Float mind = 0.95
        Int max_indel = 3
        Float bwr = 0.16
        String docker_image = "alesssia/yampdocker"
    }
    #if reads2 is not null i think we need to concatonate reads1 and reads2 maybe handle this as Array[File] and do in command block
    #where are step, stem and label defined in Nextflow?
    call yamp.QualityAssessment as RawQC{
        input:
            reads = reads1_file,
            step = "1",
            stem = "_R1",
            label = "_rawreads",
            prefix = prefix,
            docker_image = docker_image
    }

    #if read2_file is not null then need to
    call yamp.QualityAssessment as RawQC{
        input:
            reads = reads2_file,
            step = "1",
            stem = "_R2",
            label = "_rawreads",
            prefix = prefix,
            docker_image = docker_image
    }

    call yamp.Dedup {
        input:
            reads1_file=reads1_file,
            reads2_file=reads2_file,
            library_layout = library_layout,
            qin = qin,
            prefix = prefix,
            docker_image = docker_image,
    }


    call yamp.Trim {
        input:
            reads1_file = Dedup.dedupped_reads1,
            reads2_file = Dedup.dedupped_reads2[0],
            library_layout = library_layout,
            prefix = prefix,
            kcontaminants = kcontaminants,
            mink = mink,
            hdist = hdist,
            phred = phred,
            min_length = min_length,
            qin = qin,
            adapters_file = adapters_file,
            artifacts_file = artifacts_file,
            phix174ill_file = phix174ill_file,
            docker_image = docker_image

    }
    call yamp.QualityAssessment as TrimQC {
        input:
            reads = Trim.trimmed_reads,
            step = "4",
            stem = "",
            label = "_trimmedreads",
            prefix = prefix,
            docker_image = docker_image
    }


    call yamp.Decontaminate {
        input:
            infile1 = Trim.trimmed_reads,
            infile2 = Trim.trimmed_reads_r1,
            infile12 = Trim.trimmed_reads_r1[0],
            library_layout = library_layout,
            prefix = prefix,
            mind = mind,
            max_indel = max_indel,
            phred = phred,
            qin = qin,
            bwr = bwr,
            ref_foreign_genome_gz = ref_foreign_genome_qz,
            docker_image = docker_image
    }

    call yamp.QualityAssessment as DecontaminateQC {
        input:
            reads = Decontaminate.cleaned_reads,
            step = "6",
            stem = "",
            label = "_decontaminatedreads",
            prefix = prefix,
            docker_image = docker_image
    }

    call yamp.ProfileTaxa {
        input:
            in_file = Decontaminate.cleaned_reads,
            mpa_pkl_file = mpa_pkl_file,
            bowtie2db_gz_file = bowtie2db_gz_file,
            bowtiedb_files_name = bowtiedb_files_name,
            bt2_options = bt2_options,
            prefix = prefix,
            docker_image = docker_image
    }

    call yamp.AlphaDiversity {
        input:
            biom_file = ProfileTaxa.biom_file,
            prefix = prefix,
            docker_image = docker_image

    }

    call yamp.ProfileFunction {
        input:
            clean_reads_file = Decontaminate.cleaned_reads,
            metaphlan_file = ProfileTaxa.metaphlan_file,
            chocophlan_gz_file = chocophlan_gz_file,
            uniref_file = uniref_file,
            prefix = prefix,
            docker_image = docker_image
    }

    output {
        File dedupped_reads1 = Dedup.dedupped_reads1
        File dedupped_reads2 = Dedup.dedupped_reads2[0]
        File trimmed_reads = Trim.trimmed_reads
        File raw_qc_html = RawQC.html_file
        File raw_qc_fastqc = RawQC.fastqc_data_file
        File trim_qc_html = TrimQC.html_file
        File trim_qc_fastqc = TrimQC.fastqc_data_file
        File decontaminate_qc_html = DecontaminateQC.html_file
        File decontaminate_qc_fastqc = DecontaminateQC.fastqc_data_file
        File cleaned_reads = Decontaminate.cleaned_reads
        File contaminated_reads = Decontaminate.contaminated_reads
        File biom_file = ProfileTaxa.biom_file
        File metaphlan_file = ProfileTaxa.metaphlan_file
        File bt2out_file = ProfileTaxa.bt2out_file
        File alpha_diversity_file = AlphaDiversity.alpha_diversity_file
        File HUMAnN2_log_file = ProfileFunction.HUMAnN2_log_file
        File gene_families_file = ProfileFunction.gene_families_file
        File path_coverage_file = ProfileFunction.path_coverage_file
        File path_abundance_file = ProfileFunction.path_abundance_file
        File bowtie2_aligned_sam_file = ProfileFunction.bowtie2_aligned_sam_file
        File bowtie2_aligned_tsv_file = ProfileFunction.bowtie2_aligned_tsv_file
        File diamond_aligned_tsv_file = ProfileFunction.diamond_aligned_tsv_file
        File bowtie2_unaligned_file = ProfileFunction.bowtie2_unaligned_file
        File diamond_unaligned_file = ProfileFunction.diamond_unaligned_file

    }


}