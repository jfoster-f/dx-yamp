version 1.0

task Dedup {
    input {
        File reads1_file
        File? reads2_file
        String library_layout = "paired"
        String prefix
        Int qin = 33
        String docker_image = "yamp:latest"

    }

    meta {
        description: "Quality Control - STEP 1. De-duplication. Only exact duplicates are removed. This step is OPTIONAL. De-duplication should be carried on if you are using PCR amplification (in this case identical reads are technical artefacts) but not otherwise (identical reads will identify natural duplicates)."
    }

    parameter_meta {
        reads1_file: {
            description: "Input reads in fastq or fastq.gz format",
            extension: ['fastq', 'fastq.gz']
        }
        reads2_file: {
            description: "Paired read mates in fastq or fastq.gz format",
            extension: ['fastq', 'fastq.gz', 'fq']
        }
        library_layout: "One of either 'paired' or 'single' describing the input reads"
        prefix: "Filename prefix"
        qin: "Input quality offset: 33 (ASCII+33) or 64 (ASCII+64)"
        docker_image:"Name of docker image to be used"
    }

    #String cpus = "4"
    String cpus = "1"
    #String memory_amount = "32 GB"
    String memory_amount = "4 GB"
    command <<<
        sysdate=$(date)
        starttime=$(date +%s.%N)
        echo "Performing Quality Control. STEP 1 [De-duplication] at $sysdate"
        echo " "
        #Sets the maximum memory to the value requested in the config file
        maxmem=$(echo "~{memory_amount}" | sed 's/ //g' | sed 's/B//g')

        #Defines command for de-duplication
        if [ ~{library_layout} = "paired" ]; then
            #CMD="clumpify.sh -Xmx"$maxmem" in1=~{reads1_file} in2=~{reads2_file} out1=~{prefix}_dedupe_R1.fq out2=~{prefix}_dedupe_R2.fq qin=~{qin} dedupe subs=0 threads=~{cpus}"
            #CMD="clumpify.sh in1=~{reads1_file} in2=~{reads2_file} out1=~{prefix}_dedupe_R1.fq out2=~{prefix}_dedupe_R2.fq qin=~{qin} dedupe subs=0 threads=~{cpus}"
            CMD="clumpify.sh in1=~{reads1_file} in2=~{reads2_file} out1=~{prefix}_dedupe_R1.fq out2=~{prefix}_dedupe_R2.fq qin=~{qin} dedupe subs=0"
        else
            #CMD="clumpify.sh -Xmx"$maxmem" in=~{reads1_file} out=~{prefix}_dedupe.fq qin=~{qin} dedupe subs=0 threads=~{cpus}"
            #CMD="clumpify.sh in=~{reads1_file} out=~{prefix}_dedupe.fq qin=~{qin} dedupe subs=0 threads=~{cpus}"
            CMD="clumpify.sh in=~{reads1_file} out=~{prefix}_dedupe.fq qin=~{qin} dedupe subs=0"
        fi
        #Logs version of the software and executed command (BBmap prints on stderr)
        version=$(clumpify.sh --version 2>&1 >/dev/null | grep "BBMap version")
        echo "Using clumpify.sh in $version "
        echo "Executing command: $CMD "
        echo " "

        #De-duplicates
        exec $CMD 2>&1 | tee tmp.log
        #Logs some figures about sequences passing de-duplication
        echo "Clumpify's de-duplication stats: "
        echo " "
        echo sed -n '/Reads In:/,/Duplicates Found:/p' tmp.log
        echo " "
        totR=$(grep "Reads In:" tmp.log | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
        echo $totR
        remR=$(grep "Duplicates Found:" tmp.log | cut -f 1 | cut -d: -f 2 | sed 's/ //g')
        echo $remR
        survivedR=$((totR-remR))
        percentage=$(echo $survivedR $totR | awk '{print $1/$2*100}' )
        echo "$survivedR out of $totR paired reads survived de-duplication ($percentage%, $remR reads removed)"
        echo " "
        #Measures and logs execution time
        endtime=$(date +%s.%N)
        exectime=$(echo "$endtime $starttime" | awk '{print $1-$2}')
        sysdate=$(date)
        echo "STEP 1 (Quality control) terminated at $sysdate ($exectime seconds)"
        echo " "
        echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        echo " "


    >>>

    runtime {
        docker: docker_image
        cpu: cpus
        memory: memory_amount
        disks: "local-disk 20 SSD"
    }

    output {
        File dedupped_reads1 = prefix + "_dedupe_R1.fq"
        Array[File] dedupped_reads2 = glob(prefix + "_dedupe_R2.fq")
    }

}

task Trim {
    input {
        File reads1_file
        File? reads2_file
        String library_layout
        String prefix
        Int kcontaminants = 23
        Int mink = 11
        Int hdist = 1
        Int phred = 10
        Int min_length = 60
        Int qin = 33
        File adapters_file
        File artifacts_file
        File phix174ill_file
        String docker_image
    }

    meta {
        descritpion: "Trimming of low quality bases and of adapter sequences. Short reads are discarded. A decontamination of synthetic sequences is also pefoermed. If dealing with paired-end reads, when either forward or reverse of a paired-read are discarded, the surviving read is saved on a file of singleton reads."
    }

    parameter_meta {
        reads1_file: {
            description: "Input reads in fastq or fastq.gz format",
            extension: ['fq']
        }
        reads2_file: {
            description: "Paired read mates in fastq or fastq.gz format",
            extension: ['fq']
        }
        library_layout: "One of either 'paired' or 'single' describing the input reads"
        prefix: "Filename prefix"
        kcontaminants: "Kmer length used for finding contaminants"
        mink: "shorter kmers at read tips to look for"
        hdist: "maximum Hamming distance for ref kmers "
        phred: "regions with average quality BELOW this will be trimmed "
        min_length: "reads shorter than this after trimming will be discarded"
        qin: "Input quality offset: 33 (ASCII+33) or 64 (ASCII+64)"
        adapters_file: {
            description: "Adapter sequences to be trimmed in fa format",
            extension: ['fa']
        }
        artifacts_file: {
            description: "Sequencing artificats to be trimmed in fa.gz format",
            extension: ['fa.gz']
        }
        phix174ill_file: {
            description: "Processing contaminants to be trimemd in fa.gz format",
            extension: ['fa.gz']
        }
        docker_image:"Name of docker image to be used"

    }

    #String cpus = "4"
    String cpus = "1"
    #String memory_amount = "32 GB"
    String memory_amount = "8 GB"


    command <<<
        #Measures execution time
        sysdate=$(date)
        starttime=$(date +%s.%N)
        echo "Performing Quality Control. STEP 2 [Trimming] at $sysdate"
        echo " "
        #Sets the maximum memory to the value requested in the config file
        maxmem=$(echo ~{memory_amount} | sed 's/ //g' | sed 's/B//g')
        #Defines command for trimming of adapters and low quality bases
        if [ ~{library_layout} = "paired" ]; then
            #CMD="bbduk.sh -Xmx"$maxmem" in=~{reads1_file} in2=~{reads2_file} out=~{prefix}_trimmed_R1_tmp.fq out2=~{prefix}_trimmed_R2_tmp.fq outs=~{prefix}_trimmed_singletons_tmp.fq ktrim=r k=~{kcontaminants} mink=~{mink} hdist=~{hdist} qtrim=rl trimq=~{phred}  minlength=~{min_length} ref=~{adapters_file} qin=~{qin} threads=~{cpus} tbo tpe ow"
            #CMD="bbduk.sh in=~{reads1_file} in2=~{reads2_file} out=~{prefix}_trimmed_R1_tmp.fq out2=~{prefix}_trimmed_R2_tmp.fq outs=~{prefix}_trimmed_singletons_tmp.fq ktrim=r k=~{kcontaminants} mink=~{mink} hdist=~{hdist} qtrim=rl trimq=~{phred}  minlength=~{min_length} ref=~{adapters_file} qin=~{qin} threads=~{cpus} tbo tpe ow"
            CMD="bbduk.sh in=~{reads1_file} in2=~{reads2_file} out=~{prefix}_trimmed_R1_tmp.fq out2=~{prefix}_trimmed_R2_tmp.fq outs=~{prefix}_trimmed_singletons_tmp.fq ktrim=r k=~{kcontaminants} mink=~{mink} hdist=~{hdist} qtrim=rl trimq=~{phred}  minlength=~{min_length} ref=~{adapters_file} qin=~{qin} tbo tpe ow"
        else
            #CMD="bbduk.sh -Xmx"$maxmem" in=$reads1 out=~{prefix}_trimmed_tmp.fq ktrim=r k=~{kcontaminants} mink=~{mink} hdist=~{hdist} qtrim=rl trimq=~{phred}  minlength=~{min_length} ref=~{adapters_file} qin=~{qin} threads=~{cpus} tbo tpe ow"
            #CMD="bbduk.sh in=$reads1 out=~{prefix}_trimmed_tmp.fq ktrim=r k=~{kcontaminants} mink=~{mink} hdist=~{hdist} qtrim=rl trimq=~{phred}  minlength=~{min_length} ref=~{adapters_file} qin=~{qin} threads=~{cpus} tbo tpe ow"
            CMD="bbduk.sh in=$reads1 out=~{prefix}_trimmed_tmp.fq ktrim=r k=~{kcontaminants} mink=~{mink} hdist=~{hdist} qtrim=rl trimq=~{phred}  minlength=~{min_length} ref=~{adapters_file} qin=~{qin} tbo tpe ow"
        fi

        #Logs version of the software and executed command (BBMap prints on stderr)
        version=$(bbduk.sh --version 2>&1 >/dev/null | grep "BBMap version")
        echo "Using bbduk.sh in $version "
        echo "Using adapters in $params.adapters "
        echo "Using synthetic contaminants in ~{phix174ill_file} and in ~{artifacts_file} "
        echo " "
        echo "Executing command: $CMD "
        echo " "

        #Trims adapters and low quality bases
        exec $CMD 2>&1 | tee tmp.log

        #Logs some figures about sequences passing trimming
        echo  "BBduk's trimming stats (trimming adapters and low quality reads): "
        echo sed -n '/Input:/,/Result:/p' tmp.log
        echo " "
        if [ ~{library_layout} = "paired" ]; then
            unpairedR=$(wc -l ~{prefix}_trimmed_singletons_tmp.fq | cut -d" " -f 1)
            unpairedR=$((unpairedR/4))
            echo  "$unpairedR singleton reads whose mate was trimmed shorter preserved"
            echo " "
        fi
        #Defines command for removing synthetic contaminants
        if [ ~{library_layout} = "paired" ]; then
            CMD="bbduk.sh -Xmx"$maxmem" in=~{prefix}_trimmed_R1_tmp.fq in2=~{prefix}_trimmed_R2_tmp.fq out=~{prefix}_trimmed_R1.fq out2=~{prefix}_trimmed_R2.fq k=31 ref=~{phix174ill_file},~{artifacts_file} qin=~{qin} threads=~{cpus} ow"
        else
            CMD="bbduk.sh -Xmx"$maxmem" in=~{prefix}_trimmed_tmp.fq out=~{prefix}_trimmed.fq k=31 ref=~{phix174ill_file},~{artifacts_file} qin=~{qin} threads=~{cpus} ow"
        fi
        #Logs executed command
        echo "Executing command: $CMD "
        echo " "

        #Removes synthetic contaminants
        exec $CMD 2>&1 | tee tmp.log
        #Logs some figures about sequences passing deletion of contaminants
        echo  "BBduk's trimming stats (synthetic contaminants): "
        echo sed -n '/Input:/,/Result:/p' tmp.log
        echo " "
        #Removes synthetic contaminants and logs some figures (singleton read file,
        #that exists iif the library layout was 'paired')
        if [ ~{library_layout} = "paired" ]; then
            #CMD="bbduk.sh -Xmx"$maxmem" in=~{prefix}_trimmed_singletons_tmp.fq out=~{prefix}_trimmed.fq k=31 ref=~{phix174ill_file},~{artifacts_file} qin=~{qin} threads=~{cpus} ow"
            #CMD="bbduk.sh in=~{prefix}_trimmed_singletons_tmp.fq out=~{prefix}_trimmed.fq k=31 ref=~{phix174ill_file},~{artifacts_file} qin=~{qin} threads=~{cpus} ow"
            CMD="bbduk.sh in=~{prefix}_trimmed_singletons_tmp.fq out=~{prefix}_trimmed.fq k=31 ref=~{phix174ill_file},~{artifacts_file} qin=~{qin} ow"

            echo "Executing command: $CMD "
            echo " "

            #Removes synthetic contaminants
        exec $CMD 2>&1 | tee tmp.log

            #Logs some figures about sequences passing deletion of contaminants
            echo  "BBduk's trimming stats (synthetic contaminants, singleton reads): "
            echo sed -n '/Input:/,/Result:/p' tmp.log
        echo " "
        fi

        #Removes tmp files. This avoids adding them to the output channels
        rm -rf ~{prefix}_trimmed*_tmp.fq
        #Measures and log execution time
        endtime=$(date +%s.%N)
        exectime=$(echo "$endtime $starttime" | awk '{print $1-$2}')
        sysdate=$(date)
        echo "STEP 2 (Quality Control) terminated at $sysdate ($exectime seconds)"
        echo " "
        echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        echo ""
    >>>

    runtime {
        docker: docker_image
        cpu: cpus
        memory: memory_amount
        disks: "local-disk 20 SSD"
    }

    output {
        File trimmed_reads = prefix + "_trimmed.fq"
        File trimmed_reads_r1 = prefix + "_trimmed_R1.fq"
        Array[File] trimmed_reads_r2 = glob(prefix + "_trimmed_R2.fq")
    }

}

task Decontaminate {
    input{
        File infile1
        File infile2
        File infile12
        File ref_foreign_genome_gz
        String library_layout
        String prefix
        Float mind
        Int max_indel
        Int phred
        Int qin
        Float bwr
        String docker_image
    }

    meta {
        description: "Decontamination. Removes external organisms' contamination, using a previously created index. When paired-end are used, decontamination is carried on idependently on paired reads and on singleton reads thanks to BBwrap, that calls BBmap once on the paired reads and once on the singleton ones, merging the results on a single output file. Two files are outputted: the FASTQ of the decontaminated reads (including both paired-reads and singletons) and that of the contaminating reads (that can be used for refinements/checks). Please note that if keepQCtmpfile is set to false, the file of the contaminating reads is discarded"
    }

    parameter_meta {
        infile1: {
            description: "Input reads in fastq or fastq.gz format",
            extension: ['fq']
        }
        infile12: {
            description: "Singleton Input reads in fastq or fastq.gz format",
            extension: ['fq']
        }
        infile2: {
            description: "Input reads in fastq or fastq.gz format",
            extension: ['fq']
        }
        ref_foreign_genome_gz: {
            description: "Reference pan-genome for contamination. It should have been indexed beforehand.",
            extension: ['fq']
        }
        library_layout: "One of either 'paired' or 'single' describing the input reads"
        prefix: "Filename prefix"
        mind: "Approximate minimum alignment identity to look for"
        max_indel: "longest indel to look for"
        phred: "regions with average quality BELOW this will be trimmed "
        qin: "Input quality offset: 33 (ASCII+33) or 64 (ASCII+64)"
        bwr: "restrict alignment band to this"
    }

    #String cpus = "4"
    String cpus = "1"
    #String memory_amount = "32 GB"
    String memory_amount = "8 GB"

    command <<<
        sysdate=$(date)
        starttime=$(date +%s.%N)
        echo "Performing Quality Control. STEP 3 [Decontamination] at $sysdate"
        echo " "
        #Sets the maximum memory to the value requested in the config file
        #maxmem=$(echo ~{memory_amount} | sed 's/ //g' | sed 's/B//g')
        tar -xvf ~{ref_foreign_genome_gz} -C ref_foreign_genome
        #Defines command for decontamination
        if [ ~{library_layout} = "paired" ]; then
            #CMD="bbwrap.sh  -Xmx"$maxmem" mapper=bbmap append=t in1=~{infile1},~{infile12} in2=~{infile2},null outu=~{prefix}_clean.fq outm=~{prefix}_cont.fq minid=~{mind} maxindel=~{max_indel} bwr=~{bwr} bw=12 minhits=2 qtrim=rl trimq=~{phred} path=ref_foreign_genome qin=~{qin} threads=~{cpus} untrim quickmatch fast ow"
            CMD="bbwrap.sh  mapper=bbmap append=t in1=~{infile1},~{infile12} in2=~{infile2},null outu=~{prefix}_clean.fq outm=~{prefix}_cont.fq minid=~{mind} maxindel=~{max_indel} bwr=~{bwr} bw=12 minhits=2 qtrim=rl trimq=~{phred} path=ref_foreign_genome qin=~{qin} threads=~{cpus} untrim quickmatch fast ow"
        else
            #CMD="bbwrap.sh  -Xmx"$maxmem" mapper=bbmap append=t in1=~{infile1} outu=~{prefix}_clean.fq outm=~{prefix}_cont.fq minid=~{mind} maxindel=~{max_indel} bwr=~{bwr} bw=12 minhits=2 qtrim=rl trimq=~{phred} path=ref_foreign_genome qin=~{qin} threads=~{cpus} untrim quickmatch fast ow"
            CMD="bbwrap.sh  mapper=bbmap append=t in1=~{infile1} outu=~{prefix}_clean.fq outm=~{prefix}_cont.fq minid=~{mind} maxindel=~{max_indel} bwr=~{bwr} bw=12 minhits=2 qtrim=rl trimq=~{phred} path=ref_foreign_genome qin=~{qin} threads=~{cpus} untrim quickmatch fast ow"
        fi

        #Logs version of the software and executed command (BBmap prints on stderr)
        version=$(bbwrap.sh --version 2>&1 >/dev/null | grep "BBMap version")
        echo "Using bbwrap.sh in $version "
        echo "Using contaminant (pan)genome indexed in ~{ref_foreign_genome_gz}"
        echo " "
        echo "Executing command: $CMD "
        echo " "

        #Decontaminates
        exec $CMD 2>&1 | tee tmp.log

        #Logs some figures about decontaminated/contaminated reads
        echo  "BBwrap's human decontamination stats (paired reads): "
        sed -n '/Read 1 data:/,/N Rate:/p' tmp.log | head -17
        echo " "
        sed -n '/Read 2 data:/,/N Rate:/p' tmp.log
        echo " "

        if [ ~{library_layout} = "paired" ]; then
            echo  "BBmap's human decontamination stats (singletons reads): "
            sed -n '/Read 1 data:/,/N Rate:/p' tmp.log | tail -17
            echo " "
        fi
    gzip -c ~{prefix}_clean.fq > ~{prefix}_clean.fq.gz
    nClean=$(wc -l ~{prefix}_clean.fq | cut -d" " -f 1)
    nClean=$((nClean/4))
    nCont=$(wc -l ~{prefix}_cont.fq | cut -d" " -f 1)
    nCont=$((nCont/4))
    echo "$nClean reads survived decontamination ($nCont reads removed)"
    echo " "
    #Measures and log execution time
    endtime=$(date +%s.%N)
    exectime=$(echo "$endtime $starttime" | awk '{print $1-$2}')
    sysdate=$(date)
    echo "STEP 3 (Quality Control) terminated at $sysdate ($exectime seconds)"
    echo " "
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo ""
    """
    >>>

    runtime {
        docker: docker_image
        cpu: cpus
        memory: memory_amount
        disks: "local-disk 20 SSD"
    }

    output {
        File cleaned_reads = prefix + "_clean.fq.gz"
        File contaminated_reads = prefix + "_cont.fq"
    }
}

task QualityAssessment {
    input {
        File reads
        String step
        String stem
        String label
        String prefix
        String docker_image

    }

    meta {
        description: "Assessment of read quality of FASTQ file, done by means of FastQC. Multiple  plots are generated to show average phred quality scores and other metrics. This step will generate (for each end, if the layout is paired) an HTML page, showing a summary of the results and a set of plots offering a visual guidance and for assessing the quality of the sample, and a zip file, that includes the HTML page, all the images and a text report. The script will save the HTML page and the text report and delete the archive. Several information on multiple QC parameters are logged."
    }

    parameter_meta {
        reads: {
            description: "Input reads in fastq or fastq.gz format",
            extension: ['fq.gz', 'fastq']
        }
        step: "Step of QC"
        stem: "R1/R2 reads"
        label: "Input reads label"
        prefix: "Filename prefix"
    }

    #String cpus = "4"
    String cpus = "1"
    #String memory_amount = "32 GB"
    String memory_amount = "8 GB"
    command <<<
        #Measures execution time
        sysdate=$(date)
        starttime=$(date +%s.%N)
        echo "Performing Quality Control. [Assessment of read quality] at $sysdate"
        echo "File being analysed: ~{reads}"
        echo " "

        #Logs version of the software and executed command
        version=$(fastqc --version)
        CMD="fastqc --quiet --noextract --format fastq --outdir=. --threads ~{cpus} ~{reads}"

        echo "Using $version "
        echo "Executing command $CMD "
        echo " "

        #Does QC, extracts relevant information, and removes temporary files
        bash fastQC.sh ~{reads} ~{prefix}~{stem}~{label} ~{cpus} ~{reads}

        #Logging QC statistics (number of sequences, Pass/warning/fail, basic statistics, duplication level, kmers)
        base=$(basename ~{reads})
        #change last arg to standard out
        bash logQC.sh $base ~{prefix}~{stem}~{label}_fastqc_data.txt .log.$step$label

        #Measures and log execution time
        endtime=$(date +%s.%N)
        exectime=$(echo "$endtime $starttime" | awk '{print $1-$2}')
        sysdate=$(date)
        echo "Quality assessment on ~{reads} terminated at $sysdate ($exectime seconds)"
        echo " " >>
        echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        echo " "
        """
    >>>

    runtime {
        docker: docker_image
        cpu: cpus
        memory: memory_amount
        disks: "local-disk 20 SSD"
    }

    output {
        File html_file = prefix + stem + label + ".html"
        File fastqc_data_file = prefix + stem + label + ".zip"
    }
}

task ProfileTaxa {
    input {
        File in_file
        File mpa_pkl_file
        #Care this looks like a Dir, maybe provide a tar.gz
        File bowtie2db_gz_file
        String bowtiedb_files_name
        String bt2_options
        String prefix
        String docker_image
    }

    meta {
        description: "Community Characterisation - STEP 1. Performs taxonomic binning and estimates the microbial relative abundancies. MetaPhlAn2 and its databases of clade-specific markers are used to infers the presence and relative abundance of the organisms (at the specie/ strain level) that are present in the sample and to estimate their relative abundance."
    }

    parameter_meta {
        in_file: {
            description: "Input reads in fastq or fastq.gz format",
            extension: ['fq.gz', 'fastq']
        }
        mpa_pkl_file: {
            description: "BowTie2 database for MetaPhlAn2",
            extension: ['pkl']
        }
        bowtie2db_gz_file: {
            description: "BowTie2 resources for MetaPhlAn2",
            extension: ['gz']
        }
        bowtiedb_files_name: "BowTie2 resource name"
        bt2_options: "presets options for BowTie2"
        prefix: "Filename prefix"

    }

    #String cpus = "4"
    String cpus = "1"
    #String memory_amount = "32 GB"
    String memory_amount = "8 GB"

    #unzip bowtie2db to a folder and pass as arg
    command <<<
        sysdate=$(date)
        starttime=$(date +%s.%N)
        echo "Performing Community Characterisation. STEP 1 [Taxonomic binning and profiling] at $sysdate"
        echo " "

        tar -xvf ~{bowtie2db_gz_file} -C bowtie2db

        #Defines command for estimating abundances
        CMD="metaphlan2.py --input_type fastq --tmp_dir=. --biom ~{prefix}.biom --bowtie2out=~{prefix}_bt2out.txt --mpa_pkl ~{mpa_pkl_file}  --bowtie2db bowtie2db/~{bowtiedb_files_name} --bt2_ps ~{bt2_options} --nproc ~{cpus} ~{in_file} ~{prefix}_metaphlan_bugs_list.tsv"
        #Logs version of the software and executed command
        #MetaPhlAn prints on stderr
        version=$(metaphlan2.py --version 2>&1 >/dev/null | grep "MetaPhlAn")
        echo "Using $version "
        echo "Using BowTie2 database in $params.bowtie2db "
        echo " "
        echo "Executing command: $CMD "

        #Estimates microbial abundances
        exec $CMD 2>&1 | tee tmp.log
        #Sets the prefix in the biom file
        sed -i 's/Metaphlan2_Analysis/~{prefix}/g' ~{prefix}.biom
        sed -i 's/Metaphlan2_Analysis/~{prefix}/g' ~{prefix}_metaphlan_bugs_list.tsv
        #Logs some info
        tree=(kingdom phylum class order family genus species)
        for i in {2..7}
        do
            c=$(sed '1d' ~{prefix}_metaphlan_bugs_list.tsv | cut -d"|" -f $i | grep -v "k__" | cut -f 1  | sort | uniq | sed '/^\\s*$/d' | wc -l | cut -d" " -f 1)
            echo "$c ${tree[((i-1))]} found"
        done
        #Measures and log execution time
        endtime=$(date +%s.%N)
        exectime=$(echo "$endtime $starttime" | awk '{print $1-$2}')
        sysdate=$(date)
        echo ""
        echo "STEP 1 (Community Characterisation) terminated at $sysdate ($exectime seconds)"
        echo " "

        echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        echo ""
    >>>

    runtime {
        docker: docker_image
        cpu: cpus
        memory: memory_amount
        disks: "local-disk 20 SSD"
    }

    output {
        File biom_file = prefix + ".biom"
        File metaphlan_file = prefix + "_metaphlan_bugs_list.tsv"
        File bt2out_file = prefix + "_bt2out.txt"
    }
}

task AlphaDiversity {
    input{
        File biom_file
        File? tree_path_file
        String prefix
        String docker_image
    }

    meta {
        description: "Community Characterisation - STEP 2. Evaluates alpha-diversity, that is the mean species diversity the given sample. Please note that the alpha diversity is the only per-sample measure, so it is the only one evaluated by this module. If a newick tree is provided as input (see QIIME documentation for details), a further and more reliable phylogenetic measure is evaluated (i.e., PD_whole_tree)."
    }

    parameter_meta {
        biom_file: {
            description: "Biom file",
            extension: ['biom']
        }
        tree_path_file: {
            description: "Newick tree filepath, required for phylogenetic alpha diversity (PD_whole_tree, QIIME)",
            extension: ['fq.gz', 'fastq']
        }
        prefix: "Filename prefix"
    }

   #String cpus = "4"
    String cpus = "1"
    #String memory_amount = "32 GB"
    String memory_amount = "8 GB"

    command <<<
        #Measures execution time
        sysdate=$(date
        starttime=$(date +%s.%N)
        echo "Performing Community Characterisation. STEP 2 [Evaluating alpha-diversity] at $sysdate"
        echo " "
        #It checks if the profiling was successful, that is if identifies at least three species
        n=$(grep -o s__ ~{biom_file} | wc -l  | cut -d" " -f 1)
        if (( n > 3 ))
        then
        #Defines command -- if the tree path is not specified, not all the alpha
        #measures can be evaluated (that is, PD_whole_tree is skipped)
            if [ ~{tree_path_file} == null ]
            then
                CMD="alpha_diversity.py -i ~{biom_file} -o ~{prefix}_alpha_diversity.tsv -m ace,berger_parker_d,brillouin_d,chao1,chao1_ci,dominance,doubles,enspie,equitability,esty_ci,fisher_alpha,gini_index,goods_coverage,heip_e,kempton_taylor_q,margalef,mcintosh_d,mcintosh_e,menhinick,michaelis_menten_fit,observed_otus,observed_species,osd,simpson_reciprocal,robbins,shannon,simpson,simpson_e,singles,strong"
            else
                CMD="alpha_diversity.py -i ~{biom_file} -o ~{prefix}_alpha_diversity.tsv -m ace,berger_parker_d,brillouin_d,chao1,chao1_ci,dominance,doubles,enspie,equitability,esty_ci,fisher_alpha,gini_index,goods_coverage,heip_e,kempton_taylor_q,margalef,mcintosh_d,mcintosh_e,menhinick,michaelis_menten_fit,observed_otus,observed_species,osd,simpson_reciprocal,robbins,shannon,simpson,simpson_e,singles,strong,PD_whole_tree -t ~{tree_path_file}"
            fi

            #Logs version of the software and executed command
            version=$(alpha_diversity.py --version)
            echo "Using $version "
            if [ ~{tree_path_file} == null ]
            then
                echo "Newick tree not used, PD_whole_tree skipped"
            else
                echo "Using Newick tree in ~{tree_path_file}"
            fi
            echo " "
            echo "Executing command: $CMD "
            echo " "

            #Evaluates alpha diversities, redirect is done here because QIIME gets it as an extra parameter
            exec $CMD 2>&1 | tee tmp.log
        else
            #Also if the alpha are not evaluated the file should be created in order to be returned
            echo "Not enough classified species detected (N=$n). Analysis skipped."
            touch ~{prefix}_alpha_diversity.tsv
        fi
        #Measures and log execution time
        endtime=$(date +%s.%N)
        exectime=$(echo "$endtime $starttime" | awk '{print $1-$2}')
        sysdate=$(date)
        echo ""
        echo "STEP 2 (Community Characterisation) terminated at $sysdate ($exectime seconds)"
        echo " "
        echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        echo ""
    >>>

    runtime {
        docker: docker_image
        cpu: cpus
        memory: memory_amount
        disks: "local-disk 20 SSD"
    }
    output {
        File alpha_diversity_file = prefix + "_alpha_diversity.tsv"
    }
}

task ProfileFunction{
    input {
        File clean_reads_file
        File metaphlan_file
        File chocophlan_gz_file
        File uniref_file
        String prefix
        String docker_image
    }

    meta {
        description: "Community Characterisation - STEP 3. Performs the functional annotation using HUMAnN2. HUMAnN2 will bypasses the taxomonic profiling step (since it has already been performed) and uses the list of specied detected on step 7. While the aligners are forced to be Bowtie2 and DIAMOND, the user can select the UniRef database to use (UniRef50, UniRef90)."
    }

    parameter_meta {
        clean_reads_file: {
            description: "Input reads in fastq format",
            extension: ['fq', 'fastq', 'fq.gz', 'fastq.gz']
        }
        metaphlan_file: {
            description: "Metaphlan2 output file",
            extension: ['tsv']
        }
        chocophlan_gz_file: {
            description: "Gzip of Chocophlan files",
            extension: ['.gz']
        }
        uniref_file: {
            description: "Uniref file",
            extension: ['.dmnd']
        }
    }

    #String cpus = "4"
    String cpus = "1"
    #String memory_amount = "32 GB"
    String memory_amount = "8 GB"
    #deal with gz files can we pas the correct files or do we pass the whole directory?
    command <<<
        #Measures execution time
        sysdate=$(date)
        starttime=$(date +%s.%N)
        echo "Performing Community Characterisation. STEP 3 [Performing functional annotation] with HUMAnN2 at $sysdate"
        echo " "
        tar -xvf ~{chocophlan_gz_file} -C chocophlan
        #Defines HUMAnN2 command taking advantages of the MetaPhlAn2's results
        CMD="humann2 --input ~{clean_reads_file} --output . --output-basename ~{prefix} --taxonomic-profile ~{metaphlan_file} --nucleotide-database chocophlan --protein-database ~{uniref_file} --pathways metacyc --threads ~{cpus} --memory-use minimum"

        #Logs version of the software and executed command
        #HUMAnN2 prints on stderr
        version=$(humann2 --version 2>&1 >/dev/null | grep "humann2")
        echo "Using $version "
        echo "Using ChocoPhlAn database in ~{chocophlan_gz_file} "
        echo "Using UniRef database in ~{uniref_file} "
        echo " "
        echo "Executing command: $CMD > ~{prefix}_HUMAnN2.log"
        echo " "

        #Performs functional annotation, redirect is done here because HUMAnN2 freaks out
        #This is  also reported in the log.
        exec $CMD 2>&1 | tee ~{prefix}_HUMAnN2.log
        #If `|| true` is not add, nextflow stops... WTF
        grep "Total species selected from prescreen:" ~{prefix}_HUMAnN2.log
        grep "Selected species explain" ~{prefix}_HUMAnN2.log
        grep "Unaligned reads after nucleotide alignment:" ~{prefix}_HUMAnN2.log
        grep "Total gene families after translated alignment:" ~{prefix}_HUMAnN2.log
        grep "Unaligned reads after translated alignment:" ~{prefix}_HUMAnN2.log
        echo "More information on HUMAnN2 run are available in the ~{prefix}_HUMAnN2.log file"
        #Some of temporary files (if they exist) may be moved in the working directory,
        #according to the keepCCtmpfile parameter. Others (such as the bowties2 indexes),
        #are always removed. Those that should be moved, but have not been created by
        #HUMAnN2, are now created by the script (they are needed as output for the channel)
        files=(~{prefix}_bowtie2_aligned.sam ~{prefix}_bowtie2_aligned.tsv ~{prefix}_diamond_aligned.tsv ~{prefix}_bowtie2_unaligned.fa ~{prefix}_diamond_unaligned.fa)
        for i in {1..5}
        do
            if [ -f ~{prefix}_humann2_temp/${files[((i-1))]} ]
            then
                mv ~{prefix}_humann2_temp/${files[((i-1))]} .
            else
                touch ${files[((i-1))]}
            fi
        done
        rm -rf ~{prefix}_humann2_temp/
        #Measures and log execution time
        endtime=$(date +%s.%N)
        exectime=$(echo "$endtime $starttime" | awk '{print $1-$2}')
        sysdate=$(date)
        echo ""
        echo "STEP 3 (Community Characterisation) terminated at $sysdate ($exectime seconds)"
        echo " "
        echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        echo ""
    >>>

    runtime {
        docker: docker_image
        cpu: cpus
        memory: memory_amount
        disks: "local-disk 20 SSD"
    }

    output {
        File HUMAnN2_log_file = prefix + "_HUMAnN2.log"
        File gene_families_file = prefix + "_genefamilies.tsv"
        File path_coverage_file = prefix + "_pathcoverage.tsv"
        File path_abundance_file = prefix + "_pathabundance.tsv"
        File bowtie2_aligned_sam_file =  prefix + "_bowtie2_aligned.sam"
        File bowtie2_aligned_tsv_file = prefix + "_bowtie2_aligned.tsv"
        File diamond_aligned_tsv_file = prefix + "_diamond_aligned.tsv"
        File bowtie2_unaligned_file = prefix + "_bowtie2_unaligned.fa"
        File diamond_unaligned_file = prefix + "_diamond_unaligned.fa"
    }
}