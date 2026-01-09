#!/bin/bash

# Pre-requisit: 
# GNU parallel (https://www.gnu.org/software/parallel/parallel_tutorial.html)
# pigz (https://github.com/madler/pigz)
# cutadapt (https://cutadapt.readthedocs.io/en/stable/index.html)
# umi_tools (https://umi-tools.readthedocs.io/en/latest/INSTALL.html)
# STAR (https://github.com/alexdobin/STAR)
# featureCounts (subreads: https://subread.sourceforge.net/subread-package.html)
# bedtools ()
# macs2
# SICER2 (https://zanglab.github.io/SICER2/)
# in-house tools:
# tag_CB_UB.py


# linker2:
# linker1:
# Adapter:
# Library:
# DNA:
# Library:
# RNA:

toolPATH=/data02/laibb/scCOS-seq/tools
refPATH=/data02/laibb/annotation

CONFIG_FILE=$1

declare -A config  # Associative array for configuration

while IFS='=' read -r key value || [[ -n "$key" ]]; do
    # Skip empty lines and comment lines
    if [[ -z "$key" || "$key" =~ ^[[:space:]]*# ]]; then
        continue
    fi

    # Remove spaces before and after key and value
    key=$(echo "$key" | tr -d '\r\n' | xargs)
    value=$(echo "$value" | tr -d '\r\n' | xargs)

    echo "key: '$key', value: '$value'"

    config["$key"]="$value"

done < "$CONFIG_FILE"

# Check required keys
required_keys=("LibName" "SampleName" "Genome" "CUTTag"  "Samples" "Keys" "Dir_rawdata")
optional_keys=("Threads" "AlnThreads")
array_keys=("LibName" "SampleName" "CUTTag"  "Samples" "Keys")

echo -e "\nChecking required keys..."

# Check for missing keys
missing_keys=()
for key in "${required_keys[@]}"; do
#    if [[ -v config["$key"] ]]; then
    if [[ -n "${config[$key]+x}" ]]; then
        echo "Found required key '$key': ${config[$key]}"
    else
        echo "Missing required key '$key'"
        missing_keys+=("$key")
    fi
done

# Check optional keys
echo -e "\nChecking optional keys..."
for key in "${optional_keys[@]}"; do
    if [[ -n "${config[$key]+x}" ]]; then
        echo "Found optional key '$key': ${config[$key]}"
    else
        echo "Missing optional key '$key'"
    fi
done

# Exit if any required keys are missing
if [ ${#missing_keys[@]} -ne 0 ]; then
    echo -e "\nError: Missing following required keys: ${missing_keys[*]}"
    echo "Please ensure the required keys exist in the config file!"
    exit 1
fi

# Convert values to arrays
declare -A config_arrays

for key in "${array_keys[@]}"; do
    if [[ -n "${config[$key]+x}" ]]; then
        echo "Processing key '$key': ${config[$key]}"
        
        # Get original value
        original_value="${config[$key]}"
        
        # Use IFS to convert comma-separated values to array
        IFS=',' read -ra value_array <<< "$original_value"
        
        # Remove spaces before and after each element
        for i in "${!value_array[@]}"; do
            value_array[$i]=$(echo "${value_array[$i]}" | xargs)
        done
        
        # Store in config_arrays associative array
        config_arrays["$key"]=${value_array[@]}
        
        # Print the converted array
        echo "Converted array:"
        for i in "${!value_array[@]}"; do
            echo "  $i: ${value_array[$i]}"
        done
        echo ""
    else
        echo "Key '$key' does not exist."
    fi
done

# Extract array variables
lib=()
if [[ -n "${config_arrays["LibName"]+x}" ]]; then
    lib=(${config_arrays["LibName"]})
    echo "LibName array: ${lib[@]}"
else
    echo "LibName not found in config_arrays"
fi

ATAC=()
if [[ -n "${config_arrays["SampleName"]+x}" ]]; then
    ATAC=(${config_arrays["SampleName"]})
    echo "SampleName array: ${ATAC[@]}"
else
    echo "SampleName not found in config_arrays"
fi



CUTTag=()
if [[ -n "${config_arrays["CUTTag"]+x}" ]]; then
    CUTTag=(${config_arrays["CUTTag"]})
    echo "CUTTag array: ${CUTTag[@]}"
else
    echo "CUTTag not found in config_arrays"
fi

#withCUTTag=()
#if [[ -n "${config_arrays["WithCUTTag"]+x}" ]]; then
#    withCUTTag=(${config_arrays["WithCUTTag"]})
#    echo "WithCUTTag array: ${withCUTTag[@]}"
#else
#    echo "WithCUTTag not found in config_arrays"
#fi

samples=()
if [[ -n "${config_arrays["Samples"]+x}" ]]; then
    samples=(${config_arrays["Samples"]})
    echo "Samples array: ${samples[@]}"
else
    echo "Samples not found in config_arrays"
fi

keys=()
if [[ -n "${config_arrays["Keys"]+x}" ]]; then
    keys=(${config_arrays["Keys"]})
    echo "Keys array: ${keys[@]}"
else
    echo "Keys not found in config_arrays"
fi

threads=${config["Threads"]:-20}
echo "Threads: $Threads"
#addthreads=$((threads-1))

AlnThreads=${config["AlnThreads"]:-5}
echo "AlnThreads: $AlnThreads"

Dir_rawdata=${config["Dir_rawdata"]}
echo "Dir_rawdata: $Dir_rawdata"

Genome=${config["Genome"]}
echo "Genome: $Genome"

if [[ $Genome =~ ^mm.* ]]; then
    starGenome=/data02/laibb/annotation/genome/refdata-gex-mm10-2020-A/star_index
    if [[ ! -s $starGenome ]]; then
        echo "error: cannot find refGenome for STAR: $starGenome "
    fi
    gtf=/data02/laibb/annotation/gtf/mm10.gencode.gtf
    msk=/data02/laibb/annotation/gtf/mm10_rmsk.gtf
    ChromSizes=/data02/laibb/annotation/genome/refdata-gex-mm10-2020-A/fasta/mm10.main.chrom.sizes
    MacsGenome=mm
 #   gtf=/data02/laibb/annotation/gtf/refGene_mm10.gtf.sorted
    if [[ ! -s $gtf ]]; then
        echo "error: cannot find GTF file for featureCounts: $gtf "
    fi
    if [[ ! -s $msk ]]; then
        echo "error: cannot find MSK GTF file for featureCounts: $msk "
    fi
elif [[ $Genome =~ ^hg.* ]]; then
    starGenome=/data02/laibb/annotation/genome/refdata-gex-GRCh38-2024-A/star_index
    if [[ ! -s $starGenome ]]; then
        echo "error: cannot find refGenome for STAR: $starGenome "
    fi
    gtf=/data02/laibb/annotation/gtf/hg38.gencode.gtf
    msk=/data02/laibb/annotation/gtf/hg38_rmsk.gtf
    ChromSizes=/data02/laibb/annotation/genome/refdata-gex-GRCh38-2024-A/fasta/hg38.main.chrom.sizes
    MacsGenome=hs
    if [[ ! -s $gtf ]]; then
        echo "error: cannot find GTF file for featureCounts: $gtf "
    fi
    if [[ ! -s $msk ]]; then
        echo "error: cannot find MSK GTF file for featureCounts: $msk "
    fi
elif [[ $Genome =~ ^combined* ]]; then
    starGenome=/data02/laibb/annotation/genome/refdata-gex-GRCh38-and-mm10-2020-A/star_index
    if [[ ! -s $starGenome ]]; then
        echo "error: cannot find refGenome for STAR: $starGenome "
    fi

fi

#exit
### Demultiplex 3 rounds barcode
Adaptor="V3" # V1: 8bp sample barcode; V2: 6bp sample barcode; V3: sample_Ad39/linker_27

if [ ! -s Round2bc.fa ]; then
    python3 $toolPATH/generateRound2bcV2.py  #build 'Round2bc.fa' file
fi
if [ ! -s Round1bc.fa ]; then
    python3 $toolPATH/generateRound1bc.py  #build 'Round1bc.fa' file
fi
if [ ! -s Samplebc.txt ]; then
    if [ $Adaptor == "V1" ]; then
        python3 $toolPATH/generateSamplebcV1.py
    elif [ $Adaptor == "V2" ]; then
        python3 $toolPATH/generateSamplebcV2.py
    elif [ $Adaptor == "V3" ]; then
        python3 $toolPATH/generateSamplebcV3.py
    fi  #build 'Samplebc_all.fa' file

    for (( i=0; i<48; i++ )); do
        if [ ${samples[$i]} != ${keys[$i]} ]; then
            sed -i "s/${keys[$i]}/${samples[$i]}/" Samplebc_all.fa
        fi
    done

    cat Samplebc_all.fa | grep -v sample | grep -A1 ">" | grep -v "-" > Samplebc.fa
    cat Samplebc.fa | grep -v "\^" | sed 's/>//' > Sample.txt
    rm Samplebc_all.fa

    checksample=()
    while read line; do
        checksample+=("$line")
    done < Sample.txt

    if [[ `echo ${ATAC[@]} ${CUTTag[@]} ${checksample[@]} | sed 's/ /\n/g' | sort | uniq -u` ]]; then
        echo "samples are not set properly, exit..." && exit
    else
        echo "all samples set properly"
    fi
fi


for (( i=0; i< ${#lib[@]}; i++ )); do
    echo ">>> Processing ${lib[$i]}"

    
    if [ -s ${lib[$i]}/${lib[$i]}-header_1.fq.gz ] && [ -s ${lib[$i]}/${lib[$i]}-header_2.fq.gz ]; then
        echo "barcodes have been appended to header"
    else
        # Label Round3 barcode
        # Usually round3 barcodes IDs were stored in the raw fastq file names, such as sr02-R3.01_1.fq.gz to sr02-R3.96_1.fq.gz,
        # Thus R3.01 to R3.96 are the round3 barcode IDs.
        if [ -s ${lib[$i]}/${lib[$i]}-R3header_1.fq.gz ] && [ -s ${lib[$i]}/${lib[$i]}-R3header_2.fq.gz ]; then
            echo "R3 barcodes have been appended to header"
        else

            ls $Dir_rawdata/${lib[$i]}-R3.*_1.fq.gz | sed 's/\///g' | sed 's/^.*-//' | sed 's/_1.fq.gz//' > ${lib[$i]}/Round3Label.txt
            nosublib=$(cat ${lib[$i]}/Round3Label.txt | wc -l)
            if [ -s ${lib[$i]}/Round3Label.txt ] && [ $nosublib -le 96 ]; then
                echo "label Round3 barcode"
                parallel -j $threads 'zcat '$Dir_rawdata'/'${lib[$i]}'-{}_1.fq.gz | sed "1~4 s/[ \t][1-2]:N:0.*/_{}/" | \
                            pigz --fast -p 1 > '${lib[$i]}'/'${lib[$i]}'-{}.header_1.fq.gz' :::: ${lib[$i]}/Round3Label.txt
                parallel -j $threads 'zcat '$Dir_rawdata'/'${lib[$i]}'-{}_2.fq.gz | sed "1~4 s/[ \t][1-2]:N:0.*/_{}/" | \
                            pigz --fast -p 1 > '${lib[$i]}'/'${lib[$i]}'-{}.header_2.fq.gz' :::: ${lib[$i]}/Round3Label.txt
            else
                echo "unexpected file detected in ${lib[$i]}, exiting..."
                exit
            fi
            cat ${lib[$i]}/${lib[$i]}-R3.*.header_1.fq.gz > ${lib[$i]}/${lib[$i]}-R3header_1.fq.gz
            cat ${lib[$i]}/${lib[$i]}-R3.*.header_2.fq.gz > ${lib[$i]}/${lib[$i]}-R3header_2.fq.gz
            rm ${lib[$i]}/${lib[$i]}-R3.*.header* 
        #    rm ${lib[$i]}/Round3Label.txt
            echo "total reads" 
            zcat ${lib[$i]}/${lib[$i]}-R3header_1.fq.gz | wc -l 
        fi
        

        # Demultiplex Round2 barcode
        if [ -s ${lib[$i]}/${lib[$i]}-R2header_1.fq.gz ] && [ -s ${lib[$i]}/${lib[$i]}-R2header_2.fq.gz ]; then
            echo "R2 barcodes have been appended to header"
        else
            echo "demultiplex Round2 barcode for ${lib[$i]}"
            cutadapt -e 1 -j 0 -g file:Round2bc.fa --no-indels \
                -o ${lib[$i]}/${lib[$i]}-{name}_2.fq.gz \
                -p ${lib[$i]}/${lib[$i]}-{name}_1.fq.gz \
                ${lib[$i]}/${lib[$i]}-R3header_2.fq.gz \
                ${lib[$i]}/${lib[$i]}-R3header_1.fq.gz
            mv ${lib[$i]}/${lib[$i]}-unknown_1.fq.gz ${lib[$i]}/${lib[$i]}-unknownR2_1.fq.gz 
            mv ${lib[$i]}/${lib[$i]}-unknown_2.fq.gz ${lib[$i]}/${lib[$i]}-unknownR2_2.fq.gz 
            # NOTE: the {name} indicates the barcode name, such as R2.01, R2.02, to R2.96.
            # There will be 96 output paired files such as ${lib[$i]}/${lib[$i]}-R2.01_.1.fq.gz

            echo "Unknown reads"
            zcat ${lib[$i]}/${lib[$i]}-unknownR2_2.fq.gz | wc -l
            
            # Processing linker2
            ls ${lib[$i]}/${lib[$i]}-R2.*_1.fq.gz | sed 's/^.*\///g' | sed 's/_1.fq.gz//' > ${lib[$i]}/Round2linker.txt
            parallel -j $threads 'cutadapt -e 0.2 -g ^AGTATGCAGCGCGCTCAAGCACGTGG \
                                -o '${lib[$i]}'/{}.trimL2_2.fq.gz -p '${lib[$i]}'/{}.trimL2_1.fq.gz \
                                '${lib[$i]}'/{}_2.fq.gz '${lib[$i]}'/{}_1.fq.gz' :::: ${lib[$i]}/Round2linker.txt
            # Label Round2 barcode
            ls ${lib[$i]}/${lib[$i]}-R2.*.trimL2_1.fq.gz | sed 's/\///g' | sed 's/^.*-//' | sed 's/.trimL2_1.fq.gz//' > ${lib[$i]}/Round2Label.txt
            echo "label Round2 barcode for ${lib[$i]}"
            parallel -j $threads 'if [ -s '${lib[$i]}'/'${lib[$i]}'-{}.trimL2_1.fq.gz ]; then \
                        zcat '${lib[$i]}'/'${lib[$i]}'-{}.trimL2_1.fq.gz | sed "1~4 s/^@.*/&,{}/" | \
                        pigz --fast -p 1 > '${lib[$i]}'/'${lib[$i]}'-{}.header_1.fq.gz; fi' :::: ${lib[$i]}/Round2Label.txt
            parallel -j $threads 'if [ -s '${lib[$i]}'/'${lib[$i]}'-{}.trimL2_2.fq.gz ]; then \
                        zcat '${lib[$i]}'/'${lib[$i]}'-{}.trimL2_2.fq.gz | sed "1~4 s/^@.*/&,{}/" | \
                        pigz --fast -p 1 > '${lib[$i]}'/'${lib[$i]}'-{}.header_2.fq.gz; fi' :::: ${lib[$i]}/Round2Label.txt
            cat ${lib[$i]}/${lib[$i]}-R2.*.header_1.fq.gz > ${lib[$i]}/${lib[$i]}-R2header_1.fq.gz
            cat ${lib[$i]}/${lib[$i]}-R2.*.header_2.fq.gz > ${lib[$i]}/${lib[$i]}-R2header_2.fq.gz
            rm ${lib[$i]}/${lib[$i]}-R2.*.fq.gz ${lib[$i]}/Round2linker.txt 
            rm ${lib[$i]}/Round2Label.txt
        fi
     #   exit

        # Demultiplex Round1 barcode
        # The process is similar with Round2 barcodes
        echo "demultiplex Round1 barcode for ${lib[$i]}"
        cutadapt -e 1 -j 0 -g file:Round1bc.fa --no-indels \
            -o ${lib[$i]}/${lib[$i]}-{name}_2.fq.gz \
            -p ${lib[$i]}/${lib[$i]}-{name}_1.fq.gz \
            ${lib[$i]}/${lib[$i]}-R2header_2.fq.gz \
            ${lib[$i]}/${lib[$i]}-R2header_1.fq.gz
        mv ${lib[$i]}/${lib[$i]}-unknown_1.fq.gz ${lib[$i]}/${lib[$i]}-unknowR1_1.fq.gz 
        mv ${lib[$i]}/${lib[$i]}-unknown_2.fq.gz ${lib[$i]}/${lib[$i]}-unknowR1_2.fq.gz 
            
        # Processing linker1
        ls ${lib[$i]}/${lib[$i]}-R1.*_1.fq.gz | sed 's/^.*\///g' | sed 's/_1.fq.gz//' > ${lib[$i]}/Round1linker.txt
        parallel -j $threads 'cutadapt -e 0.2 -g ^TCGTACGCCGATGCGAAACATCG \
                            -o '${lib[$i]}'/{}.trimL1_2.fq.gz -p '${lib[$i]}'/{}.trimL1_1.fq.gz \
                            '${lib[$i]}'/{}_2.fq.gz '${lib[$i]}'/{}_1.fq.gz' :::: ${lib[$i]}/Round1linker.txt
        # Label Round1 and P5 barcode
        P5Label="P1."${lib[$i]:2:2}
        ls ${lib[$i]}/${lib[$i]}-R1.*.trimL1_1.fq.gz | sed 's/\///g' | sed 's/^.*-//' | sed "s/.trimL1_1.fq.gz/\t${P5Label}/" > ${lib[$i]}/Round1P5Label.txt
        echo "label Round1 and P5 barcode"
        parallel -j $threads --colsep '\t' 'if [ -s '${lib[$i]}'/'${lib[$i]}'-{1}.trimL1_1.fq.gz ]; then \
                    zcat '${lib[$i]}'/'${lib[$i]}'-{1}.trimL1_1.fq.gz | sed "1~4 s/^@.*/&,{1},{2}/" | \
                    pigz --fast -p 1 > '${lib[$i]}'/'${lib[$i]}'-{1}.header_1.fq.gz; fi' :::: ${lib[$i]}/Round1P5Label.txt
        parallel -j $threads --colsep '\t' 'if [ -s '${lib[$i]}'/'${lib[$i]}'-{1}.trimL1_2.fq.gz ]; then \
                    zcat '${lib[$i]}'/'${lib[$i]}'-{1}.trimL1_2.fq.gz | sed "1~4 s/^@.*/&,{1},{2}/" | \
                    pigz --fast -p 1 > '${lib[$i]}'/'${lib[$i]}'-{1}.header_2.fq.gz; fi' :::: ${lib[$i]}/Round1P5Label.txt
        cat ${lib[$i]}/${lib[$i]}-R1.*.header_1.fq.gz > ${lib[$i]}/${lib[$i]}-header_1.fq.gz
        cat ${lib[$i]}/${lib[$i]}-R1.*.header_2.fq.gz > ${lib[$i]}/${lib[$i]}-header_2.fq.gz
        rm ${lib[$i]}/${lib[$i]}-R1.*.fq.gz ${lib[$i]}/${lib[$i]}-unknown*
        rm ${lib[$i]}/Round1linker.txt ${lib[$i]}/Round1P5Label.txt
        rm ${lib[$i]}/${lib[$i]}-R*header*

        echo "zcat ${lib[$i]}/${lib[$i]}-header_2.fq.gz | wc -l"
        zcat ${lib[$i]}/${lib[$i]}-header_2.fq.gz | wc -l
    fi


    if [ ! -d ${lib[$i]}/trimmed ]; then mkdir -p ${lib[$i]}/trimmed; fi
    
    if [[ ${lib[$i]} =~ ^sr.* ]]; then
        
        # Extract UMIs for RNA-seq
        
        if [ -s ${lib[$i]}/trimmed/${lib[$i]}_1.fq.gz ] && [ -s ${lib[$i]}/trimmed/${lib[$i]}_2.fq.gz ]; then
            echo "UMIs have been extracted, skip"
        else
            echo "extract UMI"
            umi_tools extract -I ${lib[$i]}/${lib[$i]}-header_2.fq.gz --bc-pattern=XXNNNNNNNNNN \
                --read2-in=${lib[$i]}/${lib[$i]}-header_1.fq.gz \
                --stdout=${lib[$i]}/trimmed/${lib[$i]}_2.umi.fq.gz \
                --read2-out=${lib[$i]}/trimmed/${lib[$i]}_1.umi.fq.gz \
                --log=${lib[$i]}/trimmed/${lib[$i]}.log
            cutadapt -e 1 -g ^GCTTTTTT --no-indels --action=none \
                --discard-untrimmed \
                -o ${lib[$i]}/trimmed/${lib[$i]}_2.fq.gz -p ${lib[$i]}/trimmed/${lib[$i]}_1.fq.gz \
                ${lib[$i]}/trimmed/${lib[$i]}_2.umi.fq.gz ${lib[$i]}/trimmed/${lib[$i]}_1.umi.fq.gz
            rm ${lib[$i]}/trimmed/${lib[$i]}_2.umi.fq.gz ${lib[$i]}/trimmed/${lib[$i]}_1.umi.fq.gz
       
        fi
    fi

    if [[ ${lib[$i]} =~ ^sa.* ]] || [[ ${lib[$i]} =~ ^sc.* ]]; then
        # Demultiplex sample 
        
        if [[ ${lib[$i]} =~ ^sa.* ]]; then
            # Process DNA (ATAC - R1. R2 depends on sample)
            echo "demultiplex DNA (ATAC - R1. R2 depends on sample)"
            if [ -s ${lib[$i]}/trimmed/${lib[$i]}-${ATAC[-1]}_1.fq.gz ] && [ -s ${lib[$i]}/trimmed/${lib[$i]}-${ATAC[-1]}_2.fq.gz ]; then
                echo "samples have been demultiplexed"
            else
                echo "demultiplex sample barcode for ${lib[$i]}"
                cutadapt -e 1 -g file:Samplebc.fa --no-indels \
                    -o ${lib[$i]}/trimmed/${lib[$i]}-{name}.tmp_2.fq.gz \
                    -p ${lib[$i]}/trimmed/${lib[$i]}-{name}.tmp_1.fq.gz \
                    ${lib[$i]}/${lib[$i]}-header_2.fq.gz \
                    ${lib[$i]}/${lib[$i]}-header_1.fq.gz
                parallel -j $threads 'cutadapt -e 0.2 -g ^ATGTGTATAAGAGACAG -G ^AGATGTGTATAAGAGACAG \
                                    -o '${lib[$i]}'/trimmed/'${lib[$i]}'-{}.trim5tmp_2.fq.gz \
                                    -p '${lib[$i]}'/trimmed/'${lib[$i]}'-{}.trim5tmp_1.fq.gz \
                                    '${lib[$i]}'/trimmed/'${lib[$i]}'-{}.tmp_2.fq.gz \
                                    '${lib[$i]}'/trimmed/'${lib[$i]}'-{}.tmp_1.fq.gz' :::: Sample.txt
                parallel -j $threads 'cutadapt -e 0.2 -m 20 --overlap=1 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
                                    -o '${lib[$i]}'/trimmed/'${lib[$i]}'-{}_2.fq.gz \
                                    -p '${lib[$i]}'/trimmed/'${lib[$i]}'-{}_1.fq.gz \
                                    '${lib[$i]}'/trimmed/'${lib[$i]}'-{}.trim5tmp_2.fq.gz \
                                    '${lib[$i]}'/trimmed/'${lib[$i]}'-{}.trim5tmp_1.fq.gz' :::: Sample.txt     
                rm -f *tmp*fq.gz 
            fi
        fi
        
        if [[ ${lib[$i]} =~ ^sc.* ]]; then
            # Process DNA (CUTTag - R1. R2 depends on sample)
            echo "demultiplex DNA (ATAC - R1. R2 depends on sample)"
            if [ -s ${lib[$i]}/trimmed/${lib[$i]}-${CUTTag[-1]}_1.fq.gz ] && [ -s ${lib[$i]}/trimmed/${lib[$i]}-${CUTTag[-1]}_2.fq.gz ]; then
                echo "samples have been demultiplexed"
            else
                echo "demultiplex sample barcode for ${lib[$i]}"
                cutadapt -e 1 -g file:Samplebc.fa --no-indels \
                    -o ${lib[$i]}/trimmed/${lib[$i]}-{name}.tmp_2.fq.gz \
                    -p ${lib[$i]}/trimmed/${lib[$i]}-{name}.tmp_1.fq.gz \
                    ${lib[$i]}/${lib[$i]}-header_2.fq.gz \
                    ${lib[$i]}/${lib[$i]}-header_1.fq.gz
                parallel -j $threads 'cutadapt -e 0.2 -g ^ATGTGTATAAGAGACAG \
                                    -o '${lib[$i]}'/trimmed/'${lib[$i]}'-{}.trim5tmp_2.fq.gz \
                                    -p '${lib[$i]}'/trimmed/'${lib[$i]}'-{}.trim5tmp_1.fq.gz \
                                    '${lib[$i]}'/trimmed/'${lib[$i]}'-{}.tmp_2.fq.gz \
                                    '${lib[$i]}'/trimmed/'${lib[$i]}'-{}.tmp_1.fq.gz' :::: Sample.txt
                parallel -j $threads 'cutadapt -e 0.2 -m 20 --overlap=1 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
                                    -o '${lib[$i]}'/trimmed/'${lib[$i]}'-{}_2.fq.gz \
                                    -p '${lib[$i]}'/trimmed/'${lib[$i]}'-{}_1.fq.gz \
                                    '${lib[$i]}'/trimmed/'${lib[$i]}'-{}.trim5tmp_2.fq.gz \
                                    '${lib[$i]}'/trimmed/'${lib[$i]}'-{}.trim5tmp_1.fq.gz' :::: Sample.txt    
                rm -f *tmp*fq.gz    
            fi
        fi
    fi

done



for (( i=0; i< ${#lib[@]}; i++ )); do

    # Proceed RNA-seq
    if [ ! -d ${lib[$i]}/results ]; then mkdir -p ${lib[$i]}/results; fi
        
    if [[ ${lib[$i]} =~ ^sr.* ]]; then
        # Align
        echo "map ${lib[$i]} to $Genome genome"
        
        if [[ -s ${lib[$i]}/results/${lib[$i]}.${Genome}.bam ]]; then
            echo "star alignment already done."
        else
            STAR --chimOutType WithinBAM \
                --runThreadN $threads \
                --genomeDir $starGenome \
                --readFilesIn ${lib[$i]}/trimmed/${lib[$i]}_1.fq.gz \
                --outFileNamePrefix ${lib[$i]}/results/${lib[$i]}.${Genome}. \
                --outFilterMultimapNmax 1 \
                --outFilterScoreMinOverLread 0.33 \
                --outFilterMatchNminOverLread 0.33 \
                --outFilterMismatchNoverLmax 0.06 \
                --outSAMattributes NH HI AS nM MD \
                --limitOutSJcollapsed 5000000 \
                --outSAMtype BAM Unsorted \
                --limitIObufferSize 400000000 400000000 \
                --outReadsUnmapped None \
                --readFilesCommand zcat

            samtools sort -@ $threads -m 4G ${lib[$i]}/results/${lib[$i]}.${Genome}.Aligned.out.bam > ${lib[$i]}/results/${lib[$i]}.${Genome}.bam
            samtools index -@ $threads ${lib[$i]}/results/${lib[$i]}.${Genome}.bam
            mv ${lib[$i]}/results/${lib[$i]}.${Genome}.Log.final.out ${lib[$i]}/results/${lib[$i]}.${Genome}.align.log
            rm ${lib[$i]}/results/${lib[$i]}.${Genome}.Aligned.out.bam
            rm ${lib[$i]}/results/${lib[$i]}.${Genome}*Log.progress.out ${lib[$i]}/results/${lib[$i]}.${Genome}*SJ.out.tab
            rm -r ${lib[$i]}/results/${lib[$i]}.${Genome}._STARtmp
        fi

        # Filter chromatins and prepare for velocyto
        if [[ ! -s ${lib[$i]}/results/${lib[$i]}.${Genome}.namesort.CBUB.bam ]]; then
            
            cd ${lib[$i]}/results #--------------------
            chrs=`samtools view -H ${lib[$i]}.${Genome}.bam | grep chr | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<6)print}'`
            samtools view -@ $threads -bS -F 256 -q 255 ${lib[$i]}.${Genome}.bam `echo $chrs` | \
                samtools sort -@ $threads -m 4G -n -o ${lib[$i]}.${Genome}.namesort.bam
            

            
            # Add CB and UB tag to bam
            python3 ${toolPATH}/tag_CB_UB.py ${lib[$i]}.${Genome}.namesort.bam ${lib[$i]}.${Genome}.namesort.CBUB.sam
            samtools view ${lib[$i]}.${Genome}.namesort.CBUB.sam -b -o ${lib[$i]}.${Genome}.namesort.CBUB.bam
            rm ${lib[$i]}.${Genome}.namesort.CBUB.sam ${lib[$i]}.${Genome}.namesort.bam
        fi

        if [[ ! -s ${lib[$i]}/results/${lib[$i]}.${Genome}.CBUB.sorted.bam ]]; then
            
        #   featureCounts -T $threads -Q 10 -M -a $gtf -t exon -g gene_name \
        #       -o ${lib[$i]}.${Genome}.exon.count.txt -R SAM ${lib[$i]}.${Genome}.namesort.CBUB.bam >>Run.log
            featureCounts -T $threads -Q 10 -M -a $gtf -t gene -g gene_name \
                -o ${lib[$i]}.${Genome}.gene.count.txt -R SAM ${lib[$i]}.${Genome}.namesort.CBUB.bam >>Run.log
            
            # calculate the stat of umi count for each cell barcodes and filter barcodes. umi_thr>=500
            ${toolPATH}/stats_cb_gene_umi ${lib[$i]}.${Genome}.namesort.CBUB.bam.featureCounts.sam ${lib[$i]}.${Genome}.namesort.CBUB.bam.cb_gc_uc.stat.txt > /dev/null 2>&1
            Rscript ${toolPATH}/filter_barcodes.r ${lib[$i]}.${Genome}.namesort.CBUB.bam.cb_gc_uc.stat.txt ${lib[$i]}.${Genome}.filtered_bc.txt ${lib[$i]}.${Genome}.cb_gc_uc.summary.txt
            
            samtools view ${lib[$i]}.${Genome}.namesort.CBUB.bam.featureCounts.sam -b -o ${lib[$i]}.${Genome}.namesort.CBUB.bam.featureCounts.bam
            samtools sort -@ $threads ${lib[$i]}.${Genome}.namesort.CBUB.bam.featureCounts.bam -o ${lib[$i]}.${Genome}.CBUB.sorted.bam
            samtools index ${lib[$i]}.${Genome}.CBUB.sorted.bam
            samtools view ${lib[$i]}.${Genome}.CBUB.sorted.bam > ${lib[$i]}.${Genome}.CBUB.sorted.sam
            rm ${lib[$i]}.${Genome}.namesort.CBUB.bam.featureCounts.sam ${lib[$i]}.${Genome}.namesort.CBUB.bam.featureCounts.bam
        fi

        if [[ ! -s ${lib[$i]}/results/${lib[$i]}.${Genome}.velocyte/${lib[$i]}.${Genome}.loom ]]; then
            # Run velocyte
            velocyto run -@ $threads \
                -b ${lib[$i]}.${Genome}.filtered_bc.txt \
                -o ${lib[$i]}.${Genome}.velocyte \
                -e ${lib[$i]}.${Genome} \
                ${lib[$i]}.${Genome}.CBUB.sorted.bam \
                $gtf
            
            # Get umi count from velocyte run
            python3 ${toolPATH}/get_umic_from_velocyteloom.py ${lib[$i]}.${Genome}.velocyte/${lib[$i]}.${Genome}.loom ${lib[$i]}.${Genome}.velocyte/${lib[$i]}.${Genome}.cellID_umic.txt

            # statistics
            cd ../
            wc -l trimmed/${lib[$i]}_1.fq.gz > ${lib[$i]}_1.fq.wc.txt 
            bamToBed -i results/${lib[$i]}.${Genome}.CBUB.sorted.bam  > ${lib[$i]}.${Genome}.bed
            ${toolPATH}/RemoveRedudantReads ${lib[$i]}.${Genome}.bed ${lib[$i]}.${Genome}.noDup.bed
            wc -l ${lib[$i]}.${Genome}.bed ${lib[$i]}.${Genome}.noDup.bed > ${lib[$i]}.${Genome}.wc_bed.txt 
            cd ../
        fi
    fi

    
    
    if [[ ${lib[$i]} =~ ^sa.* ]] || [[ ${lib[$i]} =~ ^sc.* ]]; then
        
        if [ -s ${lib[$i]}/ATAC.${Genome}.txt ]; then
            echo "${lib[$i]}/ATAC.${Genome}.txt exists"
        else
            for (( j=0; j< ${#ATAC[@]}; j++ )); do
                echo ${ATAC[$j]} >> ${lib[$i]}/ATAC.${Genome}.txt
            done
        fi
        if [ -s ${lib[$i]}/CUTTag.${Genome}.txt ]; then
            echo "${lib[$i]}/CUTTag.${Genome}.txt exists"
        else
            for (( j=0; j< ${#CUTTag[@]}; j++ )); do
                echo ${CUTTag[$j]} >> ${lib[$i]}/CUTTag.${Genome}.txt
            done
        fi
    
        if [[ $Genome == mm10 ]]; then
            export BOWTIE2_INDEXES=/data02/laibb/annotation/genome/refdata-gex-mm10-2020-A/bowtie2_index
        elif [[ $Genome == mm39 ]]; then
            export BOWTIE2_INDEXES=/data02/laibb/annotation/genome/refdata-gex-GRCm39-2024-A/bowtie2_index
        elif [[ $Genome == hg38 ]]; then
            export BOWTIE2_INDEXES=/data02/laibb/annotation/genome/refdata-gex-GRCh38-2024-A/bowtie2_index
        else
            echo "genome not found"
            exit
        fi

        # Align
        echo "map ${lib[$i]}-ATAC(R2) to ${Genome} genome"
        for (( j=0; j< ${#ATAC[@]}; j++ )); do
            
            if [[ ! -s ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.${Genome}.sorted.bam ]]; then
                echo "map ${ATAC[$j]}"
                bowtie2 -p $threads -X 5000 -x ${Genome} --rg-id ${lib[$i]}-${ATAC[$j]} \
                     -1 ${lib[$i]}/trimmed/${lib[$i]}-${ATAC[$j]}_1.fq.gz -2 ${lib[$i]}/trimmed/${lib[$i]}-${ATAC[$j]}_2.fq.gz \
                     2> ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.${Genome}.bowtie2.log | \
                    samtools view -O bam -o ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.${Genome}.bam
                samtools sort -@ $threads ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.${Genome}.bam \
                    -o ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.${Genome}.sorted.bam
                samtools index ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.${Genome}.sorted.bam
            fi
        done
    
        echo "map ${lib[$i]}-CUTTag(R2) to ${Genome} genome"
        for (( j=0; j< ${#CUTTag[@]}; j++ )); do
            if [[ ! -s ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.${Genome}.sorted.bam ]]; then
                echo "map ${CUTTag[$j]}"
                bowtie2 -p $threads -X 5000 -x ${Genome} --rg-id ${lib[$i]}-${CUTTag[$j]} \
                        -1 ${lib[$i]}/trimmed/${lib[$i]}-${CUTTag[$j]}_1.fq.gz -2 ${lib[$i]}/trimmed/${lib[$i]}-${CUTTag[$j]}_2.fq.gz \
                        2> ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.${Genome}.bowtie2.log | \
                        samtools view -O bam -o ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.${Genome}.bam
                samtools sort -@ $threads ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.${Genome}.bam \
                    -o ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.${Genome}.sorted.bam
                samtools index ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.${Genome}.sorted.bam
            fi
        done

        echo "filter chr, sort by name, separate bedpe, R2 are ATAC beds, R1 depends on sc* or sa* lib"
        
        for (( j=0; j< ${#ATAC[@]}; j++ )); do
            echo ${ATAC[$j]}
            if [[ -s "${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R2.${Genome}.noDup.bed" ]]; then
                continue
            fi
            chrs=`samtools view -H ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.${Genome}.sorted.bam | \
                    grep chr | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<6)print}'`
            samtools view -@ $threads -b -q 10 -f 2 ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.${Genome}.sorted.bam `echo $chrs` | \
                    samtools sort -@ $threads -m 4G -n -o ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.${Genome}.namesort.bam
            
            bamToBed -bedpe -mate1 -i ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.${Genome}.namesort.bam \
                    > ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.${Genome}.namesort.bedpe
            cat ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.${Genome}.namesort.bedpe | awk -v OFS="\t" '{print $1,$2+4,$3-5,$7,$8,$9}' \
                    > ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R1.${Genome}.bed
            ${toolPATH}/RemoveRedudantReads ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R1.${Genome}.bed \
                    ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R1.${Genome}.noDup.bed
            cat ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.${Genome}.namesort.bedpe | awk -v OFS="\t" '{print $4,$5+4,$6-5,$7,$8,$10}' \
                    > ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R2.${Genome}.bed
            ${toolPATH}/RemoveRedudantReads ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R2.${Genome}.bed \
                    ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R2.${Genome}.noDup.bed

            generatewig_from_bed ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R2.${Genome}.noDup.bed 1 50 ${lib[$i]}-${ATAC[$j]} ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R2.${Genome}.noDup.wig  
            wigToBigWig -clip ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R2.${Genome}.noDup.wig $ChromSizes ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R2.${Genome}.noDup.bw
            generatewig_from_bed ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R1.${Genome}.noDup.bed 1 50 ${lib[$i]}-${ATAC[$j]} ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R1.${Genome}.noDup.wig  
            wigToBigWig -clip ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R1.${Genome}.noDup.wig $ChromSizes ${lib[$i]}/results/${lib[$i]}-${ATAC[$j]}.R1.${Genome}.noDup.bw
        done
        echo "filter chr, sort by name, separate bedpe, R2 are CUTTag beds, R1 depends on sc* or sa* lib"
        for (( j=0; j< ${#CUTTag[@]}; j++ )); do
            echo ${CUTTag[$j]}
            if [[ -s "${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R2.${Genome}.noDup.bed" ]]; then
                continue
            fi
            chrs=`samtools view -H ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.${Genome}.sorted.bam | \
                    grep chr  | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<6)print}'`
            samtools view -@ $threads -b -q 10 -f 2 ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.${Genome}.sorted.bam `echo $chrs` | \
                    samtools sort -@ $threads -m 4G -n -o ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.${Genome}.namesort.bam
            
            bamToBed -bedpe -mate1 -i ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.${Genome}.namesort.bam \
                    > ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.${Genome}.namesort.bedpe
            cat ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.${Genome}.namesort.bedpe | awk -v OFS="\t" '{print $1,$2+4,$3-5,$7,$8,$9}' \
                    > ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R1.${Genome}.bed
            ${toolPATH}/RemoveRedudantReads ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R1.${Genome}.bed \
                    ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R1.${Genome}.noDup.bed
            cat ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.${Genome}.namesort.bedpe | awk -v OFS="\t" '{print $4,$5+4,$6-5,$7,$8,$10}' \
                    > ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R2.${Genome}.bed
            ${toolPATH}/RemoveRedudantReads ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R2.${Genome}.bed \
                    ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R2.${Genome}.noDup.bed
            generatewig_from_bed ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R2.${Genome}.noDup.bed 1 50 ${lib[$i]}-${CUTTag[$j]} ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R2.${Genome}.noDup.wig    
            wigToBigWig -clip ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R2.${Genome}.noDup.wig $ChromSizes ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R2.${Genome}.noDup.bw
            generatewig_from_bed ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R1.${Genome}.noDup.bed 150 200 ${lib[$i]}-${CUTTag[$j]} ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R1.${Genome}.noDup.wig    
            wigToBigWig -clip ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R1.${Genome}.noDup.wig $ChromSizes ${lib[$i]}/results/${lib[$i]}-${CUTTag[$j]}.R1.${Genome}.noDup.bw
    
        done

    #    exit

    #    echo "MACS2 call peaks"
    #    parallel -j $threads 'macs2 callpeak -t '${lib[$i]}'/results/'${lib[$i]}'-{}.R2.'${Genome}'.noDup.bed \
     #       --outdir '${lib[$i]}'/results/ -n '${lib[$i]}'-{}.R2.'${Genome}'.noDup --nomodel --shift 100 --extsize 200 -g mm' :::: sa02/ATAC.${Genome}.txt
   
    fi

done    


# combine ATAC and CUTTAG beds
declare -A sa_items
declare -A sc_items
declare -A sr_items
declare -A sublibs

# 
for item in "${lib[@]}"; do
    # 
    prefix="${item:0:2}"  # 
    suffix="${item:2}"    # 
    
    # 
    case "$prefix" in
        "sa")
            sa_items["$suffix"]="$item"
            ;;
        "sc")
            sc_items["$suffix"]="$item"
            ;;
        "sr")
            sr_items["$suffix"]="$item"
            ;;
    esac
    
    sublibs["$suffix"]=1
done

tagbedfiles="combined.bedfiles.${Genome}.txt"
bcsamfile=${Genome}.bc_samfile
> $tagbedfiles
> $bcsamfile

for suffix in $(echo "${!sublibs[@]}" | tr ' ' '\n' | sort); do
    echo "sublib $suffix:"

    sa_lib=""
    sc_lib=""
    sr_lib=""
    if [ -n "${sa_items[$suffix]}" ]; then
        echo "  sa: ${sa_items[$suffix]}"
        sa_lib=${sa_items[$suffix]}
    fi
    
    if [ -n "${sc_items[$suffix]}" ]; then
        echo "  sc: ${sc_items[$suffix]}"
        sc_lib=${sc_items[$suffix]}
    fi
    
    if [ -n "${sr_items[$suffix]}" ]; then
        echo "  sr: ${sr_items[$suffix]}"
        sr_lib=${sr_items[$suffix]}
    fi
    

    if [[ -n "$sa_lib" && -n "$sc_lib" ]]; then 
        echo "Assign barcode sample tag for $sa_lib and $sc_lib"
        if [[ ! -d ${sa_lib}/results ]]; then
            echo "cannot find dir ${sa_lib}/results"
            exit
        fi
        if [[ ! -d ${sc_lib}/results ]]; then
            echo "cannot find dir ${sc_lib}/results"
            exit
        fi

        if [[ -s ${sa_lib}_${sc_lib}.${Genome}.bedfile.txt ]]; then
            echo "The file ${sa_lib}_${sc_lib}.bedfile.txt exists."
        else
            for (( j=0; j< ${#ATAC[@]}; j++ )); do
                if [[ ! -s ${sa_lib}/results/${sa_lib}-${ATAC[$j]}.R2.${Genome}.noDup.bed ]]; then
                    echo "the file ${sa_lib}/results/${sa_lib}-${ATAC[$j]}.R2.${Genome}.noDup.bed not exist"
                    exit
                fi
                printf "sa_sample_R2\t%s\t%s\n" \
                    "${ATAC[$j]}" \
                    "${sa_lib}/results/${sa_lib}-${ATAC[$j]}.R2.${Genome}.noDup.bed" \
                    >> "${sa_lib}_${sc_lib}.${Genome}.bedfile.txt"
            #   echo "sa_sample_R2\t${ATAC[$j]}\t${sa_lib}/results/${sa_lib}-${ATAC[$j]}.R2.${Genome}.noDup.bed" >> ${sa_lib}_${sc_lib}.bedfile.txt
                if [[ ! -s ${sa_lib}/results/${sa_lib}-${ATAC[$j]}.R1.${Genome}.noDup.bed ]]; then
                    echo "the file ${sa_lib}/results/${sa_lib}-${ATAC[$j]}.R1.${Genome}.noDup.bed not exist"
                    exit
                fi
                printf "sa_sample_R1\t%s\t%s\n" \
                    "${ATAC[$j]}" \
                    "${sa_lib}/results/${sa_lib}-${ATAC[$j]}.R1.${Genome}.noDup.bed" \
                    >> "${sa_lib}_${sc_lib}.${Genome}.bedfile.txt"
            #   echo "sa_sample_R1\t${ATAC[$j]}\t${sa_lib}/results/${sa_lib}-${ATAC[$j]}.R1.${Genome}.noDup.bed" >> ${sa_lib}_${sc_lib}.bedfile.txt
                if [[ ! -s ${sc_lib}/results/${sc_lib}-${ATAC[$j]}.R2.${Genome}.noDup.bed ]]; then
                    echo "the file ${sc_lib}/results/${sc_lib}-${ATAC[$j]}.R2.${Genome}.noDup.bed not exist"
                    exit
                fi
                printf "sc_sample_R2\t%s\t%s\n" \
                    "${ATAC[$j]}" \
                    "${sc_lib}/results/${sc_lib}-${ATAC[$j]}.R2.${Genome}.noDup.bed" \
                    >> "${sa_lib}_${sc_lib}.${Genome}.bedfile.txt"
            #    echo "sc_sample_R2\t${ATAC[$j]}\t${sc_lib}/results/${sc_lib}-${ATAC[$j]}.R2.${Genome}.noDup.bed" >> ${sa_lib}_${sc_lib}.bedfile.txt
                if [[ ! -s ${sc_lib}/results/${sc_lib}-${ATAC[$j]}.R1.${Genome}.noDup.bed ]]; then
                    echo "the file ${sc_lib}/results/${sc_lib}-${ATAC[$j]}.R1.${Genome}.noDup.bed not exist"
                    exit
                fi
                printf "sc_sample_R1\t%s\t%s\n" \
                    "${ATAC[$j]}" \
                    "${sc_lib}/results/${sc_lib}-${ATAC[$j]}.R1.${Genome}.noDup.bed" \
                    >> "${sa_lib}_${sc_lib}.${Genome}.bedfile.txt"
            #    echo "sc_sample_R1\t${ATAC[$j]}\t${sc_lib}/results/${sc_lib}-${ATAC[$j]}.R1.${Genome}.noDup.bed" >> ${sa_lib}_${sc_lib}.bedfile.txt
            
            done


            for (( j=0; j< ${#CUTTag[@]}; j++ )); do
                if [[ ! -s ${sa_lib}/results/${sa_lib}-${CUTTag[$j]}.R2.${Genome}.noDup.bed ]]; then
                    echo "the file ${sa_lib}/results/${sa_lib}-${CUTTag[$j]}.R2.${Genome}.noDup.bed not exist"
                    exit
                fi
                printf "sa_tag_R2\t%s\t%s\n" \
                    "${CUTTag[$j]}" \
                    "${sa_lib}/results/${sa_lib}-${CUTTag[$j]}.R2.${Genome}.noDup.bed" \
                    >> "${sa_lib}_${sc_lib}.${Genome}.bedfile.txt"
            #    echo "sa_sample_R2\t${CUTTag[$j]}\t${sa_lib}/results/${sa_lib}-${CUTTag[$j]}.R2.${Genome}.noDup.bed" >> ${sa_lib}_${sc_lib}.bedfile.txt
                if [[ ! -s ${sa_lib}/results/${sa_lib}-${CUTTag[$j]}.R1.${Genome}.noDup.bed ]]; then
                    echo "the file ${sa_lib}/results/${sa_lib}-${CUTTag[$j]}.R1.${Genome}.noDup.bed not exist"
                    exit
                fi
                printf "sa_tag_R1\t%s\t%s\n" \
                    "${CUTTag[$j]}" \
                    "${sa_lib}/results/${sa_lib}-${CUTTag[$j]}.R1.${Genome}.noDup.bed" \
                    >> "${sa_lib}_${sc_lib}.${Genome}.bedfile.txt"
            #    echo "sa_sample_R1\t${CUTTag[$j]}\t${sa_lib}/results/${sa_lib}-${CUTTag[$j]}.R1.${Genome}.noDup.bed" >> ${sa_lib}_${sc_lib}.bedfile.txt
                if [[ ! -s ${sc_lib}/results/${sc_lib}-${CUTTag[$j]}.R2.${Genome}.noDup.bed ]]; then
                    echo "the file ${sc_lib}/results/${sc_lib}-${CUTTag[$j]}.R2.${Genome}.noDup.bed not exist"
                    exit
                fi
                printf "sc_tag_R2\t%s\t%s\n" \
                    "${CUTTag[$j]}" \
                    "${sc_lib}/results/${sc_lib}-${CUTTag[$j]}.R2.${Genome}.noDup.bed" \
                    >> "${sa_lib}_${sc_lib}.${Genome}.bedfile.txt"
            #   echo "sc_sample_R2\t${CUTTag[$j]}\t${sc_lib}/results/${sc_lib}-${CUTTag[$j]}.R2.${Genome}.noDup.bed" >> ${sa_lib}_${sc_lib}.bedfile.txt
                if [[ ! -s ${sc_lib}/results/${sc_lib}-${CUTTag[$j]}.R1.${Genome}.noDup.bed ]]; then
                    echo "the file ${sc_lib}/results/${sc_lib}-${CUTTag[$j]}.R1.${Genome}.noDup.bed not exist"
                    exit
                fi
                printf "sc_tag_R1\t%s\t%s\n" \
                    "${CUTTag[$j]}" \
                    "${sc_lib}/results/${sc_lib}-${CUTTag[$j]}.R1.${Genome}.noDup.bed" \
                    >> "${sa_lib}_${sc_lib}.${Genome}.bedfile.txt"
            #   echo "sc_sample_R1\t${CUTTag[$j]}\t${sc_lib}/results/${sc_lib}-${ATAC[$j]}.R1.${Genome}.noDup.bed" >> ${sa_lib}_${sc_lib}.bedfile.txt
            
            done
        fi

        echo "Assign_barcode_sample_tag"
        if [[ ! -d ${sa_lib}_${sc_lib}.${Genome} ]]; then
            mkdir ${sa_lib}_${sc_lib}.${Genome}
        fi

        if [[ ! -s "${sa_lib}_${sc_lib}.${Genome}/sample_barcodes_mapping.txt" ]]; then
            ${toolPATH}/Assign_barcode_sample_tag -i ${sa_lib}_${sc_lib}.${Genome}.bedfile.txt -o ${sa_lib}_${sc_lib}.${Genome}/
        fi

        # Align DNA_RNA
        if [[ -n "$sr_lib" ]]; then

            echo "Align dna and rna"
            if [[ -s ${sr_lib}_${sa_lib}_${sc_lib}.${Genome}.aligned_dna_rna_bc.txt ]]; then
                echo "Align file already done"
            else

                if [[ -s ${sa_lib}_${sc_lib}.${Genome}.integrated.bedfile.txt ]]; then
                    echo "The file ${sa_lib}_${sc_lib}.${Genome}.integrated.bedfile.txt exists."
                    
                else
                    for (( j=0; j< ${#ATAC[@]}; j++ )); do
                        if [[ -s ${sa_lib}_${sc_lib}.${Genome}/${ATAC[$j]}.atac.bed ]]; then
                            printf "%s\t%s\t%s/%s.atac.bed\n" \
                                ${ATAC[$j]} \
                                "atac" \
                                "${sa_lib}_${sc_lib}.${Genome}" \
                                ${ATAC[$j]} \
                                >> "${sa_lib}_${sc_lib}.${Genome}.integrated.bedfile.txt"
                        fi

                        for (( k=0; k< ${#CUTTag[@]}; k++ )); do
                            if [[ -s ${sa_lib}_${sc_lib}.${Genome}/${ATAC[$j]}.${CUTTag[$k]}.bed ]]; then
                                printf "%s\t%s\t%s/%s.%s.bed\n" \
                                    ${ATAC[$j]} \
                                    "${CUTTag[$k]}" \
                                    "${sa_lib}_${sc_lib}.${Genome}" \
                                    "${ATAC[$j]}" \
                                    "${CUTTag[$k]}" \
                                    >> "${sa_lib}_${sc_lib}.${Genome}.integrated.bedfile.txt"
                            fi
                        done
                    done
                fi

                ${toolPATH}/Align_rna_dna -d ${sa_lib}_${sc_lib}.${Genome}.integrated.bedfile.txt \
                    -r ${sr_lib}/results/${sr_lib}.${Genome}.namesort.CBUB.bam.cb_gc_uc.stat.txt \
                    -R ${sr_lib}/results/${sr_lib}.${Genome}.velocyte/${sr_lib}.${Genome}.cellID_umic.txt \
                    -o ${sr_lib}_${sa_lib}_${sc_lib}.${Genome}.aligned_dna_rna_bc.txt
                    
            fi

            echo "${sa_lib}_${sc_lib}.${Genome}.integrated.bedfile.txt" >> $tagbedfiles
        fi



    elif [[ -n "$sa_lib" && -z "$sc_lib" ]]; then
        echo "only sa $sa_lib exists, no sc lib"

        echo "Assign barcode sample tag for $sa_lib "
        if [[ ! -d ${sa_lib}/results ]]; then
            echo "cannot find dir ${sa_lib}/results"
            exit
        fi

        if [[ -s ${sa_lib}.${Genome}.bedfile.txt ]]; then
            echo "The file ${sa_lib}.bedfile.txt exists."
        else
            for (( j=0; j< ${#ATAC[@]}; j++ )); do
                if [[ ! -s ${sa_lib}/results/${sa_lib}-${ATAC[$j]}.R2.${Genome}.noDup.bed ]]; then
                    echo "the file ${sa_lib}/results/${sa_lib}-${ATAC[$j]}.R2.${Genome}.noDup.bed not exist"
                    exit
                fi
                printf "sa_sample_R2\t%s\t%s\n" \
                    "${ATAC[$j]}" \
                    "${sa_lib}/results/${sa_lib}-${ATAC[$j]}.R2.${Genome}.noDup.bed" \
                    >> "${sa_lib}.${Genome}.bedfile.txt"
            #    echo "sa_sample_R2\t${ATAC[$j]}\t${sa_lib}/results/${sa_lib}-${ATAC[$j]}.R2.${Genome}.noDup.bed" >> ${sa_lib}_${sc_lib}.bedfile.txt
            done

            for (( j=0; j< ${#CUTTag[@]}; j++ )); do
                if [[ ! -s ${sa_lib}/results/${sa_lib}-${CUTTag[$j]}.R2.${Genome}.noDup.bed ]]; then
                    echo "the file ${sa_lib}/results/${sa_lib}-${CUTTag[$j]}.R2.${Genome}.noDup.bed not exist"
                    exit
                fi
                printf "sa_tag_R2\t%s\t%s\n" \
                    "${CUTTag[$j]}" \
                    "${sa_lib}/results/${sa_lib}-${CUTTag[$j]}.R2.${Genome}.noDup.bed" \
                    >> "${sa_lib}.${Genome}.bedfile.txt"
            done
        fi

        echo "Assign_barcode_sample_tag"
        if [[ ! -d ${sa_lib}.${Genome} ]]; then
            mkdir ${sa_lib}.${Genome}
        fi

        if [[ ! -s "${sa_lib}.${Genome}/sample_barcodes_mapping.txt" ]]; then
            ${toolPATH}/Assign_barcode_sample_tag -i ${sa_lib}.${Genome}.bedfile.txt -o ${sa_lib}.${Genome}/ -s 1
        fi


        if [[ -n "$sr_lib" ]]; then

            echo "Align dna and rna"
            if [[ -s ${sr_lib}_${sa_lib}.${Genome}.aligned_dna_rna_bc.txt ]]; then
                echo "Align file already done"
            else
                if [[ -s ${sa_lib}.${Genome}.integrated.bedfile.txt ]]; then
                    echo "The file ${sa_lib}.${Genome}.integrated.bedfile.txt exists."
                else
                    for (( j=0; j< ${#ATAC[@]}; j++ )); do
                        if [[ -s ${sa_lib}.${Genome}/${ATAC[$j]}.atac.bed ]]; then
                            printf "%s\t%s\t%s/%s.atac.bed\n" \
                                ${ATAC[$j]} \
                                "atac" \
                                "${sa_lib}.${Genome}" \
                                ${ATAC[$j]} \
                                >> "${sa_lib}.${Genome}.integrated.bedfile.txt"
                        fi

                        for (( k=0; k< ${#CUTTag[@]}; k++ )); do
                            if [[ -s ${sa_lib}.${Genome}/${ATAC[$j]}.${CUTTag[$k]}.bed ]]; then
                                printf "%s\t%s\t%s/%s.%s.bed\n" \
                                    ${ATAC[$j]} \
                                    "${CUTTag[$k]}" \
                                    "${sa_lib}.${Genome}" \
                                    "${ATAC[$j]}" \
                                    "${CUTTag[$k]}" \
                                    >> "${sa_lib}.${Genome}.integrated.bedfile.txt"
                            fi
                        done
                    done
                    
                fi

                ${toolPATH}/Align_rna_dna -d ${sa_lib}.${Genome}.integrated.bedfile.txt \
                    -r ${sr_lib}/results/${sr_lib}.${Genome}.namesort.CBUB.bam.cb_gc_uc.stat.txt \
                    -R ${sr_lib}/results/${sr_lib}.${Genome}.velocyte/${sr_lib}.${Genome}.cellID_umic.txt \
                    -o ${sr_lib}_${sa_lib}.${Genome}.aligned_dna_rna_bc.txt

            fi

            echo "${sa_lib}.${Genome}.integrated.bedfile.txt" >> $tagbedfiles
            echo "${sr_lib}_${sa_lib}.${Genome}.aligned_dna_rna_bc.txt ${sr_lib}/results/${sr_lib}.${Genome}.CBUB.sorted.sam" >> $bcsamfile
        fi
  

    elif [[ -z "$sa_lib" && -n "$sc_lib" ]]; then
        echo "only sc $sc_lib exist, no sa lib"

        
        
    else
        echo "no sa and sc libs"
    fi

done

# Merge all bedfiles for tag (including ATAC) to compute unit peaks 

bedfiles=${Genome}.tagbedfile

peakfiles=${Genome}.peakfile

if [[ -s $bedfile && -s $peakfile ]]; then
    echo "peakfile and bedfile exist"
else
    > $bedfiles
    > $peakfiles

    if [[ -s $tagbedfiles ]]; then
        if [[ ! -d combined.peaks.${Genome} ]]; then
            mkdir combined.peaks.${Genome}
        fi
        
        readarray -t file_list < "$tagbedfiles"
        if [[ ${#file_list[@]} -eq 0 ]]; then
            echo "error: $tagbedfiles has no file list"
            exit 1
        fi

        echo "found ${#file_list[@]} input files"

    
        all_tags=$(cat "${file_list[@]}" 2>/dev/null | awk -F'\t' '{print $2}' | sort -u)

        for tag in $all_tags; do
            merged_bed="combined.peaks.${Genome}/combined.${tag}.${Genome}.bed"
            # clear or setup empty outfile
            > $merged_bed
            echo "process tag $tag"
            # find tag-related bedfiles and merge and call peaks.
            for input_file in "${file_list[@]}"; do
                if [[ ! -f "$input_file" ]]; then
                    echo "   warning:  $input_file not exist, skip"
                    continue
                fi
                echo "$input_file"
                # substract
                awk -F'\t' -v t="$tag" '$2 == t {print $3}' "$input_file" | while read -r bed_file; do
                    if [[ -f "$bed_file" ]]; then
                        cat "$bed_file" >> "$merged_bed"
                        
                        echo "    add: $bed_file"
                    else
                        echo "    warning:  $bed_file not exist skip"
                    fi
                done
            done

            if [[ ! -s "$merged_bed" ]]; then
                echo "  warning: merged file empty, skip $tag"
                continue
            fi

            printf "%s\t%s\n" \
                    "$tag" \
                    "$merged_bed" \
                    >> "$bedfiles"

            # set sicer2 parameters
            if [[ "$tag" == "atac" ]]; then
                sicer_params="-w 100 -g 200"
                peak_suffix="-W100-G200"
            elif [[ "$tag" == *"K27me3"* ]] || [[ "$tag" == *"K9me3"* ]]; then
                sicer_params="-w 200 -g 600"
                peak_suffix="-W200-G600"
            else
                # default, applied to K4me1, K4me3, K27ac, K36me3, etc..
                sicer_params="-w 200 -g 200"
                peak_suffix="-W200-G200"
            fi

            peak_file="combined.peaks.${Genome}/combined.${tag}.${Genome}${peak_suffix}.scoreisland"
                
            if [[ ! -s "$peak_file" ]]; then
                echo "    run sicer2: $sicer_params"
                sicer -t "$merged_bed" -f 150 $sicer_params -s "$Genome" -o "combined.peaks.${Genome}/" > /dev/null 2>&1
            else
                echo "    sicer2 peaks exist"
            fi

            printf "%s\t%s\n" \
                    "$tag" \
                    "$peak_file" \
                    >> "$peakfiles"
            
            echo "  Genarate: $merged_bed"
            echo "  Genarate peak file: $peak_file"
            echo ""
        done

        echo "Merge bed and generate peaks done"
        echo "collect all merged bedfile list: $bedfiles"
        echo "collect all peakfile list: $peakfiles"
    else
        echo " warning :  $tagbedfiles empty"
        
    fi
fi

# Generate Mtx
if [[ -n "$bcsamfile" && -n "$bedfiles" && -n "$peakfiles" ]]; then  
    # Integrate multiomics and generate matrix.
    echo "Integrate multiomics and generate matrix"
    if [[ ! -d "mtx.${Genome}.out/" ]]; then
        mkdir mtx.${Genome}.out/
    fi
    ${toolPATH}/integrate_and_generate_mtx $bcsamfile $peakfiles $bedfiles mtx.${Genome}.out/ $threads 

   
fi

echo "Complete!"

