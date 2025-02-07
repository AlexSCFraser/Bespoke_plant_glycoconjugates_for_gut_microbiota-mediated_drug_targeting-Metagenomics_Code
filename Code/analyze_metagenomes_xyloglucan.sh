#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=502G
#SBATCH --time=0-12:00           # time (DD-HH:MM)

if [ ${listname+x} ]; then
    echo "Loading accessions from $listname via variable"
elif [[ $# -gt 0 ]]; then
    listname=$1
    echo "Loading accessions from $listname via argument"
else
    echo ERROR: No argument given! Please give a file with accessions as input.
    exit
fi

if [ ${alignment_perc_identity+x} ]; then
    echo "Loaded alignment_perc_identity=$alignment_perc_identity via variable"
elif [[ $# -gt 1 ]]; then
    alignment_perc_identity=$2
    echo "Loaded alignment_perc_identity=$alignment_perc_identity via argument"
else
    alignment_perc_identity="95.0"
    echo "No alignment percent identity variable given, defaulting to 95.0"
fi

if [ ${out_dir+x} ]; then
    echo "Loaded out_dir=$out_dir via variable"
elif [[ $# -gt 2 ]]; then
    out_dir=$3
    echo "Loaded out_dir=$out_dir via argument"
else
    out_dir="output"
    echo "No out_dir variable given, defaulting to 'output'"
fi

# jobs * threads should be less than or equal to available cpu cores on local machine
# Local machine settings
jobs=2
threads=2

# sets "unofficial bash strict mode" http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

if [ ${SLURM_TMPDIR+x} ]; then
    # load modules only on cluster, check for $SLURM_TMPDIR to verify if on cluster or local
    # Alliance Canada Cluster settings
    module load StdEnv/2020
    module load gcc/9.3.0
    module load sra-toolkit/3.0.0
    module load samtools/1.17
    module load bedtools/2.30.0
    module load blast+/2.12.0

    # jobs * threads should be equal to available cpu cores on cluster

    # NOTE: using 32 jobs requires over a terabyte of local hard drive space and is not recommended with the current alliance canada cluster settings, stick to 16 jobs
    # Use more jobs for smaller accession files which need less memory
    jobs=16
    threads=2

    # Use medium # of jobs for slightly larger accession files which need less memory
    # jobs=8
    # threads=4

    # Use fewer jobs for yet larger accession files which need more memory
    # jobs=4
    # threads=8
fi

file=(`cat "$listname" | sed 's/\r//g'`)
if [ ${SLURM_ARRAY_TASK_ID+x} ];then
    accessions=(${file[@]:SLURM_ARRAY_TASK_ID:100})
else
    accessions=(${file[@]})
fi

out_folder="${PWD}/${out_dir}/$alignment_perc_identity"
if [[ ! -d $out_folder ]]; then
    mkdir -p "$out_folder"
fi

if [ ! ${SLURM_TMPDIR+x} ]; then
    SLURM_TMPDIR=${out_folder}/temp #sets SLURM_TMPDIR to $out_folder/temp if unset
    if [[ ! -d $SLURM_TMPDIR ]]; then
        mkdir -p "$SLURM_TMPDIR"
    fi
    echo "setting temp directory to $SLURM_TMPDIR"
fi

if [ ! ${fasta_array+x} ]; then
  fasta_array=("7_xyguls.fasta")
else
  echo "Loading array of fasta databases from fasta_array variable..."
fi

db_array=("${fasta_array[@]/%.f*a/}")

for db_fasta in "${fasta_array[@]}"; do
    if [[ $db_fasta =~ [[:space:]]+ ]]; then
        mv "$db_fasta" ${db_fasta//\ /_}
        db_fasta=${db_fasta//\ /_}
    fi

    echo "$db_fasta"
    if [[ ! -f ${db_fasta/%f*a/ndb} ]]; then
        echo "Making database from fasta file: ${db_fasta}"
        makeblastdb -in "$db_fasta" -dbtype nucl -parse_seqids -out "${db_fasta/%.f*a}"
    fi
done

calculate_coverage () {
    # sets "unofficial bash strict mode" http://redsymbol.net/articles/unofficial-bash-strict-mode/
    # We set strict mode again inside this function because the earlier strict mode does not apply inside functions called by GNU parallel
    # Strict mode prevents the pipeline from continuing when any utility in the chain fails
    set -euo pipefail
    IFS=$'\n\t'

    local sra_num=$1
    local db_name=$2
    local threads=$3
    local out_folder=$4
    local percent_identity=$5


    echo "Calculate coverage, sra_num: $sra_num  db_name: $db_name  threads: $threads  out_folder: $out_folder percent_identity: $percent_identity"

    if [ ${SLURM_TMPDIR+x} ]; then
        # set higher memory usage on cluster
        mem=15
    else
        # set lower memory usage locally
        mem=1
        #sets SLURM_TMPDIR to $out_folder/temp if unset, which only happens locally
        SLURM_TMPDIR=${out_folder}/temp
    fi

    # when sra_num has a semicolon, treat it as a pair of URLs to download fastq files from
    if [[ $sra_num == *";"* ]]; then
      IFS=';' read -ra url_arrary <<< "$sra_num"
      url1="${url_arrary[0]}"
      url2="${url_arrary[1]}"
      filename="$(basename "$url1")"
      file1="${SLURM_TMPDIR}/${filename}"
      file2="${SLURM_TMPDIR}/${filename}"
      sra_num="${filename%_R.*}"
    fi

    local outbamname_filtered="${sra_num}-${db_name}.filtered_sorted.bam"
    local outbamindex="$outbamname_filtered.bai"
    local outfilename="${sra_num}-${db_name}.out"
    local outflagstatsname="${sra_num}-${db_name}.out.flagstats"
    local outbasestatsname="${sra_num}-${db_name}.out.basestats"
    local samfilename="${sra_num}-${db_name}.sam"

    start="$(date +%s)"
    magic_blast_done=$start
    samtools_done=$start
    if [[ ! -e ${out_folder}/${outfilename} ]] || [[ ! -e ${out_folder}/${outflagstatsname} ]] || [[ ! -e ${out_folder}/${outbasestatsname} ]] || [[ ! -e  ${out_folder}/$outbamname_filtered ]]; then

        if [ ${file1+x} ] && [ ${file2+x} ] && [[ ! -e $file1  || ! -e $file2 ]]; then
          wget "$url1" -O "$file1" --no-check-certificate --quiet
          wget "$url2" -O "$file2" --no-check-certificate --quiet
          echo "Downloaded paired files for identifier: ${sra_num}"
        fi

        if [[ ! -e ${SLURM_TMPDIR}/${samfilename} ]]; then

            if [ ${file1+x} ]; then
              inputarg=(-paired -query "${file1}" -query_mate "${file2}" -infmt fastq)
              echo "Running magicblast on $file1 paired with $file2"
              echo "${inputarg[@]}"
            else
              inputarg=(-sra "${sra_num}")
              echo "Running magicblast on $sra_num"
            fi
            if magicblast -db "$db_name" -num_threads "$threads" "${inputarg[@]}" -perc_identity "$percent_identity" > "${SLURM_TMPDIR}/${samfilename}";
            then
              echo "magicblast complete"
            else
              echo "magicblast failed"
              return
            fi

#            if [[ $? == 0 ]]; then
#                echo "magicblast complete"
#            else
#                echo "magicblast failed"
#                return
#            fi
            
        else
            echo "${samfilename} present, skipping magicblast"
        fi

        magic_blast_done="$(date +%s)"

        if [[ -e ${SLURM_TMPDIR}/${samfilename} ]] && [[ ! -e ${SLURM_TMPDIR}/$outbamname_filtered ]]; then
            echo "converting ${samfilename} to ${sra_num}_${db_name}.bam"
            samtools flagstat "${SLURM_TMPDIR}/${samfilename}" > "${out_folder}/${outflagstatsname}"
            samtools stats "${SLURM_TMPDIR}/${samfilename}" | grep '^SN' > "${out_folder}/${outbasestatsname}"
            samtools view -b "${SLURM_TMPDIR}/${samfilename}" -F 4 -@ "$threads" | samtools sort -o "${SLURM_TMPDIR}/$outbamname_filtered" -m ${mem}G -@ "$threads" -T "$SLURM_TMPDIR"
            samtools index "${SLURM_TMPDIR}/$outbamname_filtered" -@ "$threads"
            if [[ -e ${out_folder}/${outflagstatsname} ]]; then
                echo "${out_folder}/${outflagstatsname} written"
            else
                echo "${out_folder}/${outflagstatsname} FAILED TO WRITE"
            fi
            echo "sam to bam conversion complete"
        else
            echo "${outbamname_filtered} present, skipping conversion"
        fi

        samtools_done="$(date +%s)"


        if [[ -e ${SLURM_TMPDIR}/$outbamname_filtered ]]; then
            bedtools genomecov -ibam "${SLURM_TMPDIR}/$outbamname_filtered" > "${out_folder}/${outfilename}"
            sleep 1 # file said to not exist even though it was written fine todo: fix this file check
            if [[ -e ${out_folder}/${outfilename} ]]; then
                echo "${out_folder}/${outfilename} written"
            else
                echo "${out_folder}/${outfilename} FAILED TO WRITE"
            fi
            cp "${SLURM_TMPDIR}/$outbamname_filtered" "${out_folder}/$outbamname_filtered"
            sleep 1 # file said to not exist even though it was written fine
            if [[ -e ${out_folder}/$outbamname_filtered ]]; then
                echo "${out_folder}/$outbamname_filtered written"
            else
                echo "${out_folder}/$outbamname_filtered FAILED TO WRITE"
            fi
            cp "${SLURM_TMPDIR}/$outbamindex" "${out_folder}/$outbamindex"
            sleep 1 # file said to not exist even though it was written fine
            if [[ -e ${out_folder}/$outbamindex ]]; then
                echo "${out_folder}/$outbamindex written"
            else
                echo "${out_folder}/$outbamname_filtered FAILED TO WRITE"
            fi
        else
            echo "Could not find sorted bam file $outbamname_filtered"
        fi

#        rm "${SLURM_TMPDIR}/${sra_num}-${db_name}*"
        rm "$SLURM_TMPDIR/$outbamname_filtered"
        rm "$SLURM_TMPDIR/$outbamindex"
        rm "$SLURM_TMPDIR/$samfilename"
        if [ ${file1+x} ];then
          rm "$file1"
        fi
        if [ ${file2+x} ];then
          rm "$file2"
        fi
    else
        echo "${sra_num}_${db_name} output files already present in output folder!"
    fi
    all_done="$(date +%s)"
    duration=$(( all_done - start ))
    printf "magicblast took %02d minutes\n" $(( (magic_blast_done - start)/60 ))
    printf "samtools took %02d minutes\n" $(( (samtools_done - magic_blast_done)/60 ))
    printf "genomecov took %02d minutes\n" $(( (all_done - samtools_done)/60 ))

    printf "Coverage analysis of %s sample with %s database complete in %d:%02d:%02d\n" "$sra_num" "$db_name" $((duration/3600)) $(( (duration/60)%60 )) $((duration%60))
}


if [[ ! -d ${out_folder} ]]; then
    mkdir -p "$out_folder"
fi

export -f calculate_coverage
parallel -j${jobs} --joblog "${out_folder}/$(date +%m-%d-%y_%H-%M)"joblog.txt calculate_coverage {} "$threads" "$out_folder" $alignment_perc_identity ::: "${accessions[@]}" ::: "${db_array[@]}"
