#!/bin/bash

usage(){
echo "
Usage: assembly.sh [options]

Options:
  -i,  --input             Path to the input sequence directory.
  -o,  --outdir            Output directory.
  -a,  --assembly-options  Assembly options. (default: --megahit)
  -b,  --binning-options   Binning options. (default: --metabat2 --maxbin2 --concoct)
  -r,  --refine-options    Bin refine options. (default: -c 50 -x 10)
  -t,  --threads           Threads. (default: 1)
  -m,  --memory            memory in GB. (default: 24)
  --clean                  Clean mode. Remove all work files. (If your hard drive has less space)
  -h,  --help              Display this help message.
"
}

# default parameters
assembly_options="--megahit"
binning_options="--metabat2 --maxbin2 --concoct"
bin_refine_options="-c 50 -x 10"
threads=1
memory=24
clean=false

# parameters
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) seq_dir="$2"; shift ;;
        -o|--outdir) outdir="$2"; shift ;;
        -a|--assembly-options) assembly_options="$2"; shift ;;
        -b|--binning-options) binning_options="$2"; shift ;;
        -r|--refine-options) bin_refine_options="$2"; shift ;;
        -t|--threads) threads="$2"; shift ;;
        -m|--memory) memory="$2"; shift ;;
        --clean) clean=true ;;
        -h|--help) usage; exit 0;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done


# set root
script_dir=$(dirname "$(realpath "$0")")
ouput_folder="${outdir}/assembly_results"
mkdir $ouput_folder


# =====assembly=====
exec > >(tee -a "${ouput_folder}/output.log") 2>&1
echo "$(date '+%Y-%m-%d %H:%M:%S') - Assembly starts..."
echo "$(date '+%Y-%m-%d %H:%M:%S') - contigs processing..."

# set folder for tmp
tmp_dir=$ouput_folder/"tmp_files"
mkdir $tmp_dir

# set folder for contigs
contigs_dir=$ouput_folder/"contigs"
if [ -d $contigs_dir ]; then
  echo "$(date '+%Y-%m-%d %H:%M:%S') - Directory contigs already exists"
else
  mkdir $contigs_dir
fi

R1files=$(ls $seq_dir | awk '/_1/ {print $NF}' | sort)
R2files=$(ls $seq_dir | awk '/_2/ {print $NF}' | sort)
R1vec=()
for x in $R1files; do   R1vec+=($x); done
R2vec=()
for x in $R2files; do   R2vec+=($x); done
Len=${#R1vec[@]}

for (( i=0; i<${Len}; i++ )); do
  R1file=${R1vec[i]}
  R2file=${R2vec[i]}
  sample=${R2file%%_2*}
  outputrun=$contigs_dir/$sample

  if [ -d "$outputrun" ]; then
    if [ -f "${outputrun}/final_assembly.fasta" ]; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') - Assembling ${sample} is done"
    else
      echo "$(date '+%Y-%m-%d %H:%M:%S') - Assembling ${sample} is incomplete"
      rm -rf $outputrun
      metawrap assembly -1 "${seq_dir}/${R1file}" -2 "${seq_dir}/${R2file}" -o $contigs_dir/$sample -t $threads -m $memory $assembly_options
    fi
  else
    metawrap assembly -1 "${seq_dir}/${R1file}" -2 "${seq_dir}/${R2file}" -o $contigs_dir/$sample -t $threads -m $memory $assembly_options
  fi
done



# =====binning=====
echo "$(date '+%Y-%m-%d %H:%M:%S') - Bins processing..."

# set folder for binning
bins_dir=$ouput_folder/"bins"
if [ -d $bins_dir ]; then
  echo "$(date '+%Y-%m-%d %H:%M:%S') - Directory bins already exists"
else
  mkdir $bins_dir
fi

# main
for (( i=0; i<${Len}; i++ )); do
  R1file=${R1vec[i]}
  R2file=${R2vec[i]}
  sample=${R2file%%_2*}

  outputrun=$bins_dir/$sample
  if [ -d "$outputrun" ]; then
    if [ -f "${outputrun}/done" ]; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') - Binning ${sample} is done"
    else
      echo "$(date '+%Y-%m-%d %H:%M:%S') - Binning ${sample} is incomplete"
      contig="${contigs_dir}/${sample}/final_assembly.fasta"
      # resort reads
      $script_dir/bbmap/repair.sh in="${seq_dir}/${R1file}" in2="${seq_dir}/${R2file}" \
      out="${tmp_dir}/${sample}_1.fastq" out2="${tmp_dir}/${sample}_2.fastq"
      # metawrap
      metawrap binning -o $bins_dir/$sample -a $contig $binning_options -t $threads -m $memory "${tmp_dir}/${sample}_1.fastq" "${tmp_dir}/${sample}_2.fastq"
      rm "${tmp_dir}/${sample}_1.fastq" "${tmp_dir}/${sample}_2.fastq"
    fi

  else
    contig="${contigs_dir}/${sample}/final_assembly.fasta"
    # resort reads
    $script_dir/bbmap/repair.sh in="${seq_dir}/${R1file}" in2="${seq_dir}/${R2file}" \
    out="${tmp_dir}/${sample}_1.fastq" out2="${tmp_dir}/${sample}_2.fastq"
    # metawrap
    metawrap binning -o $bins_dir/$sample -a $contig $binning_options -t $threads -m $memory "${tmp_dir}/${sample}_1.fastq" "${tmp_dir}/${sample}_2.fastq"
    rm "${tmp_dir}/${sample}_1.fastq" "${tmp_dir}/${sample}_2.fastq"
  fi

  # done file
  num_options=$(echo $binning_options | wc -w)
  num_bins=$(ls -d $outputrun/*_bins 2>/dev/null | wc -l)
  # check bin folders are not empty
  non_empty_bins=0
  for dir in $outputrun/*_bins; do
    if [ "$(find "$dir" -type f | wc -l)" -gt 0 ]; then
      non_empty_bins=$((non_empty_bins + 1))
    fi
  done
  if [ $num_bins -ne $num_options ] || [ $non_empty_bins -ne $num_options ]; then
    touch $outputrun/"done"
  fi

  # clean mode
  if [ $clean = true ]; then
    if [ -f "${outputrun}/work_files" ]; then
      rm -rf "${outputrun}/work_files"
    fi
  fi
done
rm -r $tmp_dir


# =====bins refine=====
echo "$(date '+%Y-%m-%d %H:%M:%S') - Bins refinement processing..."

# set folder for binning
bins_refine_dir=$ouput_folder/"bins_refine"
if [ -d $bins_refine_dir ]; then
  echo "$(date '+%Y-%m-%d %H:%M:%S') - Directory bins_refine already exists"
else
  mkdir $bins_refine_dir
fi

samples=$(ls -l $bins_dir |awk '/^d/ {print $NF}')
for sample in $samples; do
  outputrun=$bins_refine_dir/$sample
  if [ -d "${outputrun}" ]; then
    if ls "${outputrun}/metawrap_*_*_bins" 1> /dev/null 2>&1; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') - Bins refinement ${sample} is done"
    else
      echo "$(date '+%Y-%m-%d %H:%M:%S') - Bins refinement ${sample} is incomplete"
      rm -rf $outputrun
      # enter into refine work
      all_bins_folder=($(find "${bins_dir}/${sample}" -type d -name "*_bins"))
      length=${#all_bins_folder[@]}

      # check number of bins and run bin_refinement
      if [ $length -eq 0 ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - No bins directory are found in ${sample}"
        exit 1
      elif [ $length -eq 1 ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Only 1 bin directory is found in ${sample}"
        exit 1
      elif [ $length -eq 2 ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - 2 bins directory are found in ${sample}"
        Adir=${all_bins_folder[0]}
        Bdir=${all_bins_folder[1]}
        metawrap bin_refinement -o $outputrun -t $threads -m $memory -A $Adir -B $Bdir $binning_options
      elif [ $length -eq 3 ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - 3 bins directory are found in ${sample}"
        Adir=${all_bins_folder[0]}
        Bdir=${all_bins_folder[1]}
        Cdir=${all_bins_folder[2]}
        metawrap bin_refinement -o $outputrun -t $threads -m $memory -A $Adir -B $Bdir -C $Cdir $bin_refine_options
      else
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Unexpected number of bins directory are found in ${sample}"
        exit 1
      fi
    fi

    # clean mode
    if [ $clean = true ]; then
      if [ -f "${outputrun}/work_files" ]; then
        rm -rf "${outputrun}/work_files"
      fi
      for dir in $outputrun/*_bins; do
        if [[ ! $(basename "$dir") =~ ^metawrap ]]; then
          rm -r "${dir}"
        fi
      done
    fi

  else
    all_bins_folder=($(find "${bins_dir}/${sample}" -type d -name "*_bins"))
    length=${#all_bins_folder[@]}

    # check number of bins and run bin_refinement
    if [ $length -eq 0 ]; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') - No bins directory are found in ${sample}"
      exit 1
    elif [ $length -eq 1 ]; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') - Only 1 bin directory is found in ${sample}"
      exit 1
    elif [ $length -eq 2 ]; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') - 2 bins directory are found in ${sample}"
      Adir=${all_bins_folder[0]}
      Bdir=${all_bins_folder[1]}
      metawrap bin_refinement -o $outputrun -t $threads -m $memory -A $Adir -B $Bdir $binning_options
    elif [ $length -eq 3 ]; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') - 3 bins directory are found in ${sample}"
      Adir=${all_bins_folder[0]}
      Bdir=${all_bins_folder[1]}
      Cdir=${all_bins_folder[2]}
      metawrap bin_refinement -o $outputrun -t $threads -m $memory -A $Adir -B $Bdir -C $Cdir $bin_refine_options
    else
      echo "$(date '+%Y-%m-%d %H:%M:%S') - Unexpected number of bins directory are found in ${sample}"
      exit 1
    fi

    # clean mode
    if [ $clean = true ]; then
      if [ -f "${outputrun}/work_files" ]; then
        rm -rf "${outputrun}/work_files"
      fi
      for dir in $outputrun/*_bins; do
        if [[ ! $(basename "$dir") =~ ^metawrap ]]; then
          rm -r "${dir}"
        fi
      done
    fi

  fi
done

echo "$(date '+%Y-%m-%d %H:%M:%S') - All done!"
