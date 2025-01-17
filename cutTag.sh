#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --job-name=CUTandTag_Analysis
#SBATCH --time=24:00:00
#SBATCH --error=job_name.err
#SBATCH --output=job_name.out
#SBATCH --mem=128G
#SBATCH --partition=compute_partition
#SBATCH --account=your_account_name

module load profile/bioinf
module load gnu/10.2.0--gcc--8.3.1
module load jre/1.8.0_111--binary
module load bowtie2
module load samtools
module load bedtools
module load picard

projPath="/path/to/project"
cores=8
ref="/path/to/reference/genome"
spikeInRef="/path/to/spike-in/reference"
chromSize="/path/to/chromosome/size/file"

sampleList="sample_names.txt"

while read histName; do
  echo "Processing sample: $histName"

  # Alignment
  bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 ${projPath}/data/raw/fastq/${histName}_R1.fastq.gz -2 ${projPath}/data/raw/fastq/${histName}_R2.fastq.gz -S ${projPath}/results/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/results/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt

  # Spike-in alignment and depth
  bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${spikeInRef} -1 ${projPath}/data/raw/fastq/${histName}_R1.fastq.gz -2 ${projPath}/data/raw/fastq/${histName}_R2.fastq.gz -S ${projPath}/results/alignment/sam/${histName}_bowtie2_spikeIn.sam &> ${projPath}/results/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.txt

  seqDepthDouble=$(samtools view -F 0x04 ${projPath}/results/alignment/sam/${histName}_bowtie2_spikeIn.sam | wc -l)
  seqDepth=$((seqDepthDouble/2))
  echo $seqDepth >${projPath}/results/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth

  # Duplicate removal
  picardCMD="java -jar /path/to/picard.jar"
  ${picardCMD} SortSam I=$projPath/results/alignment/sam/${histName}_bowtie2.sam O=$projPath/results/alignment/sam/${histName}_bowtie2.sorted.sam SORT_ORDER=coordinate
  ${picardCMD} MarkDuplicates I=$projPath/results/alignment/sam/${histName}_bowtie2.sorted.sam O=$projPath/results/alignment/removeDuplicate/${histName}_bowtie2.sorted.rmDup.sam REMOVE_DUPLICATES=true METRICS_FILE=$projPath/results/alignment/removeDuplicate/picard_summary/${histName}_picard.rmDup.txt

  # Spike-in calibration
  seqDepthFile="${projPath}/results/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth"
  seqDepth=$(cat "$seqDepthFile")
  if [[ "$seqDepth" -gt "1" ]]; then
    scale_factor=$(echo "10000 / ${seqDepth}" | bc -l)
    bedtools genomecov -bg -scale $scale_factor -i $projPath/results/alignment/bed/${histName}_bowtie2.fragments.bed -g $chromSize > ${projPath}/results/alignment/bedgraph/${histName}.fragments.normalized.bedgraph
  fi

  # Peak calling
  seacr="/path/to/SEACR.sh"
  bash ${seacr} ${projPath}/results/alignment/bedgraph/${histName}.fragments.normalized.bedgraph 0.05 non stringent $projPath/results/peakCalling/SEACR/${histName}_seacr_top0.05.peaks
done < $sampleList
