# High-Throughput Sequencing Analysis Script for CUT&Tag

## Overview

This script is designed for processing and analyzing CUT&Tag data, specifically for histone modifications, using HPC systems. It automates the alignment, duplicate removal, fragment size calculation, reproducibility assessment, spike-in calibration, and peak calling steps. The script is written in Bash and utilizes widely used bioinformatics tools.

---

## Features

1. **Alignment**
   - Aligns reads to the reference genome and spike-in genome using Bowtie2.
   - Outputs alignment statistics and depth information.

2. **Duplicate Removal**
   - Uses Picard to mark and remove duplicates.

3. **Fragment Analysis**
   - Calculates fragment size distribution.
   - Filters and formats fragments for downstream analysis.

4. **Spike-In Normalization**
   - Calculates scaling factors based on spike-in alignment.
   - Generates normalized BedGraph files.

5. **Reproducibility Assessment**
   - Binning of fragments to assess replicate reproducibility.

6. **Peak Calling**
   - Uses SEACR for peak calling with normalization and stringent thresholds.

---

## Dependencies

The script requires the following modules/tools to be available in the HPC environment:

- Bowtie2
- Samtools
- Bedtools
- Picard
- SEACR
- Java Runtime Environment (JRE)

Make sure the modules are available and can be loaded using the `module load` command.

---

## File Structure

- **Input files**:
  - FASTQ files: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`
  - `sample_names.txt`: A text file containing the sample names (one per line).
  - Reference genome: Path to the reference genome.
  - Spike-in reference genome: Path to the spike-in genome.

- **Output directories**:
  - `results/alignment/sam/`: Contains SAM files from Bowtie2.
  - `results/alignment/bam/`: Contains BAM files for mapped reads.
  - `results/alignment/bed/`: Contains BED files for fragments.
  - `results/alignment/bedgraph/`: Contains normalized BedGraph files.
  - `results/peakCalling/SEACR/`: Contains peak files from SEACR.

---

## Usage

1. **Prepare Input Files**:
   - Place raw FASTQ files in the `data/raw/fastq/` directory.
   - Create a `sample_names.txt` file with sample names (one per line).

2. **Edit Paths**:
   - Update the `projPath`, `ref`, `spikeInRef`, and `chromSize` variables in the script to match your project directory and reference genome paths.

3. **Submit Job**:
   - Submit the script to the HPC scheduler using:
     ```bash
     sbatch script_name.sh
     ```

4. **Monitor Job**:
   - Check job logs (`job_name.out` and `job_name.err`) for progress and errors.

---

## Script

```bash
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
