# Quality control on raw FASTQ files
fastqc SO_4933_C_CHR1_R1.fastq SO_4933_C_INR1_R1.fastq

# Deduplicate treatment sample (CHR1) using fastp
fastp -i SO_4933_C_CHR1_R1.fastq \
      -o SO_4933_C_CHR1_R1_dedup.fastq \
      --dedup \
      --dup_calc_accuracy 3
 
# Deduplicate control sample (INR1) using fastp
fastp -i SO_4933_C_INR1_R1.fastq \
      -o SO_4933_C_INR1_R1_dedup.fastq \
      --dedup \
      --dup_calc_accuracy 3
      
# Quality control on deduplicated FASTQ files
fastqc SO_4933_C_CHR1_R1_dedup.fastq SO_4933_C_INR1_R1_dedup.fastq

# Build bowtie2 index from reference genome
bowtie2-build 'E.coli_BW25113.fasta' genome_index

# Align treatment sample to reference genome
bowtie2 -x genome_index -U SO_4933_C_CHR1_R1_dedup.fastq -S C_CHR1_R1.sam

# Align control sample to reference genome
bowtie2 -x genome_index -U SO_4933_C_INR1_R1_dedup.fastq -S C_INR1_R1.sam

# Convert SAM to BAM format for treatment sample
samtools view -bS C_CHR1_R1.sam > C_CHR1_R1.bam

# Convert SAM to BAM format for control sample
samtools view -bS C_INR1_R1.sam > C_INR1_R1.bam

# Alternative de-duplication method in case we didnt de-duplicate with fastp before
# samtools markdup -r SO_4933_C_CHR1_sorted.bam SO_4933_C_CHR1_dedup.bam
# samtools markdup -r SO_4933_C_INR1_sorted.bam SO_4933_C_INR1_dedup.bam

# Peak calling using MACS3 with treatment vs control comparison
# g denotes E.coli genome size ~4.6 Mb. Do with conda activate conda activate macs2-env (from home dir) for macs 2.
/home/ibab/miniconda3/bin/macs3 callpeak \
    -t C_CHR1_R1.bam \
    -c C_INR1_R1.bam \
    -f BAM \
    -g 4.6e6 \
    -n ecoli_basic \
    --outdir macs3_results \
    --bdg \
    -q 0.05
    
# View in IGV (input files -> .bed, .narrowpeak, lambda.bdg, pileup.bdg.) Autoscale the two .bdg files to see peak.

# Merge nearby peaks within 500bp using HOMER
mergePeaks -d 500 macs3_results/ecoli_basic_peaks.narrowPeak > macs3_results/C_mergePeak_out.txt

# Count number of merged peaks
wc -l macs3_results/C_mergePeak_out.txt

# Convert merged peaks from txt to BED format
awk '!/^#/ && !/^Total/ && !/^name/ && /^Merged/ {print $2"\t"$3"\t"$4"\t"$1}' macs3_results/C_mergePeak_out.txt > macs3_results/C_mergePeak_out.bed

# Annotate peaks with genomic features using HOMER
annotatePeaks.pl macs3_results/ecoli_basic_peaks.narrowPeak E.coli_BW25113.fasta -gff E.coli_BW25113_annotation.gff > C_annotatedPeak.tsv

# Find closest genes to peaks using BEDTools
closestBed -a macs3_results/ecoli_basic_peaks.narrowPeak -b E.coli_BW25113_annotation.bed | awk -F"\t" '{if($18=="gene"){print}}' > C_closestGenes.bed

# Extract DNA sequences from peak regions
# Use only first 3 columns from narrowPeak file
# Extract sequences using columns 1-3 (chrom, start, end)
bedtools getfasta -fi E.coli_BW25113.fasta \
                  -bed <(cut -f1-3 macs3_results/ecoli_basic_peaks.narrowPeak) \
                  -fo C_peak_sequences.fasta
                  
# Discover enriched motifs in peak sequences using MEME
meme C_peak_sequences.fasta -dna -oc meme_results -mod zoops -nmotifs 5 -minw 6 -maxw 10 -objfun classic
