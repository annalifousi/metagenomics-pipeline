#!/bin/bash

# Set up environment variables
PROJECT_DIR="/home/projects/co_23260/data/groups/group_6"
TRIM_DIR="$PROJECT_DIR/Trim"
KMA_DIR="$PROJECT_DIR/kma"
SILVA_DIR="$PROJECT_DIR/silva"
ASSEMBLY_DIR="$PROJECT_DIR/assemblies"
MAPPED_DIR="$PROJECT_DIR/mapped"
BINS_DIR="$PROJECT_DIR/bins"
CHECKM2_DIR="$PROJECT_DIR/checkm2"
DREP_DIR="$PROJECT_DIR/drep2"
GTDB_DIR="$PROJECT_DIR/gtdb"

# Samples
SAMPLES=("1001-17-001" "6009-17-001" "8020-17-001" "9020-17-001")

# Step 1: Pre-processing with FASTP
echo "Step 1: Quality trimming and adapter removal"
for SAMPLE in "${SAMPLES[@]}"; do
    fastp -i "$PROJECT_DIR/data/${SAMPLE}_R1_001.fastq.gz" \
          -I "$PROJECT_DIR/data/${SAMPLE}_R2_001.fastq.gz" \
          -o "$TRIM_DIR/${SAMPLE}_R1.trim.fq.gz" \
          -O "$TRIM_DIR/${SAMPLE}_R2.trim.fq.gz" \
          --html "$TRIM_DIR/${SAMPLE}_fastp.html" \
          --json "$TRIM_DIR/${SAMPLE}_fastp.json"
done

# Step 2: Quality control with FASTQC
echo "Step 2: Running FASTQC for quality check"
for SAMPLE in "${SAMPLES[@]}"; do
    fastqc -o "$PROJECT_DIR/fastqc" "$TRIM_DIR/${SAMPLE}_R1.trim.fq.gz" "$TRIM_DIR/${SAMPLE}_R2.trim.fq.gz"
done

# Step 3: MultiQC report generation
echo "Step 3: Generating MultiQC report"
multiqc -o "$PROJECT_DIR/multiqc" "$PROJECT_DIR/fastqc"

# Step 4: Alignment and mapping with KMA
echo "Step 4: Mapping with KMA (Resfinder and Silva)"
for SAMPLE in "${SAMPLES[@]}"; do
    # Resfinder
    kma -ipe "$TRIM_DIR/${SAMPLE}_R1.trim.fq.gz" "$TRIM_DIR/${SAMPLE}_R2.trim.fq.gz" \
        -o "$KMA_DIR/kma_${SAMPLE}_resfinder" -ef -1t1 -cge -nf \
        -t_db "$PROJECT_DIR/databases/resfinder_db/ResFinder"
    
    # Silva
    kma -ipe "$TRIM_DIR/${SAMPLE}_R1.trim.fq.gz" "$TRIM_DIR/${SAMPLE}_R2.trim.fq.gz" \
        -o "$SILVA_DIR/silva_${SAMPLE}" -t 40 -ef -mem_mode -cge \
        -t_db "/home/projects/co_23260/data/databases/kma_silva_v138/Silva_20200116"
done

# Step 5: Genome assembly with SPAdes
echo "Step 5: De novo assembly with SPAdes"
for SAMPLE in "${SAMPLES[@]}"; do
    spades.py -1 "$TRIM_DIR/${SAMPLE}_R1.trim.fq.gz" \
              -2 "$TRIM_DIR/${SAMPLE}_R2.trim.fq.gz" \
              -o "$ASSEMBLY_DIR/$SAMPLE"
done

# Step 6: Alignment and binning
echo "Step 6: Alignment and binning"
for SAMPLE in "${SAMPLES[@]}"; do
    # Alignment with bbmap
    bbmap.sh in1="$TRIM_DIR/${SAMPLE}_R1.trim.fq.gz" \
             in2="$TRIM_DIR/${SAMPLE}_R2.trim.fq.gz" \
             ref="$ASSEMBLY_DIR/$SAMPLE/scaffolds.fasta" \
             out="$MAPPED_DIR/${SAMPLE}.sam"
    
    # Convert SAM to BAM and sort
    samtools view -bS "$MAPPED_DIR/${SAMPLE}.sam" | samtools sort -o "$MAPPED_DIR/${SAMPLE}.bam"
    rm "$MAPPED_DIR/${SAMPLE}.sam" # Clean up SAM file
    
    # Binning with Metabat2
    metabat2 -i "$ASSEMBLY_DIR/$SAMPLE/scaffolds.fasta" \
             -a "$MAPPED_DIR/${SAMPLE}.bam" \
             -o "$BINS_DIR/$SAMPLE/bin"
done

# Step 7: Genome completeness check with CheckM2
echo "Step 7: Genome completeness check with CheckM2"
checkm2 analyze --input_dir "$BINS_DIR" --output_dir "$CHECKM2_DIR"

# Step 8: Genome comparison with dRep
echo "Step 8: Genome comparison with dRep"
dRep dereplicate "$DREP_DIR" -g "$BINS_DIR/*/*.fa" --processors 40

# Step 9: Taxonomic identification with GTDB-tk
echo "Step 9: Taxonomic identification with GTDB-tk"
for SAMPLE in "${SAMPLES[@]}"; do
    gtdbtk classify_wf --genome_dir "$BINS_DIR/$SAMPLE" --out_dir "$GTDB_DIR/$SAMPLE" --cpus 40
done

# Step 10: Calculate ALR for resistance gene coverage
echo "Step 10: Calculating ALR for resistance gene coverage"
FPKMS=(28269031 19030724 25191519 18648205)
for i in "${!SAMPLES[@]}"; do
    ALR=$(awk "BEGIN {print log(${FPKMS[$i]}/1000000)}")
    echo "Sample ${SAMPLES[$i]}: ALR = $ALR"
done

# Completion message
echo "Pipeline execution completed successfully!"
