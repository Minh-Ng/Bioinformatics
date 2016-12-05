#!/bin/bash
FASTQC_OUT_DIR='fastqc_results'
CUFF_OUT_DIR='cufflinks_results'
GENOME_OUT_DIR='genome_results/'

#Make output directory if it does not exist
echo "Cleaning previous output..."
rm -f *.fasta *.sam *.bam *.out *.fpkm_tracking *.txt Genome SA SAindex *.tab *.gtf
rm -rf $FASTQC_OUT_DIR $CUFF_OUT_DIR $GENOME_OUT_DIR *_STARtmp

while getopts ":c" opt; do
  case $opt in
        c)
          echo "Done"
          exit 0
          ;;
        \?)
          echo "Invalid option: -$OPTARG" >&2
          exit 1
          ;;
  esac
done

mkdir -p $FASTQC_OUT_DIR
mkdir -p $CUFF_OUT_DIR
mkdir -p $GENOME_OUT_DIR

fastqc -o $FASTQC_OUT_DIR ./inputs/*.fastq.gz
if [ $? -ne 0 ]
then
  echo "fastqc failed"
  exit 1
else
  echo "fastqc results generated"
fi

space=" "
fastqz=""
for i in ./inputs/*.fastq.gz; do
  if [ -z "$fastqz" ]; then
    fastqz=$i
  else
    fastqz=$fastqz$space$i
  fi
done

fasta=""
for i in ./inputs/*.fasta; do
  if [ -z "$fasta" ]; then
    fasta=$i
  else
    fasta=$fasta$space$i
  fi
done

echo "Generating genome from fasta files..."
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ./$GENOME_OUT_DIR --genomeFastaFiles $fasta --genomeSAsparseD 2 --genomeSAindexNbases 14 --limitGenomeGenerateRAM 15000000000

if [ $? -ne 0 ]
then
  echo "STAR genome generation failed, most likely due to lack of ram"
  exit 1
else
  echo "Running star alignment..."
fi

prefix="experiment"
STAR --genomeDir ./$GENOME_OUT_DIR --runThreadN 8 --readFilesIn $fastqz --readFilesCommand gunzip -c --outFileNamePrefix $prefix

if [ $? -ne 0 ]
then
  echo "STAR alignment failed"
  exit 1
else
  echo "Sorting sam file from STAR alignment..."
fi

sortedSam="Aligned.out.sorted.sam"
aligned="Aligned.out.sam"
samtools sort -O sam -T $prefix$aligned -o $prefix$sortedSam $prefix$aligned

if [ $? -ne 0 ]
then
  echo "sam file sort failed"
  exit 1
else
  echo "Analyzing sam file..."
fi

cufflinks -o $CUFF_OUT_DIR $sortedSam

if [ $? -ne 0 ]
then
  echo "cufflinks run failed"
  exit 1
else
  echo "Results generated in $CUFF_OUT_DIR"
fi
