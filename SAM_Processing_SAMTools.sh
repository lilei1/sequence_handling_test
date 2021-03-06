#!/bin/bash

#PBS -l mem=12gb,nodes=1:ppn=8,walltime=24:00:00
#PBS -m abe
#PBS -M llei@umn.edu
#PBS -q lab

set -e
set -u
set -o pipefail

module load parallel

#   This script is a QSub submission script for converting a batch of SAM files to BAM files
#   To use, on line 5, change the 'user@example.com' to your own email address
#       to get notificaitons on start and completion of this script
#   Define a path to SAMTools in the 'SAMTOOLS' field on line 44
#       If using MSI, leave the definition as is to use their installation of SAMTools
#       Otherwise, it should look like this:
#           SAMTOOLS=${HOME}/software/samtools
#       Please be sure to comment out (put a '#' symbol in front of) the 'module load samtools' on line 43
#       And to uncomment (remove the '#' symbol) from lines 44 and 45
#   Add the full file path to list of samples on the 'SAMPLE_INFO' field on line 48
#       This should look like:
#           SAMPLE_INFO=${HOME}/Directory/list.txt
#       Use ${HOME}, as it is a link that the shell understands as your home directory
#           and the rest is the full path to the actual list of samples
#   Define a path to a reference genome on line 51
#       This should look like:
   #        REF_GEN=${HOME}/Directory/reference_genome.fa
#   Put the full directory path for the output in the 'SCRATCH' field on line 54
#       This should look like:
#           SCRATCH="${HOME}/Out_Directory"
#       Adjust for your own out directory.
#   Name the project in the 'PROJECT' field on line 57
#       This should look lke:
#           PROJECT=Barley
#   Run this script with the qsub command
#       qsub SAM_Processing_SAMTools.sh
#   This script outputs sorted and deduplicated BAM files for each sample

#   Load the SAMTools module for MSI, else define path to SAMTools installation
module load samtools/0.1.16
#SAMTOOLS=
#export PATH=$PATH:${SAMTOOLS}

#   List of SAM files for conversion
SAMPLE_INFO=${HOME}/Deleterious_mutation_project/mapping_sample/wild_accs_bam/sam_file_list

#   Reference genome to help base the conversion off of
REF_GEN=${HOME}/Shared/shared/References/Reference_Sequences/Barley/Morex/Morex_Reference.fasta

#   Scratch directory, for output
SCRATCH=${HOME}/Deleterious_mutation_project/mapping_sample/wild_accs_bam/

#   Name of project
PROJECT=Barley

#   Make the outdirectories
OUT=${SCRATCH}/${PROJECT}/SAM_Processing
mkdir -p ${OUT}/stats ${OUT}/deduped ${OUT}/sorted ${OUT}/raw ${OUT}/fixed_header

#   Check to see if SAMTools is installed
if `command -v samtools > /dev/null 2> /dev/null`
then
    echo "SAMTools is installed"
else
    echo "Please install SAMTools and place in your path"
    exit 1
fi

#   Define a function to process SAM files generated by BWA mem
#       This will create the following files for each input SAM file:
#           a sorted BAM file
#           alignment statistics for the sorted BAM file
#           a deduplicated BAM file
#           alignment statistics for the deduplicated BAM file
function process_sam() {
    #   Today's date
    YMD=`date +%Y-%m-%d`
    SAMFILE="$1"
    REF_SEQ="$2"
    OUTDIR="$3"
    #   Sample name, taken from full name of SAM file
    SAMPLE_NAME=`basename "${SAMFILE}" .sam`
    #   Remove unnecessary information from @PG line
    #   Could use sed's in-place option, but that fails on some systems
    #   This method bypasses that
    sed 's/-R.*$//' "${SAMFILE}" > "${OUTDIR}"/fixed_header/"${SAMPLE_NAME}"_FixedHeader.sam
    #   Generate a sorted BAM file
    samtools view -bhT "${REF_SEQ}" "${OUTDIR}"/fixed_header/"${SAMPLE_NAME}"_FixedHeader.sam > "${OUTDIR}/raw/${SAMPLE_NAME}_${YMD}_raw.bam"
    #   Create alignment statistics for the raw BAM file
    samtools flagstat "${OUTDIR}/raw/${SAMPLE_NAME}_${YMD}_raw.bam" > "${OUTDIR}/stats/${SAMPLE_NAME}_${YMD}_raw_stats.out"
    #   Sort the raw BAM file
    samtools sort "${OUTDIR}/raw/${SAMPLE_NAME}_${YMD}_raw.bam" "${OUTDIR}/sorted/${SAMPLE_NAME}_${YMD}_sorted"
    #   Deduplicate the sorted BAM file
    samtools rmdup "${OUTDIR}/sorted/${SAMPLE_NAME}_${YMD}_sorted.bam" "${OUTDIR}/deduped/${SAMPLE_NAME}_${YMD}_deduped.bam"
    #   Create alignment statistics using SAMTools
    samtools flagstat "${OUTDIR}/deduped/${SAMPLE_NAME}_${YMD}_deduped.bam" > "${OUTDIR}/stats/${SAMPLE_NAME}_${YMD}_deduped_stats.out"
}

#   Export the SAM file processing function to be used by GNU Parallel
export -f process_sam

#   Run the SAM file processing function in parallel for all input SAM files
cat ${SAMPLE_INFO} | parallel "process_sam {} ${REF_GEN} ${OUT}"

#   Add read groups and merge the BAM files
samtools merge -r ${OUT}/${PROJECT}_Merged.bam `find "${OUT}"/deduped -name "*_deduped.bam"`

#   Create a list of deduped BAM files
find ${OUT}/deduped -name "*_deduped.bam" | sort > ${OUT}/${PROJECT}_finished_BAM_files.txt
echo "List of BAM files can be found at"
echo "${OUT}/${PROJECT}_finished_BAM_files.txt"
echo
echo "Merged BAM file can be found at"
echo "${OUT}/${PROJECT}_Merged.bam"
