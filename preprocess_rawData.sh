#!/bin/bash 
#preprocess_rawData.sh
#SBATCH --job-name=preprocess_rawData
#SBATCH --partition=long
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=32 
#SBATCH --mem=100gb 
#SBATCH --time=10-00:00:00 
#SBATCH --output=epimut_raw2vc%j.log 

<< "epimut_raw2vc.sh"
For each sample:
- Clean the raw fq.gz with bbduck.sh, 
- align to TuMV reference with bwa mem
- mark duplicates with gatk4 (dedup)
- Generates a file without duplicates (nodup)
- Lofreq on complete alignment (dedup)
- LoFreq on no duplicates alignment (nodup)
 
MJ Olmo-Uceda
2022/06/02
epimut_raw2vc.sh

# paths to the samples file to analyse: one samples per line
samples=$1

### Setting directories
dir="/storage/evsysvir/EvolucionExperimentalEnEpimutantes"
REF="/storage/evsysvir/REFERENCES"
adapters=$REF"/novogen_adapters_epiMut.fasta"
ref_vir=$REF"/TuMV_Anc/TuMV_Anc.fasta"

dir_fq=$dir"/raw_data"
dir_clean_fq=$dir"/clean_data"
dir_discarded=$dir"/discarded_data"

dir_alignment=$dir"/alignments"
dir_sorted=$dir_alignment"/sorted"
dir_dedup=$dir_alignment"/dedup"
dir_metrics=$dir_alignment"/metrics"
dir_nodup=$dir_alignment"/nodup"

dir_lofreq_results=$dir"/variants_lofreq"
dir_lofreq_vcf=$dir_lofreq_results"/vcf"
dir_lofreq_nodup=$dir_lofreq_results"/vcf_nodup"

dir_qc1=$dir"/QC1"

# path to gatk4
gatk="/home/maolu/programs/gatk-4.2.2.0/gatk"

# parameters 
minlength=80
threads=20

################################################################################
# MAIN PROGRAM
################################################################################

start=`date +%s`

# PROCESS WITH EACH SAMPLE CONTAINED IN $SAMPLES
while read sample
do
    echo "****************** PROCESSING SAMPLE $sample ******************"
    echo

    echo "········································"
    echo "@--> 1. Preprocessing raw_fastqs"
    echo "········································"

    bbduk.sh \
    in=${dir_fq}/${sample}/${sample}_1.fq.gz \
    in2=${dir_fq}/${sample}/${sample}_2.fq.gz \
    out=${dir_clean_fq}/${sample}_clean_1.fq.gz \
    out2=${dir_clean_fq}/${sample}_clean_2.fq.gz \
    outm=${dir_discarded}/${sample}_dis.fq \
    ref=${adapters} \
    ktrim=r k=31 mink=11 \
    qtrim=r trimq=10 maq=5 minlength=$minlength \
    forcetrimleft=10
    echo

    echo "········································"
    echo "@--> 2. Mapping with BWA MEM"
    echo "········································"

    bwa mem -o ${dir_alignment}/${sample}_tumv.sam \
    ${ref_vir} \
    ${dir_clean_fq}/${sample}_clean_1.fq.gz \
    ${dir_clean_fq}/${sample}_clean_2.fq.gz

    echo "········································"
    echo "@--> Sorting"
    echo "········································"
    # sam to sorted bam & remove non esential files
    samtools view -bS ${dir_alignment}/${sample}_tumv.sam > ${dir_alignment}/${sample}_tumv.bam
    samtools sort ${dir_alignment}/${sample}_tumv.bam -o ${dir_sorted}/${sample}_tumv.sorted.bam
    rm ${dir_alignment}/${sample}_tumv.sam ${dir_alignment}/${sample}_tumv.bam
    
    echo "········································"
    echo "@--> Mark duplicates"
    echo "········································"
    $gatk MarkDuplicates \
    -I ${dir_sorted}/${sample}_tumv.sorted.bam \
    -O ${dir_dedup}/${sample}_tumv.dedup.bam \
    -M ${dir_metrics}/${sample}_metrics.txt

    echo "········································"
    echo "@--> Remove duplicates"
    echo "········································"
    samtools view -hbF0x400 ${dir_dedup}/${sample}_tumv.dedup.bam > ${dir_nodup}/${sample}_tumv.nodup.bam
  
    echo "········································"
    echo "@--> Variant Calling with LoFreq "
    echo "········································"    
    lofreq call -f ${ref_vir} -o ${dir_lofreq_vcf}/${sample}.vcf ${dir_sorted}/${sample}_tumv.sorted.bam 
    lofreq call -f ${ref_vir} -o ${dir_lofreq_nodup}/${sample}.nodup.vcf ${dir_nodup}/${sample}_tumv.nodup.bam 

done < $samples

# Tiempo total de ejecución
end=`date +%s`
runtime=$((end-start))
echo "Total time: $runtime s"
