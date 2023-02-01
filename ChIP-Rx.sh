###file_path
input_path=$1
file_path=`pwd`
raw_dir=${file_path}/00_rawdata
logs_dir=${file_path}/logs
fastqc_dir=${file_path}/01_fastqc
trimmedFastq_dir=${file_path}/02_trimmeddata
trimmedFastq_log_dir=${file_path}/logs/trimmeddata
trimmed_fastqc_dir=${file_path}/02_trimmed_fastqc
align_exp_dir=${file_path}/03_bam_hg19
alignexp_log_dir=${file_path}/logs/align_hg19
align_spike_dir=${file_path}/03_spikebam_mm10
alignspike_log_dir=${file_path}/logs/align_mm10
exp_bam_rmdup=${file_path}/03_bam_hg19_rmdup
rmdup_exp_log=${file_path}/logs/rmdup_state_hg19
spike_bam_rmdup=${file_path}/03_spikebam_mm10_rmdup
rmdup_spike_log=${file_path}/logs/rmdup_state_mm10
bw_dir=${file_path}/04_bw_rmdup
sampleinfo=${file_path}/sample_info.txt
sample_spikeinfo_dir=${file_path}/Pair_input_IP


#reference genome
GENOME_EXP="/share/home/Lemon/index/bowtie2-hg19/hg19"
GENOME_SPIKE="/share/home/Lemon/index/bowtie2-mm10/mm10"
SPIKE_PREFIX="mm10"

MAPQ=30


#step 0
####check files 

#step 1.1
####change file name#####
mkdir -p ${raw_dir}
ln -s `ls ${input_path}/*.gz` ${raw_dir}/ >/dev/null 2>&1
#
cat $sampleinfo| while read id;
do
arr=($id)
sample1=${arr[0]}
sample2=${arr[1]}
echo "sample1:$sample1 sample2:$sample2"
fq1=$(ls ${raw_dir}/*R1.f*q.gz |grep "$sample1") 
fq2=$(ls ${raw_dir}/*R2.f*q.gz |grep "$sample1") 
mv $fq1 ${raw_dir}/${sample2}_R1.fastq.gz
mv $fq2 ${raw_dir}/${sample2}_R2.fastq.gz
done
#
#
##step 1.2
#####fastqc of raw data ####
mkdir -p ${fastqc_dir}
#
in_path=${raw_dir}
out_path=${fastqc_dir}
#
for FILE in `ls ${in_path}/*.gz`
do 
	if [ ! -s ${out_path}/"$(basename ${FILE/.fastq.gz/_fastqc.zip})" ]
	then
	   fastqc $FILE -t 2 -o ${out_path}/ &
	fi
done
wait
#
##step 1.3
##check fq
cd ${fastqc_dir}
ls *fastqc.zip|sed 's/_R[1-2]_fastqc.zip//g'|sort|uniq -c
#
#### merge reports of fastqc
~/software/anaconda3/envs/py36/bin/multiqc ${fastqc_dir}/ -n rawdata_multiqc -o ${fastqc_dir}/
#
#
##step 2.1
#### Trimming adapters  (trim_galore)
#####remember to change --length for proseq !!!!!! #####
#
mkdir -p ${trimmedFastq_dir}
mkdir -p ${trimmedFastq_log_dir}
#
for fq1 in `ls ${raw_dir}/*R1.fastq.gz`

do
fq2=${fq1/R1.fastq.gz/R2.fastq.gz}
    if [ ! -s ${trimmedFastq_dir}/"$(basename ${fq1/.fastq.gz/_val_1.fq.gz})" ]
    then
    ~/software/anaconda3/envs/py36/bin/trim_galore --paired -o ${trimmedFastq_dir} $fq1 $fq2 --path_to_cutadapt ~/software/anaconda3/envs/py36/bin/cutadapt \
    > ${trimmedFastq_log_dir}/"$(basename ${fq1/_R1.fq.gz/_trimmed.log})" 2>&1 & 
    fi
done
wait
#
##step 2.4
######qc for trimmed data####
mkdir -p ${trimmed_fastqc_dir}
#
for FILE in `ls ${trimmedFastq_dir}/*fq.gz`
do
if [ ! -s ${trimmed_fastqc_dir}/"$(basename ${FILE/.fq.gz/_fastqc.html})" ]
then
 fastqc $FILE -t 2 -o ${trimmed_fastqc_dir} &
fi
done
wait
#### merge reports of fastqc
~/software/anaconda3/envs/py36/bin/multiqc ${trimmed_fastqc_dir}/ -n trimmed_multiqc -o ${trimmed_fastqc_dir}/


#step 3.1
###Aligning to experimental genome#####
mkdir ${align_exp_dir}
mkdir ${alignexp_log_dir}

for PAIR in $(ls ${trimmedFastq_dir} | sed 's/_R[1-2].*//' | sort | uniq )
do
if [ ! -s "${align_exp_dir}/${PAIR}_hg19.bam" ]
then
    echo "aligning ${PAIR} to experimental genome"
    (bowtie2 \
    -N 1 \
    --threads 20 \
    -x "$GENOME_EXP" \
    -1 "${trimmedFastq_dir}/${PAIR}_R1_val_1.fq.gz" \
    -2 "${trimmedFastq_dir}/${PAIR}_R2_val_2.fq.gz" \
    2> ${alignexp_log_dir}/${PAIR}_align.log) |
    samtools view -bS -F 3844 -q ${MAPQ} |
    samtools sort -@ 20 -o ${align_exp_dir}/${PAIR}_hg19.bam
    samtools index ${align_exp_dir}/${PAIR}_hg19.bam
fi
done




#step 3.2
### Aligning to spike in genome to get normalization factors ###
mkdir -p ${align_spike_dir}
mkdir -p ${alignspike_log_dir}

for PAIR in $(ls ${trimmedFastq_dir} | sed 's/_R[1-2].*//' | sort | uniq )
do
if [ ! -s "${align_spike_dir}/${PAIR}_onlymm10.bam" ]
    then
    echo "aligning ${PAIR} to spike-in genome"
    (bowtie2 \
    -N 1 \
    --threads 20 \
    -x "$GENOME_SPIKE" \
    -1 "${trimmedFastq_dir}/${PAIR}_R1_val_1.fq.gz" \
    -2 "${trimmedFastq_dir}/${PAIR}_R2_val_2.fq.gz" \
    2> ${alignspike_log_dir}/${PAIR}_onlymm10Align.log) |
    samtools view -bS -F 3844 -q ${MAPQ} |
    samtools sort -@ 20 -o ${align_spike_dir}/${PAIR}_onlymm10.bam
    samtools index ${align_spike_dir}/${PAIR}_onlymm10.bam
fi
done





# step 3.3
# remove dulipate
##remove duplicates of exp genome
mkdir -p ${exp_bam_rmdup}
mkdir -p ${rmdup_exp_log}

cd ${exp_bam_rmdup}
ls ${align_exp_dir}/*.bam|while read id;
do
sample=$(basename ${id/.bam/})
if [ ! -s "${exp_bam_rmdup}/${sample}.rmdup.bam" ]
    then
        java -jar ~/software/picard.jar MarkDuplicates -REMOVE_DUPLICATES True \
        -I $id \
        -O ${exp_bam_rmdup}/${sample}.rmdup.bam \
        -M ${exp_bam_rmdup}/${sample}.rmdup.metrics
        samtools index ${exp_bam_rmdup}/${sample}.rmdup.bam
        samtools flagstat ${exp_bam_rmdup}/${sample}.rmdup.bam > ${rmdup_exp_log}/${sample}.rmdup.stat
    fi
done

##remove duplicates of spike-in genome
mkdir -p ${spike_bam_rmdup}
mkdir -p ${rmdup_spike_log}

cd ${spike_bam_rmdup}
ls ${align_spike_dir}/*.bam|while read id;
do
sample=$(basename ${id/.bam/})
if [ ! -s "${spike_bam_rmdup}/${sample}.rmdup.bam" ]
    then
        java -jar ~/software/picard.jar MarkDuplicates -REMOVE_DUPLICATES True \
        -I $id \
        -O ${spike_bam_rmdup}/${sample}.rmdup.bam \
        -M ${spike_bam_rmdup}/${sample}.rmdup.metrics
        samtools index ${spike_bam_rmdup}/${sample}.rmdup.bam
        samtools flagstat ${spike_bam_rmdup}/${sample}.rmdup.bam > ${rmdup_spike_log}/${sample}.rmdup.stat
    fi
done




#step 3.4
### calculate normalization factors ###
hg19_path=${rmdup_exp_log}
mm10_path=${rmdup_spike_log}
align_path=${alignexp_log_dir}

spike_sample=`head -n 1 ${sample_spikeinfo_dir}/all.input`

spike_hgRead=`sed -n '1p' $hg19_path/$spike_sample"_hg19.rmdup.stat" | cut -d ' ' -f 1`
spike_mmRead=`sed -n '1p' $mm10_path/$spike_sample"_onlymm10.rmdup.stat" | cut -d ' ' -f 1`
spike_product=$(echo $spike_hgRead/$spike_mmRead | bc -l)

echo -e "sample\tALLREADS\tHg19_READS\tHg19_mapping_RATIO\tHg19_qc_READS\tHg19_qc_RATIO\tMM10_qc_READS\tMM10_qc_RATIO_intotal\tMM10_qc_RATIO_inqc\tspike_factor" >${sample_spikeinfo_dir}/spike-input.txt
cat ${sample_spikeinfo_dir}/all.input |while read sample;
do
    ALLREADS=$(cat ${align_path}/$sample"_align.log"|grep "were paired; of these:$"|cut -d "(" -f 1|awk '{print $1*2}')
    Hg19_READS=$(cat ${align_path}/$sample"_align.log"| sed 's/%//g' | awk '{printf $0"\t"}'  |cut -f 4,5,8,13,14 | sed 's/\t/\n/g' | awk '{print $1}' | awk '{printf $0"\t"}'|awk '{print 2*($1+$2+$3)+$4+$5}')
    Hg19_mapping_RATIO=$(cat ${align_path}/$sample"_align.log"|grep "overall alignment rate"|cut -d "%" -f 1)
    Hg19_qc_READS=$(cat $hg19_path/${sample}"_hg19.rmdup.stat"|grep "total (QC-passed reads"|cut -d " " -f 1)
    Hg19_qc_RATIO=$(echo "${Hg19_qc_READS}/${ALLREADS}"|bc -l)
    MM10_qc_READS=$(cat ${mm10_path}/${sample}"_onlymm10.rmdup.stat"|grep "total (QC-passed reads"|cut -d "+" -f 1)
    QC_reads=$(echo "${MM10_qc_READS}+${Hg19_qc_READS}"|bc )
    MM10_qc_RATIO_intotal=$(echo "${MM10_qc_READS}/${ALLREADS}"|bc -l)
    MM10_qc_RATIO_inqc=$(echo "${MM10_qc_READS}/${QC_reads}"|bc -l)
    spikefactor=$(echo "($Hg19_qc_READS/$MM10_qc_READS)/$spike_product"|bc -l)
    echo -e $sample"\t"$ALLREADS"\t"$Hg19_READS"\t"$Hg19_mapping_RATIO"\t"$Hg19_qc_READS"\t"$Hg19_qc_RATIO"\t"$MM10_qc_READS"\t"$MM10_qc_RATIO_intotal"\t"$MM10_qc_RATIO_inqc"\t"$spikefactor >> ${sample_spikeinfo_dir}/spike-input.txt
done


echo -e "sample\tALLREADS\tHg19_READS\tHg19_mapping_RATIO\tHg19_qc_READS\tHg19_qc_RATIO\tMM10_qc_READS\tMM10_qc_RATIO_intotal\tMM10_qc_RATIO_inqc\tScaleFactor" > ${sample_spikeinfo_dir}/scalefactor.txt
for file in ${sample_spikeinfo_dir}/ip.pair*;
do
    input=`head -n 1 $file`
    spike_factor=$(grep $input ${sample_spikeinfo_dir}/spike-input.txt | cut -f 10)
    cat $file|while read sample;do

    ALLREADS=$(cat ${align_path}/$sample"_align.log"|grep "were paired; of these:$"|cut -d "(" -f 1|awk '{print $1*2}')
    Hg19_READS=$(cat ${align_path}/$sample"_align.log"| sed 's/%//g' | awk '{printf $0"\t"}'  |cut -f 4,5,8,13,14 | sed 's/\t/\n/g' | awk '{print $1}' | awk '{printf $0"\t"}'|awk '{print 2*($1+$2+$3)+$4+$5}')
    Hg19_mapping_RATIO=$(cat ${align_path}/$sample"_align.log"|grep "overall alignment rate"|cut -d "%" -f 1)
    Hg19_qc_READS=$(cat $hg19_path/${sample}"_hg19.rmdup.stat"|grep "total (QC-passed reads"|cut -d " " -f 1)
    Hg19_qc_RATIO=$(echo "${Hg19_qc_READS}/${ALLREADS}"|bc -l)
    MM10_qc_READS=$(cat ${mm10_path}/${sample}"_onlymm10.rmdup.stat"|grep "total (QC-passed reads"|cut -d "+" -f 1)
    QC_reads=$(echo "${MM10_qc_READS}+${Hg19_qc_READS}"|bc )
    MM10_qc_RATIO_intotal=$(echo "${MM10_qc_READS}/${ALLREADS}"|bc -l)
    MM10_qc_RATIO_inqc=$(echo "${MM10_qc_READS}/${QC_reads}"|bc -l)
    SCALEFACTOR=$(echo "1/((${MM10_qc_READS}/1000000)*$spike_factor)" | bc -l )

    echo -e $sample"\t"$ALLREADS"\t"$Hg19_READS"\t"$Hg19_mapping_RATIO"\t"$Hg19_qc_READS"\t"$Hg19_qc_RATIO"\t"$MM10_qc_READS"\t"$MM10_qc_RATIO_intotal"\t"$MM10_qc_RATIO_inqc"\t"$SCALEFACTOR >> ${sample_spikeinfo_dir}/scalefactor.txt
    done
done

echo "scale factor is done"



#step 4.1
### Making CPM-normalized bigWig files with reads adjusted by spike-in ###
mkdir -p ${bw_dir}

cat  ${sample_spikeinfo_dir}/scalefactor.txt | sed '1d' | while read id;
do
arr=($id)
sample=${arr[0]}
scalefactor=${arr[9]}
bam_file=${exp_bam_rmdup}/${sample}_hg19.rmdup.bam
    if [ ! -s "${bw_dir}/${sample}_rmdup.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --blackListFileName ~/index/genomes/hg19-blacklist.v2.bed \
        --outFileName ${bw_dir}/${sample}_rmdup.bw \
        --binSize 1 \
        --scaleFactor $scalefactor \
        --numberOfProcessors 20 \
        --normalizeUsing None
    fi
done

#step 4.2
### Making CPM-normalized bigWig files without spike-in ###
cat ${sample_spikeinfo_dir}/scalefactor.txt | sed '1d' | while read id;
do
arr=($id)
sample=${arr[0]}
bam_file=${exp_bam_rmdup}/${sample}_hg19.rmdup.bam
    if [ ! -s "${bw_dir}/${sample}_rmdup_CPM.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --blackListFileName ~/index/genomes/hg19-blacklist.v2.bed \
        --outFileName ${bw_dir}/${sample}_rmdup_CPM.bw \
        --binSize 1 \
	--scaleFactor 1 \
        --numberOfProcessors 20 \
        --normalizeUsing CPM
    fi
done

#step 5.1
### Call Peaks without input ###
peak_dir=${file_path}/05_peaks
peak_log=${file_path}/logs/peaks
mkdir -p ${peak_dir}/peak_alone
ls ${exp_bam_rmdup}/*_hg19.rmdup.bam | grep -v 'nput' | while read id;
do
    sample=`basename $id | sed 's/_hg19.rmdup.bam//g'`
    bam_file=${exp_bam_rmdup}/${sample}_hg19.rmdup.bam
    blacklist="/share/home/Lemon/index/genomes/hg19-blacklist.v2.bed"
    mkdir -p ${peak_log}/${sample}
    if [ ! -f ${peak_dir}/peak_alone/${sample}_peaks.narrowPeak ]
    then
        macs2 callpeak -t ${bam_file} -f BAMPE -n ${peak_dir}/peak_alone/$sample -g hs --nomodel 2>${peak_log}/${sample}/${sample}_alone.log
        bedtools intersect -a ${peak_dir}/peak_alone/${sample}_peaks.narrowPeak -b $blacklist -f 0.25 -v > ${peak_dir}/peak_alone/${sample}_peaks.final.narrowPeak &
    fi
done
wait








