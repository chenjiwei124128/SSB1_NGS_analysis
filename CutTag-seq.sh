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
align_spike_dir=${file_path}/03_spikebam_dm3
alignspike_log_dir=${file_path}/logs/align_dm3
exp_bam_rmdup=${file_path}/03_bam_hg19_rmdup
rmdup_exp_log=${file_path}/logs/rmdup_state_hg19
spike_bam_rmdup=${file_path}/03_spikebam_dm3_rmdup
rmdup_spike_log=${file_path}/logs/rmdup_state_dm3
bw_dir=${file_path}/04_bw_rmdup
sampleinfo=${file_path}/sample_info.txt
sample_spikeinfo_dir=${file_path}/Pair_input_IP


#reference genome
GENOME_EXP="/share/home/Blueberry/reference/index/bowtie2/hg19/hg19"
GENOME_SPIKE="/share/home/Lemon/index/bowtie2-dm3/dm3"
SPIKE_PREFIX="dm3"

MAPQ=10

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
##step 2.2
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
    -t \
    -L 25 \
    -X 700 \
    --no-mixed \
    --no-discordant \
    --threads 25 \
    -x "$GENOME_EXP" \
    -1 "${trimmedFastq_dir}/${PAIR}_R1_val_1.fq.gz" \
    -2 "${trimmedFastq_dir}/${PAIR}_R2_val_2.fq.gz" \
    2> ${alignexp_log_dir}/${PAIR}_align.log) |
    samtools view -bS -f 2 -q ${MAPQ} |
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
if [ ! -s "${align_spike_dir}/${PAIR}_dm3.bam" ]
    then
    echo "aligning ${PAIR} to spike-in genome"
    (bowtie2 \
    -N 1 \
    --threads 25 \
    -x "$GENOME_SPIKE" \
    --end-to-end \
    --very-sensitive \
    --no-overlap \
    --no-dovetail \
    --no-mixed \
    --no-discordant \
    -L 25 \
    -I 10 -X 700 \
    -1 "${trimmedFastq_dir}/${PAIR}_R1_val_1.fq.gz" \
    -2 "${trimmedFastq_dir}/${PAIR}_R2_val_2.fq.gz" \
    2> ${alignspike_log_dir}/${PAIR}_dm3Align.log) |
    samtools view -bS -f 2 -q ${MAPQ} |
    samtools sort -@ 20 -o ${align_spike_dir}/${PAIR}_dm3.bam
    samtools index ${align_spike_dir}/${PAIR}_dm3.bam
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
	~/software/anaconda3/envs/py36/bin/alignmentSieve --numberOfProcessors 25 --ATACshift -b ${exp_bam_rmdup}/${sample}.rmdup.bam -o ${exp_bam_rmdup}/${sample}.rmdup_Shift.bam
	samtools sort ${exp_bam_rmdup}/${sample}.rmdup_Shift.bam -O bam -@ 25 -o ${exp_bam_rmdup}/${sample}.rmdup_Shift.sorted.bam
    samtools index ${exp_bam_rmdup}/${sample}.rmdup_Shift.sorted.bam
    samtools flagstat ${exp_bam_rmdup}/${sample}.rmdup_Shift.sorted.bam > ${rmdup_exp_log}/${sample}.rmdup_Shift.sorted.stat
    fi
done

# step 3.4
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

#step 3.5
### calculate normalization factors ###
hg19_path=${rmdup_exp_log}
dm3_path=${rmdup_spike_log}
align_path=${alignexp_log_dir}

echo -e "sample\tALLREADS\tHg19_READS\tHg19_mapping_RATIO\tHg19_qc_READS\tHg19_qc_RATIO\tMm10_qc_READS\tMm10_qc_RATIO_intotal\tMm10_qc_RATIO_inqc\tSCALEFACTOR\t" >${logs_dir}/scalefactor.txt
ls $align_path|while read file;
do
    sample=${file/_align.log/}
    ALLREADS=$(cat ${align_path}/$file|grep "were paired; of these:$"|cut -d "(" -f 1|awk '{print $1*2}')
    Hg19_READS=$(cat ${align_path}/$sample"_align.log"| sed 's/%//g' | awk '{printf $0"\t"}'  |cut -f 8,9 | sed 's/\t/\n/g' | awk '{print $1}' | awk '{printf $0"\t"}'|awk '{print 2*($1+$2)}')
    Hg19_mapping_RATIO=$(cat ${align_path}/$file|grep "overall alignment rate"|cut -d "%" -f 1)
    Hg19_qc_READS=$(cat $hg19_path/${sample}_hg19.rmdup.stat|grep "total (QC-passed reads"|cut -d " " -f 1)
    Hg19_qc_RATIO=`printf "%.2f\n" $(echo "${Hg19_qc_READS}/${ALLREADS}*100"|bc -l)`
    dm3_qc_READS=$(cat ${dm3_path}/${sample}_dm3.rmdup.stat|grep "total (QC-passed reads"|cut -d "+" -f 1)
    QC_reads=$(echo "${dm3_qc_READS}+${Hg19_qc_READS}"|bc )
    dm3_qc_RATIO_intotal=`printf "%.2f\n" $(echo "${dm3_qc_READS}/${ALLREADS}*100"|bc -l)`
    dm3_qc_RATIO_inqc=`printf "%.2f\n" $(echo "${dm3_qc_READS}/${QC_reads}*100"|bc -l)`
    SCALEFACTOR=$(echo "1000000/${dm3_qc_READS}"|bc -l)
    echo -e $sample"\t"$ALLREADS"\t"$Hg19_READS"\t"$Hg19_mapping_RATIO"\t"$Hg19_qc_READS"\t"$Hg19_qc_RATIO"\t"$dm3_qc_READS"\t"$dm3_qc_RATIO_intotal"\t"$dm3_qc_RATIO_inqc"\t"$SCALEFACTOR >> ${logs_dir}/scalefactor.txt
done

echo "scale factor is done"


#step 4.1
### Making CPM-normalized bigWig files with reads adjusted by spike-in ###
mkdir -p ${bw_dir}

cat  ${logs_dir}/scalefactor.txt | sed '1d' | while read id;
do
arr=($id)
sample=${arr[0]}
scalefactor=${arr[9]}
bam_file=${exp_bam_rmdup}/${sample}_hg19.rmdup_Shift.sorted.bam
    if [ ! -s "${bw_dir}/${sample}_rmdup_Shift.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --blackListFileName ~/index/genomes/hg19-blacklist.v2.bed \
        --outFileName ${bw_dir}/${sample}_rmdup_Shift.bw \
        --binSize 1 \
        --scaleFactor $scalefactor \
        --numberOfProcessors 23 \
        --normalizeUsing None
    fi
done

#step 4.2
### Making CPM-normalized bigWig files without spike-in ###
mkdir -p ${bw_dir}
cat ${logs_dir}/scalefactor.txt | sed '1d' | while read id;
do
arr=($id)
sample=${arr[0]}
bam_file=${exp_bam_rmdup}/${sample}_hg19.rmdup_Shift.sorted.bam
    if [ ! -s "${bw_dir}/${sample}_rmdup_Shift_CPM.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --blackListFileName ~/index/genomes/hg19-blacklist.v2.bed \
        --outFileName ${bw_dir}/${sample}_rmdup_Shift_CPM.bw \
        --binSize 1 \
	--scaleFactor 1 \
        --numberOfProcessors 23 \
        --normalizeUsing CPM
    fi
done

#step 5.1
### Call Peaks without input ###
peak_dir=${file_path}/05_peaks
peak_log=${file_path}/logs/peaks
mkdir -p ${peak_dir}/peak_alone
ls ${exp_bam_rmdup}/*_hg19.rmdup_Shift.sorted.bam | grep -v 'nput' | while read id;
do
    sample=`basename $id | sed 's/_hg19.rmdup_Shift.sorted.bam//g'`
    bam_file=${exp_bam_rmdup}/${sample}_hg19.rmdup_Shift.sorted.bam
    blacklist="/share/home/Lemon/index/genomes/hg19-blacklist.v2.bed"
    mkdir -p ${peak_log}/${sample}
    if [ ! -f ${peak_dir}/peak_alone/${sample}_peaks.narrowPeak ]
    then
        macs2 callpeak -t ${bam_file} -f BAMPE -n ${peak_dir}/peak_alone/$sample -g hs --keep-dup all --nolambda 2>${peak_log}/${sample}/${sample}_alone.log
        bedtools intersect -a ${peak_dir}/peak_alone/${sample}_peaks.narrowPeak -b $blacklist -f 0.25 -v > ${peak_dir}/peak_alone/${sample}_peaks.final.narrowPeak &
    fi
done
wait


