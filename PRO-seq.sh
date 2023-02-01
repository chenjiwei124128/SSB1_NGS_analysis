###file_path
input_path=$1
file_path=`pwd`
raw_dir=${file_path}/00_rawdata
logs_dir=${file_path}/logs
fastqc_dir=${file_path}/01_fastqc
trimmedFastq_dir=${file_path}/02_trimmeddata
trimmedFastq_log_dir=${file_path}/logs/trimmeddata
rmrRNAdata_dir=${file_path}/02_rmrRNAdata
rmrRNAdata_log_dir=${file_path}/logs/rRNA
rmrRNAdata_fastqc_dir=${file_path}/02_rmrRNAdata_fastqc
align_exp_dir=${file_path}/03_bam_hg19
alignexp_log_dir=${file_path}/logs/align_hg19
align_spike_dir=${file_path}/03_spikebam_dm3
alignspike_log_dir=${file_path}/logs/align_dm3
exp_bam_rmdup=${file_path}/03_bam_hg19_rmdup
rmdup_exp_log=${file_path}/logs/rmdup_state_hg19
spike_bam_rmdup=${file_path}/03_spikebam_dm3_rmdup
rmdup_spike_log=${file_path}/logs/rmdup_state_dm3
bw_fulllength_dir=${file_path}/04_bw_fulllength_rmdup
bw_singlebase_dir=${file_path}/04_bw_singlebase_rmdup
plot_dir=${file_path}/05_analysis
sampleinfo=${file_path}/sample_info.txt

#reference genome
GENOME_EXP="/share/home/Lemon/index/bowtie2-hg19/hg19"
GENOME_SPIKE="/share/home/Lemon/index/bowtie2-dm3/dm3"
RDNA="/share/home/Lemon/index/bowtie2-human_rDNA/human_rDNA"
SPIKE_PREFIX="dm3"

MAPQ=10


#step 0
####check files 

#step 1.1
####change file name#####
mkdir -p ${raw_dir}
ln -s `ls ${input_path}/*.gz` ${raw_dir}/ >/dev/null 2>&1

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


#step 1.2
####fastqc of raw data ####
mkdir -p ${fastqc_dir}

in_path=${raw_dir}
out_path=${fastqc_dir}

for FILE in `ls ${in_path}/*.gz`
do 
	if [ ! -s ${out_path}/"$(basename ${FILE/.fastq.gz/_fastqc.zip})" ]
	then
	   fastqc $FILE -t 2 -o ${out_path}/ &
	fi
done
wait

#step 1.3
#check fq
cd ${fastqc_dir}
ls *fastqc.zip|sed 's/_R[1-2]_fastqc.zip//g'|sort|uniq -c

### merge reports of fastqc
~/software/anaconda3/envs/py36/bin/multiqc ${fastqc_dir}/ -n rawdata_multiqc -o ${fastqc_dir}/


#step 2.1
### Trimming adapters  (trim_galore)
####remember to change --length for proseq !!!!!! #####

mkdir -p ${trimmedFastq_dir}
mkdir -p ${trimmedFastq_log_dir}

for fq1 in `ls ${raw_dir}/*R1.fastq.gz`
do
fq2=${fq1/R1.fastq.gz/R2.fastq.gz}
    if [ ! -s ${trimmedFastq_dir}/"$(basename ${fq1/.fastq.gz/_val_1.fq.gz})" ]
    then
    ~/software/anaconda3/envs/py36/bin/trim_galore -j 10 --length 15 --paired -o ${trimmedFastq_dir} $fq1 $fq2 --path_to_cutadapt ~/software/anaconda3/envs/py36/bin/cutadapt \
    > ${trimmedFastq_log_dir}/"$(basename ${fq1/_R1.fq.gz/_trimmed.log})" 2>&1 & 
    fi
done
wait

#step 2.2
####filtering rRNA reads#####
mkdir -p ${rmrRNAdata_dir}
mkdir -p ${rmrRNAdata_log_dir}

for fq1 in `ls ${trimmedFastq_dir}/*_1.fq.gz`
do
fq2=${fq1/R1_val_1.fq.gz/R2_val_2.fq.gz}
sample=$(basename ${fq1/_R1_val_1.fq.gz/})
bowtie2 \
        --fast-local \
        -x ${RDNA} \
        -1 $fq1 \
        -2 $fq2 \
        --un-conc-gz  ${rmrRNAdata_dir}/${sample}_rmrRNA.fq.gz \
        --interleaved - \
        --threads 23 2> ${rmrRNAdata_log_dir}/${sample}_rRNA_bowtie.log > /dev/null
done

#step 2.3
#####rename rmrRNAdata####
for FILE in ${rmrRNAdata_dir}/*.fq.1.gz
do
    if [ ! -s ${FILE/_rmrRNA.fq.1.gz/_R1.fq.gz} ]
    then
        mv "$FILE" ${FILE/_rmrRNA.fq.1.gz/_R1.fq.gz}
    fi
done

for FILE in ${rmrRNAdata_dir}/*.fq.2.gz
do
    if [ ! -s ${FILE/_rmrRNA.fq.2.gz/_R2.fq.gz} ]
    then
        mv "$FILE" ${FILE/_rmrRNA.fq.2.gz/_R2.fq.gz}
    fi
done


#step 2.4
#####qc for trimmed data####
mkdir -p ${rmrRNAdata_fastqc_dir}

for FILE in `ls ${rmrRNAdata_dir}/*fq.gz`
do
if [ ! -s ${rmrRNAdata_fastqc_dir}/"$(basename ${FILE/.fq.gz/_fastqc.html})" ]
then
 fastqc $FILE -t 2 -o ${rmrRNAdata_fastqc_dir} &
fi
done
wait
### merge reports of fastqc
~/software/anaconda3/envs/py36/bin/multiqc ${rmrRNAdata_fastqc_dir}/ -n rmrRNA_multiqc -o ${rmrRNAdata_fastqc_dir}/


#step 3.1
###Aligning to experimental genome#####
mkdir ${align_exp_dir}
mkdir ${alignexp_log_dir}

for PAIR in $(ls ${rmrRNAdata_dir} | sed 's/_R[1-2].fq.gz//' | uniq )
do
if [ ! -s "${align_exp_dir}/${PAIR}_hg19.bam" ]
then
    echo "aligning ${PAIR} to experimental genome"
    (bowtie2 \
    --local \
    --sensitive-local \
    --threads 25 \
    -x "$GENOME_EXP" \
    -1 "${rmrRNAdata_dir}/${PAIR}_R1.fq.gz" \
    -2 "${rmrRNAdata_dir}/${PAIR}_R2.fq.gz" \
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

for PAIR in $(ls ${rmrRNAdata_dir} | sed 's/_R[1-2].fq.gz//' | uniq )
do
if [ ! -s "${align_spike_dir}/${PAIR}_onlydm3.bam" ]
    then
    echo "aligning ${PAIR} to spike-in genome"
    (bowtie2 \
    --local \
    --very-sensitive-local \
    --threads 25 \
    --no-unal \
    --no-mixed \
    --no-discordant \
    -x "$GENOME_SPIKE" \
    -1 "${rmrRNAdata_dir}/${PAIR}_R1.fq.gz" \
    -2 "${rmrRNAdata_dir}/${PAIR}_R2.fq.gz" \
    2> ${alignspike_log_dir}/${PAIR}_onlydm3Align.log) |
    samtools view -bS -f 2 -q ${MAPQ} |
    samtools sort -@ 20 -o ${align_spike_dir}/${PAIR}_onlydm3.bam
    samtools index ${align_spike_dir}/${PAIR}_onlydm3.bam
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
dm3_path=${rmdup_spike_log}
align_path=${alignexp_log_dir}

echo -e "sample\tALLREADS\tHG19_READS\tHG19_mapping_RATIO\tHG19_qc_READS\tHG19_qc_RATIO\tDM3_qc_READS\tDM3_qc_RATIO_intotal\tDM3_qc_RATIO_inqc\tSCALEFACTOR" >${logs_dir}/scalefactor.txt
ls $align_path|while read file;
do
    sample=${file/_align.log/}
    ALLREADS=$(cat ${align_path}/$file|grep "were paired; of these:$"|cut -d "(" -f 1|awk '{print $1*2}')
    HG19_READS=$(cat ${align_path}/$file| sed 's/%//g' | awk '{printf $0"\t"}'  |cut -f 4,5,8,13,14 |
 sed 's/\t/\n/g' | awk '{print $1}' | awk '{printf $0"\t"}'|awk '{print 2*($1+$2+$3)+$4+$5}')
    HG19_mapping_RATIO=$(cat ${align_path}/$file|grep "overall alignment rate"|cut -d "%" -f 1)
    HG19_qc_READS=$(cat $hg19_path/${sample}_hg19.rmdup.stat|grep "total (QC-passed reads"|cut -d " " -f 1)
    HG19_qc_RATIO=`printf "%.2f\n" $(echo "100*${HG19_qc_READS}/${ALLREADS}"|bc -l)`
    DM3_qc_READS=$(cat ${dm3_path}/${sample}_onlydm3.rmdup.stat|grep "total (QC-passed reads"|cut -d "+" -f 1)
    QC_reads=$(echo "${DM3_qc_READS}+${HG19_qc_READS}"|bc )
    DM3_qc_RATIO_intotal=`printf "%.2f\n" $(echo "100*${DM3_qc_READS}/${ALLREADS}"|bc -l)`
    DM3_qc_RATIO_inqc=`printf "%.2f\n" $(echo "100*${DM3_qc_READS}/${QC_reads}"|bc -l)`
    SCALEFACTOR=$(echo "1000000/${DM3_qc_READS}"|bc -l)
    echo -e $sample"\t"$ALLREADS"\t"$HG19_READS"\t"$HG19_mapping_RATIO"\t"$HG19_qc_READS"\t"$HG19_qc_RATIO"\t"$DM3_qc_READS"\t"$DM3_qc_RATIO_intotal"\t"$DM3_qc_RATIO_inqc"\t"$SCALEFACTOR >> ${logs_dir}/scalefactor.txt
done



#step 4.1
### Making CPM-normalized bigWig files with full-length reads adjusted by spike-in ###
mkdir -p ${bw_fulllength_dir}

cat  ${logs_dir}/scalefactor.txt | sed '1d' | while read id;
do
arr=($id)
sample=${arr[0]}
scalefactor=${arr[9]}
bam_file=${exp_bam_rmdup}/${sample}_hg19.rmdup.bam
    if [ ! -s "${bw_fulllength_dir}/${sample}_fulllength_spike_fwd.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_fulllength_dir}/${sample}_fulllength_spike_fwd.bw \
        --binSize 1 \
        --scaleFactor  $scalefactor \
        --numberOfProcessors 23 \
        --normalizeUsing None \
        --samFlagInclude 82
    fi
    if [ ! -s "${bw_fulllength_dir}/${sample}_fulllength_spike_rev.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_fulllength_dir}/${sample}_fulllength_spike_rev.bw \
        --binSize 1 \
        --scaleFactor  $scalefactor \
        --numberOfProcessors 23 \
        --normalizeUsing None \
        --samFlagInclude 98
    fi
    if [ ! -s "${bw_fulllength_dir}/${sample}_fulllength_spike_rev_minus.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_fulllength_dir}/${sample}_fulllength_spike_rev_minus.bw \
        --binSize 1 \
        --scaleFactor  -"${scalefactor}" \
        --numberOfProcessors 23 \
        --normalizeUsing None \
        --samFlagInclude 98
    fi
done


#step 4.2
### Making CPM-normalized bigWig files with single base adjusted by spike-in for forward and reverse strands###
mkdir -p ${bw_singlebase_dir}
cat  ${logs_dir}/scalefactor.txt | sed '1d' | while read id;
do
arr=($id)
sample=${arr[0]}
scalefactor=${arr[9]}
bam_file=${exp_bam_rmdup}/${sample}_hg19.rmdup.bam
    if [ ! -s "${bw_singlebase_dir}/${sample}_singlebase_spike_fwd.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_singlebase_dir}/${sample}_singlebase_spike_fwd.bw \
        --binSize 1 \
        --scaleFactor  $scalefactor \
        --numberOfProcessors 23 \
        --normalizeUsing None \
        --Offset 1 \
        --samFlagInclude 82
    fi
    if [ ! -s "${bw_singlebase_dir}/${sample}_singlebase_spike_rev.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_singlebase_dir}/${sample}_singlebase_spike_rev.bw \
        --binSize 1 \
        --scaleFactor  $scalefactor \
        --numberOfProcessors 23 \
        --normalizeUsing None \
        --Offset 1 \
        --samFlagInclude 98
    fi
    if [ ! -s "${bw_singlebase_dir}/${sample}_singlebase_spike_rev_minus.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_singlebase_dir}/${sample}_singlebase_spike_rev_minus.bw \
        --binSize 1 \
        --scaleFactor  -"${scalefactor}" \
        --numberOfProcessors 23 \
        --normalizeUsing None \
        --Offset 1 \
        --samFlagInclude 98
    fi
done

#step 5.1
### Making CPM-normalized bigWig files with full-length reads without spike-in ###
mkdir -p ${bw_fulllength_dir}


cat  ${logs_dir}/scalefactor.txt | sed '1d' | while read id;
do
arr=($id)
sample=${arr[0]}
scalefactor=${arr[9]}
bam_file=${exp_bam_rmdup}/${sample}_hg19.rmdup.bam
    if [ ! -s "${bw_fulllength_dir}/${sample}_fulllength_CPM_fwd.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_fulllength_dir}/${sample}_fulllength_CPM_fwd.bw \
        --binSize 1 \
		--scaleFactor 1	\
        --numberOfProcessors 23 \
        --normalizeUsing CPM \
        --samFlagInclude 82
    fi
    if [ ! -s "${bw_fulllength_dir}/${sample}_fulllength_CPM_rev.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_fulllength_dir}/${sample}_fulllength_CPM_rev.bw \
        --binSize 1 \
		--scaleFactor 1	\
        --numberOfProcessors 23 \
        --normalizeUsing CPM \
        --samFlagInclude 98
    fi
    if [ ! -s "${bw_fulllength_dir}/${sample}_fulllength_CPM_rev_minus.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_fulllength_dir}/${sample}_fulllength_CPM_rev_minus.bw \
        --binSize 1 \
        --scaleFactor -1 \
        --numberOfProcessors 23 \
        --normalizeUsing CPM \
        --samFlagInclude 98
    fi
done


#step 5.2
### Making CPM-normalized bigWig files with single base without spike-in for forward and reverse strands###
mkdir -p ${bw_singlebase_dir}
cat  ${logs_dir}/scalefactor.txt | sed '1d' | while read id;
do
arr=($id)
sample=${arr[0]}
scalefactor=${arr[9]}
bam_file=${exp_bam_rmdup}/${sample}_hg19.rmdup.bam
    if [ ! -s "${bw_singlebase_dir}/${sample}_singlebase_CPM_fwd.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_singlebase_dir}/${sample}_singlebase_CPM_fwd.bw \
        --binSize 1 \
        --numberOfProcessors 23 \
        --scaleFactor 1 \
        --normalizeUsing CPM \
        --Offset 1 \
        --samFlagInclude 82
    fi
    if [ ! -s "${bw_singlebase_dir}/${sample}_singlebase_CPM_rev.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_singlebase_dir}/${sample}_singlebase_CPM_rev.bw \
        --binSize 1 \
        --numberOfProcessors 23 \
        --scaleFactor 1 \
        --normalizeUsing CPM \
        --Offset 1 \
        --samFlagInclude 98
    fi
    if [ ! -s "${bw_singlebase_dir}/${sample}_singlebase_CPM_rev_minus.bw" ]
    then
        ~/software/anaconda3/envs/py36/bin/bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${bw_singlebase_dir}/${sample}_singlebase_CPM_rev_minus.bw \
        --binSize 1 \
        --scaleFactor -1 \
        --numberOfProcessors 23 \
        --normalizeUsing CPM \
        --Offset 1 \
        --samFlagInclude 98
    fi
done




