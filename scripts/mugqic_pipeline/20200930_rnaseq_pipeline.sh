#!/bin/bash
# Exit immediately on error

set -eu -o pipefail

#-------------------------------------------------------------------------------
# RnaSeq SLURM Job Submission Bash script
# Version: 3.1.5-beta
# Created on: 2020-09-30T21:18:45
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 2 jobs
#   merge_trimmomatic_stats: 1 job
#   star: 6 jobs
#   picard_merge_sam_files: 0 job... skipping
#   picard_sort_sam: 2 jobs
#   picard_mark_duplicates: 2 jobs
#   picard_rna_metrics: 2 jobs
#   estimate_ribosomal_rna: 2 jobs
#   bam_hard_clip: 2 jobs
#   rnaseqc: 2 jobs
#   wiggle: 4 jobs
#   raw_counts: 2 jobs
#   raw_counts_metrics: 4 jobs
#   TOTAL: 31 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/lustre03/project/6001942/chris11/JD_project/output/rna-pipeline-hg38
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/RnaSeq_job_list_$TIMESTAMP
export CONFIG_FILES="/cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.5/pipelines/rnaseq/rnaseq.base.ini,/cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.5/pipelines/rnaseq/rnaseq.beluga.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.UMSCC1_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.UMSCC1_NTsh
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.UMSCC1_NTsh.87828c39c78ae886bf8ca05899a991aa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.UMSCC1_NTsh.87828c39c78ae886bf8ca05899a991aa.mugqic.done' > $COMMAND
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/UMSCC1_NTsh && \
`cat > trim/UMSCC1_NTsh/UMSCC1_NTsh.trim.adapters.fa << END
>Single
CTGTCTCTTATACACATCT
END
` && \
java -XX:ParallelGCThreads=5 -Dsamjdk.buffer_size=1048576 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 5 \
  -phred33 \
  /lustre03/project/6001942/chris11/JD_project/raw/SRR2084807.fastq \
  trim/UMSCC1_NTsh/UMSCC1_NTsh.trim.single.fastq.gz \
  ILLUMINACLIP:trim/UMSCC1_NTsh/UMSCC1_NTsh.trim.adapters.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/UMSCC1_NTsh/UMSCC1_NTsh.trim.log
trimmomatic.UMSCC1_NTsh.87828c39c78ae886bf8ca05899a991aa.mugqic.done
chmod 755 $COMMAND

trimmomatic_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem-per-cpu=2G -N 1 -n 5 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.UMSCC47_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.UMSCC47_NTsh
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.UMSCC47_NTsh.87a55883081e282900a0692e530e237b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.UMSCC47_NTsh.87a55883081e282900a0692e530e237b.mugqic.done' > $COMMAND
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/UMSCC47_NTsh && \
`cat > trim/UMSCC47_NTsh/UMSCC47_NTsh.trim.adapters.fa << END
>Single
CTGTCTCTTATACACATCT
END
` && \
java -XX:ParallelGCThreads=5 -Dsamjdk.buffer_size=1048576 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 5 \
  -phred33 \
  /lustre03/project/6001942/chris11/JD_project/raw/SRR2084809.fastq \
  trim/UMSCC47_NTsh/UMSCC47_NTsh.trim.single.fastq.gz \
  ILLUMINACLIP:trim/UMSCC47_NTsh/UMSCC47_NTsh.trim.adapters.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/UMSCC47_NTsh/UMSCC47_NTsh.trim.log
trimmomatic.UMSCC47_NTsh.87a55883081e282900a0692e530e237b.mugqic.done
chmod 755 $COMMAND

trimmomatic_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem-per-cpu=2G -N 1 -n 5 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.8f1f886ef09663ac4928819921afb342.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'merge_trimmomatic_stats.8f1f886ef09663ac4928819921afb342.mugqic.done' > $COMMAND
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/UMSCC1_NTsh/UMSCC1_NTsh.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/UMSCC1_NTsh	UMSCC1_NTsh	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/UMSCC47_NTsh/UMSCC47_NTsh.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/UMSCC47_NTsh	UMSCC47_NTsh	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"	" '{OFS="	"; if (NR==1) {if ($2=="Raw Paired Reads #") {paired=1};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$2=$2*2; $3=$3*2}; raw[$1]+=$2; surviving[$1]+=$3}}END{for (sample in raw){print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/ && \
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"} else {print $1, $2, sprintf("%\47d", $3), sprintf("%\47d", $4), sprintf("%.1f", $5)}}' metrics/trimReadsetTable.tsv` && \
pandoc \
  /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.5/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --template /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.5/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --variable trailing_min_quality=30 \
  --variable min_length=32 \
  --variable read_type=Single \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
merge_trimmomatic_stats.8f1f886ef09663ac4928819921afb342.mugqic.done
chmod 755 $COMMAND

merge_trimmomatic_stats_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# STEP: star
#-------------------------------------------------------------------------------
STEP=star
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: star_1_JOB_ID: star_align.1.UMSCC1_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.UMSCC1_NTsh
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/star/star_align.1.UMSCC1_NTsh.b8c6081246dbbee8a56a9f231bb2eaa4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'star_align.1.UMSCC1_NTsh.b8c6081246dbbee8a56a9f231bb2eaa4.mugqic.done' > $COMMAND
module load mugqic/star/2.5.3a && \
mkdir -p alignment_1stPass/UMSCC1_NTsh/UMSCC1_NTsh && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/UMSCC1_NTsh/UMSCC1_NTsh.trim.single.fastq.gz \
  --runThreadN 20 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/UMSCC1_NTsh/UMSCC1_NTsh/ \
  --outSAMattrRGline ID:"UMSCC1_NTsh" 	PL:"ILLUMINA" 			SM:"UMSCC1_NTsh" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 120000000000 \
  --limitIObufferSize 4000000000
star_align.1.UMSCC1_NTsh.b8c6081246dbbee8a56a9f231bb2eaa4.mugqic.done
chmod 755 $COMMAND

star_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 20 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: star_2_JOB_ID: star_align.1.UMSCC47_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.UMSCC47_NTsh
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/star/star_align.1.UMSCC47_NTsh.63df80998b304e8f9646d4b446b56112.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'star_align.1.UMSCC47_NTsh.63df80998b304e8f9646d4b446b56112.mugqic.done' > $COMMAND
module load mugqic/star/2.5.3a && \
mkdir -p alignment_1stPass/UMSCC47_NTsh/UMSCC47_NTsh && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/UMSCC47_NTsh/UMSCC47_NTsh.trim.single.fastq.gz \
  --runThreadN 20 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/UMSCC47_NTsh/UMSCC47_NTsh/ \
  --outSAMattrRGline ID:"UMSCC47_NTsh" 	PL:"ILLUMINA" 			SM:"UMSCC47_NTsh" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 120000000000 \
  --limitIObufferSize 4000000000
star_align.1.UMSCC47_NTsh.63df80998b304e8f9646d4b446b56112.mugqic.done
chmod 755 $COMMAND

star_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 20 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: star_3_JOB_ID: star_index.AllSamples
#-------------------------------------------------------------------------------
JOB_NAME=star_index.AllSamples
JOB_DEPENDENCIES=$star_1_JOB_ID:$star_2_JOB_ID
JOB_DONE=job_output/star/star_index.AllSamples.d8b922f40171c8aaa70e4b3799c5a4b6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'star_index.AllSamples.d8b922f40171c8aaa70e4b3799c5a4b6.mugqic.done' > $COMMAND
module load mugqic/star/2.5.3a && \
cat \
  alignment_1stPass/UMSCC1_NTsh/UMSCC1_NTsh/SJ.out.tab \
  alignment_1stPass/UMSCC47_NTsh/UMSCC47_NTsh/SJ.out.tab | \
awk 'BEGIN {OFS="	"; strChar[0]="."; strChar[1]="+"; strChar[2]="-"} {if($5>0){print $1,$2,$3,strChar[$4]}}' | sort -k1,1h -k2,2n > alignment_1stPass/AllSamples.SJ.out.tab && \
mkdir -p reference.Merged && \
STAR --runMode genomeGenerate \
  --genomeDir reference.Merged \
  --genomeFastaFiles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --runThreadN 20 \
  --limitGenomeGenerateRAM 120000000000 \
  --sjdbFileChrStartEnd alignment_1stPass/AllSamples.SJ.out.tab \
  --limitIObufferSize 1000000000 \
  --sjdbOverhang 99
star_index.AllSamples.d8b922f40171c8aaa70e4b3799c5a4b6.mugqic.done
chmod 755 $COMMAND

star_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=15:00:0 --mem=128G -N 1 -n 20 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: star_4_JOB_ID: star_align.2.UMSCC1_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.UMSCC1_NTsh
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$star_3_JOB_ID
JOB_DONE=job_output/star/star_align.2.UMSCC1_NTsh.ec910a52073ef73757e790def6bd3ec6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'star_align.2.UMSCC1_NTsh.ec910a52073ef73757e790def6bd3ec6.mugqic.done' > $COMMAND
module load mugqic/star/2.5.3a && \
mkdir -p alignment/UMSCC1_NTsh/UMSCC1_NTsh && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/UMSCC1_NTsh/UMSCC1_NTsh.trim.single.fastq.gz \
  --runThreadN 20 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/UMSCC1_NTsh/UMSCC1_NTsh/ \
  --outSAMattrRGline ID:"UMSCC1_NTsh" 	PL:"ILLUMINA" 			SM:"UMSCC1_NTsh" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 120000000000 \
  --limitBAMsortRAM 120000000000 \
  --limitIObufferSize 4000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f UMSCC1_NTsh/Aligned.sortedByCoord.out.bam alignment/UMSCC1_NTsh/UMSCC1_NTsh.sorted.bam
star_align.2.UMSCC1_NTsh.ec910a52073ef73757e790def6bd3ec6.mugqic.done
chmod 755 $COMMAND

star_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 20 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: star_5_JOB_ID: star_align.2.UMSCC47_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.UMSCC47_NTsh
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID:$star_3_JOB_ID
JOB_DONE=job_output/star/star_align.2.UMSCC47_NTsh.01f31135ff50e8844f05cabb4a61d7f6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'star_align.2.UMSCC47_NTsh.01f31135ff50e8844f05cabb4a61d7f6.mugqic.done' > $COMMAND
module load mugqic/star/2.5.3a && \
mkdir -p alignment/UMSCC47_NTsh/UMSCC47_NTsh && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/UMSCC47_NTsh/UMSCC47_NTsh.trim.single.fastq.gz \
  --runThreadN 20 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/UMSCC47_NTsh/UMSCC47_NTsh/ \
  --outSAMattrRGline ID:"UMSCC47_NTsh" 	PL:"ILLUMINA" 			SM:"UMSCC47_NTsh" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 120000000000 \
  --limitBAMsortRAM 120000000000 \
  --limitIObufferSize 4000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f UMSCC47_NTsh/Aligned.sortedByCoord.out.bam alignment/UMSCC47_NTsh/UMSCC47_NTsh.sorted.bam
star_align.2.UMSCC47_NTsh.01f31135ff50e8844f05cabb4a61d7f6.mugqic.done
chmod 755 $COMMAND

star_5_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 20 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: star_6_JOB_ID: star_report
#-------------------------------------------------------------------------------
JOB_NAME=star_report
JOB_DEPENDENCIES=$star_4_JOB_ID:$star_5_JOB_ID
JOB_DONE=job_output/star/star_report.88820a8d38f1043dbf0cc863dd909f5f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'star_report.88820a8d38f1043dbf0cc863dd909f5f.mugqic.done' > $COMMAND
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.5/bfx/report/RnaSeq.star.md \
  --variable scientific_name="Homo_sapiens" \
  --variable assembly="GRCh38" \
  /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.5/bfx/report/RnaSeq.star.md \
  > report/RnaSeq.star.md
star_report.88820a8d38f1043dbf0cc863dd909f5f.mugqic.done
chmod 755 $COMMAND

star_6_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# STEP: picard_sort_sam
#-------------------------------------------------------------------------------
STEP=picard_sort_sam
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_1_JOB_ID: picard_sort_sam.UMSCC1_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.UMSCC1_NTsh
JOB_DEPENDENCIES=$star_4_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.UMSCC1_NTsh.f1bd1d73c55b9c7dcd43ba9468c7202f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_sort_sam.UMSCC1_NTsh.f1bd1d73c55b9c7dcd43ba9468c7202f.mugqic.done' > $COMMAND
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.9.0 && \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Xmx40G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/UMSCC1_NTsh/UMSCC1_NTsh.sorted.bam \
 OUTPUT=alignment/UMSCC1_NTsh/UMSCC1_NTsh.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.UMSCC1_NTsh.f1bd1d73c55b9c7dcd43ba9468c7202f.mugqic.done
chmod 755 $COMMAND

picard_sort_sam_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_2_JOB_ID: picard_sort_sam.UMSCC47_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.UMSCC47_NTsh
JOB_DEPENDENCIES=$star_5_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.UMSCC47_NTsh.cdfa2d8f917629cb4b72032d3f6e7aae.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_sort_sam.UMSCC47_NTsh.cdfa2d8f917629cb4b72032d3f6e7aae.mugqic.done' > $COMMAND
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.9.0 && \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Xmx40G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/UMSCC47_NTsh/UMSCC47_NTsh.sorted.bam \
 OUTPUT=alignment/UMSCC47_NTsh/UMSCC47_NTsh.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.UMSCC47_NTsh.cdfa2d8f917629cb4b72032d3f6e7aae.mugqic.done
chmod 755 $COMMAND

picard_sort_sam_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.UMSCC1_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.UMSCC1_NTsh
JOB_DEPENDENCIES=$star_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.UMSCC1_NTsh.f1ba0b4465c4383c135121e1590e10b5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_mark_duplicates.UMSCC1_NTsh.f1ba0b4465c4383c135121e1590e10b5.mugqic.done' > $COMMAND
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.9.0 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/UMSCC1_NTsh/UMSCC1_NTsh.sorted.bam \
 OUTPUT=alignment/UMSCC1_NTsh/UMSCC1_NTsh.sorted.mdup.bam \
 METRICS_FILE=alignment/UMSCC1_NTsh/UMSCC1_NTsh.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000 
picard_mark_duplicates.UMSCC1_NTsh.f1ba0b4465c4383c135121e1590e10b5.mugqic.done
chmod 755 $COMMAND

picard_mark_duplicates_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=20G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.UMSCC47_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.UMSCC47_NTsh
JOB_DEPENDENCIES=$star_5_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.UMSCC47_NTsh.362c33103226a9831fea9936fe707ddb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_mark_duplicates.UMSCC47_NTsh.362c33103226a9831fea9936fe707ddb.mugqic.done' > $COMMAND
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.9.0 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/UMSCC47_NTsh/UMSCC47_NTsh.sorted.bam \
 OUTPUT=alignment/UMSCC47_NTsh/UMSCC47_NTsh.sorted.mdup.bam \
 METRICS_FILE=alignment/UMSCC47_NTsh/UMSCC47_NTsh.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000 
picard_mark_duplicates.UMSCC47_NTsh.362c33103226a9831fea9936fe707ddb.mugqic.done
chmod 755 $COMMAND

picard_mark_duplicates_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=20G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# STEP: picard_rna_metrics
#-------------------------------------------------------------------------------
STEP=picard_rna_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_1_JOB_ID: picard_rna_metrics.UMSCC1_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.UMSCC1_NTsh
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.UMSCC1_NTsh.c95ded989f8e5e46c5aaa013da54b18a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_rna_metrics.UMSCC1_NTsh.c95ded989f8e5e46c5aaa013da54b18a.mugqic.done' > $COMMAND
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.9.0 mugqic/R_Bioconductor/3.5.0_3.7 && \
mkdir -p metrics/UMSCC1_NTsh && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Xmx40G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/UMSCC1_NTsh/UMSCC1_NTsh.sorted.mdup.bam \
 OUTPUT=metrics/UMSCC1_NTsh/UMSCC1_NTsh \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Xmx40G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/UMSCC1_NTsh/UMSCC1_NTsh.sorted.mdup.bam \
 OUTPUT=metrics/UMSCC1_NTsh/UMSCC1_NTsh.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.UMSCC1_NTsh.c95ded989f8e5e46c5aaa013da54b18a.mugqic.done
chmod 755 $COMMAND

picard_rna_metrics_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_rna_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_2_JOB_ID: picard_rna_metrics.UMSCC47_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.UMSCC47_NTsh
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.UMSCC47_NTsh.12524bcd116a7af844ff0fcc157aee93.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_rna_metrics.UMSCC47_NTsh.12524bcd116a7af844ff0fcc157aee93.mugqic.done' > $COMMAND
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.9.0 mugqic/R_Bioconductor/3.5.0_3.7 && \
mkdir -p metrics/UMSCC47_NTsh && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Xmx40G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/UMSCC47_NTsh/UMSCC47_NTsh.sorted.mdup.bam \
 OUTPUT=metrics/UMSCC47_NTsh/UMSCC47_NTsh \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Xmx40G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/UMSCC47_NTsh/UMSCC47_NTsh.sorted.mdup.bam \
 OUTPUT=metrics/UMSCC47_NTsh/UMSCC47_NTsh.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.UMSCC47_NTsh.12524bcd116a7af844ff0fcc157aee93.mugqic.done
chmod 755 $COMMAND

picard_rna_metrics_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_rna_metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# STEP: estimate_ribosomal_rna
#-------------------------------------------------------------------------------
STEP=estimate_ribosomal_rna
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_1_JOB_ID: bwa_mem_rRNA.UMSCC1_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.UMSCC1_NTsh
JOB_DEPENDENCIES=$star_4_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.UMSCC1_NTsh.b857ae349274529846f89616386c1888.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'bwa_mem_rRNA.UMSCC1_NTsh.b857ae349274529846f89616386c1888.mugqic.done' > $COMMAND
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.9.0 mugqic/mugqic_tools/2.2.3 mugqic/python/2.7.13 && \
mkdir -p alignment/UMSCC1_NTsh/UMSCC1_NTsh metrics/UMSCC1_NTsh/UMSCC1_NTsh && \
java -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/UMSCC1_NTsh/UMSCC1_NTsh/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:UMSCC1_NTsh	SM:UMSCC1_NTsh	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1  -Dsamjdk.buffer_size=1048576 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=metrics/UMSCC1_NTsh/UMSCC1_NTsh/UMSCC1_NTshrRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/UMSCC1_NTsh/UMSCC1_NTsh/UMSCC1_NTshrRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/UMSCC1_NTsh/UMSCC1_NTsh/UMSCC1_NTshrRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.UMSCC1_NTsh.b857ae349274529846f89616386c1888.mugqic.done
chmod 755 $COMMAND

estimate_ribosomal_rna_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=64G -N 1 -n 20 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$estimate_ribosomal_rna_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_2_JOB_ID: bwa_mem_rRNA.UMSCC47_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.UMSCC47_NTsh
JOB_DEPENDENCIES=$star_5_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.UMSCC47_NTsh.5cdd475d38b1019d672de630f1a9216c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'bwa_mem_rRNA.UMSCC47_NTsh.5cdd475d38b1019d672de630f1a9216c.mugqic.done' > $COMMAND
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.9.0 mugqic/mugqic_tools/2.2.3 mugqic/python/2.7.13 && \
mkdir -p alignment/UMSCC47_NTsh/UMSCC47_NTsh metrics/UMSCC47_NTsh/UMSCC47_NTsh && \
java -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/UMSCC47_NTsh/UMSCC47_NTsh/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:UMSCC47_NTsh	SM:UMSCC47_NTsh	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1  -Dsamjdk.buffer_size=1048576 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=metrics/UMSCC47_NTsh/UMSCC47_NTsh/UMSCC47_NTshrRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/UMSCC47_NTsh/UMSCC47_NTsh/UMSCC47_NTshrRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/UMSCC47_NTsh/UMSCC47_NTsh/UMSCC47_NTshrRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.UMSCC47_NTsh.5cdd475d38b1019d672de630f1a9216c.mugqic.done
chmod 755 $COMMAND

estimate_ribosomal_rna_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=64G -N 1 -n 20 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$estimate_ribosomal_rna_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# STEP: bam_hard_clip
#-------------------------------------------------------------------------------
STEP=bam_hard_clip
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_1_JOB_ID: tuxedo_hard_clip.UMSCC1_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.UMSCC1_NTsh
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.UMSCC1_NTsh.adac0fe47ec8e9f5613e588c494feb60.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'tuxedo_hard_clip.UMSCC1_NTsh.adac0fe47ec8e9f5613e588c494feb60.mugqic.done' > $COMMAND
module load mugqic/samtools/1.4 && \
samtools view -h \
  alignment/UMSCC1_NTsh/UMSCC1_NTsh.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/UMSCC1_NTsh/UMSCC1_NTsh.sorted.mdup.hardClip.bam
tuxedo_hard_clip.UMSCC1_NTsh.adac0fe47ec8e9f5613e588c494feb60.mugqic.done
chmod 755 $COMMAND

bam_hard_clip_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=32G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bam_hard_clip_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_2_JOB_ID: tuxedo_hard_clip.UMSCC47_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.UMSCC47_NTsh
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.UMSCC47_NTsh.f8bdaf393564a4989fd2bb5def318540.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'tuxedo_hard_clip.UMSCC47_NTsh.f8bdaf393564a4989fd2bb5def318540.mugqic.done' > $COMMAND
module load mugqic/samtools/1.4 && \
samtools view -h \
  alignment/UMSCC47_NTsh/UMSCC47_NTsh.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/UMSCC47_NTsh/UMSCC47_NTsh.sorted.mdup.hardClip.bam
tuxedo_hard_clip.UMSCC47_NTsh.f8bdaf393564a4989fd2bb5def318540.mugqic.done
chmod 755 $COMMAND

bam_hard_clip_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=32G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bam_hard_clip_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# STEP: rnaseqc
#-------------------------------------------------------------------------------
STEP=rnaseqc
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: rnaseqc_1_JOB_ID: rnaseqc
#-------------------------------------------------------------------------------
JOB_NAME=rnaseqc
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/rnaseqc/rnaseqc.6ae9fcd907c26d7725f9f2e18b24d36c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'rnaseqc.6ae9fcd907c26d7725f9f2e18b24d36c.mugqic.done' > $COMMAND
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/bwa/0.7.12 mugqic/rnaseqc/1.1.8 && \
mkdir -p metrics/rnaseqRep && \
echo "Sample	BamFile	Note
UMSCC1_NTsh	alignment/UMSCC1_NTsh/UMSCC1_NTsh.sorted.mdup.bam	RNAseq
UMSCC47_NTsh	alignment/UMSCC47_NTsh/UMSCC47_NTsh.sorted.mdup.bam	RNAseq" \
  > alignment/rnaseqc.samples.txt && \
touch dummy_rRNA.fa && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Xmx40G -jar $RNASEQC_JAR \
  -n 1000 \
  -o metrics/rnaseqRep \
  -r /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  -s alignment/rnaseqc.samples.txt \
  -t /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.transcript_id.gtf \
  -ttype 2 \
  -singleEnd \
  -rRNA /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa && \
zip -r metrics/rnaseqRep.zip metrics/rnaseqRep
rnaseqc.6ae9fcd907c26d7725f9f2e18b24d36c.mugqic.done
chmod 755 $COMMAND

rnaseqc_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=72:00:0 --mem=48G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$rnaseqc_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: rnaseqc_2_JOB_ID: rnaseqc_report
#-------------------------------------------------------------------------------
JOB_NAME=rnaseqc_report
JOB_DEPENDENCIES=$rnaseqc_1_JOB_ID
JOB_DONE=job_output/rnaseqc/rnaseqc_report.c5351a35213b84e9cb9b916fab2d2292.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'rnaseqc_report.c5351a35213b84e9cb9b916fab2d2292.mugqic.done' > $COMMAND
module load mugqic/python/2.7.13 mugqic/pandoc/1.15.2 && \
mkdir -p report && \
cp metrics/rnaseqRep.zip report/reportRNAseqQC.zip && \
python -c 'import csv; csv_in = csv.DictReader(open("metrics/rnaseqRep/metrics.tsv"), delimiter="	")
print "	".join(["Sample", "Aligned Reads", "Alternative Alignments", "%", "rRNA Reads", "Coverage", "Exonic Rate", "Genes"])
print "\n".join(["	".join([
    line["Sample"],
    line["Mapped"],
    line["Alternative Aligments"],
    str(float(line["Alternative Aligments"]) / float(line["Mapped"]) * 100),
    line["rRNA"],
    line["Mean Per Base Cov."],
    line["Exonic Rate"],
    line["Genes Detected"]
]) for line in csv_in])' \
  > report/trimAlignmentTable.tsv.tmp && \
if [[ -f metrics/trimSampleTable.tsv ]]
then
  awk -F"	" 'FNR==NR{raw_reads[$1]=$2; surviving_reads[$1]=$3; surviving_pct[$1]=$4; next}{OFS="	"; if ($2=="Aligned Reads"){surviving_pct[$1]="%"; aligned_pct="%"; rrna_pct="%"} else {aligned_pct=($2 / surviving_reads[$1] * 100); rrna_pct=($5 / surviving_reads[$1] * 100)}; printf $1"	"raw_reads[$1]"	"surviving_reads[$1]"	"surviving_pct[$1]"	"$2"	"aligned_pct"	"$3"	"$4"	"$5"	"rrna_pct; for (i = 6; i<= NF; i++) {printf "	"$i}; print ""}' \
  metrics/trimSampleTable.tsv \
  report/trimAlignmentTable.tsv.tmp \
  > report/trimAlignmentTable.tsv
else
  cp report/trimAlignmentTable.tsv.tmp report/trimAlignmentTable.tsv
fi && \
rm report/trimAlignmentTable.tsv.tmp && \
trim_alignment_table_md=`if [[ -f metrics/trimSampleTable.tsv ]] ; then cut -f1-13 report/trimAlignmentTable.tsv | LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:"} else {print $1, sprintf("%\47d", $2), sprintf("%\47d", $3), sprintf("%.1f", $4), sprintf("%\47d", $5), sprintf("%.1f", $6), sprintf("%\47d", $7), sprintf("%.1f", $8), sprintf("%\47d", $9), sprintf("%.1f", $10), sprintf("%.2f", $11), sprintf("%.2f", $12), sprintf("%\47d", $13)}}' ; else cat report/trimAlignmentTable.tsv | LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:"} else {print $1, sprintf("%\47d", $2), sprintf("%\47d", $3), sprintf("%.1f", $4), sprintf("%\47d", $5), sprintf("%.2f", $6), sprintf("%.2f", $7), $8}}' ; fi`
pandoc \
  /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.5/bfx/report/RnaSeq.rnaseqc.md \
  --template /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.5/bfx/report/RnaSeq.rnaseqc.md \
  --variable trim_alignment_table="$trim_alignment_table_md" \
  --to markdown \
  > report/RnaSeq.rnaseqc.md
rnaseqc_report.c5351a35213b84e9cb9b916fab2d2292.mugqic.done
chmod 755 $COMMAND

rnaseqc_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$rnaseqc_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# STEP: wiggle
#-------------------------------------------------------------------------------
STEP=wiggle
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: wiggle_1_JOB_ID: bed_graph.UMSCC1_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=bed_graph.UMSCC1_NTsh
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/wiggle/bed_graph.UMSCC1_NTsh.1e5d37923409f6e76e41d335cebd2af3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'bed_graph.UMSCC1_NTsh.1e5d37923409f6e76e41d335cebd2af3.mugqic.done' > $COMMAND
module load mugqic/samtools/1.4 mugqic/bedtools/2.26.0 && \
mkdir -p tracks/UMSCC1_NTsh  && \
nmblines=$(samtools view -F 256 alignment/UMSCC1_NTsh/UMSCC1_NTsh.sorted.mdup.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed  -bg -split -scale $scalefactor \
  -ibam alignment/UMSCC1_NTsh/UMSCC1_NTsh.sorted.mdup.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/UMSCC1_NTsh/UMSCC1_NTsh.bedGraph
bed_graph.UMSCC1_NTsh.1e5d37923409f6e76e41d335cebd2af3.mugqic.done
chmod 755 $COMMAND

wiggle_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=38G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: wiggle_2_JOB_ID: wiggle.UMSCC1_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.UMSCC1_NTsh
JOB_DEPENDENCIES=$wiggle_1_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.UMSCC1_NTsh.f529b629cdbbaf92e64c9049bde6e817.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'wiggle.UMSCC1_NTsh.f529b629cdbbaf92e64c9049bde6e817.mugqic.done' > $COMMAND
module load mugqic/ucsc/v346 && \
mkdir -p tracks/bigWig && \
cat tracks/UMSCC1_NTsh/UMSCC1_NTsh.bedGraph | sort --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -v "GL\|lambda\|pUC19\|KI\|\KN\|random"  | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/UMSCC1_NTsh/UMSCC1_NTsh.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/UMSCC1_NTsh/UMSCC1_NTsh.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/UMSCC1_NTsh.bw
wiggle.UMSCC1_NTsh.f529b629cdbbaf92e64c9049bde6e817.mugqic.done
chmod 755 $COMMAND

wiggle_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem-per-cpu=4775M -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: wiggle_3_JOB_ID: bed_graph.UMSCC47_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=bed_graph.UMSCC47_NTsh
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/wiggle/bed_graph.UMSCC47_NTsh.dd04423bc025eb8b442dfba9d7e26d02.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'bed_graph.UMSCC47_NTsh.dd04423bc025eb8b442dfba9d7e26d02.mugqic.done' > $COMMAND
module load mugqic/samtools/1.4 mugqic/bedtools/2.26.0 && \
mkdir -p tracks/UMSCC47_NTsh  && \
nmblines=$(samtools view -F 256 alignment/UMSCC47_NTsh/UMSCC47_NTsh.sorted.mdup.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed  -bg -split -scale $scalefactor \
  -ibam alignment/UMSCC47_NTsh/UMSCC47_NTsh.sorted.mdup.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/UMSCC47_NTsh/UMSCC47_NTsh.bedGraph
bed_graph.UMSCC47_NTsh.dd04423bc025eb8b442dfba9d7e26d02.mugqic.done
chmod 755 $COMMAND

wiggle_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=38G -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: wiggle_4_JOB_ID: wiggle.UMSCC47_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.UMSCC47_NTsh
JOB_DEPENDENCIES=$wiggle_3_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.UMSCC47_NTsh.f1c4462bd398cf327935f336aa741ab4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'wiggle.UMSCC47_NTsh.f1c4462bd398cf327935f336aa741ab4.mugqic.done' > $COMMAND
module load mugqic/ucsc/v346 && \
mkdir -p tracks/bigWig && \
cat tracks/UMSCC47_NTsh/UMSCC47_NTsh.bedGraph | sort --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -v "GL\|lambda\|pUC19\|KI\|\KN\|random"  | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/UMSCC47_NTsh/UMSCC47_NTsh.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/UMSCC47_NTsh/UMSCC47_NTsh.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/UMSCC47_NTsh.bw
wiggle.UMSCC47_NTsh.f1c4462bd398cf327935f336aa741ab4.mugqic.done
chmod 755 $COMMAND

wiggle_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem-per-cpu=4775M -N 1 -n 10 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# STEP: raw_counts
#-------------------------------------------------------------------------------
STEP=raw_counts
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: raw_counts_1_JOB_ID: htseq_count.UMSCC1_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.UMSCC1_NTsh
JOB_DEPENDENCIES=$picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.UMSCC1_NTsh.d4284a11a7c46475aa58c85920b1bf73.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'htseq_count.UMSCC1_NTsh.d4284a11a7c46475aa58c85920b1bf73.mugqic.done' > $COMMAND
module load mugqic/samtools/1.4 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/UMSCC1_NTsh/UMSCC1_NTsh.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/UMSCC1_NTsh.readcounts.csv
htseq_count.UMSCC1_NTsh.d4284a11a7c46475aa58c85920b1bf73.mugqic.done
chmod 755 $COMMAND

raw_counts_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: raw_counts_2_JOB_ID: htseq_count.UMSCC47_NTsh
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.UMSCC47_NTsh
JOB_DEPENDENCIES=$picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.UMSCC47_NTsh.570ba6eb75ea7c09c7c8ac94197b04aa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'htseq_count.UMSCC47_NTsh.570ba6eb75ea7c09c7c8ac94197b04aa.mugqic.done' > $COMMAND
module load mugqic/samtools/1.4 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/UMSCC47_NTsh/UMSCC47_NTsh.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/UMSCC47_NTsh.readcounts.csv
htseq_count.UMSCC47_NTsh.570ba6eb75ea7c09c7c8ac94197b04aa.mugqic.done
chmod 755 $COMMAND

raw_counts_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# STEP: raw_counts_metrics
#-------------------------------------------------------------------------------
STEP=raw_counts_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_1_JOB_ID: metrics.matrix
#-------------------------------------------------------------------------------
JOB_NAME=metrics.matrix
JOB_DEPENDENCIES=$raw_counts_1_JOB_ID:$raw_counts_2_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/metrics.matrix.60c335909c5a5adb7af8b422868bc975.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'metrics.matrix.60c335909c5a5adb7af8b422868bc975.mugqic.done' > $COMMAND
module load mugqic/mugqic_tools/2.2.3 && \
mkdir -p DGE && \
gtf2tmpMatrix.awk \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  DGE/tmpMatrix.txt && \
HEAD='Gene\tSymbol' && \
for read_count_file in \
  raw_counts/UMSCC1_NTsh.readcounts.csv \
  raw_counts/UMSCC47_NTsh.readcounts.csv
do
  sort -k1,1 $read_count_file > DGE/tmpSort.txt && \
  join -1 1 -2 1 <(sort -k1,1 DGE/tmpMatrix.txt) DGE/tmpSort.txt > DGE/tmpMatrix.2.txt && \
  mv DGE/tmpMatrix.2.txt DGE/tmpMatrix.txt && \
  na=$(basename $read_count_file | rev | cut -d. -f3- | rev) && \
  HEAD="$HEAD\t$na"
done && \
echo -e $HEAD | cat - DGE/tmpMatrix.txt | tr ' ' '\t' > DGE/rawCountMatrix.csv && \
rm DGE/tmpSort.txt DGE/tmpMatrix.txt
metrics.matrix.60c335909c5a5adb7af8b422868bc975.mugqic.done
chmod 755 $COMMAND

raw_counts_metrics_1_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=5:00:0 --mem=4G -N 1 -n 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_2_JOB_ID: metrics.wigzip
#-------------------------------------------------------------------------------
JOB_NAME=metrics.wigzip
JOB_DEPENDENCIES=$wiggle_2_JOB_ID:$wiggle_4_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/metrics.wigzip.edc4e268c60b94d072db60163e113f9c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'metrics.wigzip.edc4e268c60b94d072db60163e113f9c.mugqic.done' > $COMMAND
zip -r tracks.zip tracks/bigWig
metrics.wigzip.edc4e268c60b94d072db60163e113f9c.mugqic.done
chmod 755 $COMMAND

raw_counts_metrics_2_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=5:00:0 --mem=4G -N 1 -n 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_3_JOB_ID: rpkm_saturation
#-------------------------------------------------------------------------------
JOB_NAME=rpkm_saturation
JOB_DEPENDENCIES=$raw_counts_metrics_1_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/rpkm_saturation.6527cb0fd54825a428fcb6f6ab1785f8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'rpkm_saturation.6527cb0fd54825a428fcb6f6ab1785f8.mugqic.done' > $COMMAND
module load mugqic/R_Bioconductor/3.5.0_3.7 mugqic/mugqic_tools/2.2.3 && \
mkdir -p metrics/saturation && \
Rscript $R_TOOLS/rpkmSaturation.R \
  DGE/rawCountMatrix.csv \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.genes.length.tsv \
  raw_counts \
  metrics/saturation \
  20 \
  1 && \
zip -r metrics/saturation.zip metrics/saturation
rpkm_saturation.6527cb0fd54825a428fcb6f6ab1785f8.mugqic.done
chmod 755 $COMMAND

raw_counts_metrics_3_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=90G -N 1 -n 20 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_metrics_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_4_JOB_ID: raw_count_metrics_report
#-------------------------------------------------------------------------------
JOB_NAME=raw_count_metrics_report
JOB_DEPENDENCIES=$rnaseqc_1_JOB_ID:$raw_counts_metrics_2_JOB_ID:$raw_counts_metrics_3_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/raw_count_metrics_report.bc9f6eaedc384899418daa2984235181.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'raw_count_metrics_report.bc9f6eaedc384899418daa2984235181.mugqic.done' > $COMMAND
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
cp metrics/rnaseqRep/corrMatrixSpearman.txt report/corrMatrixSpearman.tsv && \
cp tracks.zip report/ && \
cp metrics/saturation.zip report/ && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.5/bfx/report/RnaSeq.raw_counts_metrics.md \
  --variable corr_matrix_spearman_table="`head -16 report/corrMatrixSpearman.tsv | cut -f-16| awk -F"	" '{OFS="	"; if (NR==1) {$0="Vs"$0; print; gsub(/[^	]/, "-"); print} else {printf $1; for (i=2; i<=NF; i++) {printf "	"sprintf("%.2f", $i)}; print ""}}' | sed 's/	/|/g'`" \
  /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.5/bfx/report/RnaSeq.raw_counts_metrics.md \
  > report/RnaSeq.raw_counts_metrics.md
raw_count_metrics_report.bc9f6eaedc384899418daa2984235181.mugqic.done
chmod 755 $COMMAND

raw_counts_metrics_4_JOB_ID=$(echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE &&    $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_metrics_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.1


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
LOG_MD5=$(echo $USER-'10.74.73.2-RnaSeq-UMSCC1_NTsh.UMSCC1_NTsh,UMSCC47_NTsh.UMSCC47_NTsh' | md5sum | awk '{ print $1 }')
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=beluga2.int.ets1.calculquebec.ca&ip=10.74.73.2&pipeline=RnaSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,star,picard_merge_sam_files,picard_sort_sam,picard_mark_duplicates,picard_rna_metrics,estimate_ribosomal_rna,bam_hard_clip,rnaseqc,wiggle,raw_counts,raw_counts_metrics&samples=2&md5=$LOG_MD5" --quiet --output-document=/dev/null

