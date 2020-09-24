export RAP_ID=def-stbil30
export JOB_MAIL="christophe.tav@gmail.com"

mkdir -p $HOME/projects/def-stbil30/chris11/JD_project/output/rna-pipeline-hg38

rnaseq.py --job-scheduler slurm -s 1-14 \
  --log debug \
  --readsets raw/readset_20200923.txt \
  -o output/rna-pipeline-hg38 \
  --no-json \
  --config $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.base.ini \
      $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.beluga.ini \
      $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini
