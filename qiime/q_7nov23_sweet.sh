#!bin/bash

#####################
# Variables to change 
#####################
#################################################################################
# Define in and out folder (have to change)
#----------------------------------------


in_folder=./data

out_folder=./out

q_env_path=./qiime2

# Optional to change (parameters have defualt values)
# General parameters 
n_threads=60

# DADA2 
p_trim_left_f=17
p_trim_left_r=21
p_trunc_len_f=245
p_trunc_len_r=240
p_max_ee_f=2
p_max_ee_r=2

# Taxanomic classification (path to classifier)
tax_classif=./silva-138-99-nb-weighted-classifier.qza


##################################################################################
##################################################################################


# Make output directory 
mkdir -p $out_folder 

mkdir -p $out_folder/summary

mkdir -p $out_folder/export

# Activate conda env 
source activate $q_env_path


# Import files in as artefacts 
#---------------------------------------------------------
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $in_folder \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path $out_folder/demux-paired-end.qza
  
# Make a summary of demultiplex files 
#---------------------------------------------------------
qiime demux summarize \
  --i-data $out_folder/demux-paired-end.qza \
  --o-visualization $out_folder/summary/reads_summary.qzv
  
# Pick ASV using DADA2 algorithm 
#---------------------------------------------------------------------------------------------------------- 
qiime dada2 denoise-paired \
--i-demultiplexed-seqs $out_folder/demux-paired-end.qza \
--p-trim-left-f $p_trim_left_f \
--p-trim-left-r $p_trim_left_r \
--p-trunc-len-f $p_trunc_len_f \
--p-trunc-len-r $p_trunc_len_r \
--p-max-ee-f $p_max_ee_f \
--p-max-ee-r $p_max_ee_r \
--p-n-threads $n_threads \
--o-table $out_folder/asv_table.qza \
--o-representative-sequences $out_folder/rep_seqs.qza \
--o-denoising-stats $out_folder/summary/denoising_stats.qza

# Visualize statistics
qiime metadata tabulate \
--m-input-file $out_folder/summary/denoising_stats.qza \
--o-visualization $out_folder/summary/denoising_stats.qzv

# Make phylogenetic tree. 
# ---------------------------------------------------------
mkdir -p $out_folder/tree/

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences $out_folder/rep_seqs.qza \
--p-n-threads $n_threads \
--o-alignment $out_folder/tree/aligned-rep-seqs.qza \
--o-masked-alignment $out_folder/tree/masked-aligned-rep-seqs.qza \
--o-tree $out_folder/tree/unrooted-tree.qza \
--o-rooted-tree $out_folder/tree/rooted-tree.qza

# Assign taxanomy 
#--------------------------------------------------------
qiime feature-classifier classify-sklearn \
--i-classifier $tax_classif \
--i-reads $out_folder/rep_seqs.qza \
--p-n-jobs $n_threads \
--p-confidence 0.7 \
--o-classification $out_folder/taxonomy_07.qza
  
qiime metadata tabulate \
--m-input-file $out_folder/taxonomy_07.qza \
--o-visualization $out_folder/summary/taxonomy_07.qzv

conda deactivate 


# Save results as an archive 
rm $out_folder/demux-paired-end.qza

cd $out_folder

tar -cvf all_qiime.tar.gz *

cd - 

mv $out_folder/all_qiime.tar.gz $out_folder.tar.gz
