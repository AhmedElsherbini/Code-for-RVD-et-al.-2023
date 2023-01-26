#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Author: Ahmed Elsherbini
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
####################################
conda activate qiime2-2021.11

###################################
#download the data the fna and tax.txv of greengenes or silva from "https://docs.qiime2.org/2021.4/data-resources/"
#convert the fasta into .gza

qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path HOMD_16S_rRNA_RefSeq_V15.22.p9.fasta \
--output-path HOMD_data.qza
###################################

#convert the tax.tsv into artifact

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path HOMD_16S_rRNA_RefSeq_V15.22.qiime.tsv \
--output-path HOMD_tax.qza

###################################
#extract read using the primers you can see them in the (Escapa et al., 2020)
#https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00841-w

qiime feature-classifier extract-reads \
--i-sequences HOMD_data.qza \
--p-f-primer CAATTACCGCGGCTGCTGG \
--p-r-primer CCGAGTTTGATCMTGGCTCAG \
--o-reads HOMD_reads.qza 
###################################
#train using the classfier

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads HOMD_reads.qza  \
--i-reference-taxonomy HOMD_tax.qza \
--o-classifier HOMD_classfier.qza 

####################################

# IMPORT REP-SEQS FROM R IN QIIME2

qiime tools import \
--input-path rep-seqs.fna \
--type 'FeatureData[Sequence]' \
--output-path rep-seqs.qza
###################################
# CONVERT AND IMPORT SEQTAB-NOCHIM FROM R AS FEATURE-TABLE IN QIIME2
echo -n "#OTU Table" | cat - seqtab.txt > biom-table.txt
biom convert -i biom-table.txt -o feature-table.biom --table-type="OTU table" --to-hdf5
###################################
qiime tools import --input-path feature-table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path feature-table.qza
###################################
# TAXONOMIC CLASSIFICATION

qiime feature-classifier classify-sklearn \
--i-classifier HOMD_classfier.qza \
--i-reads rep-seqs.qza \
--o-classification taxonomy.qza
###################################

qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv
###################################

qiime taxa barplot \
  --i-table feature-table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file new_sample_metadata.tsv \
  --o-visualization h_taxa-bar-plots.qzv

###############################################
#let's do a phylogentic trees######### 

qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza

################################masking the tree 
##Masking helps to eliminate alignment columns that are phylogenetically uninformative or misleading before phylogenetic analysis#

qiime alignment mask \
--i-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza

#########################darwing the unrooted tree

qiime phylogeny fasttree \
--i-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza
######################### we need to root it !

qiime phylogeny  midpoint-root \
--i-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

############to determine the depth of the sequnceing#################################################### 
#first we should summarize our features tables

qiime feature-table summarize \
--i-table feature-table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file sample_metadata.tsv 
############to see the table###########################################
qiime tools view table.qzv
###############to choose the sampling depth#####!!!!!!!!!!!!!!!!!!!
####Typically you want to choose a value high enough that you capture the diversity present in samples with high counts, but low enough that you don’t get rid of a ton of your samples
#############let's do some alpha referaction plot to calculate the best depth ##takes time if you set wide range############

qiime diversity alpha-rarefaction --i-table table.qza --m-metadata-file sample_metadata.tsv --o-visualization alpha_rarefaction_curves.qzv  --p-min-depth 10 --p-max-depth 20000
#############################################################

qiime tools view alpha_rarefaction_curves.qzv
#############let's calculate the alpha and beta diversity#############

qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table feature-table.qza --p-sampling-depth 800 --m-metadata-file new_sample_metadata.tsv --output-dir metrics
#PS1:In our work, we used the beta diversity from qiime2 instead of Phyloseq due to this issue, which we also witness
#https://github.com/joey711/phyloseq/issues/956
#PS2: When we run after that for non-sorted sequences we will use 1100  as sequnce depth
###########################
#Thank you and see you soon

