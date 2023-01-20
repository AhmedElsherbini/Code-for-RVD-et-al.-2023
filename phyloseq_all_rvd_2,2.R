############################
#authour: Ahmed Elsherbini 
#ref:https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
#ref:
#update : 22-03-2022
########################
library(phyloseq) #for alpha and beta diversity measurements
library(qiime2R)
library("ape")
library(DESeq2)
library(ggplot2)
library(dplyr)  
library(microbial)
library(phyloseqCompanion)
library(ggpubr)
#########################################################################
setwd("/media/ahmed/CC69-620B6/00_Ph.D/DATA_results/3_Rv.D/4_all_RVD_MaxEE2,2")
#########################################################################
"""
SVs<-read_qza("feature-table.qza")
metadata<-read_q2metadata("New_metadata.tsv")
taxonomy<-read_qza("taxonomy_ncct.qza")
taxonomy<-parse_taxonomy(taxonomy$data)
shannon<-read_qza("shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID")
gplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))
"""
###################################################read########
physeq <- qza_to_phyloseq(
  features="feature-table.qza",
  tree="rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "new_sample_metadata.tsv")
physeq

####################################subset the samples#############################################

physeq = subset_samples(physeq, type != "0")
physeq = subset_samples(physeq, donor.id != "0")
physeq = subset_samples(physeq, type != "control")
#physeq = subset_samples(physeq, cst.manual != "")

####################################taxwork#############################################
#using micobiobal 

phy <- normalize(physeq, method = "relative")


plotbar(phy,level="Phylum")
plotbar(phy,level="Genus", group="type",fontsize.x = 10)

##################################rarefication#######################################
#refraction will be done on the result of alpha refraction curves of qiime2
rarecurve(t(otu_table(ps)), step=50, cex=0.5)

physeq.rarefied = rarefy_even_depth(physeq, rngseed=1, sample.size=800, replace=F)
##################################alpha#########################
plotalpha(physeq.rarefied, group = "type")

###########################################################
#using satistics
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

a_my_comparisons <- list(c("pre-sort", "neg"),c("pre-sort","pos"),c("pos","neg"))
plot_richness(physeq.rarefied, x="type")+
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

#####################################pvalue############################################
rich = estimate_richness(physeq.rarefied, measures = c("invsimpson"))
write.csv(rich,"estimated_simpson.csv")
wilcox.Shannon <- pairwise.wilcox.test(rich$Shannon, 
                                       sample_data(physeq.rarefied)$site, 
                                       p.adjust.method = "BH")


wilcox.Shannon <- pairwise.wilcox.test(rich$Shannon, 
                                       sample_data(physeq.rarefied)$site, 
                                       p.adjust.method = "BH")

wilcox.Observed <- pairwise.wilcox.test(rich$Observed, 
                                        sample_data(physeq.rarefied)$site, 
                                        p.adjust.method = "BH")

wilcox.Observed <- pairwise.wilcox.test(rich$Observed, 
                                        sample_data(physeq.rarefied)$site, 
                                        p.adjust.method = "BH")



wilcox.Shannon[["p.value"]]
wilcox.Observed[["p.value"]]

#######################################Beta##############################
#sample_data(physeq.rarefied)[ , 2] <- sample_data(physeq.rarefied)[ ,1]
dist = phyloseq::distance(physeq.rarefied, method="bray")
ordination = ordinate(physeq.rarefied, method="PCoA", distance=dist)
plot_ordination(physeq.rarefied, ordination, color="type") + 
  theme_classic() +
  theme(strip.background = element_blank())


dist = phyloseq::distance(physeq.rarefied, method="wunifrac", binary = TRUE)
ordination = ordinate(physeq.rarefied, method="PCoA", distance=dist)
plot_ordination(physeq.rarefied, ordination, color="type") + 
  theme_classic() +
  theme(strip.background = element_blank())

plotbeta(phy, group="type")
beta <-betatest(phy,group="type")

GP.ord = ordinate(physeq.rarefied, "PCoA", "unifrac", weighted=TRUE)
p2 = plot_ordination(physeq.rarefied, GP.ord, type="samples", color="CST", shape="sex") 
p2 + geom_point(size=4) 
##########################################diffexpersion_deseq2##################
#important diff abundacne is pairwise analysis only dummy 
physeq = subset_samples(physeq, type != "pre-sort")
physeq = subset_samples(physeq, type != "control")

res <- difftest(physeq,group="type")
plotdiff(res,level="Genus",padj=0.001)
#######################################lefese##################
phyloseq2lefse(
  physeq,
  covars = c("type"),
  file.name = "lefse_data.txt",
  taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  transpose.otus = TRUE)
################################################################
#https://www.yanh.org/2021/01/01/microbiome-r/

ps = physeq
ps = subset_samples(ps, type != "control")
ps = subset_samples(ps, type != "pre-sort")
sample_data(ps)$type <- as.factor(sample_data(ps)$type) # factorize for DESeq2
ps.taxa <- tax_glom(ps, taxrank = 'Family', NArm = FALSE)
ps.taxa.pse <- ps.taxa
#the most important step
otu_table(ps.taxa.pse) <- otu_table(ps.taxa) + 1

ps.taxa.pse.sub <- subset_samples(ps.taxa.pse, type %in% c("pos", "neg"))
ds = phyloseq_to_deseq2(ps.taxa.pse.sub, ~ type)
ds = DESeq(ds, test="Wald", fitType="parametric")
alpha = 0.05 
res = results(ds, alpha=alpha)
res = res[order(res$padj, na.last=NA), ]

taxa_sig = rownames(res[1:20, ]) # select bottom 20 with lowest p.adj values

ps.taxa.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)

ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)

ps.taxa.rel.sig <- prune_samples(colnames(otu_table(ps.taxa.pse.sub)), ps.taxa.rel.sig)

matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))

rownames(matrix) <- as.character(tax_table(ps.taxa.rel.sig)[, "Family"])

metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))
annotation_col = data.frame( 
  `type` = as.factor(metadata_sub$type), 
  check.names = FALSE)


rownames(annotation_col) = rownames(metadata_sub)

rownames(metadata_sub)

annotation_row = data.frame(
  Phylum = as.factor(tax_table(ps.taxa.rel.sig)[, "Phylum"])
)

rownames(annotation_row) = make.names(rownames(matrix), unique = TRUE)


# ann_color should be named vectors
phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
names(phylum_col) = levels(annotation_row$Phylum)
ann_colors = list(
  `type` = c(pos = "purple", neg = "yellow"),
  Phylum = phylum_col
)

ComplexHeatmap::pheatmap(matrix, scale= "row", 
                         annotation_col = annotation_col, 
                         annotation_row = annotation_row, 
                         annotation_colors = ann_colors)

#####################################################
#https://micca.readthedocs.io/en/latest/phyloseq.html#otu-differential-abundance-testing-with-deseq2
sample_data(ps)$type <- as.factor(sample_data(ps)$type) # factorize for DESeq2
otu_table(ps.taxa.pse) <- otu_table(ps.taxa) + 1

ds = phyloseq_to_deseq2(ps.taxa.pse.sub, ~ type)
ds = DESeq(ds, test="Wald", fitType="parametric")
ds = DESeq(ds)
res = results(ds, contrast=c("type", "pos", "neg"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig
res_sig = cbind(as(res_sig, "data.frame"), as(tax_table(ps)[rownames(res_sig), ], "matrix"))
ggplot(res_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
  geom_jitter(size=3, width = 0.2) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))