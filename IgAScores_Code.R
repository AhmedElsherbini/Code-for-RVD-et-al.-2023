################################################
#original author: Dr. Rob Van Dalin
#received on 26.11.2021
#modified by : Ahmed Elsherbini - 30-11-202
#last update : 03-06-2022
#based on: https://github.com/microbialman/IgAScores
#https://environmentalcomputing.net/graphics/multivariate-vis/heatmaps/
##############################################
library(IgAScores)
library(stringr)
library(xlsx)
library(tidyverse)
library(pheatmap)
library(gplots)
library(dendextend)
##############################################
############
BLUE <- "#076fa2"
RED <- "#E3120B"
BLACK <- "#202020"
GREY <- "grey50"
############
getwd()
setwd("/media/ahmed/CC69-620B5/00_Ph.D/DATA_results/3_Rv.D/0_IgAScores_AE")
#data frames with counts for the bacterial taxa in the IgA+ and IgA- fractions
igapos =  read.xlsx2("level-7_IgAScores.xlsx",  sheetName = "POS_transposed")
#convert the whole data frame to str except the first column but the names are really looooong :))
igapos[,-1] = as.data.frame(lapply(igapos[,-1], as.numeric))
#here to go to the spp level
#igapos$index = sapply(strsplit(as.character(igapos$index),';g__'), "[", 2)
###
iganeg =  read.xlsx2("level-7_IgAScores.xlsx",  sheetName = "NEG_transposed")
iganeg[,-1] = as.data.frame(lapply(iganeg[,-1], as.numeric))

#convert the counts to relative abundances using the included helper function #just skip the first column
igapos[,-1] = relabund(igapos[,-1])
iganeg[,-1] = relabund(iganeg[,-1])

#iga+ and iga- fraction sizes per sample
possize =  read.xlsx2("level-7_IgAScores.xlsx",  sheetName = "fraction_POS")
possize = as.numeric(possize[1,])
negsize =  read.xlsx2("level-7_IgAScores.xlsx",  sheetName = "fraction_NEG")
negsize = as.numeric(negsize[1,])
#set a pseudo count for handling zero values in some scoring methods
#default method is the probability ratio "probratio"
prscores = igascores(posabunds = igapos[,-1], negabunds = iganeg[,-1], 
                      possizes = possize, negsizes = negsize, 
                      pseudo =  0.00001)
#how to add the first column as a index
prscores = cbind(index = igapos$index,prscores)
#delete sample X251 of Colleague C.beck as 
prscores = subset(prscores, select = -c(X251) )
#Write as the csv
write.csv(prscores,"prscores_ahmed.csv", row.names = FALSE )
prscores$index
#End of part 1

##########################################################
#part 2
#let show the interactive spp.
#please remove every thing before K__Bacteria;p__
#please remove every thing before *;g__
#The cleaning step #step 1 keep genus spp like stpah epid #step2 remove what is #step 3 remove raws with no spp
prscores$index <- gsub("^.*;g__","", prscores$index)
prscores = prscores[-grep("k__Bacteria;p__", prscores$index),]
prscores = prscores[grep(";s__", prscores$index),]

#################################################################
#please remove every thing before K__Bacteria;p__
#can you count how many NA in each raw and add rgen to the end of the dataframe as seprate column named rawSums
prscores$count_na <- rowSums(is.na(prscores))
#Based on this count of NA can you order my dataframe so I can know the most interacting spp. The less The na the more interacting 
#you cab find them in the bottum of the dataframe
prscores <- prscores[with(prscores,order(-count_na)),]
#Print this map
write.csv(prscores,"interactive_spp.csv", row.names = FALSE )
#I can see that the most frequent 20 contain show in > 50% of profiles I can call them frequent spp
prscores_intreactive =  tail(prscores,20)
#we will convert the NA to prevalence that reflect how many time each frequent spp is showing
#why 49 ? 49  in the number of samples (number of columns in the dataframe - index spp and percentage)
prscores_intreactive$count_na = (((48 - prscores_intreactive$count_na)/48)*100) 
#now you can see that the buttom of the dataframe contains the most interactive spp 
#mak a dataframe for plot
count = prscores_intreactive$count_na
name = factor(prscores_intreactive$index, levels = prscores_intreactive$index)
#y = seq(length(name)) * 0.9
data <- data.frame(
  count , 
  name)

#plot
plt <- ggplot(data)+
  labs(y = "spp.", x = "% of samples")+
  geom_col(aes(count, name), fill = BLACK , width = 0.6)
plt
#End of part 2 
##############################################################################
#Part 3

#let make a heat map to link human profiles col with Bactria spp. raws
#delete count of NA
################################################
#prscores = igascores(posabunds = igapos[,-1], negabunds = iganeg[,-1], 
#                     possizes = possize, negsizes = negsize, 
#                     pseudo =  0.00001)
#how to add the first column as a index
#prscores = cbind(index = igapos$index,prscores)
#delete sample X251 of Colleague C.beck as 
#prscores = subset(prscores, select = -c(X251) )


############################################
prscores_intreactive
#let's delete the count na 
prscores_intreactive$count_na = NULL
#MATRIX without the first column which is spp names
m = as.matrix(prscores_intreactive[,-1])
#replace na by zero 
m[is.na(m)] <- 0
rownames(m) = prscores_intreactive$index
# Some haet mapping is good !
#https://mycolor.space/
#Make a dataframe table of the 
#data <- as.data.frame(read_csv("normalized_wetlab.csv",show_col_types = FALSE))

#data2 <- data[,-1]
#rownames(data2) <- data[,1]
#m = as.matrix(data2)
#m[is.na(m)] <- 0
#m = t(m)
#m = log(m)
metadata <- as.data.frame(read_csv("holistic_metadata.tsv",show_col_types = FALSE))
metadata2 <- metadata[,-1]
rownames(metadata2) <- metadata[,1]

#here we have to make sure that first column is the rawnames not single column
#metadata 2 we delete the first column
#metadata2 <- metadata[,-1]
#then we use the first column to name our metadata
#rownames(metadata2) <- metadata[,1]
#rownames(metadata) <- colnames(m)
#colnames(m) = rownames(metadata2)
#make a list of colors as you prefer
#https://jokergoo.github.io/ComplexHeatmap-reference/book/integrate-with-other-packages.html
annotation_colors = list(Gender = c(Male ="#6ef88a",Female= "#d357fe"),CST =c(CST1='#E1ED15',CST3='#77C25E',CST4='#ED1525',CST5='#3315ED',Unclassified='#0C0C0B'),Age =c("grey", "firebrick"))
#let's draw
pheatmap(m,show_rownames = T,show_colnames = T,fontsize_row = 8, fontsize_col = 8,annotation_col = metadata2,annotation_colors=annotation_colors)
#End of part 3
##########################################################
#let's show the postive and negative spp among our frequnts spp.
#convert NAs to zeros
prscores_intreactive[is.na(prscores_intreactive)] <- 0
#add 10 to all of you analysis to convert them postive
prscores_intreactive = data.frame(t(prscores_intreactive))
colnames(prscores_intreactive) <- prscores_intreactive[1,]

prscores_intreactive = prscores_intreactive[-1,]
rownames(prscores_intreactive) <- NULL
#To make sure no values will remove other values, I add 10 
prscores_intreactive[,-1] = prscores_intreactive[,-1] + 10
#then make a mean on the left of the 
prscores_intreactive$mean <-apply(prscores_intreactive[,-1],1,mean)
prscores_intreactive$mean
#order based on the mean value
prscores_intreactive <- prscores_intreactive[with(prscores_intreactive,order(-mean)),]

prscores_intreactive[,-1]  = prscores_intreactive[,-1] - 10

data <- data.frame(
  count = prscores_intreactive$mean, 
  name = factor(prscores_intreactive$index, levels = prscores_intreactive$index),
  y = seq(length(names)) * 0.9)

plt <- ggplot(data)+
  labs(y = "spp.", x = "IgA score")+
  geom_col(aes(count, name), fill = BLACK , width = 0.6)
plt

write.csv(prscores_intreactive,"ordered_prscores.csv", row.names = FALSE )
########the high guys##############################
prscores_high = head(prscores_intreactive, 10)
write.csv(prscores_high,"high_prscores.csv", row.names = FALSE )
prscores_high$mean = NULL
m_h = as.matrix(prscores_high[,-1])
#names please
rownames(m_h) = prscores_high$index
# Some haet mapping is good !
pheatmap(m_h,show_rownames = T,show_colnames = T,fontsize_row = 8, fontsize_col = 8)

#plot
plt <- ggplot(data)+
  labs(y = "spp.", x = "%")+
  geom_col(aes(count, name), fill = BLACK , width = 0.6)
plt

####the low spp####################################################################.


#####################################################################
#seperate top ten
prscores_high = head(prscores, 10)
#se[erate lowest ten
prscores_low =  tail(prscores,10)
#merge them 
prscores = rbind(prscores_high, prscores_low)
#remove means 
prscores$Means = NULL
#prscores1 = prscores1[-1,]
#prscores2 <- data.frame(apply(prscores1[,-1], 2, function(x) as.numeric(as.character(x))))
# make your data as a matrix
m = as.matrix(prscores[,-1])
#name it please
rownames(m) = prscores$index
#heatmap(m)
pheatmap(m,show_rownames = T,show_colnames = T,fontsize_row = 5, fontsize_col = 5)
#############################################################
#the cherry picking analysis
prscores2 <- read.csv("sel_input_very_biased.csv")
# NA TO ZERO
prscores2[is.na(prscores2)] <- 0
# REMOVE ALL RAW OF ZERO
prscores2[apply(prscores2[,-1], 1, function(x) !all(x==0)),]
#KEEP THE spp level only
prscores2$index <- gsub("^.*;g__","", prscores2$index)
# make your data as a matrix
m = as.matrix(prscores2[,-1])
#name it please
rownames(m) = prscores2$index
heatmap(m)
pheatmap(m,show_rownames = T,show_colnames = T,fontsize_row = 5, fontsize_col = 7)

############################################################



