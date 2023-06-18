#ver expresion diferencial en oncogenesis de papiloma

#no se si deberiamos usar vph 18 o 16?



script1 <-read.csv("01_RawData/ed1cvs16.tsv",dec=",",sep="\t")

View(script1)



#examdifexpbueno <- subset(examdifexp, examdifexp$logFC <0 | examdifexp$P.Value < 0.05) #EXTRAEMOS AQUELLOS SUBEXPRESADOS
#examdifexpbueno #GENES sub expresados

#genessobreexpresados <- subset(examdifexp, examdifexp$logFC >0 | examdifexp$P.Value < 0.05) SOBREEXPRESADOS



################SUBEXPRESADOS############

subexp <- subset(script1, script1$logFC <0 | script1$P.Value < 0.05)
subexp

#masexpre111 <- sort(as.numeric(pvalores11), decreasing = F)[1:10]

top10sub <- sort(subexp$logFC, decreasing = T)[1:10]
top10sub

library(dplyr)
elgensub <- subexp %>% filter_all(any_vars(. %in% c(top10sub)[1:10])) 
elgensub[1]



##############SOBREEXPRESADOS##############

sobreexp <- subset(script1, script1$logFC > 0 | script1$P.Value < 0.05)
sobreexp


top10sob <- sort(sobreexp$logFC, decreasing = T)[1:10]
top10sob

elgensobre <- sobreexp %>% filter_all(any_vars(. %in% c(top10sob)[1:10])) 
elgensobre[1]




#########################################################################
############################VPH18################
#########################################################################


vph18 <-read.csv("01_RawData/ed1cvs18.tsv",dec=",",sep="\t")

View(vph18)


subexp18 <- subset(vph18, vph18$logFC <0 | vph18$P.Value < 0.05)
subexp18

#masexpre111 <- sort(as.numeric(pvalores11), decreasing = F)[1:10]

top10sub18 <- sort(subexp18$logFC, decreasing = T)[1:10]
top10sub18

library(dplyr)

elgensub18 <- subexp18 %>% filter_all(any_vars(. %in% c(top10sub18)[1:10])) 
elgensub18[1]

#
#
#
#
#
#
#
#
#
#


##############SOBREEXPRESADOS##############

sobreexp18 <- subset(vph18, vph18$logFC > 0 | vph18$P.Value < 0.05)
sobreexp18


top10sob18 <- sort(sobreexp18$logFC, decreasing = T)[1:10]
top10sob18

elgensobre18 <- sobreexp18 %>% filter_all(any_vars(. %in% c(top10sob18)[1:10])) 
elgensobre18[1]

#5 CUST_13113_PI425136549	ENST00000421323	ENSG00000229839
#LA 10 es p10101
















######################################################
######################################################NUEVAS BASES DE DATOS ################################################
######################################################












##############GSE138081  Identification of deregulated pathways, key regulators, and novel miRNA-mRNA interactions in HPV-mediated transformation.
###################################################################

vph18_8081 <-read.csv("01_RawData/GSE138081_HPV18vsNormal.tsv",dec=",",sep="\t")
View(vph18_8081)




################NO TIENE UN LOGFC, TONS NO SE PUEDE###################





########################################################
########################GSE67522 Genome-wide analysis of gene expression to identify the probably functionally relevant pathways in cervical cancer progression
########################################################


vph16_7522 <-read.csv("01_RawData/GSE67522_HPV16vsControl.tsv",dec=",",sep="\t")

View(vph16_7522)


#subexp18 <- subset(vph18, script1$logFC <0 | script1$P.Value < 0.05)
sub16exp_7522 <- subset(vph16_7522, vph16_7522$logFC <0 | vph16_7522$P.Value < 0.05)

sub16exp_7522

#masexpre111 <- sort(as.numeric(pvalores11), decreasing = F)[1:10]

top10sub16_7522 <- sort(sub16exp_7522$logFC, decreasing = T)[1:10]
top10sub16_7522 ###########SACAMOS LOS VALORES MAS BAJOS########

library(dplyr)

elgensub16_7522 <- sub16exp_7522 %>% filter_all(any_vars(. %in% c(top10sub16_7522)[1:10])) 
elgensub16_7522[8:9]

##########SOBREEXPRESADOS----------------------------

sobr16_7522 <- subset(vph16_7522, vph16_7522$logFC > 0 | vph16_7522$P.Value < 0.05)
sobr16_7522


top10sob16_7522 <- sort(sobr16_7522$logFC, decreasing = T)[1:10]
top10sob16_7522

elgensobre16_7522 <- sobr16_7522 %>% filter_all(any_vars(. %in% c(top10sob16_7522)[1:10])) 
elgensobre16_7522[8:9]













############################################
#################### GSE151666 RNA sequencing of pre-treatment primary cervical cancer
###########################
vph16_1666 <-read.csv("01_RawData/GSE151666_VPH16vsCncer_noVPH.tsv",dec=",",sep="\t")

View(vph16_1666)


#subexp18 <- subset(vph18, script1$logFC <0 | script1$P.Value < 0.05)
sub16exp_1666 <- subset(vph16_1666, vph16_1666$log2FoldChange <0 | vph16_1666$pvalue < 0.05)

sub16exp_1666

#masexpre111 <- sort(as.numeric(pvalores11), decreasing = F)[1:10]

top10sub16_1666 <- sort(sub16exp_1666$log2FoldChange, decreasing = T)[1:10]
top10sub16_1666 ###########SACAMOS LOS VALORES MAS BAJOS########

library(dplyr)

elgensub16_1666 <- sub16exp_1666 %>% filter_all(any_vars(. %in% c(top10sub16_1666)[1:10])) 
elgensub16_1666[8:9]

##########SOBREEXPRESADOS----------------------------

sobr16_1666 <- subset(vph16_1666, vph16_1666$log2FoldChange > 0 | vph16_1666$pvalue < 0.05)
sobr16_1666


top10sob16_1666 <- sort(sobr16_1666$log2FoldChange, decreasing = T)[1:10]
top10sob16_1666

elgensobre16_1666 <- sobr16_1666 %>% filter_all(any_vars(. %in% c(top10sob16_1666)[1:10])) 
elgensobre16_1666[8:9]









################################ VPH18 VS CANCER NORMAL######################################

vph18_1666 <-read.csv("01_RawData/GSE151666_VPH18vsCancernormal.tsv",dec=",",sep="\t")

View(vph18_1666)


#subexp18 <- subset(vph18, script1$logFC <0 | script1$P.Value < 0.05)
sub18exp_1666 <- subset(vph18_1666, vph18_1666$log2FoldChange <0 | vph18_1666$pvalue < 0.05)

sub18exp_1666

#masexpre111 <- sort(as.numeric(pvalores11), decreasing = F)[1:10]

top10sub18_1666 <- sort(sub18exp_1666$log2FoldChange, decreasing = T)[1:10]
top10sub18_1666 ###########SACAMOS LOS VALORES MAS BAJOS########

library(dplyr)

elgensub18_1666 <- sub18exp_1666 %>% filter_all(any_vars(. %in% c(top10sub18_1666)[1:10])) 
elgensub18_1666[8:9]

##########SOBREEXPRESADOS----------------------------

sobr18_1666 <- subset(vph18_1666, vph18_1666$log2FoldChange > 0 | vph18_1666$pvalue < 0.05)
sobr18_1666


top10sob18_1666 <- sort(sobr18_1666$log2FoldChange, decreasing = T)[1:10]
top10sob18_1666

elgensobre18_1666 <- sobr18_1666 %>% filter_all(any_vars(. %in% c(top10sob18_1666)[1:10])) 
elgensobre18_1666[8:9]
