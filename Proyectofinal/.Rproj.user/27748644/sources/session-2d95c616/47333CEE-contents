BiocManager::install("phyloseq")
library(phyloseq)
install.packages("devtools")
library(devtools)

devtools::install_github("beadyallen/MGnifyR")
library(MGnifyR)



############################################ KOREA ###########################################
  mgclnt <- mgnify_client(usecache = T, cache_dir = '/tmp/MGnify_cache') #esto se conecta a base de dato
mgclnt
#Retrieve the list of analyses associated with a study
accession_list <- mgnify_analyses_from_studies(mgclnt, "MGYS00001168", usecache = T)
#nos permite conectarnos a esa base de datos en especifico, el numerito en verde

#Download all associated study/sample and analysis metadata
meta_dataframe <- mgnify_get_analyses_metadata(mgclnt, accession_list, usecache = T )

meta_dataframe #guarda los datos que pedimos

#Convert analyses outputs to a single `phyloseq` object
psobj1 <- mgnify_get_analyses_phyloseq(mgclnt, meta_dataframe$analysis_accession, usecache = T)
psobj1 #lo convierte ya en un objeto phyloseq

save(psobj1, file = "04_saved/psobj.RData")
#save(taxa,file="04_Processed_Data/taxa.RData")
load("04_saved/psobj.RData")

rank_names(psobj1)

sample_names(psobj1)

sample_variables(psobj1)




#
plot_net(psobj1, type = "taxa", point_label = "Class", point_size = 10, point_alpha = 0.5, maxdist = 0.5, color = "Species", distance = "bray", laymeth = "auto") 

######################################
pdf("03_results/BAR_MGYS00001168.pdf")


plot_bar(psobj1, fill = "Genus")

# Closing the graphical device
dev.off() 


###################################### MGYS00001373 FLORENCIA ##############################################################

mgclnt <- mgnify_client(usecache = T, cache_dir = '/tmp/MGnify_cache') #esto se conecta a base de dato
mgclnt
#Retrieve the list of analyses associated with a study
accession_list <- mgnify_analyses_from_studies(mgclnt, "MGYS00001373", usecache = T)
#nos permite conectarnos a esa base de datos en especifico, el numerito en verde

#Download all associated study/sample and analysis metadata
meta_dataframe <- mgnify_get_analyses_metadata(mgclnt, accession_list, usecache = T )

meta_dataframe #guarda los datos que pedimos

#Convert analyses outputs to a single `phyloseq` object
psobj2 <- mgnify_get_analyses_phyloseq(mgclnt, meta_dataframe$analysis_accession, usecache = T)
psobj2 #lo convierte ya en un objeto phyloseq

save(psobj2, file = "04_saved/psobj2.RData")
#save(taxa,file="04_Processed_Data/taxa.RData")
load("04_saved/psobj2.RData")

rank_names(psobj2)

sample_names(psobj2)

sample_variables(psobj2)


plot_bar(psobj2, fill = "Genus")

#

######################################
pdf("03_results/BAR_MGYS00001168.pdf")


plot_bar(psobj1, fill = "Genus")

# Closing the graphical device
dev.off() 