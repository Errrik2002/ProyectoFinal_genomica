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

#################analisis de frecuencia##############
total1168 = median(sample_sums(psobj1))


ps1168filtrado <- filter_taxa(psobj1, function(x) sum(x > total1168*0.20) > 0, TRUE)
ps1168filtrado



plot_heatmap(ps1168filtrado, method = "NMDS", distance = "bray")

plot_heatmap(ps1168filtrado, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Genus", taxa.order = "Genus", 
             trans=NULL)


psobj1fixed <- psobj1 %>% tax_fix()

BaseMGYS00001168 <- psobj1fixed %>%
       tax_agg("Species") %>%
       ps_arrange(.target = "otu_table") %>%
       otu_get() 


#save(psobj2, file = "04_saved/psobj2.RData")
save(BaseMGYS00001168, file = "04_saved/BaseMGYS00001168.RData")

View(tax_table(psobj1))









































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
pdf("03_results/BAR_MGYS00001373.pdf")


plot_bar(psobj2, fill = "Genus")

# Closing the graphical device
dev.off() 

#################analisis de frecuencia##############
total1373 = median(sample_sums(psobj2))


ps1373filtrado <- filter_taxa(psobj2, function(x) sum(x > total1373*0.20) > 0, TRUE)
ps1373filtrado



plot_heatmap(ps1373filtrado, method = "NMDS", distance = "bray")

plot_heatmap(ps1373filtrado, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Species", taxa.order = "Species", 
             trans=NULL)


psobj2fixed <- psobj2 %>% tax_fix()

BaseMGYS00001373 <-psobj2fixed %>%
       tax_agg("Species") %>%
       ps_arrange(.target = "otu_table") %>%
       otu_get()


#save(psobj2, file = "04_saved/psobj2.RData")
save(BaseMGYS00001373, file = "04_saved/BaseMGYS00001373.RData")


View(tax_table(psobj2))




































############3Codigo de May##################

library (phyloseq)
library (MGnifyR)

# cargamos los dstos desde mgnifyR 
# comentar esta parte del script 
bioma <- mgnify_client(usecache = T, cache_dir = '/tmp/MGnify_cache')
lista_acc <- mgnify_analyses_from_studies(bioma, "MGYS00000943", usecache = T)
meta_dataframe <- mgnify_get_analyses_metadata(bioma, lista_acc, usecache = T )
vbioma <- mgnify_get_analyses_phyloseq(bioma, meta_dataframe$analysis_accession, usecache = T)
vbioma
save(vbioma, file= "04_saved/vbioma.RData")
#save(psobj2, file = "04_saved/psobj2.RData")
load("04_saved/vbioma.RData")
# vamos a explorar las opciones dentro del bioma 
vbioma
otu_table (vbioma)
# Vemos que existen 35 muestras, se formaron 608 ASV(que van a corresponder a especies 
# diferentes dentro de cada una de las muestras, con esto ya acomodamos los organismos 
# que se encontraron dentro del muestreo)

View(sample_data (vbioma))
View(sample_variables (vbioma))
# tenemos 57 variables dentro de las que encontramos el identificador para cada una de las muestras 
# que son 35, el tipo de experimentacion que se realizo (amplicones), el tiempo de analisis 
# la plataforma que se utilizo para el analisisi que fue ilumina miseq, el tamaño de las seucencias 
# que fueron utilizadas, el tamaño que adoptaron las secuencias antes y despues del filtraje de informacion
# tambien hay una seccion en donde se explican los procedimientos utilizados para la recoleccion de muestras 


tax_table (vbioma)
rank_names (vbioma)
# donde encontramos la clasificacion taxonomica de los organismos encontrados dentro de la muestra 
# tenemos hasta especies, pero tambien hay que señalar que en este analisis tal vez no es tan conveniente
# utilizar especies o generos, porque no todos los organismos son identificados hasta este nivel taxonomico 
# aunque podemos señalar que en genero si hay bastantes organismos bien identificados. 

####################################################################
plot_bar(vbioma, fill = "Genus")

library(ggplot2)


#It is better to only consider the most abundant OTUs for heatmaps. 
#For example one can only take OTUs that represent at least 20% of reads
#in at least one sample. Remember we normalized all the sampples 
#to median number of reads (total)

total = median(sample_sums(vbioma))


vbioma_filtrado <- filter_taxa(vbioma, function(x) sum(x > total*0.20) > 0, TRUE)
vbioma_filtrado

plot_heatmap(vbioma_filtrado, method = "NMDS", distance = "bray")

  plot_heatmap(vbioma_filtrado, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Genus", taxa.order = "Genus", 
             trans=NULL)





########################
  BiocManager::install(c("microbiome", "ComplexHeatmap"), update = FALSE)
  
  install.packages(
    "microViz",
    repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
  )
  
  
  library(microViz)
  library(microbiome) 
  library(ComplexHeatmap)
    
  # Arranging by decreasing Bacteroides abundance
  vbiomafixed <- vbioma %>% tax_fix()
  
  Basevbioma <-vbiomafixed %>%
    tax_agg("Species") %>%
    ps_arrange(.target = "otu_table") %>%
    otu_get() 

  #save(psobj2, file = "04_saved/psobj2.RData")
save(Basevbioma, file = "04_saved/Basevbioma.RData")




