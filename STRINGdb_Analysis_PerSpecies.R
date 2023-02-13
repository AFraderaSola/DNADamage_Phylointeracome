###########################
####### Libraries #########
###########################

library(STRINGdb)
library(readr)
library(tidyverse)
library(ggplot2)
library(ggnet)
library(igraph)
library(ggpubr)
library(viridis)
library(ggnetwork)
set.seed(666)

###########################
#### STRINGdb Analysis ####
###########################

#### Included species ####

tx_id <- c(3702, 224308, 6239, 511145, 64091, 9606, 9606, 4932, 4896, 5691, 5911, 4577) #ncbi taxonmy ID

Species <- c("A. thaliana",
             "B. subtilis",
             "C. elegans",
             "E. coli",
             "H. salinarum",
             "H. sapiens (HEK293)",
             "H. sapiens (HeLa)",
             "S. cerevisiae",
             "S. pombe",
             "T. brucei",
             "T. thermophila",
             "Z. mays")

##### Input Files ####

# Enriched proteins

proteins <- read.csv(file = "00_Rawdata/ProteinsPerSpecies.csv",header = T)

proteins <- proteins[,c(1,3,4,10,12,13)]

proteins$Protein.IDs <- gsub(pattern = ";.*",replacement = "",x = proteins$Protein.IDs)

proteins$additional_protein_name <- gsub(pattern = ";.*",replacement = "",x = proteins$additional_protein_name)

proteins[proteins$repair == "repair",]$repair <- "Yes"

proteins[proteins$repair == "no",]$repair <- "No"

colnames(proteins)[grep("repair",colnames(proteins))] <- "Repair protein"

species <- unique(proteins$Species)

# STRINGdb (v11.5) interactions per species

files <- list.files(path = "00_Rawdata/",pattern = ".*pro.*txt$")

files <- c(files[1:6],files[6],files[7:11]) # Repeat H. sapiens; one for each cell line.

# Analysis loop per species.

plot_list <- c()

nodes_df <- c()

edges_df <- c()

scores_df <- c()

for (i in 1:length(species)) {
  
  # Get species STRINGdb
  
  string_db <- STRINGdb$new( version="11.5", species=tx_id[i],score_threshold = 150, network_type= "full", input_directory="")
  
  loop_proteins <- proteins[proteins$Species == species[i],]
  
  loop_proteins <- string_db$map(loop_proteins,"Protein.IDs",removeUnmappedRows=F)
  
  # Get STRINGdb mapping (first to protein IDs and, if NA, to the additional_protein_name)
  
  if (any(is.na(loop_proteins$STRING_id))) {
    
    loop_proteins_NA <- loop_proteins[is.na(loop_proteins$STRING_id),]
    
    loop_proteins_NA <- loop_proteins_NA[,1:6]
    
    loop_proteins_NA <- string_db$map(loop_proteins_NA,"additional_protein_name",removeUnmappedRows=F)
    
    loop_proteins <- proteins[proteins$Species == species[i],]
    
    loop_proteins <- string_db$map(loop_proteins,"Protein.IDs",removeUnmappedRows=T)
    
    loop_proteins <- rbind(loop_proteins, loop_proteins_NA)
    
    loop_proteins <- loop_proteins[!is.na(loop_proteins$STRING_id),]
    
    loop_map_df <- loop_proteins
    
    loop_map_df <- loop_map_df[,c(2,3,1,4,7)]
    
  }else{
    
    loop_map_df <- loop_proteins
    
    loop_map_df <- loop_map_df[,c(2,3,1,4,7)]
    
  }
  
  # Draw network from interactions
  
  nodes <- loop_proteins[,c(4,1,7,3,5,6,2)]
  
  nodes <- nodes[,c(1,2,3,5)]
  
  nodes <- unique(nodes)
  
  nodes$Experiment <- "Some"
  
  nodes[nodes$additional_protein_name %in% loop_proteins[loop_proteins$experiment == "8oxoG",]$additional_protein_name,]$Experiment <- "8oxoG"
  
  nodes[nodes$additional_protein_name %in% loop_proteins[loop_proteins$experiment == "abasic",]$additional_protein_name,]$Experiment <- "abasic"
  
  nodes[nodes$additional_protein_name %in% loop_proteins[loop_proteins$experiment == "RNAbase",]$additional_protein_name,]$Experiment <- "RNAbase"
  
  if (length(intersect(loop_proteins[loop_proteins$experiment == "8oxoG",]$additional_protein_name,
                       loop_proteins[loop_proteins$experiment == "abasic",]$additional_protein_name)) > 0) {
    
    nodes[nodes$additional_protein_name %in% intersect(loop_proteins[loop_proteins$experiment == "8oxoG",]$additional_protein_name,
                                                   loop_proteins[loop_proteins$experiment == "abasic",]$additional_protein_name),]$Experiment <- "8oxoG & abasic"
  }
  
  if (length(intersect(loop_proteins[loop_proteins$experiment == "8oxoG",]$additional_protein_name,
                       loop_proteins[loop_proteins$experiment == "RNAbase",]$additional_protein_name)) > 0) {
    
    nodes[nodes$additional_protein_name %in% intersect(loop_proteins[loop_proteins$experiment == "8oxoG",]$additional_protein_name,
                                                   loop_proteins[loop_proteins$experiment == "RNAbase",]$additional_protein_name),]$Experiment <- "8oxoG & RNAbase"
  }
  
  if (length(intersect(loop_proteins[loop_proteins$experiment == "abasic",]$additional_protein_name,
                       loop_proteins[loop_proteins$experiment == "RNAbase",]$additional_protein_name)) > 0) {
    
    nodes[nodes$additional_protein_name %in% intersect(loop_proteins[loop_proteins$experiment == "abasic",]$additional_protein_name,
                                                   loop_proteins[loop_proteins$experiment == "RNAbase",]$additional_protein_name),]$Experiment <- "abasic & RNAbase"
  }
  
  if (length(intersect(loop_proteins[loop_proteins$experiment == "8oxoG",]$additional_protein_name,
                       intersect(loop_proteins[loop_proteins$experiment == "abasic",]$additional_protein_name,
                                 loop_proteins[loop_proteins$experiment == "RNAbase",]$additional_protein_name))) > 0) {
    nodes[nodes$additional_protein_name %in% intersect(loop_proteins[loop_proteins$experiment == "8oxoG",]$additional_protein_name,
                                                   intersect(loop_proteins[loop_proteins$experiment == "abasic",]$additional_protein_name,
                                                             loop_proteins[loop_proteins$experiment == "RNAbase",]$additional_protein_name)),]$Experiment <- "8oxoG, abasic & RNAbase"
  }
  
  interactions <- read.delim(file = paste0("00_Rawdata/", files[i]),header = T,sep = " ")
  
  cols_tokeep <- colnames(interactions)[-grep(colnames(interactions),pattern = "cooccurence|textmining|combined")]

  cols_tofilter <- cols_tokeep[-grep(pattern = "protein",x = cols_tokeep)]

  filtered_interactions <- interactions %>%
    filter_at(vars(cols_tofilter), any_vars(.>150))
  
  cols_tofilter <- colnames(filtered_interactions)[grep(pattern = "protein",x = colnames(filtered_interactions))]
  
  scores <- filtered_interactions %>%  
    filter_at(vars(cols_tofilter), all_vars(.%in% nodes$STRING_id))
  
  if (nrow(scores) > 0) {
    
    scores$Species <- Species[i]
    
    scores <- scores[,c(11,1:10)]
    
    scores_df <- rbind(scores_df,scores)
    
  }
  
  edges <- filtered_interactions %>%  
    filter_at(vars(cols_tofilter), all_vars(.%in% nodes$STRING_id))
  
  edges <- left_join(edges, nodes, by = c("protein1" = "STRING_id"))

  edges <- left_join(edges, nodes, by = c("protein2" = "STRING_id"))
  
  edges <- edges[,grep(colnames(edges),pattern = "additional")]
  
  colnames(edges) <- c("from", "to")
  
  if (any(grepl(pattern = "9606.ENSP00000406581", x = nodes$STRING_id))) {
    
    nodes <- nodes[-grep(pattern = "9606.ENSP00000406581", x = nodes$STRING_id),] #CHD2 has two string identifiers. Keep interactions for both, exclude multiple node.
    
    
  }
  
  edges <- unique(edges)
  
  nodes <- nodes[nodes$additional_protein_name %in% c(edges$from, edges$to),]
  
  nodes$Border <- as.character(as.numeric(factor(nodes$Experiment, levels = c("8oxoG", 
                                                                                "abasic", 
                                                                                "RNAbase", 
                                                                                "8oxoG & abasic", 
                                                                                "8oxoG & RNAbase", 
                                                                                "abasic & RNAbase", 
                                                                                "8oxoG, abasic & RNAbase"))))
  if (nrow(edges) > 0) {
    
    net <- graph_from_data_frame(d = edges, vertices = nodes)
    
    colors <- c("#000000", "#ffffff","#939598", viridis(n = 6,option = "G")[3:6])
    
    names(colors) <- c("8oxoG", "abasic", "RNAbase", "8oxoG & abasic", "8oxoG & RNAbase", "abasic & RNAbase", "8oxoG, abasic & RNAbase")
    
    colorsborder <- c("#000000", "#000000","#939598", viridis(n = 6,option = "G")[3:6])
    
    names(colorsborder) <- c("1", "2", "3", "4", "5", "6", "7")
    
    colors <- c(colors, colorsborder)
    
    n <- ggnetwork(net)
    
    colors <- colors[names(colors) %in% c(nodes$Experiment, nodes$Border)]
    
    plot <- ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(color = "grey50", size = 2) +
      geom_nodes(aes(color = Border, shape = Repair.protein), size = 20, show.legend = F) +
      geom_nodes(aes(color = Experiment, shape = Repair.protein), size = 19) +
      geom_nodetext_repel(aes(label = name), size = 8,
                          fontface = "bold", box.padding = unit(2, "lines"))+
      scale_color_manual(values = colors)+
      labs(color = "Experiment", shape = "Repair protein")+
      theme_void()+
      ggtitle(Species[i])+
      theme(legend.title = element_text(size = 16, face="bold"))+
      theme(legend.text = element_text(size = 15),
            plot.title = element_text(size=20, face="bold.italic", hjust = 0.5))+
      theme(legend.position = "top", legend.box="vertical")
    
    ggsave(filename = paste0("02_OutputFiles/AdditionalANDProteinID/",
                             gsub(pattern = "-",replacement = "",x = Sys.Date()), "_03_", 
                             gsub(pattern = "\\. | ",replacement = "",x = Species[i]),"_Network.pdf"),plot = plot,width = 15,height = 15)

    loop_plot_list <- list(plot)
    
    plot_list <- c(plot_list, loop_plot_list)
    
    nodes$Species <- Species[i]
    
    nodes <- nodes[,c(7,1:5)]
    
    nodes_df <- rbind(nodes_df,nodes)
    
    if (nrow(edges) > 0) {
      
      edges$Species <- Species[i]
      
      edges <- edges[,c(3,1,2)]
      
      edges_df <- rbind(edges_df, edges)
      
    }
    
  }

}

#### Write result tables ####

write.csv(x = edges_df,
          file = paste0("02_OutputFiles/AdditionalANDProteinID/",
                        gsub(pattern = "-",replacement = "",x = Sys.Date()),
                        "_03_All_NetworkEdges.csv"),
          quote = F,row.names = F)

nodes_df$Experiment <- gsub(pattern = ",", replacement = "",x = nodes_df$Experiment)

write.csv(x = nodes_df,
          file = paste0("02_OutputFiles/AdditionalANDProteinID/",
                        gsub(pattern = "-",replacement = "",x = Sys.Date()),
                        "_03_All_NetworkNodes.csv"),
          quote = F,row.names = F)

scores_df <- scores_df[,c(1:5,7:9)]

write.csv(x = scores_df,
          file = paste0("02_OutputFiles/AdditionalANDProteinID/",
                        gsub(pattern = "-",replacement = "",x = Sys.Date()),
                        "_03_All_InteractionScores.csv"),
          quote = F,row.names = F)