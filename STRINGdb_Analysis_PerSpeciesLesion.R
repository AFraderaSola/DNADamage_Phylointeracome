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
library(scales)
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

lesion <- unique(proteins$experiment)

# STRINGdb (v11.5) interactions per species

files <- list.files(path = "00_Rawdata/",pattern = ".*pro.*txt$")

files <- c(files[1:6],files[6],files[7:11]) # Repeat H. sapiens; one for each cell line.

# Analysis loop per species/lesion.

plot_list <- c()

map_df <- c()

nodes_df <- c()

edges_df <- c()

scores_df <- c()

enrichment_df <- c()

plot_list_species <- c()

for (i in 1:length(species)) {
  
  plot_list_species <- c()
  
  for (j in 1:length(lesion)) {
    
    # Get species STRINGdb
    
    string_db <- STRINGdb$new( version="11.5", species=tx_id[i],score_threshold = 150, network_type= "full", input_directory="")
    
    loop_proteins <- proteins[proteins$Species == species[i],]
    
    loop_proteins <- loop_proteins[loop_proteins$experiment == lesion[j],]
    
    loop_proteins <- string_db$map(loop_proteins,"Protein.IDs",removeUnmappedRows=F)
    
    # Get STRINGdb mapping (first to protein IDs and, if NA, to the additional_protein_name)
    
    if (any(is.na(loop_proteins$STRING_id))) {
      
      loop_proteins_NA <- loop_proteins[is.na(loop_proteins$STRING_id),]
      
      loop_proteins_NA <- loop_proteins_NA[,1:6]
      
      loop_proteins_NA <- string_db$map(loop_proteins_NA,"additional_protein_name",removeUnmappedRows=F)
      
      loop_proteins <- proteins[proteins$Species == species[i],]
      
      loop_proteins <- loop_proteins[loop_proteins$experiment == lesion[j],]
      
      loop_proteins <- string_db$map(loop_proteins,"Protein.IDs",removeUnmappedRows=T)
      
      loop_proteins <- rbind(loop_proteins, loop_proteins_NA)
      
      loop_map_df <- loop_proteins
      
      loop_map_df <- loop_map_df[,c(2,3,1,4,7)]
      
      map_df <- rbind(map_df, loop_map_df)
      
      loop_proteins <- loop_proteins[!is.na(loop_proteins$STRING_id),]
      
      loop_proteins_species <- proteins[proteins$Species == species[i],]
      
    }else{
      
      loop_map_df <- loop_proteins
      
      loop_map_df <- loop_map_df[,c(2,3,1,4,7)]
      
      map_df <- rbind(map_df, loop_map_df)
      
      loop_proteins_species <- proteins[proteins$Species == species[i],]
      
    }
    
    # Draw network from interactions
    
    nodes <- loop_proteins[,c(4,1,7,3,5,6,2)]
    
    nodes <- nodes[,c(1,2,3,5)]
    
    nodes <- unique(nodes)
    
    nodes$Experiment <- lesion[j]
    
    if (length(intersect(loop_proteins[loop_proteins$experiment == lesion[j],]$additional_protein_name,
                         loop_proteins_species[loop_proteins_species$experiment == lesion[-j][1],]$additional_protein_name)) > 0) {
      
      nodes[nodes$additional_protein_name %in% intersect(loop_proteins[loop_proteins$experiment == lesion[j],]$additional_protein_name,
                                                         loop_proteins_species[loop_proteins_species$experiment == lesion[-j][1],]$additional_protein_name),]$Experiment <- paste0(lesion[j], " & ", lesion[-j][1])
    }
    
    if (length(intersect(loop_proteins[loop_proteins$experiment == lesion[j],]$additional_protein_name,
                         loop_proteins_species[loop_proteins_species$experiment == lesion[-j][2],]$additional_protein_name)) > 0) {
      
      nodes[nodes$additional_protein_name %in% intersect(loop_proteins[loop_proteins$experiment == lesion[j],]$additional_protein_name,
                                                         loop_proteins_species[loop_proteins_species$experiment == lesion[-j][2],]$additional_protein_name),]$Experiment <- paste0(lesion[j], " & ", lesion[-j][2])
    }
    
    if (length(intersect(loop_proteins[loop_proteins$experiment == lesion[j],]$additional_protein_name,
                         intersect(loop_proteins_species[loop_proteins_species$experiment == lesion[-j][1],]$additional_protein_name,
                                   loop_proteins_species[loop_proteins_species$experiment == lesion[-j][2],]$additional_protein_name))) > 0) {
      nodes[nodes$additional_protein_name %in% intersect(loop_proteins[loop_proteins$experiment == lesion[j],]$additional_protein_name,
                                                         intersect(loop_proteins_species[loop_proteins_species$experiment == lesion[-j][1],]$additional_protein_name,
                                                                   loop_proteins_species[loop_proteins_species$experiment == lesion[-j][2],]$additional_protein_name)),]$Experiment <- paste0(lesion[j], ", ", lesion[-j][1], " & ", lesion[-j][2])
    }
    
    interactions <- read.delim(file = paste0("00_Rawdata/", files[i]),sep = " ",header = T)
    
    cols_tokeep <- colnames(interactions)[-grep(colnames(interactions),pattern = "cooccurence|textmining|combined")]
    
    cols_tofilter <- cols_tokeep[-grep(pattern = "protein",x = cols_tokeep)]

    filtered_interactions <- interactions %>%
      filter_at(vars(all_of(cols_tofilter)), any_vars(.>150))
    
    cols_tofilter <- colnames(filtered_interactions)[grep(pattern = "protein",x = colnames(filtered_interactions))]
    
    scores <- filtered_interactions %>%  
      filter_at(vars(all_of(cols_tofilter)), all_vars(.%in% nodes$STRING_id))
    
    if (nrow(scores) > 0) {
      
      scores$Species <- Species[i]
      
      scores$Experiment <- lesion[j]
      
      scores <- scores[,c(11,12,1:10)]
      
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
    
    nodes$Border <- as.character(as.numeric(factor(nodes$Experiment, levels = c("8oxoG", "abasic", "RNAbase",
                                                                                paste0(lesion[j], " & ", lesion[-j][1]), 
                                                                                paste0(lesion[j], " & ", lesion[-j][2]),
                                                                                paste0(lesion[j], ", ", lesion[-j][1], " & ", lesion[-j][2])))))
    
    if (nrow(edges) > 0) {
      
      net <- graph_from_data_frame(d = edges, vertices = nodes)
      
      colors <- c("#000000", "#ffffff","#939598", viridis(n = 6,option = "G")[4:6])
      
      names(colors) <- c("8oxoG", "abasic", "RNAbase",
                         paste0(lesion[j], " & ", lesion[-j][1]), 
                         paste0(lesion[j], " & ", lesion[-j][2]),
                         paste0(lesion[j], ", ", lesion[-j][1], " & ", lesion[-j][2]))
      
      colorsborder <- c("#000000", "#000000","#939598", viridis(n = 6,option = "G")[4:6])
      
      names(colorsborder) <- c("1", "2", "3", "4", "5", "6")
      
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
      
      ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                               gsub(pattern = "-",replacement = "",x = Sys.Date()), "_04_", 
                               gsub(pattern = "\\. | ",replacement = "",x = Species[i]),"_",lesion[j],"_Network.pdf"),plot = plot,width = 15,height = 15)

      loop_plot_list <- list(plot)
      
      plot_list <- c(plot_list, loop_plot_list)
      
      plot_list_species <- c(plot_list_species, loop_plot_list)
      
    }
    
    if (nrow(nodes) > 0) {
      
      nodes$Species <- Species[i]
      
      nodes$Lesion <- lesion[j]
      
      nodes <- nodes[,c(7,8,1:5)]
      
      nodes_df <- rbind(nodes_df,nodes)
      
    }
    
    if (nrow(edges) > 0) {
      
      edges$Species <- Species[i]
      
      edges$Lesion <- lesion[j]
      
      edges <- edges[,c(3,4,1,2)]
      
      edges_df <- rbind(edges_df, edges)
      
      
    }
    
    # Get STRINGdb functional analysis
    
    loop_enrichment <- string_db$get_enrichment(loop_proteins$STRING_id)
    
    if (nrow(loop_enrichment) > 0) {
      
      loop_enrichment$Species <- Species[i]
      
      loop_enrichment$Lesion <- lesion[j]
      
      loop_enrichment <- loop_enrichment[,c(11,12,1:10)]
      
      enrichment_df <- rbind(enrichment_df, loop_enrichment)
      
      }
    
  }
  
  if (length(plot_list_species) > 0) {
    
    plot <- ggarrange(plotlist=plot_list_species, ncol = length(plot_list_species), nrow = 1)
    
    ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                             gsub(pattern = "-",replacement = "",x = Sys.Date()),
                             "_04_", 
                             gsub(pattern = "\\. | ",replacement = "",x = Species[i]),
                             "_AllperLesion_Network.pdf"),plot = plot,width = 15*length(plot_list_species),height = 15)
  
  }

}

#### Write result tables ####

write.csv(x = edges_df,
          file = paste0("03_OutputFiles/AdditionalANDProteinID/",
                        gsub(pattern = "-",replacement = "",x = Sys.Date()),
                        "_04_All_perLesion_NetworkEdges.csv"),
          quote = F,row.names = F)

nodes_df$Experiment <- gsub(pattern = ",", replacement = "",x = nodes_df$Experiment)

write.csv(x = nodes_df,
          file = paste0("03_OutputFiles/AdditionalANDProteinID/",
                        gsub(pattern = "-",replacement = "",x = Sys.Date()),
                        "_04_All_perLesion_NetworkNodes.csv"),
          quote = F,row.names = F)

write.csv(x = map_df,
          file = paste0("03_OutputFiles/AdditionalANDProteinID/",
                        gsub(pattern = "-",replacement = "",x = Sys.Date()),
                        "_01_All_perLesion_MappingSTRINGdb.csv"),
          quote = F,row.names = F)

scores_df <- scores_df[,c(1:6,8:10)]

write.csv(x = scores_df,
          file = paste0("03_OutputFiles/AdditionalANDProteinID/",
                        gsub(pattern = "-",replacement = "",x = Sys.Date()),
                        "_04_All_perLesion_InteractionScores.csv"),
          quote = F,row.names = F)

enrichment_df <- enrichment_df[enrichment_df$category %in% c("Process", "Function", "Component", "KEGG", "Pfam"),]

enrichment_df$preferredNames <- gsub(pattern = ",",replacement = ";",x = enrichment_df$preferredNames)

enrichment_df$inputGenes <- gsub(pattern = ",",replacement = ";",x = enrichment_df$inputGenes)

enrichment_df$description <- gsub(pattern = ",",replacement = " ",x = enrichment_df$description)

enrichment_df_tosave <- enrichment_df[enrichment_df$category %in% c("Process","KEGG"),]

write.csv(x = enrichment_df_tosave,
          file = paste0("03_OutputFiles/AdditionalANDProteinID/",
                        gsub(pattern = "-",replacement = "",x = Sys.Date()),
                        "_02_All_Enrichment.csv"),
          quote = F,row.names = F)

#### STRINGdb Mapping and number of interactions ####

# Mapping

map_df$Species <- factor(x = map_df$Species, levels = c("E. coli",
                                                          "B. subtilis",
                                                          "H. salinarum",
                                                          "T. brucei",
                                                          "T. thermophila",
                                                          "S. cerevisiae",
                                                          "S. pombe",
                                                          "C. elegans",
                                                          "H. sapiens (HeLa)",
                                                          "H. sapiens (HEK293)",
                                                          "Z. mays",
                                                          "A. thaliana"))

map_df$experiment <- factor(x = map_df$experiment, levels = c("8oxoG", "abasic", "RNAbase"))

map_df$Fill <- is.na(map_df$STRING_id)

map_df <- map_df %>% group_by(Species,experiment, Fill) %>% tally()

map_df$Fill <- as.character(map_df$Fill)

map_df[map_df$Fill == "TRUE",]$Fill <- "Not included"

map_df[map_df$Fill == "FALSE",]$Fill <- "Included"

plot <- ggplot(map_df, aes(fill = Fill,y = n , x= Species)) +
  facet_wrap(~experiment)+
  geom_bar(position="stack", stat="identity", size = 1.5)+
  labs(fill = "STRINGdb")+
  scale_fill_brewer(palette = "Paired",direction = -1)+
  theme_bw()+
  theme(axis.text.x=element_text(size=16,hjust = 1,angle = 45),
        axis.text.y=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=20, face="bold", hjust = 0.5))+
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 16,face="bold"))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  theme(strip.text = element_text(size=16))+
  theme(axis.text.x = element_text(face = "italic"))+
  ylab("Number of proteins")

ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                         gsub(pattern = "-",replacement = "",x = Sys.Date()),
                         "_01_All_perLesion_STRINGdbInclusion.pdf"),plot = plot,width = 15,height = 7.5)

# Number of interactions

n_interactions <- edges_df %>% group_by(Species,Lesion) %>% tally()
n_interactions$n <- n_interactions$n / 2
colnames(n_interactions)[2] <- "Experiment"



n_interactions_sal <- data.frame(Species = rep("H. salinarum",3),
                                 Experiment = c("8oxoG", "abasic", "RNAbase"),
                                 n = rep(0, 3))

n_interactions <- rbind(n_interactions, n_interactions_sal)

n_interactions_species <- unique(n_interactions$Species)

n_interactions_experiment <- as.character(unique(n_interactions$Experiment))

for (i in 1:length(n_interactions_species)) {
  
  loop_df <- n_interactions[n_interactions$Species == n_interactions_species[i],]
  
  if (length(grep(pattern = paste0(loop_df$Experiment,collapse = "|"), n_interactions_experiment)) < 3) {
    
    new_loop_df_experiment <- n_interactions_experiment[-grep(pattern = paste0(loop_df$Experiment,collapse = "|"), n_interactions_experiment)]
    
    new_loop_df <- data.frame(Species = rep(n_interactions_species[i], length(new_loop_df_experiment)),
                              Experiment = new_loop_df_experiment,
                              n = rep(0, length(new_loop_df_experiment)))
    
    n_interactions <- rbind(n_interactions, new_loop_df)
    
  }
}

n_interactions <- n_interactions %>% arrange(Species, Experiment)

n_interactions$Experiment <- factor(x = n_interactions$Experiment, levels = c("8oxoG", "abasic", "RNAbase"))

n_interactions$Species <- factor(x = n_interactions$Species, levels = c("E. coli",
                                                                        "B. subtilis",
                                                                        "H. salinarum",
                                                                        "T. brucei",
                                                                        "T. thermophila",
                                                                        "S. cerevisiae",
                                                                        "S. pombe",
                                                                        "C. elegans",
                                                                        "H. sapiens (HeLa)",
                                                                        "H. sapiens (HEK293)",
                                                                        "Z. mays",
                                                                        "A. thaliana"))

colors_fill <- c("#000000", "#ffffff","#939598")

colors_colour <- c("#000000", "#000000", "#939598")

write.csv(x = n_interactions,
          file = paste0("03_OutputFiles/AdditionalANDProteinID/",
                        gsub(pattern = "-",replacement = "",x = Sys.Date()),
                        "_01_All_perLesion_STRINGdbInteractions.csv"),
          quote = F,row.names = F)

plot <- ggplot(n_interactions, aes(y = n , x= Species, fill = Experiment, colour = Experiment)) +
  geom_bar(position=position_dodge(.7), stat="identity", width=.6)+
  labs(fill = "Experiment")+
  scale_fill_manual(values = colors_fill)+
  scale_color_manual(values = colors_colour, guide = "none")+
  theme_minimal()+
  theme(axis.text.x=element_text(size=16,hjust = 1,angle = 45),
        axis.text.y=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=20, face="bold", hjust = 0.5))+
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 16,face="bold"))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  theme(strip.text = element_text(size=16))+
  theme(axis.text.x = element_text(face = "italic"))+
  ylab("Number of interactions")

ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                         gsub(pattern = "-",replacement = "",x = Sys.Date()),
                         "_01_All_perLesion_STRINGdbInteractions.pdf"),plot = plot,width = 15,height = 7.5)


########################################
#### STRINGdb Functional Enrichment ####
########################################

#### Barplots per species ####

species <- unique(enrichment_df$Species)

for (i in 1:length(species)) {
  
  mf <- enrichment_df[enrichment_df$category == "Function",]
  
  mf <- mf[mf$Species == species[i],]
  
  mf <- mf %>% group_by(Species) %>% slice_max(fdr,n = 5)
  
  if (nrow(mf)>0) {
    
    plot <- ggplot(mf, aes(fill = fdr,y = description  , x= number_of_genes)) +
      facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
      geom_bar(position="stack", stat="identity", size = 1.5)+
      # facet_wrap(~Species, scales = "free_y", ncol = 1)+
      labs(fill = "FDR")+
      scale_fill_viridis_c(option = "G")+
      theme_bw()+
      ggtitle(species[i])+
      theme(axis.text.x=element_text(size=16,hjust = 1),
            axis.text.y=element_text(size=16),
            axis.title=element_text(size=18,face="bold"),
            plot.title = element_text(size=20, face="bold.italic"))+
      theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
            legend.title = element_text(size = 16,face="bold"))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position = "top")+
      theme(strip.text = element_text(size=16))+
      xlab("Number of genes")+
      ylab("GO Molecular Function term")
    
    ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                             gsub(pattern = "-",replacement = "",x = Sys.Date()),
                             "_02_", 
                             gsub(pattern = "\\. | ",replacement = "",x = species[i]),
                             "_GOMF_Barplot.pdf"),plot = plot,width = 15,height = 15)
    
  }
  
  bp <- enrichment_df[enrichment_df$category == "Process",]
  
  bp <- bp[bp$Species == species[i],]
  
  bp <- bp %>% group_by(Species) %>% slice_max(fdr,n = 5)
  
  if (nrow(bp)>0) {
    
    plot <- ggplot(bp, aes(fill = fdr,y = description  , x= number_of_genes)) +
      geom_bar(position="stack", stat="identity", size = 1.5)+
      facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
      # facet_wrap(~Species, scales = "free_y", ncol = 1)+
      labs(fill = "FDR")+
      scale_fill_viridis_c(option = "G")+
      theme_bw()+
      ggtitle(species[i])+
      theme(axis.text.x=element_text(size=16,hjust = 1),
            axis.text.y=element_text(size=16),
            axis.title=element_text(size=18,face="bold"),
            plot.title = element_text(size=20, face="bold.italic"))+
      theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
            legend.title = element_text(size = 16,face="bold"))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position = "top")+
      theme(strip.text = element_text(size=16))+
      xlab("Number of genes")+
      ylab("GO Biological Process term")
    
    ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                             gsub(pattern = "-",replacement = "",x = Sys.Date()),
                             "_02_", 
                             gsub(pattern = "\\. | ",replacement = "",x = species[i]),
                             "_GOBP_Barplot.pdf"),plot = plot,width = 15,height = 15)
  }
  
  
  
  cc <- enrichment_df[enrichment_df$category == "Component",]
  
  cc <- cc[cc$Species == species[i],]
  
  cc <- cc %>% group_by(Species) %>% slice_max(fdr,n = 5)
  
  if (nrow(cc)>0) {
    
    plot <- ggplot(cc, aes(fill = fdr,y = description  , x= number_of_genes)) +
      geom_bar(position="stack", stat="identity", size = 1.5)+
      facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
      # facet_wrap(~Species, scales = "free_y", ncol = 1)+
      labs(fill = "FDR")+
      scale_fill_viridis_c(option = "G")+
      theme_bw()+
      ggtitle(species[i])+
      theme(axis.text.x=element_text(size=16,hjust = 1),
            axis.text.y=element_text(size=16),
            axis.title=element_text(size=18,face="bold"),
            plot.title = element_text(size=20, face="bold.italic"))+
      theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
            legend.title = element_text(size = 16,face="bold"))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position = "top")+
      theme(strip.text = element_text(size=16))+
      xlab("Number of genes")+
      ylab("GO Cellular Component term")
    
    ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                             gsub(pattern = "-",replacement = "",x = Sys.Date()),
                             "_02_", 
                             gsub(pattern = "\\. | ",replacement = "",x = species[i]),
                             "_GOCC_Barplot.pdf"),plot = plot,width = 15,height = 15)
    
  }
  
  
  
  kegg <- enrichment_df[enrichment_df$category == "KEGG",]
  
  kegg <- kegg[kegg$Species == species[i],]
  
  if (nrow(kegg)>0) {
   
    plot <- ggplot(kegg, aes(fill = fdr,y = description  , x= number_of_genes)) +
      geom_bar(position="stack", stat="identity", size = 1.5)+
      facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
      # facet_wrap(~Species, scales = "free_y", ncol = 1)+
      labs(fill = "FDR")+
      scale_fill_viridis_c(option = "G")+
      theme_bw()+
      ggtitle(species[i])+
      theme(axis.text.x=element_text(size=16,hjust = 1),
            axis.text.y=element_text(size=16),
            axis.title=element_text(size=18,face="bold"),
            plot.title = element_text(size=20, face="bold.italic"))+
      theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
            legend.title = element_text(size = 16,face="bold"))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position = "top")+
      theme(strip.text = element_text(size=16))+
      xlab("Number of genes")+
      ylab("KEGG term")
    
    ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                             gsub(pattern = "-",replacement = "",x = Sys.Date()),
                             "_02_", 
                             gsub(pattern = "\\. | ",replacement = "",x = species[i]),
                             "_KEGG_Barplot.pdf"),plot = plot,width = 15,height = 15)
    
  }
  
  
  
  pfam <- enrichment_df[enrichment_df$category == "Pfam",]
  
  pfam <- pfam[pfam$Species == species[i],]
  
  if (nrow(pfam)>0) {
    
    plot <- ggplot(pfam, aes(fill = fdr,y = description  , x= number_of_genes)) +
      geom_bar(position="stack", stat="identity", size = 1.5)+
      facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
      # facet_wrap(~Species, scales = "free_y", ncol = 1)+
      labs(fill = "FDR")+
      scale_fill_viridis_c(option = "G")+
      theme_bw()+
      ggtitle(species[i])+
      theme(axis.text.x=element_text(size=16,hjust = 1),
            axis.text.y=element_text(size=16),
            axis.title=element_text(size=18,face="bold"),
            plot.title = element_text(size=20, face="bold.italic"))+
      theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
            legend.title = element_text(size = 16,face="bold"))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position = "top")+
      theme(strip.text = element_text(size=16))+
      xlab("Number of genes")+
      ylab("Pfam term")
    
    ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                             gsub(pattern = "-",replacement = "",x = Sys.Date()),
                             "_02_", 
                             gsub(pattern = "\\. | ",replacement = "",x = species[i]),
                             "_Pfam_Barplot.pdf"),plot = plot,width = 15,height = 15)
    
  }
  
}

#### TOPN DotPlots and HeatMaps ####

top <- c(5,3)

for (j in 1:length(top)) {
  
  bp_dp <- enrichment_df[enrichment_df$category == "Process",]
  
  bp_dp <- bp_dp %>% group_by(Species, Lesion) %>% slice_max(fdr,n = top[j])
  
  bp_dp$Ratio <- bp_dp$number_of_genes/bp_dp$number_of_genes_in_background 
  
  bp_dp$Lesion <- factor(x = bp_dp$Lesion, levels = c("8oxoG", "abasic", "RNAbase"))
  
  bp_dp$Species <- factor(x = bp_dp$Species, levels = c("E. coli",
                                                        "B. subtilis",
                                                        "H. salinarum",
                                                        "T. brucei",
                                                        "T. thermophila",
                                                        "S. cerevisiae",
                                                        "S. pombe",
                                                        "C. elegans",
                                                        "H. sapiens (HeLa)",
                                                        "H. sapiens (HEK293)",
                                                        "Z. mays",
                                                        "A. thaliana"))
  
  plot <- ggplot(bp_dp, aes(fill = fdr,y = description, size = Ratio, x= Species)) +
    geom_point()+
    facet_wrap(~Lesion, scales = "free_y", ncol = 1)+
    # facet_wrap(~Species, scales = "free_y", ncol = 1)+
    labs(fill = "FDR",size = "Gene Ratio")+
    scale_fill_viridis_c(option = "G")+
    scale_size_binned()+
    theme_bw()+
    # ggtitle(species[i])+
    theme(axis.text.x=element_text(size=16,hjust = 1, angle = 45),
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=18,face="bold"),
          plot.title = element_text(size=20, face="bold"))+
    theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
          legend.title = element_text(size = 16,face="bold"))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.position = "top")+
    theme(strip.text = element_text(size=16))+
    theme(axis.text.x = element_text(face = "italic"))+
    xlab("Species")+
    ylab("GO Biological Process term")
  
  ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                           gsub(pattern = "-",replacement = "",x = Sys.Date()),
                           "_02_All_top",
                           top[j], "_GOBP_Dotplot.pdf"),plot = plot,width = 15,height = 15)
  
  bp_dp$SpeciesandLesion <- paste0(bp_dp$Species, " - ", bp_dp$Lesion)
  
  plot <- ggplot(bp_dp, aes(fill = fdr,y = description, size = Ratio, x= SpeciesandLesion)) +
    geom_point()+
    # facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
    # facet_wrap(~Species, scales = "free_y", ncol = 1)+
    labs(fill = "FDR",size = "Gene Ratio")+
    scale_fill_viridis_c(option = "G")+
    scale_size_binned()+
    theme_bw()+
    # ggtitle(species[i])+
    theme(axis.text.x=element_text(size=16,hjust = 1, angle = 45),
          axis.text.y=element_text(size=16),
          axis.title=element_text(size=18,face="bold"),
          plot.title = element_text(size=20, face="bold"))+
    theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
          legend.title = element_text(size = 16,face="bold"))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.position = "top")+
    theme(strip.text = element_text(size=16))+
    theme(axis.text.x = element_text(face = "italic"))+
    xlab("Species and Lesion")+
    ylab("GO Biological Process term")
  
  ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                           gsub(pattern = "-",replacement = "",x = Sys.Date()),
                           "_02_AllnoFacet_top",
                           top[j], "_GOBP_Dotplot.pdf"),plot = plot,width = 15,height = 15)
  
  bp_hm <- enrichment_df[enrichment_df$category == "Process",]
  
  bp_hm <- bp_hm %>% group_by(Species, Lesion) %>% slice_max(fdr,n = top[j])
  
  bp_hm$Ratio <- bp_hm$number_of_genes/bp_hm$number_of_genes_in_background 
  
  bp_hm$Lesion <- factor(x = bp_hm$Lesion, levels = c("8oxoG", "abasic", "RNAbase"))
  
  bp_hm$Species <- factor(x = bp_hm$Species, levels = c("E. coli",
                                                        "B. subtilis",
                                                        "H. salinarum",
                                                        "T. brucei",
                                                        "T. thermophila",
                                                        "S. cerevisiae",
                                                        "S. pombe",
                                                        "C. elegans",
                                                        "H. sapiens (HeLa)",
                                                        "H. sapiens (HEK293)",
                                                        "Z. mays",
                                                        "A. thaliana"))
  
  bp_hm$SpeciesandLesion <- paste0(bp_dp$Species, " - ", bp_dp$Lesion)
  
  bp_hm$SpeciesandLesion <- factor(x = bp_hm$SpeciesandLesion, levels = c("E. coli - 8oxoG",
                                                                          "B. subtilis - 8oxoG",
                                                                          "H. salinarum - 8oxoG",
                                                                          "T. brucei - 8oxoG",
                                                                          "T. thermophila - 8oxoG",
                                                                          "S. cerevisiae - 8oxoG",
                                                                          "S. pombe - 8oxoG",
                                                                          "C. elegans - 8oxoG",
                                                                          "H. sapiens (HeLa) - 8oxoG",
                                                                          "H. sapiens (HEK293) - 8oxoG",
                                                                          "Z. mays - 8oxoG",
                                                                          "A. thaliana - 8oxoG",
                                                                          "E. coli - abasic",
                                                                          "B. subtilis - abasic",
                                                                          "H. salinarum - abasic",
                                                                          "T. brucei - abasic",
                                                                          "T. thermophila - abasic",
                                                                          "S. cerevisiae - abasic",
                                                                          "S. pombe - abasic",
                                                                          "C. elegans - abasic",
                                                                          "H. sapiens (HeLa) - abasic",
                                                                          "H. sapiens (HEK293) - abasic",
                                                                          "Z. mays - abasic",
                                                                          "A. thaliana - abasic",
                                                                          "E. coli - RNAbase",
                                                                          "B. subtilis - RNAbase",
                                                                          "H. salinarum - RNAbase",
                                                                          "T. brucei - RNAbase",
                                                                          "T. thermophila - RNAbase",
                                                                          "S. cerevisiae - RNAbase",
                                                                          "S. pombe - RNAbase",
                                                                          "C. elegans - RNAbase",
                                                                          "H. sapiens (HeLa) - RNAbase",
                                                                          "H. sapiens (HEK293) - RNAbase",
                                                                          "Z. mays - RNAbase",
                                                                          "A. thaliana - RNAbase"))
  
  bp_hm <- bp_hm %>% arrange(SpeciesandLesion)
  
  bp_hm <- bp_hm[,c(14,12,13)]
  
  filter <- unique(as.character(bp_hm$SpeciesandLesion))
  
  bp_df_hm <- bp_hm[bp_hm$SpeciesandLesion == filter[1],]
  
  colnames(bp_df_hm)[3] <- paste0(unique(bp_df_hm$SpeciesandLesion), "_GeneRatio")
  
  bp_df_hm <- bp_df_hm[,2:3]
  
  for (i in 2:length(filter)) {
    
    loop_bp_hm <- bp_hm[bp_hm$SpeciesandLesion == filter[i],]
    
    colnames(loop_bp_hm)[3] <- paste0(unique(loop_bp_hm$SpeciesandLesion), "_GeneRatio")
    
    loop_bp_hm <- loop_bp_hm[,2:3]
    
    bp_df_hm <- full_join(x = bp_df_hm,y = loop_bp_hm, by = "description")
    
    
  }
  
  bp_df_hm <- as.data.frame(bp_df_hm)
  
  rownames(bp_df_hm) <- bp_df_hm[,1]
  
  bp_df_hm <- bp_df_hm[,c(2:ncol(bp_df_hm))]
  
  bp_df_hm <- as.matrix(bp_df_hm)
  
  colours <- c(viridis(n = 100, option = "G",direction = -1))
  
  Experiment <- c(rep("8oxoG", length(grep(pattern = "8oxoG",x = colnames(bp_df_hm)))),
                  rep("abasic", length(grep(pattern = "abasic",x = colnames(bp_df_hm)))),
                  rep("RNAbase", length(grep(pattern = "RNAbase",x = colnames(bp_df_hm)))))
  
  df_col <- data.frame(Experiment)
  
  rownames(df_col) <- colnames(bp_df_hm)
  
  Experiment <- c("#000000", "#ffffff","#939598")
  
  names(Experiment) <- c("8oxoG", "abasic", "RNAbase")
  
  anno_colors <- list(Experiment = Experiment)
  
  filename <- paste0("03_OutputFiles/AdditionalANDProteinID/",
                     gsub(pattern = "-",replacement = "",x = Sys.Date()),
                     "_02_AllnoFacet_top",
                     top[j], "_GOBP_HeatMap.pdf")
  
  pheatmap(bp_df_hm,
           col=colours,
           cluster_cols = F,
           cluster_rows = F,
           show_rownames = T,
           show_colnames = T,
           labels_col = str_to_title(gsub(pattern = " - .*",replacement = "",x = colnames(bp_df_hm))),
           annotation_col = df_col,
           # annotation_row = df_row,
           annotation_colors = anno_colors,
           fontsize = 16,border_color = "#b4b4b4",
           height = 15,width = 15,filename = filename
  )
  
  mf_dp <- enrichment_df[enrichment_df$category == "Function",]
  
  mf_dp <- mf_dp %>% group_by(Species, Lesion) %>% slice_max(fdr,n = top[j])
  
  mf_dp$Ratio <- mf_dp$number_of_genes/mf_dp$number_of_genes_in_background 
  
  mf_dp$Lesion <- factor(x = mf_dp$Lesion, levels = c("8oxoG", "abasic", "RNAbase"))
  
  mf_dp$Species <- factor(x = mf_dp$Species, levels = c("E. coli",
                                                        "B. subtilis",
                                                        "H. salinarum",
                                                        "T. brucei",
                                                        "T. thermophila",
                                                        "S. cerevisiae",
                                                        "S. pombe",
                                                        "C. elegans",
                                                        "H. sapiens (HeLa)",
                                                        "H. sapiens (HEK293)",
                                                        "Z. mays",
                                                        "A. thaliana"))
  
  plot <- ggplot(mf_dp, aes(fill = fdr,y = description, size = Ratio, x= Species)) +
    geom_point()+
    facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
    # facet_wrap(~Species, scales = "free_y", ncol = 1)+
    labs(fill = "FDR",size = "Gene Ratio")+
    scale_fill_viridis_c(option = "G")+
    scale_size_binned()+
    theme_bw()+
    # ggtitle(species[i])+
    theme(axis.text.x=element_text(size=16,hjust = 1, angle = 45),
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=18,face="bold"),
          plot.title = element_text(size=20, face="bold"))+
    theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
          legend.title = element_text(size = 16,face="bold"))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.position = "top")+
    theme(strip.text = element_text(size=16))+
    theme(axis.text.x = element_text(face = "italic"))+
    xlab("Species")+
    ylab("GO Molecular Function term")
  
  ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                           gsub(pattern = "-",replacement = "",x = Sys.Date()),
                           "_02_All_top",
                           top[j], "_GOMF_Dotplot.pdf"),plot = plot,width = 15,height = 15)
  
  mf_dp$SpeciesandLesion <- paste0(mf_dp$Species, " - ", mf_dp$Lesion)
  
  plot <- ggplot(mf_dp, aes(fill = fdr,y = description, size = Ratio, x= SpeciesandLesion)) +
    geom_point()+
    # facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
    # facet_wrap(~Species, scales = "free_y", ncol = 1)+
    labs(fill = "FDR",size = "Gene Ratio")+
    scale_fill_viridis_c(option = "G")+
    scale_size_binned()+
    theme_bw()+
    # ggtitle(species[i])+
    theme(axis.text.x=element_text(size=16,hjust = 1, angle = 45),
          axis.text.y=element_text(size=16),
          axis.title=element_text(size=18,face="bold"),
          plot.title = element_text(size=20, face="bold"))+
    theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
          legend.title = element_text(size = 16,face="bold"))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.position = "top")+
    theme(strip.text = element_text(size=16))+
    theme(axis.text.x = element_text(face = "italic"))+
    xlab("Species and Lesion")+
    ylab("GO Molecular Function term")
  
  ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                           gsub(pattern = "-",replacement = "",x = Sys.Date()),
                           "_02_AllnoFacet_top",
                           top[j], "_GOMF_Dotplot.pdf"),plot = plot,width = 15,height = 15)
  
  mf_hm <- enrichment_df[enrichment_df$category == "Function",]
  
  mf_hm <- mf_hm %>% group_by(Species, Lesion) %>% slice_max(fdr,n = top[j])
  
  mf_hm$Ratio <- mf_hm$number_of_genes/mf_hm$number_of_genes_in_background 
  
  mf_hm$Lesion <- factor(x = mf_hm$Lesion, levels = c("8oxoG", "abasic", "RNAbase"))
  
  mf_hm$Species <- factor(x = mf_hm$Species, levels = c("E. coli",
                                                        "B. subtilis",
                                                        "H. salinarum",
                                                        "T. brucei",
                                                        "T. thermophila",
                                                        "S. cerevisiae",
                                                        "S. pombe",
                                                        "C. elegans",
                                                        "H. sapiens (HeLa)",
                                                        "H. sapiens (HEK293)",
                                                        "Z. mays",
                                                        "A. thaliana"))
  
  mf_hm$SpeciesandLesion <- paste0(mf_dp$Species, " - ", mf_dp$Lesion)
  
  mf_hm$SpeciesandLesion <- factor(x = mf_hm$SpeciesandLesion, levels = c("E. coli - 8oxoG",
                                                                          "B. subtilis - 8oxoG",
                                                                          "H. salinarum - 8oxoG",
                                                                          "T. brucei - 8oxoG",
                                                                          "T. thermophila - 8oxoG",
                                                                          "S. cerevisiae - 8oxoG",
                                                                          "S. pombe - 8oxoG",
                                                                          "C. elegans - 8oxoG",
                                                                          "H. sapiens (HeLa) - 8oxoG",
                                                                          "H. sapiens (HEK293) - 8oxoG",
                                                                          "Z. mays - 8oxoG",
                                                                          "A. thaliana - 8oxoG",
                                                                          "E. coli - abasic",
                                                                          "B. subtilis - abasic",
                                                                          "H. salinarum - abasic",
                                                                          "T. brucei - abasic",
                                                                          "T. thermophila - abasic",
                                                                          "S. cerevisiae - abasic",
                                                                          "S. pombe - abasic",
                                                                          "C. elegans - abasic",
                                                                          "H. sapiens (HeLa) - abasic",
                                                                          "H. sapiens (HEK293) - abasic",
                                                                          "Z. mays - abasic",
                                                                          "A. thaliana - abasic",
                                                                          "E. coli - RNAbase",
                                                                          "B. subtilis - RNAbase",
                                                                          "H. salinarum - RNAbase",
                                                                          "T. brucei - RNAbase",
                                                                          "T. thermophila - RNAbase",
                                                                          "S. cerevisiae - RNAbase",
                                                                          "S. pombe - RNAbase",
                                                                          "C. elegans - RNAbase",
                                                                          "H. sapiens (HeLa) - RNAbase",
                                                                          "H. sapiens (HEK293) - RNAbase",
                                                                          "Z. mays - RNAbase",
                                                                          "A. thaliana - RNAbase"))
  
  mf_hm <- mf_hm %>% arrange(SpeciesandLesion)
  
  mf_hm <- mf_hm[,c(14,12,13)]
  
  filter <- unique(as.character(mf_hm$SpeciesandLesion))
  
  mf_df_hm <- mf_hm[mf_hm$SpeciesandLesion == filter[1],]
  
  colnames(mf_df_hm)[3] <- paste0(unique(mf_df_hm$SpeciesandLesion), "_GeneRatio")
  
  mf_df_hm <- mf_df_hm[,2:3]
  
  i <- 3
  
  for (i in 2:length(filter)) {
    
    loop_mf_hm <- mf_hm[mf_hm$SpeciesandLesion == filter[i],]
    
    colnames(loop_mf_hm)[3] <- paste0(unique(loop_mf_hm$SpeciesandLesion), "_GeneRatio")
    
    loop_mf_hm <- loop_mf_hm[,2:3]
    
    mf_df_hm <- full_join(x = mf_df_hm,y = loop_mf_hm, by = "description")
    
    
  }
  
  mf_df_hm <- as.data.frame(mf_df_hm)
  
  rownames(mf_df_hm) <- mf_df_hm[,1]
  
  mf_df_hm <- mf_df_hm[,c(2:ncol(mf_df_hm))]
  
  mf_df_hm <- as.matrix(mf_df_hm)
  
  colours <- c(viridis(n = 100, option = "G",direction = -1))
  
  Experiment <- c(rep("8oxoG", length(grep(pattern = "8oxoG",x = colnames(mf_df_hm)))),
                  rep("abasic", length(grep(pattern = "abasic",x = colnames(mf_df_hm)))),
                  rep("RNAbase", length(grep(pattern = "RNAbase",x = colnames(mf_df_hm)))))
  
  df_col <- data.frame(Experiment)
  
  rownames(df_col) <- colnames(mf_df_hm)
  
  Experiment <- c("#000000", "#ffffff","#939598")
  
  names(Experiment) <- c("8oxoG", "abasic", "RNAbase")
  
  anno_colors <- list(Experiment = Experiment)
  
  filename <- paste0("03_OutputFiles/AdditionalANDProteinID/",
                     gsub(pattern = "-",replacement = "",x = Sys.Date()),
                     "_02_AllnoFacet_top",
                     top[j], "_GOMF_HeatMap.pdf")
  
  pheatmap(mf_df_hm,
           col=colours,
           cluster_cols = F,
           cluster_rows = F,
           show_rownames = T,
           show_colnames = T,
           labels_col = str_to_title(gsub(pattern = " - .*",replacement = "",x = colnames(mf_df_hm))),
           annotation_col = df_col,
           # annotation_row = df_row,
           annotation_colors = anno_colors,
           fontsize = 16,border_color = "#b4b4b4",
           height = 15,width = 15,filename = filename
  )
  
  cc_dp <- enrichment_df[enrichment_df$category == "Component",]
  
  cc_dp <- cc_dp %>% group_by(Species, Lesion) %>% slice_max(fdr,n = top[j])
  
  cc_dp$Ratio <- cc_dp$number_of_genes/cc_dp$number_of_genes_in_background 
  
  cc_dp$Lesion <- factor(x = cc_dp$Lesion, levels = c("8oxoG", "abasic", "RNAbase"))
  
  cc_dp$Species <- factor(x = cc_dp$Species, levels = c("E. coli",
                                                        "B. subtilis",
                                                        "H. salinarum",
                                                        "T. brucei",
                                                        "T. thermophila",
                                                        "S. cerevisiae",
                                                        "S. pombe",
                                                        "C. elegans",
                                                        "H. sapiens (HeLa)",
                                                        "H. sapiens (HEK293)",
                                                        "Z. mays",
                                                        "A. thaliana"))
  
  plot <- ggplot(cc_dp, aes(fill = fdr,y = description, size = Ratio, x= Species)) +
    geom_point()+
    facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
    # facet_wrap(~Species, scales = "free_y", ncol = 1)+
    labs(fill = "FDR",size = "Gene Ratio")+
    scale_fill_viridis_c(option = "G")+
    scale_size_binned()+
    theme_bw()+
    # ggtitle(species[i])+
    theme(axis.text.x=element_text(size=16,hjust = 1, angle = 45),
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=18,face="bold"),
          plot.title = element_text(size=20, face="bold"))+
    theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
          legend.title = element_text(size = 16,face="bold"))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.position = "top")+
    theme(strip.text = element_text(size=16))+
    theme(axis.text.x = element_text(face = "italic"))+
    xlab("Species")+
    ylab("GO Cellular Component term")
  
  ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                           gsub(pattern = "-",replacement = "",x = Sys.Date()),
                           "_02_All_top",
                           top[j], "_GOCC_Dotplot.pdf"),plot = plot,width = 15,height = 15)
  
  cc_dp$SpeciesandLesion <- paste0(cc_dp$Species, " - ", cc_dp$Lesion)
  
  plot <- ggplot(cc_dp, aes(fill = fdr,y = description, size = Ratio, x= SpeciesandLesion)) +
    geom_point()+
    # facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
    # facet_wrap(~Species, scales = "free_y", ncol = 1)+
    labs(fill = "FDR",size = "Gene Ratio")+
    scale_fill_viridis_c(option = "G")+
    scale_size_binned()+
    theme_bw()+
    # ggtitle(species[i])+
    theme(axis.text.x=element_text(size=16,hjust = 1, angle = 45),
          axis.text.y=element_text(size=16),
          axis.title=element_text(size=18,face="bold"),
          plot.title = element_text(size=20, face="bold"))+
    theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
          legend.title = element_text(size = 16,face="bold"))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.position = "top")+
    theme(strip.text = element_text(size=16))+
    theme(axis.text.x = element_text(face = "italic"))+
    xlab("Species and Lesion")+
    ylab("GO Cellular Component term")
  
  ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                           gsub(pattern = "-",replacement = "",x = Sys.Date()),
                           "_02_AllnoFacet_top",
                           top[j], "_GOCC_Dotplot.pdf"),plot = plot,width = 15,height = 15)
  
  cc_hm <- enrichment_df[enrichment_df$category == "Component",]
  
  cc_hm <- cc_hm %>% group_by(Species, Lesion) %>% slice_max(fdr,n = top[j])
  
  cc_hm$Ratio <- cc_hm$number_of_genes/cc_hm$number_of_genes_in_background 
  
  cc_hm$Lesion <- factor(x = cc_hm$Lesion, levels = c("8oxoG", "abasic", "RNAbase"))
  
  cc_hm$Species <- factor(x = cc_hm$Species, levels = c("E. coli",
                                                        "B. subtilis",
                                                        "H. salinarum",
                                                        "T. brucei",
                                                        "T. thermophila",
                                                        "S. cerevisiae",
                                                        "S. pombe",
                                                        "C. elegans",
                                                        "H. sapiens (HeLa)",
                                                        "H. sapiens (HEK293)",
                                                        "Z. mays",
                                                        "A. thaliana"))
  
  cc_hm$SpeciesandLesion <- paste0(cc_dp$Species, " - ", cc_dp$Lesion)
  
  cc_hm$SpeciesandLesion <- factor(x = cc_hm$SpeciesandLesion, levels = c("E. coli - 8oxoG",
                                                                          "B. subtilis - 8oxoG",
                                                                          "H. salinarum - 8oxoG",
                                                                          "T. brucei - 8oxoG",
                                                                          "T. thermophila - 8oxoG",
                                                                          "S. cerevisiae - 8oxoG",
                                                                          "S. pombe - 8oxoG",
                                                                          "C. elegans - 8oxoG",
                                                                          "H. sapiens (HeLa) - 8oxoG",
                                                                          "H. sapiens (HEK293) - 8oxoG",
                                                                          "Z. mays - 8oxoG",
                                                                          "A. thaliana - 8oxoG",
                                                                          "E. coli - abasic",
                                                                          "B. subtilis - abasic",
                                                                          "H. salinarum - abasic",
                                                                          "T. brucei - abasic",
                                                                          "T. thermophila - abasic",
                                                                          "S. cerevisiae - abasic",
                                                                          "S. pombe - abasic",
                                                                          "C. elegans - abasic",
                                                                          "H. sapiens (HeLa) - abasic",
                                                                          "H. sapiens (HEK293) - abasic",
                                                                          "Z. mays - abasic",
                                                                          "A. thaliana - abasic",
                                                                          "E. coli - RNAbase",
                                                                          "B. subtilis - RNAbase",
                                                                          "H. salinarum - RNAbase",
                                                                          "T. brucei - RNAbase",
                                                                          "T. thermophila - RNAbase",
                                                                          "S. cerevisiae - RNAbase",
                                                                          "S. pombe - RNAbase",
                                                                          "C. elegans - RNAbase",
                                                                          "H. sapiens (HeLa) - RNAbase",
                                                                          "H. sapiens (HEK293) - RNAbase",
                                                                          "Z. mays - RNAbase",
                                                                          "A. thaliana - RNAbase"))
  
  cc_hm <- cc_hm %>% arrange(SpeciesandLesion)
  
  cc_hm <- cc_hm[,c(14,12,13)]
  
  filter <- unique(as.character(cc_hm$SpeciesandLesion))
  
  cc_df_hm <- cc_hm[cc_hm$SpeciesandLesion == filter[1],]
  
  colnames(cc_df_hm)[3] <- paste0(unique(cc_df_hm$SpeciesandLesion), "_GeneRatio")
  
  cc_df_hm <- cc_df_hm[,2:3]
  
  i <- 3
  
  for (i in 2:length(filter)) {
    
    loop_cc_hm <- cc_hm[cc_hm$SpeciesandLesion == filter[i],]
    
    colnames(loop_cc_hm)[3] <- paste0(unique(loop_cc_hm$SpeciesandLesion), "_GeneRatio")
    
    loop_cc_hm <- loop_cc_hm[,2:3]
    
    cc_df_hm <- full_join(x = cc_df_hm,y = loop_cc_hm, by = "description")
    
    
  }
  
  cc_df_hm <- as.data.frame(cc_df_hm)
  
  rownames(cc_df_hm) <- cc_df_hm[,1]
  
  cc_df_hm <- cc_df_hm[,c(2:ncol(cc_df_hm))]
  
  cc_df_hm <- as.matrix(cc_df_hm)
  
  colours <- c(viridis(n = 100, option = "G",direction = -1))
  
  Experiment <- c(rep("8oxoG", length(grep(pattern = "8oxoG",x = colnames(cc_df_hm)))),
                  rep("abasic", length(grep(pattern = "abasic",x = colnames(cc_df_hm)))),
                  rep("RNAbase", length(grep(pattern = "RNAbase",x = colnames(cc_df_hm)))))
  
  df_col <- data.frame(Experiment)
  
  rownames(df_col) <- colnames(cc_df_hm)
  
  Experiment <- c("#000000", "#ffffff","#939598")
  
  names(Experiment) <- c("8oxoG", "abasic", "RNAbase")
  
  anno_colors <- list(Experiment = Experiment)
  
  filename <- paste0("03_OutputFiles/AdditionalANDProteinID/",
                     gsub(pattern = "-",replacement = "",x = Sys.Date()),
                     "_02_AllnoFacet_top",
                     top[j], "_GOCC_HeatMap.pdf")
  
  pheatmap(cc_df_hm,
           col=colours,
           cluster_cols = F,
           cluster_rows = F,
           show_rownames = T,
           show_colnames = T,
           labels_col = str_to_title(gsub(pattern = " - .*",replacement = "",x = colnames(cc_df_hm))),
           annotation_col = df_col,
           # annotation_row = df_row,
           annotation_colors = anno_colors,
           fontsize = 16,border_color = "#b4b4b4",
           height = 15,width = 15,filename = filename
  )
  
  kegg_dp <- enrichment_df[enrichment_df$category == "KEGG",]
  
  kegg_dp <- kegg_dp %>% group_by(Species, Lesion) %>% slice_max(fdr,n = top[j])
  
  kegg_dp$Ratio <- kegg_dp$number_of_genes/kegg_dp$number_of_genes_in_background 
  
  kegg_dp$Lesion <- factor(x = kegg_dp$Lesion, levels = c("8oxoG", "abasic", "RNAbase"))
  
  kegg_dp$Species <- factor(x = kegg_dp$Species, levels = c("E. coli",
                                                        "B. subtilis",
                                                        "H. salinarum",
                                                        "T. brucei",
                                                        "T. thermophila",
                                                        "S. cerevisiae",
                                                        "S. pombe",
                                                        "C. elegans",
                                                        "H. sapiens (HeLa)",
                                                        "H. sapiens (HEK293)",
                                                        "Z. mays",
                                                        "A. thaliana"))
  
  plot <- ggplot(kegg_dp, aes(fill = fdr,y = description, size = Ratio, x= Species)) +
    geom_point()+
    facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
    # facet_wrap(~Species, scales = "free_y", ncol = 1)+
    labs(fill = "FDR",size = "Gene Ratio")+
    scale_fill_viridis_c(option = "G")+
    scale_size_binned()+
    theme_bw()+
    # ggtitle(species[i])+
    theme(axis.text.x=element_text(size=16,hjust = 1, angle = 45),
          axis.text.y=element_text(size=16),
          axis.title=element_text(size=18,face="bold"),
          plot.title = element_text(size=20, face="bold"))+
    theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
          legend.title = element_text(size = 16,face="bold"))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.position = "top")+
    theme(strip.text = element_text(size=16))+
    theme(axis.text.x = element_text(face = "italic"))+
    xlab("Species")+
    ylab("KEGG term")
  
  kegg_dp$SpeciesandLesion <- paste0(kegg_dp$Species, " - ", kegg_dp$Lesion)
  
  ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                           gsub(pattern = "-",replacement = "",x = Sys.Date()),
                           "_02_All_top",
                           top[j], "_KEGG_Dotplot.pdf"),plot = plot,width = 15,height = 15)
  
  plot <- ggplot(kegg_dp, aes(fill = fdr,y = description, size = Ratio, x= SpeciesandLesion)) +
    geom_point()+
    # facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
    # facet_wrap(~Species, scales = "free_y", ncol = 1)+
    labs(fill = "FDR",size = "Gene Ratio")+
    scale_fill_viridis_c(option = "G")+
    scale_size_binned()+
    theme_bw()+
    # ggtitle(species[i])+
    theme(axis.text.x=element_text(size=16,hjust = 1, angle = 45),
          axis.text.y=element_text(size=16),
          axis.title=element_text(size=18,face="bold"),
          plot.title = element_text(size=20, face="bold"))+
    theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
          legend.title = element_text(size = 16,face="bold"))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.position = "top")+
    theme(strip.text = element_text(size=16))+
    theme(axis.text.x = element_text(face = "italic"))+
    xlab("Species and Lesion")+
    ylab("KEGG term")
  
  ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                           gsub(pattern = "-",replacement = "",x = Sys.Date()),
                           "_02_AllnoFacet_top",
                           top[j], "_KEGG_Dotplot.pdf"),plot = plot,width = 15,height = 15)
  
  kegg_hm <- enrichment_df[enrichment_df$category == "KEGG",]
  
  kegg_hm <- kegg_hm %>% group_by(Species, Lesion) %>% slice_max(fdr,n = top[j])
  
  kegg_hm$Ratio <- kegg_hm$number_of_genes/kegg_hm$number_of_genes_in_background 
  
  kegg_hm$Lesion <- factor(x = kegg_hm$Lesion, levels = c("8oxoG", "abasic", "RNAbase"))
  
  kegg_hm$Species <- factor(x = kegg_hm$Species, levels = c("E. coli",
                                                        "B. subtilis",
                                                        "H. salinarum",
                                                        "T. brucei",
                                                        "T. thermophila",
                                                        "S. cerevisiae",
                                                        "S. pombe",
                                                        "C. elegans",
                                                        "H. sapiens (HeLa)",
                                                        "H. sapiens (HEK293)",
                                                        "Z. mays",
                                                        "A. thaliana"))
  
  kegg_hm$SpeciesandLesion <- paste0(kegg_dp$Species, " - ", kegg_dp$Lesion)
  
  kegg_hm$SpeciesandLesion <- factor(x = kegg_hm$SpeciesandLesion, levels = c("E. coli - 8oxoG",
                                                                          "B. subtilis - 8oxoG",
                                                                          "H. salinarum - 8oxoG",
                                                                          "T. brucei - 8oxoG",
                                                                          "T. thermophila - 8oxoG",
                                                                          "S. cerevisiae - 8oxoG",
                                                                          "S. pombe - 8oxoG",
                                                                          "C. elegans - 8oxoG",
                                                                          "H. sapiens (HeLa) - 8oxoG",
                                                                          "H. sapiens (HEK293) - 8oxoG",
                                                                          "Z. mays - 8oxoG",
                                                                          "A. thaliana - 8oxoG",
                                                                          "E. coli - abasic",
                                                                          "B. subtilis - abasic",
                                                                          "H. salinarum - abasic",
                                                                          "T. brucei - abasic",
                                                                          "T. thermophila - abasic",
                                                                          "S. cerevisiae - abasic",
                                                                          "S. pombe - abasic",
                                                                          "C. elegans - abasic",
                                                                          "H. sapiens (HeLa) - abasic",
                                                                          "H. sapiens (HEK293) - abasic",
                                                                          "Z. mays - abasic",
                                                                          "A. thaliana - abasic",
                                                                          "E. coli - RNAbase",
                                                                          "B. subtilis - RNAbase",
                                                                          "H. salinarum - RNAbase",
                                                                          "T. brucei - RNAbase",
                                                                          "T. thermophila - RNAbase",
                                                                          "S. cerevisiae - RNAbase",
                                                                          "S. pombe - RNAbase",
                                                                          "C. elegans - RNAbase",
                                                                          "H. sapiens (HeLa) - RNAbase",
                                                                          "H. sapiens (HEK293) - RNAbase",
                                                                          "Z. mays - RNAbase",
                                                                          "A. thaliana - RNAbase"))
  
  kegg_hm <- kegg_hm %>% arrange(SpeciesandLesion)
  
  kegg_hm <- kegg_hm[,c(14,12,13)]
  
  filter <- unique(as.character(kegg_hm$SpeciesandLesion))
  
  kegg_df_hm <- kegg_hm[kegg_hm$SpeciesandLesion == filter[1],]
  
  colnames(kegg_df_hm)[3] <- paste0(unique(kegg_df_hm$SpeciesandLesion), "_GeneRatio")
  
  kegg_df_hm <- kegg_df_hm[,2:3]
  
  i <- 3
  
  for (i in 2:length(filter)) {
    
    loop_kegg_hm <- kegg_hm[kegg_hm$SpeciesandLesion == filter[i],]
    
    colnames(loop_kegg_hm)[3] <- paste0(unique(loop_kegg_hm$SpeciesandLesion), "_GeneRatio")
    
    loop_kegg_hm <- loop_kegg_hm[,2:3]
    
    kegg_df_hm <- full_join(x = kegg_df_hm,y = loop_kegg_hm, by = "description")
    
    
  }
  
  kegg_df_hm <- as.data.frame(kegg_df_hm)
  
  rownames(kegg_df_hm) <- kegg_df_hm[,1]
  
  kegg_df_hm <- kegg_df_hm[,c(2:ncol(kegg_df_hm))]
  
  kegg_df_hm <- as.matrix(kegg_df_hm)
  
  colours <- c(viridis(n = 100, option = "G",direction = -1))
  
  Experiment <- c(rep("8oxoG", length(grep(pattern = "8oxoG",x = colnames(kegg_df_hm)))),
                  rep("abasic", length(grep(pattern = "abasic",x = colnames(kegg_df_hm)))),
                  rep("RNAbase", length(grep(pattern = "RNAbase",x = colnames(kegg_df_hm)))))
  
  df_col <- data.frame(Experiment)
  
  rownames(df_col) <- colnames(kegg_df_hm)
  
  Experiment <- c("#000000", "#ffffff","#939598")
  
  names(Experiment) <- c("8oxoG", "abasic", "RNAbase")
  
  anno_colors <- list(Experiment = Experiment)
  
  filename <- paste0("03_OutputFiles/AdditionalANDProteinID/",
                     gsub(pattern = "-",replacement = "",x = Sys.Date()),
                     "_02_AllnoFacet_top",
                     top[j], "_KEGG_HeatMap.pdf")
  
  pheatmap(kegg_df_hm,
           col=colours,
           cluster_cols = F,
           cluster_rows = F,
           show_rownames = T,
           show_colnames = T,
           labels_col = str_to_title(gsub(pattern = " - .*",replacement = "",x = colnames(kegg_df_hm))),
           annotation_col = df_col,
           # annotation_row = df_row,
           annotation_colors = anno_colors,
           fontsize = 16,border_color = "#b4b4b4",
           height = 15,width = 15,filename = filename
  )

  pfam_dp <- enrichment_df[enrichment_df$category == "Pfam",]
  
  pfam_dp <- pfam_dp %>% group_by(Species, Lesion) %>% slice_max(fdr,n = top[j])
  
  pfam_dp$Ratio <- pfam_dp$number_of_genes/pfam_dp$number_of_genes_in_background 
  
  pfam_dp$Lesion <- factor(x = pfam_dp$Lesion, levels = c("8oxoG", "abasic", "RNAbase"))
  
  pfam_dp$Species <- factor(x = pfam_dp$Species, levels = c("E. coli",
                                                            "B. subtilis",
                                                            "H. salinarum",
                                                            "T. brucei",
                                                            "T. thermophila",
                                                            "S. cerevisiae",
                                                            "S. pombe",
                                                            "C. elegans",
                                                            "H. sapiens (HeLa)",
                                                            "H. sapiens (HEK293)",
                                                            "Z. mays",
                                                            "A. thaliana"))
  
  plot <- ggplot(pfam_dp, aes(fill = fdr,y = description, size = Ratio, x= Species)) +
    geom_point()+
    facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
    # facet_wrap(~Species, scales = "free_y", ncol = 1)+
    labs(fill = "FDR",size = "Gene Ratio")+
    scale_fill_viridis_c(option = "G")+
    scale_size_binned()+
    theme_bw()+
    # ggtitle(species[i])+
    theme(axis.text.x=element_text(size=16,hjust = 1, angle = 45),
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=18,face="bold"),
          plot.title = element_text(size=20, face="bold"))+
    theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
          legend.title = element_text(size = 16,face="bold"))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.position = "top")+
    theme(strip.text = element_text(size=16))+
    theme(axis.text.x = element_text(face = "italic"))+
    xlab("Species")+
    ylab("Pfam term")
  
  ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                           gsub(pattern = "-",replacement = "",x = Sys.Date()),
                           "_02_All_top",
                           top[j], "_Pfam_Dotplot.pdf"),plot = plot,width = 15,height = 15)
  
  pfam_dp$SpeciesandLesion <- paste0(pfam_dp$Species, " - ", pfam_dp$Lesion)
  
  plot <- ggplot(pfam_dp, aes(fill = fdr,y = description, size = Ratio, x= SpeciesandLesion)) +
    geom_point()+
    # facet_wrap(~Lesion,scales = "free_y", ncol = 1)+
    # facet_wrap(~Species, scales = "free_y", ncol = 1)+
    labs(fill = "FDR",size = "Gene Ratio")+
    scale_fill_viridis_c(option = "G")+
    scale_size_binned()+
    theme_bw()+
    # ggtitle(species[i])+
    theme(axis.text.x=element_text(size=16,hjust = 1, angle = 45),
          axis.text.y=element_text(size=16),
          axis.title=element_text(size=18,face="bold"),
          plot.title = element_text(size=20, face="bold"))+
    theme(legend.text = element_text(size = 14,hjust = 1, angle = 45),
          legend.title = element_text(size = 16,face="bold"))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.position = "top")+
    theme(strip.text = element_text(size=16))+
    theme(axis.text.x = element_text(face = "italic"))+
    xlab("Species and Lesion")+
    ylab("Pfam term")
  
  ggsave(filename = paste0("03_OutputFiles/AdditionalANDProteinID/",
                           gsub(pattern = "-",replacement = "",x = Sys.Date()),
                           "_02_AllnoFacet_top",
                           top[j], "_Pfam_Dotplot.pdf"),plot = plot,width = 15,height = 15)
  
  pfam_hm <- enrichment_df[enrichment_df$category == "Pfam",]
  
  pfam_hm <- pfam_hm %>% group_by(Species, Lesion) %>% slice_max(fdr,n = top[j])
  
  pfam_hm$Ratio <- pfam_hm$number_of_genes/pfam_hm$number_of_genes_in_background 
  
  pfam_hm$Lesion <- factor(x = pfam_hm$Lesion, levels = c("8oxoG", "abasic", "RNAbase"))
  
  pfam_hm$Species <- factor(x = pfam_hm$Species, levels = c("E. coli",
                                                        "B. subtilis",
                                                        "H. salinarum",
                                                        "T. brucei",
                                                        "T. thermophila",
                                                        "S. cerevisiae",
                                                        "S. pombe",
                                                        "C. elegans",
                                                        "H. sapiens (HeLa)",
                                                        "H. sapiens (HEK293)",
                                                        "Z. mays",
                                                        "A. thaliana"))
  
  pfam_hm$SpeciesandLesion <- paste0(pfam_dp$Species, " - ", pfam_dp$Lesion)
  
  pfam_hm$SpeciesandLesion <- factor(x = pfam_hm$SpeciesandLesion, levels = c("E. coli - 8oxoG",
                                                                          "B. subtilis - 8oxoG",
                                                                          "H. salinarum - 8oxoG",
                                                                          "T. brucei - 8oxoG",
                                                                          "T. thermophila - 8oxoG",
                                                                          "S. cerevisiae - 8oxoG",
                                                                          "S. pombe - 8oxoG",
                                                                          "C. elegans - 8oxoG",
                                                                          "H. sapiens (HeLa) - 8oxoG",
                                                                          "H. sapiens (HEK293) - 8oxoG",
                                                                          "Z. mays - 8oxoG",
                                                                          "A. thaliana - 8oxoG",
                                                                          "E. coli - abasic",
                                                                          "B. subtilis - abasic",
                                                                          "H. salinarum - abasic",
                                                                          "T. brucei - abasic",
                                                                          "T. thermophila - abasic",
                                                                          "S. cerevisiae - abasic",
                                                                          "S. pombe - abasic",
                                                                          "C. elegans - abasic",
                                                                          "H. sapiens (HeLa) - abasic",
                                                                          "H. sapiens (HEK293) - abasic",
                                                                          "Z. mays - abasic",
                                                                          "A. thaliana - abasic",
                                                                          "E. coli - RNAbase",
                                                                          "B. subtilis - RNAbase",
                                                                          "H. salinarum - RNAbase",
                                                                          "T. brucei - RNAbase",
                                                                          "T. thermophila - RNAbase",
                                                                          "S. cerevisiae - RNAbase",
                                                                          "S. pombe - RNAbase",
                                                                          "C. elegans - RNAbase",
                                                                          "H. sapiens (HeLa) - RNAbase",
                                                                          "H. sapiens (HEK293) - RNAbase",
                                                                          "Z. mays - RNAbase",
                                                                          "A. thaliana - RNAbase"))
  
  pfam_hm <- pfam_hm %>% arrange(SpeciesandLesion)
  
  pfam_hm <- pfam_hm[,c(14,12,13)]
  
  filter <- unique(as.character(pfam_hm$SpeciesandLesion))
  
  pfam_df_hm <- pfam_hm[pfam_hm$SpeciesandLesion == filter[1],]
  
  colnames(pfam_df_hm)[3] <- paste0(unique(pfam_df_hm$SpeciesandLesion), "_GeneRatio")
  
  pfam_df_hm <- pfam_df_hm[,2:3]
  
  i <- 3
  
  for (i in 2:length(filter)) {
    
    loop_pfam_hm <- pfam_hm[pfam_hm$SpeciesandLesion == filter[i],]
    
    colnames(loop_pfam_hm)[3] <- paste0(unique(loop_pfam_hm$SpeciesandLesion), "_GeneRatio")
    
    loop_pfam_hm <- loop_pfam_hm[,2:3]
    
    pfam_df_hm <- full_join(x = pfam_df_hm,y = loop_pfam_hm, by = "description")
    
    
  }
  
  pfam_df_hm <- as.data.frame(pfam_df_hm)
  
  rownames(pfam_df_hm) <- pfam_df_hm[,1]
  
  pfam_df_hm <- pfam_df_hm[,c(2:ncol(pfam_df_hm))]
  
  pfam_df_hm <- as.matrix(pfam_df_hm)
  
  colours <- c(viridis(n = 100, option = "G",direction = -1))
  
  Experiment <- c(rep("8oxoG", length(grep(pattern = "8oxoG",x = colnames(pfam_df_hm)))),
                  rep("abasic", length(grep(pattern = "abasic",x = colnames(pfam_df_hm)))),
                  rep("RNAbase", length(grep(pattern = "RNAbase",x = colnames(pfam_df_hm)))))
  
  df_col <- data.frame(Experiment)
  
  rownames(df_col) <- colnames(pfam_df_hm)
  
  Experiment <- c("#000000", "#ffffff","#939598")
  
  names(Experiment) <- c("8oxoG", "abasic", "RNAbase")
  
  anno_colors <- list(Experiment = Experiment)
  
  filename <- paste0("03_OutputFiles/AdditionalANDProteinID/",
                     gsub(pattern = "-",replacement = "",x = Sys.Date()),
                     "_02_AllnoFacet_top",
                     top[j], "_Pfam_HeatMap.pdf")
  
  pheatmap(pfam_df_hm,
           col=colours,
           cluster_cols = F,
           cluster_rows = F,
           show_rownames = T,
           show_colnames = T,
           labels_col = str_to_title(gsub(pattern = " - .*",replacement = "",x = colnames(pfam_df_hm))),
           annotation_col = df_col,
           # annotation_row = df_row,
           annotation_colors = anno_colors,
           fontsize = 16,border_color = "#b4b4b4",
           height = 15,width = 15,filename = filename
  )
  
}
