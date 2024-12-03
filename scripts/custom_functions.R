#exclude unclassified and homo sapiens
get_top <- function(df, number_of_hits=5){
  
  top_df <- df %>% dplyr::select(-rel_abun) %>%
    dplyr::filter(!k2_taxon %in% c("Homo sapiens", "unclassified"))
  
  #Recalculate relative abundance while excluding Human and unclassified reads
  top_df$rel_abun <- (top_df$paired_reads / sum(top_df$paired_reads))*100
  
  top_df <- top_df  %>%
    slice_max(rel_abun, n=number_of_hits) %>% dplyr::select(rel_abun, k2_taxon, LIBID)
  
  other_rel_abun <- 100 - sum(top_df$rel_abun)
  library <- top_df$LIBID %>% unique()
  
  other_df <- data.frame(k2_taxon="others", 
                         rel_abun=other_rel_abun,
                         LIBID=library)
  output <- rbind(top_df, other_df)
  
  taxa_to_add <- top_df$k2_taxon
  
  output$k2_taxon <- factor(output$k2_taxon)
  
  output$k2_taxon <- relevel(output$k2_taxon, "others")
  
  return(output)
}


get_simple_profile <- function(input_list, input_metadata, n) {
  
  profile <- lapply(input_list, function(x){get_top(x,n)})
  
  profile <- do.call("rbind", profile)
  
  taxa_to_add <- profile %>% pull(k2_taxon) %>% unique() %>% sort()
  taxa_to_add <- taxa_to_add[!taxa_to_add %in% c("others" )] %>% as.character()
  
  #Now this ENTIRE VECTOR of taxa_to_add will be used to re-separate the species of EACH library into "others" vs the species to be plotted on a bar chart
  
  final_profile <- lapply(input_list, function(df){
    
    taxa_df <- df  %>% dplyr::filter(k2_taxon %in% taxa_to_add) %>%
      dplyr::select(paired_reads, k2_taxon, LIBID)
    
    others_df <- df %>% 
      dplyr::filter(!k2_taxon %in% c("unclassified", "Homo sapiens",
                                     taxa_to_add)) %>%
      dplyr::select(paired_reads, k2_taxon, LIBID)
    
    #The denominator for the new relative abundance
    total_reads <- sum(taxa_df$paired_reads) + sum(others_df$paired_reads)
    
    taxa_df$rel_abun <- (taxa_df$paired_reads/total_reads)*100
    
    other_rel_abun <- (sum(others_df$paired_reads) / total_reads)*100
    
    id <- taxa_df %>% pull(LIBID) %>% unique()
    
    other_df <- data.frame(k2_taxon="others", 
                           rel_abun=other_rel_abun,
                           LIBID=id)
    
    df_out <- rbind(taxa_df %>% dplyr::select(-paired_reads), other_df)
    
    
    return(df_out)
    
  })
  
  final_profile <- do.call("rbind", final_profile)
  
  output <- merge (final_profile, input_metadata, by = "LIBID", all.x=TRUE)
  
  #re-arrange factor levels so that "others" always come first :)
  #Unclassified read fraction should have been removed prior to this.
  output$k2_taxon <- factor(output$k2_taxon)
  
  output$k2_taxon <- relevel(output$k2_taxon, "others")
  
  return(output)
  
} 


###function to get top n relative abundance of any species. Anything else outside the top n will be under "others"
#unclassified and homo_sapiens were already filtered out from the inputs
#sample_identifier can be a string like "CRAM_ID", "LIBRARY_ID", "LIBID" etc
##Here,we only pick the top species that have passed the k2 minimizer check for false positives.
get_top_filtered <- function(df, number_of_hits=5, sample_identifier){
  
  top_df <- df %>% dplyr::filter(taxa_pass_filter==TRUE) %>%
    slice_max(rel_abun, n=number_of_hits) %>% dplyr::select(rel_abun, k2_taxon, all_of(sample_identifier))
  
  other_rel_abun <- 100 - sum(top_df$rel_abun)
  identifier_vec <- top_df %>% pull(get(sample_identifier)) %>% unique()
  
  other_df <- data.frame(k2_taxon="others", 
                         rel_abun=other_rel_abun,
                         PLACEHOLDER=identifier_vec)
  
  colnames(other_df)[3] <- sample_identifier
  
  output <- rbind(top_df, other_df)
  
  taxa_to_add <- top_df$k2_taxon
  
  output$k2_taxon <- factor(output$k2_taxon)
  
  output$k2_taxon <- relevel(output$k2_taxon, "others")
  
  return(output)
}

#subset_key, identifier and merging_key can be a string like "Library_ID"
get_simple_profile_filtered <- function(input_list, input_metadata, subset_key, n, identifier, merging_key) {
  
  subset_vector <- input_metadata %>% pull(get(subset_key))
  
  input_list_subset <- input_list[subset_vector]
  
  profile <- lapply(input_list_subset, function(x){get_top_filtered(x,n, sample_identifier = identifier)})
  
  profile <- do.call("rbind", profile)
  
  taxa_to_add <- profile %>% pull(k2_taxon) %>% unique() %>% sort()
  taxa_to_add <- taxa_to_add[!taxa_to_add %in% c("others" )] %>% as.character()
  
  #Now this ENTIRE VECTOR of taxa_to_add will be used to re-separate the species of EACH library into "others" vs the species to   be plotted on a bar chart
  
  final_profile <- lapply(input_list_subset, function(df){
    
    taxa_df <- df %>% dplyr::filter(taxa_pass_filter==TRUE) %>% dplyr::filter(k2_taxon %in% taxa_to_add) %>%
      dplyr::select(rel_abun, k2_taxon, all_of(identifier))
    
    #This works because we are not removing any features from "taxa_df" = not changing the denominator for rel abundance  
    other_rel_abun <- 100 - sum(taxa_df$rel_abun)
    id <- taxa_df %>% pull(get(identifier)) %>% unique()
    
    other_df <- data.frame(k2_taxon="others", 
                           rel_abun=other_rel_abun,
                           PLACEHOLDER=id)
    
    colnames(other_df)[3] <- identifier
    
    df_out <- rbind(taxa_df, other_df)
    
    
    return(df_out)
    
  })
  
  final_profile <- do.call("rbind", final_profile)
  
  output <- merge (final_profile, input_metadata, by = merging_key, all.x=TRUE)
  
  #re-arrange factor levels so that "others" always come first :)
  #Unclassified read fraction should have been removed prior to this.
  output$k2_taxon <- factor(output$k2_taxon)
  
  output$k2_taxon <- relevel(output$k2_taxon, "others")
  
  return(output)
  
} 


#colorblind friendly palettes

c23 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2",
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)



c24 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2",
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)



rel_abun_plot_fn <- function(df, x_var, title){
  
  x_var_to_plot <- enquo(x_var)
  
  df %>%
    ggplot(aes(fill=k2_taxon, y=rel_abun, x=!!x_var_to_plot)) +   #change x to SAMPLE_ID for larger datasets
    geom_col()+
    theme_classic() +
    scale_fill_manual(values = c24 ) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(legend.text=element_text(size=8), axis.title.x=element_blank()) +
    ggtitle(title) + 
    ylab("relative abundance") +guides(fill=guide_legend(ncol=2))}



rel_abun_filt_plot_fn <- function(df, x_var, title){
  
  x_var_to_plot <- enquo(x_var)
  
  df %>%
    ggplot(aes(fill=k2_taxon, y=rel_abun, x=!!x_var_to_plot)) + 
    geom_col()+
    theme_classic() +
    scale_fill_manual(values = c24 ) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.text=element_text(size=8), axis.title.x=element_blank()) +
    ggtitle(title) + 
    ylab("relative abundance") +guides(fill=guide_legend(ncol=2))}


rel_abun_plot_fn2 <- function(df, x_var, title){
  
  x_var_to_plot <- enquo(x_var)
  
  df %>%
    ggplot(aes(fill=k2_taxon, y=rel_abun, x=!!x_var_to_plot)) +
    geom_col()+
    theme_classic() +
    scale_fill_manual(values = Glasbey26 ) +  #Use Glasbey26 to expand the color palette if necessary
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(legend.text=element_text(size=8), axis.title.x=element_blank()) +
    ggtitle(title) + 
    ylab("relative abundance") }


#For each site for each species, compute the geometric mean and median transcriptional activity.
geom_mean <- function(x){
  exp(mean(log(x)))
}



compute_site_prevalence <- function(k2_df, LIB_vec, site){
  
  site_df <- k2_df %>% dplyr::filter(LIBID %in% LIB_vec)
  
  #number of libraries for a skin site out of the 102
  n_site <- length(LIB_vec)
  
  #Organism prevalence per site
  freq_measured_df <- plyr::count(site_df,vars = "k2_taxon")
  freq_measured_df$site_size <- n_site
  
  freq_measured_df$percent_detected_at_site <- (freq_measured_df$freq/freq_measured_df$site_size)*100
  
  freq_measured_df$region <- site
  
  return(freq_measured_df)
  
}

#mol_type can be "DNA" or "RNA"
get_species_prevalence_all_sites <- function(k2_df, mol_type){
  
  input <- k2_df
  
  if (mol_type=="RNA"){
    Sc_prevalence <- compute_site_prevalence(input, mtx_Sc_ids ,site="Sc")
    Ch_prevalence <- compute_site_prevalence(input, mtx_Ch_ids ,site="Ch")
    Ac_prevalence <- compute_site_prevalence(input, mtx_Ac_ids ,site="Ac")
    Vf_prevalence <- compute_site_prevalence(input, mtx_Vf_ids ,site="Vf")
    Tw_prevalence <- compute_site_prevalence(input, mtx_Tw_ids ,site="Tw")
    
  } else if (mol_type == "DNA"){
    Sc_prevalence <- compute_site_prevalence(input, mgx_Sc_ids ,site="Sc")
    Ch_prevalence <- compute_site_prevalence(input, mgx_Ch_ids ,site="Ch")
    Ac_prevalence <- compute_site_prevalence(input, mgx_Ac_ids ,site="Ac")
    Vf_prevalence <- compute_site_prevalence(input, mgx_Vf_ids ,site="Vf")
    Tw_prevalence <- compute_site_prevalence(input, mgx_Tw_ids ,site="Tw")
    
  }
  
  
  output <- do.call("rbind", list(Sc_prevalence,
                                  Ch_prevalence,
                                  Ac_prevalence,
                                  Vf_prevalence,
                                  Tw_prevalence))
  
  return(output)
  
}

define_core_variable <- function(df, min_detect=3){
  
  #Showing the core and variable taxa which are detected in at least THREE samples
  df_filt <- df %>% dplyr::filter(freq>=min_detect)
  
  df_filt$category_within_site <- ifelse(df_filt$percent_detected_at_site > 75, "core",
                                         ifelse(df_filt$percent_detected_at_site > 50, "common",
                                                "variable")) 
  
  return(df_filt)
}

#species contribution to a specific pathway on top, mtx relative abundance at bottom
facet_pathway_vs_community_abundance_plot_fn <- function(mtx_list=rna_k2_minimizer_renorm,
                                                         pathway_df=mtx_pathabun_stratified_TSS,
                                                         reaction_module, #e.g. "M00861"
                                                         metadata = mtx_mgx_stats_chosen,
                                                         chosen_species,
                                                         chosen_region,
                                                         manual_sample_order=FALSE,
                                                         subj_region_order #a vector
){
  
  options(dplyr.summarise.inform = FALSE)
  
  mtx_libs_to_choose <- metadata %>% dplyr::filter(region==chosen_region) %>% pull(mtx_LIBID)
  
  
  mtx_df <- lapply(mtx_list[mtx_libs_to_choose], function(df){
    
    df <- merge(df %>% dplyr::rename(LIBID=mtx_LIBID), 
                metadata %>% dplyr::rename(LIBID=mtx_LIBID) %>% dplyr::select(LIBID, subj_region), by="LIBID")
    
    df$species_fmt<- ifelse(df$k2_taxon %in% chosen_species, df$k2_taxon, "others")
    
    df_summarized <- df %>% group_by(LIBID, subj_region, species_fmt) %>% summarise(rel_abun_sum=sum(rel_abun))
    
    df_summarized$type <- "community"
    
    return(df_summarized)
    
  }) %>% do.call("rbind",.) %>% dplyr::select(species_fmt, subj_region, rel_abun_sum, type)
  
  
  pathway_df_subset <- mtx_pathabun_stratified_TSS %>% dplyr::filter(module == reaction_module & 
                                                                       region == chosen_region) %>% 
    dplyr::rename(rel_abun_sum=rel_abun_CoPM) %>% dplyr::select(species_fmt, subj_region, rel_abun_sum)
  
  
  pathway_df_subset$type <- reaction_module
  
  #Genus level and above classifications were abstracted to "unclassified", which means "unclassified at species level"
  pathway_df_subset$species_fmt <- gsub(pattern="unclassified", replacement="unclassified at species level", x=pathway_df_subset$species_fmt)
  
  pathway_df_subset$species_fmt <- ifelse(pathway_df_subset$species_fmt %in% c(chosen_species,"unclassified at species level"), 
                                          pathway_df_subset$species_fmt, "others")
  
  
  combined_df <- rbind(pathway_df_subset, mtx_df)  
  
  #Re-order factor levels to put "others" first.
  #renormalized RNA K2 abundances do not have "unclassified" any more
  taxa_to_add <- chosen_species[!chosen_species %in% c("unclassified" )]
  #combined_df$species_relabelled <- factor(combined_df$species_relabelled)
  
  combined_df<- combined_df %>% mutate(species_fmt= factor(species_fmt, 
                                                           levels = c(taxa_to_add,"others", "unclassified at species level")))
  
  if (manual_sample_order==TRUE){
    combined_df$subj_region <- factor(combined_df$subj_region, levels=subj_region_order)
  }
  
  combined_df$type <- factor(combined_df$type, levels=c( reaction_module,"community"))
  
  
  ggplot(combined_df, 
         aes(x=subj_region, y=rel_abun_sum, fill=species_fmt)) + geom_col() + #coord_flip() #+ 
    facet_wrap(~type, ncol=1, nrow=2)+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x=element_blank()) + 
    ggtitle(paste0("Metatranscriptomes, ",chosen_region)) + scale_fill_manual(values = c24)
  
  
}

##This function creates a custom term2gene dataframe for specific Malassezia species, based on Eggnog annotated malassezia transcriptomes.

get_Malassezia_term2gene_df <- function(tx2gene=Malassezia_tx2gene, species, annotation_used="KEGG"){
  
  if (species== "Malassezia_restricta"){
    
    loci<-"DNF"
    
  } else if (species=="Malassezia_globosa"){
    
    loci<-"MGL"
    
  } else if (species=="Malassezia_sympodialis"){
    
    loci<-"MSYG"
    
  } else if (species=="Malassezia_furfur"){
    
    loci <-"CBS"
    
  } else if (species=="Malassezia_arunalokei"){
    
    loci <-"MARU"
  }
  
  #Get the species specific tx2gene df
  
  tx2gene_specific <- tx2gene %>% dplyr::filter(str_detect(tx2gene$gene, pattern = loci))
  #Add  KEGG pathway terms
  tx2gene_specific <- merge(tx2gene_specific, Malassezia_OG_annot, by = c("query",
                                                                          "gene"), all.x=TRUE)
  
  tx2gene_specific[is.na( tx2gene_specific)] <- "-"
  
  if (annotation_used=="KEGG"){
    Malassezia_anno_KO <-  tx2gene_specific %>% dplyr::select(c("gene", "KEGG_Pathway")) %>% 
      dplyr::filter(KEGG_Pathway != "-")
    
    Malassezia_anno_KO <- Malassezia_anno_KO %>% tidyr::separate_rows(KEGG_Pathway, sep = ",")
    
    Malassezia_anno_map <- Malassezia_anno_KO %>% dplyr::filter(str_detect(string=Malassezia_anno_KO$KEGG_Pathway,
                                                                           pattern="map"))
    
    #Create the species specific term2gene dataframe for KEGG pathways
    
    term2gene <- data.frame(from=Malassezia_anno_map$KEGG_Pathway, to=Malassezia_anno_map$gene)
    
  } else if (annotation_used == "GO"){
    
    Malassezia_anno_GO <- tx2gene_specific %>% dplyr::select(c("gene", "GOs")) %>% 
      dplyr::filter(GOs != "-")
    
    
    Malassezia_anno_GO <- Malassezia_anno_GO %>% tidyr::separate_rows(GOs, sep = ",") %>% unique()
    
    term2gene <- Malassezia_anno_GO %>% dplyr::rename(from="GOs", to="gene") %>% dplyr::select(from, to)
  }
  
  return(term2gene)
  
}


calculate_species_TPM <- function(df, species){
  
  df_filt <- df %>% dplyr::filter(str_detect(df$pangene, pattern = species))
  
  df_filt$RPK <- (df_filt$unpaired_read_count/df_filt$length) *1000
  
  #does not count unmapped reads
  per_mill_scaling_factor <- sum(df_filt$RPK)/1E6
  
  df_filt$TPM <- df_filt$RPK/per_mill_scaling_factor
  
  df_filt$propan_clusterID <- str_match(string=df_filt$pangene, pattern="cluster_[0-9]*")[,1]
  
  #merge clusters together to fit the propan definitions
  df_filt_merge <- df_filt %>% group_by(propan_clusterID) %>% 
    summarise(unpaired_read_count_sum=sum(unpaired_read_count),
              RPK_sum=sum(RPK),
              TPM_sum=sum(TPM)) %>% ungroup()
  
  df_out <- df_filt_merge %>% dplyr::select(propan_clusterID, unpaired_read_count_sum,RPK_sum, TPM_sum) %>% unique()
  df_out$species <- species
  
  library_id <- df_filt %>% pull(LIBID) %>% unique()
  df_out$LIBID <- library_id
  
  return(df_out)
  
  
}

multi_facet_pathway_vs_community_abundance_plot_fn <- function(mtx_list=rna_k2_minimizer_renorm,
                                                               pathway_df=mtx_pathabun_stratified_TSS,
                                                               reaction_module, #e.g. "M00861"
                                                               metadata = mtx_mgx_stats_chosen,
                                                               chosen_species,
                                                               chosen_region, #can take multiple regions
                                                               manual_sample_order=FALSE,
                                                               subj_region_order #a vector
){
  
  options(dplyr.summarise.inform = FALSE)
  
  mtx_libs_to_choose <- metadata %>% dplyr::filter(region %in% chosen_region) %>% pull(mtx_LIBID)
  
  
  mtx_df <- lapply(mtx_list[mtx_libs_to_choose], function(df){
    
    df <- merge(df %>% dplyr::rename(LIBID=mtx_LIBID), 
                metadata %>% dplyr::rename(LIBID=mtx_LIBID) %>% dplyr::select(LIBID, region, subj_region), by="LIBID")
    
    df$species_fmt<- ifelse(df$k2_taxon %in% chosen_species, df$k2_taxon, "others")
    
    df_summarized <- df %>% group_by(LIBID, region, subj_region, species_fmt) %>% summarise(rel_abun_sum=sum(rel_abun))
    
    df_summarized$type <- "community"
    
    return(df_summarized)
    
  }) %>% do.call("rbind",.) %>% dplyr::select(species_fmt, region, subj_region, rel_abun_sum, type)
  
  
  pathway_df_subset <- mtx_pathabun_stratified_TSS %>% dplyr::filter(module == reaction_module & 
                                                                       region %in% chosen_region) %>% 
    dplyr::rename(rel_abun_sum=rel_abun_CoPM) %>% 
    dplyr::select(species_fmt, region, subj_region, rel_abun_sum)
  
  
  pathway_df_subset$type <- reaction_module
  
  #Genus level and above classifications were abstracted to "unclassified", which means "unclassified at species level"
  pathway_df_subset$species_fmt <- gsub(pattern="unclassified", replacement="unclassified at species level", x=pathway_df_subset$species_fmt)
  
  pathway_df_subset$species_fmt <- ifelse(pathway_df_subset$species_fmt %in% c(chosen_species,"unclassified at species level"), 
                                          pathway_df_subset$species_fmt, "others")
  
  
  combined_df <- rbind(pathway_df_subset, mtx_df)  
  
  #Re-order factor levels to put "others" first.
  #renormalized RNA K2 abundances do not have "unclassified" any more
  taxa_to_add <- chosen_species[!chosen_species %in% c("unclassified" )]
  #combined_df$species_relabelled <- factor(combined_df$species_relabelled)
  
  combined_df<- combined_df %>% mutate(species_fmt= factor(species_fmt, 
                                                           levels = c(taxa_to_add,"others", "unclassified at species level")))
  
  if (manual_sample_order==TRUE){
    combined_df$subj_region <- factor(combined_df$subj_region, levels=subj_region_order)
  }
  
  combined_df$type <- factor(combined_df$type, levels=c( reaction_module,"community"))
  combined_df$region<- factor(combined_df$region, levels=chosen_region)
  
  
  ggplot(combined_df, 
         aes(x=subj_region, y=rel_abun_sum, fill=species_fmt)) + geom_col() + #coord_flip() #+ 
    facet_wrap(type~region, nrow=2, scales = "free")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x=element_blank()) + 
    ggtitle(paste0("Metatranscriptomes")) + scale_fill_manual(values = c24)
  
  
}

