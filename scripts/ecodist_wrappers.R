#wrappers for vegan and ecodist, to calculate contributional diversity

#calculate geometric mean of non-zero elements in a vector

geom_mean_non_zero <- function(x){
  
  exp(mean(log(x[x>0])))
}

counts_to_rel_abun <- function(vec){
  
  total_counts <- sum(vec)
  output <- vec/total_counts
  
  return(output)
  
}

#to calculate alpha diversity
get_contributional_alpha_diversity <- function(input_counts=mtx_KO_table,
                                               core_module_df = mtx_pathcounts_core_filt,
                                               diversity_index="simpson"){
  
  #Format column names of the mtx_KO_table
  colnames(input_counts) <- gsub("_genefamilies_Abundance","",colnames(input_counts))
  
  all_site_alpha_div_list <- lapply(c("Sc","Ch","Ac","Vf","Tw"), function(i){
    
    #vector of library IDs
    lib_for_site <- mtx_stats_chosen %>% dplyr::filter(region==i) %>% pull(LIBID)
    
    #get a vector of site specific core modules
    site_core_modules <- mtx_pathcounts_core_filt %>% 
      dplyr::filter(region==i) %>% pull(module)
    
    site_alpha_div_df <- lapply(site_core_modules, function(module_name, site=i){
      
      #Counts per million (from Humann3) for a module/pathway, excluding contributions from unclassified
      site_pathway_counts_df <- input_counts %>% 
        dplyr::filter(str_detect(string = .$module,pattern=paste0(module_name,"\\|")) &
                        str_detect(string = .$module, pattern="\\|unclassified",negate=TRUE)) %>%
        dplyr::select(all_of(lib_for_site))
      
      
      #before transpose, should remove any libraries for which all the species CPMs are all 0.
      site_pathway_counts_df <- site_pathway_counts_df[rowSums(site_pathway_counts_df)>0,colSums(site_pathway_counts_df)>0]
      
      alpha_diversity <- vegan::diversity(t(site_pathway_counts_df), index=diversity_index)
      
      alpha_diversity_df <- data.frame(LIBID=names(alpha_diversity),alpha_diversity_index=alpha_diversity)
      
      #Arithmetic mean of the Simpson index for a pathway in a skin site
      #We use arithmetic mean because Simpson score of 0 is valid
      #The Simpson index (1-D) is equals the probability that the two entities represent different types
      alpha_diversity_df$site_mean_alpha_diversity <- mean(alpha_diversity_df$alpha_diversity_index)
      
      alpha_diversity_df$module <- module_name
      alpha_diversity_df$region <- site
      
      return(alpha_diversity_df)
      
    }) %>% do.call("rbind",.)
    
    return(site_alpha_div_df)
  })
  
  output <- do.call("rbind",all_site_alpha_div_list)
  
  return(output)
}

#To calculate beta diversity (BC dissimilarity)
get_contributional_beta_diversity <- function(input_counts=mtx_KO_table,
                                              core_module_df = mtx_pathcounts_core_filt){
  
  #Format column names of the mtx_KO_table
  colnames(input_counts) <- gsub("_genefamilies_Abundance","",colnames(input_counts))
  
  all_site_beta_div_list <- lapply(c("Sc","Ch","Ac","Vf","Tw"), function(i){
    
    #vector of library IDs
    lib_for_site <- mtx_stats_chosen %>% dplyr::filter(region==i) %>% pull(LIBID)
    
    #get a vector of site specific core modules
    site_core_modules <- mtx_pathcounts_core_filt %>% 
      dplyr::filter(region==i) %>% pull(module)
    
    site_beta_div_df <- lapply(site_core_modules, function(module_name, site=i){
      
      #Counts per million (from Humann3) for a module/pathway, excluding contributions from unclassified
      site_pathway_counts_df <- input_counts %>% 
        dplyr::filter(str_detect(string = .$module,pattern=paste0(module_name,"\\|")) &
                        str_detect(string = .$module, pattern="\\|unclassified",negate=TRUE)) %>%
        dplyr::select(all_of(lib_for_site))
      
      
      #remove any rows or columns for which all values are 0.
      site_pathway_counts_df <- site_pathway_counts_df[rowSums(site_pathway_counts_df)>0,
                                                       colSums(site_pathway_counts_df)>0]
      
      #Convert to relative abundance for beta diversity calculations
      
      site_pathway_ra_df <- apply(site_pathway_counts_df, MARGIN=2, 
                                  counts_to_rel_abun, simplify = TRUE)
      
      #BC dissimilarity
      #rows as samples and columns as species
      beta_diversity <- ecodist::bcdist(t(site_pathway_ra_df), rmzero = TRUE) %>% as.matrix()
      
      beta_diversity_df  <- reshape2::melt(as.matrix(beta_diversity), varnames = c("row", "col"))
      
      #Arithmetic mean of the BC distance per module per site
      beta_diversity_df <- beta_diversity_df %>% dplyr::filter(col!=row)
      
      mean_beta_div <- mean(beta_diversity_df$value)
      
      mean_beta_div_df <- data.frame(site_mean_beta_diversity=mean_beta_div,
                                     module=module_name,
                                     region=site)
      
      
      return(mean_beta_div_df)
      
    }) %>% do.call("rbind",.)
    
    return(site_beta_div_df)
  })
  
  output <- do.call("rbind",all_site_beta_div_list)
  
  return(output)
}

####

facet_abundance_plot_fn_b <- function(mtx_list=rna_k2_minimizer_renorm, 
                                      mgx_list=k2_minimizer_renorm, 
                                      metadata = mtx_mgx_stats_chosen,
                                      chosen_species,
                                      chosen_region,
                                      manual_sample_order=FALSE,
                                      subj_region_order #a vector
){
  
  options(dplyr.summarise.inform = FALSE)
  
  mtx_libs_to_choose <- metadata %>% dplyr::filter(region==chosen_region) %>% pull(mtx_LIBID)
  mgx_libs_to_choose <- metadata %>% dplyr::filter(region==chosen_region) %>% pull(mgx_LIBID)
  
  
  mtx_df <- lapply(mtx_list[mtx_libs_to_choose], function(df){
    
    df <- merge(df %>% dplyr::rename(LIBID=mtx_LIBID), 
                metadata %>% dplyr::rename(LIBID=mtx_LIBID) %>% dplyr::select(LIBID, subj_region), by="LIBID")
    
    df$species_relabelled <- ifelse(df$k2_taxon %in% chosen_species, df$k2_taxon, "others")
    
    df_summarized <- df %>% group_by(LIBID, subj_region, species_relabelled) %>% summarise(rel_abun_sum=sum(rel_abun))
    
    df_summarized$type <- "MTX"
    
    return(df_summarized)
    
  }) %>% do.call("rbind",.)
  
  mgx_df <- lapply(mgx_list[mgx_libs_to_choose], function(df){
    
    df <- merge(df %>% dplyr::rename(LIBID=mgx_LIBID), 
                metadata %>% dplyr::rename(LIBID=mgx_LIBID) %>% dplyr::select(LIBID, subj_region), by="LIBID")
    
    df$species_relabelled <- ifelse(df$k2_taxon %in% chosen_species, df$k2_taxon, "others")
    
    df_summarized <- df %>% group_by(LIBID, subj_region, species_relabelled) %>% summarise(rel_abun_sum=sum(rel_abun))
    
    df_summarized$type <- "MGX"
    
    return(df_summarized)
    
  }) %>% do.call("rbind",.)
  
  combined_df <- rbind(mtx_df, mgx_df)  
  
  #Re-order factor levels to put "others" first.
  #renormalized RNA K2 abundances do not have "unclassified" any more
  taxa_to_add <- chosen_species[!chosen_species %in% c("unclassified" )]
  #combined_df$species_relabelled <- factor(combined_df$species_relabelled)
  
  combined_df<- combined_df %>% mutate(species_relabelled = factor(species_relabelled, 
                                                                   levels = c(taxa_to_add,"others")))
  
  if (manual_sample_order==TRUE){
    combined_df$subj_region <- factor(combined_df$subj_region, levels=subj_region_order)
  }
  
  
  ggplot(combined_df, 
         aes(x=subj_region, y=rel_abun_sum, fill=species_relabelled)) + geom_col() + #coord_flip() #+ 
    facet_wrap(~type, ncol=1, nrow=2)+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x=element_blank()) + 
    ggtitle(paste0("Taxa on ",chosen_region)) + scale_fill_manual(values = c24)
  
  
}
