##This function can either perform ORA or GSEA(preferred)
## Wrapper for clusterprofiler
tally_features_and_enrichments <- function(site_comparison, 
                                           term2gene_arg, term2name_arg, analysis_mode="ORA"){
  
  bact_DE_l2fc <- read_tsv(file=paste0("../data/DESeq2_out/bact_DESeq_",site_comparison,"_l2fc.tsv"), show_col_types = FALSE)
  
  fungi_DE_l2fc <- read_tsv(file=paste0("../data/DESeq2_out/fungi_DESeq_",site_comparison,"_l2fc.tsv"), show_col_types = FALSE)
  
  bact_n_upreg <- bact_DE_l2fc %>% dplyr::filter(padj < 0.05 & log2FoldChange > 0) %>% nrow()
  
  bact_n_downreg <- bact_DE_l2fc %>% dplyr::filter(padj < 0.05 & log2FoldChange < 0) %>% nrow()
  
  fungi_n_upreg <- fungi_DE_l2fc %>% dplyr::filter(padj < 0.05 & log2FoldChange > 0) %>% nrow()
  
  fungi_n_downreg <- fungi_DE_l2fc %>% dplyr::filter(padj < 0.05 & log2FoldChange < 0) %>% nrow()
  
  DE_tally_df <- data.frame(DE_comparison=site_comparison,
                            bact_n_upreg = bact_n_upreg,
                            bact_n_downreg = bact_n_downreg,
                            fungi_n_upreg = fungi_n_upreg,
                            fungi_n_downreg = fungi_n_downreg)
  
  DE_tally_df$total_n_upreg <- DE_tally_df$bact_n_upreg + DE_tally_df$fungi_n_upreg
  
  DE_tally_df$total_n_downreg <- DE_tally_df$bact_n_downreg + DE_tally_df$fungi_n_downreg
  
  ###
  bact_DE_l2fc_labelled <- bact_DE_l2fc %>%
    dplyr::select(feature,baseMean,log2FoldChange,lfcSE,pvalue,padj)
  
  bact_DE_l2fc_labelled$feature_from <- "bacteria"
  
  fungi_DE_l2fc_labelled <- fungi_DE_l2fc %>%
    dplyr::select(feature,baseMean,log2FoldChange,lfcSE,pvalue,padj)
  
  fungi_DE_l2fc_labelled$feature_from <- "fungi"
  
  ###
  #combine fungal and bacterial DE results
  combined_DE_res <- rbind(bact_DE_l2fc_labelled, 
                           fungi_DE_l2fc_labelled)
  
  combined_DE_res$DE_comparison <- site_comparison
  
  
  
  ###Clusterprofiler for upregulated genes
  
  #To do over-representation analysis 
  if (analysis_mode=="ORA"){
    bacteria_upreg_clusterprofiler_out <- run_clusterprofiler(changing_df = bact_DE_l2fc %>% dplyr::filter(padj < 0.05 & log2FoldChange > 0),
                                                              background_df=bact_DE_l2fc,
                                                              term2gene = term2gene_arg,
                                                              term2name = term2name_arg)
    
    bacteria_upreg_clusterprofiler_df <- as.data.frame(bacteria_upreg_clusterprofiler_out)
    
    if(nrow(bacteria_upreg_clusterprofiler_df >0)){
      bacteria_upreg_clusterprofiler_df <- bacteria_upreg_clusterprofiler_df %>% dplyr::filter(Description != "No_annotation")
    }
    
    
    
    fungi_upreg_clusterprofiler_out <- run_clusterprofiler(changing_df = fungi_DE_l2fc %>% 
                                                             dplyr::filter(padj < 0.05 & log2FoldChange > 0),
                                                           background_df=fungi_DE_l2fc, 
                                                           term2gene = term2gene_arg, term2name = term2name_arg)
    
    fungi_upreg_clusterprofiler_df <- as.data.frame(fungi_upreg_clusterprofiler_out)
    
    if(nrow(fungi_upreg_clusterprofiler_df >0)){
      fungi_upreg_clusterprofiler_df <- fungi_upreg_clusterprofiler_df %>% dplyr::filter(Description != "No_annotation")
      
    }
    
    
    ###Clusterprofiler for downregulated genes
    
    
    bacteria_downreg_clusterprofiler_out <- run_clusterprofiler(changing_df = bact_DE_l2fc %>% dplyr::filter(padj < 0.05 & log2FoldChange < 0),
                                                                background_df=bact_DE_l2fc,
                                                                term2gene = term2gene_arg,
                                                                term2name = term2name_arg)
    
    bacteria_downreg_clusterprofiler_df <- as.data.frame(bacteria_downreg_clusterprofiler_out)
    
    if(nrow(bacteria_downreg_clusterprofiler_df >0)){
      bacteria_downreg_clusterprofiler_df <- bacteria_downreg_clusterprofiler_df %>%  dplyr::filter(Description != "No_annotation")
    }
    
    
    
    fungi_downreg_clusterprofiler_out <- run_clusterprofiler(changing_df = fungi_DE_l2fc %>% dplyr::filter(padj < 0.05 & log2FoldChange < 0),
                                                             background_df=fungi_DE_l2fc,
                                                             term2gene = term2gene_arg,
                                                             term2name = term2name_arg)
    
    fungi_downreg_clusterprofiler_df <- as.data.frame(fungi_downreg_clusterprofiler_out)
    
    if(nrow(fungi_downreg_clusterprofiler_df >0)){
      fungi_downreg_clusterprofiler_df <- fungi_downreg_clusterprofiler_df %>%  dplyr::filter(Description != "No_annotation")
    }
    
    
    #For GO enrichment summary
    enrichment_tally_df <- data.frame(bact_upreg_n_enriched=nrow(bacteria_upreg_clusterprofiler_df),
                                      fungi_upreg_n_enriched=nrow(fungi_upreg_clusterprofiler_df),
                                      bact_downreg_n_enriched=nrow(bacteria_downreg_clusterprofiler_df),
                                      fungi_downreg_n_enriched=nrow(fungi_downreg_clusterprofiler_df),
                                      DE_comparison=site_comparison)
    
    
    output <- tibble::lst(DE_tally_df, combined_DE_res,
                          bacteria_upreg_clusterprofiler_out,
                          bacteria_upreg_clusterprofiler_df,
                          fungi_upreg_clusterprofiler_out,
                          fungi_upreg_clusterprofiler_df,
                          bacteria_downreg_clusterprofiler_out,
                          bacteria_downreg_clusterprofiler_df,
                          fungi_downreg_clusterprofiler_out,
                          fungi_downreg_clusterprofiler_df,
                          enrichment_tally_df)
  } else if (analysis_mode=="GSEA") {
    
    bact_DESeq_GSEA_results <- GSEA_from_DESeq(genelist=bact_DE_l2fc$log2FoldChange, 
                                               gene_names=bact_DE_l2fc$feature, term2gene = term2gene_arg, term2name = term2name_arg)
    
    fungi_DESeq_GSEA_results <- GSEA_from_DESeq(genelist=fungi_DE_l2fc$log2FoldChange, 
                                                gene_names=fungi_DE_l2fc$feature, term2gene = term2gene_arg, term2name = term2name_arg)
    
    output <- tibble::lst(DE_tally_df, 
                          combined_DE_res,
                          bact_DESeq_GSEA_results,
                          fungi_DESeq_GSEA_results)
  }
  
  
  return(output)
  
}

GSEA_from_DESeq <- function(genelist, gene_names, term2gene, term2name = KEGG_term2name){
  
  names(genelist) <- as.character(gene_names)
  # omit any NA values 
  genelist<-na.omit(genelist)
  #sort the genelist in decreasing order.
  genelist <- sort(genelist, decreasing = TRUE)
  
  set.seed(123)
  GSEA_out <- GSEA(geneList = genelist, TERM2GENE = term2gene, 
                   TERM2NAME = term2name, eps= 0, nPermSimple=10000, seed=TRUE)
  
  return(GSEA_out)
  
}

#https://www.biostars.org/p/467197/ for plotting pathway GSEA normalized enrichment scores in a bar plot

make_NES_barplot <- function(GSEA_df, title, x_axis_lab="Function/Process"){
  
  GSEA_df$direction <- ifelse(GSEA_df$NES > 0, "upregulated", "downregulated")
  cols <- c("downregulated" = "darkblue", "upregulated" = "red")
  
  GSEA_df$ID_and_Desc <- paste0(GSEA_df$ID," ",GSEA_df$Description)
  
  
  ggplot(GSEA_df, aes(reorder(ID_and_Desc, NES), NES, fill = direction)) +
    geom_col() +
    scale_fill_manual(values = cols) +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1)) +
    coord_flip() +
    labs(x=x_axis_lab, y="Normalized Enrichment Score") + ggtitle(title) + theme_classic()
  
}

#We want to find a set of GO slim metagenome terms that cover all/as much data as possible and represent them in the NES plot
#adapted from https://support.bioconductor.org/p/128407/
#go_ids is a vector of GO terms of interest
#Mode can be "MF" or "BP"
associate_GO_to_slim <- function(go_ids, 
                                 obo_path="../metadata/goslim_metagenomics.obo",
                                 mode="MF"){
  
  myCollection <- GOCollection(go_ids)
  
  # Retrieve GOslims from GO OBO file set
  slim <- getOBOCollection(obo_path)
  
  
  if (mode=="BP"){
    # Retrieve Biological Process (BP) GOslims
    slimdf <- goSlim(myCollection, slim, "BP", verbose)
    # List of GOslims and all GO IDs from `go_ids`
    gomap <- as.list(GOBPOFFSPRING[rownames(slimdf)])
    
  } else if (mode=="MF"){
    # Retrieve Biological Process (MF) GOslims
    slimdf <- goSlim(myCollection, slim, "MF", verbose)
    # List of GOslims and all GO IDs from `go_ids`
    gomap <- as.list(GOMFOFFSPRING[rownames(slimdf)])
  }
  
  
  # Maps `go_ids` to matching GOslims: lapply(gomap, intersect, ids(myCollection))
  
  # Append all mapped GO IDs to `slimdf`
  # `sapply` needed to apply paste() to create semi-colon delimited values
  slimdf$ids <- sapply(lapply(gomap, intersect, ids(myCollection)), paste, collapse=";")
  
  # Remove "character(0) string from "ids" column
  slimdf$ids[slimdf$ids == "character(0)"] <- ""
  
  slimdf[slimdf==""]<-NA
  
  slimdf$GO_slim_parent <- rownames(slimdf)
  
  slimdf$annotation_type <- mode
  
  # Add self-matching GOIDs to "ids" column, if not present
  #i.e the slimmed parent can be itself if it is a slimmed term also
  self_matching_parent_ids  <- intersect(go_ids,slimdf$GO_slim_parent)
  
  slimdf_nonself_subset <- slimdf %>% dplyr::filter(!GO_slim_parent %in% self_matching_parent_ids)
  slimdf_self_subset <- slimdf %>% dplyr::filter(GO_slim_parent %in% self_matching_parent_ids)
  
  
  
  slimdf_self_subset$ids <- ifelse(slimdf_self_subset$Count<=1, 
                                   slimdf_self_subset$GO_slim_parent,
                                   paste0(slimdf_self_subset$ids,";",slimdf_self_subset$GO_slim_parent))
  
  slimdf_self_subset$Count <- str_count(slimdf_self_subset$ids, "GO")
  
  
  slimdf_out <- rbind(slimdf_nonself_subset, slimdf_self_subset)
  slimdf_out <- slimdf_out %>% dplyr::rename(ID=ids) %>% dplyr::select(-Percent)
  
  return(slimdf_out)
  
}

#Function to color code NES bar plots by GO slim (parent categories)
#chosen_parent_ids can be a vector like c("transport","generation of precursor metabolites and energy"...)

prepare_NES_plot_input <- function(df, slimdf, chosen_parent_ids){
  
  slimdf_subset <- slimdf %>% tidyr::separate_rows(ID, sep = ";") %>% 
    dplyr::select(c("ID","Term")) %>% 
    dplyr::filter(Term %in% chosen_parent_ids) %>%
    dplyr::rename(GO_slim_parent_term=Term)
  
  
  df_out <- merge(df, slimdf_subset, by="ID", all.x=TRUE)
  df_out$GO_slim_parent_term <- as.factor(df_out$GO_slim_parent_term)
  
  return(df_out)
  
}


