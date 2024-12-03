
tally_signif_DE <- function(input_df){
  
  #upregulated 
  upreg_df <- input_df %>% dplyr::filter(log2FoldChange > 0 & padj < 0.05)
  upreg_tally <- length(unique(upreg_df$pangene))
  
  #Downregulated
  downreg_df <- input_df %>% dplyr::filter(log2FoldChange < 0 & padj < 0.05)
  downreg_tally <- length(unique(downreg_df$pangene))
  
  
  output <- data.frame(comparison=unique(input_df$comparison),
                       num_upreg_features=upreg_tally, 
                       num_downreg_features=downreg_tally)
  
  return(output)
  
} 

#For pan-microbial comparison
run_DESeq_with_inputs <- function(mtx_count_list, mgx_count_list, cond, ref_cond, 
                                  size_factor_function = "poscounts", annot){
  
  libs_to_pull <- eligible_mtx_metadata %>% dplyr::filter(region %in% c(cond, ref_cond)) %>% pull(LIBID)
  
  ##############################################
  ## prepare input read count matrix for DESeq2#
  ##############################################
  
  #Obtain count data for eligible mtx and mgx
  rawcounts <- c(mtx_count_list[libs_to_pull], mgx_count_list[libs_to_pull] )%>% do.call("rbind",.)
  
  count_matrix <- pivot_wider(rawcounts %>% 
                                dplyr::select(feature, read_count_sum, LIBID), names_from = LIBID, values_from = read_count_sum) %>%
    as.data.frame()
  
  rownames(count_matrix) <- count_matrix$feature
  
  count_matrix <- count_matrix %>% dplyr::select(-feature) %>% as.matrix()
  
  #This matrix operation is much faster than mutate_all(~replace(.,is.na(.),0)) for tibbles/df
  count_matrix[is.na(count_matrix)] <- 0 
  
  ##############################################
  ## prepare metadata for DESeq2               #
  ##############################################
  
  #the row names of metadata should automatically be ordered in a way that matches the input matrix
  selected_metadata <- eligible_mtx_metadata %>% dplyr::filter(region %in% c(cond, ref_cond)) %>% 
    dplyr::select(LIBID, subject, region, mol_type) %>% as.data.frame()
  
  rownames(selected_metadata) <- selected_metadata$LIBID
  
  selected_metadata <- selected_metadata %>% dplyr::select(-LIBID)
  
  selected_metadata$subject <- as.factor(selected_metadata$subject)
  selected_metadata$region <- as.factor(selected_metadata$region)
  selected_metadata$mol_type <- as.factor(selected_metadata$mol_type)
  
  #reorder rows of the metadata
  idx <- match(colnames(count_matrix), rownames(selected_metadata))
  selected_metadata <- selected_metadata[idx,]
  
  #a vector of LIBIDs
  RNA_libs <- selected_metadata %>% dplyr::filter(mol_type=="RNA") %>% row.names()
  DNA_libs <- selected_metadata %>% dplyr::filter(mol_type=="DNA") %>% row.names()
  
  ##############################################
  ## Run DESeq2                                #
  ##############################################
  
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                                colData = selected_metadata, 
                                design = ~subject+region+mol_type+region:mol_type)
  
  #Pre-filtering to keep rows with a minimum median read count for transcripts AND gene counts.
  #Caution: You cannot naively run "keep <- rowMedians(counts(dds)) >= 10", because the DNA read counts in the matrix might obscure this 
  #There are cases where the gene counts are extremely low but transcript counts are high e.g. Sc vs Vf MET_02106687
  
  
  mat <- counts(dds)
  rna_mat <- mat[, colnames(mat) %in% RNA_libs]
  dna_mat <- mat[, colnames(mat) %in% DNA_libs]
  
  keep <- (rowMedians(rna_mat) >= 10) & (rowMedians(dna_mat) >= 10)
  dds <- dds[keep,]
  
  #####Estimate size factors separately for mtx and mgx data (subsets of the dds) https://support.bioconductor.org/p/67455/#67498
  sf <- numeric(ncol(dds))  #a numeric vector
  idx1 <- dds$mol_type == "RNA"
  sf[ idx1 ] <- estimateSizeFactorsForMatrix(counts(dds)[ , idx1], type=size_factor_function) #poscounts by default for zero inflated metagenomic data
  idx2 <- dds$mol_type == "DNA"
  sf [idx2 ] <- estimateSizeFactorsForMatrix(counts(dds)[ , idx2], type=size_factor_function)
  sizeFactors(dds) <- sf
  
  
  #dds <- estimateSizeFactors(dds, type = size_factor_function) #only if you want to estimate size factors using combined data for the two assays
  
  #factor relevel such that ref_cond and "DNA" are always the reference conditions
  
  dds$region <- relevel(dds$region, ref_cond)
  dds$mol_type <- relevel(dds$mol_type, "DNA")
  
  normalized_counts <- counts(dds, normalized=TRUE)
  
  #See error and fix:https://support.bioconductor.org/p/98634/
  #See DESeq2 vignette about blind dispersion estimation
  
  vsd <- vst(dds, blind=FALSE, nsub=sum( rowMeans( counts(dds, normalized=TRUE)) > 5 ))  
  
  vsd_mat <- assay(vsd)
  
  vsd_cor <- cor(vsd_mat)
  
  #vsd_cor_heatmap <- pheatmap(vsd_cor, annotation=select(selected_metadata, region))
  
  #plotPCA(vsd, intgroup = "region")
  
  
  #dds <- DESeq(dds, reduced = ~subject+mol_type+region, test="LRT")   #Only if you want the likelihood ratio test
  #ddsClean <- dds[which(mcols(dds)$fullBetaConv),]  #Removing rows which do not converge in beta 
  
  
  dds <- DESeq(dds)   #wald test
  ddsClean <- dds[which(mcols(dds)$betaConv),]  #Removing rows which do not converge in beta 
  
  
  #plotDispEsts(ddsClean)
  
  DE_results <- results(dds,name = paste0("region",cond,".mol_typeRNA"), alpha=0.05)
  
  DE_results_clean <- results(ddsClean,name = paste0("region",cond,".mol_typeRNA"), alpha=0.05)
  
  #For ashr, if res is provided, then coef and contrast are ignored.
  DE_results_shrunken <- lfcShrink(dds, type="ashr", res=DE_results)
  
  DE_results_clean_shrunken <- lfcShrink(ddsClean, type="ashr", res=DE_results_clean)
  
  #plotMA(DE_results_shrunken, ylim = c(-20,20))
  
  #summary(DE_results_shrunken)
  
  DE_results_shrunken_df <- data.frame(DE_results_shrunken)
  
  DE_results_clean_shrunken_df <- data.frame(DE_results_clean_shrunken)
  
  #Add annotations like bact_og_metadata
  
  DE_results_shrunken_df$feature <- rownames(DE_results_shrunken_df)
  
  colnames(annot)[1] <- "feature"
  
  DE_results_shrunken_df <- merge(DE_results_shrunken_df, annot, by = "feature", all.x=TRUE)
  
  DE_results_clean_shrunken_df$feature <- rownames(DE_results_clean_shrunken_df)
  
  DE_results_clean_shrunken_df <- merge(DE_results_clean_shrunken_df, annot, by = "feature", all.x=TRUE)
  
  ##############################################
  ## Save outputs                              #
  ##############################################
  #results are presented as cond vs ref_cond
  outputs <- tibble::lst(skin_site=cond, ref_skin_site=ref_cond, 
                         count_matrix, keep, selected_metadata,
                         dds, ddsClean,
                         normalized_counts,
                         vsd, vsd_mat, vsd_cor,
                         DE_results, DE_results_shrunken,
                         DE_results_clean, DE_results_clean_shrunken,
                         DE_results_shrunken_df,
                         DE_results_clean_shrunken_df)
  
  return(outputs)
}

###

run_DESeq_comparisons <- function(site_to_compare, site_reference){
  
  #bacterial OGs
  
  bact_OG_DE <- run_DESeq_with_inputs(mtx_count_list=mtx_bact_counts, 
                                      mgx_count_list=mgx_bact_counts, 
                                      cond=site_to_compare, ref_cond=site_reference,
                                      annot=bact_og_metadata)
  
  #fungal OGs
  
  fungal_OG_DE <- run_DESeq_with_inputs(mtx_count_list=mtx_fungi_counts, 
                                        mgx_count_list=mgx_fungi_counts, 
                                        cond=site_to_compare, ref_cond=site_reference,
                                        annot=fungi_og_metadata)
  
  
  #Other features with no match to bacterial or fungal OG
  
  #pangene_uniref_DE <- run_DESeq_with_inputs(mtx_count_list=mtx_uniref_pangene_counts, 
  #                      mgx_count_list=mgx_uniref_pangene_counts, 
  #                      cond=site_to_compare, ref_cond=site_reference,
  #                      annot=mtx_uniref_pangene_metadata)
  
  #list of lists
  outputs <- tibble::lst(bact_OG_DE, fungal_OG_DE)
  return(outputs)
  
}


run_DESeq_with_inputs_species_specific <- function(input_list, 
                                                   input_metadata=DESeq_metadata,
                                                   lib_vec,
                                                   run_test=FALSE,
                                                   batch_correct=TRUE,
                                                   size_factor_function = "poscounts"){
  
  ##############################################
  ## prepare input read count matrix for DESeq2#
  ##############################################
  
  #Obtain count data for eligible mtx and mgx
  rawcounts <- input_list[lib_vec] %>% do.call("rbind",.)
  
  count_matrix <- pivot_wider(rawcounts %>% 
                                dplyr::select(pangene, unpaired_read_count, LIBID), names_from = LIBID, 
                              values_from = unpaired_read_count) %>%
    as.data.frame()
  
  rownames(count_matrix) <- count_matrix$pangene
  
  count_matrix <- count_matrix %>% dplyr::select(-pangene) %>% as.matrix()
  
  #This matrix operation is much faster than mutate_all(~replace(.,is.na(.),0)) for tibbles/df
  count_matrix[is.na(count_matrix)] <- 0 
  
  ##############################################
  ## prepare metadata for DESeq2               #
  ##############################################
  
  #the row names of metadata should automatically be ordered in a way that matches the input matrix
  selected_metadata <- input_metadata[lib_vec,]
  selected_metadata$condition <- as.factor(selected_metadata$condition)
  selected_metadata$expt_class <- as.factor(selected_metadata$expt_class)
  
  tested_conditions <- unique(selected_metadata$condition)
  
  #reorder rows of the metadata
  idx <- match(colnames(count_matrix), rownames(selected_metadata))
  selected_metadata <- selected_metadata[idx,]
  
  
  ##############################################
  ## Run DESeq2                                #
  ##############################################
  
  if(batch_correct==TRUE & run_test==FALSE){
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                                  colData = selected_metadata, 
                                  design = ~study+condition)
    
  }else if(batch_correct==FALSE & run_test==FALSE){
    
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                                  colData = selected_metadata, 
                                  design = ~condition)
    
  }else if(batch_correct==TRUE & run_test==TRUE){  #Comparing in vivo vs in vitro
    
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                                  colData = selected_metadata, 
                                  design = ~study+expt_class)
    
  }else if(batch_correct==FALSE & run_test==TRUE){
    
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                                  colData = selected_metadata, 
                                  design = ~expt_class)
  }
  
  
  #Pre-filtering to keep rows with a minimum median read count for transcripts.
  
  mat <- counts(dds)
  
  keep <- (rowMedians(mat) >= 10)
  dds <- dds[keep,]
  
  #####Estimate size factors
  dds <- estimateSizeFactors(dds,  type = size_factor_function)
  
  #in_vitro will be the reference condition for the actual DESeq test
  dds$expt_class <- relevel(dds$expt_class, "in_vitro")
  
  normalized_counts <- counts(dds, normalized=TRUE)
  
  #
  
  vsd <- vst(dds, blind=TRUE)  
  
  vsd_mat <- assay(vsd)
  
  vsd_cor <- cor(vsd_mat)
  
  #vsd_cor_heatmap <- pheatmap(vsd_cor, annotation=select(selected_metadata, region))
  
  #plotPCA(vsd, intgroup = "region")
  
  if (run_test==FALSE){
    outputs <- tibble::lst(count_matrix, keep, selected_metadata,
                           dds,
                           normalized_counts,
                           vsd, vsd_mat, vsd_cor)
  }else if (run_test==TRUE){
    dds <- DESeq(dds)
    
    DE_results <- results(dds,
                          contrast=c("expt_class","in_vivo","in_vitro"),
                          alpha=0.05)
    
    DE_results_shrunken <- lfcShrink(dds,
                                     type="ashr",
                                     res=DE_results)
    
    DE_results_shrunken_df <- data.frame(DE_results_shrunken)
    
    DE_results_shrunken_df$pangene <- rownames(DE_results_shrunken_df)
    
    outputs <- tibble::lst(lib_vec, count_matrix, keep, selected_metadata,
                           tested_conditions,
                           dds,
                           normalized_counts,
                           vsd, vsd_mat, vsd_cor,
                           DE_results, DE_results_shrunken,
                           DE_results_shrunken_df)
    
  }
  
  return(outputs)
}


#This function selects for mtx libraries corresponding to sites of choice.
#Assumes input is already a count matrix from Salmon quant
#Species can be "Malassezia_restricta", "Malassezia_globosa" etc
run_DESeq_with_salmon_matrix <- function(input_matrix, 
                                         input_metadata=mtx_stats_chosen,
                                         lib_vec,
                                         cond, 
                                         ref_cond, 
                                         run_test=FALSE,
                                         control_for_individuals=TRUE,
                                         size_factor_function = "poscounts",
                                         species="Malassezia_restricta"){
  
  libs_to_pull <- input_metadata %>% dplyr::filter(LIBID %in% lib_vec) %>% pull(LIBID)
  
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
  
  
  ##############################################
  ## prepare input read count matrix for DESeq2#
  ##############################################
  
  #Subset read count matrix with the selected libraries AND species.
  
  chosen_matrix <- input_matrix[grep(loci, row.names(input_matrix)),libs_to_pull]
  
  #Convert to matrix of integers
  
  chosen_matrix <- round(chosen_matrix)
  mode(chosen_matrix) <- "integer"
  
  ##############################################
  ## prepare metadata for DESeq2               #
  ##############################################
  
  
  selected_metadata <- input_metadata %>% dplyr::filter(LIBID %in% lib_vec) %>% 
    dplyr::select(LIBID, subject, region) %>% as.data.frame()
  
  rownames(selected_metadata) <- selected_metadata$LIBID
  
  selected_metadata <- selected_metadata %>% dplyr::select(-LIBID)
  
  selected_metadata$subject <- as.factor(selected_metadata$subject)
  selected_metadata$region <- as.factor(selected_metadata$region)
  
  #the row names of metadata should automatically be ordered in a way that matches the input matrix
  #reorder rows of the metadata
  idx <- match(colnames(chosen_matrix), rownames(selected_metadata))
  selected_metadata <- selected_metadata[idx,]
  
  
  #a vector of LIBIDs
  RNA_libs <- selected_metadata %>% row.names()
  
  ##############################################
  ## Run DESeq2                                #
  ##############################################
  #https://www.biostars.org/p/412320/ DESeqDataSetFromTximport rounds salmon matrix to integers, unlike DESeqDataSetFromMatrix
  #Here I am using the salmon read matrix straight off, so I have manually done the rounding
  
  #https://support.bioconductor.org/p/9137115/
  #If your design is not confounded (e.g. you have balanced your conditions of interest within batches or other nuisance variables) then you can do, e.g. ~batch + condition or       #~batch + nuisance1 + nuisance2 + condition.
  
  if(control_for_individuals==TRUE){
    dds <- DESeqDataSetFromMatrix(countData = chosen_matrix, 
                                  colData = selected_metadata, 
                                  design = ~subject+region)
  } else if (control_for_individuals==FALSE){
    dds <- DESeqDataSetFromMatrix(countData = chosen_matrix, 
                                  colData = selected_metadata, 
                                  design = ~region)  
    
  }
  
  #Pre-filtering to keep rows with a minimum median read count <= 10 for genes
  
  mat <- counts(dds)
  
  keep <- (rowMedians(mat) >= 10)
  dds <- dds[keep,]
  
  #####Estimate size factors
  dds <- estimateSizeFactors(dds,  type = size_factor_function)
  
  dds$region <- relevel(dds$region, ref_cond)
  
  normalized_counts <- counts(dds, normalized=TRUE)
  
  #See error and fix:https://support.bioconductor.org/p/98634/
  #See DESeq2 vignette about blind dispersion estimation
  #https://www.biostars.org/p/428369/
  
  vsd <- vst(dds, blind=FALSE)  
  
  vsd_mat <- assay(vsd)
  
  vsd_cor <- cor(vsd_mat)
  
  
  if (run_test==FALSE){
    outputs <- tibble::lst(chosen_matrix, keep, selected_metadata,
                           dds,
                           normalized_counts,
                           vsd, vsd_mat, vsd_cor)
  } else if (run_test==TRUE){
    
    dds <- DESeq(dds)   #wald test
    
    DE_results <- results(dds,contrast=c("region",cond,ref_cond),
                          alpha=0.05)
    
    DE_results_shrunken <- lfcShrink(dds,
                                     type="ashr",
                                     res=DE_results)
    
    DE_results_shrunken_df <- data.frame(DE_results_shrunken)
    
    DE_results_shrunken_df$gene <- rownames(DE_results_shrunken_df)
    
    outputs <- tibble::lst(lib_vec, chosen_matrix, keep, selected_metadata,
                           skin_site=cond, ref_skin_site=ref_cond,
                           dds,
                           normalized_counts,
                           vsd, vsd_mat, vsd_cor,
                           DE_results, DE_results_shrunken,
                           DE_results_shrunken_df)
  }
  return(outputs)
}


# GSEA wrappers


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


make_NES_barplot_rotate <- function(GSEA_df, title){
  
  GSEA_df$direction <- ifelse(GSEA_df$NES > 0, "upregulated", "downregulated")
  cols <- c("downregulated" = "darkblue", "upregulated" = "red")
  
  
  ggplot(GSEA_df, aes(reorder(Description, NES), NES, fill = direction)) +
    geom_col() +
    scale_fill_manual(values = cols) +
    labs(x="KEGG Pathway", y="Normalized Enrichment Score") + ggtitle(title) + theme_classic() +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1))
  
}





make_custom_term2gene_df <- function(input_annot, pangene_names, annotation_used="KEGG"){
  
  if(annotation_used=="KEGG"){
    colnames(input_annot)[7] <- "KEGG_Pathway"
    
    annot_KO <- input_annot %>% dplyr::select(propan_clusterID, KEGG_Pathway) %>% dplyr::filter(KEGG_Pathway != "-")
    
    annot_KO <- annot_KO %>% tidyr::separate_rows(KEGG_Pathway, sep = ";")
    
    annot_KO <- annot_KO %>% dplyr::filter(str_detect(string=annot_KO$KEGG_Pathway, pattern="map"))
    
    annot_KO$KEGG_Pathway_term <- str_match(string=annot_KO$KEGG_Pathway, pattern="map[0-9]*")[,1]
    
    pangene_names$propan_clusterID <- str_match(string=pangene_names$pangene, pattern="cluster_[0-9]*")[,1]
    
    pangene_names_filt <- merge(pangene_names, annot_KO, by ="propan_clusterID")
    
    #to=pangene_names_filt$pangene is also possible
    output_term2gene <- data.frame(from=pangene_names_filt$KEGG_Pathway_term, 
                                   to=pangene_names_filt$propan_clusterID) %>% unique()
    
    #remove rows corresponding to pathways for organismal systems or human diseases
    
    output_term2gene <- output_term2gene %>% dplyr::filter(!from %in% c(KEGG_hum_disease_pathways$map_ID,
                                                                        KEGG_org_sys_pathways$map_ID))
    
  } else if (annotation_used=="GO"){
    
    annot_GO <- input_annot %>% dplyr::select(propan_clusterID, GOs) %>% 
      dplyr::filter(GOs != "-")
    
    annot_GO <- annot_GO %>% tidyr::separate_rows(GOs, sep = ",") %>% unique()
    
    output_term2gene <- annot_GO %>% dplyr::rename(from="GOs", to="propan_clusterID") %>% dplyr::select(from, to)
    
  }
  
  
  return(output_term2gene)
  
}


gmean_of_non_NA <- function(x, threshold=3){
  
  row_wise_entries <- sum(!is.na(x))
  #The geometric mean is more robust to outliers
  if(row_wise_entries >= threshold){
    vector_no_NA <- x[!is.na(x)]
    output <- exp(mean(log(vector_no_NA)))
  } else {
    output <- NA
  }
  return(output)
  
}

#function to format dataframe for rowwise geometric mean calculation
#For ECnum

fmt_df_for_gmean <- function(df){
  
  df_wide <- pivot_wider(df %>% 
                           dplyr::select(ECnum, TPM, LIBID), names_from = LIBID, 
                         values_from = TPM) %>%
    as.data.frame()
  
  #E.C numbers in row names
  df_wide_fmt <- df_wide
  
  rownames(df_wide_fmt) <- df_wide_fmt$ECnum
  df_wide_fmt <- df_wide_fmt %>% dplyr::select(-ECnum)
  
  return(df_wide_fmt)
  
}

#For model_gene_id

fmt_df_for_gmean_alt <- function(df){
  
  df_wide <- pivot_wider(df %>% 
                           dplyr::select(model_gene_id, TPM, LIBID), names_from = LIBID, 
                         values_from = TPM) %>%
    as.data.frame()
  
  #E.C numbers in row names
  df_wide_fmt <- df_wide
  
  rownames(df_wide_fmt) <- df_wide_fmt$model_gene_id
  df_wide_fmt <- df_wide_fmt %>% dplyr::select(-model_gene_id)
  
  return(df_wide_fmt)
  
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




make_NES_barplot_with_slim <- function(GSEA_df, title, x_axis_lab="Function/Process"){
  
  GSEA_df$ID_and_Desc <- paste0(GSEA_df$ID," ",GSEA_df$Description)
  
  
  ggplot(GSEA_df, aes(reorder(ID_and_Desc, NES), NES, fill = GO_slim_parent_term)) +
    geom_col(colour="black",position="dodge") +
    scale_fill_manual(values = c23) +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1)) +
    coord_flip() +
    labs(x=x_axis_lab, y="Normalized Enrichment Score") + ggtitle(title) + theme_classic()
  
}


#####Wrappers for QC or analyzing different BACTERIAL (not fungi) species and running site specific DE on them.

#The following function is for running DESeq from a input list object derived from pileup.sh's outputs
#Updated to analyze read counts for features that are aggregated into original propan_clusterIDs
run_DESeq_with_input_list <- function(input_list, 
                                      input_metadata=mtx_stats_chosen,
                                      cond, 
                                      ref_cond,
                                      lib_vec,
                                      run_test=FALSE, 
                                      control_for_individuals=TRUE,
                                      manual_vst_nsub = FALSE,
                                      size_factor_function = "poscounts"
){
  
  ##############################################
  ## prepare input read count matrix for DESeq2#
  ##############################################
  
  #Obtain count data for eligible mtx and mgx
  rawcounts <- input_list[lib_vec] %>% do.call("rbind",.)
  
  count_matrix <- pivot_wider(rawcounts %>% 
                                dplyr::select(propan_clusterID, unpaired_read_count_sum, LIBID), names_from = LIBID, 
                              values_from = unpaired_read_count_sum) %>%
    as.data.frame()
  
  rownames(count_matrix) <- count_matrix$propan_clusterID
  
  count_matrix <- count_matrix %>% dplyr::select(-propan_clusterID) %>% as.matrix()
  
  #This matrix operation is much faster than mutate_all(~replace(.,is.na(.),0)) for tibbles/df
  count_matrix[is.na(count_matrix)] <- 0 
  
  ##############################################
  ## prepare metadata for DESeq2               #
  ##############################################
  
  
  selected_metadata <- input_metadata %>% dplyr::filter(LIBID %in% lib_vec) %>% 
    dplyr::select(LIBID, subject, region) %>% as.data.frame()
  
  rownames(selected_metadata) <- selected_metadata$LIBID
  
  selected_metadata <- selected_metadata %>% dplyr::select(-LIBID)
  
  selected_metadata$subject <- as.factor(selected_metadata$subject)
  selected_metadata$region <- as.factor(selected_metadata$region)
  
  #reorder rows of the metadata
  idx <- match(colnames(count_matrix), rownames(selected_metadata))
  selected_metadata <- selected_metadata[idx,]
  
  
  ##############################################
  ## Run DESeq2                                #
  ##############################################
  
  if (control_for_individuals==TRUE){
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                                  colData = selected_metadata, 
                                  design = ~ subject + region)
  } else if (control_for_individuals==FALSE){
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                                  colData = selected_metadata, 
                                  design = ~ region)
    
  }
  
  
  #Pre-filtering to keep rows with a minimum median read count for transcripts/propan clusters.
  
  mat <- counts(dds)
  
  keep <- (rowMedians(mat) >= 10)
  dds <- dds[keep,]
  
  #####Estimate size factors
  dds <- estimateSizeFactors(dds,  type = size_factor_function)
  
  dds$region <- relevel(dds$region, ref_cond)
  
  normalized_counts <- counts(dds, normalized=TRUE)
  
  # set manual_vst_nsub == TRUE if number of rows in dataset is lower than the DESeq2 default of 1000
  if (manual_vst_nsub == FALSE){
    vsd <- vst(dds, blind=FALSE) 
  } else if (manual_vst_nsub == TRUE){
    vsd <- vst(dds, blind=FALSE, nsub=sum( rowMedians(normalized_counts) >= 10 )) 
  }
  
  
  vsd_mat <- assay(vsd)
  
  vsd_cor <- cor(vsd_mat)
  
  
  if (run_test==FALSE){
    outputs <- tibble::lst(count_matrix, keep, selected_metadata,
                           dds,
                           normalized_counts,
                           vsd, vsd_mat, vsd_cor)
  }else if (run_test==TRUE){
    dds <- DESeq(dds)
    
    ##Removing rows which do not converge in beta 
    dds <- dds[which(mcols(dds)$betaConv),]
    
    DE_results <- results(dds,
                          contrast=c("region",cond,ref_cond),
                          alpha=0.05)
    
    DE_results_shrunken <- lfcShrink(dds,
                                     type="ashr",
                                     res=DE_results)
    
    DE_results_shrunken_df <- data.frame(DE_results_shrunken)
    
    DE_results_shrunken_df$propan_clusterID <- rownames(DE_results_shrunken_df)
    
    outputs <- tibble::lst(lib_vec, count_matrix, keep, 
                           cond, ref_cond,
                           selected_metadata,
                           dds,
                           normalized_counts,
                           vsd, vsd_mat, vsd_cor,
                           DE_results, DE_results_shrunken,
                           DE_results_shrunken_df)
    
  }
  
  return(outputs)
}



#The updated wrapper function now simply accepts a pre-processed input list of pangene counts
#Make sure the list of counts are species specific
#The input list comprises of dfs with three columns: LIBID, propan_clusterID, unpaired_read_count_sum
#The input list is a named list, with each Library ID being a name

QC_bacteria_RNA_sites <- function(list_of_counts,
                                  species,
                                  site, ref_site,
                                  count_threshold,
                                  specify_vst_nsub = FALSE,
                                  shared_individuals_bet_sites=TRUE){
  
  site_libs <- mtx_stats_chosen %>% dplyr::filter(region %in% c(site,ref_site)) %>% pull(LIBID)
  
  species_site_metadata <- lapply(rna_k2_minimizer[site_libs], function(df){
    
    metadata_out <- df %>% dplyr::filter(k2_taxon == species) %>% dplyr::rename(LIBID=mtx_LIBID)
    
    metadata_out <- merge(metadata_out, mtx_stats_chosen %>% 
                            dplyr::select(c("LIBID","subj_region", "subject","region","low_conc")),
                          by="LIBID")
    
    return(metadata_out)
  }) %>% do.call("rbind",.)
  
  species_site_mtx_counts <- list_of_counts[site_libs]
  
  control_for_individuals_boolean <- shared_individuals_bet_sites
  
  DE_PCA_input <- run_DESeq_with_input_list(input_list = species_site_mtx_counts,
                                            lib_vec=species_site_metadata %>% 
                                              dplyr::filter(paired_counts >= count_threshold) %>%  pull(LIBID),
                                            cond=site,
                                            ref_cond=ref_site,
                                            control_for_individuals = control_for_individuals_boolean,
                                            manual_vst_nsub = specify_vst_nsub,
                                            run_test = FALSE)
  
  PCA_plot <- plotPCA(DE_PCA_input[["vsd"]], intgroup = "region", returnData=FALSE) + geom_text_repel(aes(label = name))
  
  output <- tibble::lst(species, 
                        site, 
                        ref_site,
                        count_threshold,
                        species_site_metadata,
                        species_site_mtx_counts,
                        PCA_plot)  
  
  return(output)
  
  
}

#The updated wrapper function now simply accepts a pre-processed input list of pangene counts
#Make sure the list of counts are species specific
#The input list comprises of dfs with three columns: LIBID, propan_clusterID, unpaired_read_count_sum
#The input list is a named list, with each Library ID being a name
#This function performs DESeq2 analysis and overrepresentation analysis (ORA) using clusterprofiler
#Aug 2024 update: GSEA is preferred over ORA.
compare_bacteria_RNA_sites <- function(list_of_counts, species_anno, species_pangene_names, 
                                       species, site, ref_site,
                                       count_threshold,
                                       specify_vst_nsub = FALSE,
                                       shared_individuals_bet_sites=TRUE){
  
  
  site_libs <- mtx_stats_chosen %>% dplyr::filter(region %in% c(site,ref_site)) %>% pull(LIBID)
  
  species_site_metadata <- lapply(rna_k2_minimizer[site_libs], function(df){
    
    metadata_out <- df %>% dplyr::filter(k2_taxon == species) %>% dplyr::rename(LIBID=mtx_LIBID)
    
    metadata_out <- merge(metadata_out, mtx_stats_chosen %>% dplyr::select(c("LIBID","subj_region", "subject","region","low_conc")),
                          by="LIBID")
    
    return(metadata_out)
  }) %>% do.call("rbind",.)
  
  
  species_site_mtx_counts <- list_of_counts[site_libs]
  
  
  control_for_individuals_boolean <- shared_individuals_bet_sites
  
  
  species_site_filt_DESeq <-  run_DESeq_with_input_list(input_list = species_site_mtx_counts,
                                                        lib_vec=species_site_metadata %>% 
                                                          dplyr::filter(paired_counts >= count_threshold) %>%
                                                          pull(LIBID),
                                                        cond=site,
                                                        ref_cond=ref_site,
                                                        control_for_individuals = control_for_individuals_boolean,
                                                        manual_vst_nsub = specify_vst_nsub,
                                                        run_test = TRUE)
  
  species_site_DE_res <- species_site_filt_DESeq[["DE_results_shrunken_df"]] %>% dplyr::rename(feature=propan_clusterID)
  
  ### Run KEGG enrichment analysis 
  
  species_term2gene <- make_custom_term2gene_df(input_annot=species_anno,
                                                pangene_names = species_pangene_names)
  
  species_site_DE_KEGG_upreg <- run_clusterprofiler(changing_df = species_site_DE_res %>% dplyr::filter(padj < 0.05 & log2FoldChange > log(1.5)/log(2)),
                                                    background_df=species_site_DE_res,
                                                    term2gene = species_term2gene)
  
  species_site_DE_KEGG_downreg <- run_clusterprofiler(changing_df = species_site_DE_res %>% dplyr::filter(padj < 0.05 & log2FoldChange < -(log(1.5)/log(2))),
                                                      background_df=species_site_DE_res,
                                                      term2gene = species_term2gene)
  
  species_site_DE_KEGG_upreg_df <- as.data.frame(species_site_DE_KEGG_upreg)
  
  species_site_DE_KEGG_downreg_df <- as.data.frame(species_site_DE_KEGG_downreg) 
  
  
  
  output <- tibble::lst(species, 
                        site, 
                        ref_site,
                        count_threshold,
                        species_site_metadata,
                        species_site_mtx_counts,
                        control_for_individuals_boolean,
                        species_site_filt_DESeq,
                        species_site_DE_res,
                        species_term2gene,
                        species_site_DE_KEGG_upreg,
                        species_site_DE_KEGG_upreg_df,
                        species_site_DE_KEGG_downreg,
                        species_site_DE_KEGG_downreg_df
  ) 
  
  return(output)
  
} 


#Changing df contains the changing genes (up or downregulated)
#background_df contains all detected genes in the DESeq2 results
run_clusterprofiler <- function(changing_df, background_df, term2gene, term2name=KEGG_term2name){
  
  #the "universe"/background of genes, is all the detected genes in the DEseq2 results and NOT all possible eggnog entries in the bacterial pangenome
  term2gene_filt <- term2gene %>% dplyr::filter(to %in% background_df$feature)
  
  #We actually want to define the universe as all detected genes, which is bigger than the intersect of the background_df and term2gene (previous line)
  other_background <- background_df %>% dplyr::filter(!feature %in% term2gene_filt$to)
  
  if (nrow(other_background) > 0){
    
    other_background_term2gene <- data.frame(from="map00000", to=other_background$feature)
    
    term2gene_universe <- rbind(term2gene_filt, other_background_term2gene)
    
  } else {
    
    term2gene_universe <- term2gene_filt
    
  }
  
  #Add dummy term to term2name
  term2name_dummy <- data.frame(from="map00000", to="No_annotation")
  term2name <- rbind(term2name, term2name_dummy)
  
  
  #changing_df contains the genes that are significantly different between conditions
  
  #clusterProfiler ORA analysis requires a ranked gene list, which contains three features
  #numeric vector: fold change or other type of numerical variable
  #named vector: every number has a name, the corresponding gene ID
  #sorted vector: number should be sorted in decreasing order
  
  #prepare the numeric vector of l2fc of significant DE genes
  geneList <- changing_df$log2FoldChange
  
  #Name the vectors
  names(geneList) <- as.character(changing_df$feature)
  
  #sort the vectors in decreasing order.
  geneList <- sort(geneList, decreasing = TRUE)
  
  #run the ORA analysis
  
  output <- enricher(gene=names(geneList),
                     pvalueCutoff = 0.05, 
                     pAdjustMethod = "BH",
                     minGSSize = 10,
                     maxGSSize = 500,
                     TERM2GENE = term2gene_universe,
                     TERM2NAME = term2name)
  
  
  return(output)
  
}

#input_lists can be a list of dfs like tibble::lst(S_epi_Seb_vs_all_log_phase_GSEA,
#S_epi_Seb_vs_all_OS_GSEA, S_epi_Seb_vs_all_stat_phase_GSEA)
plot_facetted_NES <- function(input_list,
                              shared_gene_sets,
                              plot_title){
  
  #subset for shared gene sets between all the comparisons
  
  combined_df_subset <- lapply(input_list, function(input_df){
    
    df_subset <- input_df %>% dplyr::filter(ID %in% shared_gene_sets)
    
    df_subset$ID_and_Desc <- paste0(df_subset$ID," ", df_subset$Description)
    
    return(df_subset)
    
  }) %>% do.call("rbind",.)
  
  
  ggplot(combined_df_subset, aes(x=ID_and_Desc, y= NES, fill = comparison)) +
    geom_col(colour="black",position="dodge") +
    facet_wrap(~set_annotation, nrow=2, ncol=2, scales='free') +
    scale_fill_manual(values = c23) +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1)) +
    coord_flip() +
    labs(x="gene sets/pathways", y="Normalized Enrichment Score") + ggtitle(plot_title) + theme_classic()
  
  
}

