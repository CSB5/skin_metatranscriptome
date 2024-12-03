
#species_transcriptome_name is a string that identifies a species of interest for their transcriptomes to be analyzed 
#secretory_prediction_df is the processed output of signalP v 6
#aldex_raw_clr_values[["Sc_bracken_clr"]] etc can be the aldex_clr_input 
#site_libs is a character vector 

get_mtx_vs_mgx_correlation <- function(mtx_count_matrix, aldex_clr_input, secretory_prediction_df, site_libs, 
                                       species_transcriptome_name, species_transcriptome_read_threshold=200000, 
                                       responding_species, 
                                       cor_test_method="spearman",
                                       manual_vst_nsub = FALSE){
  
  ######
  #raw counts, species specific and site specific
  count_matrix <- mtx_count_matrix[str_detect(rownames(mtx_count_matrix), pattern=species_transcriptome_name), site_libs]
  
  #pre-filtering
  #keep if at least 10 reads in 10 samples across the metatranscriptomes?
  # rowSums(count_matrix >= 10) >= 10
  
  keep_transcript <- (rowMedians(count_matrix) >= 10)
  keep_sample <- (colSums(count_matrix)>=species_transcriptome_read_threshold)
  
  count_matrix <- count_matrix[keep_transcript, keep_sample]
  
  
  #vst() can be called on a matrix of counts and it will also call estimateSizeFactors within the function if you have not already supplied size factors: https://www.biostars.org/p/434003/
  
  if (manual_vst_nsub == FALSE){
    count_matrix_vst <- vst(count_matrix, blind=FALSE)
  } else if (manual_vst_nsub == TRUE){
    count_matrix_vst <- vst(count_matrix, blind=FALSE, nsub=sum( rowMedians(count_matrix) >= 10 ))
  }
  
  #subset for proteins in the secretory pathway
  
  count_matrix_vst_secretory <- count_matrix_vst[rownames(count_matrix_vst) %in% 
                                                   secretory_prediction_df$pangene, ]
  
  pangenes_to_test <- rownames(count_matrix_vst_secretory)
  
  
  
  #prepare the species abundance inputs (clr transformed metagenomic abundances with matching mtx library ID label)
  #load all Monte Carlo instances from aldex2 and do correlation tests in a way similar to the aldexcorr function,
  #where correlation tests are done for each MC instance and an average correlation statistic and p-value are reported
  #See ?propr::aldex.cor
  #https://rdrr.io/github/tpq/propr/src/R/7-aldex2propr.R
  
  aldex_clr_MC <- getMonteCarloInstances(aldex_clr_input)
  
  #extract the species specific clr MC instances
  aldex_clr_MC_species_specific <- lapply(names(aldex_clr_MC), function(LIBID){
    
    mat_out <- aldex_clr_MC[[LIBID]][responding_species,]
    
  })
  
  names(aldex_clr_MC_species_specific) <- names(aldex_clr_MC)
  
  aldex_clr_MC_species_specific_mat <- do.call("rbind", aldex_clr_MC_species_specific)
  
  
  #Make sure that the libraries in the clr input match the filtered libraries from the mtx data
  matching_mtx_ids <- colnames(count_matrix_vst_secretory)
  
  matching_mgx_ids <- mtx_mgx_stats_chosen %>% dplyr::filter(mtx_LIBID %in% matching_mtx_ids)
  matching_mgx_ids <- matching_mgx_ids[match(matching_mtx_ids, matching_mgx_ids$mtx_LIBID), ] %>% pull(mgx_LIBID)
  
  #subset the aldex clr matrix in the order of the matching mgx_ids
  #each column is a monte carlo instance
  aldex_clr_MC_species_specific_mat_filt <- aldex_clr_MC_species_specific_mat[matching_mgx_ids,]
  
  
  if (cor_test_method=="pearson"){
    
    cor_test_coef <- "cor"
    
  } else if (cor_test_method=="spearman"){
    
    cor_test_coef <- "rho"
    
  }
  
  
  #initialize vectors to hold test results
  MC_corr_vec <- vector(mode="numeric", length=1000)
  MC_pval_vec <- vector(mode="numeric", length=1000)
  #Initialize vectors to hold final outputs
  corr_p_vals <- vector(mode="numeric", length=length(pangenes_to_test))
  corr_coeff <- vector(mode="numeric", length=length(pangenes_to_test))
  
  #For each pangene (in the secretory pathway...)
  #perform the correlation tests for each MC instance (1000 of them)
  for (i in 1:length(pangenes_to_test)){
    
    #vector of vst counts for a given pangene
    mtx_vst_vector <- count_matrix_vst_secretory[pangenes_to_test[i], matching_mtx_ids]
    
    for (j in 1:1000){
      MC_corr_results <- cor.test(aldex_clr_MC_species_specific_mat_filt[,j], mtx_vst_vector,  
                                  method=cor_test_method, exact=FALSE )
      
      MC_pval_vec[j] <- MC_corr_results$p.value
      MC_corr_vec[j] <- MC_corr_results$estimate[[cor_test_coef]]
      
    }
    
    corr_p_vals[i] <- mean(MC_pval_vec)
    corr_coeff[i] <- mean(MC_corr_vec)
    
  } 
  
  corr_df <- data.frame(pangene=pangenes_to_test,
                        corr_p_vals,
                        corr_coeff,
                        species_tested_for=responding_species)
  
  corr_df$p_adj <-  p.adjust(corr_df$corr_p_vals, method="fdr")
  
  
  output <- tibble::lst(count_matrix_vst, count_matrix_vst_secretory, corr_df)
  
  return(output)
  
}

#
transcript_vs_clr_scatterplot <- function(vst_mat, mean_clr, transcript, species){
  
  
  vst_df <- t(vst_mat) %>% as.data.frame()
  
  vst_df$mtx_LIBID <- rownames(vst_df)
  
  species_mean_clr <- mean_clr %>% dplyr::select(species,"mtx_LIBID")
  
  transcript_and_clr_df <- merge(vst_df,species_mean_clr, by ="mtx_LIBID")
  
  ggscatter(transcript_and_clr_df, x = species, 
            y = transcript, add = "reg.line") + stat_cor(method="spearman") + theme_classic()
  
  
}
