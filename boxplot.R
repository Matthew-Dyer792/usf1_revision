library(TFregulomeR)

# tf's in question
tfs <- c("MM1_HSA_H1-hESC_USF1", "MM1_HSA_H1-hESC_USF2", "MM1_HSA_H1-hESC_ATF3")

exclusivePeak_result_output <- c()
commonPeak_result_output <- c()
exclusive_names <- c()
common_names <- c()

# a loop to generate the exclusive peaks and then common peaks
for (target_tf in tfs) {
  # target_tf is the tf in the x angle
  # remaining_tfs is a list of tfs without the current target_tf
  remaining_tfs <- tfs[!(tfs %in% target_tf)]
  
  # generates the peaks for your target_tf that do not include the peaks of the specified tfs 
  exclusivePeak_output <- exclusivePeaks(target_peak_id = target_tf,
                                         motif_only_for_target_peak = TRUE,
                                         excluded_peak_id = remaining_tfs,
                                         TFregulome_url = 'http://localhost/MethMotif/2021/')
  
  exclusivePeak_result <- exclusivePeakResult(exclusivePeaks = exclusivePeak_output,
                                              return_exclusive_peak_sites = TRUE)
  
  # generates a name for the list below
  exclusive_names <- c(exclusive_names, paste0("exclusivePeaks_", target_tf))
  # a list variable to store the results from this loop
  exclusivePeak_result_output <- c(exclusivePeak_result_output, exclusivePeak_result[[1]])
  
   
  # a loop to generate the common peaks for all possible combinations
  for (tf_y in remaining_tfs) {
    # common peaks function returns all shared peaks between target_tf (angle x) and tf_y (angle y)
    commonPeak_output <- commonPeaks(target_peak_id = target_tf,
                                     compared_peak_id = tf_y,
                                     TFregulome_url = 'http://localhost/MethMotif/2021/')
    
    commonPeak_result <- commonPeakResult(commonPeaks = commonPeak_output,
                                          return_common_peak_sites = TRUE)
    
    # generates a name for the list below
    common_names <- c(common_names, paste0("commonPeaks_", target_tf, "_", tf_y))
    # a list variable to store the results from this loop
    commonPeak_result_output <- c(commonPeak_result_output, commonPeak_result[[1]])
    
    # the tf that is not tf_x or tf_y
    missing_tf <- remaining_tfs[!(remaining_tfs %in% tf_y)]
    customPeaks <- commonPeak_result$common_peak_list[[1]]
    
    # generates the peaks for your target_tf that do not include the peaks of the specified tfs 
    exclusivePeak_output <- exclusivePeaks(user_target_peak_list = list(customPeaks),
                                           user_target_peak_id = c(target_tf),
                                           motif_only_for_target_peak = TRUE,
                                           excluded_peak_id = missing_tf,
                                           TFregulome_url = 'http://localhost/MethMotif/2021/')
    
    exclusivePeak_result <- exclusivePeakResult(exclusivePeaks = exclusivePeak_output,
                                                return_exclusive_peak_sites = TRUE)
    
    # generates a name for the list below
    exclusive_names <- c(exclusive_names, paste0("exclusivePeaks_", target_tf, "_", tf_y))
    # a list variable to store the results from this loop
    exclusivePeak_result_output <- c(exclusivePeak_result_output, exclusivePeak_result[[1]])
  }
}

# applying the names generated to the results lists 
names(exclusivePeak_result_output) <- exclusive_names
names(commonPeak_result_output) <- common_names

# example of the plotting process
boxplot(exclusivePeak_result_output$`exclusivePeaks_MM1_HSA_H1-hESC_USF1`$tag_fold_change,
        exclusivePeak_result_output$`exclusivePeaks_MM1_HSA_H1-hESC_USF1_MM1_HSA_H1-hESC_USF2`$tag_fold_change,
        commonPeak_result_output$`commonPeaks_MM1_HSA_H1-hESC_USF1_MM1_HSA_H1-hESC_USF2`$tag_fold_change,
        exclusivePeak_result_output$`exclusivePeaks_MM1_HSA_H1-hESC_USF1_MM1_HSA_H1-hESC_ATF3`$tag_fold_change,
        commonPeak_result_output$`commonPeaks_MM1_HSA_H1-hESC_USF1_MM1_HSA_H1-hESC_ATF3`$tag_fold_change,
        main = "Multiple boxplots for comparision with USF1 as TF x",
        names = c("only USF1", "USF1 x USF2 - ATF3", "USF1 x USF2", "USF1 x ATF3 - USF2", "USF1 x ATF3"),
        ylab = "TAG Fold Change"
        # las = 2,
        # col = c("orange","red"),
        # border = "brown",
        # horizontal = TRUE,
        # notch = TRUE
)

boxplot(exclusivePeak_result_output$`exclusivePeaks_MM1_HSA_H1-hESC_USF2`$tag_fold_change,
        exclusivePeak_result_output$`exclusivePeaks_MM1_HSA_H1-hESC_USF2_MM1_HSA_H1-hESC_USF1`$tag_fold_change,
        commonPeak_result_output$`commonPeaks_MM1_HSA_H1-hESC_USF2_MM1_HSA_H1-hESC_USF1`$tag_fold_change,
        exclusivePeak_result_output$`exclusivePeaks_MM1_HSA_H1-hESC_USF2_MM1_HSA_H1-hESC_ATF3`$tag_fold_change,
        commonPeak_result_output$`commonPeaks_MM1_HSA_H1-hESC_USF2_MM1_HSA_H1-hESC_ATF3`$tag_fold_change,
        main = "Multiple boxplots for comparision with USF2 as TF x",
        names = c("only USF2", "USF2 x USF1 - ATF3", "USF2 x USF1", "USF2 x ATF3 - USF1", "USF2 x ATF3"),
        ylab = "TAG Fold Change"
)

boxplot(exclusivePeak_result_output$`exclusivePeaks_MM1_HSA_H1-hESC_ATF3`$tag_fold_change,
        exclusivePeak_result_output$`exclusivePeaks_MM1_HSA_H1-hESC_ATF3_MM1_HSA_H1-hESC_USF1`$tag_fold_change,
        commonPeak_result_output$`commonPeaks_MM1_HSA_H1-hESC_ATF3_MM1_HSA_H1-hESC_USF1`$tag_fold_change,
        exclusivePeak_result_output$`exclusivePeaks_MM1_HSA_H1-hESC_ATF3_MM1_HSA_H1-hESC_USF2`$tag_fold_change,
        commonPeak_result_output$`commonPeaks_MM1_HSA_H1-hESC_ATF3_MM1_HSA_H1-hESC_USF2`$tag_fold_change,
        main = "Multiple boxplots for comparision with ATF3 as TF x",
        names = c("only ATF3", "ATF3 x USF1 - USF2", "ATF3 x USF1", "ATF3 x USF2 - USF1", "ATF3 x USF2"),
        ylab = "TAG Fold Change"
)

save.image(file = "~/Downloads/revision_boxplot.RData")
