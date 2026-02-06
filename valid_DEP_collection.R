# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Explicitly ensure dplyr functions are prioritized
select <- dplyr::select
filter <- dplyr::filter

# ============================================================================
# CONFIGURATION: Specify comparisons of interest and their labels
# ============================================================================
comparisons_of_interest <- c(
  "diabetic_empty_42-nondiabetic_empty_42" = "Diabetic Untreated vs Non-diabetic Untreated",
  "diabetic_PCL_42-nondiabetic_PCL_42" = "Diabetic PCL Scaffold vs Non-diabetic PCL Scaffold"
)

############################################
# Read data files
# ---
output_dir <- "valid_DEPs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read protein abundance data
abundance_file_AI_unscaled <- "normalization_results_AI/RobNorm_normalized_data.csv"
abundance_file_AS_unscaled <- "normalization_results_AS/RobNorm_normalized_data.csv"
abundance_file_AI <- "normalization_results_AI/scaled/RobNorm_scaled_normalized_data.csv"
abundance_file_AS <- "normalization_results_AS/scaled/RobNorm_scaled_normalized_data.csv"
metadata_file_AI <- "normalization_results_AI/scaled/meta_data.csv"
metadata_file_AS <- "normalization_results_AS/scaled/meta_data.csv"
ai_data_unscaled <- read.csv(abundance_file_AI_unscaled, stringsAsFactors = FALSE, check.names = FALSE)
as_data_unscaled <- read.csv(abundance_file_AS_unscaled, stringsAsFactors = FALSE, check.names = FALSE)
ai_data <- read.csv(abundance_file_AI, stringsAsFactors = FALSE, check.names = FALSE)
as_data <- read.csv(abundance_file_AS, stringsAsFactors = FALSE, check.names = FALSE)
meta_ai <- read.csv(metadata_file_AI, stringsAsFactors = FALSE)
meta_as <- read.csv(metadata_file_AS, stringsAsFactors = FALSE)

# Read meta-analysis proteins
meta_proteins <- read.csv("input/bone_caps_meta_analysis.csv", stringsAsFactors = FALSE)
meta_genes <- unique(toupper(meta_proteins$Gene))

# Read DE results
de_results_AI <- read.csv("de_analysis_results_AI_RobNorm_scaled/de_results_raw_with_orthologs.csv", stringsAsFactors = FALSE)
de_results_AS <- read.csv("de_analysis_results_AS_RobNorm_scaled/de_results_raw_with_orthologs.csv", stringsAsFactors = FALSE)
# Filter DE results to only include specified comparisons
de_results_AI <- de_results_AI[de_results_AI$Comparison %in% names(comparisons_of_interest), ]
de_results_AS <- de_results_AS[de_results_AS$Comparison %in% names(comparisons_of_interest), ]
############################################

############################################
# Create Sample_ID to match AI and AS (same animal, condition, scaffold, intervention, timepoint, tissue) 
# ---
meta_ai$Sample_ID <- paste(meta_ai$Animal, meta_ai$Condition, meta_ai$Scaffold, 
                           meta_ai$Intervention, meta_ai$Timepoint, meta_ai$Tissue, sep="_")
meta_as$Sample_ID <- paste(meta_as$Animal, meta_as$Condition, meta_as$Scaffold, 
                           meta_as$Intervention, meta_as$Timepoint, meta_as$Tissue, sep="_")

# Match AI and AS samples
matched_samples <- inner_join(
  meta_ai %>% select(Sample_ID, Column_AI = Column, 
                     Condition, Scaffold, Intervention, Timepoint, Tissue, Animal),
  meta_as %>% select(Sample_ID, Column_AS = Column),
  by = "Sample_ID"
)

cat("Matched", nrow(matched_samples), "AI/AS sample pairs\n\n")
############################################

############################################
# HELPER FUNCTIONS
# Function to parse group names (e.g., "diabetic_PCL_42") from the comparison name
parse_group <- function(group_str) {
  parts <- strsplit(group_str, "_")[[1]]
  condition <- parts[1] # Get the very first list element
  timepoint <- as.integer(parts[length(parts)])
  scaffold <- parts[(length(parts)-1)]
  return(list(condition = condition, scaffold = scaffold, timepoint = timepoint))
}

get_fraction_abundances <- function(protein, samples, ai_data, as_data) {
  # Get AI abundances
  ai_abundances <- numeric()
  for (i in 1:nrow(samples)) { # The number of samples
    ai_col <- samples$Column_AI[i] # Get column name for the i-th sample
    # Select ai abundance value for this protein
    ai_value <- as.numeric(ai_data[ai_data$Protein.IDs == protein, ai_col])
    if (length(ai_value) > 0 && !is.na(ai_value)) {
      # There is a value measured
      ai_abundances <- c(ai_abundances, 2^ai_value)
    }
  }
  
  # Get AS abundances
  as_abundances <- numeric()
  for (i in 1:nrow(samples)) {
    as_col <- samples$Column_AS[i]
    as_value <- as.numeric(as_data[as_data$Protein.IDs == protein, as_col])
    if (length(as_value) > 0 && !is.na(as_value)) {
      as_abundances <- c(as_abundances, 2^as_value)
    }
  }
  return (list(ai_abundances=ai_abundances, 
               as_abundances=as_abundances))
}

# Function to calculate AI/AS ratios for a group of samples and get fraction-exclusive proteins
calculate_group_ratios_and_exclusive_proteins <- function(samples, ai_data, as_data, 
                                                          ai_data_unscaled, as_data_unscaled) {
  protein_ratios <- data.frame()
  ai_exclusive_proteins <- character()
  as_exclusive_proteins <- character()
  
  # Get all unique proteins from both fractions
  all_proteins <- unique(c(ai_data$Protein.IDs, as_data$Protein.IDs))
  
  for (protein_id in all_proteins) {
    abundances_unscaled = get_fraction_abundances(protein_id, samples, ai_data_unscaled, as_data_unscaled)
    ai_abundances_unscaled = abundances_unscaled$ai_abundances
    as_abundances_unscaled = abundances_unscaled$as_abundances
    
    abundances = get_fraction_abundances(protein_id, samples, ai_data, as_data)
    ai_abundances = abundances$ai_abundances
    as_abundances = abundances$as_abundances
    
    if (length(as_abundances) > 0 && length(ai_abundances) > 0) {
      # PROTEIN IS MEASURED IN AT LEAST ONE SAMPLE IN BOTH FRACTIONS - calculate ratio
      mean_ai <- mean(ai_abundances, na.rm = TRUE)
      mean_as <- mean(as_abundances, na.rm = TRUE)
      # THIS LINE NAIVELY ASSUMES THAT THE SUM OF BOTH FRACTIONS IS EQUAL TO THE TOTAL INTENSITY
      # THIS IS TECHNICALLY NOT CORRECT
      total <- mean_ai + mean_as
      
      protein_ratios <- rbind(protein_ratios, data.frame(
        Protein.IDs = protein_id,
        Mean_AI = mean_ai,
        Mean_AS = mean_as,
        Pct_AI = (mean_ai / total) * 100,
        Fraction_Status = "Both",
        stringsAsFactors = FALSE
      ))
    } else if(length(as_abundances_unscaled) > 0 && !(length(ai_abundances_unscaled) > 0)) {
        # The protein is exclusively measured in AS 
        as_exclusive_proteins <- c(as_exclusive_proteins, protein_id)
        
        protein_ratios <- rbind(protein_ratios, data.frame(
          Protein.IDs = protein_id,
          Mean_AI = NA,
          Mean_AS = NA,
          Pct_AI = 0,
          Fraction_Status = "AS_Exclusive",
          stringsAsFactors = FALSE
        ))
    } else if(length(ai_abundances_unscaled) > 0 && !(length(as_abundances_unscaled) > 0)) {
        ai_exclusive_proteins <- c(ai_exclusive_proteins, protein_id)
        
        protein_ratios <- rbind(protein_ratios, data.frame(
          Protein.IDs = protein_id,
          Mean_AI = NA,
          Mean_AS = NA,
          Pct_AI = 100,
          Fraction_Status = "AI_Exclusive",
          stringsAsFactors = FALSE
        ))
    }
  }
  
  return(list(
    ratios = protein_ratios,
    ai_exclusive = ai_exclusive_proteins,
    as_exclusive = as_exclusive_proteins
  ))
}

############################################
# Calculate fold changes across fractions for every protein in every sample (pair)
# ---
all_fc <- list()
cat("Processing", nrow(matched_samples), "sample pairs...\n")

for (i in 1:nrow(matched_samples)) {
  ai_col <- matched_samples$Column_AI[i]
  as_col <- matched_samples$Column_AS[i]
  
  # Get AI data for this sample
  ai_sample <- data.frame(
    Protein.IDs = ai_data$Protein.IDs,
    AI_abundance = as.numeric(ai_data[[ai_col]]),
    stringsAsFactors = FALSE
  )
  
  # Get AS data for this sample
  as_sample <- data.frame(
    Protein.IDs = as_data$Protein.IDs,
    AS_abundance = as.numeric(as_data[[as_col]]),
    stringsAsFactors = FALSE
  )
  
  # Merge by Protein.IDs to match correctly
  # Only keeps proteins detected in BOTH fractions
  merged <- inner_join(ai_sample, as_sample, by = "Protein.IDs")
  
  # Calculate fold change: larger / smaller (always >= 1)
  merged$AI_linear <- 2^merged$AI_abundance
  merged$AS_linear <- 2^merged$AS_abundance
  merged$Fold_Change <- pmax(merged$AI_linear, merged$AS_linear, na.rm = FALSE) / 
    pmin(merged$AI_linear, merged$AS_linear, na.rm = FALSE)
  merged$Dominant_Fraction <- ifelse(merged$AI_linear > merged$AS_linear, "AI", "AS")
  
  # Store result
  all_fc[[i]] <- data.frame(
    Protein.IDs = merged$Protein.IDs,
    Sample_ID = matched_samples$Sample_ID[i],
    Fold_Change = merged$Fold_Change,
    AI_linear = merged$AI_linear, 
    AS_linear = merged$AS_linear,
    Dominant_Fraction = merged$Dominant_Fraction,
    stringsAsFactors = FALSE
  )
}

# Combine all
# Each row = one protein in one sample pair (AI + AS)
# Same protein (e.g., P12345) appears in multiple rows, one for each sample
complete_fcs_proteinwise <- bind_rows(all_fc) %>%
  filter(!is.na(Fold_Change) & is.finite(Fold_Change))

cat("Total data points:", nrow(complete_fcs_proteinwise), "\n\n")
############################################

############################################
# Calculate percentiles and plot
# ---
median_fc <- median(complete_fcs_proteinwise$Fold_Change)
p70 <- quantile(complete_fcs_proteinwise$Fold_Change, 0.70)
p80 <- quantile(complete_fcs_proteinwise$Fold_Change, 0.80)
p90 <- quantile(complete_fcs_proteinwise$Fold_Change, 0.90)

cat("Percentiles:\n")
cat("  Median:", round(median_fc, 2), "\n")
cat("  70th:", round(p70, 2), "\n")
cat("  80th:", round(p80, 2), "\n")
cat("  90th:", round(p90, 2), "\n")
cat("  Max:", round(max(complete_fcs_proteinwise$Fold_Change), 2), "\n\n")

percentile_data <- data.frame(
  label = c("Median", "70th", "80th", "90th"),
  value = c(median_fc, p70, p80, p90),
  text = sprintf("%.2f", c(median_fc, p70, p80, p90))
)

plot <- ggplot(complete_fcs_proteinwise, aes(x = Fold_Change, fill = Dominant_Fraction)) +
  geom_histogram(bins = 100, alpha = 0.7, color = "white", position = "dodge") +
  geom_vline(data = percentile_data, aes(xintercept = value, color = label), 
             linewidth = 1.2, linetype = "dashed") +
  scale_fill_manual(values = c("AI" = "#2E86AB", "AS" = "#E63946")) +
  scale_color_manual(
    values = c("Median" = "#FFB703", "70th" = "#06A77D", "80th" = "#D62828", "90th" = "#7209B7"),
    breaks = c("Median", "70th", "80th", "90th"),
    labels = setNames(
      paste0(percentile_data$label, ": ", percentile_data$text),
      percentile_data$label
    )
  ) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 500, 1000)) +
  labs(
    title = "AI vs AS Fold Change Distribution",
    subtitle = paste0("Every protein in every sample (n = ", nrow(complete_fcs_proteinwise), ")"),
    x = "Fold Change (Larger / Smaller)",
    y = "Count",
    color = "Percentiles",
    fill = "Dominant Fraction"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "FC_distribution_separated.png"), plot, width = 12, height = 7, dpi = 300)

# Find proteins with extreme fold changes
extreme_fc <- complete_fcs_proteinwise %>%
  filter(Fold_Change > 100) %>%  # or whatever threshold seems reasonable
  arrange(desc(Fold_Change))

# Check the actual intensity values
cat("\nTop 10 extreme fold changes:\n")
for(i in 1:min(10, nrow(extreme_fc))) {
  cat("\nProtein:", extreme_fc$Protein.IDs[i], "\n")
  cat("  Sample:", extreme_fc$Sample_ID[i], "\n")
  cat("  AI intensity:", extreme_fc$AI_linear[i], "\n")
  cat("  AS intensity:", extreme_fc$AS_linear[i], "\n")
  cat("  Fold Change:", extreme_fc$Fold_Change[i], "\n")
  cat("  Dominant:", extreme_fc$Dominant_Fraction[i], "\n")
}
###################################################################

###################################################################
# ============================================================================
# MAIN THRESHOLD TESTING LOOPS
# ============================================================================

fc_thresholds <- c(2.07, 2.49, 3.38)
stability_thresholds <- c(3, 5, 7)

# Store all results across all threshold combinations
all_threshold_results <- list()

for (fc_thresh in fc_thresholds) {
  for (stab_thresh in stability_thresholds) {
    
    cat("\n========================================\n")
    cat("Testing FC threshold:", fc_thresh, "| Stability:", stab_thresh, "%\n")
    cat("========================================\n")
    
    # Create subdirectory for this combination
    thresh_dir <- file.path(output_dir, paste0("FC", fc_thresh, "_Stab", stab_thresh))
    dir.create(thresh_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Initialize validated deps list for THIS threshold combination
    all_validated_deps <- list()
    
    # ========================================================================
    # Loop through comparisons
    # ========================================================================
    for(comp in names(comparisons_of_interest)) {
      cat("\nProcessing:", comparisons_of_interest[comp], "\n")
      
      # Parse comparison to get the two groups
      comp_parts <- strsplit(comp, "-")[[1]]
      group1 <- parse_group(comp_parts[1])
      group2 <- parse_group(comp_parts[2])
      
      # Get samples for each group
      group1_samples <- matched_samples %>%
        filter(Condition == group1$condition, 
               Scaffold == group1$scaffold, 
               Timepoint == group1$timepoint)
      
      group2_samples <- matched_samples %>%
        filter(Condition == group2$condition, 
               Scaffold == group2$scaffold, 
               Timepoint == group2$timepoint)
      
      comparison_samples <- rbind(group1_samples, group2_samples)
      
      # ======================================================================
      # Calculate per-protein fold changes FOR THIS COMPARISON
      # For every protein we take the mean of the intensities of all samples in each fraction 
      # ======================================================================
      comparison_fc <- data.frame()
      
      all_proteins <- unique(c(ai_data$Protein.IDs, as_data$Protein.IDs))
      
      for (protein_id in all_proteins) {
        # Get abundances using scaled data
        abundances <- get_fraction_abundances(protein_id, comparison_samples, ai_data, as_data)
        ai_abundances <- abundances$ai_abundances
        as_abundances <- abundances$as_abundances
        
        # Calculate fold change if protein is in both fractions
        if (length(ai_abundances) > 0 && length(as_abundances) > 0) {
          mean_ai <- mean(ai_abundances, na.rm = TRUE)
          mean_as <- mean(as_abundances, na.rm = TRUE)
          
          comparison_fc <- rbind(comparison_fc, data.frame(
            Protein.IDs = protein_id,
            mean_AI = mean_ai,
            mean_AS = mean_as,
            FC_between_fractions = pmax(mean_ai, mean_as) / pmin(mean_ai, mean_as),
            stringsAsFactors = FALSE
          ))
        }
      }      
      # ======================================================================
      # Calculate ratios for both groups in this comparison
      # ======================================================================
      group1_ratios <- calculate_group_ratios_and_exclusive_proteins(group1_samples, ai_data, as_data,
                                                                     ai_data_unscaled, as_data_unscaled)
      group2_ratios <- calculate_group_ratios_and_exclusive_proteins(group2_samples, ai_data, as_data,
                                                                     ai_data_unscaled, as_data_unscaled)
      
      # Compare ratios
      ratio_comparison <- inner_join(
        group1_ratios$ratios %>% select(Protein.IDs, Pct_AI_Group1 = Pct_AI),
        group2_ratios$ratios %>% select(Protein.IDs, Pct_AI_Group2 = Pct_AI),
        by = "Protein.IDs"
      )
      
      ratio_comparison$Pct_AI_Diff <- abs(ratio_comparison$Pct_AI_Group1 - ratio_comparison$Pct_AI_Group2)
      
      # ======================================================================
      # Get First-level and Second-level validated DEPs
      # ======================================================================
      deps_ai <- de_results_AI[de_results_AI$Comparison == comp & de_results_AI$Change != "No Change", ]
      deps_as <- de_results_AS[de_results_AS$Comparison == comp & de_results_AS$Change != "No Change", ]
      all_deps <- unique(c(deps_ai$Protein.IDs, deps_as$Protein.IDs))
      
      validated_deps <- data.frame()
      
      for(protein in all_deps) {
        # Check if protein is fraction-exclusive
        is_ai_exclusive <- protein %in% group1_ratios$ai_exclusive && 
          protein %in% group2_ratios$ai_exclusive
        is_as_exclusive <- protein %in% group1_ratios$as_exclusive && 
          protein %in% group2_ratios$as_exclusive
        
        # Check if protein is DEP in each fraction
        in_ai <- protein %in% deps_ai$Protein.IDs
        in_as <- protein %in% deps_as$Protein.IDs
        
        # Get direction of change
        direction_ai <- ifelse(in_ai, deps_ai[deps_ai$Protein.IDs == protein, "Change"], NA)
        direction_as <- ifelse(in_as, deps_as[deps_as$Protein.IDs == protein, "Change"], NA)
        
        # Get fold change between fractions
        fc_info <- comparison_fc[comparison_fc$Protein.IDs == protein, ]
        fc_between <- ifelse(nrow(fc_info) > 0, fc_info$FC_between_fractions[1], NA)
        # fc_info$FC_between_fractions[1] to make sure we get NA if there is no value
        
        # Validation logic
        is_valid <- FALSE
        validation_reason <- ""
        
        if(is_ai_exclusive && in_ai) {
          is_valid <- TRUE
          validation_reason <- "Exclusive to AI fraction"
        } else if(is_as_exclusive && in_as) {
          is_valid <- TRUE
          validation_reason <- "Exclusive to AS fraction"
        } else if(in_ai & in_as) {
          if(!is.na(direction_ai) & !is.na(direction_as) & direction_ai == direction_as) {
            is_valid <- TRUE
            validation_reason <- "Significant in both fractions (same direction)"
          }
        } else if((in_ai | in_as) & !(in_ai & in_as) & !is.na(fc_between) & fc_between >= fc_thresh) {
          is_valid <- TRUE
          validation_reason <- paste0("Dominant in ", ifelse(in_ai, "AI", "AS"), 
                                      " fraction (FC=", round(fc_between, 2), ")")
        }
        
        if(is_valid) {
          gene_name <- ifelse(in_ai, 
                              deps_ai[deps_ai$Protein.IDs == protein, "Gene.Names"],
                              deps_as[deps_as$Protein.IDs == protein, "Gene.Names"])
          
          ortholog_name <- ifelse(in_ai, 
                                  deps_ai[deps_ai$Protein.IDs == protein, "Orthologs"],
                                  deps_as[deps_as$Protein.IDs == protein, "Orthologs"])
          
          validated_deps <- rbind(validated_deps, data.frame(
            Protein.IDs = protein,
            Gene.Names = gene_name,
            Orthologs = ortholog_name,
            Direction = ifelse(!is.na(direction_ai), direction_ai, direction_as),
            DEP_in_AI = in_ai,
            DEP_in_AS = in_as,
            FC_between_fractions = fc_between,
            Validation = validation_reason,
            stringsAsFactors = FALSE
          ))
        }
      }
      
      # ======================================================================
      # Filter second-level validated DEPs by ratio consistency
      # ======================================================================
      if(nrow(validated_deps) > 0) {
        validated_deps <- validated_deps %>%
          left_join(ratio_comparison %>% 
                      select(Protein.IDs, Pct_AI_Group1, Pct_AI_Group2, Pct_AI_Diff),
                    by = "Protein.IDs")
        
        n_before <- nrow(validated_deps)
        
        validated_deps <- validated_deps %>%
          filter(
            grepl("both fractions|Exclusive", Validation) |
              (grepl("Dominant", Validation) & Pct_AI_Diff <= stab_thresh)
          )
        
        n_after <- nrow(validated_deps)
        cat("  Filtered by ratio consistency (≤", stab_thresh, "%):", n_before, "->", n_after, "\n")
      }
      
      # Store with comparison info
      if(nrow(validated_deps) > 0) {
        validated_deps$Comparison <- comp
        validated_deps$Comparison_Label <- comparisons_of_interest[comp]
      }
      all_validated_deps[[comp]] <- validated_deps
      
      cat("  Total unique DEPs (AI+AS):", length(all_deps), "\n")
      cat("  Validated DEPs:", nrow(validated_deps), "\n")
    }
    
    # ========================================================================
    # Save results for this threshold combination
    # ========================================================================
    validated_deps_df <- bind_rows(all_validated_deps)
    write.csv(validated_deps_df, 
              file.path(thresh_dir, "validated_DEPs.csv"), 
              row.names = FALSE)
    cat("\nSaved validated DEPs to:", file.path(thresh_dir, "validated_DEPs.csv"), "\n")
    
    # ========================================================================
    # Generate summary statistics and plots for this threshold combination
    # ========================================================================
    summary_stats <- data.frame()
    
    for(comp in names(comparisons_of_interest)) {
      # Get orthologs for tested proteins in this comparison
      tested_orthologs_ai <- de_results_AI[de_results_AI$Comparison == comp, "Orthologs"]
      tested_orthologs_as <- de_results_AS[de_results_AS$Comparison == comp, "Orthologs"]
      all_tested_orthologs <- unique(c(tested_orthologs_ai, tested_orthologs_as))
      tested_genes <- toupper(sapply(strsplit(as.character(all_tested_orthologs), ";"), `[`, 1))
      tested_genes <- unique(tested_genes)
      bone_caps_tested <- sum(tested_genes %in% meta_genes, na.rm = TRUE)
      
      deps_ai_count <- sum(de_results_AI$Comparison == comp & de_results_AI$Change != "No Change")
      deps_as_count <- sum(de_results_AS$Comparison == comp & de_results_AS$Change != "No Change")
      total_unique_deps <- length(unique(c(
        de_results_AI[de_results_AI$Comparison == comp & de_results_AI$Change != "No Change", "Protein.IDs"],
        de_results_AS[de_results_AS$Comparison == comp & de_results_AS$Change != "No Change", "Protein.IDs"]
      )))
      
      all_dep_orthologs <- unique(c(
        de_results_AI[de_results_AI$Comparison == comp & de_results_AI$Change != "No Change", "Orthologs"],
        de_results_AS[de_results_AS$Comparison == comp & de_results_AS$Change != "No Change", "Orthologs"]
      ))
      all_dep_genes <- toupper(sapply(strsplit(as.character(all_dep_orthologs), ";"), `[`, 1))
      all_dep_genes <- unique(all_dep_genes)
      overlap_total <- sum(all_dep_genes %in% meta_genes, na.rm = TRUE)
      
      validated_count <- ifelse(comp %in% names(all_validated_deps), 
                                nrow(all_validated_deps[[comp]]), 0)
      
      if(validated_count > 0) {
        validated <- all_validated_deps[[comp]]
        n_both <- sum(validated$Validation == "Significant in both fractions (same direction)")
        n_dominant <- sum(grepl("Dominant in", validated$Validation))
        n_exclusive <- sum(grepl("Exclusive to", validated$Validation))
        
        first_level_proteins <- validated[validated$Validation == "Significant in both fractions (same direction)" | 
                                            grepl("Exclusive to", validated$Validation), ]
        first_level_genes <- toupper(sapply(strsplit(as.character(first_level_proteins$Orthologs), ";"), `[`, 1))
        first_level_genes <- unique(first_level_genes)
        overlap_first_level <- sum(first_level_genes %in% meta_genes, na.rm = TRUE)
        
        second_level_proteins <- validated[grepl("Dominant in", validated$Validation), ]
        second_level_genes <- toupper(sapply(strsplit(as.character(second_level_proteins$Orthologs), ";"), `[`, 1))
        second_level_genes <- unique(second_level_genes)
        overlap_second_level <- sum(second_level_genes %in% meta_genes, na.rm = TRUE)
      } else {
        n_both <- 0
        n_dominant <- 0
        n_exclusive <- 0
        overlap_first_level <- 0
        overlap_second_level <- 0
      }
      
      summary_stats <- rbind(summary_stats, data.frame(
        Comparison = comp,
        Comparison_Label = comparisons_of_interest[comp],
        DEPs_AI = deps_ai_count,
        DEPs_AS = deps_as_count,
        Bone_Caps_Tested = bone_caps_tested,
        Total_Unique_DEPs = total_unique_deps,
        total_unique_orthologs_count = length(all_dep_genes),
        Overlap_Total = overlap_total,
        Validated_Both_Fractions = n_both,
        Validated_Exclusive = n_exclusive,
        Validated_Dominant = n_dominant,
        Validated_Total = validated_count,
        Overlap_First_Level = overlap_first_level,
        Overlap_Second_Level = overlap_second_level,
        FC_Threshold = fc_thresh,
        Stability_Threshold = stab_thresh,
        stringsAsFactors = FALSE
      ))
    }
    
    # Save summary stats
    write.csv(summary_stats, 
              file.path(thresh_dir, "summary_statistics.csv"), 
              row.names = FALSE)
    
    # Create plot
    plot_data <- summary_stats %>%
      mutate(
        First_Level_Total = Validated_Both_Fractions + Validated_Exclusive,
        Second_Level_Total = Validated_Dominant
      ) %>%
      select(Comparison_Label, Total_Unique_DEPs, Overlap_Total,
             First_Level_Total, Overlap_First_Level,
             Second_Level_Total, Overlap_Second_Level) %>%
      pivot_longer(cols = -Comparison_Label,
                   names_to = "Category",
                   values_to = "Count")
    
    plot_data$Category <- factor(plot_data$Category, 
                                 levels = c("Total_Unique_DEPs", "Overlap_Total",
                                            "First_Level_Total", "Overlap_First_Level",
                                            "Second_Level_Total", "Overlap_Second_Level"))
    
    p <- ggplot(plot_data, aes(x = Comparison_Label, y = Count, fill = Category)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_text(aes(label = Count),
                position = position_dodge(width = 0.9),
                vjust = -0.5, size = 3, fontface = "bold") +
      scale_fill_manual(
        values = c("Total_Unique_DEPs" = "#E63946",
                   "Overlap_Total" = "#F5A3AD",
                   "First_Level_Total" = "#2E86AB",
                   "Overlap_First_Level" = "#A8DADC",
                   "Second_Level_Total" = "#F77F00",
                   "Overlap_Second_Level" = "#FCBF49"),
        labels = c("Total Unique DEPs (AI + AS)",
                   "Overlap Total with Meta-analysis",
                   "First-Level Validated DEPs",
                   "Overlap First-Level with Meta-analysis",
                   "Second-Level Validated DEPs",
                   "Overlap Second-Level with Meta-analysis")
      ) +
      theme_minimal() +
      labs(
        title = "DEPs with Two-Level Validation Strategy",
        subtitle = paste0("FC ≥ ", fc_thresh, " | Stability ≤ ", stab_thresh, "%"),
        y = "Number of Proteins",
        x = "",
        fill = ""
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        plot.margin = margin(20, 20, 20, 60)
      )
    
    ggsave(file.path(thresh_dir, "validated_deps_summary.png"), 
           plot = p, width = 14, height = 8, dpi = 300)
    
    # Store this threshold's results
    all_threshold_results[[paste0("FC", fc_thresh, "_Stab", stab_thresh)]] <- list(
      validated_deps = validated_deps_df,
      summary_stats = summary_stats
    )
    
    cat("\n=== Completed FC", fc_thresh, "Stab", stab_thresh, "===\n")
  }
}

cat("\n\n========================================")
cat("\nALL THRESHOLD COMBINATIONS COMPLETE")
cat("\n========================================\n")