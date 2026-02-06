# Performs GO and KEGG enrichment analysis on validated DEPs across multiple 
# FC/stability threshold combinations. Outputs: (1) complete enrichment results, 
# (2) significant pathways only (p.adj < 0.05), and (3) summary statistics.

# Enrichment Analysis on Validated DEPs
library(clusterProfiler)
library(org.Rn.eg.db)
library(proharmed)
library(dplyr)

setwd("/home/ole/symbod_proteomics")

# Setup paths
base_output_dir <- "valid_DEPs/"

# Define the threshold combinations to process
fc_thresholds <- c(2.07, 2.49, 3.38)
stability_thresholds <- c(3, 5, 7)

# Paths to full DE results for background universe
ai_base_path <- "de_analysis_results_AI_RobNorm_scaled/"
as_base_path <- "de_analysis_results_AS_RobNorm_scaled/"

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
# Function to convert gene names to Entrez IDs
convert_to_entrez <- function(gene_names) {
  entrez_ids <- bitr(gene_names, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Rn.eg.db)
  return(entrez_ids$ENTREZID)
}

# Function to convert human gene symbols to rat gene symbols via proharmed
convert_human_to_rat_symbols <- function(human_gene_names) {
  cat("\n--- Converting human genes to rat orthologs ---\n")
  
  # Create a minimal dataframe for proharmed
  temp_df <- data.frame(HumanGenes = human_gene_names)
  
  # Use proharmed to map orthologs from human to rat
  orthologs_result <- proharmed::map_orthologs(
    data = temp_df, 
    gene_column = "HumanGenes", 
    organism = "human",
    tar_organism = "rat",
    res_column = "RatOrthologs",
    keep_empty = FALSE  # Remove genes without orthologs
  )
  
  orthologs_data <- orthologs_result$Modified_Data
  
  if (nrow(orthologs_data) == 0) {
    cat("No rat orthologs found\n")
    return(character(0))
  }
  
  # Extract rat gene symbols (they come as lists, need to unlist)
  rat_symbols <- unlist(orthologs_data$RatOrthologs)
  rat_symbols <- unique(rat_symbols[!is.na(rat_symbols)])
  
  cat("Found", length(rat_symbols), "rat ortholog symbols from", 
      length(human_gene_names), "human genes\n")
  
  return(rat_symbols)
}

# Function to run enrichment analysis with custom background
run_enrichment_with_background <- function(gene_names, dataset_name, universe) {
  # Clean gene names
  gene_names <- gene_names[!is.na(gene_names) & gene_names != "" & gene_names != "NA"]
  
  if (length(gene_names) == 0) {
    cat("No valid gene names found for", dataset_name, "\n")
    return(NULL)
  }
  
  # Convert to Entrez IDs
  gene_entrez <- convert_to_entrez(gene_names)
  
  cat("Converted", length(gene_entrez), "genes to Entrez IDs from", 
      length(gene_names), "input genes (", 
      round(100 * length(gene_entrez) / length(gene_names), 1), "% success)\n")
  
  if (length(gene_entrez) == 0) {
    cat("No genes could be converted to Entrez IDs for", dataset_name, "\n")
    return(NULL)
  }
  
  cat("Testing", length(gene_entrez), "genes against universe of", 
      length(universe), "genes for", dataset_name, "\n")
  
  # Skip enrichment for very small gene sets
  if (length(gene_entrez) < 10) {
    cat("WARNING: Gene set too small (", length(gene_entrez), 
        " genes) for reliable enrichment. Skipping.\n")
    return(NULL)
  }
  
  # Run enrichments with custom background
  bp <- enrichGO(gene = gene_entrez,
                 universe = universe,
                 OrgDb = org.Rn.eg.db,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1,  # Changed to 1 to get all results
                 qvalueCutoff = 1,  # Changed to 1 to get all results
                 readable = TRUE)
  
  mf <- enrichGO(gene = gene_entrez,
                 universe = universe,
                 OrgDb = org.Rn.eg.db, 
                 ont = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1,
                 qvalueCutoff = 1,
                 readable = TRUE)
  
  cc <- enrichGO(gene = gene_entrez,
                 universe = universe,
                 OrgDb = org.Rn.eg.db,
                 ont = "CC", 
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1,
                 qvalueCutoff = 1,
                 readable = TRUE)
  
  kegg <- enrichKEGG(gene = gene_entrez,
                     universe = universe,
                     organism = 'rno',
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1,
                     qvalueCutoff = 1)
  
  return(list(BP = bp, MF = mf, CC = cc, KEGG = kegg))
}
# ============================================================================

# ============================================================================
# LOOP THROUGH ALL THRESHOLD AND TOGGLE COMBINATIONS
# ============================================================================
for (fc_thresh in fc_thresholds) {
  for (stab_thresh in stability_thresholds) {
    
    cat("\n\n========================================\n")
    cat("PROCESSING: FC", fc_thresh, "| Stability", stab_thresh, "%\n")
    cat("========================================\n\n")
    
    # Define paths for this threshold combination
    thresh_dir <- file.path(base_output_dir, paste0("FC", fc_thresh, "_Stab", stab_thresh))
    validated_deps_file <- file.path(thresh_dir, "validated_DEPs.csv")
    
    # Check if validated DEPs file exists
    if (!file.exists(validated_deps_file)) {
      cat("WARNING: No validated_DEPs.csv found for FC", fc_thresh, "Stab", stab_thresh, "\n")
      cat("Skipping this combination...\n")
      next
    }
    
    # Read validated DEPs
    validated_deps <- read.csv(validated_deps_file, stringsAsFactors = FALSE)

    cat("Loaded", nrow(validated_deps), "validated DEPs across", 
    length(unique(validated_deps$Comparison)), "comparisons\n\n")
      
    # Process each comparison
    for (comparison in unique(validated_deps$Comparison)) {
      cat("\n========================\n")
      cat("Processing comparison:", comparison, "\n")
      
      # Get validated DEPs for this comparison
      comp_deps <- validated_deps[validated_deps$Comparison == comparison, ]
      
      # Separate first-level and second-level
      first_level_deps <- comp_deps[grepl("both fractions|Exclusive", comp_deps$Validation), ]
      second_level_deps <- comp_deps[grepl("Dominant in", comp_deps$Validation), ]
      
      # Load comparison-specific exception proteins
      exception_file <- file.path("network_enrichment_results", comparison, "exception_proteins.csv")
      exception_genes <- NULL
      exception_entrez_count <- 0
      if (file.exists(exception_file)) {
        exceptions <- read.csv(exception_file)
        exception_genes <- exceptions$Gene
        
        # Convert human to rat and get Entrez IDs
        exception_rat <- convert_human_to_rat_symbols(exception_genes)
        exception_entrez <- convert_to_entrez(exception_rat)
        exception_entrez_count <- length(exception_entrez)
      }

      # CREATE 3 SETS:
      # Set 1: First-level genes only
      set1_genes <- first_level_deps$Gene.Names
      
      # Set 2: First-level + exception proteins (converted to rat)
      set2_genes <- set1_genes
      if (!is.null(exception_genes)) {
        set2_genes <- unique(c(set1_genes, exception_rat))
      }
      
      # Set 3: First-level + second-level
      set3_genes <- unique(c(first_level_deps$Gene.Names, second_level_deps$Gene.Names))
      
      # Define the 3 sets
      sets_to_analyze <- list(
        list(name = "set1_first_level", genes = set1_genes),
        list(name = "set2_first_plus_network", genes = set2_genes),
        list(name = "set3_first_plus_second", genes = set3_genes)
      )
      
      # Read FULL DE results to get background universe
      ai_file_full <- file.path(ai_base_path, "de_results_raw.csv") 
      as_file_full <- file.path(as_base_path, "de_results_raw.csv")
      
      # Read full DE results
      full_ai <- read.csv(ai_file_full)
      full_as <- read.csv(as_file_full)
      
      # Filter for this comparison
      full_ai_comp <- full_ai[full_ai$Comparison == comparison, ]
      full_as_comp <- full_as[full_as$Comparison == comparison, ]
      
      if (nrow(full_ai_comp) == 0 || nrow(full_as_comp) == 0) {
        cat("Warning: No data found for comparison", comparison, "\n")
        next
      }
      
      # Get all gene names for background universe
      full_ai_genes <- unique(full_ai_comp$Gene.Names)
      full_as_genes <- unique(full_as_comp$Gene.Names)
      
      # Clean gene names
      valid_ai_genes <- full_ai_genes[!is.na(full_ai_genes) & 
                                        full_ai_genes != "" & 
                                        full_ai_genes != "NA"]
      valid_as_genes <- full_as_genes[!is.na(full_as_genes) & 
                                        full_as_genes != "" & 
                                        full_as_genes != "NA"]
      
      # Combine both backgrounds for universe
      all_background_genes <- unique(c(valid_ai_genes, valid_as_genes))
      # Convert background to Entrez first
      background_entrez <- convert_to_entrez(all_background_genes)
      
      cat("Background universe:", length(all_background_genes), "genes →", 
          length(background_entrez), "Entrez IDs\n")
      
      enrichment_summary <- data.frame()
      for (set_info in sets_to_analyze) {
        # Collect all enrichment results
        all_enrichment_results <- list()
        
        set_name <- set_info$name
        genes <- set_info$genes
        
        output_dir <- file.path(thresh_dir, "enrichment", set_name)
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        
        # For Set 2 expand background with exception proteins
        current_background <- background_entrez
        if (set_name == "set2_first_plus_network" && !is.null(exception_genes)) {
          current_background <- unique(c(background_entrez, exception_entrez))
        }
        
        enrichment_results <- run_enrichment_with_background(
          genes, 
          paste0(comparison, "_", set_name), 
          current_background
        )
        
        # Collect enrichment results (ALL pathways, not just significant)
        if (!is.null(enrichment_results)) {
          n_significant <- 0
          for (category in c("BP", "MF", "CC", "KEGG")) {
            if (!is.null(enrichment_results[[category]])) {
              result_data <- enrichment_results[[category]]@result
              
              if (nrow(result_data) > 0) {
                # Add metadata columns
                result_data$Category <- category
                result_data$Comparison <- comparison
                
                # Add significance flag
                result_data$Significant <- result_data$p.adjust < 0.05
                
                all_enrichment_results[[length(all_enrichment_results) + 1]] <- result_data
                
                # Count significant
                n_sig <- sum(result_data$p.adjust < 0.05)
                n_significant <- n_significant + n_sig
                
                cat("Found", n_sig, "significant /", nrow(result_data), 
                    "tested", category, "terms for", comparison, "\n")
              }
            }
          }
          
          # Save results for THIS set only
          if (length(all_enrichment_results) > 0) {
            combined_enrichment <- dplyr::bind_rows(all_enrichment_results)
            write.csv(combined_enrichment, 
                      file.path(output_dir, paste0(comparison, "_", set_name, "_enrichment_complete.csv")))
          }
          
          # Get the actual number of genes tested (Entrez IDs)
          genes_tested <- if (!is.null(enrichment_results)) {
            gene_entrez <- convert_to_entrez(genes)
            length(gene_entrez)
          } else {
            0
          }
          
          summary_row <- data.frame(
            Comparison = comparison,
            Set = set_name,
            Total_Genes_Input = length(genes),  # Gene symbols input
            Total_Genes_Tested = genes_tested,  # Entrez IDs actually tested
            Exception_Proteins = ifelse(set_name == "set2_first_plus_network", 
                                        length(exception_genes %||% character(0)), 0),
            Exception_Proteins_Tested = exception_entrez_count,
            Background_Universe = length(current_background),
            Enrichment_Performed = !is.null(enrichment_results),
            Significant_Pathways = n_significant
          )
          enrichment_summary <- rbind(enrichment_summary, summary_row)
        }
      } # End input gene set loop
      
      # Save summary ONCE per comparison (after all 3 sets)
      comparison_summary_dir <- file.path(thresh_dir, "enrichment")
      dir.create(comparison_summary_dir, recursive = TRUE, showWarnings = FALSE)
      write.csv(enrichment_summary, 
                file.path(comparison_summary_dir, paste0(comparison, "_enrichment_summary.csv")), 
                row.names = FALSE)
      cat("Summary saved for comparison:", comparison, "\n")
    } # END comparison loop
  } # END stability threshold loop
} # END FC threshold loop

cat("\n\n========================================")
cat("\n ALL ENRICHMENT ANALYSES COMPLETE")
cat("\n========================================\n")