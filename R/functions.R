#Function to check if files provided in main.R arguments exist.
arg_check = funtion(file) {
    if (!file.exists(file)) {
    stop(paste("This file does not exist:", file))
  }
}

# Function to conduct DESeq2 differential expression analysis
# This function performs DESeq2 analysis on the given counts matrix and model data
deseq = function(counts, model, condition) {
  # Ensure that the sample names in the counts matrix match those in the model matrix
  if (all(colnames(counts) == model$sample)) {
    print("MATCH")  # Print "MATCH" if the column names and model sample names align
  } else {
    stop("MISMATCH: STOP ANALYSIS RIIIIIIIIGHT NOW")  # Stop the analysis if they do not match
  }
  
  # Prepare the metadata for DESeq2 and run the differential expression analysis
  colData = data.frame(row.names = colnames(counts), group = condition)  # Create metadata frame
  dds = DESeqDataSetFromMatrix(countData = counts_filtered, colData = colData, design = ~ group)  # Create DESeq2 object
  data.frame(dds$group)  # View the group assignments
  dds = DESeq(dds)  # Perform DESeq2 analysis
}

# Function to conduct GSEA on differentially expressed genes
# This function conducts Gene Set Enrichment Analysis (GSEA) on DESeq2 results
gsea <- function(df, set) {
  tryCatch({
    # Convert ENSEMBL IDs to gene symbols
    ensembl_ids = row.names(df)  # Extract ENSEMBL IDs from DE results
    ensembl_ids = gsub("\\..*", "", ensembl_ids)  # Remove version numbers from ENSEMBL IDs
    gene_symbols = AnnotationDbi::select(org.Hs.eg.db, keys = ensembl_ids,  # Map ENSEMBL IDs to gene symbols
                                         columns = c("SYMBOL"),
                                         keytype = "ENSEMBL")
    
    # Prepare the dataframe for GSEA by merging gene symbols and DESeq2 results
    df$ENSEMBL = row.names(df)  # Add ENSEMBL IDs to the dataframe
    df$ENSEMBL = gsub("\\..*", "", df$ENSEMBL)  # Remove version numbers again
    gene_symbols$Gene = gene_symbols$SYMBOL  # Rename the SYMBOL column to Gene
    df = na.omit(merge(gene_symbols, df, by = "ENSEMBL"))  # Merge gene symbols and DE results, removing NAs
    df = df[ ,-(1:2)]  # Remove ENSEMBL IDs and SYMBOL columns
    
    # Prepare the gene list for GSEA
    msigdb = msigdbr(species = "Homo sapiens", category = set)  # Load the MSigDB gene sets
    gene_list = df$log2FoldChange  # Use the log2 fold change for GSEA
    names(gene_list) = df$Gene  # Assign gene names to the gene list
    gene_list = sort(gene_list, decreasing = TRUE)  # Sort the gene list in descending order
    
    # Run the GSEA analysis using the MSigDB gene sets
    gsea_res = GSEA(geneList = gene_list,
                    TERM2GENE = msigdb %>% dplyr::select(gs_name, gene_symbol),  # Use MSigDB gene sets
                    pvalueCutoff = 0.1,  # Set p-value cutoff for GSEA results
                    verbose = FALSE)  # Suppress verbose output
    gsea_res = data.frame(gsea_res)  # Convert GSEA results to a dataframe
    gsea_res$gene_set = set  # Add the gene set category to the results
    
    return(gsea_res)  # Return the GSEA results dataframe
    
  }, error = function(e) {
    # Handle the error (e.g., log the error message and continue)
    print(paste("Error occurred during GSEA analysis for set:", set))
    print(paste("Error message:", conditionMessage(e)))
    
    # Optionally, return NULL or some placeholder value if you want to continue
    return(NULL)
  })
}

# Function to combine GSEA results
bind_gsea = function(df_list) {
  result = data.frame()
  for (i in 1:length(df_list)) {
    df = df_list[[as.numeric(i)]]
    
    # Try to add a column and rbind, catching any errors
    tryCatch({
      # Check if the data frame is empty
      if (nrow(df) > 0) {
        # Bind to the result
        result = rbind(result, df)
      } else {
        print(paste("Skipping empty data frame at index", i))
      }
    }, error = function(e) {
      # Handle the error (e.g., skip this data frame)
      print(paste("Error in data frame at index", i, ":", conditionMessage(e)))
    })
  }
  return(result)
}

# Function to perform cell deconvolution
# This function uses the immunedeconv package to perform cell deconvolution
deconv = function(matrix, tool, name) {
  name = immunedeconv::deconvolute(matrix, tool)  # Perform cell deconvolution with the specified tool
  name = data.frame(name)  # Convert the result to a dataframe
  row.names(name) = name$cell_type  # Set row names to the cell type
  name = name[, -1]  # Remove the cell type column from the dataframe
  colnames(name) = gsub("\\.", "-", colnames(name))  # Replace dots with dashes in the column names (sample IDs)
  return(name)  # Return the cell deconvolution results
}

# Function to filter cell deconvolution results by subtype
# This function filters the cell deconvolution results by the predefined experimental groups
filter_subtypes = function(tool, name) {
  # Create variable names for different subtypes and tools
  bnn_pos_var = paste0("bnn_pos_", tool)
  bnn_neg_var = paste0("bnn_neg_", tool)
  twt_pos_var = paste0("twt_pos_", tool)
  twt_neg_var = paste0("twt_neg_", tool)
  all_pos_var = paste0("all_pos_", tool)
  all_neg_var = paste0("all_neg_", tool)
  
  # Assign the filtered data to the appropriate variables
  assign(bnn_pos_var, data.frame(name[, colnames(name) %in% bnn_pos]))
  assign(bnn_neg_var, data.frame(name[, colnames(name) %in% bnn_neg]))
  assign(twt_pos_var, data.frame(name[, colnames(name) %in% twt_pos]))
  assign(twt_neg_var, data.frame(name[, colnames(name) %in% twt_neg]))
  assign(all_pos_var, data.frame(name[, colnames(name) %in% ecdna_pos]))
  assign(all_neg_var, data.frame(name[, colnames(name) %in% ecdna_neg]))
  
  # Return a list of the filtered results for each subtype
  return(list(
    bnn_pos = get(bnn_pos_var),
    bnn_neg = get(bnn_neg_var),
    twt_pos = get(twt_pos_var),
    twt_neg = get(twt_neg_var), 
    all_pos = get(all_pos_var),
    all_neg = get(all_neg_var)
  ))
}

# Function to perform the Mann-Whitney U test (Wilcoxon rank-sum test)
# This function compares two groups of cell deconvolution data using the Wilcoxon test
mannU = function(list, comp1, comp2, tool) {
  comp1_df = data.frame(list[comp1])  # Get the first comparison group data
  comp2_df = data.frame(list[comp2])  # Get the second comparison group data
  df = data.frame(row.names(comp1_df))  # Create an empty dataframe to store results
  
  # Loop through each cell type and perform the Wilcoxon rank-sum test
  for (i in 1:nrow(comp1_df)) {
    posValues = as.numeric(comp1_df[i, ])  # Extract positive group values for the cell type
    negValues = as.numeric(comp2_df[i, ])  # Extract negative group values for the cell type
    combined = c(posValues, negValues)  # Combine the two groups for the test
    
    # Perform the Wilcoxon test
    test_result = wilcox.test(posValues, negValues)
    df$pval[i] = test_result$p.value  # Store the p-value for the test
    df$padj = p.adjust(df$pval, method = "fdr")  # Adjust p-values using the false discovery rate (FDR) method
    df$logfc[i] = log(mean(posValues)/mean(negValues), 2)  # Calculate log2 fold change between the groups
  }
  
  # Add metadata columns to the results dataframe
  df$tool = tool  # Add the deconvolution tool used
  colnames(df) = c("cell_type", "pval", "padj", "logfc", "tool")  # Set the column names for the dataframe
  
  return(df)  # Return the dataframe containing Mann-Whitney U test results
}

