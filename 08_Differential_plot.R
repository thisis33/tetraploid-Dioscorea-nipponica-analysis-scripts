# --- 1. Configuration ---

# Load the required libraries
library(ggplot2)
library(dplyr)
library(ggrepel) 

# a. Define the CSV filenames and corresponding tissue names for the three tissues
files_and_tissues <- list(
  list(file = "leaf_deg.csv", name = "Leaf"),
  list(file = "rhizome_deg.csv", name = "Rhizome"),
  list(file = "stem_deg.csv", name = "Stem")
)

# b. Define the filename containing the list of highlighted genes
highlight_file <- "highlight_genes.csv"

# c. Custom colors
color_upregulated   <- "#FD9C99"
color_downregulated <- "#BCCF90"
color_highlighted   <- "#ff7f0e"

# d. Define the output image filename
output_filename <- "Differential_expression_plot.pdf"

plot_data_output_file <- "plot_data_output.csv"

# --- 2. Read and merge data  

list_of_dfs <- list()
for (item in files_and_tissues) {
  tryCatch({
    df <- read.csv(item$file)
    df$tissue <- item$name
    list_of_dfs[[item$name]] <- df
  }, error = function(e) {
    stop(paste("Error: Unable to read tissue file '", item$file, "'.", sep=""))
  })
}
combined_df <- bind_rows(list_of_dfs)

tryCatch({
  highlight_df <- read.csv(highlight_file)
}, error = function(e) {
  stop(paste("Error: Unable to read the highlighted-gene file '", highlight_file, "'.\n", e))
})

# --- 3. Automatically identify and rename columns --

lfc_col_actual <- grep("^log2FoldChange$", names(combined_df), ignore.case = TRUE, value = TRUE)
if (length(lfc_col_actual) != 1) stop("A unique 'log2FoldChange' column could not be found in the DEG files.")
gene_col_actual <- grep("^gene$", names(combined_df), ignore.case = TRUE, value = TRUE)
if (length(gene_col_actual) != 1) stop("A unique 'gene' column could not be found in the DEG files.")
combined_df_renamed <- combined_df %>%
  rename(internal_log2fc = !!lfc_col_actual, internal_gene = !!gene_col_actual)

gene_id_col_actual <- grep("^gene.id$", gsub(" ", ".", names(highlight_df)), ignore.case = TRUE, value = TRUE)
if (length(gene_id_col_actual) != 1) stop("A unique 'Gene id' or 'Gene_id' column could not be found in the highlighted-gene file.")
gene_name_col_actual <- grep("^gene.name$", gsub(" ", ".", names(highlight_df)), ignore.case = TRUE, value = TRUE)
if (length(gene_name_col_actual) != 1) stop("A unique 'Gene name' or 'Gene_name' column could not be found in the highlighted-gene file.")
tissue_col_actual <- grep("^tissue$", names(highlight_df), ignore.case = TRUE, value = TRUE)
if (length(tissue_col_actual) != 1) stop("A unique 'Tissue' column could not be found in the highlighted-gene file.")
highlight_df_renamed <- highlight_df %>%
  rename(internal_gene_id = !!gene_id_col_actual, internal_gene_name = !!gene_name_col_actual, internal_tissue = !!tissue_col_actual)

cat("\n--- 3a. First six rows after standardizing DEG column names (combined_df_renamed) ---\n")
print(head(combined_df_renamed))

plot_data <- combined_df_renamed %>%
  left_join(highlight_df_renamed, by = c("internal_gene" = "internal_gene_id", "tissue" = "internal_tissue")) %>%
  mutate(
    regulation = ifelse(internal_log2fc > 0, "Upregulated", "Downregulated"),
    highlight_status = ifelse(!is.na(internal_gene_name), "Highlighted", regulation)
  )

# Convert the 'tissue' column to a factor and arrange it in the specified order
plot_data$tissue <- factor(plot_data$tissue, levels = c("Leaf", "Stem", "Rhizome"))

# ---  Print statements ---
cat("\n--- 3c. First six rows of the final data used for plotting (plot_data) ---\n")
cat("Please verify that the 'internal_gene_name', 'regulation', and 'highlight_status' columns were generated correctly.\n")
print(head(plot_data))

# --- Save plot_data to a file ---
write.csv(plot_data, file = plot_data_output_file, row.names = FALSE)
cat(paste("\n--- The final plotting data have been saved successfully to:", plot_data_output_file, "---\n"))

cat("\n--- End of printed output; beginning plotting ---\n")
# --- End of printed output ---

color_map <- c(
  "Upregulated" = color_upregulated,
  "Downregulated" = color_downregulated,
  "Highlighted" = color_highlighted
)

# --- 4. Draw the plot using ggplot2 ---

p <- ggplot(plot_data, aes(x = tissue, y = internal_log2fc, color = highlight_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  
 
  
  # 1. Draw all non-highlighted points with the size set to 2
  geom_jitter(
    data = filter(plot_data, highlight_status != "Highlighted"), 
    width = 0.25, 
    size = 0.5,     # <---  Change this value to adjust the size of ordinary points
    alpha = 0.6
  ) +
  
  # 2. Draw highlighted points on the upper layer with the size set to 4 for greater emphasis
  geom_jitter(
    data = filter(plot_data, highlight_status == "Highlighted"), 
    width = 0.25, 
    size = 2,     # <---  Change this value to adjust the size of highlighted points
    alpha = 1
  ) +
  
 
  
  scale_color_manual(values = color_map, name = "Regulation Status") +
  
  geom_text_repel(
    data = filter(plot_data, highlight_status == "Highlighted"),
    aes(label = internal_gene_name),
    size = 3,   # <---  Change this value to adjust the size of highlighted text
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey50'
  ) +
  
  labs(
    title = "Gene Expression across Tissues",
    x = "Tissue",
    y = "log2 Fold Change"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )

print(p)

# --- 5. Save the plot to a file ---
ggsave(output_filename, plot = p, width = 10, height = 8)
cat(paste("\nThe plot has been saved successfully as '", output_filename, "'\n", sep=""))
