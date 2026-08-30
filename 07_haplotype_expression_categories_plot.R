# ==============================================================================
# Draw a stacked bar chart of haplotype expression-category proportions (with adjustable legend size)
# ==============================================================================

# Step 0: Install and load the required R packages
# If 'tidyverse' is not installed, uncomment and run the following line first
# install.packages("tidyverse")

library(tidyverse) # Load core packages including ggplot2, dplyr, and readr
library(grid)      # Load grid to use the unit() function

# ==============================================================================
# User-defined settings
# ==============================================================================

# 1. List of file paths
file_paths <- c(
  "leaf_classification_count.txt",
  "stem_Classification_count.txt",
  "rhizome_classification_count.txt"
)

# 2. List of tissue names
tissue_names <- c(
  "Leaf",
  "Stem",
  "Rhizome"
)

# 3. Custom colors
custom_colors <- c(
  "Balanced" = "#D3D3D3",
  "AD" = "#F08080",
  "BD" = "#FFFACD",
  "CD" = "#87CEEB",
  "DD" = "#FFDAB9",
  "AS" = "#E6E6FA",
  "BS" = "#40E0D0",
  "CS" = "#AFEEEE",
  "DS" = "#FFBBFF"
)

# ==============================================================================
# Data processing and plotting
# ==============================================================================

all_data <- data.frame()
for (i in seq_along(file_paths)) {
  temp_data <- read.table(file_paths[i], header = FALSE, sep = "\t", skip = 1,
                          col.names = c("Category", "Count"))
  temp_data$Tissue <- tissue_names[i]
  all_data <- rbind(all_data, temp_data)
}

plot_data <- all_data %>%
  group_by(Tissue) %>%
  mutate(
    TotalCount = sum(Count),
    Proportion = Count / TotalCount
  ) %>%
  ungroup() %>%
  mutate(Category_simplified = case_when(
    Category == "Balanced" ~ "Balanced",
    grepl("_dominant", Category) ~ paste0(sub("_dominant", "", Category), "D"),
    grepl("_suppressed", Category) ~ paste0(sub("_suppressed", "", Category), "S"),
    TRUE ~ as.character(Category)
  ))

category_order <- c("Balanced", "AD", "BD", "CD", "DD", "AS", "BS", "CS", "DS")
tissue_order <- rev(tissue_names)

plot_data$Category_simplified <- factor(plot_data$Category_simplified, levels = category_order)
plot_data$Tissue <- factor(plot_data$Tissue, levels = tissue_order)

# Begin plotting
stacked_bar_plot <- ggplot(plot_data, aes(x = Proportion, y = Tissue, fill = Category_simplified)) +
  geom_col(width = 0.6, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = custom_colors, name = NULL, drop = FALSE) +
  scale_x_continuous(
    labels = scales::percent_format(),
    expand = c(0, 0)
  ) +
  labs(
    title = " Expression categories Proportions",
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  # Add parameters in theme() to adjust the legend size
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold",size = 9),
    # 1. Adjust the size of legend keys
    legend.key.size = unit(0.3, 'cm'), # Try values such as 0.3, 0.4, or 0.5
    
    # 2. Adjust the legend text size
    legend.text = element_text(size = 7), # Try values such as 7, 8, or 10
    # Adjust the size of Y-axis labels such as "Leaf" and "Rhizome"
    axis.text.y = element_text(size = 7), # Try values such as 10, 12, or 14
    axis.text.x = element_text(size = 7), # Try values such as 10, 12, or 14
    # margin(top, right, bottom, left, unit)
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    # =====================================================================
  ) +
  guides(fill = guide_legend(nrow = 1))

# Display the plot
print(stacked_bar_plot)

# Save the plot
ggsave("haplotype_expression_plot.pdf", plot = stacked_bar_plot, width = 7, height = 2, dpi = 300)
