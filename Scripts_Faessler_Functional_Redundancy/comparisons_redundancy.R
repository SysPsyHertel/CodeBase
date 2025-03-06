library(SYNCSA)
library(PNWColors)
library(FunRed)
library(MCMCpack)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
set.seed(333)
color_palette <- pnw_palette("Bay", 3)

## This file simulates the scenario of a keystone species and compares it with traditional approaches of functional redundancy
# Define the species set and functional values
species_count <- 100  # Set species count to 100
f_values <- list("100" = c(1, rep(0, 99)))

# Create a function to compute functional redundancy or Rao's diversity
compute_results <- function(species_count, f_values, is_rao = FALSE) {
  abundance_list <- list()
  fr_results <- list()
  
  for (i in 1:species_count) {
    abundance <- rep(0, species_count)
    abundance[1:i] <- 1 / i
    abundance_list[[i]] <- abundance
    
    if (is_rao) {
      community_matrix <- matrix(
        abundance,
        nrow = 1, byrow = TRUE,
        dimnames = list(c("Site_1"), paste0("Sp_", 1:species_count))
      )
      traits_matrix <- matrix(
        f_values[[as.character(species_count)]],
        ncol = 1,
        dimnames = list(paste0("Sp_", 1:species_count), c("Trait_1"))
      )
      fr_results[[i]] <- rao.diversity(community_matrix, traits = traits_matrix)$FunRedundancy
    } else {
      fr_results[[i]] <- FunRed::fredundancy(f_values[[as.character(species_count)]], abundance)$abundance_based
    }
  }
  
  # Prepare data for plotting with species from 1 to species_count
  plot_data <- data.frame(
    species = 1:species_count,  # Changed to range from 1 to species_count
    fr = sapply(fr_results, identity),
    group = paste0(species_count, "_species")
  )
  
  return(plot_data)
}


# Run the compute_results function for 100 species
results <- compute_results(species_count, f_values, is_rao = FALSE)

# Calculate -log(100) once for reuse
y_intercept_value <- -log(100)

# Create a new data frame for reference-based redundancy (-log(100) for each x value)
reference_data <- data.frame(species = seq(0, 100, by = 1), 
                             fr = rep(-log(100), length.out = 101))  # -log(100) for all x values

fr_plot <- ggplot(results, aes(x = species, y = fr)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.4) +  
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "black", size = 0.4) + 
  # Black dots with legend label for Sample/abundance-based
  geom_point(aes(color = "Sample/abundance-based"), size = 1.5, alpha = 0.5) +  
  # Blue points with alpha=0.5 for reference-based redundancy
  geom_point(data = reference_data, aes(x = species, y = fr, color = "Reference-based"), 
             size = 1.5, alpha = 0.5) +  
  # Horizontal lines
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 9, color = "black"),  
    axis.title.y = element_text(size = 9, color = "black", vjust = -10),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(size = 9, hjust = 0.5),
    legend.position = c(0.77, 0.8),
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    legend.key.size = unit(2.2, "mm"),
    legend.background = element_rect(color = "black", size = 0.2, linetype = "solid", fill = "white"),
    legend.spacing.y = unit(1.5, "mm")
  ) +
  labs(
    x = "Species richness", 
    y = "Functional redundancy",
    title = "Trait-based functional redundancy"
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 10),
    labels = seq(0, 100, by = 10)
  ) +
  scale_y_continuous(
    breaks = c(seq(-5, 0, by = 1), -log(100)),
    labels = c(seq(-5, 0, by = 1), "-log(100)")
  ) +
  scale_color_manual(values = c(
    "Sample/abundance-based" = "black", 
    "Reference-based" = color_palette[1]  # Ensure color_palette[1] is defined correctly
  )) +  # Legend colors
  guides(color = guide_legend(title = NULL, ncol = 1))  # Arrange legend in one column

# Print the plot
print(fr_plot)

# Run the compute_results function for 100 species
results_rao <- compute_results(species_count, f_values, is_rao = TRUE)



# Create the rao_plot with reduced font sizes
rao_plot <- ggplot(results_rao, aes(x = species, y = fr)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.4) +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.4) + 
  geom_point(color = "black", size = 1.5, alpha=0.5) +  # Blue dots
  theme_bw() +
  ylim(0,1) + 
  theme(
    axis.title.x = element_text(size = 9, color = "black"),  # Decrease size of X-axis title
    axis.title.y = element_text(size = 9, color = "black"),  # Decrease size of Y-axis title
    axis.text.x = element_text(size = 9, color = "black"),   # Decrease size of X-axis text
    axis.text.y = element_text(size = 9, color = "black"),   # Decrease size of Y-axis text
    axis.line = element_line(color = "black"),                 # X and Y axis lines color
    panel.grid = element_blank(),                             # No grid
    plot.title = element_text(size = 9, hjust = 0.5)         # Set title size here
  ) +
  labs(
    x = "Species richness", 
    y = bquote("Functional redundancy (FR) \n (based on Rao's quadratic entropy)"),
    title = "Community-based functional redundancy"
  )+
  scale_x_continuous(
    breaks = seq(0, 100, by = 10),    # Set x-axis breaks to 0, 10, 20, ..., 100
    labels = seq(0, 100, by = 10)      # Set x-axis labels to match the breaks
  )

#rao_fr_plots <- plot_grid(fr_plot, rao_plot, ncol = 2)

# Create a common title as a plot
title_top_plot <- ggdraw() + 
  draw_label("Simulating trait-based- and community-based functional redundancy of a single-species contribution", 
             size = 10, fontface = 'bold') + 
  theme_void() + 
  theme(plot.background = element_rect(fill = "white", color = NA))  # Set the background to white

# Combine the title with the combined plot
#final_plot_top <- plot_grid(title_top_plot, rao_fr_plots, ncol = 1, rel_heights = c(0.1, 1))

rao_fr_plots <- plot_grid(fr_plot, rao_plot, ncol = 2, labels = "AUTO")