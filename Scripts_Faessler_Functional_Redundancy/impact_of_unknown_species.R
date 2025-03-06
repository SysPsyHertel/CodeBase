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

## This file simulates functional contributions and abundances based
# on dirichlet parametrization
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




############## remove species
# Define function to generate abundance with different evenness levels
generate_abundance <- function(num_species, alpha_value) {
  alpha <- rep(alpha_value, num_species)  # Dirichlet parameters
  abundance <- rdirichlet(1, alpha)  # Generate relative frequencies
  return(abundance)
}

# Define function for calculating FR values
calculate_FR <- function(f, a) {
  FR <- FunRed::fredundancy(f, a, n_reference=100)
  return(list(FR_sample = FR$sample_based,
              FR_reference = FR$reference_based,
              FR_abundance = FR$abundance_based))
}

# Simulation function
simulate_FR <- function(num_species = 100, num_steps = 9, num_simulations = 1000) {
  results <- data.frame()
  
  for (sim in 1:num_simulations) {
    # Assign even or uneven distribution randomly
    if (runif(1) < 0.5) {
      alpha_value <- 1  # Uneven distribution (alpha = 1)
      distribution_type <- "alpha = 1"
    } else {
      alpha_value <- 5  # More even distribution (alpha = 5)
      distribution_type <- "alpha = 5"
    }
    
    f <- generate_abundance(num_species, alpha_value)
    a <- generate_abundance(num_species, alpha_value)
    
    for (step in 1:num_steps) {
      # Calculate FR values
      FR <- calculate_FR(f, a)
      
      # Store results
      results <- rbind(results, data.frame(Simulation = sim,
                                           Step = step,
                                           FR_sample = FR$FR_sample,
                                           FR_reference = FR$FR_reference,
                                           FR_abundance = FR$FR_abundance,
                                           Parameter = distribution_type))
      
      # Remove 10 species and renormalize
      non_zero_indices <- which(f > 0)
      if (length(non_zero_indices) > 10) {
        set_zero_indices <- sample(non_zero_indices, 10)
        f[set_zero_indices] <- 0
        a[set_zero_indices] <- 0
      }
      
      f <- f / sum(f)
      a <- a / sum(a)
    }
  }
  return(results)
}

# Run simulation
simulation_results <- simulate_FR()

# Reshape data to long format for ggplot
results_long <- simulation_results %>%
  pivot_longer(cols = starts_with("FR_"), names_to = "Measure", values_to = "FR_value")

#results_long$Parameter <- factor(results_long$Parameter, levels = c("alpha = 5", "alpha = 1"))
plot_FR <- function(data, measure) {
  
  # Check unique values in the Parameter column to confirm their presence
  unique_parameters <- unique(data$Parameter)
  print(unique_parameters)  # To see what the unique values are
  
  # Manually define color mapping
  color_mapping <- c("alpha = 1" = color_palette[3], "alpha = 5" = color_palette[1])
  
  # Create the ggplot
  ggplot(data %>% filter(Measure == measure),
         aes(x = factor(Step), y = FR_value, fill = Parameter)) +
    geom_jitter(aes(color = Parameter), size = 0.3, alpha = 0.5, 
                position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8), 
                shape = 21, stroke = 0.1) +  
    geom_violin(alpha = 0.5, trim = TRUE, 
                position = position_dodge(width = 0.8), 
                size = 0.2) + 
    labs(title = paste0(measure, " FR"), x = "Unknown species", y = "Functional redundancy") +
    scale_x_discrete(labels = seq(0, 80, by = 10)) +
    scale_y_continuous(limits = c(-3, 0)) +
    
    # Define fill and color scales with appropriate labels
    scale_fill_manual(
      values = color_mapping,  # Use the color mapping for the fill
      labels = c("alpha = 1" = "alpha = 1 (uneven distributions)", "alpha = 5" = "alpha = 5 (even distributions)")  # Descriptive labels
    ) +
    scale_color_manual(
      values = color_mapping,  # Use the color mapping for the color
      labels = c("alpha = 1" = "alpha = 1 (uneven distributions)", "alpha = 5" = "alpha = 5 (even distributions)")  # Descriptive labels
    ) +
    ylim(-2.5, 0) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 9, color = "black"),
      axis.title.y = element_text(size = 9, color = "black"),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      strip.text = element_text(size = 9, color = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
      legend.position = "right",  # Position the legend on the right
      legend.direction = "horizontal",  # Make legend items horizontal
      legend.title = element_text(size = 7),  # Adjust the legend title size
      legend.text = element_text(size = 7),   # Adjust the legend text size
      legend.key.size = unit(0.5, "cm"),  # Adjust the size of the legend keys (color boxes)
      legend.key.height = unit(0.5, "cm"),  # Adjust the height of the legend keys
      legend.key.width = unit(0.8, "cm")  # Adjust the width of the legend keys
    )
}

# Generate plots with descriptive titles (including the legend)
plot_sample <- plot_FR(results_long, "FR_sample") + labs(title = "Sample taxon-based functional redundancy")
plot_reference <- plot_FR(results_long, "FR_reference") + labs(title = "Reference taxon-based functional redundancy")
plot_abundance <- plot_FR(results_long, "FR_abundance") + labs(title = "Abundance-based functional redundancy")

# Ensure legend is visible by printing one of the plots (e.g., plot_sample)
print(plot_sample)

# Extract the legend and make its background white
legend_bottom <- cowplot::get_legend(plot_sample) + 
  theme(
    # Set the legend's background to white
    legend.background = element_rect(fill = "white", color = NA),
    
    # Ensure the legend keys (colored boxes) are visible
    legend.key = element_rect(fill = "white", color = NA),  # Ensures keys stay visible with a white background
    
    # Optional: Set legend title and text size to ensure readability
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  )

# Display the legend with a white background
legend_bottom


# Display the legend with a white background
legend_bottom

# Remove the legends from individual plots
plot_sample <- plot_sample + theme(legend.position = "none")
plot_reference <- plot_reference + theme(legend.position = "none")
plot_abundance <- plot_abundance + theme(legend.position = "none")

# Combine the plots without the legends
combined_plot_bottom <- plot_grid(plot_sample, plot_reference, plot_abundance, ncol = 3)

# Add a title above the plots
title_plot_bottom <- ggdraw() + 
  draw_label("Impact of unknown species on measures of functional redundancy", 
             size = 10, fontface = 'bold') + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Combine everything (title, plots, and legend)
final_plot_bottom <- plot_grid(
  title_plot_bottom, 
  combined_plot_bottom, 
  legend_bottom, 
  ncol = 1, 
  rel_heights = c(0.1, 1, 0.1)  # Adjust the relative height for the legend
)


F2 <- plot_grid(rao_fr_plots, final_plot_bottom, nrow=2,  rel_heights=c(1,1.2), labels=c("", "C"))


# Save the final figure
ggsave("F:/Functional_Redundancy_CSBJ/ReSubmission/Figures/Figure_2/Figure_2.svg", 
       plot = F2, 
       width = 190, height = 150, units = "mm", dpi=600, scale=1)




