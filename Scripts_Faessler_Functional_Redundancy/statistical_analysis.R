library(xlsx)
library(ggplot2)
library(data.table)
library(gridExtra)
library(reshape2)
library(tidyr)
library(dplyr)
library(cowplot)

setwd(r'(C:\Users\faesslerd\Documents\Projects\CRC_Redundancy\Final_Vers)')
dsecr <- read.csv('Microbe_Secretion.csv')
d <- read.csv('RedundancyMeasures.csv')

## VMH identifiers from the cobratoolbox
VMH_metabolites <- fread('VMH_Metabolites.tsv')

mets <- names(d)[grep("^KL_Sample", names(d))]
mets <- gsub(".*KL_Sample_", "", mets)

d$healthstat <- as.factor(d$healthstat)
measures <- c("KL_Sample", "KL_Reference", "KL_Abundance")


##############################################
#### Interrelations redundancy measures   ####
##############################################

dc <- d
Rsqtable <- data.frame(VMHId=character(),
                       Producing.Species=double(),
                       Sample_Reference=double(),
                       Sample_Abundance=double(),
                       Reference_Abundance=double(),
                       stringsAsFactors=FALSE)

plot_list <- vector("list", length(mets))
cnt <- 0
# Extract intercepts from frm12
coefs <- c()
for(i in mets){
  cnt <- cnt+1
  frm12 <- as.formula(paste(paste0("KL_Sample_",i), paste(paste0("KL_Reference_",i), collapse=" + "), sep=" ~ "))
  frm13 <- as.formula(paste(paste0("KL_Sample_",i), paste(paste0("KL_Abundance_",i), collapse=" + "), sep=" ~ "))
  frm23 <- as.formula(paste(paste0("KL_Reference_",i), paste(paste0("KL_Abundance_",i), collapse=" + "), sep=" ~ "))
  
  fit1 <- lm(frm12, data = d)
  coefs <- c(coefs, as.numeric(fit1$coefficients[1]))
  fit2 <- lm(frm13, data = d)
  fit3 <- lm(frm23, data = d)

  Rsqtable[nrow(Rsqtable)+1, 1] <- i
  Rsqtable[nrow(Rsqtable), 2:ncol(Rsqtable)] <- c(mean(d[, paste0("N_",i)]), summary(fit1)$r.squared, summary(fit2)$r.squared, summary(fit3)$r.squared)
  print(i)
}

RSQT <- Rsqtable
RSQT$Metabolite <- VMH_metabolites$fullName[match(RSQT$VMHId, VMH_metabolites$abbreviation)]
  

Mapping <- RSQT[, c("Metabolite", "VMHId")]
ref <- c()

for(i in mets){
  metexch <- dsecr$X[grepl(paste0("_", i, "$"), dsecr$X)]
  print(length(metexch))
  ref <- c(ref, length(metexch))
}

RSQT$NoReference <- ref
colnames(RSQT)
new_order <- c("Metabolite", "VMHId", "NoReference", "Producing.Species", "Sample_Reference", "Sample_Abundance", "Reference_Abundance") # Define the new order
RSQT <- RSQT[, new_order]
#write.xlsx(RSQT, "Tables/TableS3.xlsx", row.names = FALSE)

combined_data <- bind_rows(
  transform(d, Category = "Glycine", X = KL_Sample_gly, Y = KL_Reference_gly),
  transform(d, Category = "Tryptamine", X = KL_Sample_trypta, Y = KL_Reference_trypta),
  transform(d, Category = "Glucosamine", X = KL_Sample_gam, Y = KL_Reference_gam),
)

# Define custom colors for each category
custom_colors <- c("Glycine" = "blue", "Tryptamine" = "purple", "Glucosamine" = "orange")
fig_3a <- ggplot(combined_data, aes(x = X, y = Y, color = Category)) +
  geom_point() +
  xlim(-5.5, 0) +
  ylim(-5.5, 0) +
  labs(x = "Sample taxon-based redundancy", y = "Reference taxon-based redundancy", color = "Metabolite") +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(x = 0, y = -1, label = "y = x", color = "black", hjust = 1, vjust = 1, size = 5) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = c(0.85, 0.15), # Adjust the legend position inside the plot
    aspect.ratio = 1,  # Ensure the aspect ratio is square
    plot.title = element_text(hjust = 0)  # Align the title to the left
  ) + coord_fixed(ratio = 1)

## Plot Number of species in the reference vs. Regression coefficients
data <- data.frame(Prod.Species = RSQT$Producing.Species, No.Reference=RSQT$NoReference, coefs=coefs)

## Average number of species that can perform the funktion
model <- lm(coefs ~ log(Prod.Species), data=data)

new_data <- data.frame(Prod.Species = seq(min(data$Prod.Species), max(data$Prod.Species), length.out = 100))
predictions <- predict(model, newdata = new_data, interval = "confidence")

predictions_df <- data.frame(Prod.Species = new_data$Prod.Species, 
                             fit = predictions[, "fit"], 
                             lwr = predictions[, "lwr"], 
                             upr = predictions[, "upr"])


p_value <- coef(summary(model))[2, 4]
r_squared <- summary(model)$r.squared
p_value_label <- paste("p-value =", format(p_value, digits = 2))
r_squared_label <- paste("R-squared =", format(r_squared, digits = 2))

# Plotting
fig_3b <- ggplot(data) +
  geom_point(mapping = aes(x = Prod.Species, y = coefs)) +
  xlab("Average number of species in the community that can secrete the metabolite") +
  ylab("Intercept") +
  theme_bw() +
  geom_line(data = predictions_df, aes(x = Prod.Species, y = fit), color = "blue", size = 1.2) +
  geom_ribbon(data = predictions_df, aes(x = Prod.Species, ymin = lwr, ymax = upr), alpha = 0.2) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black")
  ) +
  annotate("text", x = max(data$Prod.Species), y = min(data$coefs)-0.15,
           label = p_value_label,
           hjust = 1, vjust = 0, size = 5, color = "black") +
  annotate("text", x = max(data$Prod.Species), y = min(data$coefs) + 0.15,
           label = r_squared_label,
           hjust = 1, vjust = 0, size = 5, color = "black")
build_plot <- ggplot_build(fig_3b)
xlim <- build_plot$layout$panel_scales_x[[1]]$range$range
ylim <- build_plot$layout$panel_scales_y[[1]]$range$range

x_range <- diff(xlim)
y_range <- diff(ylim)
ratio <- x_range / y_range
fig_3b <- fig_3b + coord_fixed(ratio = ratio)

## No Reference
model <- lm(coefs ~ log(No.Reference), data=data)

new_data <- data.frame(No.Reference = seq(min(data$No.Reference), max(data$No.Reference), length.out = 100))
predictions <- predict(model, newdata = new_data, interval = "confidence")

predictions_df <- data.frame(No.Reference = new_data$No.Reference, 
                             fit = predictions[, "fit"], 
                             lwr = predictions[, "lwr"], 
                             upr = predictions[, "upr"])

p_value <- coef(summary(model))[2, 4]
r_squared <- summary(model)$r.squared
p_value_label <- paste("p-value =", format(p_value, digits = 2))
r_squared_label <- paste("R-squared =", format(r_squared, digits = 2))

# Plotting
fig_3c <- ggplot(data) +
  geom_point(mapping = aes(x = No.Reference, y = coefs)) +
  xlab("Number of species in the reference that can secrete the metabolite") +
  ylab("Intercept") +
  theme_bw() +
  geom_line(data = predictions_df, aes(x = No.Reference, y = fit), color = "blue", size = 1.2) +
  geom_ribbon(data = predictions_df, aes(x = No.Reference, ymin = lwr, ymax = upr), alpha = 0.2) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black")
  ) +
  annotate("text", x = max(data$No.Reference), y = min(data$coefs)-0.15,
           label = p_value_label,
           hjust = 1, vjust = 0, size = 5, color = "black") +
  annotate("text", x = max(data$No.Reference), y = min(data$coefs) + 0.15,
           label = r_squared_label,
           hjust = 1, vjust = 0, size = 5, color = "black")

build_plot <- ggplot_build(fig_3c)
xlim <- build_plot$layout$panel_scales_x[[1]]$range$range
ylim <- build_plot$layout$panel_scales_y[[1]]$range$range

x_range <- diff(xlim)
y_range <- diff(ylim)
ratio <- x_range / y_range
fig_3c <- fig_3c + coord_fixed(ratio = ratio)

agcgamform <-  as.formula(paste(paste0("KL_Reference_acgam"), paste(paste0("KL_Abundance_acgam"), collapse=" + "), sep=" ~ "))
dmetmeas <- dc[, colnames(dc) %in% c(paste0("KL_Sample_acgam"), paste0("KL_Reference_acgam"), paste0("KL_Abundance_acgam"))]
dmetmeas <- na.omit(dmetmeas)

fitacgam <- lm(agcgamform, data = dmetmeas)

predictions <- predict(fitacgam, newdata = dmetmeas, interval = "confidence")

dmetmeas <- dmetmeas %>%
  mutate(predacgam = predictions[, "fit"],
         lower_conf = predictions[, "lwr"],
         upper_conf = predictions[, "upr"])

p_value <- coef(summary(fitacgam))[2,4]
r_squared <- summary(fitacgam)$r.squared


fig_3d <- ggplot(dmetmeas) +
  geom_point(mapping = aes(x = !!sym("KL_Abundance_acgam"), y = !!sym("KL_Reference_acgam"))) +
  geom_line(mapping = aes(x = !!sym("KL_Abundance_acgam"), y = predacgam), color = "blue", size = 1.2) +
  geom_ribbon(mapping = aes(x = !!sym("KL_Abundance_acgam"), ymin = lower_conf, ymax = upper_conf), 
              fill = "darkgrey", alpha = 0.4) +  # Use "darkgrey" and adjust alpha for better visibility
  xlab("Reference-based redundancy of N-Acetyl-D-glucosamine") +  
  ylab("Abundance-based redundancy of N-Acetyl-D-glucosamine") + 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 12, color = "black"),  # Set axis titles size and color
    axis.text = element_text(size = 12, color = "black")    # Set axis text size and color
  ) +
  labs(title = "N-Acetyl-D-glucosamine") +  # Specify title within labs()
  theme(panel.grid = element_blank()) +
  annotate("text", x = min(dmetmeas$KL_Abundance_acgam), y = min(dmetmeas$KL_Reference_acgam),
           label = paste("p-value =", format(p_value, digits = 2)),
           hjust = 0, vjust = 0, size = 5, color = "black") +
  annotate("text", x = min(dmetmeas$KL_Abundance_acgam), y = min(dmetmeas$KL_Reference_acgam) + 0.17,
           label = paste("R-squared =", format(r_squared, digits = 2)),
           hjust = 0, vjust = 0, size = 5, color = "black")

build_plot <- ggplot_build(fig_3d)
xlim <- build_plot$layout$panel_scales_x[[1]]$range$range
ylim <- build_plot$layout$panel_scales_y[[1]]$range$range

x_range <- diff(xlim)
y_range <- diff(ylim)
ratio <- x_range / y_range
fig_3d <- fig_3d + coord_fixed(ratio = ratio)



#########################################
#### species species interrelations  ####
#########################################
dn <- d[, grep("^KLI", colnames(d))]

dn_long <- gather(dn, key = "met", value = "value")
range_values <- aggregate(value ~ met, data = dn_long, FUN = function(x) diff(range(x)))
dn_long$met <- factor(dn_long$met, levels = range_values$met[order(range_values$value)])

range_values$met[order(range_values$value)]
custom_labels <- gsub("KLI_", "", range_values$met[order(range_values$value)])
custom_labels <- RSQT$Metabolite[match(custom_labels, RSQT$VMHId)]

fig_3e <- ggplot(dn_long, aes(x = met, y = value)) +
  geom_boxplot() +
  scale_x_discrete(name = "Metabolite", labels = custom_labels) +
  scale_y_continuous(name = "Interdependency index") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, color = "black"),  # Adjust x-axis text size and color
    axis.text.y = element_text(size = 12, color = "black"),  # Adjust y-axis text size and color
    axis.title = element_text(size = 12, color = "black"),   # Adjust axis title size and color
  )



fig_3f <- ggplot() +
  geom_point(data = d, aes(x = sNB_lac_L, y = KL_Abundance_lac_L, color = "purple"), alpha = 0.6) +
  geom_point(data = d, aes(x = sNB_but, y = KL_Abundance_but, color = "blue"), alpha = 0.6) +
  geom_point(data = d, aes(x = sNB_succ, y = KL_Abundance_succ, color = "orange"), alpha = 0.6) +
  labs(x = "Total abundance of species that can produce the metabolite", y = "Abundance-based redundancy") +
  scale_color_manual(
    values = c("purple" = "purple", "blue" = "blue", "orange" = "orange"),
    labels = c("purple" = "L-Lactic acid", "blue" = "Butyric acid", "orange" = "Succinic acid"),
    breaks = c("blue", "purple", "orange")
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    legend.position = c(1, 0.5), 
    legend.justification = c(1, 0),
    axis.text.x = element_text(size = 12, color = "black"),  
    axis.text.y = element_text(size = 12, color = "black"),  
    axis.title = element_text(size = 12, color = "black"),   
    legend.key = element_rect(fill = "white"),  # Set legend key background to white
  ) +
  labs(color = "Metabolite")

build_plot <- ggplot_build(fig_3f)
xlim <- build_plot$layout$panel_scales_x[[1]]$range$range
ylim <- build_plot$layout$panel_scales_y[[1]]$range$range

x_range <- diff(xlim)
y_range <- diff(ylim)
ratio <- x_range / y_range
fig_3f <- fig_3f + coord_fixed(ratio = ratio)


combined_plot_fig3 <- plot_grid(
  fig_3a + labs(title = "A)"),
  fig_3b + labs(title = "B)"),
  fig_3c + labs(title = "C)"),
  fig_3d + labs(title = "D)"),
  fig_3e + labs(title = "E)"),
  fig_3f + labs(title = "F)"),
  ncol = 3
)

combined_plot_fig3 <- combined_plot_fig3 + theme(panel.grid = element_blank(),
                                                 plot.background = element_rect(fill = "white", color =NA))

ggsave("Figures/Figure_3.tiff", combined_plot_fig3, width = 18, height = 12, dpi = 300)

##########################################################
#### ### Species diversity <-> functional redundancy  ####
##########################################################

div_red <- data.frame(Metabolite=character(),
                            VMHId=character(),
                            No.Reference=double(),
                            Producing.Species=double(),
                            Total.Sum.Prod=double(),
                            median.Interdependencies=double(),
                            Estimate.Sample.CI=double(),
                            pvalue.KL_Sample=double(),
                            Rsq.Sample=double(),
                            Estimate.Reference.CI=double(),
                            pvalue.KL_Reference=double(),
                            Rsq.Reference=double(),
                            Estimate.Abundance.CI=double(),
                            pvalue.KL_Abundance=double(),
                            Rsq.Abundance=double(),
                            stringsAsFactors=FALSE)

for(i in mets){
  modelSamp <- lm(d[, paste0("KL_Sample_", i)] ~ d[, "shannon"])
  rsquared_samp <- summary(modelSamp)$r.squared
  
  modelRef <- lm(d[, paste0("KL_Reference_", i)] ~ d[, "shannon"])
  rsquared_ref <- summary(modelRef)$r.squared
  
  modelAbu <- lm(d[, paste0("KL_Abundance_", i)] ~ d[, "shannon"])
  rsquared_abu <- summary(modelAbu)$r.squared
  
  modelSampv <- paste0(format(summary(modelSamp)$coef[2, 1], scientific = FALSE, digits = 2)," (",
                format(confint(modelSamp)[2,1], scientific = FALSE, digits = 2), ",", format(confint(modelSamp)[2,2], scientific = FALSE, digits = 2),")")
  
  modelRefv <- paste0(format(summary(modelRef)$coef[2, 1], scientific = FALSE, digits = 2)," (",
                format(confint(modelRef)[2,1], scientific = FALSE, digits = 2), ",", format(confint(modelRef)[2,2], scientific = FALSE, digits = 2),")")
  
  modelAbuv <- paste0(format(summary(modelAbu)$coef[2, 1], scientific = FALSE, digits = 2)," (",
                format(confint(modelAbu)[2,1], scientific = FALSE, digits = 2), ",", format(confint(modelAbu)[2,2], scientific = FALSE, digits = 2),")")
  
  NoReference <- RSQT$NoReference[match(i, RSQT$VMHId)]
  ProSpecies <- RSQT$Producing.Species[match(i, RSQT$VMHId)]
  snB <- mean(d[, paste0("sNB_", i)])
  med <- median(d[, paste0("KLI_", i)], na.rm=TRUE)
  
  
  div_red[nrow(div_red)+1, ] <- c(Mapping$Metabolite[match(i, Mapping$VMHId)], i, NoReference,ProSpecies, snB, med,
                                  modelSampv, summary(modelSamp)$coef[2, 4], rsquared_samp,
                                  modelRefv,summary(modelRef)$coef[2, 4], rsquared_ref,
                                  modelAbuv,summary(modelAbu)$coef[2, 4], rsquared_abu)
}

pvaluesdiv_red <- c(div_red$pvalue.KL_Sample, div_red$pvalue.KL_Reference, div_red$pvalue.KL_Abundance)
pvaluesdiv_red <- as.numeric(pvaluesdiv_red)
Bonfpvaldiv_red <- p.adjust(pvaluesdiv_red, method="bonferroni")

nonsigSample <- div_red$Metabolite[which(div_red$BonfKL_Sample>=0.05)]
nonsigRef <- div_red$Metabolite[which(div_red$BonfKL_Reference>=0.05)]
div_red$Metabolite[div_red$BonfKL_Sample<0.05 & sign(as.numeric(div_red$Estimate.KL_Sample))==1]

div_red$Metabolite[div_red$BonfKL_Abundance<0.05 & sign(as.numeric(div_red$Estimate.KL_Abundance))==1]
div_red$Metabolite[div_red$BonfKL_Abundance<0.05]

div_red$BonfKL_Sample <- Bonfpvaldiv_red[1:45]
div_red$BonfKL_Reference <- Bonfpvaldiv_red[46:90]
div_red$BonfKL_Abundance <- Bonfpvaldiv_red[91:135]

columns_to_numeric <- c(3:6, 9, 11, 12, 14, 15)
div_red[, columns_to_numeric] <- lapply(div_red[, columns_to_numeric], as.numeric)

write.xlsx(div_red, "Tables/TableS4.xlsx", row.names = FALSE)

##########################
## Supplementary plot
div_red_long <- melt(div_red, measure.vars = c("Rsq.Sample", "Rsq.Reference"))
violinplots <- ggplot(div_red_long, aes(x = "", y = value, fill = variable)) +
  geom_violin() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("Rsq.Sample" = "skyblue", "Rsq.Reference" = "lightgreen")) +
  facet_wrap(~ variable, scales = "free", nrow = 1, strip.position = "bottom", labeller = labeller(variable = c("Rsq.Sample" = "Using sample taxon-based functional \n redundancy as response", "Rsq.Reference" = "Using reference taxon-based functional \n redundancy as response"))) +
  ylab("R-squared Value") +
  ylim(0,1) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, color = "black"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        axis.title = element_text(size = 12)) + 
  guides(fill = FALSE) +
  labs(x = NULL)


Rsq_scatterplot <- ggplot(div_red, aes(x = Rsq.Sample, y = Rsq.Reference)) +  
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +  
  geom_text(x = 0.6, y = 0.65, label = "y = x", color = "black", size = 4) +  
  xlab("R-squared using the sample taxon-based measure") +
  ylab("R-squared using the reference taxon-based measure") +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid = element_blank(),
    axis.title = element_text(size = 12)
  )

build_plot <- ggplot_build(Rsq_scatterplot)
xlim <- build_plot$layout$panel_scales_x[[1]]$range$range
ylim <- build_plot$layout$panel_scales_y[[1]]$range$range

x_range <- diff(xlim)
y_range <- diff(ylim)
ratio <- x_range / y_range
Rsq_scatterplot <- Rsq_scatterplot + coord_fixed(ratio = ratio)


# Set white background for individual plots
violinplots <- violinplots + theme(plot.background = element_rect(fill = "white", color = NA))
Rsq_scatterplot <- Rsq_scatterplot + theme(plot.background = element_rect(fill = "white", color = NA))

# Combine the plots with a white background
combined_plot_supp_f1 <- plot_grid(
  violinplots + labs(title = "A)"),
  Rsq_scatterplot + labs(title = "B)"),
  ncol = 2
)

# Save the combined plot with a white background
ggsave("Figures/Supplementary_Figure_1.tiff", combined_plot_supp_f1, width = 18, height = 5, dpi = 300, bg = "white")


create_plot <- function(data, x, y, formula, xlab, ylab, title, y_lim, annotation_shift = 0, annotate_metabolites = FALSE) {
  model <- lm(as.formula(paste(y, "~", x)), data = data)
  p_value <- coef(summary(model))[2,4]
  r_squared <- summary(model)$r.squared
  
  # Base plot
  plot <- ggplot(data, aes_string(x = x, y = y, color = "Total.Sum.Prod")) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", fill = "grey", formula = formula, size = 1) +
    labs(x = xlab, y = ylab, color = "Average total \n species abundance \n that can perform \n the function") +
    scale_color_viridis_c(direction = -1) +
    theme_bw() +
    ylim(y_lim) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 12, color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      legend.position = "none"
    ) +
    # R-squared and p-value annotations
    annotate("text", x = min(data[[x]]) + 0.05 * diff(range(data[[x]])), y = y_lim[2] - annotation_shift, 
             label = paste("R-squared =", format(r_squared, digits = 2)),
             hjust = 0, vjust = 1, size = 4, color = "black") +
    annotate("text", x = min(data[[x]]) + 0.05 * diff(range(data[[x]])), y = y_lim[2] - 0.05 - annotation_shift,
             label = paste("p-value =", format(p_value, digits = 2)),
             hjust = 0, vjust = 1, size = 4, color = "black")
  
  # Add metabolite names for Producing.Species > 60
  if (annotate_metabolites) {
    plot <- plot + geom_text(data = subset(data, Producing.Species > 60), aes(label = Metabolite), 
                             vjust = -1.2, hjust = 1.1, size = 3, color = "black")
  }
  
  build_plot <- ggplot_build(plot)
  xlim <- build_plot$layout$panel_scales_x[[1]]$range$range
  x_range <- diff(xlim)
  ratio <- x_range
  plot
  #plot + coord_fixed(ratio = ratio)
}



# Create plots with adjusted annotation positions and add metabolite names for plot A
# Create plots with adjusted annotation positions and add metabolite names for plot A
prod_spec_plot_samp <- create_plot(div_red, "Producing.Species", "Rsq.Sample", y ~ x, 
                                   "Average producing species", "R-squared", "A)", c(0, 1), annotation_shift = 0, annotate_metabolites = TRUE)
prod_spec_plot_ref <- create_plot(div_red, "Producing.Species", "Rsq.Reference", y ~ x, 
                                  "Average producing species", "R-squared", "C)", c(0, 1), annotation_shift = 0)

prod_spec_plot_abundance <- create_plot(div_red, "Producing.Species", "Rsq.Abundance", y ~ x, 
                                        "Average producing species", "R-squared", "E)", c(0, 1), annotation_shift = 0)
interdependencies_plot_sample <- create_plot(div_red, "median.Interdependencies", "Rsq.Sample", y ~ x, 
                                             "Median interdependency", "R-squared", "B)", c(0, 1), annotation_shift = 0)
interdependencies_plot_ref <- create_plot(div_red, "median.Interdependencies", "Rsq.Reference", y ~ x, 
                                          "Median interdependency", "R-squared", "D)", c(0, 1), annotation_shift = 0)
interdependencies_plot_abundance <- create_plot(div_red, "median.Interdependencies", "Rsq.Abundance", y ~ x, 
                                                "Median interdependency", "R-squared", "F)", c(0, 1), annotation_shift = 0)


# Create text elements for titles
top_text <- ggdraw() + 
  draw_label("Sample taxon-based functional redundancy", fontface = 'bold', size = 15, hjust = 0.5) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

middle_text <- ggdraw() + 
  draw_label("Reference taxon-based functional redundancy", fontface = 'bold', size = 15, hjust = 0.5) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

bottom_text <- ggdraw() + 
  draw_label("Abundance-based functional redundancy", fontface = 'bold', size = 15, hjust = 0.5) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

# Combine the plots into grids
top_plots <- plot_grid(prod_spec_plot_samp + labs(title = "A)"), interdependencies_plot_sample + labs(title = "B)"), ncol = 2)
middle_plots <- plot_grid(prod_spec_plot_ref + labs(title = "C)"), interdependencies_plot_ref + labs(title = "D)"), ncol = 2)
bottom_plots <- plot_grid(prod_spec_plot_abundance + labs(title = "E)"), interdependencies_plot_abundance + labs(title = "F)"), ncol = 2)

# Define a plot just for extracting the legend
legend_plot <- ggplot(div_red, aes(x = Producing.Species, y = Rsq.Sample, color = Total.Sum.Prod)) +
  geom_point(size = 3) +
  scale_color_viridis_c(direction = -1) +
  theme_void() +
  guides(color = guide_colorbar(
    title = "Average total species abundance that can perform the function",
    title.position = "left", 
    title.hjust = 0.5,
    title.vjust = 0.3,  # Adjust the title position
    label.position = "bottom",
    bar.width = 1,
    bar.height = 15,
    title.theme = element_text(vjust = 0.5),  # Adjust title text vertical position
    label.theme = element_text(vjust = 0.5),  # Adjust label text vertical position
    bar.background = element_rect(fill = "white", color = "black")  # White background for color bar with black frame
  )) +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = "black"),  # White background with black frame
    legend.box.background = element_rect(fill = "white", color = "black"), # White background with black frame
    legend.margin = margin(t = 10, r = 10, b = 10, l = 10),  # Add space around the legend
    legend.box.spacing = unit(10, "pt"),  # Add space between legend and the plot area
    legend.title = element_text(vjust = 0.5),  # Lower title text
    legend.text = element_text(vjust = 0.5)   # Lower label text
  )

# Extract the legend as a grob
legend_grob <- ggplotGrob(legend_plot)$grobs[[which(sapply(ggplotGrob(legend_plot)$grobs, function(x) x$name) == "guide-box")]]

# Combine everything with the legend at the bottom
Figure_4 <- plot_grid(
  top_text, top_plots, 
  middle_text, middle_plots, 
  bottom_text, bottom_plots, 
  legend_grob,
  ncol = 1, 
  rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1, 0.2)  # Adjust height for the legend
)

# Print the final figure
print(Figure_4)

# Save the figure
ggsave("Figures/Figure_4.tiff", Figure_4, width = 10, height = 15, dpi = 300, bg = "white")

#########################################################
#### faecal concentraions <-> functional redundancy  ####
#########################################################
# Initialize the list to store valid metabolites
metsana <- c()

# Check each metabolite and filter based on the log transformation condition
for (i in mets) {
  j <- ifelse(grepl("^[[:digit:]]+", i), paste0("X", i), i)
  if (!j %in% names(d)) next  # Skip if 'j' is not a column in 'd'
  
  logvi <- log(d[[j]])
  logvi[is.infinite(logvi)] <- NA
  if (sum(is.na(logvi)) < 0.5 * length(logvi)) {
    metsana <- c(metsana, j)
  }
}

# Remove the prefixed 'X' from valid metabolites
metsana <- gsub("^X", "", metsana)

# Initialize the parameter table
paramtab_faec <- data.frame(
  Metabolite = character(),
  VMHId = character(),
  Estimate.Sample.CI = character(),
  pvalue.KL_Sample = double(),
  Estimate.Reference.CI = character(),
  pvalue.KL_Reference = double(),
  Estimate.Abundance.CI = character(),
  pvalue.KL_Abundance = double(),
  stringsAsFactors = FALSE
)

# Populate the parameter table
for (i in metsana) {
  j <- ifelse(grepl("^[[:digit:]]+", i), paste0("X", i), i)
  
  logvi <- log(d[[j]])
  logvi[is.infinite(logvi)] <- NA
  
  modelKL_Sample <- lm(logvi ~ d[[paste0("KL_Sample_", i)]] + d$Age + d$BMI + d$Gender + d$healthstat)
  modelKL_Reference <- lm(logvi ~ d[[paste0("KL_Reference_", i)]] + d$Age + d$BMI + d$Gender + d$healthstat)
  modelKL_Abundance <- lm(logvi ~ d[[paste0("KL_Abundance_", i)]] + d$Age + d$BMI + d$Gender + d$healthstat)
  
  modelSampf <- sprintf("%.5f(%.5f,%.5f)", coef(summary(modelKL_Sample))[2, 1], confint(modelKL_Sample)[2, 1], confint(modelKL_Sample)[2, 2])
  modelReff <- sprintf("%.5f(%.5f,%.5f)", coef(summary(modelKL_Reference))[2, 1], confint(modelKL_Reference)[2, 1], confint(modelKL_Reference)[2, 2])
  modelAbuf <- sprintf("%.5f(%.5f,%.5f)", coef(summary(modelKL_Abundance))[2, 1], confint(modelKL_Abundance)[2, 1], confint(modelKL_Abundance)[2, 2])
  
  paramtab_faec <- rbind(paramtab_faec, data.frame(
    Metabolite = Mapping$Metabolite[match(i, Mapping$VMHId)],
    VMHId = i,
    Estimate.Sample.CI = modelSampf,
    pvalue.KL_Sample = summary(modelKL_Sample)$coef[2, 4],
    Estimate.Reference.CI = modelReff,
    pvalue.KL_Reference = summary(modelKL_Reference)$coef[2, 4],
    Estimate.Abundance.CI = modelAbuf,
    pvalue.KL_Abundance = summary(modelKL_Abundance)$coef[2, 4],
    stringsAsFactors = FALSE
  ))
}

# Adjust p-values for multiple testing
pvalues <- c(paramtab_faec$pvalue.KL_Sample, paramtab_faec$pvalue.KL_Reference, paramtab_faec$pvalue.KL_Abundance)
Bonfpval <- p.adjust(pvalues, method = "bonferroni")

paramtab_faec$BonfKL_Sample <- Bonfpval[1:nrow(paramtab_faec)]
paramtab_faec$BonfKL_Reference <- Bonfpval[(nrow(paramtab_faec) + 1):(2 * nrow(paramtab_faec))]
paramtab_faec$BonfKL_Abundance <- Bonfpval[(2 * nrow(paramtab_faec) + 1):(3 * nrow(paramtab_faec))]

# Filter significant metabolites
KL_Samplesig <- paramtab_faec$VMHId[paramtab_faec$BonfKL_Sample < 0.05]
KL_Referencesig <- paramtab_faec$VMHId[paramtab_faec$BonfKL_Reference < 0.05]
KL_Abundancesig <- paramtab_faec$VMHId[paramtab_faec$BonfKL_Abundance < 0.05]

# Check the structure of 'dt'
## Plot
dt <- d

# Prepare columns for plotting
cols <- unique(c(KL_Samplesig, KL_Referencesig, KL_Abundancesig))
cols <- gsub("^(\\d+)", "X\\1", cols)  # Ensure variables starting with digits are prefixed with "X"

# Log-transform the relevant columns
dt[cols] <- lapply(dt[cols], function(x) {
  logx <- log(x)
  logx[is.infinite(logx)] <- NA
  return(logx)
})

plot_faecal_concentration <- function(data, x, y, xlab, ylab, title) {
  data[[y]][is.infinite(data[[y]])] <- NA
  filtered_data <- data[!is.na(data[[x]]) & !is.na(data[[y]]), ]
  model <- lm(reformulate(x, y), data = filtered_data)
  p_value <- coef(summary(model))[2, 4]
  r_squared <- summary(model)$r.squared
  
  # Build the initial plot
  plot <- ggplot(filtered_data, aes_string(x = x, y = y)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", color = "blue", se = TRUE, size = 1.2) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title) +
    theme_bw() +
    annotate("text", x = min(filtered_data[[x]]), y = min(filtered_data[[y]]) + 0.4,
             label = paste("R-squared =", format(r_squared, digits = 2)),
             hjust = 0, vjust = 1, size = 4, color = "black") +
    annotate("text", x = min(filtered_data[[x]]), y = min(filtered_data[[y]]) + 0.04,
             label = paste("p-value =", format(p_value, digits = 2)),
             hjust = 0, vjust = 1, size = 4, color = "black") +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 12, color = "black"),  
          axis.text.y = element_text(size = 12, color = "black"),  
          axis.title = element_text(size = 12, color = "black"))
  
  # Build the plot to get the x and y range
  build_plot <- ggplot_build(plot)
  xlim <- build_plot$layout$panel_scales_x[[1]]$range$range
  ylim <- build_plot$layout$panel_scales_y[[1]]$range$range
  x_range <- diff(xlim)
  y_range <- diff(ylim)
  ratio <- x_range / y_range
  
  # Add coord_fixed() with the computed ratio
  plot + coord_fixed(ratio = ratio)
}


# Create plots
f5_a <- plot_faecal_concentration(dt, "KL_Reference_4abut", "X4abut", 
                                  "Reference taxon-based redundancy of gamma-Aminobutyric acid",
                                  "Log faecal concentration of gamma-Aminobutyric acid [nmol/g]", "A)")

f5_b <- plot_faecal_concentration(dt, "KL_Reference_isoval", "isoval", 
                                  "Reference taxon-based redundancy of isovaleric acid",
                                  "Log faecal concentration of isovaleric acid [nmol/g]", "B)")

f5_c <- plot_faecal_concentration(dt, "KL_Reference_taur", "taur", 
                                  "Reference taxon-based redundancy of taurine",
                                  "Log faecal concentration of taurine [nmol/g]", "C)")
# D) by hand
x <- "but"
y <- "KL_Abundance_but"
xlab <- "Log faecal concentration of butyric acid [nmol/g]"
ylab <- "Abundance-based redundancy of butyric acid"
title <- "D)"

# Filter out infinite values
dt[[y]][is.infinite(dt[[y]])] <- NA
filtered_data <- dt[!is.na(dt[[x]]) & !is.na(dt[[y]]), ]

# Create the model
model <- lm(reformulate(x, y), data = filtered_data)
p_value <- coef(summary(model))[2, 4]
r_squared <- summary(model)$r.squared

# Build the plot
f5_d <- ggplot(filtered_data, aes_string(x = x, y = y)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", color = "blue", se = TRUE, size = 1.2) +
  xlab(xlab) +
  ylab(ylab) +
  ggtitle(title) +
  theme_bw() +
  annotate("text", x = min(filtered_data[[x]]), y = min(filtered_data[[y]]) + 0.3,  # Adjust y to move text up
           label = paste("R-squared =", format(r_squared, digits = 2)),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  annotate("text", x = min(filtered_data[[x]]), y = min(filtered_data[[y]])+0.05,  # Adjust y to move text up
           label = paste("p-value =", format.pval(p_value, digits = 2)),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"),  
        axis.text.y = element_text(size = 12, color = "black"),  
        axis.title = element_text(size = 12, color = "black"))

# Build the plot to get the x and y range
build_plot <- ggplot_build(f5_d)
xlim <- build_plot$layout$panel_scales_x[[1]]$range$range
ylim <- build_plot$layout$panel_scales_y[[1]]$range$range
x_range <- diff(xlim)
y_range <- diff(ylim)
ratio <- x_range / y_range

# Add coord_fixed() with the computed ratio
f5_d <- f5_d + coord_fixed(ratio = ratio)

f5_e <- plot_faecal_concentration(dt, "KL_Sample_15dap", "X15dap", 
                                  "Sample taxon-based redundancy",
                                  "Log faecal concentration of 1,5-Diaminopentane [nmol/g]", "E)")

f5_f <- plot_faecal_concentration(dt, "KL_Reference_15dap", "X15dap", 
                                  "Reference taxon-based redundancy",
                                  "Log faecal concentration of 1,5-Diaminopentane [nmol/g]", "F)")

f5_g <- plot_faecal_concentration(dt, "KL_Abundance_15dap", "X15dap", 
                                  "Abundance-based redundancy",
                                  "Log faecal concentration of 1,5-Diaminopentane [nmol/g]", "F)")

# Define the top row with 4 plots
top_row <- plot_grid(
  f5_a + labs(title = "A)"),
  f5_b + labs(title = "B)"),
  f5_c + labs(title = "C)"),
  f5_d + labs(title = "D)"),
  ncol = 4,
  align = "h",  # Horizontal alignment
  rel_widths = c(1, 1, 1, 1)  # Equal width for all plots
)


empty_plot <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "white", color = NA))


bottom_row <- plot_grid(
  empty_plot,
  f5_e + labs(title = "E)"),
  f5_f + labs(title = "F)"),
  f5_g + labs(title = "G)"),
  empty_plot,
  ncol = 5,
  align = "h",  # Horizontal alignment
  rel_widths = c(0.5,1, 1, 1,0.5)  # Equal width for all plots
)

# Define the bottom title
bottom_title <- ggdraw() + 
  draw_label("1,5-Diaminopentane (cadaverine)", 
             fontface = "bold", 
             size = 20, 
             color = "black", 
             hjust = 0.5,  # Center horizontally
             vjust = 1)   # Position title at the top of the `ggdraw` area




# Combine the bottom title and the bottom row of plots
bottom_combined <- plot_grid(
  bottom_title,  # Title above the bottom plots
  bottom_row,    # The bottom row with plots
  ncol = 1,
  rel_heights = c(0.1, 1)  # Height of title vs. the row with plots
)

# Combine the top row with the bottom title and bottom row of plots
combined_plot_f5 <- plot_grid(
  top_row,        # The top plots
  bottom_combined,  # The bottom title and plots
  ncol = 1,
  rel_heights = c(1, 1.1)  # Equal heights for the top and bottom sections
)

combined_plot_f5 <- combined_plot_f5 + theme(plot.background = element_rect(fill = "white", color = NA))

print(combined_plot_f5)

ggsave("Figures/Figure_5.tiff", combined_plot_f5, width = 24, height = 12, dpi = 300)
length(which(paramtab_faec$pvalue.KL_Sample < 0.05))
length(which(paramtab_faec$pvalue.KL_Reference < 0.05))
length(which(paramtab_faec$pvalue.KL_Abundance < 0.05))


# Round p-values to export
columns_to_numeric<- c(4, 6, 8:10)
paramtab_faec[, columns_to_numeric] <- lapply(paramtab_faec[, columns_to_numeric], as.numeric)
write.xlsx(paramtab_faec, "Tables/TableS5.xlsx", row.names = FALSE)

#################################################
#### Functional redundancy <-> Healthstatus  ####
#################################################
paramtabHealthstat <- data.frame(Metabolite=character(),
                                 VMHId=character(),
                                 Estimate.Sample.CI=double(),
                                 pvalue.KL_Sample=double(),
                                 Estimate.Reference.CI=double(),
                                 pvalue.KL_Reference=double(),
                                 Estimate.Abundance.CI=double(),
                                 pvalue.KL_Abundance=double(),
                                 Estimate.Interdependency.CI=double(),
                                 pvalue.Interdependency=double(),
                                 stringsAsFactors=FALSE)

for(i in mets){
  modelSamp <- lm(d[, paste0("KL_Sample_", i)] ~ d[, "healthstat"] + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
  modelRef <- lm(d[, paste0("KL_Reference_", i)] ~ d[, "healthstat"]  + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
  modelAbu <- lm(d[, paste0("KL_Abundance_", i)] ~ d[, "healthstat"] + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
  logvi <- log(d[, paste0("KLI_", i)])
  logvi[which(logvi==-Inf)] <- NaN
  modelI <- lm(logvi ~ d[, "healthstat"] + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
  
  modelSampv <- paste0(format(summary(modelSamp)$coef[2, 1], scientific = FALSE, digits = 5),"(",
                       format(confint(modelSamp)[2,1], scientific = FALSE, digits = 5), ",", format(confint(modelSamp)[2,2], scientific = FALSE, digits = 5),")")
  
  modelRefv <- paste0(format(summary(modelRef)$coef[2, 1], scientific = FALSE, digits = 5),"(",
                      format(confint(modelRef)[2,1], scientific = FALSE, digits = 5), ",", format(confint(modelRef)[2,2], scientific = FALSE, digits = 5),")")
  
  modelAbuv <- paste0(format(summary(modelAbu)$coef[2, 1], scientific = FALSE, digits = 5),"(",
                      format(confint(modelAbu)[2,1], scientific = FALSE, digits = 5), ",", format(confint(modelAbu)[2,2], scientific = FALSE, digits = 5),")")
  
  modelIv <- paste0(format(summary(modelI)$coef[2, 1], scientific = FALSE, digits = 5),"(",
                    format(confint(modelI)[2,1], scientific = FALSE, digits = 5), ",", format(confint(modelI)[2,2], scientific = FALSE, digits = 5),")")
  
  paramtabHealthstat[nrow(paramtabHealthstat)+1, ] <- c(Mapping$Metabolite[match(i, Mapping$VMHId)], i, modelSampv, summary(modelSamp)$coef[2, 4],
                                                        modelRefv,summary(modelRef)$coef[2, 4], modelAbuv,summary(modelAbu)$coef[2, 4],
                                                        modelIv, summary(modelI)$coef[2,4])
}

pvalues <- c(paramtabHealthstat$pvalue.KL_Sample, paramtabHealthstat$pvalue.KL_Reference, paramtabHealthstat$pvalue.KL_Abundance, paramtabHealthstat$pvalue.Interdependency)
pvalues <- as.numeric(pvalues)
Bonfpval<- p.adjust(pvalues, method="bonferroni")


paramtabHealthstat$BonfKL_Sample <- Bonfpval[1:45]
paramtabHealthstat$BonfKL_Reference <- Bonfpval[46:90]
paramtabHealthstat$BonfKL_Abundance <- Bonfpval[91:135]
paramtabHealthstat$BonfKL_I <- Bonfpval[136:180]


modelI_g_spec <- lm(d[, "spec"] ~ d[, "healthstat"] + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
summary(modelI_g_spec)
confint(modelI_g_spec)


modelI_g_shannon <- lm(d[, "shannon"] ~ d[, "healthstat"] + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
summary(modelI_g_shannon)
confint(modelI_g_shannon)

columns_to_numeric <- c(4,6,8)

paramtabHealthstat[, columns_to_numeric] <- lapply(paramtabHealthstat[, columns_to_numeric], as.numeric)
write.xlsx(paramtabHealthstat, "Tables/TableS8.xlsx")

## Summary stats of pyruvate
modelRef <- lm(d[, paste0("KL_Reference_", "pyr")] ~ d[, "healthstat"]  + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
coef(summary(modelRef))[2,4]
confint(modelRef)

## Global redundancy measure
mets_I <- gsub("^", "KLI_", mets)

d$median_mets_I <- log(apply(d[, mets_I], 1, median))

modelI_g_sh <- lm(d[, "median_mets_I"] ~ d[, "healthstat"] + d[, "Age"] + d[, "BMI"] + d[, "Gender"])
summary(modelI_g_sh)
confint(modelI_g_sh)
