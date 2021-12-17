# Plots for iNovo paper

# Read in packages

library(tidyverse)
library(cowplot)
library(reshape2)

# Load data
cometabolism1 <- read_csv("/Users/Alex/Desktop/iNovo/Model_results/SA_VA_pHBA_dFBA_results.csv")
cometabolism2 <- read_csv("/Users/Alex/Desktop/iNovo/Model_results/glu_SA_VA_pHBA_dFBA_results.csv")
glucose_ratios <- read_csv("/Users/Alex/Desktop/iNovo/Model_results/aromatic_glucose_ratios.csv")

####### Co-metabolism - dFBA results
# Panel 1 - no glucose
# Select columns to plot
to_plot <- c("exSA", "exVA", "expHBA")

# Add Time columns for hours and days
cometabolism1$Hours <- cometabolism1$Time / 60
cometabolism1$Days <- cometabolism1$Hours / 24
cometabolism1 <- cometabolism1[which(cometabolism1$Days < 1.2), ]

# Pull out and reshape the variables for plotting
reshaped1 <- cometabolism1[, c("Time", to_plot)]
reshaped1 <- melt(reshaped1, id = "Time")
reshaped1$Hours <- reshaped1$Time / 60
reshaped1$Days <- reshaped1$Hours / 24

p1 <- ggplot(reshaped1, aes(x = Days, y = value, color = variable)) + geom_line(size = 1) + labs(x = NULL, y = "mmol/L") + scale_color_manual(values = c("gold2", "darkturquoise", "olivedrab3")) + scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01)) + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none")

# add another panel of biomass
biomass1 <- cometabolism1[, c("Time", "Biomass")]
biomass1 <- melt(biomass1, id = "Time")
biomass1$Hours <- biomass1$Time / 60
biomass1$Days <- biomass1$Hours / 24

b1 <- ggplot(biomass1, aes(x = Days, y = value)) + labs(x = "Days", y = "Biomass (g)") + geom_line(size = 1) + scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 0.45)) + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none")


# Panel 2 - with glucose
to_plot <- c("exSA", "exVA", "expHBA", "exC00031")
cometabolism2$Hours <- cometabolism2$Time / 60
cometabolism2$Days <- cometabolism2$Hours / 24
cometabolism2 <- cometabolism2[which(cometabolism2$Days < 1.2), ]

reshaped2 <- cometabolism2[, c("Time", to_plot)]
reshaped2 <- melt(reshaped2, id = "Time")
reshaped2$Hours <- reshaped2$Time / 60
reshaped2$Days <- reshaped2$Hours / 24

p2 <- ggplot(reshaped2, aes(x = Days, y = value, color = variable)) + geom_line(size = 1) + labs(x = NULL, y = "mmol/L", color = "Substrate") + scale_color_manual(values = c("gold2", "darkturquoise", "olivedrab3", "hotpink3"), labels = c("Syringic acid", "Vanillic acid", "p-HBA", "Glucose")) + theme_bw() + theme(panel.grid.minor = element_blank())

biomass2 <- cometabolism2[, c("Time", "Biomass")]
biomass2 <- melt(biomass2, id = "Time")
biomass2$Hours <- biomass2$Time / 60
biomass2$Days <- biomass2$Hours / 24

b2 <- ggplot(biomass2, aes(x = Days, y = value)) + labs(x = "Days", y = "Biomass (g)") + geom_line(size = 1) + scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 0.45)) + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none")


legend <- get_legend(p2)
p2 <- ggplot(reshaped2, aes(x = Days, y = value, color = variable)) + geom_line(size = 1) + labs(x = NULL, y = "mmol/L") + scale_color_manual(values = c("gold2", "darkturquoise", "olivedrab3", "hotpink3")) + scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01)) + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none")

cometabolism_plot <- plot_grid(p1, p2, legend, b1, b2, nrow = 2, ncol = 3, labels = c("A", "B", "", "C", "D"), rel_widths = c(1, 1, 0.5))

save_plot("/Users/Alex/Desktop/iNovo/Plots_and_Tables/cometabolism_plot.pdf", cometabolism_plot)


####### Plot aromatic:glucose ratios
# Pull out needed data and convert to long format

colnames(glucose_ratios) <- c("Substrate", "Ratio", "Rate")
glucose_ratios$Ratio <- factor(glucose_ratios$Ratio, levels = c("1:4", "2:3", "1:1", "3:2", "4:1", "9:1"))
glucose_ratios$Substrate <- factor(glucose_ratios$Substrate, levels = c("VA", "p-HBA", "SA"))

ratios_plot <- ggplot(glucose_ratios, aes(x = Ratio, y = Rate, fill = Substrate)) + geom_col(position = "dodge", color = "black") + scale_fill_manual(values = c("darkturquoise", "olivedrab3", "gold"), labels = c("Vanillic acid", "p-HBA", "Syringic acid")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(expand = c(0, 0), limits = c(0, 0.0225)) + labs(x = "Aromatic:Glucose", y = "g/L/hr PDC produced")

save_plot("/Users/Alex/Desktop/iNovo/Plots_and_Tables/ratios_plot.pdf", ratios_plot)
