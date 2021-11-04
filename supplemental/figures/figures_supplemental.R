# Supplemental Materials
# Add figures comparing ISP to Fixed (code for these sims in 'FixedvsISPsims.Rmd')
# Add Google Drive > Independent Segments Procedure > Simulations > Figure 1 (Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing). Make ISP solid, fixed dash. Code in 'Figures.R'
# Add Google Drive > Independent Segments Procedure > Simulations > Figure 4 (Bias in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing). Code in 'Figures.Rmd'
# Add Google Drive > Independent Segments Procedure > Simulations > Figure 7 (Bias-Variance in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing).  Code in 'Figures.Rmd'
# Add Google Drive > Independent Segments Procedure > Simulations > Figure 6 (Estimated Cohen's d using Fixed versus Independent Segments Hypothesis Testing).  Code in 'Figures Appendix.Rmd'


# Set up ------------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(ggtext)
library(grid)

df = read.csv("supplemental/simulations/data_supplemental.csv")

# Figure S1 ----------------------------------------------------------------
# Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing

plt = function(f){
  D = mean(f$d) #get separate plot for each combo of power and d
  target_power = mean(f$power)
  t = sprintf("\u03b4 = %0.1f, 1 - *\u03b2* = %0.1f", D, target_power)
  
  ggplot(data = f, aes(x = ES, group = proc, linetype = proc)) +
    geom_density() +
    theme_classic() +
    theme(
      legend.position = "none", 
      plot.title = ggtext::element_markdown(size = 11, hjust = 0.5, face = "bold"),
      axis.title = element_blank(), 
      text = element_text(size = 11), 
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1), 
      axis.ticks.y = element_blank()
      ) +
    scale_x_continuous(breaks = c(0, 0.5*D, D, 1.5*D),
                       limits = c(-D, 2.5*D)) +
    scale_linetype_manual(values=c("dashed", "solid"))+
    labs(title = t)
}

p = lapply(split(df, list(df$power, df$d_actual), drop = TRUE), plt)

p[[1]] = p[[1]] + #manually set legend in first facet
  annotate("segment", x = -0.2, y = 4.5, xend = -.1, yend = 4.5) +
  annotate("segment", x = -0.2, y = 3.5, xend = -0.1, yend = 3.5, linetype = "dashed") +
  annotate(geom = "text", x = -.08, y = 4.5, label = "ISP", hjust = 0) +
  annotate(geom = "text", x = -.08, y = 3.5, label = "Fixed", hjust = 0) 

#tiff(file="Figure1.tiff",width=2100,height=2000, units = "px", res = 300)
Fig1 <- grid.arrange(grobs = p, 
                     ncol = 3, 
                     top = textGrob("Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing 
                     \n Number of simulations = 100,000", 
                                    gp = gpar(fontsize = 12)),
                     bottom = textGrob("Observed Cohen's d", 
                                       gp = gpar(fontsize = 12)),
                     left = textGrob("Density", 
                                     gp=gpar(fontsize = 12), 
                                     rot = 90, 
                                     just = "right"))
#dev.off()

# Figure S2 ----------------------------------------------------------------
# Bias in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing


# Figure S3 ---------------------------------------------------------------
# Bias-Variance in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing


# Figure S4 ---------------------------------------------------------------
# Exaggeration of Significant Results using Fixed versus Independent Segments Hypothesis Testing




