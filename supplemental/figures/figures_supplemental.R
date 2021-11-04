# Supplemental Materials
# Add figures comparing ISP to Fixed (code for these sims in 'FixedvsISPsims.Rmd')
# Add Google Drive > Independent Segments Procedure > Simulations > Figure 7 (Bias-Variance in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing).  Code in 'Figures.Rmd'
# Add Google Drive > Independent Segments Procedure > Simulations > Figure 6 (Estimated Cohen's d using Fixed versus Independent Segments Hypothesis Testing).  Code in 'Figures Appendix.Rmd'

# Set up ------------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(ggtext)
library(grid)
library(gridExtra)

zipped_df = gzfile("supplemental/simulations/data_supplemental.csv.gz", "rt")
df = read.csv(zipped_df, header = T) 

# Figure S1 ----------------------------------------------------------------
# Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing
# Expected effect sizes only (i.e., d powered for == actual d)
# Google Drive > Independent Segments Procedure > Simulations > Figure 1 (Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing). Code in 'Figures.R'

df_expected <- df %>% 
  filter(d_forpower == d_actual)

plt = function(f){
  D = mean(f$d_actual) #get separate plot for each combo of power and d
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

p = lapply(split(df_expected, list(df_expected$power, df_expected$d_actual), drop = TRUE), plt)

p[[1]] = p[[1]] + #manually set legend in first facet
  annotate("segment", x = -0.2, y = 4.5, xend = -.1, yend = 4.5) +
  annotate("segment", x = -0.2, y = 3.5, xend = -.1, yend = 3.5, linetype = "dashed") +
  annotate(geom = "text", x = -.08, y = 4.3, label = "ISP", hjust = 0) +
  annotate(geom = "text", x = -.08, y = 3.5, label = "Fixed", hjust = 0) 

tiff(file="supplemental/figures/figureS1.tiff",width=2100,height=2700, units = "px", res = 300)
figS1 <- grid.arrange(
  grobs = p, 
  ncol = 3, 
  top = textGrob(
    "Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing", 
    gp = gpar(fontsize = 12)
    ),
  bottom = textGrob(
    "Observed Cohen's d", 
    gp = gpar(fontsize = 12)
    ),
  left = textGrob(
    "Density", 
    gp=gpar(fontsize = 12), 
    rot = 90, 
    just = "right")
  )
dev.off()

# Figure S2 ----------------------------------------------------------------
# Bias in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing
# Includes unexpected effect sizes
# Focus is only on power = 0.8
# Add Google Drive > Independent Segments Procedure > Simulations > Figure 4 (Bias in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing). Code in 'Figures.Rmd'

df_unexpected <- df %>% 
  filter(power == 0.8)

df_summary = df_unexpected %>% 
  group_by(d_forpower, proc, d_actual) %>%  
  summarize(ES = mean(ES)) %>% 
  mutate(bias = ES - d_actual)

bplot = function(f) {
  ggplot(data = f, 
         mapping = aes(x = d_actual, 
                       y = bias, 
                       group = proc, 
                       linetype = proc)) +
    geom_smooth(color = "black", se = F, size = .6, method = "gam", formula = y ~ s(x)) + 
    #there is *some* bias in the fixed procedure, but it's imperceptibly small, especially when using a smoother  
    theme_classic() +
    #geom_smooth(color = "black", 
     #           se = F, 
      #          size = .6, 
       #         method = "gam", 
        #        formula = y ~ s(x, k = 11)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       labels = paste0("\u03B4 = ", seq(0, 1, 0.2))) +
    scale_y_continuous(breaks = c(-max(abs(f$bias)), 0, max(abs(f$bias))),
                       labels = c(round(-max(abs(f$bias)), 2), 0, round(max(abs(f$bias)), 2)),
                       limits = c(-max(abs(f$bias)), max(abs(f$bias)))) +
    labs(title = paste0("H\u00b2: \u03B4 = ", mean(f$d_forpower), ", 1 - *\u03b2* = 0.8")) +
    theme(
      axis.title = element_blank(),
      plot.title = ggtext::element_markdown(hjust = 0.5),
      axis.text.x = ggtext::element_markdown(size = 10),
      axis.text.y = ggtext::element_markdown(size = 12),
      legend.position = "none"
      )
}

bp = lapply(split(df_summary, df_summary$d_forpower), bplot) #create 4 plots, one for each d_forpower
bp[[1]] = bp[[1]]  + annotate("segment", x = 0.6, y = .02, xend = 0.7, yend = .02, linetype = "dashed") +
  annotate("segment", x = 0.6, y = .01, xend = 0.7, yend = .01) +
  annotate(geom = "text", x = .7, y = .02, label = "ISP", hjust = 0) +
  annotate(geom = "text", x = .7, y = .01, label = "Fixed", hjust = 0) #manually set legend in first plot

tiff(file="supplemental/figures/figureS2.tiff",
     width=3000,
     height=3000, 
     units = "px", 
     res = 300)
figS2 <- grid.arrange(
  grobs = bp, 
  ncol=2, 
  top = textGrob("Bias in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing", 
                 gp = gpar(fontsize = 12)),
  bottom = textGrob("True Effect Size (Cohen's d)", 
                    gp = gpar(fontsize = 12)),
  left = textGrob("Average Estimated Effect Size - True Effect Size", 
                  gp=gpar(fontsize = 12), 
                  rot = 90, 
                  just = "center"))

dev.off()

quartz(type = 'pdf', file = 'Figure4.pdf', width=15)
grid.arrange(grobs = bp, ncol=4, top = textGrob("Bias in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing \n Number of Simulations = 10,000", gp = gpar(fontsize = 12)),
             bottom = textGrob("True Effect Size (Cohen's d)", gp = gpar(fontsize = 12)),
             left = textGrob("Average Estimated Effect Size - True Effect Size", gp=gpar(fontsize = 12), rot = 90, just = "center"))



grid.arrange(grobs = bp, ncol =4, top = textGrob("Bias in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing \n Number of Simulations = 10,000", gp = gpar(fontsize = 12)),
             bottom = textGrob("True Effect Size (Cohen's d)", gp = gpar(fontsize = 12)),
             left = textGrob("Average Significant Effect - True Effect Size", gp=gpar(fontsize = 12), rot = 90, just = "center"))



# Figure S3 ---------------------------------------------------------------
# Bias-Variance in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing


# Figure S4 ---------------------------------------------------------------
# Exaggeration of Significant Results using Fixed versus Independent Segments Hypothesis Testing




