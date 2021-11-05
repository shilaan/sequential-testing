# Supplemental Materials

# Set up ------------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(ggtext)
library(grid)
library(gridExtra)


# Read data ---------------------------------------------------------------
zipped_df = gzfile("supplemental/simulations/data_supplemental.csv.gz", "rt")
df = read.csv(zipped_df, header = T) 

df_summary = df %>% 
  group_by(d_forpower, power, proc, d_actual) %>%  
  summarize(mean_ES = mean(ES),
            mse = mean((ES - d_actual)^2),
            var = mean((ES - mean(ES))^2)) %>% 
  mutate(bias = mean_ES - d_actual,
         bias2 = (mean_ES - d_actual)^2) 

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
  annotate(geom = "text", x = -.08, y = 4.5, label = "ISP", hjust = 0) +
  annotate(geom = "text", x = -.08, y = 3.5, label = "Fixed", hjust = 0) 

#tiff(file="supplemental/figures/figureS1.tiff",width=2100,height=2700, units = "px", res = 300)
quartz(type = "pdf", file = "supplemental/figures/figureS1.pdf", width=9, height = 8)
grid.arrange(
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

bplot = function(f) {
  ggplot(data = f, 
         mapping = aes(x = d_actual, 
                       y = bias, 
                       group = proc, 
                       linetype = proc)) +
    geom_smooth(color = "black", 
                se = F, 
                size = .6, 
                method = "gam", 
                formula = y ~ s(x, k = 8)) + 
    theme_classic() +
    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       labels = paste0("\u03B4 = ", seq(0, 1, 0.2))) +
    scale_y_continuous(breaks = c(-max(abs(f$bias)), 0, max(abs(f$bias))),
                       labels = c(round(-max(abs(f$bias)), 2), 0, round(max(abs(f$bias)), 2)),
                       limits = c(-max(abs(f$bias)), max(abs(f$bias)))) +
    scale_linetype_manual(values=c("dashed", "solid")) + 
    labs(title = paste0("H\u00b2: \u03B4 = ", 
                        mean(f$d_forpower), 
                        ", 1 - *\u03b2* = ", 
                        mean(f$power))) +
    theme(
      axis.title = element_blank(),
      plot.title = ggtext::element_markdown(hjust = 0.5),
      axis.text.x = ggtext::element_markdown(size = 10),
      axis.text.y = ggtext::element_markdown(size = 12),
      legend.position = "none"
      )
}

bp = lapply(split(df_summary, list(df_summary$d_forpower, df_summary$power)), bplot) #create 4 plots, one for each d_forpower
bp[[1]] = bp[[1]]  + 
  annotate("segment", x = 0.6, y = .05, xend = 0.75, yend = .05) +
  annotate("segment", x = 0.6, y = .035, xend = 0.75, yend = .035, linetype = "dashed") +
  annotate(geom = "text", x = .8, y = .05, label = "ISP", hjust = 0) +
  annotate(geom = "text", x = .8, y = .035, label = "Fixed", hjust = 0) #manually set legend in first plot

#tiff(file="supplemental/figures/figureS2.tiff", width=4000, height=3000, units = "px", res = 300)
quartz(type = "pdf", file = "supplemental/figures/figureS2.pdf", width=13, height = 9)
grid.arrange(
  grobs = bp, 
  ncol = 4, 
  top = textGrob("Bias in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing", 
             gp = gpar(fontsize = 12)),
  bottom = textGrob("True Effect Size (Cohen's d)", 
                    gp = gpar(fontsize = 12)),
  left = textGrob("Average Estimated Effect Size - True Effect Size", 
                  gp=gpar(fontsize = 12), 
                  rot = 90, 
                  just = "center"))

dev.off()

# Figure S3 ---------------------------------------------------------------
# Bias-Variance in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing
# Add Google Drive > Independent Segments Procedure > Simulations > Figure 7 (Bias-Variance in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing).  Code in 'Figures.Rmd'

proc.label = c("Fixed Sample Hypothesis Test", "Independent Segments Procedure")
names(proc.label) = c("Fixed", "ISP")
plt = function(f) {
  D = mean(f$d_forpower) #get separate plot for each d powered for
  t = paste0(sprintf("H\u00b2: \u03b4 = %0.1f, 1 - *\u03b2* = ", D), mean(f$power))
  
  ggplot(data = f, mapping = aes(x = d_actual, y = bias2)) +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 12), se = F, color = "black") +
    #geom_smooth(color = "black", se = F, size = .6, method = "gam") + #bias2
    geom_smooth(mapping = aes(y = mse), se = F, size = .6, method = "gam",
                linetype = "dashed", color = "red") + #MSE
    geom_smooth(mapping = aes(y = var), linetype = "dotted", se = F, 
                size = .6, method = "gam") + #variance
    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       labels = paste0("\u03B4 = ", seq(0, 1, 0.2))) +
    scale_y_continuous(breaks = c(0, max(f$mse)/2,  max(f$mse)),
                       labels = c(0, round(max(f$mse)/2, 2), round(max(f$mse), 2)),
                       limits = c(0, max(f$mse))) +
    theme_bw() +
    labs(title = t) +
    facet_wrap(vars(proc),labeller = labeller(proc = proc.label)) +
    theme(axis.title = element_blank(),
          strip.background =element_rect(fill="white"),
          plot.title = ggtext::element_markdown(size = 12, hjust = 0.5),
          axis.text.x = ggtext::element_markdown(size = 6, angle = 45, hjust = 1))
}

bvp = lapply(split(df_summary, list(df_summary$d_forpower, df_summary$power)), plt) #create 12 plots, one for each combination of d_forpower and power

facet1 = data.frame(proc = "Fixed", 
                    label1 = "Mean Squared Error",
                    label2 = "Variance",
                    label3 = "Squared Bias",
                    x1 = 0.1, y1 = .025, y2 = .020, y3 = .015, x0 = 0)

bvp[[1]] = bvp[[1]] +
  geom_text(data = facet1, aes(x = x1, y = y1, label = label1, hjust = 0), size = 3,
            color = "red") +
  geom_text(data = facet1, aes(x = x1, y = y2, label = label2, hjust = 0), size = 3,
            color = "blue") +
  geom_text(data = facet1, aes(x = x1, y = y3, label = label3, hjust = 0), size = 3) +
  geom_segment(data = facet1, aes(x = x0, xend = x1, y = y1, yend = y1),
               linetype = "dashed", color = "red") +
  geom_segment(data = facet1, aes(x = x0, xend = x1, y = y2, yend = y2),
               linetype = "dotted", color = "blue") +
  geom_segment(data = facet1, aes(x = x0, xend = x1, y = y3, yend = y3))

quartz(type = 'pdf', file = 'supplemental/figures/figureS3.pdf', width=20, height = 10.5)
grid.arrange(
  grobs = bvp, 
  ncol = 4, 
  top = textGrob(
    "Bias-Variance in Obtained Effect Sizes using Fixed versus Independent Segments Hypothesis Testing",
    gp = gpar(fontsize = 12)
    ), 
  bottom = textGrob(
    "True Effect Size (Cohen's \u03B4)", 
    gp = gpar(fontsize = 12)
    )
  )
dev.off()
# Figure S4 ---------------------------------------------------------------
# Exaggeration of Significant Results using Fixed versus Independent Segments Hypothesis Testing
# Add Google Drive > Independent Segments Procedure > Simulations > Figure 6 (Estimated Cohen's d using Fixed versus Independent Segments Hypothesis Testing).  Code in 'Figures Appendix.Rmd'


