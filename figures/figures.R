#.rs.restartR() #restart when corrupt
library(beepr)
library(ggplot2)
library(tidyverse)
library(ggtext)
library(grid)
library(papaja)
library(gridExtra)
library(metafor) #For an overview and introduction to the package please type: help(metafor).
library(rpact)
library(BayesFactor)
library(metafor)

#Insert code section with Ctrl+Shift+R

# TO DO -------------------------------------------------------------------
# Figure 6: Add legend for corrected ES estimate 
# Add line for correct ES estimate BFDA?
# Check why citations don't render correctly (perhaps add re-add those papers to Zotero?)
# Change figures after taking out MSPRT (height/width/labels with \n)
    # Set Fixed Sample Figure to left, and the Sequential Figures to right in a 2x2
    # OR: Add extra panel for Fixed Sample Test only?
# Change group-sequential designs figures (binding futility bounds + different rates for OBF changed the N). # Change this^ in manuscript too
# Change table 2a and 2b code in manuscript

# MAYBE -------------------------------------------------------------------
# Run meta-analysis and include meta-analytic effect size figure 
# Take figure 9 out of paper
# Rerun simulations on Sherlock?
# Run with 90% power (ideally, both?)
# Make density by segment plot
# Also vary number of looks 
# Add a table that show results for all effect sizes (maybe for supplemental?)
# Change subscript in the Bayes plots (fig 2a and 2b) (don't remember what I meant by this.. maybe +0?)

#=============================================================================#
########################## FIGURE 1A: POCOCK EXAMPLE ########################## 
#=============================================================================#

alpha_gs = function(proc, return = "n") {
  bs = ifelse(proc == "asOF", "bsOF", "bsP")
  
  # Specify the design
  design = getDesignGroupSequential(
    sided = 1,
    alpha = 0.05, 
    beta = 0.2,
    kMax = 3, 
    typeOfDesign = proc,
    typeBetaSpending= bs) #also specify information rates
  
  # Get parameters
  parameters = getSampleSizeMeans(design = design, groups = 2, alternative = 0.5)
  n_gs = ceiling(parameters$maxNumberOfSubjects/3/2) # n per look per group
  alpha = parameters$criticalValuesPValueScale
  futility = parameters$futilityBoundsPValueScale
  
  if(return == "alpha") {return(alpha)}
  if(return == "n") {return(n_gs)}
  if(return == "futility") {return(futility)}
}

n_p = alpha_gs(proc = "asP")
alpha_p = alpha_gs(proc = "asP", return = "alpha")
futility_p = alpha_gs(proc = "asP", return = "futility") #futility on critical p-value scale
futility_pd = round(
  c(qt(1 - futility_p[1], df = (n_p*2)-2)/sqrt(n_p/2), 
    qt(1 - futility_p[2], df = 4*n_p-2)/sqrt(2*n_p/2)), 
  2) #futility for pococock on critical d-value

n_of = alpha_gs(proc = "asOF")
alpha_of = alpha_gs(proc = "asOF", return = "alpha")
futility_of = alpha_gs(proc = "asOF", return = "futility")  
futility_ofd = round(
  c(qt(1 - futility_of[1], df = (n_of*2)-2)/sqrt(n_of/2), 
    qt(1 - futility_of[2], df = 4*n_of-2)/sqrt(2*n_of/2)), 
  2) #futility for OBF on critical d-value scale

d.plot = function(n, alpha, t) {
  df = 2*n-2
  critical_d = round(qt(1-alpha, df = 2*n-2) * sqrt(2*(1/n)),2)
  p <- ggplot(data.frame(x = c(-4, 4)), aes(x)) +
    stat_function(fun = dt, args = list(df = df)) +
    stat_function(fun = dt, args = list(df = df), 
                  xlim = c(qt(p = 1 - alpha, df = df), 4),
                  geom = "area", fill = "red", alpha = 0.6) +
    labs(title = t, y = "Density", subtitle = "*p*-value") +
    theme_bw() +
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
          plot.subtitle = ggtext::element_markdown(size = 13, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 11)) + 
    geom_vline(xintercept = qt(1-alpha, df = df), linetype = 2, color = "grey") 
  return(p)
}


p1 = d.plot(n = n_p, alpha = alpha_p[1], t = "Pocock (First Look)") +
  scale_x_continuous(
    limits = c(0, 1*sqrt(n_p/2)), 
    breaks = c(0, .2*sqrt(n_p/2), .4*sqrt(n_p/2), .6*sqrt(n_p/2), .8*sqrt(n_p/2)),
    labels = c(0, .2, .4, .6, .8),
    sec.axis = dup_axis( #create secondary p-value axis
      labels = c(
        signif(
          sapply(
            c(0, .2*sqrt(n_p/2), 
              .4*sqrt(n_p/2), 
              .6*sqrt(n_p/2), 
              .8*sqrt(n_p/2)), 
            function(x) pt(-abs(x), df = 2*n_p-2)), 1)))) +
  annotate(geom = "text", x = 1.525, y = 0.29, fontface = "bold", size = 2.3,
           label = paste0(
             "CONTINUE \n", 
             futility_pd[1],
             " \u2264 d < ",
             round(qt(1-alpha_p[1], df = 2*n_p-2) * sqrt(2*(1/n_p)),2),
             "\n",
             signif(alpha_p[1], 2),
             " < p \u2264 ",
             signif(futility_p[1], 2))) +
  annotate(geom = "text", x = 3, y = 0.29,   
           label = paste0(
             "REJECT \n d \u2265 ", 
             round(qt(1-alpha_p[1], df = 2*n_p-2) * sqrt(2*(1/n_p)),2),
             "\n p \u2264 ", 
             gsub("0.", ".", as.character(signif(alpha_p[1], 2)), fixed = T)), 
           fontface = "bold", size = 2.3) +
  geom_vline(xintercept = futility_pd[1] *sqrt(n_p/2), linetype = "dashed", color = "grey") +
  stat_function(
    fun = dt, args = list(df = 2*n_p -2), 
    xlim = c(0, qt(1 - futility_p[1], df = 2*n_p -2)),
    geom = "area", fill = "grey", alpha = 0.3) +
  annotate(geom = "text", x = 0.1, y = 0.29, label = "FAIL \n TO \n REJECT", fontface = "bold", size = 2.3) +
  annotate(geom = "text", x = 3.2, y = 0.07,   label = paste0("N = ", n_p), fontface = "italic", size = 2.3)

p2 = d.plot(n = 2*n_p, alpha = alpha_p[2], t = "Pocock (Second Look)") +
  scale_x_continuous(
    limits = c(0.05*sqrt(2*n_p/2), 0.6*sqrt(2*n_p/2)),
    breaks = c(.2*sqrt(2*n_p/2), .4*sqrt(2*n_p/2)),
    labels = c(.2, .4),
    sec.axis = dup_axis(
      labels = c(signif(
        sapply(
          c(.2*sqrt(2*n_p/2), 
            .4*sqrt(2*n_p/2)), 
          function(x) pt(-abs(x), df = 4*n_p-2)), 1)))) +
  geom_vline(xintercept = futility_pd[2] *sqrt(2*n_p/2), linetype = "dashed", color = "grey") +
  stat_function(
    fun = dt, args = list(df = 4*n_p -2), 
    xlim = c(0, qt(1 - futility_p[2], df = 4*n_p -2)),
    geom = "area", fill = "grey", alpha = 0.3) +
  annotate(geom = "text", x = 1.625, y = 0.29, fontface = "bold", size = 2.3,
           label = paste0(
             "CONTINUE \n", 
             futility_pd[2],
             " \u2264 d < ",
             round(qt(1-alpha_p[2], df = 4*n_p-2) / sqrt(2*n_p/2),2),
             "\n",
             signif(alpha_p[2], 2),
             " < p \u2264 ",
             signif(futility_p[2], 2))) +
  annotate(geom = "text", x = 2.6, y = 0.29,   
           label = paste0(
             "REJECT \n d \u2265 ", 
             round(qt(1-alpha_p[2], df = 4*n_p-2) / sqrt(2*n_p/2),2),
             "\n p \u2264 ", 
             gsub("0.", ".", as.character(signif(alpha_p[2], 2)), fixed = T)), 
           fontface = "bold", size = 2.3) +
  annotate(geom = "text", x = 0.43, 
           y = 0.29,   
           label = "FAIL \n TO \n REJECT", 
           fontface = "bold", 
           size = 2.3) +
  annotate(geom = "text", 
           x = 2.75, 
           y = 0.07,
           label = paste0("N = ", 2*n_p), #check if this is still correct
           fontface = "italic", 
           size = 2.3)

p3 = d.plot(n = 3*n_p, alpha = alpha_p[3], t = "Pocock (Third Look)") +
  scale_x_continuous(
    limits = c(0, 0.6*sqrt(3*n_p/2)),
    breaks = c(0, .2*sqrt(3*n_p/2), .4*sqrt(3*n_p/2), .6*sqrt(3*n_p/2)),
    labels = c(0, .2, .4, .6),
    sec.axis = dup_axis(
      labels = c(signif(
        sapply(c(
          0, .2*sqrt(3*n_p/2), 
          .4*sqrt(3*n_p/2), 
          .6*sqrt(3*n_p/2)), 
          function(x) pt(-abs(x), df = 6*n_p-2)), 1)))) +
  stat_function(
    fun = dt, 
    args = list(df = 6*n_p -2), 
    xlim = c(0, qt(1 - alpha_p[3], df = 6*n_p-2)),
    geom = "area", 
    fill = "grey", 
    alpha = 0.3) +
  annotate(geom = "text", 
           x = 0.23, 
           y = 0.29,   
           label = "FAIL \n TO \n REJECT", 
           fontface = "bold", 
           size = 2.3) +
  annotate(geom = "text", 
           x = 3, 
           y = 0.29,   
           label = paste0(
             "REJECT \n d \u2265 ", 
             round(qt(1-alpha_p[3], df = 6*n_p-2) / sqrt(3*n_p/2),2),
             "\n p \u2264 ", 
             gsub("0.", ".", as.character(signif(alpha_p[3], 2)), fixed = T)), 
           fontface = "bold", size = 2.3)  +
  annotate(geom = "text", 
           x = 3.3, 
           y = 0.07,   
           label = paste0("N = ", 3*n_p), 
           fontface = "italic", 
           size = 2.3)

tiff(file="figures/figure1a.tiff",width=2500,height=800, units = "px", res = 300)
grid.arrange(p1, p2, p3, nrow = 1, bottom = "Cohen's d", 
             left = textGrob("Density under H0: \u03b4 = 0", rot = 90, hjust = 0.57,
                             gp = gpar(fontsize = 10)))
dev.off()


#======================================================================================#
########################## FIGURE 1B: O'BRIEN-FLEMING EXAMPLE ########################## 
#======================================================================================#

of1 = d.plot(n = n_of, alpha = alpha_of[1], t = "O'Brien-Fleming (First Look)") +
  scale_x_continuous(
    limits = c(0, 1.5*sqrt(n_of/2)), 
    breaks = c(0, .4*sqrt(n_of/2), .8*sqrt(n_of/2), 1.2*sqrt(n_of/2)),
    labels = c(0, .4, .8, 1.2),
    sec.axis = dup_axis(
      labels = c(
        signif(
          sapply(
            c(0, .4*sqrt(n_of/2), 
              .8*sqrt(n_of/2), 
              1.2*sqrt(n_of/2)), 
            function(x) pt(-abs(x), df = 2*n_of-2)), 1)))) +
  annotate(geom = "text", x = 2, y = 0.27,   fontface = "bold", size = 2.3,
           label = paste0(
             "CONTINUE \n", 
             futility_ofd[1],
             " \u2264 d < ",
             round(qt(1-alpha_of[1], df = 2*n_of-2) * sqrt(2*(1/n_of)),2),
             "\n",
             signif(alpha_of[1], 2),
             " < p \u2264 ",
             signif(futility_of[1], 2))) +
  annotate(geom = "text", x = 4.1, y = 0.27,   
           label = paste0(
             "REJECT \n d \u2265 ", 
             round(qt(1-alpha_of[1], df = 2*n_of-2) * sqrt(2*(1/n_of)),2),
             "\n p \u2264 ", 
             gsub("0.", ".", as.character(signif(alpha_of[1], 2)), fixed = T)), 
           fontface = "bold", size = 2.3) +
  annotate(geom = "text", 
           x = 4.3, 
           y = 0.07,   
           label = paste0("N = ", n_of), 
           fontface = "italic", 
           size = 2.3)

of2 = d.plot(n = 2*n_of, alpha = alpha_of[2], t = "O'Brien-Fleming (Second Look)") +
  scale_x_continuous(
    limits = c(0, 0.8*sqrt(2*n_of/2)),
    breaks = c(0, .2*sqrt(2*n_of/2), .4*sqrt(2*n_of/2), .6*sqrt(2*n_of/2)),
    labels = c(0, .2, .4, .6),
    sec.axis = dup_axis(labels = c(signif(
      sapply(c(0, .2*sqrt(2*n_of/2), 
               .4*sqrt(2*n_of/2), 
               .6*sqrt(2*n_of/2)), 
             function(x) pt(-abs(x), df = 4*n_of-2)), 1)))) +
  geom_vline(xintercept = futility_ofd[2] *sqrt(2*n_of/2), linetype = "dashed", color = "grey") +
  stat_function(
    fun = dt, 
    args = list(df = 4*n_of -2), 
    xlim = c(0, qt(1 - futility_of[2], df = 4*n_of -2)),
    geom = "area", 
    fill = "grey", 
    alpha = 0.3) +
  annotate(geom = "text", x = 1.65, y = 0.27,   fontface = "bold", size = 2.3,
           label = paste0(
             "CONTINUE \n", 
             futility_ofd[2],
             " \u2264 d < ",
             round(qt(1-alpha_of[2], df = 4*n_of-2) / sqrt(2*n_of/2),2),
             "\n",
             signif(alpha_of[2], 2),
             " < p \u2264 ",
             signif(futility_of[2], 2))) +
  annotate(geom = "text", x = 2.8, y = 0.27,   
           label = paste0(
             "REJECT \n d \u2265 ", 
             round(qt(1-alpha_of[2], df = 4*n_of-2) / sqrt(2*n_of/2),2),
             "\n p \u2264 ", 
             gsub("0.", ".", as.character(signif(alpha_of[2], 2)), fixed = T)), 
           fontface = "bold", size = 2.3) +
  annotate(geom = "text", 
           x = 0.25, 
           y = 0.27,   
           label = "FAIL \n TO \n REJECT", 
           fontface = "bold", 
           size = 2.3) +
  annotate(geom = "text", 
           x = 3.2, 
           y = 0.07,   
           label = paste0("N = ", 2*n_of), 
           fontface = "italic", 
           size = 2.3)

of3 = d.plot(n = 3*n_of, alpha = alpha_of[3], t = "O'Brien-Fleming (Third Look)") +
  scale_x_continuous(
    limits = c(0, 0.6*sqrt(3*n_of/2)),
    breaks = c(0, .2*sqrt(3*n_of/2), .4*sqrt(3*n_of/2), .6*sqrt(3*n_of/2)),
    labels = c(0, .2, .4, .6),
    sec.axis = dup_axis(labels = c(signif(
      sapply(c(0, 
               .2*sqrt(3*n_of/2), 
               .4*sqrt(3*n_of/2), 
               .6*sqrt(3*n_of/2)), 
             function(x) pt(-abs(x), df = 6*n_of-2)), 1)))) +
  stat_function(
    fun = dt, 
    args = list(df = 6*n_of -2), 
    xlim = c(0, qt(1 - alpha_of[3], df = 6*n_of-2)),
    geom = "area", 
    fill = "grey", 
    alpha = 0.3) +
  annotate(geom = "text", 
           x = 0.25, 
           y = 0.27,   
           label = "FAIL \n TO \n REJECT", 
           fontface = "bold", 
           size = 2.3) +
  annotate(geom = "text", 
           x = 2.8, 
           y = 0.27,   
           label = paste0(
             "REJECT \n d \u2265 ", 
             round(qt(1-alpha_of[3], df = 6*n_of-2) / sqrt(3*n_of/2),2),
             "\n p \u2264 ", 
             gsub("0.", ".", as.character(signif(alpha_of[3], 2)), fixed = T)), 
           fontface = "bold", 
           size = 2.3) +
  annotate(geom = "text", 
           x = 3, 
           y = 0.07,   
           label = paste0("N = ", 3*n_of), 
           fontface = "italic", 
           size = 2.3)

tiff(file="figures/figure1b.tiff",width=2500,height=800, units = "px", res = 300)
grid.arrange(of1, of2, of3, nrow = 1, bottom = "Cohen's d", 
             left = textGrob("Density under H0: \u03b4 = 0", rot = 90, hjust = 0.57, 
                             gp = gpar(fontsize = 10)))
dev.off()

#============================================================================#
########################## FIGURE 2AB: BAYES EXAMPLE ########################## 
#============================================================================#
# Adapted from https://www.shinyapps.org/apps/RGraphCompendium/index.php#prior-and-posterior
# Big thank you to Eric-Jan Wagenmakers & Quentin Gronau for maintaining this site!
# I took the code from the side above and recoded it from base R to ggplot

.likelihoodShiftedT <- function(par, data) {
  -sum(log(dt((data - par[1])/par[2], par[3])/par[2]))
}

.dposteriorShiftedT <- function(x, parameters) {
  # Function that returns the posterior density
  ifelse(x >= 0, (dt((x - parameters[1])/parameters[2], parameters[3])/parameters[2])/
           pt((0 - parameters[1])/parameters[2], parameters[3], lower.tail = FALSE), 
         0)
}

.dprior <- function(x, r) {
  # Function that returns the prior density
  y <- ifelse(x < 0, 0, 2/(pi * r * (1 + (x/r)^2)))
  return(y)
}

.plotPosterior.ttest.ggplot <- function(x = NULL, 
                                        BF, 
                                        BFH1H0 = TRUE, 
                                        iterations = 10000, 
                                        rscale = "medium") {
  
  if (rscale == "medium") {r <- sqrt(2)/2}
  if (rscale == "wide") {r <- 1}
  if (rscale == "ultrawide") {r <- sqrt(2)}
  if (mode(rscale) == "numeric") {r <- rscale}
  
  nullInterval <- c(0, Inf)
  
  # Sample from delta posterior
  samples <- BayesFactor::ttestBF(x = x, 
                                  posterior = TRUE, 
                                  iterations = iterations, 
                                  rscale = r)
  
  delta <- samples[, "delta"]
  
  # Fit shifted t distribution
  deltaHat <- mean(x)/sd(x)
  N <- length(x)
  df <- N - 1
  sigmaStart <- 1/N
  
  if (sigmaStart < 0.01) 
    sigmaStart <- 0.01
  
  parameters <- try(
    silent = TRUE, 
    expr = optim(par = c(deltaHat, sigmaStart, df), 
                 fn = .likelihoodShiftedT, 
                 data = delta, 
                 method = "BFGS")$par)
  
  if (class(parameters) == "try-error") {
    parameters <- try(
      silent = TRUE, 
      expr = optim(par = c(deltaHat, sigmaStart, df), 
                   fn = .likelihoodShiftedT, 
                   data = delta, 
                   method = "Nelder-Mead")$par)
  }
  
  if (BFH1H0) {
    BF10 <- BF
    BF01 <- 1/BF10
  } else {
    BF01 <- BF
    BF10 <- 1/BF01
  }
  
  # Set limits plot
  xlim <- vector("numeric", 2)
  xlim[1] <- min(-1, quantile(delta[delta >= 0], probs = 0.01)[[1]])
  xlim[2] <- max(2, quantile(delta[delta >= 0], probs = 0.99)[[1]])
  
  stretch <- 1.32
  xticks <- pretty(xlim)
  ylim <- vector("numeric", 2)
  
  ylim[1] <- 0
  dmax <- optimize(function(x) 
    .dposteriorShiftedT(x, 
                        parameters = parameters), 
    interval = range(xticks), 
    maximum = TRUE)$objective
  
  ylim[2] <- max(stretch * .dprior(0, r), stretch * dmax)  # get maximum density
  
  # Calculate position of 'nice' tick marks and create labels
  yticks <- pretty(ylim)
  xlabels <- formatC(xticks, 1, format = "f")
  ylabels <- formatC(yticks, 1, format = "f")
  
  # compute 95% credible interval & median:
  CIlow <- quantile(delta[delta >= 0], probs = 0.025)[[1]]
  CIhigh <- quantile(delta[delta >= 0], probs = 0.975)[[1]]
  medianPosterior <- median(delta[delta >= 0])
  
  posteriorLine <- .dposteriorShiftedT(
    x = seq(min(xticks), max(xticks), length.out = 1000), 
    parameters = parameters)
  
  priorLine <- .dprior(seq(min(xticks), max(xticks), length.out = 1000), r = r)
  
  xlim <- c(min(CIlow, range(xticks)[1]), max(range(xticks)[2], CIhigh))
  
  ##ggplot 
  ggplot(mapping = aes(x = seq(xlim[1], xlim[2], length.out = 1000), 
                       y = posteriorLine)) +
    geom_line(size = 1) +
    geom_line(mapping = aes(y = priorLine), linetype = "dashed", size = 1) + 
    geom_point(mapping = aes(x = 0, y = .dprior(0, r)), 
               size = 5, fill = "grey", shape = 21) + 
    geom_point(mapping = aes(x = 0, y = posteriorLine[posteriorLine >0][1]),
               size = 5, fill = "grey", shape = 21) +
    theme_classic() +
    scale_x_continuous(breaks = xticks,
                       labels = xlabels) +
    scale_y_continuous(breaks = yticks,
                       labels = ylabels,
                       limits = c(min(yticks), max(yticks))) +
    labs(x = expression(paste("Effect size ", delta)),
         y = "Density") +
    theme(text = element_text(size = 25),
          axis.ticks.length = unit(.45, "cm")) +
    annotate("text", x = xticks[2], y = ylim[2] - 0.3, size = 6,
             label = paste0("BF0+ = ", format(round(BF01, 2), nsmall = 2))) +
    annotate("text", x = xticks[2], y = ylim[2], size = 6,
             label = paste0("BF+0 = ", format(round(BF10, 2), nsmall = 2))) +
    geom_line(mapping = aes(x = c(xticks[length(xticks) - 2],
                                  xticks[length(xticks) - 2] + 0.4),
                            y = c(ylim[2] - 0.3, ylim[2] - 0.3)),
              linetype = "dashed", size = 1) +
    geom_line(mapping = aes(x = c(xticks[length(xticks) - 2],
                                  xticks[length(xticks) - 2] + 0.4),
                            y = c(ylim[2], ylim[2])),
              size = 1) + 
    annotate("text", x = xticks[length(xticks) -1], y = ylim[2] - 0.3,
             label = "Prior", size = 6, hjust = 0) + 
    annotate("text", x = xticks[length(xticks) - 1], y = ylim[2],
             label = "Posterior", size = 6, hjust = 0)
  
}

### SUPPORT FOR THE ALTERNATIVE H+: DELTA > 0

# Generate data
set.seed(2)
x <- rnorm(22, 0.30355)

# Calculate Bayes factor
BF = extractBF(ttestBF(x, rscale = "medium", nullInterval = c(0, Inf)))$bf[1]

# Plot
BF10_plot = .plotPosterior.ttest.ggplot(x = x, BF = BF)
ggsave(filename = "figures/figure2a.tiff", device = "tiff", width = 6, height = 6, dpi = 300)

### SUPPORT FOR THE NULL 

# Generate data
set.seed(2)
x <- rnorm(22, -0.1031)

# Calculate Bayes factor
BF = extractBF(ttestBF(x, rscale = "medium", nullInterval = c(0, Inf)))$bf[1]

# Plot
BF01_plot = .plotPosterior.ttest.ggplot(x = x, BF = BF)
ggsave(filename = "figures/figure2b.tiff", device = "tiff", width = 6, height = 6, dpi = 300)

#===================================================================#
########################## FIGURE 3A: INTRO ########################## 
#===================================================================#
n_fixed = 51
n_s = 25
alpha_total = 0.05
alpha_strong = 0.025
alpha_weak = 0.281784510170711 #See Simulations.R for function alpha_weak
max_n_segments = 3

critical.d = function(proc, segment) {
  alpha = ifelse(proc == "Fixed", alpha_total,
                 ifelse(segment == max_n_segments, alpha_weak, alpha_strong))
  n = ifelse(proc == "Fixed", n_fixed, n_s)
  t = ifelse(proc == "Fixed", "Fixed Sample Hypothesis Test",
             paste0("ISP ", ifelse(segment == 3, "(Last Segment)", "(First Two Segments)")))
  df = 2*n-2
  
  p <- ggplot(data.frame(x = c(-4, 4)), aes(x)) +
    stat_function(fun = dt, args = list(df = df)) +
    stat_function(fun = dt, args = list(df = df), 
                  xlim = c(qt(p = 1 - alpha, df = df), 4),
                  geom = "area", fill = "red", alpha = 0.6) +
    labs(title = t, y = "Density") +
    theme_classic() +
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
          plot.subtitle = ggtext::element_markdown(size = 13, hjust = 0.5),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_x_continuous(limits = c(0, .8*sqrt(n/2))) +
    geom_vline(xintercept = qt(1-alpha, df = df), linetype = 2, color = "grey") +
    stat_function(fun = dt, args = list(df = df), 
                  xlim = c(-1, qt(1 - ifelse(proc == "Fixed", alpha, alpha_weak), df = df)),
                  geom = "area", fill = "grey", alpha = 0.3) +
    annotate(geom = "text", x = 0.26, y = 0.29,   
             label = "FAIL \n TO \n REJECT", fontface = "bold", size = 2.3)
  return(p)
}

# Create-plots
p1 = critical.d(proc = "Fixed") + 
  annotate(geom = "text", x = 3, y = 0.29, size = 2.3, 
           label = paste0("REJECT \n p \u2264 \u03b1"), fontface = "bold") +
  scale_x_continuous(limits = c(0, .8*sqrt(n_fixed/2)),
                     breaks = qt(1-alpha_total, df = 2*n_fixed-2),
                     labels = "\u03b1") +
  theme(axis.text.x = element_text(size = 12, face = "bold")) 

p2 = critical.d(proc = "ISP", segment = 1) +
  annotate(geom = "text", x = 2.5, y = 0.29, size = 2.3, 
           label = paste0("REJECT \n p \u2264 \u03b1 strong"), fontface = "bold") +
  geom_vline(xintercept = qt(1-alpha_weak, df = 2*n_s-2), linetype = 2, color = "grey") +
  annotate(geom = "text", x = 1.45, y = 0.29,   size = 2.3,
           label = paste0("CONTINUE \n \u03b1 strong < p \u2264 \u03b1 weak"),
           fontface = "bold") +
  scale_x_continuous(limits = c(0, .8*sqrt(n_s/2)),
                     breaks = c(qt(1 - alpha_weak, df = 2*n_s-2), 
                                qt(1 - alpha_strong, df = 2*n_s-2)),
                     labels = c("\u03b1 weak", "\u03b1 strong" )) +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

p3 = critical.d(proc = "ISP", segment = 3) +
  annotate(geom = "text", x = 2.5,y = 0.29, size = 2.3, 
           label = paste0("REJECT \n p \u2264 \u03b1 weak"), fontface = "bold") +
  scale_x_continuous(limits = c(0, .8*sqrt(n_s/2)),
                     breaks = qt(1 - alpha_weak, df = 2*n_s-2),
                     labels = "\u03b1 weak") +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

tiff(file="figures/figure3a.tiff",width=2800,height=700, units = "px", res = 300)
grid.arrange(p1, p2, p3, nrow = 1, bottom = "Cohen's d", 
             left = textGrob("Density under H0: \u03b4 = 0", rot = 90, hjust = 0.4,
                             gp = gpar(fontsize = 10)))
dev.off()

#=================================================================================#
########################## FIGURE 3B: FIXED + ISP EXAMPLE ########################## 
#=================================================================================#

critical.d = function(proc, segment) {
  alpha = ifelse(proc == "Fixed", alpha_total,
                 ifelse(segment == max_n_segments, alpha_weak, alpha_strong))
  n = ifelse(proc == "Fixed", n_fixed, n_s)
  t = ifelse(proc == "Fixed", "Fixed Sample Hypothesis Test",
             paste0("ISP ", ifelse(segment == 3, "(Last Segment)", "(First Two Segments)")))
  df = 2*n-2
  critical_d = round(qt(1-alpha, df = 2*n-2) * sqrt(2*(1/n)),2)
  
  p <- ggplot(data.frame(x = c(-4, 4)), aes(x)) +
    stat_function(fun = dt, args = list(df = df)) +
    stat_function(fun = dt, args = list(df = df), 
                  xlim = c(qt(p = 1 - alpha, df = df), 4),
                  geom = "area", fill = "red", alpha = 0.6) +
    labs(title = t, y = "Density", subtitle = "*p*-value") +
    theme_bw() +
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
          plot.subtitle = ggtext::element_markdown(size = 13, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 11)) + 
    scale_x_continuous(
      limits = c(0, .8*sqrt(n/2)), 
      breaks = c(0, .2*sqrt(n/2), .4*sqrt(n/2), .6*sqrt(n/2)),
      labels = c(0, .20, .40, .60),
      sec.axis = dup_axis(labels = c(signif(
        sapply(c(0, 
                 .2*sqrt(n/2), 
                 .4*sqrt(n/2), 
                 .6*sqrt(n/2)), 
               function(x) pt(-abs(x), df = df)), c(2, 2, 1, 1))))) +
    geom_vline(xintercept = qt(1-alpha, df = df), linetype = 2, color = "grey") +
    stat_function(fun = dt, args = list(df = df), 
                  xlim = c(0, qt(1 - ifelse(proc == "Fixed", alpha, alpha_weak), df = df)),
                  geom = "area", fill = "grey", alpha = 0.3) +
    annotate(geom = "text", 
             x = 0.25, 
             y = 0.27,   
             label = "FAIL \n TO \n REJECT", 
             fontface = "bold", 
             size = 2.3) +
    annotate(geom = "text", 
             x = 2.5, 
             y = 0.27,   
             label = paste0(
               "REJECT \n d \u2265 ", 
               critical_d, 
               "\n p \u2264 ", 
               gsub("0.", ".", as.character(signif(alpha, 2)), fixed = T)), 
             fontface = "bold", size = 2.3)
  return(p)
}

p1 = critical.d(proc = "Fixed") +
  annotate(geom = "text", 
           x = 3.8, 
           y = 0.07,   
           label = paste0("N = ", n_fixed), 
           fontface = "italic", 
           size = 2.3)

p2 = critical.d(proc = "ISP", segment = 1) +
  geom_vline(xintercept = qt(1-alpha_weak, df = 2*n_s-2), linetype = 2, color = "grey") +
  annotate(geom = "text", 
           x = 1.5, 
           y = 0.27,   
           size = 2.3,
           label = paste0(
             "CONTINUE \n ", 
             round(qt(1-alpha_weak, df = 2*n_s-2) / sqrt(n_s/2),2),
             " \u2264 d < ", 
             round(qt(1 - alpha_strong, df = 2*n_s-2) / sqrt(n_s/2), 2), "\n",
             alpha_strong, 
             " < p \u2264 ", gsub("0.", ".", as.character(round(alpha_weak, 2)), fixed = T)),
           fontface = "bold") +
  annotate(geom = "text", 
           x = 2.65, 
           y = 0.07,   
           label = paste0("N = ", n_s), 
           fontface = "italic", 
           size = 2.3)

p3 = critical.d(proc = "ISP", segment = 3) +
  annotate(geom = "text", 
           x = 2.65, 
           y = 0.07,   
           label = paste0("N = ", n_s), 
           fontface = "italic", 
           size = 2.3)

tiff(file="figures/figure3b.tiff",width=2500,height=800, units = "px", res = 300)
grid.arrange(p1, p2, p3, nrow = 1, bottom = "Cohen's d", 
             left = textGrob("Density under H0: \u03b4 = 0", rot = 90, hjust = 0.57,
                             gp = gpar(fontsize = 10)))
dev.off()
#=========================================================================#
########################## FIGURE 4: ERROR RATES ########################## 
#=========================================================================#
# WITHOUT THE MSPRT
df = read.csv("simulations/data.csv")

df = df %>% 
  mutate(
    proc = fct_relevel(proc, "Fixed", after = 0),
    facet = case_when(
      proc == "Fixed" ~ 1,
      proc == "Bayes" ~ 2,
      proc == "ISP" ~ 3,
      proc == "asP" ~ 4,
      proc == "asOF"~ 5)
    )

facet.label = c("Fixed Hypothesis and Equivalence Test", 
                "Sequential Bayes Factor", 
                "Independent Segments Procedure", 
                "Pocock-like GS design",
                "O'Brien-Fleming-like GS design")
names(facet.label) = c(1, 2, 3, 4, 5)

ER = df %>% 
  filter(d_actual != -0.2) %>% 
  mutate(
    correct_inference = ifelse(
      d_actual > 0 & decision == "reject.null", 1,
      ifelse(d_actual == 0 & decision == "reject.alt", 1, 
             0)),
    incorrect_inference = ifelse(
      d_actual == 0 & decision == "reject.null", 1,
      ifelse(d_actual > 0 & decision == "reject.alt", 1,
             0)),
    inconclusive = ifelse(decision == "inconclusive", 1, 0))

ER_summary = ER %>% 
  group_by(proc, facet, d_actual) %>% 
  summarize(
    correct_rate = mean(correct_inference),
    incorrect_rate = mean(incorrect_inference),
    inconclusive_rate = mean(inconclusive))

tiff(file="figures/figure4.tiff",width=2350,height=1500, units = "px", res = 300)
ggplot(data = ER_summary %>% filter(d_actual > 0), mapping = aes(x = d_actual, y = correct_rate)) +
  facet_wrap(vars(facet), labeller = labeller(facet = facet.label)) +
  geom_line(color = "green", alpha = 0.7) +
  geom_point(shape = "v", size = 3, color = "green") +
  geom_line(mapping = aes(y = incorrect_rate), color = "red", alpha = 0.7) +
  geom_line(mapping = aes(y = inconclusive_rate), color = "grey", linetype = "dashed", alpha = 0.7) +
  geom_point(data = ER_summary %>% filter(d_actual < 0.8),
             mapping = aes(y = inconclusive_rate), 
             shape = "?", 
             size = 3, 
             color = "grey") +
  geom_point(mapping = aes(y = incorrect_rate), 
             shape = "x", 
             size = 3, 
             color = "red") +
  theme_bw() +
  scale_x_continuous(limits = c(0.2, 1)) +
  scale_y_continuous(labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme(strip.background =element_rect(fill="white")) +
  labs(x = "True Effect Size (Cohen's d)", 
       y = "True Positives, False Negatives, Inconclusive Evidence")
dev.off()

# Error Rates Table 2a [see code in manuscript.Rmd] 
t1 = df %>% filter(d_forpower == d_actual) %>%
  group_by(proc) %>% 
  summarize(true_positive = mean(decision == "reject.null"),
            inconclusive = mean(decision == "inconclusive"),
            false_negative = mean(decision == "reject.alt"))

# Error Rates Table 2b [see code in manuscript.Rmd]
t2 = df %>% filter(d_actual == 0) %>% 
  group_by(proc) %>% 
  summarize(true_negative = mean(decision == "reject.alt"),
            inconclusive = mean(decision == "inconclusive"),
            false_positive = mean(decision == "reject.null"))


#=========================================================================#
########################## FIGURE 4: ERROR RATES ########################## 
#=========================================================================#
# WITH MSPRT
# 
# df = df %>% 
#   mutate(
#     proc = fct_relevel(proc, "Fixed", after = 0),
#     facet = case_when(proc == "Fixed" ~ 1,
#                       proc == "Bayes" ~ 2,
#                       proc == "ISP" ~ 3,
#                       proc == "asP" ~ 4,
#                       proc == "asOF"~ 5,
#                       proc == "SPRT" ~ 6))
# 
# facet.label = c("Fixed Hypothesis and Equivalence Test", "Sequential Bayes Factor", 
#                 "Independent Segments Procedure", "Pocock-like GS design",
#                 "O'Brien-Fleming-like GS design", "Modified Sequential Probability Ratio Test")
# names(facet.label) = c(1, 2, 3, 4, 5, 6)
# 
# ER = df %>% 
#   filter(d_actual != -0.2) %>% 
#   mutate(correct_inference = ifelse(d_actual > 0 & decision == "reject.null", 1,
#                                     ifelse(d_actual == 0 & decision == "reject.alt", 1, 
#                                            0)),
#          incorrect_inference = ifelse(d_actual == 0 & decision == "reject.null", 1,
#                                       ifelse(d_actual > 0 & decision == "reject.alt", 1,
#                                              0)),
#          inconclusive = ifelse(decision == "inconclusive", 1, 0))
# 
# ER_summary = ER %>% 
#   group_by(proc, facet, d_actual) %>% 
#   summarize(correct_rate = mean(correct_inference),
#             incorrect_rate = mean(incorrect_inference),
#             inconclusive_rate = mean(inconclusive))
# 
# tiff(file="figures/figure4.tiff",width=2350,height=1500, units = "px", res = 300)
# ggplot(data = ER_summary %>% filter(d_actual > 0), mapping = aes(x = d_actual, y = correct_rate)) +
#   facet_wrap(vars(facet), labeller = labeller(facet = facet.label)) +
#   geom_line(color = "green", alpha = 0.7) +
#   geom_point(shape = "v", size = 3, color = "green") +
#   geom_line(mapping = aes(y = incorrect_rate), color = "red", alpha = 0.7) +
#   geom_line(mapping = aes(y = inconclusive_rate), color = "grey", linetype = "dashed", alpha = 0.7) +
#   geom_point(data = ER_summary %>% filter(d_actual < 0.8),
#              mapping = aes(y = inconclusive_rate), shape = "?", size = 3, color = "grey") +
#   geom_point(mapping = aes(y = incorrect_rate), shape = "x", size = 3, color = "red") +
#   theme_bw() +
#   scale_x_continuous(limits = c(0.2, 1)) +
#   scale_y_continuous(labels = c("0%", "25%", "50%", "75%", "100%")) +
#   theme(strip.background =element_rect(fill="white")) +
#   labs(x = "True Effect Size (Cohen's d)", y = "True Positives, False Negatives, Inconclusive Evidence")
# dev.off()
# 
# # Error Rates Table 2a [see code in manuscript.Rmd] #change code in manuscript
# t1 = df %>% filter(d_forpower == d_actual) %>%
#   group_by(proc) %>% 
#   summarize(true_positive = mean(decision == "reject.null"),
#             inconclusive = mean(decision == "inconclusive"),
#             false_negative = mean(decision == "reject.alt"))
# 
# # Error Rates Table 2b [see code in manuscript.Rmd]
# t2 = df %>% filter(d_actual == 0) %>% 
#   group_by(proc) %>% 
#   summarize(true_negative = mean(decision == "reject.alt"),
#             inconclusive = mean(decision == "inconclusive"),
#             false_positive = mean(decision == "reject.null"))

#========================================================================#
########################## FIGURE 5: Efficiency ########################## 
#========================================================================#
# WITHOUT MSPRT
#facet.label = c("Fixed Sample Hypothesis Test", "Sequential Bayes Factor", 
#                "Independent Segments Procedure", "Group-Sequential Procedure")
#names(facet.label) = c(1, 2, 3, 4)

E_n = df %>% group_by(proc, d_actual) %>% 
  #mutate(proc = fct_relevel(proc, "Fixed", after = 0)) %>% 
  summarize(E_n = mean(n_cumulative)) # %>% 
#mutate(proc = case_when(proc == "Fixed" ~ "Fixed",
#                        proc == "asP" ~ "Pocock",
#                        proc == "asOF" ~ "O'Brien-Fleming",
#                        proc == "Bayes" ~ "Bayes",
#                        proc == "ISP" ~ "ISP")) %>% 
#rename(Procedure = proc)

tiff(file="figures/figure5.tiff",width=1800,height=1500, units = "px", res = 300)
ggplot(data = E_n, mapping = aes(x = d_actual, y = E_n, group = proc, linetype = proc, color = proc)) +
  geom_line(size = 1.2) +
  annotate(geom = "text", x = 0.3, y = 47.7, label = "Bayes", color = "red") +
  #annotate(geom = "segment", x = 0.54, y = 48, xend = 0.66, yend = 48, linetype = "longdash", color = "red", size = 1) +
  annotate(geom = "text", x = 0.85, y = 50, label = "Fixed") +
  annotate(geom = "text", x = 0.87, y = 35, label = "O'Brien \n Fleming", color = "green") +
  annotate(geom = "text", x = 0.23, y = 43.2, label = "ISP", color = "blue") +
  annotate(geom = "text", x = 0.75, y = 28, label = "Pocock", color = "grey") +
  #annotate(geom = "text", x = -0.1, y = 47, label = "mSPRT", color = "purple") + 
  theme_bw() +
  scale_linetype_manual(values=c("solid", "twodash", "solid", "longdash", "dotted")) +
  scale_color_manual(values=c("black", "green", "grey",  "red", "blue"))+
  labs(x = "True Effect Size (Cohen's d)", y = "Average Sample Size") +
  scale_x_continuous(breaks = seq(-0.2, 1, 0.2),
                     limits = c(-0.2, 1)) +
  theme(legend.position = "none")
dev.off()


#t3 = E_n %>% mutate(E_n = E_n %>% round()) %>% 
#  pivot_wider(names_from = d_actual, names_prefix = "d = ", values_from = E_n) 

#knitr::kable(t3)

#========================================================================#
########################## FIGURE 5: Efficiency ########################## 
#========================================================================#
# WITH MSPRT 

#facet.label = c("Fixed Sample Hypothesis Test", "Sequential Bayes Factor", 
#                "Independent Segments Procedure", "Group-Sequential Procedure")
#names(facet.label) = c(1, 2, 3, 4)

E_n = df %>% group_by(proc, d_actual) %>% 
  #mutate(proc = fct_relevel(proc, "Fixed", after = 0)) %>% 
  summarize(E_n = mean(n_cumulative)) # %>% 
#mutate(proc = case_when(proc == "Fixed" ~ "Fixed",
#                        proc == "asP" ~ "Pocock",
#                        proc == "asOF" ~ "O'Brien-Fleming",
#                        proc == "Bayes" ~ "Bayes",
#                        proc == "ISP" ~ "ISP")) %>% 
#rename(Procedure = proc)

# tiff(file="figures/figure5.tiff",width=1800,height=1500, units = "px", res = 300)
# ggplot(data = E_n, mapping = aes(x = d_actual, y = E_n, group = proc, linetype = proc, color = proc)) +
#   geom_line(size = 1.2) +
#   annotate(geom = "text", x = 0.3, y = 47.7, label = "Bayes", color = "red") +
#   #annotate(geom = "segment", x = 0.54, y = 48, xend = 0.66, yend = 48, linetype = "longdash", color = "red", size = 1) +
#   annotate(geom = "text", x = 0.85, y = 50, label = "Fixed") +
#   annotate(geom = "text", x = 0.65, y = 41, label = "O'Brien \n Fleming", color = "green") +
#   annotate(geom = "text", x = -0.11, y = 28.2, label = "ISP", color = "blue") +
#   annotate(geom = "text", x = 0.75, y = 28, label = "Pocock", color = "grey") +
#   annotate(geom = "text", x = -0.1, y = 47, label = "mSPRT", color = "purple") + 
#   theme_bw() +
#   scale_linetype_manual(values=c("solid", "twodash", "solid", "longdash", "dotted", "dotdash")) +
#   scale_color_manual(values=c("black", "green", "grey",  "red", "blue", "purple"))+
#   labs(x = "True Effect Size (Cohen's d)", y = "Average Sample Size") +
#   scale_x_continuous(breaks = seq(-0.2, 1, 0.2),
#                      limits = c(-0.2, 1)) +
#   theme(legend.position = "none")
# dev.off()


#t3 = E_n %>% mutate(E_n = E_n %>% round()) %>% 
#  pivot_wider(names_from = d_actual, names_prefix = "d = ", values_from = E_n) 

#knitr::kable(t3)

#=================================================================================#
#=================================================================================#
########################## FIGURE 6: EFFECT SIZE DENSITY ########################## 
#=================================================================================#
# Unconditional Bias
# Density of Obtained Effect Size Estimates
# Excepted effect size ---> d = 0.5
# True effect size == Hypothesized effect size (i.e., effect size powered for)
# Add legend for corrected ES estimate 

median_ES = df %>% 
  filter(d_actual == d_forpower) %>% 
  group_by(facet) %>% 
  summarize(ES_corrected = median(ES_corrected)) %>% 
  pull(ES_corrected) %>% round(2)


median_ES_uncorrected = df %>% 
  filter(d_actual == d_forpower) %>% 
  group_by(facet) %>% 
  summarize(ES = median(ES)) %>% 
  pull(ES) %>% round(2)

text_uncorrected = data.frame(facet = 1:5,
                              label = paste0("Median ES = ", median_ES_uncorrected),
                              x = -0.5, y = 2.35)

text_corrected = data.frame(facet = 1:5,
                  label = paste0("Corrected ES = ", median_ES),
                  x = -0.5, y = 2.1) %>% 
  filter(facet %in% c(2, 4, 5))

density = function(dat) {
  ggplot(data = dat, mapping = aes(x = ES)) +
    geom_density() +
    geom_density(mapping = aes(x = ES_corrected), linetype = "dashed") + 
    facet_wrap(vars(facet), labeller = labeller(facet = facet.label)) +
    theme_bw() +
    labs(x = "Empirical Effect Size Estimate", y = "Density") +
    scale_x_continuous(limits = c(-0.5, 1.5)) +
    #scale_y_continuous(limits=c(0,3)) +
    theme(legend.position = "none", strip.background = element_rect(fill = "white")) +
    geom_text(data = text_corrected, 
              aes(x = x, y = y, label = label, hjust = 0), 
              size = 2.7) +
    geom_text(data = text_uncorrected, 
              aes(x = x, y = y, label = label, hjust = 0), 
              size = 2.7)
}

tiff(file="figures/figure6.tiff",width=2350,height=1200, units = "px", res = 300)
density(df %>% filter(d_forpower == d_actual))
dev.off()

#========================================================================================#
########################## FIGURE 7: BIAS ACROSS UNEXPECTED D'S ########################## 
#========================================================================================#
# How does the procedure perform when d is unexpected?

df_summary = df %>% 
  group_by(proc, facet, d_actual) %>% 
  summarize(
    mse = mean((ES_corrected - d_actual) ^2),
    var = mean((ES_corrected - mean(ES_corrected))^2),
    median.ES = median(ES_corrected),
    ES = mean(ES_corrected)
    ) %>%
  mutate(
    bias = ES - d_actual,
    median.bias = median.ES - d_actual,
    bias.sq = (ES - d_actual)^2,
    median.bias.sq = (median.ES - d_actual)^2
    )

# Bias
tiff(file="figures/figure7a.tiff",width=2300,height=1200, units = "px", res = 300)
ggplot(data = df_summary, 
       mapping = aes(x = d_actual, y = bias)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
  geom_segment(aes(x = d_actual, xend = d_actual, y = 0, yend = bias),
               color = "red", linetype = "dashed") + 
  geom_point() +
  facet_wrap(vars(facet), labeller = labeller(facet = facet.label)) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1),
                     limits = c(-0.2, 1)) +
  labs(x = "True Effect Size (Cohen's d)", y = "Mean Estimated Effect Size - True Effect Size")
dev.off()

# Median bias
tiff(file="figures/figure7b.tiff",width=2300,height=1200, units = "px", res = 300)
ggplot(data = df_summary, 
       mapping = aes(x = d_actual, y = median.bias)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
  geom_segment(aes(x = d_actual, xend = d_actual, y = 0, yend = median.bias),
               color = "red", linetype = "dashed") + 
  geom_point() +
  facet_wrap(vars(facet), labeller = labeller(facet = facet.label)) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1),
                     limits = c(-0.2, 1)) +
  labs(x = "True Effect Size (Cohen's d)", y = "Median Estimated Effect Size - True Effect Size")
dev.off()

#===========================================================================#
########################## FIGURE 8: BIAS-VARIANCE ########################## 
#===========================================================================#
#df_summary = df %>% 
#  group_by(proc, facet, line, d_actual) %>% 
#  summarize(mse = mean((ES - d_actual) ^2),
#            var = mean((ES - mean(ES))^2),
#            ES = mean(ES)) %>% 
#  mutate(bias = ES - d_actual,
#         bias.sq = (ES - d_actual)^2)

#df_summary = df_summary %>% 
  #mutate(proc = fct_relevel(proc, "Fixed", after = 0)) %>% 
#  filter(proc != "asOF")

#facet.label[4] = "Group-Sequential (Pocock)"
facet1 = data.frame(
  facet = 1,
  label1 = "Mean Squared Error",
  label2 = "Variance",
  label3 = "Squared Bias",
  x1 = -0.2, 
  x2 = -0.04, 
  y1 = .13, 
  y2 = .11, 
  y3 = .09
  )

tiff(file="figures/figure8.tiff",width=2300,height=1300, units = "px", res = 300)
ggplot(data = df_summary, mapping = aes(x = d_actual, y = bias.sq)) +
  facet_wrap(vars(facet), labeller = labeller(facet = facet.label)) +
  theme_bw() +
  #geom_line(mapping = aes(y = median.bias.sq), color = "grey", linetype = "dotdash") + 
  geom_smooth(color = "black", se = F, size = .6) +
  geom_smooth(mapping= aes(y = mse), se = F, size = .6,
              color = "red", linetype = "dashed") + #MSE
  geom_smooth(mapping = aes(y = var), se = F, size = .6, linetype = "dotted") + #Variance
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  labs(x = "True Effect Size (Cohen's d)", y = "Bias-Variance") +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1),
                     limits = c(-0.2, 1)) +
  geom_segment(data = facet1, 
               aes(x = x1, xend = x2, y = y1, yend = y1), 
               linetype = "dashed", 
               color = "red") +
  geom_segment(data = facet1, 
               aes(x = x1, xend = x2, y = y2, yend = y2), 
               linetype = "dotted", 
               color = "blue") +
  geom_segment(data = facet1, aes(x = x1, xend = x2, y = y3, yend = y3)) +
  geom_text(data = facet1, 
            aes(x = x2, y = y1, label = label1, hjust = 0), 
            size = 2.7, 
            color = "red") +
  geom_text(data = facet1, 
            aes(x = x2, y = y2, label = label2, hjust = 0), 
            size = 2.7, 
            color = "blue") +
  geom_text(data = facet1, 
            aes(x = x2, y = y3, label = label3, hjust = 0), 
            size = 2.7)
dev.off()

#==============================================================================================#
########################## FIGURE 9: MSE CONDITIONAL ON STOPPING TIME ########################## 
#==============================================================================================#
#Without MSPR
#Take this figure out of paper?

sub = df %>% filter(proc != "Fixed" & d_forpower == d_actual) %>% 
  group_by(proc, facet, segment, d_actual) %>% 
  summarize(mse = mean((ES_corrected - d_actual)^2)) %>% ungroup() %>% 
  select(proc, facet, segment, mse) %>% 
  bind_rows(
    df %>% filter(proc != "Fixed" & d_forpower == d_actual) %>%
      group_by(proc, facet) %>% 
      summarize(mse = mean((ES_corrected - d_actual)^2)) %>% 
      mutate(segment = rep(4))) #%>% 
#mutate(proc = fct_relevel(proc, "Bayes", "asP", "ISP"))

facet.label = c("Fixed Hypothesis \n and Equivalence Test", 
                "Sequential Bayes Factor", 
                "Independent Segments \n Procedure", 
                "Pocock-like GS design",
                "O'Brien-Fleming-like GS design")
names(facet.label) = c(1, 2, 3, 4, 5)

tiff(file="figures/figure9.tiff",width=2000,height=1500, units = "px", res = 300)
ggplot(data = sub, mapping = aes(x = as.factor(segment), y = mse)) +
  facet_wrap(vars(facet), labeller = labeller(facet = facet.label), nrow = 2) +
  geom_segment(aes(x = segment, xend = segment, y = 0, yend = mse),
               color = "red", linetype = "dashed") +
  geom_point() +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  labs(x = "Segment", y = "Mean Squared Error") +
  scale_x_discrete(breaks = c(1, 2, 3, 4),
                   labels = c(1, 2, 3, "Overall"))
dev.off()

#tiff(file="Figure10.tiff",width=2200,height=1200, units = "px", res = 300)
#ggplot(data = sub, mapping = aes(x = as.factor(segment), y = mse)) +
#  facet_wrap(vars(facet), labeller = labeller(facet = facet.label)) +
#  geom_bar(stat = "identity", fill = "lightgrey", color = "black", width = .8) +
#  theme_bw() +
#  theme(strip.background =element_rect(fill="white")) +
#  labs(x = "Segment", y = "Mean Squared Error") +
#  scale_x_discrete(breaks = c(1, 2, 3, 4),
#                   labels = c(1, 2, 3, "Overall"))
#dev.off()


#==============================================================================================#
########################## FIGURE 9: MSE CONDITIONAL ON STOPPING TIME ########################## 
#==============================================================================================#
# WITH MSPRT
# sub = df %>% filter(proc != "Fixed" & d_forpower == d_actual) %>% 
#   group_by(proc, facet, segment, d_actual) %>% 
#   summarize(mse = mean((ES_corrected - d_actual)^2)) %>% ungroup() %>% 
#   select(proc, facet, segment, mse) %>% 
#   bind_rows(
#     df %>% filter(proc != "Fixed" & d_forpower == d_actual) %>%
#       group_by(proc, facet) %>% 
#       summarize(mse = mean((ES_corrected - d_actual)^2)) %>% 
#       mutate(segment = rep(4))) #%>% 
#   #mutate(proc = fct_relevel(proc, "Bayes", "asP", "ISP"))
# 
# facet.label = c("Fixed Hypothesis \n and Equivalence Test", 
#                 "Sequential Bayes Factor", 
#                 "Independent Segments \n Procedure", 
#                 "Pocock-like GS design",
#                 "O'Brien-Fleming-like GS design", 
#                 "Modified Sequential \n Probability Ratio Test")
# names(facet.label) = c(1, 2, 3, 4, 5, 6)
# 
# tiff(file="figures/figure9.tiff",width=2900,height=750, units = "px", res = 300)
# ggplot(data = sub, mapping = aes(x = as.factor(segment), y = mse)) +
#   facet_wrap(vars(facet), labeller = labeller(facet = facet.label), nrow = 1) +
#   geom_segment(aes(x = segment, xend = segment, y = 0, yend = mse),
#                color = "red", linetype = "dashed") +
#   geom_point() +
#   theme_bw() +
#   theme(strip.background =element_rect(fill="white")) +
#   labs(x = "Segment", y = "Mean Squared Error") +
#   scale_x_discrete(breaks = c(1, 2, 3, 4),
#                    labels = c(1, 2, 3, "Overall"))
# dev.off()

#tiff(file="Figure10.tiff",width=2200,height=1200, units = "px", res = 300)
#ggplot(data = sub, mapping = aes(x = as.factor(segment), y = mse)) +
#  facet_wrap(vars(facet), labeller = labeller(facet = facet.label)) +
#  geom_bar(stat = "identity", fill = "lightgrey", color = "black", width = .8) +
#  theme_bw() +
#  theme(strip.background =element_rect(fill="white")) +
#  labs(x = "Segment", y = "Mean Squared Error") +
#  scale_x_discrete(breaks = c(1, 2, 3, 4),
#                   labels = c(1, 2, 3, "Overall"))
#dev.off()

#==========================================================================================#
########################## FIGUREXI: MEDIAN EFFECT SIZE ESTIIMATE ########################## 
#==========================================================================================#
med = df %>% 
  group_by(proc, facet, d_actual) %>% 
  summarize(medianES = median(ES_corrected)) 

tiff(file="figures/figurexi.tiff",width=2000,height=1400, units = "px", res = 300)
ggplot(data = med, mapping = aes(x = d_actual, y = medianES)) +
  geom_point() +
  facet_wrap(vars(facet), labeller = labeller(facet = facet.label)) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  theme(legend.position = "none", strip.background = element_rect(fill = "white")) +
  labs(x = "True Effect Size (Cohen's d)", y = "Median Effect Size Estimate") +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1),
                     limits = c(-0.2, 1)) 
dev.off()

#t4 = med %>% mutate(medianES = medianES %>% round(2)) %>%  
#  pivot_wider(names_from = d_actual, names_prefix = "*d* = ", values_from = medianES)

#==============================================================================================#
################# FIGUREXII: EFFECT SIZE DENSITY CONDITIONAL ON REJECTION ###################### 
#==============================================================================================#
# Conditional Bias (Bias conditional on Rejection of the Null)
# Density of Obtained Significant Effects
# True effect size == Hypothesized effect size 
# Taking this fig out of the paper

#facet4 = facet4 %>% mutate(y1 = 2.5, y2 = 2)
#tiff(file="figurexii.tiff",width=1500,height=1200, units = "px", res = 300)
#density(df %>% filter(d_forpower == d_actual & decision == "reject.null"))
#dev.off()

#========================================================================================#
##################### FIGUREXIII: META-ANALYTIC EFFECT SIZE ESTIMATE ##################### 
#========================================================================================#
#rm(list = ls())
#df = read.csv("simulations/data.csv")
# Recode this later
df = df %>% 
  mutate(proc = fct_relevel(proc, "Fixed", "ISP", "asP", "asOF", "Bayes"),
         n = ifelse(proc == "ISP", 25, n)) 

# set.seed(42)
# df = df %>% # meta-analysis on all 4 million observations exhausts vector memory
#   group_by(proc, d_forpower, d_actual, power) %>% 
#   slice_sample(n = 10000) # randomly select 400,000 obs (10,000 for each unique specification)

# Compute sampling variance for meta-analysis
models = df %>% split(list(.$proc, .$d_actual)) %>% 
  map(~ 1/.$n + 1/.$n + .$ES^2/(4*.$n)) # compute sampling variance

df = df %>% cbind(
  models %>% 
    tibble() %>% 
    unnest(cols = c(.)) %>% 
    rename(., "var" = .) %>% 
    select(var))

# Run meta-analysis
meta = df %>% split(list(.$proc, .$d_actual)) %>% 
  map(~ rma(yi = .$ES, vi = .$var, method = "FE")$b) #takes 5 mins to run
beep()
meta_analysis = data.frame(
  proc = rep(c("Fixed Sample Hypothesis Test", 
               "Independent Segments Procedure", 
               "Group-Sequential (Pocock)",
               "O'Brien-Fleming", 
               "Sequential Bayes Factor"), 8),
  d_actual = rep(c(-0.2, 0, 0.2, 0.4, 0.5, 0.6, 0.8, 1), each = 5),
  meta.ES = meta %>% 
    tibble() %>% 
    unnest(cols = c(.)) %>% 
    rename(., "meta" = .) %>% 
    select(meta)) %>% 
  mutate(proc = fct_relevel(
    proc, "Fixed Sample Hypothesis Test", 
    "Sequential Bayes Factor", 
    "Group-Sequential (Pocock)", 
    "O'Brien-Fleming", 
    "Independent Segments Procedure"))

tiff(file="figures/figurexiii.tiff",width=2900,height=1400, units = "px", res = 300)
ggplot(data = meta_analysis, mapping = aes(x = d_actual, y = meta)) +
  facet_wrap(vars(proc)) +
  geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed") +
  geom_point() +
  theme_bw() +
  labs(x = "True effect size (Cohen's d)", y = "Meta-analytic effect size estimate") +
  scale_x_continuous(limits = c(-0.2, 1),
                     breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_y_continuous(breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme(strip.background = element_rect(fill = "white"))
dev.off()

# For Bayes, I use the mean of the posterior distribution as the ES estimator, 
# which shrinks the estimate in early terminations (which happen more often if the true ES is greater)
# ----> Explains why the SBF meta-analytic effect size estimate underestimates the true effect as the true population effect size increases
# No longer applies?^

## Figure 10
## *Meta-analytic Effect Size Estimates* 
#\noindent
#```{r, fig.heigth = 4.72, fig.width = 5}
#knitr::include_graphics("Figure9.tiff", dpi = 300)
#```
#==============================================================================================#
