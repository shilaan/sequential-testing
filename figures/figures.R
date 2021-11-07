# TO DO -------------------------------------------------------------------
# Include both bias-corrected and uncorrected estimates in the manuscript
# Write up final results
# Update Notes for figures
# Write translational abstract

# Set up ------------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(ggtext)
library(grid)
library(papaja)
library(gridExtra)
library(rpact)
library(BayesFactor)

# Read data ---------------------------------------------------------------
df <- read.csv("simulations/data.csv")

df <- df %>% 
  mutate(
    proc = fct_relevel(proc, "Fixed", after = 0),
    facet = case_when(
      proc == "Fixed" ~ 1,
      proc == "Bayes" ~ 2,
      proc == "ISP" ~ 3,
      proc == "asP" ~ 4,
      proc == "asOF"~ 5)
  )

facet.label <- c("Fixed Hypothesis and Equivalence Test", 
                "Sequential Bayes Factor", 
                "Independent Segments Procedure", 
                "Pocock-like GS design",
                "O'Brien-Fleming-like GS design")
names(facet.label) <- c(1, 2, 3, 4, 5)

#=============================================================================#
########################## FIGURE 1A: POCOCK EXAMPLE ########################## 
#=============================================================================#

gs <- function(proc, return = "n") {
  bs <- ifelse(proc == "asOF", "bsOF", "bsP")
  if(proc == "asOF") {rates = c(0.50, 0.75, 1)} else{rates = c(1/3, 2/3, 1)} #specify information rates for OBF vs Pocock
  
  
  # Specify the design
  design <- getDesignGroupSequential(
    sided = 1,
    alpha = 0.05, 
    beta = 0.2,
    kMax = 3, 
    typeOfDesign = proc,
    typeBetaSpending= bs,
    informationRates = rates,
    bindingFutility = TRUE
    ) 
  
  # Get parameters
  parameters <- getSampleSizeMeans(design = design, groups = 2, alternative = 0.5)
  n_gs = ceiling(c(parameters$numberOfSubjects[1], diff(parameters$numberOfSubjects))/2) # n per look per group
  alpha <- parameters$criticalValuesPValueScale
  futility <- parameters$futilityBoundsPValueScale
  
  if(return == "alpha") {return(alpha)}
  if(return == "n") {return(n_gs)}
  if(return == "futility") {return(futility)}
}

n_p <- gs(proc = "asP")[1]
alpha_p <- gs(proc = "asP", return = "alpha") #p-value scale
alpha_pd <- round(c( #critical d-value scale
  qt(1-alpha_p[1], df = 2*n_p-2) / sqrt(n_p/2),
  qt(1-alpha_p[2], df = 4*n_p-2) / sqrt(2*n_p/2),
  qt(1-alpha_p[3], df = 6*n_p-2) / sqrt(3*n_p/2)
  )
  ,2)

futility_p <- gs(proc = "asP", return = "futility") #futility on critical p-value scale
futility_pd <- round(
  c(qt(1 - futility_p[1], df = (n_p*2)-2)/sqrt(n_p/2), 
    qt(1 - futility_p[2], df = 4*n_p-2)/sqrt(2*n_p/2)), 
  2) #futility for pocock on critical d-value

n_of <- gs(proc = "asOF")
alpha_of <- gs(proc = "asOF", return = "alpha")#p-value scale
alpha_ofd <- round(c( #critical d-value scale
  qt(1-alpha_of[1], df = 2*n_of[1]-2) / sqrt(n_of[1]/2),
  qt(1-alpha_of[2], df = 2*(n_of[1]+n_of[2])-2) / sqrt((n_of[1]+n_of[2])/2),
  qt(1-alpha_of[3], df = 2*(sum(n_of))-2) / sqrt(sum(n_of)/2)
  ),2)

futility_of <- gs(proc = "asOF", return = "futility")  #futility on critical p-value scale
futility_ofd <- round(
  c(qt(1 - futility_of[1], df = (n_of[1]*2)-2)/sqrt(n_of[1]/2), 
    qt(1 - futility_of[2], df = 2*(n_of[1]+n_of[2])-2)/sqrt((n_of[1]+n_of[2])/2)), 
  2) #futility for OBF on critical d-value scale

d.plot <- function(n, alpha, t) {
  df <- 2*n-2
  critical_d <- round(qt(1-alpha, df = 2*n-2) * sqrt(2*(1/n)),2)
  p <- ggplot(data.frame(x = c(-4, 4)), aes(x)) +
    stat_function(
      fun = dt, 
      args = list(df = df)
      ) +
    stat_function(
      fun = dt, 
      args = list(df = df), 
      xlim = c(qt(p = 1 - alpha, df = df), 4),
      geom = "area", 
      fill = "red", 
      alpha = 0.3
      ) +
    labs(title = t, y = "Density") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      plot.subtitle = ggtext::element_markdown(size = 13, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), 
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.text.x.top = ggtext::element_markdown(),
      axis.text.x.bottom = ggtext::element_markdown(),
      panel.grid = element_blank()
      ) + 
    geom_vline(xintercept = qt(1-alpha, df = df), linetype = 2, color = "grey") 
  return(p)
}


p1 <- d.plot(n = n_p, alpha = alpha_p[1], t = "Pocock (First Look)") +
  scale_x_continuous(
    limits = c(-0.44, 1*sqrt(n_p/2)), #create secondary p-value axis
    breaks = c(futility_pd[1]*sqrt(n_p/2), alpha_pd[1]*sqrt(n_p/2)),
    labels = paste0("*d* = ", round(c(futility_pd[1], alpha_pd[1]), 2)),
    sec.axis = dup_axis(
      labels = paste0("*p* = ", round(c(futility_p[1], alpha_p[1]), 2)) 
    )
  ) +
  annotate(
    geom = "text", 
    x = 1.525, 
    y = 0.27, 
    fontface = "bold", 
    size = 2.2,
    label = paste0(
      "CONTINUE \n", 
      futility_pd[1],
      " < d < ",
      alpha_pd[1],
      "\n",
      signif(alpha_p[1], 1),
      " < p < ",
      signif(futility_p[1], 2))
    ) +
  annotate(
    geom = "text", 
    x = 3, 
    y = 0.27, 
    label = paste0(
      "REJECT H0 \n d > ", 
      alpha_pd[1],
      "\n p < ", 
      signif(alpha_p[1], 1)
      ), 
    fontface = "bold", 
    size = 2.2
    ) +
  geom_vline(
    xintercept = futility_pd[1] *sqrt(n_p/2), 
    linetype = "dashed", 
    color = "grey"
    ) +
  stat_function(
    fun = dt, 
    args = list(df = 2*n_p -2), 
    xlim = c(-0.44, qt(1 - futility_p[1], df = 2*n_p -2)),
    geom = "area", 
    fill = "red", 
    alpha = 0.3
    ) +
  stat_function(
    fun = dt, 
    args = list(df = 2*n_p -2), 
    xlim = c(
      qt(1 - futility_p[1], df = 2*n_p -2), 
      qt(p = 1 - alpha_p[1], df = 2*n_p-2)
      ),
    geom = "area", 
    fill = "grey", 
    alpha = 0.3
  ) +
  annotate(
    geom = "text", 
    x = -0.035, 
    y = 0.27, 
    label = paste0(
      "REJECT H1 \n d < ",
      futility_pd[1],
      "\n p > ", 
      round(futility_p[1], 2)
      ), 
    fontface = "bold", 
    size = 2.2
    ) +
  annotate(
    geom = "text", 
    x = 3.2, 
    y = 0.07,   
    label = paste0("N = ", n_p), 
    fontface = "italic", 
    size = 2.2)

p2 <- d.plot(
  n = 2*n_p, 
  alpha = alpha_p[2], 
  t = "Pocock (Second Look)"
  ) +
  scale_x_continuous(
    limits = c(0.02*sqrt(2*n_p/2), 0.6*sqrt(2*n_p/2)),
    breaks = c(futility_pd[2]*sqrt(2*n_p/2), alpha_pd[2]*sqrt(2*n_p/2)),
    labels = paste0("*d* = ", round(c(futility_pd[2], alpha_pd[2]), 2)),
    sec.axis = dup_axis(
      labels = paste0("*p* = ", round(c(futility_p[2], alpha_p[2]), 2)) 
    )
  ) +
  geom_vline(
    xintercept = futility_pd[2] *sqrt(2*n_p/2), 
    linetype = "dashed", 
    color = "grey"
    ) +
  stat_function(
    fun = dt, 
    args = list(df = 4*n_p -2), 
    xlim = c(0, qt(1 - futility_p[2], df = 4*n_p -2)),
    geom = "area", 
    fill = "red", 
    alpha = 0.3
    ) +
  stat_function(
    fun = dt, 
    args = list(df = 2*n_p -2), 
    xlim = c(
      qt(1 - futility_p[2], df = 4*n_p -2), 
      qt(p = 1 - alpha_p[2], df = 4*n_p-2)
    ),
    geom = "area", 
    fill = "grey", 
    alpha = 0.3
  ) +
  annotate(
    geom = "text", 
    x = 1.6, 
    y = 0.27, 
    fontface = "bold", 
    size = 2.2,
    label = paste0(
      "CONTINUE \n", 
      futility_pd[2],
      " < d < ",
      alpha_pd[2],
      "\n",
      signif(alpha_p[2], 1),
      " < p < ",
      signif(futility_p[2], 2))
    ) +
  annotate(
    geom = "text", 
    x = 2.6, 
    y = 0.27, 
    label = paste0(
      "REJECT H0 \n d > ",
      alpha_pd[2],
      "\n p < ", 
      signif(alpha_p[2], 1)
      ),
    fontface = "bold", 
    size = 2.2
    ) +
  annotate(
    geom = "text", 
    x = 0.41, 
    y = 0.27,  
    label = paste0(
      "REJECT H1 \n d < ",
      futility_pd[2],
      "\n p > ",
      round(futility_p[2], 2)
      ), 
    fontface = "bold", 
    size = 2.2
    ) +
  annotate(
    geom = "text", 
    x = 2.75, 
    y = 0.07,
    label = paste0("N = ", 2*n_p), 
    fontface = "italic", 
    size = 2.2
    )

p3 <- d.plot(
  n = 3*n_p, 
  alpha = alpha_p[3], 
  t = "Pocock (Third Look)") +
  scale_x_continuous(
    limits = c(0, 0.6*sqrt(3*n_p/2)),
    breaks = alpha_pd[3]*sqrt(3*n_p/2),
    labels = paste0("*d* = ", round(alpha_pd[3], 2)),
    sec.axis = dup_axis(
      labels = paste0("*p* = ", round(alpha_p[3], 2))
      )
    ) +
  stat_function(
    fun = dt, 
    args = list(df = 6*n_p -2), 
    xlim = c(0, qt(1 - alpha_p[3], df = 6*n_p-2)),
    geom = "area", 
    fill = "grey", 
    alpha = 0.3
    ) +
  annotate(
    geom = "text", 
    x = 0.23, 
    y = 0.27,  
    label = "FAIL \n TO \n REJECT", 
    fontface = "bold", 
    size = 2.2
    ) +
  annotate(
    geom = "text",
    x = 3, 
    y = 0.27, 
    label = paste0(
      "REJECT H0 \n d > ",
      alpha_pd[3],
      "\n p < ", 
      signif(alpha_p[3], 1)
      ), 
    fontface = "bold", 
    size = 2.2
    )  +
  annotate(
    geom = "text",
    x = 3.3, 
    y = 0.07, 
    label = paste0("N = ", 3*n_p), 
    fontface = "italic",
    size = 2.2
    )

tiff(file="figures/figure1a.tiff",width=2500,height=800, units = "px", res = 300)
grid.arrange(p1, p2, p3, nrow = 1, 
             left = textGrob("Density under H0: \u03b4 = 0", 
                             rot = 90, hjust = 0.57,
                             gp = gpar(fontsize = 10)))
dev.off()


#======================================================================================#
########################## FIGURE 1B: O'BRIEN-FLEMING EXAMPLE ########################## 
#======================================================================================#

of1 <- d.plot(n = n_of[1], alpha = alpha_of[1], t = "O'Brien-Fleming (First Look)") +
  scale_x_continuous(
    limits = c(-0.46, 1*sqrt(n_of[1]/2)), 
    breaks = c(futility_ofd[1]*sqrt(n_of[1]/2), alpha_ofd[1]*sqrt(n_of[1]/2)),
    labels = paste0("*d* = ", round(c(futility_ofd[1], alpha_ofd[1]), 2)),
    sec.axis = dup_axis( 
      labels = paste0("*p* = ", round(c(futility_of[1], alpha_of[1]), 2)) 
    ) #create secondary p-value axis
  ) +
  annotate(
    geom = "text", 
    x = 1.7, 
    y = 0.27, 
    fontface = "bold", 
    size = 2.2,
    label = paste0(
      "CONTINUE \n", 
      futility_ofd[1],
      " < d < ",
      alpha_ofd[1],
      "\n",
      signif(alpha_of[1], 1),
      " < p < ",
      signif(futility_of[1], 2))
  ) +
  annotate(
    geom = "text", 
    x = 3.2, 
    y = 0.27, 
    label = paste0(
      "REJECT H0 \n d > ", 
      alpha_ofd[1],
      "\n p < ", 
      signif(alpha_of[1], 1)
      ), 
    fontface = "bold", 
    size = 2.2
  ) +
  geom_vline(
    xintercept = futility_ofd[1] *sqrt(n_of[1]/2), 
    linetype = "dashed", 
    color = "grey"
  ) +
  stat_function(
    fun = dt, 
    args = list(df = 2*n_of[1]-2), 
    xlim = c(-0.46, qt(1 - futility_of[1], df = 2*n_of[1] -2)),
    geom = "area", 
    fill = "red", 
    alpha = 0.3
  ) +
  stat_function(
    fun = dt, 
    args = list(df = 2*n_of[1] -2), 
    xlim = c(
      qt(p = 1 - futility_of[1], df = 2*n_of[1] -2), 
      qt(p = 1 - alpha_of[1], df = 2*n_of[1] -2)
    ),
    geom = "area", 
    fill = "grey", 
    alpha = 0.3
  ) +
  annotate(
    geom = "text", 
    x = 0.005, 
    y = 0.27, 
    label = paste0(
      "REJECT H1 \n d < ",
      futility_ofd[1],
      "\n p > ", 
      round(futility_of[1], 2)
    ), 
    fontface = "bold", 
    size = 2.2
  ) +
  annotate(
    geom = "text", 
    x = 3.2, 
    y = 0.07,   
    label = paste0("N = ", n_of[1]), 
    fontface = "italic", 
    size = 2.2)

of2 <- d.plot(n = n_of[1] + n_of[2], alpha = alpha_of[2], t = "O'Brien-Fleming (Second Look)") +
  scale_x_continuous(
    limits = c(-0.15, 0.8*sqrt((n_of[1]+n_of[2])/2)), 
    breaks = c(futility_ofd[2]*sqrt((n_of[1]+n_of[2])/2), alpha_ofd[2]*sqrt((n_of[1]+n_of[2])/2)),
    labels = paste0("*d* = ", round(c(futility_ofd[2], alpha_ofd[2]), 2)),
    sec.axis = dup_axis( 
      labels = paste0("*p* = ", round(c(futility_of[2], alpha_of[2]), 2)) 
    ) #create secondary p-value axis
  ) +
  annotate(
    geom = "text", 
    x = 1.56, 
    y = 0.27, 
    fontface = "bold", 
    size = 2.2,
    label = paste0(
      "CONTINUE \n", 
      futility_ofd[2],
      " < d < ",
      alpha_ofd[2],
      "\n",
      signif(alpha_of[2], 1),
      " < p < ",
      signif(futility_of[2], 2))
  ) +
  annotate(
    geom = "text", 
    x = 3, 
    y = 0.27, 
    label = paste0(
      "REJECT H0 \n d > ", 
      alpha_ofd[2],
      "\n p < ", 
      signif(alpha_of[2], 1)
      ), 
    fontface = "bold", 
    size = 2.2
  ) +
  geom_vline(
    xintercept = futility_ofd[2] *sqrt((n_of[1]+n_of[2])/2), 
    linetype = "dashed", 
    color = "grey"
  ) +
  stat_function(
    fun = dt, 
    args = list(df = 2*(n_of[1]+n_of[2])-2), 
    xlim = c(-0.15, qt(1 - futility_of[2], df = 2*(n_of[1]+n_of[2])-2)),
    geom = "area", 
    fill = "red", 
    alpha = 0.3
  ) +
  stat_function(
    fun = dt, 
    args = list(df = 2*(n_of[1]+n_of[2])-2), 
    xlim = c(
      qt(p = 1 - futility_of[2], df = 2*(n_of[1]+n_of[2])-2), 
      qt(p = 1 - alpha_of[2], df = 2*(n_of[1]+n_of[2])-2)
    ),
    geom = "area", 
    fill = "grey", 
    alpha = 0.3
  ) +
  annotate(
    geom = "text", 
    x = 0.255, 
    y = 0.27, 
    label = paste0(
      "REJECT H1 \n d < ",
      futility_ofd[2],
      "\n p > ", 
      round(futility_of[2], 2)
    ), 
    fontface = "bold", 
    size = 2.2
  ) +
  annotate(
    geom = "text", 
    x = 3.2, 
    y = 0.07,   
    label = paste0("N = ", n_of[1]+n_of[2]), 
    fontface = "italic", 
    size = 2.2)

of3 <- d.plot(
  n = sum(n_of), 
  alpha = alpha_of[3], 
  t = "O'Brien-Fleming (Third Look)") +
  scale_x_continuous(
    limits = c(0, 0.6*sqrt(sum(n_of)/2)),
    breaks = alpha_ofd[3]*sqrt(sum(n_of)/2),
    labels = paste0("*d* = ", round(alpha_ofd[3], 2)),
    sec.axis = dup_axis(
      labels = paste0("*p* = ", round(alpha_of[3], 2))
    )
  ) +
  stat_function(
    fun = dt, 
    args = list(df = 2*sum(n_of)-2), 
    xlim = c(0, qt(1 - alpha_of[3], df = 2*sum(n_of)-2)),
    geom = "area", 
    fill = "grey", 
    alpha = 0.3
  ) +
  annotate(
    geom = "text", 
    x = 0.25, 
    y = 0.27,  
    label = "FAIL \n TO \n REJECT", 
    fontface = "bold", 
    size = 2.2
  ) +
  annotate(
    geom = "text",
    x = 2.8, 
    y = 0.27, 
    label = paste0(
      "REJECT H0 \n d > ",
      alpha_ofd[3],
      "\n p < ", 
      signif(alpha_of[3], 1)
    ), 
    fontface = "bold", 
    size = 2.2
  )  +
  annotate(
    geom = "text",
    x = 3, 
    y = 0.07, 
    label = paste0("N = ", sum(n_of)), 
    fontface = "italic",
    size = 2.2
  )


tiff(file="figures/figure1b.tiff",width=2500,height=800, units = "px", res = 300)
grid.arrange(of1, of2, of3, nrow = 1,
             left = textGrob("Density under H0: \u03b4 = 0", rot = 90, hjust = 0.57, 
                             gp = gpar(fontsize = 10)))
dev.off()

#============================================================================#
########################## FIGURE 2AB: BAYES EXAMPLE ########################## 
#============================================================================#
# Adapted from https://www.shinyapps.org/apps/RGraphCompendium/index.php#prior-and-posterior
# Big thank you to Eric-Jan Wagenmakers & Quentin Gronau for maintaining this site!
# I took the code from the site above and recoded it from base R to ggplot

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
  label1 <- sprintf("BF[1] == '%0.2f'", BF10)
  label2 <- sprintf("BF[0] == '%0.2f'", BF01)

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
    annotate("text", x = xticks[2], y=ylim[2] - 0.3, size = 6,
             label=label1, parse = TRUE) +
    annotate("text", x = xticks[2], y=ylim[2], size = 6,
             label=label2, parse = TRUE) +
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
    stat_function(
      fun = dt,
      args = list(df = df)
      ) +
    stat_function(
      fun = dt, 
      args = list(df = df), 
      xlim = c(qt(p = 1 - alpha, df = df), 4),
      geom = "area", 
      fill = "red", 
      alpha = 0.6
      ) +
    labs(
      title = t, 
      y = "Density"
      ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      plot.subtitle = ggtext::element_markdown(size = 13, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), 
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 11),
      axis.text.x.top = ggtext::element_markdown(),
      axis.text.x.bottom = ggtext::element_markdown(),
      panel.grid = element_blank()
      ) + 
    scale_x_continuous(
      limits = c(0, .8*sqrt(n/2)),
      breaks = qt(1 - alpha, df = df),
      labels = paste0("*d* = ", critical_d),
      sec.axis = dup_axis(
        labels = paste0("*p* = ", round(alpha, 2))
      )
    ) + 
    geom_vline(
      xintercept = qt(1-alpha, df = df), 
      linetype = 2, 
      color = "grey"
      ) +
    stat_function(
      fun = dt, 
      args = list(df = df), 
      xlim = c(0, qt(1 - ifelse(proc == "Fixed", alpha, alpha_weak), df = df)),
      geom = "area", 
      fill = "grey", 
      alpha = 0.3
      ) +
    annotate(
      geom = "text", 
      x = 0.29, 
      y = 0.27,  
      label = "FAIL \n TO \n REJECT",
      fontface = "bold", 
      size = 2.3
      ) +
    annotate(
      geom = "text", 
      x = 2.5, 
      y = 0.27,  
      label = paste0(
        "REJECT \n d > ", 
        critical_d,
        "\n p < ", 
        round(alpha, 2)
        ), 
      fontface = "bold", 
      size = 2.3
      )
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
  geom_vline(
    xintercept = qt(1-alpha_weak, df = 2*n_s-2), 
    linetype = 2, 
    color = "grey") +
  annotate(geom = "text", 
           x = 1.5, 
           y = 0.27,   
           size = 2.3,
           label = paste0(
             "CONTINUE \n ", 
             round(qt(1-alpha_weak, df = 2*n_s-2) / sqrt(n_s/2),2),
             " < d < ", 
             round(qt(1 - alpha_strong, df = 2*n_s-2) / sqrt(n_s/2), 2), "\n",
             alpha_strong, 
             " < p < ", 
             round(alpha_weak, 2)
             ),
           fontface = "bold") +
  annotate(geom = "text", 
           x = 2.65, 
           y = 0.07,   
           label = paste0("N = ", n_s), 
           fontface = "italic", 
           size = 2.3) + 
  scale_x_continuous(
    limits = c(0, .8*sqrt(n_s/2)),
    breaks = c(qt(1-alpha_weak, df = 2*n_s-2), qt(1 - alpha_strong, df = 2*n_s-2)),
    labels = paste0("*d* = ", c(round(qt(1-alpha_weak, df = 2*n_s-2) * sqrt(2*(1/n_s)),2), 
                                round(qt(1-alpha_strong, df = 2*n_s-2) * sqrt(2*(1/n_s)),2))),
    sec.axis = dup_axis(
      labels = paste0("*p* = ", c(round(alpha_weak, 2), round(alpha_strong, 3)))
    )
  ) 

p3 = critical.d(proc = "ISP", segment = 3) +
  annotate(geom = "text", 
           x = 2.65, 
           y = 0.07,   
           label = paste0("N = ", n_s), 
           fontface = "italic", 
           size = 2.3)

tiff(file="figures/figure3b.tiff",width=2500,height=800, units = "px", res = 300)
grid.arrange(p1, p2, p3, nrow = 1, 
             left = textGrob("Density under H0: \u03b4 = 0", rot = 90, hjust = 0.57,
                             gp = gpar(fontsize = 10)))
dev.off()
#=========================================================================#
########################## FIGURE 4: ERROR RATES ########################## 
#=========================================================================#

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

#========================================================================#
########################## FIGURE 5: Efficiency ########################## 
#========================================================================#

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

# t3 = E_n %>% mutate(E_n = E_n %>% round()) %>% 
#   pivot_wider(names_from = d_actual, names_prefix = "d = ", values_from = E_n) 
# 
# knitr::kable(t3)


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

text_uncorrected = data.frame(
  facet = 1:5,
  label = paste0("Median ES = ", median_ES_uncorrected),
  x = -0.5, 
  y = 2.35
  )

text_corrected = data.frame(
  facet = 1:5,
  label = paste0("Corrected ES = ", median_ES),
  x = -0.5, 
  y = 2.1
  ) %>% 
  filter(facet %in% c(2, 4, 5))

legend = data.frame(
  facet = rep(1),
  label = c("Uncorrected", "Corrected"),
  y = c(2.35, 2.1),
  x = rep(1)
)

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
              size = 2.7) +
    geom_text(data = legend,
              aes(x = x, y = y, label = label, hjust = 0),
              size = 2.7) +
    geom_line(data = data.frame(x = c(0.7, .95),
                                y = rep(2.35),
                                facet = 1),
              mapping = aes(x = x, y = y)) +
    geom_line(data = data.frame(x = c(0.7, .95),
                                y = rep(2.1),
                                facet = 1),
              mapping = aes(x = x, y = y),
              linetype = "dashed")
              
}


tiff(file="figures/figure6.tiff",width=2350,height=1200, units = "px", res = 300)
density(df %>% filter(d_forpower == d_actual))
dev.off()

#========================================================================================#
########################## FIGURE 7: BIAS ACROSS UNEXPECTED D'S ########################## 
#========================================================================================#
# How does the procedure perform when d is unexpected?

df_long <- df %>% 
  pivot_longer(
    cols = c(ES, ES_corrected),
    names_to = c("ES_type"),
    values_to = "ES"
  ) 

df_summary <- df_long %>% 
  group_by(proc, facet, d_actual, ES_type) %>% 
  summarize(
    mse = mean((ES - d_actual) ^2),
    var = mean((ES - mean(ES))^2),
    median.ES = median(ES),
    ES = mean(ES)
  ) %>% 
  mutate(
    bias = ES - d_actual,
    median.bias = median.ES - d_actual,
    bias.sq = (ES - d_actual)^2,
    median.bias.sq = (median.ES - d_actual)^2
  ) %>% 
  filter(!((ES_type == "ES_corrected") & (proc == "Fixed" | proc == "ISP"))) #%>% 
  

# Mean Bias

legend <- data.frame(
  facet = rep(1),
  label = c("Uncorrected", "Corrected"),
  y = c(0.08, 0.05),
  x = rep(0.1),
  ES_type = c("ES", "ES_corrected")
)

tiff(file="figures/figure7a.tiff",width=2300,height=1200, units = "px", res = 300)
ggplot(data = df_summary, 
       mapping = aes(x = d_actual, y = bias, color = ES_type)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
  geom_segment(
    aes(
      x = d_actual, 
      xend = d_actual, 
      y = 0, 
      yend = bias, 
      color = ES_type
      ),
    linetype = "dashed"
    ) + 
  geom_point() +
  facet_wrap(vars(facet), labeller = labeller(facet = facet.label)) +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1),
    limits = c(-0.2, 1)
    ) +
  labs(x = "True Effect Size (Cohen's d)", 
       y = "Mean Estimated Effect Size - True Effect Size") +
  scale_color_manual(values = c("red", "blue")) +
  geom_text(data = legend, 
            aes(x = x, y = y, label = label, hjust = 0)) + 
  geom_line(
    data = data.frame(
      x = c(-0.2, 0.05),
      y = rep(0.08),
      facet = 1
    ),
    mapping = aes(x = x, y = y),
    color = "red",
    linetype = "dashed"
  ) +
  geom_line(
    data = data.frame(
      x = c(-0.2, 0.05),
      y = rep(0.05),
      facet = 1
    ),
    mapping = aes(x = x, y = y),
    color = "blue",
    linetype = "dashed"
  )
  
dev.off()

# Median bias

legend <- data.frame(
  facet = rep(1),
  label = c("Uncorrected", "Corrected"),
  y = c(0.11, 0.06),
  x = rep(0.1),
  ES_type = c("ES", "ES_corrected")
)

tiff(file="figures/figure7b.tiff",width=2300,height=1200, units = "px", res = 300)
ggplot(data = df_summary, 
       mapping = aes(x = d_actual, y = median.bias, color = ES_type)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
  geom_segment(
    aes(
      x = d_actual, 
      xend = d_actual, 
      y = 0, 
      yend = median.bias, 
      color = ES_type
    ),
    linetype = "dashed"
  ) + 
  geom_point() +
  facet_wrap(vars(facet), labeller = labeller(facet = facet.label)) +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1),
    limits = c(-0.2, 1)
  ) +
  labs(x = "True Effect Size (Cohen's d)", 
       y = "Median Estimated Effect Size - True Effect Size") +
  scale_color_manual(values = c("red", "blue")) +
  geom_text(data = legend, 
            aes(x = x, y = y, label = label, hjust = 0)) + 
  geom_line(
    data = data.frame(
      x = c(-0.2, 0.05),
      y = rep(0.11),
      facet = 1
    ),
    mapping = aes(x = x, y = y),
    color = "red",
    linetype = "dashed"
  ) +
  geom_line(
    data = data.frame(
      x = c(-0.2, 0.05),
      y = rep(0.06),
      facet = 1
    ),
    mapping = aes(x = x, y = y),
    color = "blue",
    linetype = "dashed"
  )
dev.off()


#===========================================================================#
########################## FIGURE 8: BIAS-VARIANCE ########################## 
#===========================================================================#
df_summary <- df_summary %>% 
  mutate(proc_ES = factor(paste0(proc, ES_type))) %>% 
  mutate(proc_ES = fct_relevel(proc_ES, c(
    "FixedES", 
    "BayesES", 
    "asPES", 
    "asOFES",
    "ISPES", 
    "BayesES_corrected", 
    "asPES_corrected", 
    "asOFES_corrected"
    )
  ))

proc_ES.label <- c("Fixed Hypothesis and \n Equivalence Test", 
                 "Sequential Bayes Factor \n Naive effect size estimate", 
                 "Pocock-like GS design  \n Naive effect size estimate",
                 "O'Brien-Fleming-like GS design  \n Naive effect size estimate",
                 "Independent Segments \n Procedure", 
                 "Sequential Bayes Factor \n Bias-adjusted estimator",
                 "Pocock-like GS design \n Bias-adjusted estimator",
                 "O'Brien-Fleming-like GS design \n Bias-adjusted estimator")
names(proc_ES.label) <- levels(df_summary$proc_ES)

facet1 = data.frame(
  proc_ES = factor("FixedES"),
  label1 = "Mean Squared Error",
  label2 = "Variance",
  label3 = "Squared Bias",
  x1 = -0.2, 
  x2 = -0.04, 
  y1 = .13, 
  y2 = .11, 
  y3 = .09
  )

tiff(file="figures/figure8b.tiff",width=2400,height=1300, units = "px", res = 300)
ggplot(data = df_summary, mapping = aes(x = d_actual, y = bias.sq)) +
  facet_wrap(
    ~ proc_ES, 
    labeller = labeller(proc_ES = proc_ES.label), 
    ncol = 4
    ) +
  theme_bw() +
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

#===============================================================================#
################ FIGURE 9: MSE CONDITIONAL ON STOPPING TIME #################### 
#==============================================================================#

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
