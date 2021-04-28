########################## SET UP ########################## 
rm(list = ls())
library(beepr)
library(rpact)
library(pracma)
library(BayesFactor)
library(coda)
library(tidyverse)
library(parallel)

########################## READ ISP FILES ##########################
source("simulations/engine.R") # Obtained from Ulrich & Miller 2020

########################## FIND ALPHA-WEAK FUNCTION ##########################

find_alpha_weak <- function(alpha_total, max_n_segments, alpha_strong) {
  # Obtained from Ulrich & Miller 2020
  # Use numerical search to find the appropriate alpha_weak value
  # that will produce the desired overall_alpha level for the indicated
  # values of n_segments and alpha_strong
  compute_error <- function(try_aw) {
    diff <- try_aw - alpha_strong
    diff_accumulated <- diff^(max_n_segments - 1)
    alpha_strong * (1-diff_accumulated)/(1-diff)+
      try_aw * diff_accumulated-alpha_total 
  }
  aw_range <- c(alpha_strong, alpha_total^(1/max_n_segments))
  final_aw <- fzero(compute_error, aw_range)
  final_aw$x
}

#find_alpha_weak(0.05, 3, .025) #returns .28

########################## SIMULATION FUNCTIONS ##########################

generate_data = function(nsims, n, d) { 
  # Generate raw data
  raw = lapply(1:nsims, function(i) {
    data.frame(nsim = i, n = n, d = d,
               m2 = rnorm(n = n, mean = 0),
               m1 = rnorm(n = n, mean = 0) + d)})
  bind_rows(raw)
}

run_tests = function(df, n, proc, segment = 1, alpha) {
  # Run t-tests on raw data
  p = lapply(1:max(df$nsim), function(i) {
    df %>% filter(nsim == i) %>% #run 1 analysis per simulation
      summarize(p = t.test(m1, m2, alternative = "greater", var.equal = T)$p.value)
    })
  
  # Create final data-frame
  df %>% group_by(nsim) %>% 
    summarize(proc = proc, segment = segment, 
              alpha = alpha, n = mean(n), 
              sd1 = sd(m1), sd2 = sd(m2),
              m1 = mean(m1), m2 = mean(m2)) %>% 
    mutate(sdpooled = sqrt((sd1^2 +sd2^2)/2)) %>% 
    mutate(ES = (m1 - m2) / sdpooled) %>% #ES estimate = Cohen's d
    cbind(tibble(p) %>% 
            unnest(cols = c(p))) #add p-value
} 

########################## SEQUENTIAL PROCEDURE FUNCTION ########################## 

#########  FUNCTION WITH BIAS ADJUSTMENT AND ADJUSTED INFORMATION RATES FOR O'BRIEN-FLEMING #########  

group_sequential = function(proc, target_power, d_forpower, d_actual) { 

bs = ifelse(proc == "asOF", "bsOF", "bsP") #Specify beta-spending function
if(proc == "asOF") {rates = c(0.50, 0.75, 1)} else{rates = c(1/3, 2/3, 1)} #specify information rates for Pocock vs OBF

# Specify the design
design = getDesignGroupSequential(sided = 1,
                                  alpha = alpha_total, 
                                  beta = 1 - target_power,
                                  kMax = max_n_segments, 
                                  typeOfDesign = proc,
                                  typeBetaSpending= bs,
                                  informationRates = rates,
                                  bindingFutility = TRUE) #binding futility bounds

# Get parameters
parameters = getSampleSizeMeans(design = design, 
                                groups = 2, 
                                alternative = d_forpower)
n_gs = ceiling(c(parameters$numberOfSubjects[1], diff(parameters$numberOfSubjects))/2) # n per look per group
alpha = parameters$criticalValuesPValueScale
futility = parameters$futilityBoundsPValueScale

# RUN LOOK 1
set.seed(1234*(d_forpower*target_power))
raw1 = generate_data(nsims = nsims, n = n_gs[1], d = d_actual) %>% 
  mutate(id = nsim, 
         n_cumulative = n_gs[1],
         look = 1) 

gs1 = run_tests(df = raw1, n = n_gs[1], proc = proc, alpha = alpha[1]) %>% 
  mutate(id = nsim, 
         n_cumulative = n_gs[1],
         decision = ifelse(p <= alpha, "Reject",
                           ifelse(p > futility[1], "SupportNull", "Continue")))
#Support null not to be taken literally; stopped for futility based on beta-spending

keep1 = gs1 %>% filter(decision == "Continue") %>% select(nsim) %>% pull()

# RUN LOOK 2
gs1keep = raw1 %>% #keep raw data of unfinished trials after look 1
  filter(nsim %in% keep1) %>% 
  mutate(id = nsim, #create experiment-id variable 
         nsim = rep(1:length(keep1), each = n_gs[1])) 

set.seed(1512*(d_forpower*target_power))
raw2 = generate_data(nsims = sum(gs1$decision=="Continue"), n = n_gs[2], d = d_actual) %>% 
  mutate(id = rep(keep1, each = n_gs[2]), 
         n_cumulative = n_gs[1] + n_gs[2],  #Do I need to adjust n here?
         look = 2) 

gs2 = bind_rows(gs1keep, raw2) %>% 
  run_tests(n = n_gs[2], proc = proc, alpha = alpha[2], segment = 2) %>% 
  mutate(id = keep1,
         n = n_gs[2],
         n_cumulative = n_gs[1] + n_gs[2],
         decision = ifelse(p <= alpha, "Reject", ifelse(p > futility[2], "SupportNull", "Continue")))
#Support null not to be taken literally; stopped for futility based on beta-spending

keep2 = gs2 %>% filter(decision == "Continue") %>% select(nsim) %>% pull()

# RUN LOOK 3
gs1keep = raw1 %>% filter(nsim %in% keep1[keep2]) %>%
  mutate(nsim = rep(1:length(keep2), each = n_gs[1])) #raw data look 1, unfinished trials
gs2keep = raw2 %>% filter(nsim %in% keep2) %>% 
  mutate(nsim = rep(1:length(keep2), each = n_gs[2])) #raw data look 2, unfinished trials

set.seed(0304*(d_forpower*target_power))
raw3 = generate_data(nsims = sum(gs2$decision=="Continue"), n = n_gs[3], d = d_actual) %>% 
  mutate(id = rep(keep1[keep2], each = n_gs[3]), #mistake here 
         n_cumulative = sum(n_gs),
         look = 3)

gs3 = bind_rows(gs1keep, gs2keep, raw3) %>% 
  run_tests(n = n_gs[3], proc = proc, alpha = alpha[3], segment = 3) %>% 
  mutate(id = keep1[keep2],
         n = n_gs[3],
         n_cumulative = sum(n_gs), 
         decision = ifelse(p <= alpha, "Reject", "FTR"))

# Bias-correction with rpact

raw_merged = bind_rows(raw1, raw2, raw3) %>% 
  group_by(id, look) %>% 
  summarize(n = mean(n),
            sd1 = sd(m1),
            sd2 = sd(m2),
            m1 = mean(m1),
            m2 = mean(m2)) %>% 
  mutate(sdpooled = sqrt((sd1^2 +sd2^2)/2),
         ES = (m1 - m2)/ sdpooled) #ES = uncorrected Cohen's d

# Get median unbiased estimate

unbiased = raw_merged %>% 
  split(.$id) %>% 
  map(~ getDataset(
    n1 = .$n,
    n2 = .$n,
    means1 = .$m1,
    means2 = .$m2,
    stDevs1 = .$sd1,
    stDevs2 = .$sd2
  )) %>% 
    map(~ getFinalConfidenceInterval(
    design = design,
    dataInput = .x
  )) %>%
  map("medianUnbiased") %>% 
  data_frame(ES_corrected = .) 

##############################
# Add to all dataframes: n, n_cumulative, ES_corrected

# Create final data-frame

df = bind_rows(gs1, gs2, gs3) %>% 
  filter(decision != "Continue") %>% 
  arrange(id) %>% 
  bind_cols(unbiased) %>% 
  mutate(d_forpower = d_forpower, d_actual = d_actual, power = target_power) %>% 
  select(id, proc, segment, n, n_cumulative, d_forpower, d_actual,  power, decision, ES, ES_corrected)

}

#########  FUNCTION WITH ADJUSTED INFORMATION RATES FOR O'BRIEN-FLEMING #########
group_sequential = function(proc, target_power, d_forpower, d_actual) { #adjusted information rates for O'Brien-Fleming
  bs = ifelse(proc == "asOF", "bsOF", "bsP")
  if(proc == "asOF") {rates = c(0.50, 0.75, 1)} else{rates = c(1/3, 2/3, 1)}
  
  # Specify the design
  design = getDesignGroupSequential(sided = 1,
                                    alpha = alpha_total, 
                                    beta = 1 - target_power,
                                    kMax = max_n_segments, 
                                    typeOfDesign = proc,
                                    typeBetaSpending= bs,
                                    informationRates = rates) 
  
  # Get parameters
  parameters = getSampleSizeMeans(design = design, 
                                  groups = 2, 
                                  alternative = d_forpower)
  n_gs = ceiling(c(parameters$numberOfSubjects[1], diff(parameters$numberOfSubjects))/2) # n per look per group
  alpha = parameters$criticalValuesPValueScale
  futility = parameters$futilityBoundsPValueScale
  
  # RUN LOOK 1
  set.seed(1234*(d_forpower*target_power))
  raw1 = generate_data(nsims = nsims, n = n_gs[1], d = d_actual)
  gs1 = run_tests(df = raw1, n = n_gs[1], proc = proc, alpha = alpha[1]) %>% 
    mutate(decision = ifelse(p <= alpha, "Reject", 
                             ifelse(p > futility[1], "SupportNull", "Continue")))
  #Support null not to be taken literally; stopped for futility based on beta-spending
  
  keep1 = gs1 %>% filter(decision == "Continue") %>% select(nsim) %>% pull()
  
  # RUN LOOK 2
  gs1keep = raw1 %>% 
    filter(nsim %in% keep1) %>% 
    mutate(nsim = rep(1:length(keep1), each = n_gs[1])) #keep raw data of unfinished trials after look 1
  
  set.seed(1512*(d_forpower*target_power))
  raw2 = generate_data(nsims = sum(gs1$decision=="Continue"), n = n_gs[2], d = d_actual)
  
  gs2 = bind_rows(gs1keep, raw2) %>% 
    run_tests(n = n_gs[2], proc = proc, alpha = alpha[2], segment = 2) %>% 
    mutate(n = n_gs[1] + n_gs[2], decision = ifelse(p <= alpha, "Reject", 
                                                 ifelse(p > futility[2], "SupportNull", "Continue")))
  #Support null not to be taken literally; stopped for futility based on beta-spending
  
  keep2 = gs2 %>% filter(decision == "Continue") %>% select(nsim) %>% pull()
  
  # RUN LOOK 3
  gs1keep = raw1 %>% filter(nsim %in% keep1[keep2]) %>%
    mutate(nsim = rep(1:length(keep2), each = n_gs[1])) #raw data look 1, unfinished trials
  gs2keep = raw2 %>% filter(nsim %in% keep2) %>% 
    mutate(nsim = rep(1:length(keep2), each = n_gs[2])) #raw data look 2, unfinished trials
  
  set.seed(0304*(d_forpower*target_power))
  raw3 = generate_data(nsims = sum(gs2$decision=="Continue"), n = n_gs[3], d = d_actual)
  gs3 = bind_rows(gs1keep, gs2keep, raw3) %>% 
    run_tests(n = n_gs[3], proc = proc, alpha = alpha[3], segment = 3) %>% 
    mutate(n = sum(n_gs), decision = ifelse(p <= alpha, "Reject", "FTR"))
  
  # Create final data-frame
  bind_rows(gs1, gs2, gs3) %>% 
    filter(decision != "Continue") %>% 
    mutate(id = row_number(), d_forpower = d_forpower, d_actual = d_actual, power = target_power) %>% 
    select(id, proc, segment, n, d_forpower, d_actual,  power, decision, ES)
}

################# FUNCTION WITH EQUAL INFORMATION RATES #################

group_sequential = function(proc, target_power, d_forpower, d_actual) {
  bs = ifelse(proc == "asOF", "bsOF", "bsP") #equal information rates
  
  # Specify the design
  design = getDesignGroupSequential(sided = 1,
                                    alpha = alpha_total, 
                                    beta = 1 - target_power,
                                    kMax = max_n_segments, 
                                    typeOfDesign = proc,
                                    typeBetaSpending= bs) 
  
  # Get parameters
  parameters = getSampleSizeMeans(design = design, 
                                  groups = 2, 
                                  alternative = d_forpower)
  n_gs = ceiling(parameters$maxNumberOfSubjects / max_n_segments/2) # n per look per group
  alpha = parameters$criticalValuesPValueScale
  futility = parameters$futilityBoundsPValueScale
  
  # RUN LOOK 1
  set.seed(1234*(d_forpower*target_power))
  raw1 = generate_data(nsims = nsims, n = n_gs, d = d_actual)
  gs1 = run_tests(df = raw1, n = n_gs, proc = proc, alpha = alpha[1]) %>% 
    mutate(decision = ifelse(p <= alpha, "Reject", 
                             ifelse(p > futility[1], "SupportNull", "Continue")))
  #Support null not to be taken literally; stopped for futility based on beta-spending
  
  keep1 = gs1 %>% filter(decision == "Continue") %>% select(nsim) %>% pull()
  
  # RUN LOOK 2
  gs1keep = raw1 %>% 
    filter(nsim %in% keep1) %>% 
    mutate(nsim = rep(1:length(keep1), each = n_gs)) #keep raw data of unfinished trials after look 1
  
  set.seed(1512*(d_forpower*target_power))
  raw2 = generate_data(nsims = sum(gs1$decision=="Continue"), n = n_gs, d = d_actual)
  
  gs2 = bind_rows(gs1keep, raw2) %>% 
    run_tests(n = n_gs, proc = proc, alpha = alpha[2], segment = 2) %>% 
    mutate(n = n_gs * segment, decision = ifelse(p <= alpha, "Reject", 
                                                 ifelse(p > futility[2], "SupportNull", "Continue")))
  #Support null not to be taken literally; stopped for futility based on beta-spending
  
  keep2 = gs2 %>% filter(decision == "Continue") %>% select(nsim) %>% pull()
  
  # RUN LOOK 3
  gs1keep = raw1 %>% filter(nsim %in% keep1[keep2]) %>%
    mutate(nsim = rep(1:length(keep2), each = n_gs)) #raw data look 1, unfinished trials
  gs2keep = raw2 %>% filter(nsim %in% keep2) %>% 
    mutate(nsim = rep(1:length(keep2), each = n_gs)) #raw data look 2, unfinished trials
  
  set.seed(0304*(d_forpower*target_power))
  raw3 = generate_data(nsims = sum(gs2$decision=="Continue"), n = n_gs, d = d_actual)
  gs3 = bind_rows(gs1keep, gs2keep, raw3) %>% 
    run_tests(n = n_gs, proc = proc, alpha = alpha[3], segment = 3) %>% 
    mutate(n = n_gs*segment, decision = ifelse(p <= alpha, "Reject", "FTR"))
  
  # Create final data-frame
  bind_rows(gs1, gs2, gs3) %>% 
    filter(decision != "Continue") %>% 
    mutate(id = row_number(), d_forpower = d_forpower, d_actual = d_actual, power = target_power) %>% 
    select(id, proc, segment, n, d_forpower, d_actual,  power, decision, ES)
}

######################### FULL PROCEDURE FUNCTION ########################## 
# All frequentist procedures use naive Cohen's d as the effect size estimate 
# Bayesian Effect Size estimate is based on MCMC samples from the posterior
# For the Bayesian hypothesis test, I obtain a Bayes Factor using the one-sided Bayesian two-sample t-test
# For the Bayesian effect size estimate, I obtain a posterior distribution using the two-sided Bayesian two-sample t-test
# (As recommended in Van Doorn et al. 2020)

procedure = function(d_forpower, d_actual) {
  
  ########################## FIXED SAMPLE ##########################
  
  # Calculate fixed n per group for the one-tailed two-sample t-test
  n_fixed = ceiling(power.t.test(delta = d_forpower, sig.level = alpha_total, power = target_power, 
                                 type = "two.sample", alternative = "one.sided")$n)
  
  # RUN FIXED SAMPLE
  set.seed(1234*(d_forpower*target_power))
  raw = generate_data(nsims = nsims, n = n_fixed, d = d_actual)
  fixed = run_tests(df = raw, n = n_fixed, proc = "Fixed", alpha = alpha_total) %>% 
    mutate(decision = ifelse(p <= alpha, "Reject", "FTR"), id = nsim, 
           d_forpower = d_forpower, d_actual = d_actual, power = target_power) %>% 
    select(id, proc, segment, n, d_forpower, d_actual, power, decision, ES)
  
  ########################## ISP ##########################
  
  # Calculate n per group per segment for the ISP
  n_s = ceiling(n_for_power(target_power, stat_procedure_name = "2t", effect_size = d_forpower, 
                            max_n_segments, alpha_total, alpha_strong))/2
  
  # RUN SEGMENT 1
  set.seed(1234*(d_forpower*target_power))
  raw = generate_data(nsims = nsims, n = n_s, d = d_actual)
  isp1 = run_tests(df = raw, n = n_s, proc = "ISP", alpha = alpha_strong) %>% 
    mutate(decision = ifelse(p <= alpha_strong, "Reject",
                             ifelse(p > alpha_weak, "SupportNull", "Continue")))
  #Support null not to be taken literally 
  
  # RUN SEGMENT 2
  set.seed(1512*(d_forpower*target_power))
  raw = generate_data(nsims = sum(isp1$decision=="Continue"), n = n_s, d = d_actual)
  isp2 = run_tests(df = raw, n = n_s, proc = "ISP", alpha = alpha_strong, segment = 2) %>% 
    mutate(decision = ifelse(p <= alpha_strong, "Reject",
                             ifelse(p > alpha_weak, "SupportNull", "Continue")))
  #Support null not to be taken literally
  
  # RUN SEGMENT 3
  set.seed(0304*(d_forpower*target_power))
  raw = generate_data(nsims = sum(isp2$decision=="Continue"), n = n_s, d = d_actual)
  isp3 = run_tests(df = raw, n = n_s, proc = "ISP", alpha = alpha_weak, segment = 3) %>% 
    mutate(decision = ifelse(p <= alpha_weak, "Reject", "FTR"))
  
  isp = bind_rows(isp1, isp2, isp3) %>% filter(decision != "Continue") %>% 
    mutate(id = row_number(), n = n_s*segment, d_forpower = d_forpower, d_actual = d_actual, power = target_power) %>% 
    select(id, proc, segment, n, d_forpower, d_actual, power, decision, ES)
  
  ########################## GROUP SEQUENTIAL ##########################
  
  pocock = group_sequential(proc = "asP", target_power = target_power, d_forpower = d_forpower, d_actual = d_actual)
  obrien = group_sequential(proc = "asOF", target_power = target_power, d_forpower = d_forpower, d_actual = d_actual)
  
  ########################## BAYES ##########################
  
  # Set parameters
  minN = 22 # minimum n before optional stopping is started
  maxN = 22*3 # maximum n - if that is reached without hitting a boundary, the run is aborted
  cumulativeNs = seq(minN, maxN, minN)
  maxBoundary = log(3) # maximal boundary for stopping, in log-units
  
  # Run all data in batches of minN
  set.seed(1234*(d_forpower*target_power))
  l = lapply(1:length(cumulativeNs), function(i){generate_data(nsims = nsims, n = minN, d = d_actual)})
  df = tibble(l) %>% unnest(cols = c(l)) %>% arrange(nsim) %>% 
    mutate(segment = rep(rep(1:length(cumulativeNs), each = minN), nsims))
  
  # Run Bayesian t-test (get log BF for H2: delta > 0; log BF for delta = 0: -logBF)
  b = lapply(1:length(cumulativeNs), function(i){df %>% 
      group_by(nsim) %>% filter(segment %in% 1:i) %>% 
      summarize(BF = ttestBF(m2, m1, nullInterval = c(0, Inf))@bayesFactor$bf[1])}) #one-sided
  
  # Keep all data until we've hit the boundary (if no boundary is hit, keep all data)
  logBF = tibble(b) %>% unnest(cols = c(b)) %>% 
    arrange(nsim) %>% mutate(hit = ifelse(abs(BF) >= maxBoundary, 1, 0),
                             segment = rep(1:length(cumulativeNs), nsims)) %>% 
    group_by(nsim) %>% mutate(hits = cumsum(hit)) %>% 
    filter(cumsum(hits) < 2) %>% slice_tail() %>% ungroup() %>% 
    mutate(decision = ifelse(BF >= maxBoundary, "Reject",
                             ifelse(BF <= -maxBoundary, "SupportNull", "FTR")),
           n = segment * minN) %>% 
    select(-hits, -hit)
  
  # Get Bayesian Effect Size Estimate
  df = df %>% mutate(keep = rep(logBF$segment, each = length(cumulativeNs)*minN)) %>% 
    group_by(nsim) %>% filter(segment <= keep)
  
  set.seed(1512*(d_forpower*target_power))
  models = df %>% split(.$nsim) %>%
    map(~ ttestBF(.$m2, .$m1)) %>% #two-sided
    map(posterior, iterations = 1000, progress = FALSE) %>% 
    map(~ mean(.x[,"delta"])) %>% 
    tibble() %>% unnest(cols = c(.)) %>% rename(ES = ".")
  
  # Combine results
  bayes = logBF %>% cbind(models) %>% 
    mutate(id = row_number(), proc = "Bayes", d_forpower = d_forpower, d_actual = d_actual, power = target_power) %>% 
    select(id, proc, segment, n, d_forpower, d_actual, power, decision, ES)
  
  ########################## MERGE DATA ##########################
  
  bind_rows(fixed, isp, pocock, obrien, bayes)
}

########################## RUN SIMULATION 01 ##########################
# One-tailed, two-sample t-tests
# 5 procedures: Fixed, ISP, GS Pocock, GS O'Brien Fleming, & Sequential Bayes Factor
# Includes expected effect sizes (d = 0.5) and unexpected effect sizes
# Returns simulation01.csv

########################## SET PARAMETERS ##########################
nsims = 50000
max_n_segments = 3
#d_actual = 0.2
#d_forpower = 0.5

d_actual = c(seq(-0.2, 0.4, 0.2), 0.5, seq(0.6, 1, 0.2))
d_forpower = rep(0.5, 8)
target_power = 0.8
alpha_total = 0.05
alpha_strong = 0.025
alpha_weak= find_alpha_weak(alpha_total, max_n_segments, alpha_strong) #calculate alpha_weak for ISP

########################## PARALELLELIZE ##########################

system.time(dt <- mclapply(1:length(d_actual), 
                           function(i) procedure(d_forpower[i], d_actual[i]),
                           mc.cores = 2)) #takes ~56mins to run when nsim = 10,000 & ncores = 2
                                          #takes ~7.9 hours to run when nsim = 50,000 & ncores = 2
beep()
df = dt %>% bind_rows()
#write.csv(df, "simulations/simulation001.csv", row.names = F)

#write in zip to compress: 
write.csv(df, file=gzfile("simulations/simulation001.csv.gz"), row.names = F)
zz = gzfile("simulations/simulation001.csv.gz", 'rt')
df = read.csv(zz, header = T)

########################## SPLIT UP PROCEDURES IN BAYES & FREQUENTIST  ########################## 

procedure = function(d_forpower, d_actual) {
  
  ########################## FIXED SAMPLE ##########################
  
  # Calculate fixed n per group for the one-tailed two-sample t-test
  n_fixed = ceiling(power.t.test(delta = d_forpower, sig.level = alpha_total, power = target_power, 
                                 type = "two.sample", alternative = "one.sided")$n)
  
  # RUN FIXED SAMPLE
  set.seed(1234*(d_forpower*target_power))
  raw = generate_data(nsims = nsims, n = n_fixed, d = d_actual)
  fixed = run_tests(df = raw, n = n_fixed, proc = "Fixed", alpha = alpha_total) %>% 
    mutate(decision = ifelse(p <= alpha, "Reject", "FTR"), id = nsim, 
           d_forpower = d_forpower, d_actual = d_actual, power = target_power) %>% 
    select(id, proc, segment, n, d_forpower, d_actual, power, decision, ES)
  
  ########################## ISP ##########################
  
  # Calculate n per group per segment for the ISP
  n_s = ceiling(n_for_power(target_power, stat_procedure_name = "2t", effect_size = d_forpower, 
                            max_n_segments, alpha_total, alpha_strong))/2
  
  # RUN SEGMENT 1
  set.seed(1234*(d_forpower*target_power))
  raw = generate_data(nsims = nsims, n = n_s, d = d_actual)
  isp1 = run_tests(df = raw, n = n_s, proc = "ISP", alpha = alpha_strong) %>% 
    mutate(decision = ifelse(p <= alpha_strong, "Reject",
                             ifelse(p > alpha_weak, "SupportNull", "Continue")))
  #Support null not to be taken literally 
  
  # RUN SEGMENT 2
  set.seed(1512*(d_forpower*target_power))
  raw = generate_data(nsims = sum(isp1$decision=="Continue"), n = n_s, d = d_actual)
  isp2 = run_tests(df = raw, n = n_s, proc = "ISP", alpha = alpha_strong, segment = 2) %>% 
    mutate(decision = ifelse(p <= alpha_strong, "Reject",
                             ifelse(p > alpha_weak, "SupportNull", "Continue")))
  #Support null not to be taken literally
  
  # RUN SEGMENT 3
  set.seed(0304*(d_forpower*target_power))
  raw = generate_data(nsims = sum(isp2$decision=="Continue"), n = n_s, d = d_actual)
  isp3 = run_tests(df = raw, n = n_s, proc = "ISP", alpha = alpha_weak, segment = 3) %>% 
    mutate(decision = ifelse(p <= alpha_weak, "Reject", "FTR"))
  
  isp = bind_rows(isp1, isp2, isp3) %>% filter(decision != "Continue") %>% 
    mutate(id = row_number(), n = n_s*segment, d_forpower = d_forpower, d_actual = d_actual, power = target_power) %>% 
    select(id, proc, segment, n, d_forpower, d_actual, power, decision, ES)
  
  ########################## GROUP SEQUENTIAL ##########################
  
  pocock = group_sequential(proc = "asP", target_power = target_power, d_forpower = d_forpower, d_actual = d_actual)
  obrien = group_sequential(proc = "asOF", target_power = target_power, d_forpower = d_forpower, d_actual = d_actual)
  
  ########################## MERGE DATA ##########################
  
  bind_rows(fixed, isp, pocock, obrien)
}

########################## RUN SIMULATION 01 ##########################
# One-tailed, two-sample t-tests
# 4 procedures: Fixed, ISP, GS Pocock, GS O'Brien Fleming
# Includes expected effect sizes (d = 0.5) and unexpected effect sizes
# Returns simulation001-freq.csv

########################## SET PARAMETERS ##########################
nsims = 10000
max_n_segments = 3
#d_actual = 0.2
#d_forpower = 0.5

d_actual = c(seq(-0.2, 0.4, 0.2), 0.5, seq(0.6, 1, 0.2))
d_forpower = rep(0.5, 8)
target_power = 0.8
alpha_total = 0.05
alpha_strong = 0.025
alpha_weak= find_alpha_weak(alpha_total, max_n_segments, alpha_strong) #calculate alpha_weak for ISP

########################## PARALELLELIZE ##########################

system.time(dt <- mclapply(1:length(d_actual), 
                           function(i) procedure(d_forpower[i], d_actual[i]),
                           mc.cores = 2)) #takes ~23 mins to run when nsim = 10,000 & ncores = 2 & bayes excluded 
beep()
#df_freq = dt %>% bind_rows()
#write.csv(df_freq, "simulations/simulation001-freq.csv", row.names = F)

#=========================================================#
########################## BAYES ##########################
#=========================================================#

bayes = function(d_forpower, d_actual) {
  
  ########################## BAYES ##########################
  
  # Set parameters
  minN = 23 # minimum n before optional stopping is started
  maxN = 23*3 # maximum n - if that is reached without hitting a boundary, the run is aborted
  cumulativeNs = seq(minN, maxN, minN)
  maxBoundary = log(3) # maximal boundary for stopping, in log-units
  
  # Run all data in batches of minN
  set.seed(1234*(d_forpower*target_power))
  l = lapply(1:length(cumulativeNs), function(i){generate_data(nsims = nsims, n = minN, d = d_actual)})
  df = tibble(l) %>% unnest(cols = c(l)) %>% arrange(nsim) %>% 
    mutate(segment = rep(rep(1:length(cumulativeNs), each = minN), nsims))
  
  # Run Bayesian t-test (get log BF for H2: delta > 0; log BF for delta = 0: -logBF)
  b = lapply(1:length(cumulativeNs), function(i){df %>% 
      group_by(nsim) %>% filter(segment %in% 1:i) %>% 
      summarize(BF = ttestBF(m2, m1, nullInterval = c(0, Inf))@bayesFactor$bf[1])}) #one-sided
  
  # Keep all data until we've hit the boundary (if no boundary is hit, keep all data)
  logBF = tibble(b) %>% unnest(cols = c(b)) %>% 
    arrange(nsim) %>% mutate(hit = ifelse(abs(BF) >= maxBoundary, 1, 0),
                             segment = rep(1:length(cumulativeNs), nsims)) %>% 
    group_by(nsim) %>% mutate(hits = cumsum(hit)) %>% 
    filter(cumsum(hits) < 2) %>% slice_tail() %>% ungroup() %>% 
    mutate(decision = ifelse(BF >= maxBoundary, "Reject",
                             ifelse(BF <= -maxBoundary, "SupportNull", "FTR")),
           n = segment * minN) %>% 
    select(-hits, -hit)
  
  # Get Bayesian Effect Size Estimate
  df = df %>% mutate(keep = rep(logBF$segment, each = length(cumulativeNs)*minN)) %>% 
    group_by(nsim) %>% filter(segment <= keep)
  
  set.seed(1512*(d_forpower*target_power))
  models = df %>% split(.$nsim) %>%
    map(~ ttestBF(.$m2, .$m1)) %>% #two-sided
    map(posterior, iterations = 1000, progress = FALSE) %>% 
    map(~ mean(.x[,"delta"])) %>% 
    tibble() %>% unnest(cols = c(.)) %>% rename(ES = ".")
  
  # Combine results
  logBF %>% cbind(models) %>% 
    mutate(id = row_number(), proc = "Bayes", d_forpower = d_forpower, d_actual = d_actual, power = target_power) %>% 
    select(id, proc, segment, n, d_forpower, d_actual, power, decision, ES)
  
}

########################## SET PARAMETERS ##########################
nsims = 10000
d_actual = c(seq(-0.2, 0.4, 0.2), 0.5, seq(0.6, 1, 0.2))
d_forpower = rep(0.5, 8)
target_power = 0.8
########################## PARALELLELIZE ##########################
system.time(dt <- mclapply(1:length(d_actual), 
                           function(i) bayes(d_forpower[i], d_actual[i]),
                           mc.cores = 2)) #takes ~31mins to run when nsim = 10,000 & ncores = 2
beep()

df_bayes = dt %>% bind_rows()
write.csv(df_bayes, "simulation001-bayes.csv", row.names = F)


#df = bind_rows(df_freq, df_bayes)
#write.csv(df, "simulation001.csv", row.names = F)