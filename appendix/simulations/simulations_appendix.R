########################## SET UP ########################## 
library(data.table)
library(rpact)
library(pracma)
library(tidyverse)
library(parallel)
library(TOSTER)
library(here)

########################## GLOBAL HELPER FUNCTIONS ##########################
generate_data = function(nsims, n, d) { 
  # Generate raw data
  raw = lapply(1:nsims, function(i) {
    data.table(nsim = i, n = n, d = d,
               m2 = rnorm(n = n, mean = 0),
               m1 = rnorm(n = n, mean = 0) + d)})
  data.table::rbindlist(raw)
}

run_tests = function(df, n, proc, segment = 1, alpha) {
  
  # Run t-tests and create data-frame
  df %>% 
    group_by(nsim) %>% 
    summarize(
      p = t.test(m1, m2, alternative = "greater", var.equal = T)$p.value,
      proc = proc, 
      segment = segment, 
      alpha = alpha, 
      n = mean(n), 
      sd1 = sd(m1), 
      sd2 = sd(m2),
      m1 = mean(m1), 
      m2 = mean(m2)
    ) %>% 
    mutate(sdpooled = sqrt((sd1^2 +sd2^2)/2)) %>% 
    mutate(ES = (m1 - m2) / sdpooled) #ES estimate = Cohen's d
  
} 

########################## GLOBAL PARAMETERS ##########################
max_n_segments = 3
alpha_total = 0.05

# Simulations for appendix data ---------------------------------------

########################## FULL PROCEDURE FUNCTION ##########################

run_all_procedures = function(nsims){
  
  ########################## FIXED SAMPLING PROCEDURE FUNCTION ########################## 
  
  run_fixed = function(d_actual, d_forpower, target_power) {
    
    ########################## FIXED SAMPLE ##########################
    
    # Calculate fixed n per group for the one-tailed two-sample t-test
    n_fixed = ceiling(
      power.t.test(
        delta = d_forpower, 
        #delta = d_actual, 
        sig.level = alpha_total, 
        power = target_power, 
        type = "two.sample", 
        alternative = "one.sided")$n)
    
    # RUN FIXED SAMPLE NHST
    set.seed(1234*(d_actual*target_power))
    fixed = generate_data(
      nsims = nsims,
      n = n_fixed, 
      d = d_actual
    ) %>% 
      run_tests(
        n = n_fixed, 
        proc = "Fixed",
        alpha = alpha_total
      ) %>% 
      mutate(
        id = nsim,
        n_cumulative = n, 
        ES_corrected = ES,
        d_forpower = d_forpower,
        #d_forpower = d_actual,
        d_actual = d_actual, 
        power = target_power
      )
    
    # RUN EQUIVALENCE TEST
    eq_bounds = powerTOSTtwo(
      alpha = alpha_total, #set bounds
      statistical_power = target_power,
      N = n_fixed)
    
    eq_tests = fixed %>% 
      split(.$nsim) %>%
      map(~ TOSTtwo(
        m1 = .$m1,
        m2 = .$m2,
        sd1 = .$sd1,
        sd2 = .$sd2,
        n1 = .$n,
        n2 = .$n,
        low_eqbound_d = eq_bounds[1],
        high_eqbound_d = eq_bounds[2],
        alpha = alpha_total,
        var.equal = T,
        verbose = F,
        plot = F)
      ) %>% 
      data.table::rbindlist() %>% 
      rowwise() %>% 
      mutate(outcome = max(TOST_p1, TOST_p2) <= alpha_total)
    
    # Create final data-frame
    fixed = fixed %>% 
      mutate(equivalence = eq_tests$outcome) %>% 
      mutate(decision = ifelse(
        p <= alpha, "reject.null", 
        ifelse(equivalence == TRUE, "reject.alt",
               "inconclusive"))) 
    
    data.table(fixed)[, .(id, 
                          proc, 
                          segment, 
                          n, 
                          n_cumulative, 
                          d_forpower, 
                          d_actual, 
                          power, 
                          decision, 
                          ES, 
                          ES_corrected)]
    
  }
  
  ########################## RUN FIXED SAMPLING PROCEDURE ##########################
  
  ########################## SET PARAMETERS ##########################
  target_power = rep(c(0.7, 0.8, 0.9), each = 84)
  d_actual = rep(rep(seq(0, 1, 0.05), 4), 3)
  d_forpower = rep(rep(seq(0.2, 0.8, 0.2), 21), 3) 
  
  ########################## PARALELLELIZE ##########################
  dt <- mclapply(1:length(d_actual), 
                 function(i) run_fixed(d_actual[i], 
                                       d_forpower[i],
                                       target_power[i]),
                 mc.cores = 2)
  print("The fixed sampling procedure has finished running!")
  
  fixed = data.table::rbindlist(dt)
  
  ############## INDEPENDENT SEGMENTS PROCEDURE FUNCTION ################ 
  
  ########################## READ ISP FUNCTIONS ##########################
  source("simulations/isp-helpers.R") # Obtained from Ulrich & Miller 2020
  
  # Find alpha weak
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
  
  run_isp = function(d_actual, d_forpower, target_power) {
    
    ########################## ISP ##########################
    
    # Calculate n per group per segment for the ISP
    ns = ceiling(
      n_for_power(
        target_power, 
        stat_procedure_name = "2t", 
        effect_size = d_forpower, 
        #effect_size = d_actual,
        max_n_segments, 
        alpha_total, 
        alpha_strong))/2
    
    # Run all data in batches of minN
    set.seed(1234*(d_actual*target_power))
    l = lapply(1:max_n_segments, 
               function(i){generate_data(nsims = nsims, 
                                         n = ns, 
                                         d = d_actual)})
    
    df = data.table::rbindlist(l)[order(nsim)] %>% 
      mutate(segment = rep(rep(1:max_n_segments, each = ns), nsims))
    
    #Run hypothesis test
    tests = lapply(1:max_n_segments, function(i){
      df %>% 
        group_by(nsim) %>% 
        filter(segment == i) %>% #get data by individual segments
        summarize(
          p = t.test(m1, m2, alternative = "greater", var.equal = T)$p.value,
          sd1 = sd(m1),
          sd2 = sd(m2),
          m1 = mean(m1),
          m2 = mean(m2),
          segment = i,
          decision = case_when(
            segment < 3 & p <= alpha_strong ~ "reject.null",
            segment < 3 & p > alpha_weak ~ "reject.alt",
            segment == 3 & p <= alpha_weak ~ "reject.null",
            is.integer(segment) ~ "inconclusive" #code all other cases as inconclusive
          )
        )
    }) %>% 
      data.table::rbindlist() %>% 
      arrange(nsim) %>% 
      group_by(nsim) %>% 
      mutate(final_segment = ifelse(decision != "inconclusive" | segment == 3, 1, 0)) %>% 
      filter(final_segment == 1) %>% 
      slice_head() %>% 
      mutate(
        id = nsim,
        proc = "ISP",
        n = ns,
        n_cumulative = n * segment,
        d_forpower = d_forpower, 
        #d_forpower = d_actual,
        d_actual = d_actual, 
        power = target_power,
        sdpooled = sqrt((sd1^2 +sd2^2)/2),
        ES = (m1 - m2) / sdpooled, #ES estimate = Cohen's d 
        ES_corrected = ES) %>% 
      ungroup() 
    
    # Get final data frame
    data.table(tests)[, .(id, 
                          proc, 
                          segment, 
                          n, 
                          n_cumulative, 
                          d_forpower, 
                          d_actual, 
                          power, 
                          decision, 
                          ES, 
                          ES_corrected)]
    
  }
  
  ########################## RUN INDEPENDENT SEGMENTS ##########################
  
  ########################## SET PARAMETERS ##########################
  alpha_strong = 0.025
  alpha_weak= find_alpha_weak(alpha_total, max_n_segments, alpha_strong) #calculate alpha_weak for ISP
  #d_actual = rep(seq(0, 1, 0.05), 4) 
  #d_forpower = rep(seq(0.2, 0.8, 0.2), 21) #powered to .8 to detect this effect size
  
  ########################## PARALELLELIZE ##########################
  dt <- mclapply(1:length(d_actual), 
                 function(i) run_isp(d_actual[i], 
                                     d_forpower[i],
                                     target_power[i]),
                 mc.cores = 2)
  print("The Independent Segments Procedure has finished running!")
  
  isp = data.table::rbindlist(dt)
  
  ########################## CREATE FINAL DATA-FRAME ##########################
  df = bind_rows(fixed, isp)
  return(df)
  
}

########################## RUN ALL PROCEDURES ##########################

system.time(dt <- mclapply(1, function(i) run_all_procedures(nsims = 10000),
                           mc.cores = 2)) #10,000 sims take 45 minutes to run

system("say All your procedures have finished running!")
df = data.table::rbindlist(dt)

#################### CREATE FINAL DATA-FRAME #####################
write.csv(df, file=gzfile("appendix/simulations/data_appendix.csv.gz"), row.names = F) #write in zip to compress
