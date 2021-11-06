########################## SET UP ########################## 
library(data.table)
library(rpact)
library(pracma)
library(BayesFactor)
library(coda)
library(tidyverse)
library(parallel)
library(TOSTER)
library(BFDA)
library(microbenchmark)
library(NCmisc)
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
target_power = 0.8 #might also add 0.9
d_forpower = 0.5

########################## FULL PROCEDURE FUNCTION ##########################

run_all_procedures = function(nsims){
  
  ########################## FIXED SAMPLING PROCEDURE FUNCTION ########################## 
  
  run_fixed = function(d_actual) {
    
    ########################## FIXED SAMPLE ##########################
    
    # Calculate fixed n per group for the one-tailed two-sample t-test
    n_fixed = ceiling(
      power.t.test(
        delta = d_forpower, 
        sig.level = alpha_total, 
        power = target_power, 
        type = "two.sample", 
        alternative = "one.sided")$n)
    
    # RUN FIXED SAMPLE NHST
    set.seed(1234*(d_forpower*target_power))
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
  d_actual = c(seq(-0.2, 0.4, 0.2), 0.5, seq(0.6, 1, 0.2))
  
  ########################## PARALELLELIZE ##########################
  
  dt <- mclapply(1:length(d_actual), 
                 function(i) run_fixed(d_actual[i]),
                 mc.cores = 2)
  print("The fixed sampling procedure has finished running!")
  
  fixed = data.table::rbindlist(dt)
  
  ########################## GROUP-SEQUENTIAL PROCEDURE FUNCTION ########################## 
  
  #  FUNCTION WITH BIAS ADJUSTMENT AND ADJUSTED INFORMATION RATES FOR O'BRIEN-FLEMING #  
  run_group_sequential = function(proc, d_actual) {
    
    if(proc == "asOF") {design = design_obf} else{design = design_pocock}
    if(proc == "asOF") {parameters = parameters_obf} else{parameters = parameters_pocock}
    
    ns = ceiling(c(parameters$numberOfSubjects[1], 
                   diff(parameters$numberOfSubjects))/2) # n per look per group
    alpha = parameters$criticalValuesPValueScale
    futility = parameters$futilityBoundsPValueScale
    
    # Run all data in batches
    set.seed(1234*(d_forpower*target_power))
    l = lapply(1:length(ns), 
               function(i){generate_data(nsims = nsims, 
                                         n = ns[i], 
                                         d = d_actual)})
    
    df = data.table::rbindlist(l)[order(nsim)] %>% 
      mutate(
        segment = rep(
          lapply(1:length(ns), 
                 function(i) {rep(i, times = ns[i])}) %>%  
            unlist(), 
        nsims)
        )
    
    # Run tests
    tests = lapply(1:length(ns), function(i){
      df %>% 
        group_by(nsim) %>% 
        filter(segment %in% 1:i) %>% #get data of all segments up to i 
        summarize(
          p = t.test(m1, m2, alternative = "greater", var.equal = T)$p.value,
          sd1 = sd(m1),
          sd2 = sd(m2),
          m1 = mean(m1),
          m2 = mean(m2),
          segment = max(segment),
          n = case_when(
            segment == 1 ~ ns[1],
            segment == 2 ~ ns[2],
            segment == 3 ~ ns[3]
            ),
          n_cumulative = case_when(
            segment == 1 ~ ns[1],
            segment == 2 ~ ns[1] + ns[2],
            segment == 3 ~ sum(ns)
            ),
          alpha = case_when(
            segment == 1 ~ alpha[1],
            segment == 2 ~ alpha[2],
            segment == 3 ~ alpha[3]
            ),
          futility = case_when(
            segment == 1 ~ futility[1],
            segment == 2 ~ futility[2],
            segment == 3 ~ 1
            ),
          decision = case_when(
            p <= alpha ~ "reject.null",
            p > futility ~ "reject.alt",
            p > alpha & p < futility ~ "inconclusive"
            ) 
        )
    }) %>% 
      data.table::rbindlist() %>% 
      arrange(nsim) %>%
      group_by(nsim) %>%
      mutate(final_segment = ifelse(decision != "inconclusive" | segment == 3, 1, 0)) %>%
      filter(final_segment == 1) %>%
      slice_head()
    
    #Get bias corrected ES estimate - first split data into individual segments
    raw_data = lapply(1:length(ns), function(i){
      df %>% 
        group_by(nsim) %>% 
        filter(segment == i) %>% #get data by individual segments
        summarize(
          sd1 = sd(m1),
          sd2 = sd(m2),
          m1 = mean(m1),
          m2 = mean(m2),
          segment = i,
          n = case_when(
            segment == 1 ~ ns[1],
            segment == 2 ~ ns[2],
            segment == 3 ~ ns[3]
            ))
    }) %>% 
      data.table::rbindlist() %>% 
      arrange(nsim) %>% 
      mutate(
        final_segment = tests %>% 
          pull(segment) %>% 
          rep(each = 3)
        ) %>% 
      filter(segment <= final_segment)
    
    # Get median unbiased estimate
    ES_corrected = raw_data %>% 
      split(.$nsim) %>% 
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
      unlist()
    
    # Create final data-frame
    tests = tests %>% 
      ungroup() %>% 
      mutate(
        id = nsim, 
        proc = proc,
        d_forpower = d_forpower, 
        d_actual = d_actual, 
        power = target_power, 
        sdpooled = sqrt((sd1^2 +sd2^2)/2),
        ES = (m1 - m2) / sdpooled,
        ES_corrected = ES_corrected
        )
    
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
  
  ########################## RUN GROUP SEQUENTIAL ##########################
  
  ########################## SET PARAMETERS ##########################
  
  get_design = function(proc) {
    bs = ifelse(proc == "asOF", "bsOF", "bsP") #Specify beta-spending function
    if(proc == "asOF") {rates = c(0.50, 0.75, 1)} else{rates = c(1/3, 2/3, 1)} #specify information rates for Pocock vs OBF
    
    # Specify the design
    getDesignGroupSequential(
      sided = 1,
      alpha = alpha_total, 
      beta = 1 - target_power,
      kMax = max_n_segments, 
      typeOfDesign = proc,
      typeBetaSpending= bs,
      informationRates = rates,
      bindingFutility = TRUE #binding futility bounds
      ) 
  }
  
  get_parameters = function(design) {
    # Get parameters
    parameters = getSampleSizeMeans(
      design = design, 
      groups = 2, 
      alternative = d_forpower)
  }
  
  # Get pocock design and parameters
  design_pocock = get_design(proc = "asP")
  parameters_pocock = get_parameters(design = design_pocock)
  
  # Get O'Brien-Fleming design and parameters
  design_obf = get_design(proc = "asOF")
  parameters_obf = get_parameters(design = design_obf)
  
  d_actual = rep(c(seq(-0.2, 0.4, 0.2), 0.5, seq(0.6, 1, 0.2)), 2)
  proc = rep(c("asP", "asOF"), each = 8) 
  
  ########################## PARALELLELIZE ##########################
  dt <- mclapply(1:length(d_actual), 
                 function(i) run_group_sequential(proc[i], d_actual[i]),
                 mc.cores = 2)
  print("The group-sequential procedures have finished running!")
  
  group_sequential = data.table::rbindlist(dt)
  
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
  
  run_isp = function(d_actual) {
    
    ########################## ISP ##########################
    
    # Calculate n per group per segment for the ISP
    ns = ceiling(
      n_for_power(
        target_power, 
        stat_procedure_name = "2t", 
        effect_size = d_forpower, 
        max_n_segments, 
        alpha_total, 
        alpha_strong))/2
    
    # Run all data in batches of minN
    set.seed(1234*(d_forpower*target_power))
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
  d_actual = c(seq(-0.2, 0.4, 0.2), 0.5, seq(0.6, 1, 0.2))
  alpha_strong = 0.025
  alpha_weak= find_alpha_weak(alpha_total, max_n_segments, alpha_strong) #calculate alpha_weak for ISP
  
  ########################## PARALELLELIZE ##########################
  dt <- mclapply(1:length(d_actual), 
                 function(i) run_isp(d_actual[i]),
                 mc.cores = 2)
  print("The Independent Segments Procedure has finished running!")
  
  isp = data.table::rbindlist(dt)
  
  ############################################################################
  ######################### SEQUENTIAL BAYES FACTOR ########################## 
  ############################################################################
  
  run_sbf = function() { 
    
    df = lapply(1:length(d_actual), function(i){
      BFDA.sim(
        expected.ES = d_actual[i], 
        type="t.between",
        prior = list("normal", list(prior.mean = 0.5, prior.variance = 0.3)),
        n.min = 25, 
        n.max = 75, 
        alternative = "greater", #directionality of the prior
        boundary = 3, 
        B = nsims,
        verbose = TRUE, 
        cores = 4, 
        stepsize = 25, 
        design = "sequential") %>% 
        BFDA.analyze(boundary = 3) 
    }) 
    
    df = df %>% 
      map("endpoint") %>% 
      data.table::rbindlist() %>% 
      group_by(true.ES) %>%  
      mutate(
        id = row_number(),
        proc = "Bayes",
        n_cumulative = n,
        n = 25,
        decision = case_when(
          logBF >= log(3) ~ "reject.null",
          logBF <= log(1/3) ~ "reject.alt",
          is.numeric(logBF) ~ "inconclusive"
          ),
        segment = n_cumulative/n,
        d_forpower = d_forpower,
        power = target_power) %>% 
      rowwise() %>% 
      mutate(
        ES_corrected = integrate( #ES_corrected = mean of posterior distribution 
          EV, #EV function in bayes-helpers.R: written by Angelika Stefan
          lower = -Inf, 
          upper = Inf,
          t = statistic,
          n1 = n_cumulative,
          n2 = n_cumulative)$value 
        ) 
    
    data.table(df)[, .(id, 
                       proc, 
                       segment, 
                       n, 
                       n_cumulative, 
                       d_forpower, 
                       d_actual = true.ES, 
                       power, 
                       decision, 
                       ES = emp.ES, #ES = Cohen's d (2*t/sqrt(degrees of freedom)) 
                       ES_corrected)]
    
  }
  
  ########################## RUN SEQUENTIAL BAYES ##########################
  
  ########################## SET PARAMETERS  ##########################
  source("simulations/bayes-helpers.R") # Helper functions obtained from the BFDA package and with very generous help from Angelika Stefan
  d_actual = c(seq(-0.2, 0.4, 0.2), 0.5, seq(0.6, 1, 0.2))
  
  dt <- mclapply(1:length(d_forpower), 
                 function(i) run_sbf(),
                 mc.cores = 2)
  print("The Sequential Bayes Factor has finished running!")
  
  bayes = data.table::rbindlist(dt)
  
  # MSPRT with adjusted information rates -----------------------------------
  
  ########################### MSPRT ###########################
  
  # MSPRT = function(d_actual) {
  #   
  #   # Run all data in batches of minN
  #   set.seed(1234*(d_forpower*target_power))
  #   l = lapply(1:length(ns), function(i){generate_data(nsims = nsims, n = ns[i], d = d_actual)})
  #   
  #   df = data.table::rbindlist(l)[order(nsim)] %>% 
  #     mutate(segment = rep(lapply(1:length(ns), function(i) {rep(i, times = ns[i])}) %>%  unlist(), nsims),
  #            n_cumulative = case_when(
  #              segment == 1 ~ ns[1],
  #              segment == 2 ~ ns[1] + ns[2],
  #              segment == 3 ~ sum(ns)
  #            ))
  #   
  #   # Run LR tests
  #   lr = df %>% 
  #     split(.$nsim) %>%
  #     map(~ implement.MSPRT(
  #       obs1 = .$m1,
  #       obs2 = .$m2,
  #       design.MSPRT.object = design_SPRT,
  #       plot.it = 0, 
  #       verbose = F
  #     ))
  #   
  #   #Keep all data until we've hit the boundary (if no boundary is hit, keep all data)
  #   df %>% 
  #     mutate(n_final = lr %>% 
  #              map("n1") %>%
  #              unlist() %>% 
  #              rep(each = maxN)) %>% 
  #     filter(n_cumulative <= n_final) %>% 
  #     group_by(nsim) %>% 
  #     summarize(sd1 = sd(m1),
  #               sd2 = sd(m2),
  #               m1 = mean(m1),
  #               m2 = mean(m2),
  #               segment = max(segment),
  #               n = case_when(
  #                 segment == 1 ~ ns[1],
  #                 segment == 2 ~ ns[2],
  #                 segment == 3 ~ ns[3]),
  #               n_cumulative = max(n_cumulative)) %>% 
  #     ungroup() %>% 
  #     mutate(id = row_number(),
  #            proc = "SPRT",
  #            d_forpower = d_forpower, 
  #            d_actual = d_actual, 
  #            power = target_power,
  #            decision = lr %>% 
  #              map("decision") %>% 
  #              unlist(),
  #            sdpooled = sqrt((sd1^2 +sd2^2)/2),
  #            ES = (m1 - m2)/ sdpooled, #uncorrected effect size estimate
  #            ES_corrected = ES) %>% 
  #     select(id, proc, segment, n, n_cumulative, d_forpower, d_actual, power, decision, ES, ES_corrected) 
  #   
  # }
  # 
  # ########################## SET PARAMETERS ##########################
  # d_actual = c(seq(-0.2, 0.4, 0.2), 0.5, seq(0.6, 1, 0.2))
  # ns = c(28, 14, 14) #matched with OBF sample sizes
  # maxN = sum(ns)
  # 
  # design_SPRT = design.MSPRT(test.type = "twoT",
  #                            theta1 = d_forpower,
  #                            Type1.target = alpha_total,
  #                            Type2.target = 1 - target_power,
  #                            N1.max = maxN,
  #                            N2.max = maxN,
  #                            batch1.size = ns,
  #                            batch2.size = ns,
  #                            verbose = F) 
  # 
  # ########################## PARALELLELIZE ##########################
  # system.time(dt <- mclapply(1:length(d_actual), 
  #                            function(i) MSPRT(d_actual[i]),
  #                            mc.cores = 2)) 
  # 
  # print("The Modified Sequential Probability Ratio Test has finished running!")
  # 
  # sprt = data.table::rbindlist(dt)
  
  ########################## CREATE FINAL DATA-FRAME ##########################
  #df = bind_rows(fixed, group_sequential, isp, bayes, sprt)
  
  df = bind_rows(fixed, group_sequential, isp, bayes)
  return(df)
  
}

########################## RUN ALL PROCEDURES ##########################
system.time(dt <- mclapply(1, function(i) run_all_procedures(nsims = 10000),
                           mc.cores = 2)) #10,000 sims take about 2.5-3 hours to run (2h40min)

system("say All your procedures have finished running!")
df = data.table::rbindlist(dt)

########################## CREATE FINAL DATA-FRAME ##########################
write.csv(df, "simulations/data.csv", row.names = F)

########################## READ IN DATA ##########################
df = read.csv("simulations/data.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~