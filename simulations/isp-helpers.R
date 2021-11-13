####ISP Engine provided by Ulrich & Miller at https://github.com/milleratotago/Independent_Segments_R

# Class SegmentedHypTestEngine

# Runs (computes the probabilistic results of) complete segmented hypothesis testing scenarios

##########################################
# Creation methods - following Wickham
##########################################
#=========================================
# new_SegmentedHypTestEngine
#=========================================

new_SegmentedHypTestEngine <- function()
{
  
  # No fields needed at creation
  new_SegmentedHypTestEngine <- list()
  
  attr(new_SegmentedHypTestEngine, "class") <- "SegmentedHypTestEngine"
  
  # Return instance
  return(new_SegmentedHypTestEngine)
  
} # new_SegmentedHypTestEngine
#=========================================

#=========================================
# validate_SegmentedHypTestEngine
#=========================================
# Checks additional properties of elements, for example, that values are in allowed ranges -- business rules
validate_SegmentedHypTestEngine <- function(segmented_hyp_test_engine)
{
  # More elaborate checking in here if needed. Stop if things go wrong.
  
} # end validate_SegmentedHypTestEngine
#=========================================

#=========================================
# SegmentedHypTestEngine - Public ctor
#=========================================
#-----------------
# Roxygen comments

#' SegmentedHypTestEngine
#'
#' \code{SegmentedHypTestEngine}  is the public-facing SegmentedHypTestEngine
#' constructor. It accepts no parameters. This object exposes the primary methods
#' for numerically determing the expected outcomes of a segmented hypothesis testing
#' study. See run_scenario.SegmentedHypTestEngine for an example.
#' @return Validated SegmentedHypTestEngine instance
#'
#' @export
SegmentedHypTestEngine <- function()
{
  # Create instance
  instance = new_SegmentedHypTestEngine()
  
  # Validate instance. Execution halts here if things go wrong
  validate_SegmentedHypTestEngine(instance)
  
  # Return validated instance
  return(instance)
  
} # end SegmentedHypTestEngine ctor
#=========================================


##########################################
# Public methods
##########################################
#' @export
#=========================================
# print.SegmentedHypTestEngine: System will dispatch to this method on print(instance)
#=========================================
print.SegmentedHypTestEngine <- function(x, ...)
{
  data_values <- unclass(x) # This grabs the original list from the struct built when attr was called
  df <- data.frame(data_values)
  print(df)
}

#=========================================
# run_scenario: Top-level wrapper for running a complete scenario.
#
# Accepts researcher, statistical procedure and list (must be a list) of TrueEffect instances. Each TrueEffect
# hold size and baserate probability. A TrueEffect with size 0 (i.e. H0 is true) may be included.
# Probabilities across all effect sizes must sum to 1. A single SegmentedResult object is returned, which is the
# average across all effect sizes, weighted by associated probability.
#=========================================
#-----------------
# Roxygen comments

#' run_scenario.SegmentedHypTestEngine
#'
#' \code{run_scenario.SegmentedHypTestEngine} numerically determines expected
#' outcomes for a complete study scenario. It wraps various internal methods of
#' class SegmentedHypTestEngine (developers can cf. source code).
#' @param segmented_hyp_test_engine A SegmentedHypTestEngine instance for method
#'   dispatch
#' @param segmented_researcher Correctly initialised instance of
#'   SegmentedResearcher
#' @param stat_procedure Correctly initialised instance of any class descended
#'   from StatProcedureBase, e.g. OneSampleT, PearsonR, etc.
#' @param effects_list A list of TrueEffect instances. This must be a list, not
#'   simply a vector. Cast with list() if necessary. Each TrueEffect instance
#'   holds fields effect_size and effect_size_probability. A TrueEffect with
#'   effect_size 0 (i.e. H0 is true) may be included. Probabilities across all effect
#'   sizes must sum to 1 or an exception is thrown.
#' @return A SegmentedHypTestResult instance which holds the average expected outcomes
#' across all effect sizes, weighted by associated probability.
#'
#' @export
run_scenario.SegmentedHypTestEngine <- function(segmented_hyp_test_engine,
                                                segmented_researcher,
                                                stat_procedure,
                                                effects_list)
{
  # Check that effect probs sum to 1. If they do not, throw an exception
  
  # Grab all the probability values from the true_effects list
  probs <- sapply(effects_list, function(te) te[["effect_size_probability"]])
  
  if (sum(probs) != 1)
  {
    stop("The sum of all effect size probabilities must be 1")
  }
  
  # Determine how many effects were passed in
  n_effects <- length(effects_list)
  
  # Make skeleton to hold results
  individual_results <- vector(mode="list", length = n_effects)
  
  # Run the experimental scenario (i.e. compute outcome probabilities) for each effect size and accrue results
  for(i in 1:n_effects)
  {
    individual_results[[i]] <- run_any_subjects_scenario(segmented_researcher, stat_procedure, effects_list[[i]])
  }
  
  # Pass results and probs to weighted averager. Returns a SegmentedHypTestResult instance
  result <- weighted_avg_hyp_test_result(individual_results, probs)
  
  # Return weighted average SegmentedHypTestResult instance
  return(result)
}

#####################################################################################
# average_power(): By computational definition, the average power across a set of
# experiments. Needs to run from the raw data, as SegHTResults do not carry their
# original true effect size(s), so needs access to the scenario methods.
# Called by user-facing methods in SegmentedHypTester.R
#####################################################################################
#' @export
average_power.SegmentedHypTestEngine <- function(segmented_hyp_test_engine,
                                                 max_n_segments,
                                                 n_per_segment,
                                                 alpha_total,
                                                 alpha_strong,
                                                 stat_procedure,
                                                 effects_list)
{
  
  # We exclude effect sizes of zero (i.e. when H0 is true) from this power calculation
  # Drop any element of effects_list with effect_size == 0
  real_effects_list <- effects_list[sapply(effects_list, function(te) te[["effect_size"]] > 0)]
  
  
  # Skeleton to hold all results, one for each provided effect size greater than 0
  n_real_effects <- length(real_effects_list)
  results_list <- list(n_real_effects)
  
  
  # For each provided true effect, run the experimental scenario, producing a SegmentedHypTestResult
  
  # Same researcher information used throughout
  segmented_researcher <- SegmentedResearcher(max_n_segments,n_per_segment, alpha_total, alpha_strong)
  
  # We call run_any_subject_scenario here to run the partial analysis with total effect prob < 1
  # Here we just want this single result, with or without baserate, not the weighted average
  # Run for each effect and accrue the results
  for (i in 1:length(real_effects_list))
  {
    curr_effect <- real_effects_list[[i]]
    result <- run_any_subjects_scenario(segmented_researcher, stat_procedure, curr_effect)
    results_list[[i]] <- result
  }
  
  # Get the average power. Computations by definition (cf. Miller, 2020)
  
  # We multiply each baserate by its associated probability of a true positive, and sum these.
  # Average power is that sum divided by the total baserate where effect size is greater than 0.
  
  pr_true_pos <- sapply(results_list, function(result) result[["pr_true_pos"]])[1:n_real_effects]
  effect_probs <- sapply(real_effects_list, function(te) te[["effect_size_probability"]])
  
  # Perform the computation as described
  power_numerator <- sum(pr_true_pos * effect_probs)
  power_denominator <- sum(effect_probs)
  avg_power <- power_numerator/power_denominator
  
  return(avg_power)
  
} # end average_power

##########################################
# Internal Methods
##########################################
#=========================================
# outcome_probabilities
#=========================================

# This method computes the probability of each outcome (reject, fail to reject, or continue) at each
# of kmax stages for a given experimental scenario. See Miller & Ulrich, 2020, esp. Figures 2A and 2B for further
# computational details.

# Note that when Ho is true, "pr_beat_strong" (causing E to reject and stop) and "pr_beat_weak" (1 - pr_beat_weak causing
# E to fail to reject and stop) are exactly equal to alpha strong and alpha weak. However, when Ho is false these
# values depend on the effect size and statistical test being used.

outcome_probabilities <- function(max_n_segments,pr_beat_strong,pr_beat_weak)
{
  # The method returns a list with the following fields:
  
  #   pr_reject_total: Total probability of rejecting at any step during the experiment
  #   pr_reject_by_segment: Probability of stopping after step i & rejecting.
  #   pr_ftr_by_segmenmt: Probability of stopping after step i & failing to reject.
  #   pr_continue_by_segment: Probability of continuing (neither reject or ftr) at step i
  #   expected_n_segments: The mean of the number of steps taken.
  #   variance_n_segments: The variance of the number of steps taken.
  
  
  # Prepare skeleton lists to hold outputs
  pr_reject_by_segment <- numeric(max_n_segments)
  pr_ftr_by_segment <- numeric(max_n_segments)
  pr_continue_by_segment <- numeric(max_n_segments)
  
  # continue if p is neither greater than weak or less than strong
  pr_continue <- pr_beat_weak - pr_beat_strong
  
  # Compute the probabilistic outcomes for each of max_n_segments steps
  for (curr_segment in 1:max_n_segments)
  {
    # To reject at this segment, you have continued curr_segment-1 previous times,
    # and now you beat (i.e. are less than or equal to) alpha strong.
    # Take the product of the probabilities of these events
    pr_reject_by_segment[curr_segment] <- pr_beat_strong * pr_continue^(curr_segment-1)
    
    # To fail to reject at this segment, you have continued curr_segment-1 previous times,
    # and now you DO NOT BEAT (i.e. are not less than or equal to; are greater than) pr_beat_weak
    pr_ftr_by_segment[curr_segment] <- (1-pr_beat_weak) * pr_continue^(curr_segment-1)
    
    # To continue at this segment, you have continued at all previous segments, and this one as well
    pr_continue_by_segment[curr_segment] <- pr_continue^curr_segment
  }
  
  # The probability of continuing at the last segments is, by definition 0. Adjust that here.
  pr_continue_by_segment[max_n_segments] <- 0;
  
  excess_prob <- 1 - sum(pr_reject_by_segment) - sum(pr_ftr_by_segment)
  pr_reject_by_segment[max_n_segments] <- pr_reject_by_segment[max_n_segments] + excess_prob
  
  # Total probabilitiy of rejecting is the sum pr(reject) across all segments
  pr_reject_total <- sum(pr_reject_by_segment)
  
  
  # By definition
  expected_n_segments <- (1 - pr_continue^max_n_segments)/(1-pr_continue)
  
  # This vector holds the probability of terminating at each segment, by either decision
  pr_termination_by_segment <- pr_reject_by_segment + pr_ftr_by_segment
  
  # By definition
  segment_ord_squared <- (1:max_n_segments)^2
  expected_n_segments_squared <- sum(pr_termination_by_segment * segment_ord_squared)
  variance_n_segments <- expected_n_segments_squared - expected_n_segments^2
  
  result_list <- list(pr_reject_total = pr_reject_total,
                      pr_reject_by_segment = pr_reject_by_segment,
                      pr_ftr_by_segment = pr_ftr_by_segment,
                      pr_continue_by_segment = pr_continue_by_segment,
                      expected_n_segments = expected_n_segments,
                      variance_n_segments = variance_n_segments)
  
  return(result_list)
  
} # end outcome_probabilities

#=========================================
# run_int_subjects_scenario: Run scenario with an integer number of subjects
#=========================================

# This method generates a description of the probabalistic results of a single segmented
# hypothesis testing scenario for a given "researcher" (alphaweak,alphastrong and kmax),
# statistical test and one or more effect sizes whose total probability sums to 1. It returns
# an instance of SegmentedHypTestResult, which is effectivley a table of the useful descriptive
# metrics (see source file).

# This method requires that the provided SegmentedResearcher's number of subjects is an integer. See wrapper
# run_any_subjects_scenario for handling of non-integer cases that can be encountered, for example, during simulations.

run_int_subjects_scenario  <- function(segmented_researcher, stat_procedure, true_effect)
{
  
  # Create a skeleton instance to hold the results. Number of segments is provided to accomodate
  # the output vectors which hold a value for each segment (e.g. pr(reject) at segment i)
  results <- SegmentedHypTestResult(segmented_researcher$max_n_segments) # all other properties default to 0
  
  
  # A SegmentedResearcher object specified the number of subjects to be used for each segment
  # totalled across all groups. Thus the actual size per group depends on the number of groups under
  # the statistical procedure being used. For example, a one-sample t places all subjects in a segment
  # into a single group while a two-sample t puts half in each group. The statistical procedure classes
  # expose method divide_total_sample to perform this calculation
  
  n_per_group <- divide_total_sample(stat_procedure, segmented_researcher$n_per_segment)
  
  # Query the statistical procedure instance to compute the probabilities of "beating" (i.e. getting a
  # smaller p observed value than) alpha_strong for the effect size(s). This is logically equivalent to power
  # when Ho is false or the Type I error rate when Ho is true, those outcomes associated with rejecting Ho
  
  # These values are required to compute the prob of rejecting, failing to reject, or continuing at each
  # segment. When H0 is true (effect size = 0), they are equal to alpha_strong and alpha_weak, respectively.
  # When Ho is false they depend on the statistical procedure, the effect size and number of subjects per group.
  
  pr_beat_strong <- compute_power(stat_procedure,
                                  segmented_researcher$alpha_strong,
                                  true_effect$effect_size,
                                  n_per_group)
  
  # Perform the equivalent computation for alpha_weak. Failure to beat alpha_weak triggers failure to reject.
  
  pr_beat_weak <- compute_power(stat_procedure,
                                segmented_researcher$alpha_weak,
                                true_effect$effect_size,
                                n_per_group)
  
  # Use these values to compute the various outcome probabilities (for reject, ftr and continue) for each
  # segment.
  
  outcomes <- outcome_probabilities(segmented_researcher$max_n_segments, pr_beat_strong, pr_beat_weak)
  
  # The outcome_probabilities method returns a list with the fields as shown. We store some of
  # these in the skeleton results object.
  
  #   pr_reject_total: Total probability of rejecting at any step during the experiment
  #   pr_reject_by_segment: Probability of stopping after step i & rejecting.
  #   pr_ftr_by_segmenmt: Probability of stopping after step i & failing to reject.
  #   pr_continue_by_segment: Probability of continuing (neither reject or ftr) at step i
  #   expected_n_segments: The mean of the number of steps taken.
  #   variance_n_segments: The variance of the number of steps taken.
  
  results$pr_reject_total <- outcomes$pr_reject_total
  results$pr_reject_by_segment <- outcomes$pr_reject_by_segment
  results$pr_ftr_by_segment <- outcomes$pr_ftr_by_segment
  results$pr_continue_by_segment <- outcomes$pr_continue_by_segment
  
  # Set up values of alpha-beta matrix and store in result
  if (true_effect$effect_size == 0) {
    alpha_beta <- c(0, results$pr_reject_total, 1 - results$pr_reject_total, 0)
  } else {
    alpha_beta <- c(results$pr_reject_total, 0, 0, 1 - results$pr_reject_total)
  }
  
  results$pr_true_pos <- alpha_beta[1] # R vectors are 1-indexed
  results$pr_false_pos <- alpha_beta[2]
  results$pr_true_neg <- alpha_beta[3]
  results$pr_false_neg <- alpha_beta[4]
  
  # With all these important values computed, the expected numbers of subjects and expected effect size
  # for cases where Ho is rejected can be computed, by the magic of maths, as follows (pers comm, Miller, 2020):
  
  exp_n_subj <- outcomes$expected_n_segments * segmented_researcher$n_per_segment
  exp_n_segments_sqr <- outcomes$variance_n_segments + outcomes$expected_n_segments^2
  exp_n_subj_sqr <- exp_n_segments_sqr * segmented_researcher$n_per_segment^2
  
  # Assign the relevant values to the results object
  results$exp_n_subj <- exp_n_subj
  results$exp_n_subj_sqr <- exp_n_subj_sqr
  
  # Expected size of significant effect, conditional on rejection. Useful for simulations onle.
  
  # Find the expected effect size (mean of effect size across infinitely many iterations) conditional on
  # a significant result, for a given researcher object (kmax, nsubjects, alphastrong and alpha weak), and
  # statistical procedure. Simulations using this computation have shown that the expected effect size is
  # generally larger than the true effect size. That is, looking only at those experiments where Ho is
  # rejectd (as is typical in the publication process) leads to a general overestimation of true effect size.
  
  # Adjust sample size for number of groups per experiment (e.g. one-sample has 1, two-sample has 2)
  n_per_group <- divide_total_sample(stat_procedure, segmented_researcher$n_per_segment)
  
  # Compute probability-weighted average of expected significant effect sizes across all possible strong cut-off stopping segments.
  # Segments preceding the strong cut-off stop have p's bounded between alpha weak and alpha strong.
  
  
  # Compute the three expected effect sizes -- against alpha weak, between alpha weak and strong, and against alpha strong
  exp_effect_size_weak <- expected_significant_effect(stat_procedure, segmented_researcher$alpha_weak, n_per_group, true_effect$effect_size)
  exp_effect_size_between <- expected_significant_effect_bounded(stat_procedure, segmented_researcher$alpha_strong, segmented_researcher$alpha_weak, n_per_group, true_effect$effect_size)
  exp_effect_size_strong <- expected_significant_effect(stat_procedure, segmented_researcher$alpha_strong, n_per_group, true_effect$effect_size)
  
  # The weak result in the baseline nonsegmented case
  if ((segmented_researcher$alpha_strong == 0) || (segmented_researcher$max_n_segments == 1)){
    exp_effect_size_overall <- exp_effect_size_weak
  } else {
    # Iteratively determine the expected effect size for rejection at each of the kmax segments
    
    # Prepare skeleton to hold results
    exp_effect_size_by_segment <- numeric(segmented_researcher$max_n_segments)
    
    for(curr_segment in 1:segmented_researcher$max_n_segments - 1)
    {
      exp_effect_size_by_segment[curr_segment] <- (exp_effect_size_strong + ((curr_segment-1) * exp_effect_size_between))/curr_segment
    }
    exp_effect_size_by_segment[segmented_researcher$max_n_segments] <- (exp_effect_size_weak + ((segmented_researcher$max_n_segments -1)* exp_effect_size_between))/segmented_researcher$max_n_segments
    
    # Compute the overall expected significant effect size as the weighted average of each segment's expected effect size and probability
    # These probabilities have been computed above.
    exp_effect_size_overall <- sum(exp_effect_size_by_segment * outcomes$pr_reject_by_segment)/ outcomes$pr_reject_total
  } # end else
  
  
  # Store the value in the results data structre
  results$exp_effect_size <- exp_effect_size_overall
  
  return(results)
  
} # end run_int_subjects_scenario

#=========================================
# run_any_subjects_scenarioL Run scenario with integer or real number of subjects. Wrapper for
# run_int_subjects_scenario
#=========================================
# This function handles cases in which Researcher.NperSegment is a real number (i.e., need not be an integer). In that
# case, the scenario is simulated once with the floor and once with the ceiling of NperSegment, and the weighted average
# of those two results is returned. Weighting based on fractional part of original NPerSegment. For example, if nPerSegment
# is 5.14, simulate with 5 and 6, then take their weighted average with weights of .86 and .14 respectively.

run_any_subjects_scenario <- function(segmented_researcher, stat_procedure, true_effect)
{
  # Get number of subjects to use
  n_per_segment <- segmented_researcher$n_per_segment
  
  # If n is an integer, simply run the integer scenarion
  if (n_per_segment %% 1 == 0)
    result <- run_int_subjects_scenario(segmented_researcher, stat_procedure, true_effect)
  else
  {
    # If n has a fractional part, do the adjustment described above
    
    # Get the integer bounds
    lower_n <- floor(n_per_segment)
    upper_n <- ceiling(n_per_segment)
    
    # Assign the fractional part and 1- fractional part as probabilities for the weighted averaging
    pr_upper_n = n_per_segment - lower_n;
    pr_lower_n = 1 - pr_upper_n;
    
    # We will run the scenario twice -- once with n = lower and once with n = upper. In each case, we
    # need to pass the researcher into the run method, so we need to be able to change its n_per_segment
    # value. To avoid risk of destroying it in some way, we'll make a copy here.
    
    researcher_data_values <- unclass(segmented_researcher)
    # SegmentedResearcher <- function(max_n_segments, n_per_segment, alpha_total, alpha_strong = 0, alpha_weak = NULL)
    temp_segmented_researcher <- SegmentedResearcher(researcher_data_values$max_n_segments,
                                                     researcher_data_values$n_per_segment,
                                                     researcher_data_values$alpha_total,
                                                     researcher_data_values$alpha_strong)
    
    # run the lower n
    temp_segmented_researcher$n_per_segment <- lower_n
    lower_result <- run_int_subjects_scenario(temp_segmented_researcher, stat_procedure, true_effect)
    
    # run the upper n
    temp_segmented_researcher$n_per_segment <- upper_n
    upper_result <- run_int_subjects_scenario(temp_segmented_researcher, stat_procedure, true_effect)
    
    # Compute the weighted average of the two results. See method definition below.
    result <- weighted_avg_hyp_test_result(list(lower_result, upper_result), list(pr_lower_n, pr_upper_n))
    
  } # end else
  return(result)
} # end run_any_subjects_scenario

#=========================================
# Accepts a list of SegmentedHypTestResults and a parallel list of probabilities.
# Returns a SegmentedHypTestResult instance whose values are the weighted (by prob)
# average of the input objects.
#=========================================
weighted_avg_hyp_test_result <- function(results_list, probs_list)
{
  if (length(results_list) != length(probs_list)){
    stop("Results and probabilities lists must be the same length")
  }
  
  # We can pull representative values from any element of results_list. We use the 1st one, as per convention
  
  # Get the number of segments
  max_n_segments <- results_list[[1]]$max_n_segments
  
  # Get the names of all the fields, for associative indexing
  field_list <- names(results_list[[1]])
  
  # Fields 1 to 10 ields are scalars. we can process them in a single loop.
  # The remaining fields are vectors, which are processed separately
  
  n_scalars = 10
  scalar_list <- field_list[1:n_scalars]
  
  
  # Prepare an empty instance to hold the summary values. Missing args all default to 0
  weighted_avg_results <- SegmentedHypTestResult(max_n_segments)
  
  # iterate over the scalar fields, taking the weighted mean of each value, and inserting into the
  # prepared results object
  for (scalar_field_name in scalar_list)
  {
    
    # Pull the values for this scalar field into a list
    field_values <- sapply(results_list, function(result) result[scalar_field_name])
    
    # stats::weighted.mean wants vectors, not lists, so we need to unlist here.
    weighted_avg <- stats::weighted.mean(unlist(field_values), unlist(probs_list))
    
    # Store in skeleton
    weighted_avg_results[scalar_field_name] <- weighted_avg
  }
  
  # Iterate over vector fields
  # Each result has three vector fields for the  probabilities of rejecting, failing to reject and
  # continuing. Each of these vectors has one element for each segment.
  # Here we want to take a weighted average of those values, using the probability of each result.
  # So, for example, assume that we have 5 results. Let kmax = 3, so each of the five result objects
  # has a 3-element vector of pr(reject) for segments 1, 2 and 3, and also has a probability of occurrence.
  # We take the pr(reject vector) of result 1 and multiple each element vector-wise by the probability of
  # occurence of result 1. We then take result 2 and multiply its vector elements by its probability, and
  # WE SUM THE VALUES FROM 1 AND 2 VECTOR-WISE, obtaining a new vector, still of length three. We continue
  # over the remaining results, accruing the sum for each. The resulting 3-element vector which contains the
  # sum of five weighted values is the result to return.
  
  n_results <- length(results_list) # How many results have been passed in
  
  probs_vector <- unlist(probs_list)
  for (i in 1:n_results)
  {
    wght_reject <- results_list[[i]]$pr_reject_by_segment * probs_vector[i]
    weighted_avg_results$pr_reject_by_segment =  weighted_avg_results$pr_reject_by_segment + wght_reject
    
    wght_ftr <- results_list[[i]]$pr_ftr_by_segment * probs_vector[i]
    weighted_avg_results$pr_ftr_by_segment = weighted_avg_results$pr_ftr_by_segment + wght_ftr
    
    wght_cont <- results_list[[i]]$pr_continue_by_segment * probs_vector[i]
    weighted_avg_results$pr_continue_by_segment = weighted_avg_results$pr_continue_by_segment + wght_cont
  }
  
  return(weighted_avg_results)
  
} # end weighted_avg_hyp_test_result


##########################################
# Generics for public-facing methods
##########################################
#-----------------
# Roxygen comments

#' run_scenario.SegmentedHypTestEngine
#'
#' \code{run_scenario.SegmentedHypTestEngine} numerically determines expected
#' outcomes for a complete study scenario. It wraps various internal methods of
#' class SegmentedHypTestEngine (developers can cf. source code).
#' @param segmented_hyp_test_engine A SegmentedHypTestEngine instance for method
#'   dispatch
#' @param segmented_researcher Correctly initialised instance of
#'   SegmentedResearcher
#' @param stat_procedure Correctly initialised instance of any class descended
#'   from StatProcedureBase, e.g. OneSampleT, PearsonR, etc.
#' @param effects_list A list of TrueEffect instances. This must be a list, not
#'   simply a vector. Cast with list() if necessary. Each TrueEffect instance
#'   holds fields effect_size and effect_size_probability. A TrueEffect with
#'   effect_size 0 (i.e. H0 is true) may be included. Probabilities across all effect
#'   sizes must sum to 1 or an exception is thrown.
#' @return A SegmentedHypTestResult instance which holds the average expected outcomes
#' across all effect sizes, weighted by associated probability.
#'
#' @export
run_scenario <- function(segmented_hyp_test_engine,
                         segmented_researcher,
                         stat_procedure,
                         effects_list)
{
  UseMethod("run_scenario", segmented_hyp_test_engine)
}


#-----------------
# Roxygen comments

#' average_power.SegmentedHypTestEngine
#'
#' \code{average_power.SegmentedHypTestEngine} computes the average power across
#' a set of studies as defined by the effects_list argument. Needs to run from
#' the raw data, not SegHTResults, as SegHTResults do not carry their original
#' true effect size(s). Called by user-facing methods in SegmentedHypTester.R.
#' See source code for logic.
#' @param segmented_hyp_test_engine A SegmentedHypTestEngine instance for method
#'   dispatch
#' @param max_n_segments Maximum number of segments to be run before termination
#' @param n_per_segment Number of subjects to be tested in each segment
#' @param alpha_total Desired total Type 1 error probability across the entire
#'   experiment
#' @param alpha_strong For decision rule: If p-observed for a segment is less
#'   than alpha strong, reject and stop.
#' @param stat_procedure Character code for hypothesis testing procedure to
#'   be used in the experiment. Currently must be one of {'1t', '2t', 'r', '1z',
#'   '2z'}
#' @param effects_list List of TrueEffect instances

#' @return Average power across all studies based on effects in effect_list
#'
#' @export
average_power <- function(segmented_hyp_test_engine,
                          max_n_segments,
                          n_per_segment,
                          alpha_total,
                          alpha_strong,
                          stat_procedure,
                          effects_list)
{
  UseMethod("average_power", segmented_hyp_test_engine)
}

# Class SegmentedHypTestResult

# This class holds information about the expected outcomes for a segmented hypothesis testing experiment for given
# SegmentedHypTestResult (kmax, alphaTotal, alphaWeak and alphaStrong) & StatType (e.g. one sample t-test or correlation r).

# Can be produced either for a single true effect size or for a probability-weighted average of
# multiple true effect sizes which may include and effect size of 0 (i.e. H0 is true) with associated baserate probability.

# Note that some values are vectors holding a result for each of the possible ending segments 1 to k.

#############################################################
# Class properties - all default to 0 at creation
#############################################################
# max_n_segments:         maximum number of segments that can be run in the described experiment before termination
# pr_true_pos:            hit - reject Ho when it is false (there is an effect)
# pr_false_pos:           false alarm - reject Ho when it is true (there is no effect)
# pr_true_neg:            correct failure to reject
# pr_false_neg:           miss - failure to reject when there is really an effect
# pr_reject_total:        probability of rejecting (correctly or incorrectly) at any of the k segments
# exp_n_subj:             expected (mean across inf iterations) number of subjects used when testing for this effect.
# exp_n_subj_sqr:         expected (mean across inf iterations) square of the number of subjects used. Useful for variance computations.
# exp_effect_size:        expected (mean across inf iterations) observed effect size (not generally equal to the true effect size)
# pr_reject_by_segment    vector giving prob of rejecting at each segment
# pr_ftr_by_segment       vector giving probability of stopping with FTR at each segment
# pr_continue_by_segment  vector giving probability of continuing at each segment
#
# n_simulations           0 for computed results, or number of simulations for simulated results (not used in this project)

# The public-facing constructor requires only max_n_segments, and sets all other fields to 0
# but if you want to pack an instance with a set of prepared result values, you can assign to l0 the
# scalar input args.

#####################################################################
# Creation methods - following Wickham
#####################################################################
new_SegmentedHypTestResult <- function(max_n_segments,
                                       pr_true_pos = 0,
                                       pr_false_pos = 0,
                                       pr_true_neg = 0,
                                       pr_false_neg = 0,
                                       pr_reject_total = 0,
                                       exp_n_subj = 0,
                                       exp_n_subj_sqr = 0,
                                       exp_effect_size = 0,
                                       n_simulations = 0)
{
  
  new_segmented_hyp_test_result <- list(max_n_segments = max_n_segments,
                                        pr_true_pos = pr_true_pos,
                                        pr_false_pos = pr_false_pos,
                                        pr_true_neg = pr_true_neg,
                                        pr_false_neg = pr_false_neg,
                                        pr_reject_total = pr_reject_total,
                                        exp_n_subj = exp_n_subj,
                                        exp_n_subj_sqr = exp_n_subj_sqr,
                                        exp_effect_size = exp_effect_size,
                                        n_simulations = n_simulations,
                                        pr_reject_by_segment = numeric(max_n_segments),
                                        pr_ftr_by_segment = numeric(max_n_segments),
                                        pr_continue_by_segment = numeric(max_n_segments))
  
  
  
  attr(new_segmented_hyp_test_result, "class") <- "SegmentedHypTestResult"
  
  # Return instance
  return(new_segmented_hyp_test_result)
}


#----------------------------------------------------------------------
validate_SegmentedHypTestResult <- function(segemented_hyp_test_result)
{
  # Perform additional checking here, if wanted.
}


#---------------------------------------------------------------------------
#-----------------
# Roxygen comments

#' SegmentedHypTestResult Constructor
#'
#' A SegmentedHypTestResult instance holds information about the expected
#' outcomes of a segmented hypothesis testing study for given maximum number of
#' segments, alpha total, alpha strong, and statistical procedure.These objects
#' are typically created with no parameter values supplied so all fields default
#' to 'empty', and initisalised during computation to hold the simulated or
#' numerically-derived expected outcomes for the scenario. See
#' segHT::run_scenario.SegmentedHypTestEngine for an example.
#'
#' @param max_n_segments Maximum number of segments that can be run in the
#'   described experiment before termination
#' @param pr_true_pos Overall probability of a hit - reject Ho when it is false
#'   (there is an effect)
#' @param pr_false_pos Overall probability of a false alarm - reject Ho when it is
#'   true (there is no effect)
#' @param pr_true_neg Overall probability of a correct failure to reject
#' @param pr_false_neg Overall probability of a miss - failure to reject when there is really
#'   an effect
#' @param pr_reject_total Overall probability of rejecting (correctly or
#'   incorrectly) at any of the k segments
#' @param exp_n_subj Expected (mean across infinite iterations) number of
#'   subjects used when testing for this effect.
#' @param exp_n_subj_sqr Expected (mean across infinite iterations) square of
#'   the number of subjects used. Useful for variance computations.
#' @param exp_effect_size Expected (mean across infinite iterations) observed
#'   effect size (not generally equal to the true effect size)
#' @param n_simulations 0 for computed results, or the number of
#'   simulations for simulated results (to allow extension to simulation)
#' @return SegmentedHypTestResult instance
#'
#' @export
SegmentedHypTestResult <- function(max_n_segments,
                                   pr_true_pos = 0,
                                   pr_false_pos = 0,
                                   pr_true_neg = 0,
                                   pr_false_neg = 0,
                                   pr_reject_total = 0,
                                   exp_n_subj = 0,
                                   exp_n_subj_sqr = 0,
                                   exp_effect_size = 0,
                                   n_simulations = 0)
{
  
  # Create instance
  instance = new_SegmentedHypTestResult(max_n_segments,
                                        pr_true_pos,
                                        pr_false_pos,
                                        pr_true_neg,
                                        pr_false_neg,
                                        pr_reject_total,
                                        exp_n_subj,
                                        exp_n_subj_sqr,
                                        exp_effect_size,
                                        n_simulations)
  
  # Validate instance. Execution halts here if things go wrong
  validate_SegmentedHypTestResult(instance)
  
  # Return validated instance
  return(instance)
}


##########################################
# Public methods
##########################################
#=========================================
# print.SegmentendHypTestResult: System will dispatch to this method on print(instance)
#=========================================
#' @export
print.SegmentedHypTestResult <- function(x, ...)
{
  data_values <- unclass(x) # This grabs the original list from the struct built when attr was called
  
  # Print only the scalars, since this is primarily for debugging anyway. If you include the
  # vectors, R replicates the scalars to fill the data frame and it is confusing.
  # The scalars are the first 10 fields.
  
  n_scalar_fields <- 10
  df <- data.frame(data_values[1:n_scalar_fields])
  print(df)
}

# For maximum ease of access, these end-user facing methods are not contained
# in a class, so they need no generics or class instance for dispatch.
# Top level methods for non-programming users.

#' @importFrom pracma fzero
#' @importFrom pracma fminbnd

# Ignore this 13-03-20
placeholder <- function()
{
}

##########################################
# alpha_weak: User-facing wrapper for alpha_weak_computation
##########################################
#-----------------
# Roxygen comments

#' alpha_weak
#'
#' \code{alpha_weak}  returns the value of alpha_weak which
#' will produce the correct overall Type 1 error probability (i.e. alpha_total)
#' for a segmented hypothesis testing design, given values for maximum
#' number of segments, alpha total and alpha strong. This function wraps the
#' alpha_weak_computation method of the SegmentedResearcher class. If alpha
#' strong has been specified as zero (boundary case) alpha weak has an
#' analytical solution alpha_total^(1/max_n_segments). Otherwise alpha weak is
#' found via numerical search using function pracma::fzero. See Miller and
#' Ulrich, 2020 for discussion.
#'
#' @param max_n_segments Maximum number of segments to be run before termination
#' @param alpha_total Desired total Type 1 error probability across the entire
#'   study
#' @param alpha_strong For decision rule: If p-observed for a segment is less
#'   than alpha strong, stop and reject Ho.
#'
#' @return Numeric. Value of alpha weak.
#'
#' @export
alpha_weak <- function(max_n_segments, alpha_total, alpha_strong)
{
  # We create a SegmentedResearcher to call its alpha_weak method.
  segmented_researcher <- SegmentedResearcher()
  alpha_weak <- alpha_weak_computation(segmented_researcher, max_n_segments, alpha_total, alpha_strong)
  
  return(alpha_weak)
}

##########################################
# n_for_power: User-facing wrapper for n_for_power_computation
##########################################

#-----------------
# Roxygen comments

#' n_for_power
#'
#' \code{n_for_power}  returns the number of subjects
#' required to achieve a target statistical power for a segmented
#' hypothesis testing design, given the other design descriptors (see Arguments,
#' below). This method wraps n_for_power_computation to prevent the end user from
#' having to create his/her own component class instances. The result is
#' determined by numerical search using pracma::fminbnd. The function minimised
#' computes the difference between the power obtained by a given value of n and
#' the target power. The power obtained for any given value of n is computed
#' numerically. See n_for_power_computation in this file and
#' seght::run_scenario.SegmentedHypTestEngine for computation logic.
#'
#' @param target_power Desired power level (probability of detecting a false H0
#'   when present) across all segments of the study
#' @param stat_procedure_name Character code for hypothesis testing procedure to
#'   be used in the study. Currently must be one of {'1t', '2t', 'r', '1z',
#'   '2z'}
#' @param effect_size Size of the true effect when it is present.
#'   For tests of correlation, this is the true correlation rho.
#'   For tests of means (1t, 2t, 1z, 2z), this Cohen's d.
#' @param max_n_segments Maximum number of segments to be run before termination
#' @param alpha_total Desired total Type 1 error probability across the entire
#'   study
#' @param alpha_strong For decision rule: If p-observed for a segment is less
#'   than alpha_strong, stop and reject Ho.
#' @param alpha_weak For decision rule: If p-observed for a segment is greater
#'   than alpha weak, stop and fail to reject Ho. This argument is not normally provided,
#'   in which case the correct value is computed in the method. If the argument is
#'   provided, it is used without being checked (for speed) so it is the user's
#'   responsibility to ensure that it is the correct value.
#'
#' @return Returns the number
#'   of subjects per segment that yields the target power. The result is not
#'   necessarily a whole number, in which case you must use the next
#'   smaller or larger integer and get a little lower or higher power.
#' @export
n_for_power <- function(target_power,
                        stat_procedure_name,
                        effect_size,
                        max_n_segments = 3,
                        alpha_total = 0.5,
                        alpha_strong = 0.025,
                        alpha_weak = NA)
{
  
  # Assign a generous starting guess for n_per_segment. Computation code will raise
  # this if necessary to gain the required power level.
  n_per_segment <- 50
  
  # Build necessary classes for computation
  segmented_researcher <- SegmentedResearcher(max_n_segments, n_per_segment, alpha_total, alpha_strong, alpha_weak)
  stat_procedure <- stat_procedure_factory(stat_procedure_name)
  
  # Call the computation method
  # This method returns a three-item vector [n per segment, exact power, computed alpha weak]
  n_for_power_outcomes <- n_for_power_computation(target_power,
                                                  segmented_researcher,
                                                  stat_procedure,
                                                  effect_size)
  
  n_for_power <- n_for_power_outcomes[1] # First item of the three returned by no_for_power_computation
  
  # return the result.
  return(n_for_power)
  
} # n_for_power wrapper function


###########################################################################################
# segmented_hyp_test_outcomes: User-facing wrapper for SegmentedHypTestEngine scenario methods.
#
# User supplies scalar descriptors for an experimental scenario (kmax, total alpha, alpha strong,
# statistical test, total n per segment, true effect size(s) and base rate(s), and receives
# a summary of the expected experimental outcomes (alpha weak, average power, expected number of subjects used,
# pr(stop & reject) and pr(stop and fail to reject) at each segment.
###########################################################################################
#-----------------
# Roxygen comments

#' segmented_hyp_test_outcomes
#'
#' \code{segmented_hyp_test_outcomes}  wraps
#' seght::run_scenario.SegmentedHypTestEngine to numerically determine the
#' probabilities of all outcomes (reject, fail to reject, or continue) for each
#' potential segment of a segmented hypothesis testing design, along with a
#' variety of useful summary descriptors of expected study performance.
#'
#' @param max_n_segments Maximum number of segments to be run before termination
#' @param n_per_segment Number of subjects to be tested in each segment
#' @param alpha_total Desired total Type 1 error probability across the entire
#'   experiment
#' @param alpha_strong For decision rule: If p-observed for a segment is less
#'   than alpha strong, stop and reject Ho.
#' @param stat_procedure_name Character code for hypothesis testing procedure to
#'   be used in the experiment. Currently must be one of {'1t', '2t', 'r', '1z',
#'   '2z'}
#' @param effect_sizes Value, or vector of values, representing the true effect when it is present.
#'   For tests of correlation, this is the true correlation rho.
#'   For tests of means (1t, 2t, 1z, 2z), this Cohen's d.
#'   Any effect size of zero is logically Ho = true.
#' @param base_rates Base rate or vector of base rates giving the probability of each effect size.
#'   Each base_rate value matches the corresponding element of the effect_sizes argument, above.
#'   Each base rate must be between 0 and 1, and their total should be less than or equal to 1.
#'   If the total of base_rates is less than 1, the
#'   remaining probability is allocated to effect size 0 (i.e. H0 is true)
#'
#' @return Returns a list (collection of named items) containing all of the method
#'   input values (for reference), plus: \itemize{ \item exp_n_sub: Expected (mean
#'   across infinite iterations) number of subjects used in the study \item
#'   avg_power: Power level of the study on average across all cases where Ho is false
#'   (weighted average according to base rates).
#'   \item pr_reject_by_segment: Probablility of rejecting H0 at each segment
#'   \item pr_ftr_by_segment: Probability of failing to reject H0 at each
#'   segment }
#'
#' @export
segmented_hyp_test_outcomes <- function(max_n_segments,
                                        n_per_segment,
                                        alpha_total,
                                        alpha_strong,
                                        stat_procedure_name,
                                        effect_sizes,
                                        base_rates = 1)
{
  
  # We are aiming for a call to run_scenario.SegmentedHypTestEngine. We prepare our arguments first.
  segmented_researcher <- SegmentedResearcher(max_n_segments, n_per_segment, alpha_total, alpha_strong)
  stat_procedure <- stat_procedure_factory(stat_procedure_name)
  
  # run_scenario takes a list of true effects whose probabilities sum to 1. If the user has passed
  # in a baserate or vector of baserates that do not sum to 1, we add an element with
  # effect size 0 (i.e. H0 true) with base rate as required to get the total to 1
  
  if (length(base_rates) != length(effect_sizes)) {
    stop("must have equal numbers of effect_sizes and base_rates")
  }
  
  total_effect_prob <- sum(base_rates)
  
  if (total_effect_prob > 1) {
    stop("sum of base_rates cannot exceed 1.0")
  }
  
  if (total_effect_prob < 1) {
    ho_base_rate <- 1 - total_effect_prob
    ho_effect_size <- 0
    effect_sizes <- c(effect_sizes, ho_effect_size)
    base_rates <- c(base_rates, ho_base_rate)
  }
  
  # Create TrueEffect objects for each effect_size/base_rate pair
  effects_list <- list(length(effect_sizes))
  for (i in 1:length(effect_sizes))
  {
    te <- TrueEffect(effect_sizes[i], base_rates[i])
    effects_list[[i]] <- te
  }
  
  # Make the call
  engine <- SegmentedHypTestEngine()
  experiment_results <- run_scenario(engine, segmented_researcher, stat_procedure, effects_list)
  
  # Compute the average power
  avg_power <-average_power(engine,
                            max_n_segments,
                            n_per_segment,
                            alpha_total,
                            alpha_strong,
                            stat_procedure,
                            effects_list)
  
  
  # Pull other necessary values to build the return list
  alpha_weak <- segmented_researcher$alpha_weak # This value initialised in the SegmentedResearcher ctor
  exp_n_subj <- experiment_results$exp_n_subj
  pr_reject_by_segment = experiment_results$pr_reject_by_segment
  pr_ftr_by_segment = experiment_results$pr_ftr_by_segment
  
  # Build it
  output_value <- list(max_n_segments = max_n_segments,
                       n_per_segment = n_per_segment,
                       alpha_total = alpha_total,
                       alpha_strong = alpha_strong,
                       alpha_weak = alpha_weak,
                       stat_procedure = stat_procedure_name,
                       effects_list = effects_list,
                       exp_n_subj = exp_n_subj,
                       avg_power = avg_power,
                       pr_reject_by_segment = pr_reject_by_segment,
                       pr_ftr_by_segment = pr_ftr_by_segment)
  
  # Return it
  return(output_value)
}


####################################################################################
# search_kmax: Accepts total alpha, alpha strong, a test name, an effect size with optional
# effect probability (baserate; defaults to 1 if not provided) and a target power level. Returns
# a dataframe with one row for each value of kmax between 2 and 10 (a practical maximum).
# Columns indicate the appropriate alpha weak, the n per segment required to achieve the target power,
# and the expected number of subjects and expected number of segments for the study
####################################################################################

#-----------------
# Roxygen comments

#' search_kmax
#'
#' \code{search_kmax} computes a set  of performance
#' descriptors (see Values, below) for a specified study scenario for each
#' value of kmax (maximum number of segments) between 2 and 10. Returned values
#' are based on the n_per_segment required to achieve the provided target statistical power.
#'
#' @param alpha_total Desired total Type 1 error probability across the entire
#'   study
#' @param alpha_strong For decision rule: If p-observed for a segment is less
#'   than alpha strong, stop and reject Ho.
#' @param stat_procedure_name Character code for hypothesis testing procedure to
#'   be used in the experiment. Currently must be one of {'1t', '2t', 'r', '1z',
#'   '2z'}
#' @param target_power Desired power level (probability of detecting a false H0
#'   when present) across all segments of the study
#' @param effect_size Size of the true effect when it is present.
#'   For tests of correlation, this is the true correlation rho.
#'   For tests of means (1t, 2t, 1z, 2z), this Cohen's d.
#' @param base_rate Base rate of the true effect_size. Must be between 0 and 1.
#'   If this value is less than 1, the remaining probability is allocated to
#'   effect size 0 (i.e. H0 is true)
#' @return Returns a data frame with one row for each value of kmax between 2
#'   and 10, and columns: \itemize{ \item kmax: Maximum number of segments in
#'   the study \item alpha_weak: For decision rule: If p-observed for a segment
#'   is greater than alpha weak, stop and fail to reject Ho. \item n per segment:
#'   Number of subjects per segment needed to achieve the target power. \item expected n subjects: Expected (mean)
#'   number of subjects required across infinite iterations of the experiment
#'   \item expected n segments: Expected (mean) number of segments that will be
#'   performed prior to termination (with any outcome) across infinite
#'   iterations of the experiment. }
#'
#' @export
search_kmax <- function(alpha_total,
                        alpha_strong,
                        stat_procedure_name,
                        target_power,
                        effect_size,
                        base_rate = 1)
{
  # Create engine instance to run scenarios
  engine <- SegmentedHypTestEngine()
  
  # Prepare skeleton vectors to hold outputs
  kmax_holder <- c()
  alpha_weak_holder <- c()
  n_per_segment_holder <- c()
  exp_n_subj_holder <- c()
  exp_n_segments_holder <- c()
  
  # Prepare required args
  stat_procedure <- stat_procedure_factory(stat_procedure_name)
  
  n_segment_upper_limit = 10 # By definition
  n_subjects_default = 50
  
  # Loop over potential values of kmax
  
  # For each potential value of kmax, we call n_for_power_computation (internal method of this class)
  # to find an appropriate n_per_segment to get the desired power.
  # We then run the scenario with that value to get expected n subjects and expected n segments.
  
  # fzero throws if max_n_segments is < 2
  for (max_n_segments in 2:n_segment_upper_limit)
  {
    # SegmentedResearcher$alpha_weak is initialised in the ctor
    segmented_researcher <- SegmentedResearcher(max_n_segments,
                                                n_subjects_default,
                                                alpha_total,
                                                alpha_strong)
    
    # n_for_power returns a vector [n, power]
    n_per_segment_output <- n_for_power_computation(target_power,
                                                    segmented_researcher,
                                                    stat_procedure,
                                                    effect_size)
    n_per_segment <- n_per_segment_output[1]
    
    # We wish to call run_scenario with this n_per_segment.
    # run_scenario takes a list of true effects whose probabilities sum to 1. If the user has passed
    # in a baserate, we create the appropriate TrueEffect instance for effect size 0 (i.e. H0 true)
    
    real_effect <- TrueEffect(effect_size, base_rate)
    if (base_rate == 1) {
      effects_list <- list(real_effect)
    } else {
      null_effect_probability <- 1 - base_rate
      null_effect <- TrueEffect(0, null_effect_probability)
      effects_list <- list(real_effect, null_effect)
    }
    
    # Remake the researcher with the n_per_segment computed by find_n_for_power
    segmented_researcher <- SegmentedResearcher(max_n_segments,
                                                n_per_segment,
                                                alpha_total,
                                                alpha_strong)
    
    experiment_results <- run_scenario(engine, segmented_researcher, stat_procedure, effects_list)
    
    # Expected number of segments is (expected total number of subjects)/(number of subjects per segment)
    exp_n_subjects <- experiment_results$exp_n_subj
    exp_n_segments <- exp_n_subjects/n_per_segment
    
    # Gather up the output values
    kmax_holder <- c(kmax_holder, max_n_segments)
    alpha_weak_holder <- c(alpha_weak_holder, segmented_researcher$alpha_weak)
    n_per_segment_holder <- c(n_per_segment_holder, n_per_segment)
    exp_n_subj_holder <- c(exp_n_subj_holder, exp_n_subjects )
    exp_n_segments_holder <- c(exp_n_segments_holder, exp_n_segments)
    
  } # end for all values of kmax
  
  # Build the output dataframe
  output_df <- data.frame(kmax_holder,
                          alpha_weak_holder,
                          n_per_segment_holder,
                          exp_n_subj_holder,
                          exp_n_segments_holder)
  
  colnames(output_df) <- c("kmax", "alpha_weak", "n_per_segment", "exp_n_subjects", "exp_n_segments")
  
  # Return it
  return(output_df)
}



#####################
# Internal methods
#####################
#=========================================
# stat_procedure_factory: Accepts string from {'1t', '2t', '1z', '2z', 'r'}. Returns initialised objects
# instance from StatProcedure family.
#=========================================
stat_procedure_factory <- function(stat_procedure_name)
{
  stat_names <- c("1t", "2t", "1z", "2z", "r")
  
  if (!any(stat_names == stat_procedure_name))
  {
    stop("stat_procedure_name must be one of {'1t', '2t', '1z', '2z', 'r'} ")
  }
  
  switch(stat_procedure_name,
         "1t" = stat_procedure <- OneSampleT(),
         "2t" = stat_procedure <- TwoSampleT(),
         "1z" = stat_procedure <- OneSampleZ(),
         "2z" = stat_procedure <- TwoSampleZ(),
         "r"  = stat_procedure <- PearsonR())
  
  return(stat_procedure)
}


#####################################################################################
# Find the N that is needed _PER SEGMENT_ to have the indicated TargetPower for the
# specified researcher, stat type and true effect
#####################################################################################
n_for_power_computation <- function(target_power, segmented_researcher, stat_procedure, effect_size, sample_max = 50)
{
  
  # Copy input researcher for use in computation
  researcher_data_values <- unclass(segmented_researcher)
  temp_segmented_researcher <- SegmentedResearcher(researcher_data_values$max_n_segments,
                                                   researcher_data_values$n_per_segment,
                                                   researcher_data_values$alpha_total,
                                                   researcher_data_values$alpha_strong)
  
  # Never a base rate in this computation, by definition
  true_effect <- TrueEffect(effect_size, effect_size_probability = 1)
  
  engine <- SegmentedHypTestEngine()
  
  
  # Find minimum sample size for this statistical procedure. Maximum is an input argument
  # that defaults to 50 (pers. comm. Miller, 2020)
  sample_min <- stat_procedure$min_sample_n
  
  # Make sure that the desired power is possible between sample min and max. If not, slide the
  # bounds up in sample_max size steps, until it is.
  found_bounds <- FALSE
  
  while (!found_bounds)
  {
    temp_segmented_researcher$n_per_segment <- sample_max
    temp_results <- run_scenario(engine, temp_segmented_researcher, stat_procedure, list(true_effect))
    if (temp_results$pr_reject_total < target_power) {
      sample_min <- sample_max
      sample_max <- sample_max * 2
    } else {
      found_bounds <- TRUE
    }
  } # not found_bounds
  
  # With your bounds in hand, use numerical search to generate the target sample size
  
  # Define the function for minimisation of distance between pr(reject) and target power.
  # Variables other than n must be in scope an initialised
  power_error_fn <- function(n)
  {
    temp_segmented_researcher$n_per_segment <- n
    result <- run_scenario(engine, temp_segmented_researcher, stat_procedure, list(true_effect))
    error <- abs(result$pr_reject_total - target_power)
  }
  
  # Here is the search. It generates messages about its algorithms, which we do not want.
  # We use capture.output to suppress these.
  capture.output(best_n <- pracma::fminbnd(power_error_fn, sample_min, sample_max, tol = 1e-2))
  
  # Compute the precise power for your guess
  temp_segmented_researcher$n_per_segment <- best_n$xmin
  final_result <- run_scenario(engine, temp_segmented_researcher, stat_procedure, list(true_effect))
  
  output_value <- c(best_n$xmin, final_result$pr_reject_total, temp_segmented_researcher$alpha_weak)
  return(output_value)
  
} # find N for power computation

# Class SegmentedResearcher: Holds collection of summary values to describe a specific research scenario



##########################################
# Class properties
##########################################
# Passed in to constructor

#   max_n_segments: Maximum number of segments to test.
#   n_per_segment:  SampleN per segment (total N across all samples if there are multiple samples per segment, e.g., 2-sample t-test).
#   alpha_total:    Total probability of Type I error across all segments.
#   alpha_strong:   Stop sampling and reject Ho if p<AlphaStrong in any segment. Defaults to 0 if passed null
#   alpha_weak:     Stop sampling Maximum p value for observed results in each segment. Computed if not provided


##########################################
# Creation methods - following Wickham
##########################################
#=========================================
# new_SegmentedResearcher: Default values determined empirically. Users are advised to override.
#=========================================
# Low-level constructor
new_SegmentedResearcher <- function(max_n_segments = 3,
                                    n_per_segment = 50,
                                    alpha_total = 0.05,
                                    alpha_strong = 0,
                                    alpha_weak = NA)
{
  
  new_seg_res <- list(max_n_segments = max_n_segments,
                      n_per_segment = n_per_segment,
                      alpha_total = alpha_total,
                      alpha_strong = alpha_strong,
                      alpha_weak = alpha_weak)
  
  
  attr(new_seg_res, "class") <- "SegmentedResearcher"
  
  # Return instance
  return(new_seg_res)
  
} # new_SegmentedResearcher
#=========================================

#=========================================
# validate_SegmentedResearcher
#=========================================
# Checks additional properties of elements, for example, that values are in allowed ranges -- business rules
validate_SegmentedResearcher <- function(segmented_researcher)
{
  # Pull field values for clarity
  max_n_segments <- segmented_researcher$max_n_segments
  alpha_total <- segmented_researcher$alpha_total
  alpha_strong <- segmented_researcher$alpha_strong
  alpha_weak <- segmented_researcher$alpha_weak
  
  # Make sure that alpha weak is initialised correctly.
  # If a value was supplied (e.g. to improve efficiency in simulations), make sure that it is
  # mathematically correct.
  # If no value was supplied, compute it here.
  
  # If supplied, check the maths (pers. comm. Miller, 2020)
  tolerance <- 0.0001
  if (!is.na(alpha_weak)) {
    
    alpha_dif <- alpha_weak - alpha_strong
    pred_alpha_total <- alpha_strong*(1 - alpha_dif^(max_n_segments-1)) / (1- alpha_dif) + alpha_weak * alpha_dif^(max_n_segments-1)
    
    if (abs(pred_alpha_total - alpha_total) > tolerance) {
      stop('Cannot proceed because provided values of alpha_strong, and alpha_weak do not produce alpha_total')
    }
  } else {
    
    # Computes appropriate alpha weak value to obtain desired alpha total. See method for computation.
    alpha_weak <- alpha_weak_computation(segmented_researcher, max_n_segments, alpha_total, alpha_strong)
    
    # Assign to instance. Objects are passed by reference in this construction, so this works.
    segmented_researcher$alpha_weak <- alpha_weak
  } # end alpha_weak passed as default NA
  # Return updated object instance
  return(segmented_researcher)
  
} # end validate_SegmentedResearcher
#=========================================


#=========================================
# SegmentedResearcher - public-facing ctor
#=========================================
#-----------------
# Roxygen comments

#' SegmentedResearcher Constructor
#'
#' A SegmentedResearcher instance holds a collection of values which
#' define a specific research scenario.
#'
#'
#' \code{SegmentedResearcher}  is the public-facing SegmentedResearcher constructor. It
#' provides typical default values for its parameters, but in most cases the user
#' will supply their own values to define the study scenario they wish to explore.
#'
#' @param max_n_segments An integer specifying the maximum number of segments to
#'   be run prior to termination
#' @param n_per_segment A \code{numeric} specifying the number of subjects to be
#'   tested in each segment. Note that this value should comprise all subjects
#'   in a segment, summed across any groups.
#' @param alpha_total Required Type 1 error probability across the entire study.
#'   Defaults to 0.05
#' @param alpha_strong For segmented hypothesis testing decision rule: If
#'   p-observed for any segment is less than alpha_strong, reject H0 and stop
#'   testing.
#' @param alpha_weak For decision rule: If p-observed for any segment is
#'   greater than alpha_weak, fail to reject H0 and stop testing. NB: The correct value
#'   of alpha_weak is determined entirely by the values of alpha_total and
#'   alpha_strong and must be derived numerically from those parameters. In most
#'   cases, a value SHOULD NOT be supplied for this argument. When no value is
#'   supplied, the correct value is automatically computed in succeeding
#'   computations.
#' @return Validated SegmentedResearcher instance
#'
#' @export
SegmentedResearcher <- function(max_n_segments = 3,
                                n_per_segment = 50,
                                alpha_total = 0.05,
                                alpha_strong = 0,
                                alpha_weak = NA)
{
  
  # Create instance
  instance <- new_SegmentedResearcher(max_n_segments, n_per_segment, alpha_total, alpha_strong, alpha_weak)
  
  # Validate instance. Initialises alpha weak if necessary.
  valid_instance <- validate_SegmentedResearcher(instance)
  
  # Return validated instance
  return(valid_instance)
  
} # end SegmentedResearcher
#=========================================


##########################################
# Public methods
##########################################
#=========================================
# print.SegmentedResearcher: System will dispatch to this method on print(instance)
#=========================================
#' @export
print.SegmentedResearcher <- function(x, ...)
{
  data_values <- unclass(x) # This grabs the original list from the struct built when attr was called
  df <- data.frame(data_values)
  print(df)
}


##########################################
# Internal methods
##########################################
#=========================================
# alpha_weaK_computation - For given values of max_n_segments and alpha_strong, compute the value of alpha_weak
# that will give the desired overall value of alpha_total (total Type 1 error rate). See Miller & Ulrich, 2020
# for detailed discussion.
#=========================================
alpha_weak_computation <- function(segmented_researcher, max_n_segments, alpha_total, alpha_strong, tolerance = 1e-8)
{
  # Logically, alpha_strong can not be larger than alpha_total
  if (alpha_strong > alpha_total)
  {
    stop("Alpha strong must be less than or equal to the overall alpha")
  }
  
  # If alpha_strong has been specified as zero (boundary case) alpha weak has an analytical solution. Otherwise
  # use numerical search function fzero to minimise function which computes the difference between alpha_total obtained
  # with a target value of alpha weak, and the desired alpha_total.
  
  if (alpha_strong == 0) {
    final_alpha_weak <- alpha_total^(1/max_n_segments)
  } else {
    
    # Error function to minimise. Computation by definition. Cf. Miller & Ulrich, 2020.
    compute_error <- function(try_alpha_weak)
    {
      diff <- try_alpha_weak - alpha_strong
      diff_accumulated  <- diff^(max_n_segments - 1)
      err <- alpha_strong * (1 - diff_accumulated )/(1 - diff )+try_alpha_weak * diff_accumulated - alpha_total
      return(err)
    }
    
    # Define sensible range for numeric search
    aw_range <- c(alpha_strong , alpha_total^(1/max_n_segments ))
    
    # Perform numeric search to minimise function
    search_result <- pracma::fzero(compute_error, aw_range )
    
    # Pull value of x that minimises error function
    final_alpha_weak <- search_result$x
  } # end if alpha_strong != 0
  
  
  return (final_alpha_weak)
}

# Class StatProcedureBase

# Base class for objects representing the statistical test being used in a research
# scenario. Children are 1 and 2 group t-tests, correlation and 1 and 2 group z-test
# following Miller & Ulrich 2020


##########################################
# Class properties: Values by definition for each child test, set in ctor
##########################################

# min_sample_n
# display_abbrev
# display_name


##########################################
# Creation methods (following Wickman)
##########################################

#=========================================
# new_StatProcedureBase: Children will, in typical cases, provide hard-coded values for all fields.
# Values for display_abbrev must match SegmentedHypTestEngine$stat_procedure_factory switch statement.
#=========================================
# Low-level constructor - following Wickman
new_StatProcedureBase <- function(min_sample_n, display_abbrev, display_name)
{
  
  new_stat_procedure_base <- list(min_sample_n = min_sample_n,
                                  display_abbrev = display_abbrev,
                                  display_name = display_name)
  
  
  attr(new_stat_procedure_base, "class") <- "StatProcedureBase"
  
  # Return instance
  return(new_stat_procedure_base)
  
} # new_StatProcedureBase
#=========================================

#=========================================
# validate_StatProcedureBase
#=========================================
# Checks additional properties of elements, for example, that values are in allowed ranges -- business rules
validate_StatProcedureBase <- function(stat_procedure_base)
{
  if (stat_procedure_base$min_sample_n < 1) {
    stop("Minimum sample size must be greater than 0")
  }
  
} # end validate_StatProcedureBase
#=========================================


#=========================================
# StatProcedureBase - Public ctor
#=========================================
#-----------------
# Roxygen comments

#' StatProcedureBase Placeholder Ctor
#'
#' StatProcedureBase is the abstract base class for all stat procedure classes.
#'
#' \code{StatProcedureBase}  should not be called directly. Descend all child classes
#' from StatProcedureBase by assigning the class attribute in the child constructor.
#' For example, the low-level constructor for class OneSampleT contains the statement:
#'
#' \code{attr(new_one_sample_t, "class") <- c("OneSampleT", "StatProcedureBase")}
#'
#' @param min_sample_n Smallest n allowed for this test. Supplied in constructor.
#' @param display_abbrev Two-character abbreviation for this test. Developers see stat_procedure_factory.SegmentedHypTester
#' @param display_name String value useful for generating tidy outputs for this test
#'
#'
#' @export
StatProcedureBase <- function(min_sample_n, display_abbrev, display_name)
{
  # Create instance
  instance = new_StatProcedureBase(min_sample_n, display_abbrev, display_name)
  
  # Validate instance. Execution halts here if things go wrong
  validate_StatProcedureBase(instance)
  
  # Return validated instance
  return(instance)
  
} # end StatProcedureBase ctor
#=========================================


##########################################
# Class methods
##########################################

#=========================================
# print.StatProcedureBase: System will dispatch to this method on print(instance)
#=========================================
#' @export
print.StatProcedureBase <- function(x, ...)
{
  data_values <- unclass(x) # This grabs the original list from the struct built when attr was called
  df <- data.frame(data_values)
  print(df)
}


# Remaining family methods are all abstract in the base type
#=========================================
# compute_power: Given a directional alpha, anticipated effect size and sample n, return power for the child test.
# These have analytical solutions for each test. See children for details
#=========================================
compute_power.StatProcedureBase <- function(stat_procedure, alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
{
  # Abstract
} # end compute_power


#=========================================
# p_level: Given an observed stat value and sample n, compute the p-level for the child test.
# These are computed by selecting from the appropriate sampling distribution
#=========================================
p_level.StatProcedureBase <- function(stat_procedure, observed_stat, sample_n)
{
  # Abstract
} # end p_level

#=========================================
# critical_value: Given alpha and sample n, compute the critical value for the child test.
# Determined by accessing the cumulative density function for the appropriate sampling distribution
#=========================================
critical_value.StatProcedureBase <- function(stat_procedure, alpha_one_tailed, sample_n)
{
  # Abstract
} # end critical_value

#=========================================
# noncentrality: Provides a measure of the extent to which the null hypothesis is false, determining the
# true distribution from which the observed statistic is drawn. These have analytical solutions for each
# test. See the child code.
#=========================================
noncentrality.StatProcedureBase <- function(stat_procedure, sample_n, effect_size)
{
  # Abstract
} # end noncentrality

#=========================================
# effect_size_from_stat: Measure of effect size derived from observed statistic value and n. These
# have analytical solutiosn. See the child code.
#=========================================
effect_size_from_stat.StatProcedureBase <- function(stat_procedure, stat_value, sample_n)
{
  # Abstract
} # end effect_size_from_stat

#=========================================
# expected_significant_effect: Determines the expected (mean) effect size, if one considers only those
# cases where the null hypothesis is reject (i.e. the observed stat is greater than the critical). For each
# test, this is computed by integrating across stat * stat_density for all values above the critical
# to obtain an expected (mean) value for that portion of the distribution, than dividing by the total area
# under the curve to obtain the correct proportional value. See child code.
#=========================================
expected_significant_effect.StatProcedureBase <- function(stat_procedure, alpha_one_tailed, sample_n, effect_size)
{
  # Abstract
} # end expected_significant_effect


#=========================================
# expected_significant_effect_bounded: As for expected_significant_effect, above, but integrated only
# over a region determined by two bounding alpha levels.
#=========================================
expected_significant_effect_bounded.StatProcedureBase <- function(stat_procedure, alpha_strong, alpha_weak, sample_n, effect_size)
{
  # Abstract
} # end expected_significant_effect


#=========================================
# divide_total_sample: Convert from n_per_segment (which is always the total subjects used per segment)
# to n_per_group, depending on the number of groups the child test applies to.
#=========================================
divide_total_sample.StatProcedureBase <- function(stat_procedure, n_per_segment)
{
  # Abstract
} # end divide_total_sample



##########################################
# Generics: Following the S3 method dispatch model. This causes each child instance to invoke its own
# polymorphic method implementation. See Wickman and others for details of the function binding process.
##########################################

#=========================================
# Generics -- inherited by child classes
#=========================================
#-----------------
# Roxygen comments

#' compute_power
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all children
#'
#' \code{compute_power} computes the statistical power of a statistical test
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param alpha_one_tailed Alpha level for one-tailed test
#' @param effect_size Presumed true effect size
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#' @param critical_value Optional to speed computation. If omitted, this value is computed in the method
#'
#' @return Computed statistical power (probability of rejecting H0 when the effect is present)
#'
#' @export
compute_power <- function(stat_procedure, alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
{
  UseMethod("compute_power", stat_procedure)
}

#-----------------
# Roxygen comments

#' p_level
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all children
#'
#' \code{p_level} computes the p_level for an observed statistical values and sample size
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param observed_stat Observed test value of interest
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#'
#' @return Computed p-level
#'
#' @export
p_level <- function(stat_procedure, observed_stat, sample_n)
{
  UseMethod("p_level", stat_procedure)
}

#-----------------
# Roxygen comments

#' critical_value
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all children
#'
#' \code{critical_value} computes the critical value of a test for given alpha and n
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param alpha_one_tailed Alpha level for one-tailed test
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#'
#' @return Critical value of test for alpha and n
#'
#' @export
critical_value <- function(stat_procedure, alpha_one_tailed, sample_n)
{
  UseMethod("critical_value", stat_procedure)
}

#-----------------
# Roxygen comments

#' noncentrality
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all children
#'
#' \code{noncentrality} computes the noncentrality parameter for a given n and effect size
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param effect_size Presumed true effect size
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#'
#' @return Computed noncentrality parameter
#'
#' @export
noncentrality <- function(stat_procedure, sample_n, effect_size)
{
  UseMethod("noncentrality", stat_procedure)
}

#-----------------
# Roxygen comments

#' effect_size_from_stat
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all
#' children
#'
#' \code{effect_size_from_stat} computes the effect size for a test given alpha,
#' n, and an observed value of the test statistic
#'
#' @param stat_procedure Instance for method dispatch. Member of
#'   StatProcedureBase family
#' @param stat_value Value of test statistic
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#'
#' @return Computed effect size
#'
#' @export
effect_size_from_stat <- function(stat_procedure,stat_value, sample_n)
{
  UseMethod("effect_size_from_stat", stat_procedure)
}

##-----------------
# Roxygen comments

#' expected_significant_effect
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all
#'children
#'
#' \code{expected_significant_effect} is used for simulations to demonstrate
#' bias in published observed effect sizes. Determines the expected (mean) effect
#' size, if one considers only those cases where the null hypothesis is rejected
#' (i.e. the observed stat is greater than the critical). For each test, this is
#' computed by integrating across stat * stat_density for all values above the
#' critical to obtain an expected (mean) value for that portion of the
#' distribution, than dividing by the total area under the curve to obtain the
#' correct proportional value. See child code.
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param alpha_one_tailed Alpha level for one-tailed test
#' @param effect_size Presumed true effect size
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#'
#' @return Expected (mean) effect size
#'
#' @export
expected_significant_effect <- function(stat_procedure, alpha_one_tailed, sample_n, effect_size)
{
  UseMethod("expected_significant_effect", stat_procedure)
}

##-----------------
## Roxygen comments

#' expected_significant_effect_bounded
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all
#'children
#'
#' \code{expected_significant_effect_bounded} as for expected_significant_effect but for
#' a specified portion of the sampling distribution
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param alpha_strong Upper bound in sampling distribution
#' @param alpha_weak Lower bound in sampling distribution
#' @param effect_size Presumed true effect size
#' @param sample_n Group sizes. Must be a vector for tests with multiple groups
#'
#' @return Expected (mean) effect size
#'
#' @export
expected_significant_effect_bounded <- function(stat_procedure, alpha_strong, alpha_weak, sample_n, effect_size)
{
  UseMethod("expected_significant_effect_bounded", stat_procedure)
}



## Roxygen comments

#' divide_total_sample
#'
#' Abstract method defined in StatProcedureBase. Should be overridden by all
#'children
#'
#' \code{divide_total_sample} returns n per group given total n
#'
#' @param stat_procedure Instance for method dispatch. Member of StatProcedureBase family
#' @param n_per_segment Total n per segment
#'
#' @return N per group for the associated statistical procedure object
#'
#' @export
divide_total_sample <- function(stat_procedure, n_per_segment)
{
  UseMethod("divide_total_sample", stat_procedure)
}

# Class OneSampleZ descends from StatProcedureBase

# See class OneSampleT for detailed explanations of method logic

# All the computations ported from original source code (Miller, 2020)

# Note that the sampling distribution for z, when Ho is true, is a normal with mean = 0 and s = 1.
# Noncentrality shifts the distribution mean depending on n and the true effect size

# Use
# dnorm(x, mean = 0, sd = 1, log = FALSE)
# pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
# qnorm(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)

#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats dnorm
#' @importFrom pracma integral
#' @importFrom utils capture.output

########################################
# Creation methods - following Wickham
########################################
#====================================================================
new_OneSampleZ <- function(min_sample_n, display_abbrev, display_name)
{
  new_one_sample_z <- list(min_sample_n = min_sample_n,
                           display_abbrev = display_abbrev,
                           display_name = display_name)
  
  
  attr(new_one_sample_z, "class") <- c("OneSampleZ", "StatProcedureBase")
  
  return(new_one_sample_z)
}
#====================================================================
validate_OneSampleZ <- function(one_sample_z)
{
  # Do any validation here. Use stop(msg) to exit if things go wrong
}
#====================================================================
#-----------------
# Roxygen comments

#' OneSampleZ Constructor
#'
#' OneSampleZ is a child of StatProcedureBase, representing a one-sample z-test
#'
#' \code{OneSampleZ}  is the public-facing OneSampleZ constructor. It creates and returns a class instance
#' with appropriately initialised values for all fields (min_sample_n, display_abbrev, and
#' display_name). It requires no parameters.
#'
#' @param min_sample_n Smallest n allowed for this test. Supplied in constructor.
#' @param display_abbrev Two-character abbreviation for this test. Developers see stat_procedure_factory.SegmentedHypTester
#' @param display_name String value useful for generating tidy outputs for this test
#'
#' @return Validated OneSampleZ instance
#'
#' @export
OneSampleZ <- function(min_sample_n = 2, display_abbrev = "1z", display_name = '1-sample z-test')
{
  instance <- new_OneSampleZ(min_sample_n, display_abbrev, display_name)
  validate_OneSampleZ(instance)
  return(instance)
}


#################################
# Instance methods
#################################
#=========================================
# compute_power(alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
#=========================================
#' @export
compute_power.OneSampleZ <- function(stat_procedure, alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
{
  z_crit_cumul_prob = 1 - alpha_one_tailed
  
  if (is.null(critical_value))
  {
    z_crit = qnorm(z_crit_cumul_prob)
  }
  else
  {
    z_crit = critical_value
  }
  
  noncentrality_parameter = noncentrality(stat_procedure, sample_n, effect_size)
  
  # When Ho is false, z-critical remains the same, but the true z-distribution shifts upward by the effect
  # size and sample size as reflected by the noncentrality parameter. The area that causes rejection is that
  # above z-critical in this new distribution, i.e. that above z-critical - the noncen parameter, reflecting the
  # shift.
  power = 1 - pnorm(z_crit - noncentrality_parameter)
  
  return(power)
  
} # end compute_power


#=========================================
# p_level(one_sample_z, observed_stat, sample_n)
#=========================================
#' @export
p_level.OneSampleZ <- function(stat_procedure, observed_stat, sample_n)
{
  
  # the observed sample exceeds this proportion of the normal when Ho is true
  cumul_prob_obs <- pnorm(observed_stat)
  
  # by definition of p-value
  p <- 1 - cumul_prob_obs
  
  return(p)
  
} # end p_level

#=========================================
# critical_value(one_sample_z, alpha_one_tailed, sample_n)
#=========================================
#' @export
critical_value.OneSampleZ <- function(stat_procedure, alpha_one_tailed, sample_n)
{
  #print(stat_procedure)
  
  
  # To obtain your critical value, find the z at 1-alpha percentile in the Ho: n(0,1) distribution
  # e.g. if alpha = 0.05, you want the z that cuts off the lower 95% of the Ho distribution.
  
  z_crit_cumul_prob = 1 - alpha_one_tailed
  
  z_critical = qnorm(z_crit_cumul_prob)
  
  return(z_critical)
  
} # end critical_value

#=========================================
# noncentrality(sample_n, effect_size)
#=========================================
#' @export
noncentrality.OneSampleZ <- function(stat_procedure, sample_n, effect_size)
{
  # By definition
  noncen = sqrt(sample_n) * effect_size
  return(noncen)
  
} # end noncentrality

#=========================================
# effect_size_from_stat(stat_value, sample_n)
#=========================================
#' @export
effect_size_from_stat.OneSampleZ <- function(stat_procedure, stat_value, sample_n)
{
  # By definition
  effect_size = stat_value / sqrt(sample_n);
  return(effect_size)
  
} # end effect_size_from_stat

#=========================================
# expected_significant_effect - See StatProcedureOneSampleT.R for explanation of logic
#=========================================
#' @export
expected_significant_effect.OneSampleZ <- function(stat_procedure, alpha_one_tailed, sample_n, effect_size)
{
  
  # Fetch the appropriate critical value and noncentrality parameter as defined for this statistic
  z_critical <- critical_value(stat_procedure, alpha_one_tailed, sample_n)
  noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)
  
  # We will integrate over the portion of the noncentral sampling distribution greater than z-critical
  lower_bound <- z_critical
  upper_bound <- Inf
  
  # Required to make a proportional adjustment in the computation of expected effect size
  total_area_above_crit <- 1 - pnorm(z_critical - noncen_parameter)
  
  # The function to integrate over. The variable noncen_parameter must be
  # in scope and initialised at the time of the call to integral(args) below
  z_prob_times_value <- function(z)
  {
    dnorm(z - noncen_parameter) * z
  }
  
  # Integral produces a message about algorithm used. We don't want it, so we wrap the
  # call in capture.output
  capture.output(expected_value_numerator <- pracma::integral(z_prob_times_value, lower_bound, upper_bound))
  expected_z <- expected_value_numerator/total_area_above_crit
  
  expected_effect_size = effect_size_from_stat(stat_procedure, expected_z, sample_n)
  
  return(expected_effect_size)
  
} # end expected_significant_effect

#=========================================
# expected_significant_effect_bounded - See StatProcedureOneSampleT.R for explanation of logic
#=========================================
#' @export
expected_significant_effect_bounded.OneSampleZ <- function(stat_procedure, alpha_strong, alpha_weak, sample_n, effect_size)
{
  z_critical_strong <- critical_value(stat_procedure, alpha_strong, sample_n)
  z_critical_weak <- critical_value(stat_procedure, alpha_weak, sample_n)
  noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)
  
  lower_bound <- z_critical_weak
  upper_bound <- z_critical_strong
  
  total_area_above_strong <- 1 - pnorm(z_critical_strong - noncen_parameter)
  total_area_above_weak <- 1 - pnorm(z_critical_weak - noncen_parameter)
  total_area_between <- total_area_above_weak - total_area_above_strong
  
  z_prob_times_value <- function(z)
  {
    dnorm(z - noncen_parameter) * z
  }
  
  # Integral produces a message about algorithm used. We don't want it, so we wrap the
  # call in capture.output
  capture.output(expected_value_numerator <- pracma::integral(z_prob_times_value, lower_bound, upper_bound))
  expected_z <- expected_value_numerator/total_area_between
  
  expected_effect_size = effect_size_from_stat(stat_procedure, expected_z, sample_n)
  
  return(expected_effect_size)
  
} # end expected_significant_effect_bounded


#=========================================
# divide_total_sample
#=========================================
#' @export
divide_total_sample.OneSampleZ <- function(stat_procedure, n_per_segment)
{
  # In a one-sample z-test, all subjects are in a single group
  return(n_per_segment)
} # end divid_total_sample

# Class TrueEffect: Holds effect size (TrueEffect$effect_size) and probablity (TrueEffect$effect_size_probability),
# information needed to describe the true state of the world.
#
# An effect size of 0 implies that the null hypothesis (H0) is true.
#
# An experimental scenario may have one or more true effects, each with a probability of occurrence,
# (more generally called a "base rate"). The most typical situations are:
#
#    1) a single effect with probability 1, meaning we consider only those cases where the effect is present
#
#    2) a single effect with probability (base rate) < 1. The researcher expects the effect to be present in this
#    proportion of studies across an infinite number of studies in the area. The researcher expects there to be no effect (i.e. H0 to be true) in the remaining 1- base rate
#    proportion of studies.


##########################################
# Class properties
##########################################
#   effect_size:              Default = 0
#   effect_size_probability:  Default = 1
##########################################


##########################################
# Creation methods (following Wickham)
##########################################
#=========================================
# new_TrueEffect - low-level ctor.
#=========================================
new_TrueEffect <- function(effect_size = 0, effect_size_probability = 1)
{
  new_true_effect <- list(effect_size = effect_size,
                          effect_size_probability = effect_size_probability)
  
  
  attr(new_true_effect, "class") <- "TrueEffect"
  
  # Return instance
  return(new_true_effect)
  
} # new_TrueEffect
#=========================================


#=========================================
# validate_TrueEffect - error checking for raw TrueEffect instance
#=========================================
validate_TrueEffect <- function(true_effect)
{
  
  if ((true_effect$effect_size < 0) || (true_effect$effect_size > 1)){
    stop("Effect size value must be between 0 and 1, inclusive")
  }
  
  
  if ((true_effect$effect_size_probability < 0) || (true_effect$effect_size_probability > 1)){
    stop("Effect size probability values must be between 0 and 1, inclusive")
  }
  
} # end validate_TrueEffect
#=========================================


#=========================================
# TrueEffect - Public ctor
#=========================================
#-----------------
# Roxygen comments

#' TrueEffect
#'
#' Class TrueEffect holds effect size (TrueEffect$effect_size) and probablity (TrueEffect$effect_size_probability),
#' information needed to describe the true state of the world.
#'
#' An effect size of 0 implies that the null hypothesis (H0) is true.
#
#' An experimental scenario may have one or more true effects, each with a probability of occurrence,
#' (more generally called a "base rate"). The most typical situations are:
#'
#'    1) a single effect with probability 1, meaning we consider only those cases where the effect is present
#'
#'    2) a single effect with probability (base rate) < 1. The researcher expects the effect to be present in this
#'    proportion of studies across an infinite number of studies in the area. The researcher expects there to be no effect (i.e. H0 to be true) in the remaining 1- base rate
#'    proportion of studies.
#'
#' \code{TrueEffect}  is the public-facing TrueEffect constructor. It creates an unvalidated class instance with
#'  new_TrueEffect, and passes it to validate_TrueEffect. An exception is thrown if either parameter value
#'  is illegal (outside range {0..1}). It returns the instance if validation does not fail.
#'
#' @param effect_size A \code{numeric} that describes the effect size on a scale from 0 to 1 (e.g. Cohen's d or eta-squared)
#' @param effect_size_probability A \code{numeric} between 0 and 1; the base rate of this effect size. Set at 1 if only interested in cases where
#'   the effect is known to exist. Set at less than 1 to allow probability of Ho true.
#' @return Validated TrueEffect instance
#' @export
TrueEffect <- function(effect_size = 0, effect_size_probability = 1)
{
  
  # Create instance
  instance = new_TrueEffect(effect_size, effect_size_probability)
  
  # Validate instance. Execution halts here if things go wrong
  validate_TrueEffect(instance)
  
  # Return validated instance
  return(instance)
  
} # end TrueEffect ctor
#=========================================


##########################################
# Public methods
##########################################
#=========================================
# print.TrueEffect: System will dispatch to this method on print(instance)
#=========================================
#' @export
print.TrueEffect <- function(x, ...)
{
  data_values <- unclass(x) # This grabs the original list from the struct built when attr was called
  df <- data.frame(data_values)
  print(df)
}

# Class TwoSampleT descends from StatProcedureBase

# See class OneSampleT (StatProcedureOneSampleT.R) for detailed explanation of method logic

# All the computations ported from original source code (Miller, 2020)

# NB: All sample n's provided to these methods must be a vector of length 2. Passing
# a scalar or vector of length one throws an exception

#' @importFrom stats qt
#' @importFrom stats pt
#' @importFrom stats dt
#' @importFrom pracma integral
#' @importFrom utils capture.output

########################################
# Creation methods - following Wickham
########################################
new_TwoSampleT <- function(min_sample_n, display_abbrev, display_name)
{
  new_two_sample_t <- list(min_sample_n = min_sample_n,
                           display_abbrev = display_abbrev,
                           display_name = display_name)
  
  
  attr(new_two_sample_t, "class") <- c("TwoSampleT", "StatProcedureBase")
  
  return(new_two_sample_t)
}

#====================================================================
validate_TwoSampleT <- function(two_sample_t)
{
  # Do any validation here. Use stop(msg) to exit if things go wrong
}

#====================================================================
#-----------------
# Roxygen comments

#' TwoSampleT Constructor
#'
#' TwoSampleT is a child of StatProcedureBase, representing a two sample t-test
#'
#' \code{TwoSampleT}  is the public-facing TwoSampleT constructor. It creates and returns a class instance
#' with appropriately initialised values for all fields (min_sample_n, display_abbrev, and
#' display_name). It requires no parameters.
#'
#' @param min_sample_n Smallest n allowed for this test. Supplied in constructor.
#' @param display_abbrev Two-character abbreviation for this test. Developers see stat_procedure_factory.SegmentedHypTester
#' @param display_name String value useful for generating tidy outputs for this test
#'
#' @return Validated TwoSampleT instance
#'
#' @export
TwoSampleT <- function(min_sample_n = 4, display_abbrev = "2t", display_name = '2-sample t-test')
{
  instance <- new_TwoSampleT(min_sample_n, display_abbrev, display_name)
  validate_TwoSampleT(instance)
  return(instance)
}
#################################
# Instance methods
#################################
#=========================================
# compute_power(alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
#=========================================
#-----------------
# Roxygen comments

#' compute_power.TwoSampleT
#'
#' Implementation of compute_power method for class TwoSampleT
#'
#' \code{TwoSampleT} computes the statistical power of a described two sample
#' t-test as the probability density of t-critical in the t-distribution for a given effect size
#'
#' @param stat_procedure TwoSampleT instance for method dispatch
#' @param alpha_one_tailed Alpha level for one-tailed test
#' @param effect_size Presumed true effect size on Cohen's d scale
#' @param sample_n Group sizes. Must be a numeric vector of length 2
#' @param critical_value Optional to speed computation. If omitted, this value is computed in the method
#'
#' @return Computed statistical power (probability of rejecting H0 when the effect is present)
#'
#' @export
compute_power.TwoSampleT <- function(stat_procedure, alpha_one_tailed, effect_size, sample_n, critical_value = NULL)
{
  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")
  
  # degrees of freedom for a two sample t-test is n1 + n2 - 2
  df <- sample_n[1] + sample_n[2] - 2
  
  t_crit_cumul_prob <- 1 - alpha_one_tailed
  
  if (is.null(critical_value))
  {
    t_crit <- qt(t_crit_cumul_prob, df)
  }
  else
  {
    t_crit <- critical_value
  }
  
  noncentrality_parameter <- noncentrality(stat_procedure, sample_n, effect_size)
  
  power <- 1 - pt(t_crit, df, ncp = noncentrality_parameter)
  
  return(power)
} # end compute_power.TwoSampleT


#=========================================
# p_level(two_sample_t, observed_stat, sample_n)
#=========================================
#' @export
p_level.TwoSampleT <- function(stat_procedure, observed_stat, sample_n)
{
  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")
  
  # degrees of freedom for a two sample t-test is n1 + n2 - 2
  df <- sample_n[1] + sample_n[2] - 2
  
  # the observed sample exceeds this proportion of t when Ho is true
  cumul_prob_obs <- pt(observed_stat, df)
  
  # by definition of p-value
  p <- 1 - cumul_prob_obs
  
  return(p)
  
} # end p_level

#=========================================
# critical_value(two_sample_t, alpha_one_tailed, sample_n)
#=========================================
#' @export
critical_value.TwoSampleT <- function(stat_procedure, alpha_one_tailed, sample_n)
{
  #print(stat_procedure)
  
  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")
  
  # degrees of freedom for a two sample t-test is n1 + n2 - 2
  df <- sample_n[1] + sample_n[2] - 2
  
  t_crit_cumul_prob <- 1 - alpha_one_tailed
  
  t_critical <- qt(t_crit_cumul_prob, df)
  
  return(t_critical)
  
} # end critical_value

#=========================================
# noncentrality(sample_n, effect_size)
#=========================================
#' @export
noncentrality.TwoSampleT <- function(stat_procedure, sample_n, effect_size)
{
  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")
  
  
  # By definition
  noncen <- sqrt((sample_n[1] * sample_n[2]) / (sample_n[1] + sample_n[2])) * effect_size
  return(noncen)
  
} # end noncentrality

#=========================================
# effect_size_from_stat(stat_value, sample_n)
#=========================================
#' @export
effect_size_from_stat.TwoSampleT <- function(stat_procedure, stat_value, sample_n)
{
  
  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")
  
  
  # By definition
  effect_size <- stat_value * sqrt( (sample_n[1] + sample_n[2])/ (sample_n[1] * sample_n[2]))
  return(effect_size)
  
} # end effect_size_from_stat

#=========================================
# expected_significant_effect
#=========================================
#' @export
expected_significant_effect.TwoSampleT <- function(stat_procedure, alpha_one_tailed, sample_n, effect_size)
{
  
  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")
  
  t_critical <- critical_value(stat_procedure, alpha_one_tailed, sample_n)
  noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)
  
  lower_bound <- t_critical
  upper_bound <- Inf
  
  df <- sample_n[1] + sample_n[2] - 2
  total_area_above_crit <- 1 - pt(t_critical, df, noncen_parameter)
  
  
  # The function to integrate over. The variables df and noncen_parameter must be
  # in scope and initialised at the time of the call to integral(args) below
  t_prob_times_value <- function(t)
  {
    dt(t, df = df, ncp = noncen_parameter) * t
  }
  
  # Integral produces a message about algorithm used. We don't want it, so we wrap the
  # call in capture.output
  capture.output(expected_value_numerator <- pracma::integral(t_prob_times_value, lower_bound, upper_bound))
  expected_t <- expected_value_numerator/total_area_above_crit
  
  expected_effect_size <- effect_size_from_stat(stat_procedure, expected_t, sample_n)
  
  return(expected_effect_size)
} # end expected_significant_effect

#=========================================
#expected_significant_effect_bounded
#=========================================
#' @export
expected_significant_effect_bounded.TwoSampleT <- function(stat_procedure, alpha_strong, alpha_weak, sample_n, effect_size)
{
  # Check that sample_n is a vector of length two. If not, throw exception
  if (length(sample_n) < 2)
    stop("Sample n values for a two sample t-test must be a vector of length 2 -- one value for each group size. Sizes may be equal or unequal.")
  
  
  t_critical_strong <- critical_value(stat_procedure, alpha_strong, sample_n)
  t_critical_weak <- critical_value(stat_procedure, alpha_weak, sample_n)
  noncen_parameter <- noncentrality(stat_procedure, sample_n, effect_size)
  
  lower_bound <- t_critical_weak
  upper_bound <- t_critical_strong
  
  df <- sample_n[1] + sample_n[2] - 2
  total_area_above_strong <- 1 - pt(t_critical_strong, df, noncen_parameter)
  total_area_above_weak <- 1 - pt(t_critical_weak, df, noncen_parameter)
  total_area_between <- total_area_above_weak - total_area_above_strong
  
  
  # The function to integrate over. The variables df and noncen_parameter must be
  # in scope and initialised at the time of the call to integral(args) below
  t_prob_times_value <- function(t)
  {
    dt(t, df = df, ncp = noncen_parameter) * t
  }
  
  # Integral produces a message about algorithm used. We don't want it, so we wrap the
  # call in capture.output
  capture.output(expected_value_numerator <- pracma::integral(t_prob_times_value, lower_bound, upper_bound))
  expected_t <- expected_value_numerator/total_area_between
  
  expected_effect_size <- effect_size_from_stat(stat_procedure, expected_t, sample_n)
  
  return(expected_effect_size)
  
  
} # end expected_significant_effect_bounded

#=========================================
# divide_total_sample
#=========================================
#' @export
divide_total_sample.TwoSampleT <- function(stat_procedure, n_per_segment)
{
  
  # split total subject into equal groups (or as close as possible, if total n is odd)
  group_sizes <- numeric(2)
  group_sizes[1] <- floor(n_per_segment/2)
  group_sizes[2]  <- n_per_segment - group_sizes[1]
  
  return(group_sizes)
} # end divid_total_sample