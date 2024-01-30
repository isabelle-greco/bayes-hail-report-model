### Script to calculate Bayesian p-values, evaluation metrics for dichotomous 
### events, and scores for probabilistic events on model output.
### Draws from work in 230222_metrics.qmd
###
### Author: Isabelle Greco
### Last updated: 2023-07-31

# par arguments
args = commandArgs(trailingOnly = TRUE)

# assert correct number of args
if (length(args) != 4) {
  stop(paste("Must give library for R package installs, directory in which to",
	     "save the results of this code, directory in which to find the",
	     "model results, and the name of the model."),
       call. = FALSE)
}


# load/install relevant packages
# note putting tidyverse later stops MASS masking the dplyr::select function
for (package in c("foreach", "doParallel", "poibin", "tidyverse", "lubridate", "bigstatsr", "loo")) {
  if (!require(package, character.only = TRUE)) {
    # if not installed, install in the directory given by first argument
    install.packages(package, lib = args[1], 
                     repos = "https://cran.csiro.au", dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# define functions 
calculate_pvalues_by_grouping_col <- function(data, group, report, 
                                              ppc = "post_pred") {
  # Calculate p-values (lower and upper) by grouping by a given column.
  # 
  # Inputs:
  #   data    tibble
  #           assumed to have columns `group` with the grouping variable, `ppc#` 
  #           with the ppc draws, and `report` with the observations. Names user
  #           supplied.
  #   group   colunm of data
  #           contains the grouping value 
  #   ppc     str
  #           name of the columns with the posterior predictive 
  #           samples from which metrics are calculated. 
  #           default is "ppc"
  #   report  column of data
  #           contains the ground truth
  #   
  # Output:
  #   _       tibble
  #           colunms corresponding to the group and lower and upper tailed 
  #           p-values
  
  # turn given column name into a format r likes
  group_name <- enquo(group)
  report_col <- enquo(report)
  
  data %>%
    # group by the groups
    group_by(!!group_name) %>%
    # take sums within groups
    summarise(across(c(!!report_col, starts_with(ppc)), sum)) %>%
    # make rowwise tibble
    rowwise() %>%
    # calculate the pvalues
    mutate(pvalue_upper = mean(c_across(starts_with(ppc)) >= !!report_col),
           pvalue_lower = mean(c_across(starts_with(ppc)) <= !!report_col)) %>%
    # gets rid of ppc for display 
    select(-starts_with(ppc), -!!report_col)
}

calculate_binary_metrics_by_grouping_col <- function(data, group, report, 
                                                     ppc = "post_pred", 
                                                     summary_quantile = 0.05) {
  # Calculate binary metrics by grouping by a given column.
  # 
  # Inputs:
  #   data              tibble
  #                     assumed to have columns `group` with the grouping 
  #                     variable `post_pred` with the ppc draws, and `report`  
  #                     with the observations. exact names user supplied
  #   group             column of data
  #                     contains the grouping value 
  #   ppc               str
  #                     name of the columns with the posterior predictive 
  #                     samples from which metrics are calculated. 
  #                     default is "post_pred"
  #   report            column of data
  #                     contains the ground truth
  #   summary_quantile  double
  #                     quantile to use to summarise the data, default 5%
  #   
  # Output:
  #   _               tibble
  #                   colunms corresponding to the group ad then each of the
  #                   desired metrics (coded in the function)
  
  # turn given column name into a format r likes
  group_name <- enquo(group)
  report_col <- enquo(report)
  
  # create contingency table 
  contingency_table_with_group <- data %>%
    # group with provided groups
    group_by(!!group_name) %>%
    # calculate components for each group
    summarise(across(starts_with(ppc), .fns = list(
      # true positives
      a = function(x) sum(x == 1 & !!report_col == 1),
      # false negatives
      b = function(x) sum(x == 0 & !!report_col == 1),
      # false positives
      c = function(x) sum(x == 1 & !!report_col == 0),
      # true negatives 
      d = function(x) sum(x == 0 & !!report_col == 0)
    ))
    ) %>%
    # longer to get columns
    pivot_longer(-!!group_name)%>%
    # separating column names
    # won't work if _a, ..., _d appears in `ppc`
    mutate(metric = str_split_i(name, "_", -1)) %>% # last part after "_"
    # first part before _a etc.
    mutate(ppc_idx = str_split_i(name, paste("_", metric, sep = ""), 1)) %>% 
    # wider to get columns and groups
    pivot_wider(names_from = metric, values_from = value, 
                id_cols = c(ppc_idx, !!group_name))
  
  # metrics of interest as a function of the four components of the contingency 
  # table
  
  # EDI
  metric_extremal_dependence_index <- function(a, b, c, d) {
    # will return NaN if at least one of hit rate or false alarm rate is zero
    # positively oriented (i.e. maximised at one is best)
    # hit rate
    hit_rate <- a / (a + b)
    # false alarm rate
    far <- c / (c + d)
    # edi
    (log(far) - log(hit_rate)) / (log(far) + log(hit_rate))
  }
  # SEDI
  metric_symmetric_extremal_dependence_index <- function(a, b, c, d) {
    # will also return NaN if at least one of the hit rate or false alarm rate
    # is zero
    # also positively oriented (i.e. maximised at one is best)
    # hit rate
    hit_rate <- a / (a + b)
    # false alarm rate
    far <- c / (c + d)
    # sedi
    (log(far) - log(hit_rate) + log(1 - hit_rate) - log(1 - far)) / 
      (log(far) + log(hit_rate) + log(1 - hit_rate) + log(1 - far))
  }
  # Peirce Skill Score
  # positively oriented (i.e. maximised at one is best)
  metric_peirce_skill_score <- function(a, b, c, d){
    (a / (a + b)) - (c / (c + d))
  }
  # Accuracy
  # positively oriented (i.e. maximised at one is best)
  metric_accuracy <- function(a, b, c, d){
    (a + d) / (a + b + c + d)
  }
  # Critical success index
  # positively oriented (i.e. maximised at one is best)
  metric_critical_success_index <- function(a, b, c, d){
    a / (a + b + c)
  }
  
  # calculating metrics
  metrics_with_group <- contingency_table_with_group %>%
    mutate(metric_extremal_dependence_index = 
             metric_extremal_dependence_index(a, b, c, d),
           metric_symmetric_extremal_dependence_index = 
             metric_symmetric_extremal_dependence_index(a, b, c, d),
           metric_peirce_skill_score = metric_peirce_skill_score(a, b, c, d),
           metric_accuracy = metric_accuracy(a, b, c, d),
           metric_critical_success_index = 
             metric_critical_success_index(a, b, c, d)) %>%
    select(-c(a, b, c, d))
  
  # summarise the desired percentile 
  # note NaN results are carried through - if NA for any of the PPC then result is NA
  # will only effect EDI and SEDI
  metrics_with_group %>%
    group_by(!!group_name) %>%
    summarise(across(starts_with("metric"),
		     ~ tryCatch(quantile(.x, prob = summary_quantile), error = function(e) NA)))
}

calculate_prob_metrics_by_grouping_col <- function(data, group, report, 
                                                   prob = "prob_report",
                                                   summary_quantile = 0.05) {
  # Calculate probabilistic metrics by grouping by a given column.
  # 
  # Inputs:
  #   data              tibble
  #                     assumed to have columns `group` with the grouping 
  #                     variable `prob#` with the prob draws, and `report` with 
  #                     the observations. exact names user supplied
  #   group             colunm of data
  #                     contains the grouping value 
  #   prob              str
  #                     name of the columns with the probability of reporting
  #                     samples from which metrics are calculated. 
  #                     default is "prob_report"
  #   report            column of data
  #                     contains the ground truth
  #   summary_quantile  double
  #                     quantile to use to summarise the data, default 5%
  #   
  # Output:
  #   _               tibble
  #                   colunms corresponding to the group ad then each of the
  #                   desired metrics (coded in the function)
  
  # turn given column name into a format r likes
  group_name <- enquo(group)
  report_col <- enquo(report)
  
  # area under receiver operating curve
  metric_auc <- function(obs, pred){
    # positively oriented (i.e. maximised at one is best)
    # only want area (A) not everything else
    # will error if only successes or failures observed - we set to return NA 
    # in this case returns NA
    tryCatch(suppressWarnings(bigstatsr::AUC(pred, obs)), 
             error = function(error) NA)
  }
  # logarithmic score
  metric_log <- function(obs, pred){
    # positively oriented (i.e. maximised at zero is best)
    # log of the probability assigned to the event occurred
    # the multiplication and summation deals with the yes/no logic 
    mean(log(obs * pred + (1 - obs) * (1 - pred)))
  }
  # discrimination statistic (DO)
  metric_discrimination <- function(obs, pred){
    # positively oriented (maximising is best, normalised to [0, 1])
    mu_x <- mean(obs)
    mu_f_given_x0 <- mean(pred[obs == 0])
    mu_f_given_x1 <- mean(pred[obs == 1])
    mu_f <- mean(pred)
    dis <- (1 - mu_x) * ((mu_f_given_x0 - mu_f)^2) + 
      mu_x * ((mu_f_given_x1 - mu_f) ^ 2)
    dis / (mu_x * (1 - mu_x))
  }
  # Poisson-Binomial calibration check
  metric_pb_calibration <- function(obs, pred){
    # positively oriented (i.e. maximised is best)
    # technically bounded at 1 but more realistically bounded just over 0.5
    #
    # supposing forecast probs are true, probability that at most the number of 
    # events observed are observed or at least the number of events observed
    # are observed.
    #
    # using the refined normal approximation (RNA) method is a fast and good
    # approach to avoid segfaults with larger data
    prob_at_most_obs_events <- poibin::ppoibin(sum(obs), pp = pred, method = "RNA")
    prob_at_least_obs_events <- 1 - poibin::ppoibin(sum(obs) - 1, pp = pred, method = "RNA")
    # takes the smallest tail probability    
    min(prob_at_most_obs_events, prob_at_least_obs_events)
  }
  
  ### the bulk of the calculations and pivoting ### 
  # passing in data
  data %>%
    # selecting relevant columns
    select(!!group_name, !!report_col, starts_with(prob)) %>% 
    # grouping by the group
    group_by(!!group_name) %>%
    # computing the various metrics and naming columns appropriately 
    summarise(across(starts_with(prob), 
                     list(metric_auc = ~ metric_auc(!!report_col, .x),
                          metric_log = ~ metric_log(!!report_col, .x),
                          metric_discrimination = 
                            ~ metric_discrimination(!!report_col, .x),
                          metric_pb_calibration = 
                            ~ metric_pb_calibration(!!report_col, .x)), 
                     .names = "{.fn}__{.col}")) %>%
    # pivoting longer
    pivot_longer(-!!group_name) %>%
    # separating column names into metric and draw number
    mutate(metric = str_split_i(name, "__", 1),  
           draw_number = str_split_i(name, "__", 2)) %>%
    # dropping now redundant name column 
    select(-name) %>%
    # wider to get columns for metrics by each group and ppc
    pivot_wider(names_from = metric, values_from = value) %>%
    # grouping by each group to...
    group_by(!!group_name) %>%
    # ...calculate the quantiles 
    # carries NAs through - if any iteration gives NA then whole result is NA
    summarise(across(starts_with("metric"), 
		     ~ tryCatch(quantile(.x, prob = summary_quantile), error = function(e) NA)))
}

aggregate_all_metrics <- function(data, group, report, folder_name) {
  # Aggregate the worst metrics calculated for a given grouping column and 
  # ground truth
  # 
  # Inputs:
  #   data        tibble
  #               assumed to have columns `group` with the grouping variable 
  #               `prob#` with the prob draws, `ppc#` with the ppc draws, and 
  #               `report` with the observations. exact names user supplied 
  #   group       column of data
  #               contains the grouping value 
  #   report      column of data
  #               contains the ground truth
  #   folder_name str
  #               folder in which to save the group-level results prior to  
  #               aggregating across levels
  #   
  # Output:
  #   _       tibble with one row
  #           columns corresponding to each of the desired metrics (coded in 
  #           the function)
  
  # turn given column name into a format r likes
  group_name <- enquo(group)
  report_col <- enquo(report)
  
  # calculating each of the metrics required
  pvalues <- data %>%
    calculate_pvalues_by_grouping_col(!!group_name, !!report_col) %>%
    # selecting smaller of the two tails
    mutate(pvalue = min(pvalue_upper, pvalue_lower)) %>%
    select(-pvalue_upper, -pvalue_lower)
  binary <- data %>%
    calculate_binary_metrics_by_grouping_col(!!group_name, !!report_col)
  prob <- data %>%
    calculate_prob_metrics_by_grouping_col(!!group_name, !!report_col)
  
  # calculating size of each level for reference
  level_size <- data %>% 
    group_by(!!group_name) %>% 
    summarise(num_obs_in_level = n(),
              num_report_in_level = sum(as.numeric(as.character(!!report_col))))
  
  all_metrics <- pvalues %>%
    inner_join(binary, by = join_by(!!group_name)) %>%
    inner_join(prob, by = join_by(!!group_name)) %>%
    inner_join(level_size, by = join_by(!!group_name))
  
  # saving the results presented by level for manual analysis if desired
  file_name <- paste(folder_name, "/metrics_", as_label(group_name), ".csv", 
                     sep = "")
  all_metrics %>%
    write_csv(file_name)
  
  # calculating worst across the groups
  all_metrics %>%
    # don't need the number in level any more
    select(-num_obs_in_level, -num_report_in_level) %>%
    # drop the grouping column
    select(-!!group_name) %>%
    # undo the grouping
    ungroup() %>%
    # get the minimum (the worst performing level) 
    # note that if ANY level is NA, this will result in an NA
    summarise(across(everything(), min))
}

cut_name_in_original_scale <- function(x, fun, round_digits = 2) {
  # Convenience function which will, using the supplied inverse function
  # take names generated by `cut` and relable to the correct scale.
  #
  # Inputs:
  #   x             str
  #                 the string in the form (BLAH, BLAH] to convert
  #                 not sensitive to ordering of input bracket but will 
  #                 always output in form (BLAH, BLAH]
  #   fun           function
  #                 the supplied inverse function (e.g. exp will 
  #                 undo the log transofrmation to `cut`)
  #   round_digits  str
  #                 number of digits to which the result should be
  # 		    rounded. default two

  # split the string using ( , ] characters	
  split_label <- str_split_1(x, "\\(|\\]|\\,")
  # put it on the new scale by first dropping empty strings...
  new_scale <- split_label[nzchar(split_label)] %>%
    # ...then making numeric...
    as.numeric %>%
    # ...then applying the function...
    fun %>%
    # ..and rounding
    round(digits = round_digits)
  # pasting together into correct format
  paste("(", new_scale[1], ",", new_scale[2], "]", sep = "")
}

# directory in which to save results
folder_name <- args[2]

# reading in data as suppplied
ppc <- readRDS(paste(args[3], "/", args[4], "_posterior_predictive_samples_cv.rds", sep = ""))
prob <- readRDS(paste(args[3], "/", args[4], "_probability_report_cv.rds", sep = ""))
join_cols <- prob %>% 
  select(-starts_with("prob_report")) %>% 
  colnames()
test_data <- inner_join(ppc, prob, by = join_cols) %>%
  mutate(report = as.integer(as.character(report)))

# create the groups used for evaluation
test_data_grouped <- test_data %>%
  mutate(
    # spatial grouping
    group_density_individual = as_factor(pop_dens),
    group_density_binned = cut(log(pop_dens), breaks = 4), # four bins in log space - ensures enough reports in each
    group_density_urban_rural = as_factor(case_when(pop_dens < 100 ~ "rural",
                                                    TRUE ~ "urban")),
    group_density_latitude = as_factor(y_bins),
    group_density_longitude = as_factor(x_bins),
    # mesh grouping
    group_mesh_fine_bins = cut(log(mesh + 1), breaks = 5),
    group_mesh_threshold1 = as_factor(case_when(mesh < 22 ~ "mesh < 22",
                                                TRUE ~ "mesh >= 22")),
    group_mesh_threshold2 = as_factor(case_when(mesh < 32 ~ "mesh < 32",
                                                TRUE ~ "mesh >= 32")),
    group_mesh_threshold3 = as_factor(case_when(mesh < 44 ~ "mesh < 44",
                                                TRUE ~ "mesh >= 44")),
    group_everything = as_factor(1)
    ) %>%
  # re-using existing information for temporal grouping
  rename(group_time_dow = dow,
	 group_time_month = month,
	 group_time_storm_season = storm_season,
	 group_time_weekend = weekend,
	 group_time_time_of_day = time_of_day) %>%
  # relabelling transformed factors
  mutate(group_density_binned = fct_relabel(
            group_density_binned, .fun = ~ map_chr(., cut_name_in_original_scale, fun = exp)
          ),
	 group_mesh_fine_bins = fct_relabel(
            group_mesh_fine_bins, .fun = ~ map_chr(., cut_name_in_original_scale, fun = function(x) exp(x) - 1)
         )
    )

# specifying grouping columns to group over
loop_over_cols <- colnames(test_data_grouped)[grep("group", colnames(test_data_grouped))]

# setting up a cluster with one worker for each group
cl <- makeCluster(length(loop_over_cols))
registerDoParallel(cl, cores = parallel::detectCores() - 4) # assumes at least 4 cores

# calculating metrics for multiple groups - save results for each group
# level within the function 
# executed in parallel
leaf_metrics <- foreach(col = loop_over_cols, .combine = rbind, .packages=c("dplyr", "tidyr", "stringr", "readr")) %dopar% {
  test_data_grouped %>%
    aggregate_all_metrics(!!sym(col), report, folder_name = folder_name)
} %>%
  add_column(name = loop_over_cols, .before = "pvalue") %>%
  # removing the group_density_individual, only wanted the saved files
  filter(name != "group_density_individual") %>%
  # same for group_mesh_fine_bins
  filter(name != "group_mesh_fine_bins") %>%
  # and for the month
  filter(name != "group_time_month")

# saving leaf metrics prior to aggregation
leaf_metrics %>%
  write_csv(paste(folder_name, "metrics_leaf.csv", sep = "/"))

# TODO: tune these numbers
# aggregating up leaves
agg_metrics <- leaf_metrics %>%
  # get the larger group name (e.g. MESH, density, time, everything)
  mutate(group_name = str_split_i(name, "_", 2)) %>%
  # get the smaller group name (leaf name)
  mutate(leaf_name = str_split_i(name, paste(group_name, "_", sep = ""), 2)) %>%
  # drop the now redundant full name 
  select(-name) %>%
  # grouping leaves by overall groups
  group_by(group_name) %>%
  # take the average score of the leaves within the groups - assumption of 
  # equally important performance between the different catagorisations
  # note we now drop the nas
  summarise(across(where(is.numeric), \(x) tryCatch(mean(x), error = function(e) NA)), 
	    .groups = "drop") %>%
  # weight the overall groups 
  summarise(across(where(is.numeric), 
                   # weighted down as will be examined separately
                   ~ (1 / 7) * .x[group_name == "everything"] 
                   + (2 / 7) * .x[group_name == "mesh"] 
                   + (2 / 7) * .x[group_name == "density"]
                   + (2 / 7) * .x[group_name == "time"]
  )) %>%
  # then aggregate the different metrics, placed at the front
  # the metrics employed are all in range [0,  1] though we note that really
  # anything past 0.5 we're pretty happy with
  # hence consider 'realistic' upper bound to be around 2/3
  # from experiements with 'idealised' models (i.e. drawing from the correct
  # probabilities and perturbations of the probabilites) believe a realistic
  # 'really good' value is 0.45 ish
  # the difference is due in part to the nature of the problem (lots of 
  # low probability events, we're not trying to build a classifier) and also
  # the fact that we're looking at the 5% quantile
  mutate(agg_metric =
           0.2 * metric_peirce_skill_score +                           		    # deterministic weighted lower
           0.8 * ((2/3) * min(metric_pb_calibration / 0.1, 1)+ (1/3) * metric_auc), # want discrimination subject to calibration hence the weighting
                                                                       		    # note that the discrimination value also considers discrimination but take AUC as more holistic
                                                                       		    # log score also useful but encapsulated in separate likelihood investigation
	 .before = everything()
  )

# save aggregation
agg_metrics %>% 
  write_csv(paste(folder_name, "metrics_aggregated.csv", sep = "/"))

# exracting pointwise log likelihood
full_log_lik <- readRDS(paste(args[3], "/", args[4], "_log_likelihood_full.rds", sep = ""))

# provide relative effective sample sizes for better performance
r_eff <- relative_eff(exp(full_log_lik), cores = parallel::detectCores() - 4)

# calculate loo
loo_obj <- loo(full_log_lik, r_eff = r_eff, cores = parallel::detectCores() - 4)

# saving 
loo_obj %>%
  saveRDS(paste(folder_name, "loo_object.rds", sep = "/"))
