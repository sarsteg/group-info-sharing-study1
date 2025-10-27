
# Common Functions


library(psych)
library(dplyr)
library(tidyr)
library(haven)
library(descr)
library(ggplot2)
library(readr)

#---------------------------------------------------------------------------




# Example
# corr_table(filtered_data, )


# Requires: psych
corr_table <- function(
    data,
    vars,
    labels = NULL,                   # either same-length vector or a *named* vector mapping var -> label
    method = c("spearman","pearson","kendall"),
    adjust = "none",                 # psych::corr.test adjust arg: "none","holm","bonferroni", etc.
    digits = 2,
    star_cutoffs = c("***"=0.001,"**"=0.01,"*"=0.05),
    blank_upper = TRUE,
    blank_diag  = TRUE,
    use = "pairwise",                # passed to psych::corr.test (pairwise or complete)
    coerce_factors = FALSE,          # set TRUE to coerce factors to numeric (levels)
    return = c("table","both")       # "both" also returns r, p, and the psych object
) {
  method <- match.arg(method)
  return <- match.arg(return)
  
  # subset & optional factor coercion
  df <- data[, vars, drop = FALSE]
  if (coerce_factors) {
    df[] <- lapply(df, function(x) if (is.factor(x)) as.numeric(x) else x)
  }
  
  # apply labels
  if (!is.null(labels)) {
    if (!is.null(names(labels))) {
      # named mapping (var -> label)
      match_idx <- match(colnames(df), names(labels))
      new_names <- ifelse(!is.na(match_idx), labels[match_idx], colnames(df))
      colnames(df) <- new_names
    } else {
      # positional labels
      if (length(labels) != ncol(df)) stop("labels must match length of vars or be a named vector.")
      colnames(df) <- labels
    }
  }
  
  # correlation (psych gives r and p nicely)
  ct <- psych::corr.test(df, use = use, method = method, adjust = adjust)
  
  r <- ct$r
  p <- ct$p
  
  # build stars matrix
  stars <- matrix("", nrow = nrow(p), ncol = ncol(p), dimnames = dimnames(p))
  for (lab in names(star_cutoffs)) {
    stars[p < star_cutoffs[lab]] <- lab
  }
  
  # formatted r with stars
  formatted <- matrix(
    paste0(format(round(r, digits), nsmall = digits), stars),
    nrow = nrow(r), ncol = ncol(r),
    dimnames = dimnames(r)
  )
  
  # blank upper/diag if requested
  if (blank_upper) {
    formatted[upper.tri(formatted, diag = blank_diag)] <- ""
  }
  
  out_tbl <- as.data.frame(formatted, stringsAsFactors = FALSE)
  
  if (return == "table") {
    return(out_tbl)
  } else {
    return(list(
      table = out_tbl,
      r = r,
      p = p,
      psych = ct
    ))
  }
}






#--------------------------------------------
  
  
  
  
  
  
# # corr_spear: example usage ----
#   
# vars1 <- c(
#     "structure_number01","motivation_number01",
#     "public_more_shared","public_less_shared",
#     "private_more_shared","private_less_shared",
#     "Manipulation_Mean","Q14.1","sex","SES_mean",
#     "NeedCog_Mean","DirtyDozen_Mean","SVO_angle"
#   )
# 
# labels1 <- c(
#   "Structure","Motivation",
#   "Common High","Common Low",
#   "Unique High","Unique Low",
#   "StateCompCoop","Age","Sex","SES",
#   "NFC","Dirty Dozen","SVO"
# )
# 
# corr_spear <- corr_table(
#   data = filtered_data,
#   vars = vars1,
#   labels = labels1,
#   method = "spearman"
# )
# 
# corr_spear  # your formatted lower-tri table












#---------------------------------------------------------------------------
  
  


# Create table for descriptives with stats listed below
summarise_numeric <- function(.data, vars) {
  present <- intersect(vars, names(.data))
  if (length(present) == 0) {
    stop("None of the requested variables are in `filtered_data`.")
  }
  .data |>
    summarise(across(all_of(present),
                     list(
                       n    = ~sum(!is.na(.)),
                       mean = ~round(mean(., na.rm = TRUE), 2),
                       sd   = ~round(sd(., na.rm = TRUE), 2),
                       min  = ~round(min(., na.rm = TRUE), 2),
                       max  = ~round(max(., na.rm = TRUE), 2)
                     ),
                     .names = "{.col}_{.fn}")) |>
    pivot_longer(everything(),
                 names_to = c("variable", ".value"),
                 names_pattern = "(.*)_(n|mean|sd|min|max)") |>
    arrange(variable)
}





#---------------------------------------------------------------------------
  
  
  
  
## --- 3) 2×2 ANOVA helper
run_two_way_anova <- function(var) {
  f <- as.formula(paste(var, "~ structure * motivation"))
  fit <- aov(f, data = filtered_data)
  
  list(
    dv          = var,
    descriptives= desc_by_cell(var),
    anova_table = car::Anova(fit, type = 2) %>% broom::tidy(),
    eta_sq      = effectsize::eta_squared(fit, partial = TRUE) %>% as_tibble(),
    shapiro     = shapiro.test(residuals(fit)),
    levene      = leveneTest(f, data = filtered_data)
  )
}





#---------------------------------------------------------------------------






# Cook's distance at cluster level for glmgee/geeglm objects
# Flags clusters above cutoff (4/n or custom)
# Pipe in a glmgee/geeglm object

# Example use:

# # geepack::geeglm fit
# No, does not work for this package

# # glmtoolbox::glmgee 
# res2 <- cooks_by_cluster(glmgee_fit)
# res2$flagged

library(glmtoolbox)

cooks_by_cluster_glmgee <- function(fit) {
  # Cook's distance at cluster level (glmtoolbox handles this)
  cd_mat <- stats::cooks.distance(
    fit,
    method  = "full",
    level   = "clusters",
    plot.it = FALSE
  )
  
  # cluster IDs are stored directly in glmgee fits
  cluster_ids <- unique(fit$id)
  
  # results table
  cooks_tbl <- data.frame(
    ResponseId = cluster_ids,
    CooksD     = as.numeric(cd_mat[, 1])
  )
  
  # cutoff rule: 4/N
  cutoff <- 4 / nrow(cooks_tbl)
  cooks_tbl$Flag <- cooks_tbl$CooksD > cutoff
  
  cooks_tbl
}






#---------------------------------------------------------------------------




# desc_tbl <- map_dfr(controls, desc_by_cell)
# desc_tbl


desc_by_cell <- function(var) {
  # Check if var is quoted
  if (!is.character(var)) {
    stop("Variable name must be in quotes.")
  }
  
  filtered_data %>%
    group_by(structure, motivation) %>%
    summarise(
      n  = sum(!is.na(.data[[var]])),
      M  = mean(.data[[var]], na.rm = TRUE),
      SD = sd(.data[[var]],   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(variable = var) %>%
    relocate(variable)
}




#---------------------------------------------------------------------------






desc_by_cell_input <- function(data, group_vars, var) {
  # data: your dataset (e.g., dat_full)
  # group_vars: vector of column names to group by (quoted or unquoted)
  # var: name of the variable to summarize (quoted or unquoted)
  
  data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n  = sum(!is.na(.data[[deparse(substitute(var))]])),
      M  = mean(.data[[deparse(substitute(var))]], na.rm = TRUE),
      SD = sd(.data[[deparse(substitute(var))]],   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(variable = deparse(substitute(var))) %>%
    relocate(variable)
}


# Example of use: 
# desc_by_cell_input(dat_full, group_vars = c("structure", "motivation"), var = "Count")




#---------------------------------------------------------------------------





library(qgraph)
library(psych)

# install.packages(c("psych","qgraph"), dependencies = TRUE)



# Example code

# # Example: all 0/1 items
# vars <- c("REASONS_1","REASONS_2","REASONS_3","REASONS_4","REASONS_5")
# 
# # 1) Best-practice for dichotomous: tetrachoric (no p by default; clean APA table)
# tbl <- apa_corr_matrix(data = your_df, vars = vars, type = "tetrachoric", digits = 2)
# tbl
# attr(tbl, "note")
# 
# # 2) If you need stars and can afford a (slow) bootstrap:
# tbl_b <- apa_corr_matrix(data = your_df, vars = vars, type = "tetrachoric",
#                          bootstrap_p = TRUE, B = 500, digits = 2)
# tbl_b$table
# tbl_b$note
# 
# # 3) If you prefer phi with classic p-values/stars:
# tbl_phi <- apa_corr_matrix(data = your_df, vars = vars, type = "phi", digits = 2)
# tbl_phi
# attr(tbl_phi, "note")
# 
# # 4) Auto mode: tetrachoric if all 0/1, otherwise Pearson
# tbl_auto <- apa_corr_matrix(data = your_df, vars = vars, type = "auto", digits = 2)





apa_corr_matrix <- function(
    data,
    vars,
    labels = NULL,                     # optional vector; positional or named mapping var -> label
    type = c("auto","tetrachoric","phi","pearson","spearman","kendall"),
    digits = 2,
    stars = TRUE,
    star_cutoffs = c("***"=0.001,"**"=0.01,"*"=0.05),
    blank_upper = TRUE,
    blank_diag  = TRUE,
    use   = "pairwise",                # for psych::corr.test
    adjust = "none",                   # p adjustment for non-tetra methods
    coerce_factors = FALSE,
    bootstrap_p = FALSE,               # only used for tetrachoric
    B = 500,                           # bootstrap reps
    seed = NULL,
    return = c("table","both")         # "both" also returns r, p, and notes
){
  type <- match.arg(type)
  return <- match.arg(return)
  
  df <- data[, vars, drop = FALSE]
  
  # helper: 0/1 binary?
  is_binary01 <- function(x) is.numeric(x) && all(na.omit(x) %in% c(0,1))
  
  # optional factor -> numeric
  if (coerce_factors) {
    df[] <- lapply(df, function(x) if (is.factor(x)) as.numeric(x) else x)
  }
  
  # resolve labels
  var_names <- colnames(df)
  if (!is.null(labels)) {
    if (!is.null(names(labels))) {
      match_idx <- match(var_names, names(labels))
      var_names <- ifelse(!is.na(match_idx), labels[match_idx], var_names)
    } else {
      if (length(labels) != ncol(df)) stop("labels must match length of vars or be a named vector.")
      var_names <- labels
    }
  }
  
  # auto type
  if (type == "auto") {
    type <- if (all(vapply(df, is_binary01, logical(1)))) "tetrachoric" else "pearson"
  }
  
  note <- NULL
  r <- p <- NULL
  
  if (type %in% c("pearson","spearman","kendall")) {
    ct <- psych::corr.test(df, use = use, method = type, adjust = adjust)
    r <- ct$r; p <- ct$p
  } else if (type == "phi") {
    if (!all(vapply(df, is_binary01, logical(1)))) {
      stop("type='phi' requires all selected variables to be numeric 0/1.")
    }
    ct <- psych::corr.test(df, use = use, method = "pearson", adjust = adjust)
    r <- ct$r; p <- ct$p
    note <- "Note. Phi correlations (Pearson computed on 0/1 items)."
  } else if (type == "tetrachoric") {
    if (!all(vapply(df, is_binary01, logical(1)))) {
      stop("type='tetrachoric' requires all selected variables to be numeric 0/1.")
    }
    tc <- psych::tetrachoric(df)
    r <- tc$rho
    # default: no p-values for tetrachoric
    p <- matrix(NA_real_, nrow = nrow(r), ncol = ncol(r), dimnames = dimnames(r))
    note <- "Note. Tetrachoric correlations among 0/1 items. p-values not computed."
    
    # optional bootstrap to approximate significance (slow)
    if (bootstrap_p) {
      if (!is.null(seed)) set.seed(seed)
      n <- nrow(df); k <- ncol(df)
      # store bootstrapped r for each pair
      r_boot <- array(NA_real_, dim = c(k, k, B))
      for (b in seq_len(B)) {
        idx <- sample.int(n, n, replace = TRUE)
        rb <- try(psych::tetrachoric(df[idx, , drop = FALSE])$rho, silent = TRUE)
        if (!inherits(rb, "try-error"))
          r_boot[,,b] <- rb
      }
      # compute a p-like measure: 2*min(prop>0, prop<0)
      p_est <- matrix(NA_real_, k, k, dimnames = dimnames(r))
      for (i in 1:k) for (j in 1:k) if (i != j) {
        samp <- r_boot[i,j,]
        samp <- samp[is.finite(samp)]
        if (length(samp) > 10) {
          prop_pos <- mean(samp > 0, na.rm = TRUE)
          prop_neg <- mean(samp < 0, na.rm = TRUE)
          p_est[i,j] <- 2 * min(prop_pos, prop_neg)
        }
      }
      p <- p_est
      note <- paste0(
        "Note. Tetrachoric correlations among 0/1 items. p-values approximated via bootstrap (B=",
        B, ")."
      )
    }
  }
  
  # round & stars
  stars_mat <- matrix("", nrow = nrow(r), ncol = ncol(r), dimnames = dimnames(r))
  if (stars && !all(is.na(p))) {
    for (lab in names(star_cutoffs)) stars_mat[p < star_cutoffs[lab]] <- lab
  }
  
  fmt <- function(x) format(round(x, digits), nsmall = digits)
  formatted <- matrix(
    paste0(fmt(r), ifelse(is.na(p), "", stars_mat)),
    nrow = nrow(r), ncol = ncol(r),
    dimnames = list(var_names, var_names)
  )
  
  # blanking rules
  if (blank_upper) {
    formatted[upper.tri(formatted, diag = blank_diag)] <- ""
  }
  
  out_tbl <- as.data.frame(formatted, stringsAsFactors = FALSE)
  attr(out_tbl, "note") <- note
  attr(out_tbl, "type") <- type
  
  if (return == "table") {
    return(out_tbl)
  } else {
    return(list(
      table = out_tbl,
      r = r,
      p = p,
      note = note,
      type = type
    ))
  }
}










#---------------------------------------------------------------------------


















#---------------------------------------------------------------------------
