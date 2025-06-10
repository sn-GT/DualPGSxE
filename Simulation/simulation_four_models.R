library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(2025)
Nsim         <- 1e5
M            <- 5000
M_causal     <- 2500
hsq          <- 0.5
desired_prev <- 0.12
scale_factor <- 1
eps          <- 1e-4   # for clamping

# Simulated genotypes, betas, PRS and environment
maf    <- runif(M, 0.1, 0.9)
X      <- sapply(maf, function(p) rbinom(Nsim, 2, p))
X_s    <- scale(X)
p_c    <- sample(M, M_causal)
b_c    <- rnorm(M_causal, 0, sqrt(hsq / M_causal))
prs_base     <- X_s[, p_c] %*% b_c
PRS    <- as.vector(scale(prs_base) * sqrt(hsq))
E_base <- rbinom(Nsim, 1, 0.44)  

# Prevalence vs percentile PGS (exposure specific bins)
prevperc_groupPRSbins <- function(df, n_bins = 100) {
  df %>%
    group_by(GROUP) %>%
    mutate(PGS = as.integer(ntile(invGRS, n_bins))) %>%
    group_by(GROUP, PGS) %>%
    summarise(
      Prev    = 100 * mean(BCA),
      meanPRS = mean(invGRS),
      n       = n(),
      .groups = "drop"
    ) %>%
    arrange(GROUP, PGS)
}

# Liability threshold disease risk modeling on prevalence vs percentile PGS relationships - additive (null) model 
add_comp <- function(PP, overallprev) {
  PP2 <- PP %>%
    mutate(
      Prev_prop = pmin(pmax(Prev/100, eps), 1 - eps),  
      ENV       = if_else(GROUP == first(GROUP[which.max(Prev_prop)]), 1, 0)
    )
  t0  <- qnorm(1 - overallprev)
  PP2 <- PP2 %>%
    mutate(
      t_PGS = qnorm(1 - Prev_prop),
      ui    = t0 - t_PGS
    )
  fit    <- lm(ui ~ meanPRS + ENV, data = PP2)
  co     <- coef(fit)
  sigma0 <- summary(fit)$sigma
  PP2 %>%
    mutate(
      Estimated_ui         = co[1] + co[2]*meanPRS + co[3]*ENV,
      ExpectedPrev_nostoch = 1 - pnorm(t0, mean = Estimated_ui, sd = 1)
    )
}

# Liability threshold disease risk modeling on prevalence vs percentile PGS relationships - interaction model 
int_comp <- function(PP, overallprev) {
  PP2 <- PP %>%
    mutate(
      Prev_prop = pmin(pmax(Prev/100, eps), 1 - eps),
      ENV       = if_else(GROUP == first(GROUP[which.max(Prev_prop)]), 1, 0)
    )
  t0  <- qnorm(1 - overallprev)
  PP2 <- PP2 %>%
    mutate(
      t_PGS = qnorm(1 - Prev_prop),
      ui    = t0 - t_PGS
    )
  fit  <- lm(ui ~ meanPRS + ENV + meanPRS * ENV, data = PP2)
  co   <- coef(fit)
  PP2 %>%
    mutate(
      Estimated_ui         = co[1] + co[2]*meanPRS + co[3]*ENV + co[4]*meanPRS*ENV,
      ExpectedPrev_nostoch = 1 - pnorm(t0, mean = Estimated_ui, sd = 1)
    )
}


# Simulate disease case control status based on Y, E and PRS
# Then compute prevalence vs percentile PRS
generate_PP <- function(y_vec, E_vec, custom_PRS = NULL) {
  p_E1 <- mean(E_vec == 1)
  p_E0 <- 1 - p_E1
  
  # Prevalence of high-risk E=1 is set to 14%
  target_prev_E1 <- 0.14
  target_prev_E0 <- (desired_prev - p_E1 * target_prev_E1) / p_E0
  
  lp0 <- uniroot(function(l) mean(plogis(y_vec[E_vec == 0] + l)) - target_prev_E0,
                 lower = -15, upper = 15)$root
  lp1 <- uniroot(function(l) mean(plogis(y_vec[E_vec == 1] + l)) - target_prev_E1,
                 lower = -15, upper = 15)$root
  
  lp_vec <- ifelse(E_vec == 1, lp1, lp0)
  D <- rbinom(length(y_vec), 1, plogis(y_vec + lp_vec))
  
  PRS_use <- if (is.null(custom_PRS)) PRS else custom_PRS
  
  PP <- prevperc_groupPRSbins(data.frame(invGRS = PRS_use, BCA = D, GROUP = as.character(E_vec)))
  list(PP = PP, D = D)
}



# Four models - liability scores Y and PGS, E
scenarios <- list(
  # Null (Additive model)	
  Null = list(y = PRS + 0.5*E_base + rnorm(Nsim), E = E_base),
  # Alternate (Interaction model)
  Alt  = list(y = PRS + 0.5*E_base + 0.3*(PRS*E_base) + rnorm(Nsim), E = E_base),
  # Heterscedasticity: residual variance vary for E=0 and E=1
  Het  = list(y = PRS + 0.5*E_base + rnorm(Nsim, 0, ifelse(E_base==1,1,1.5)), E = E_base),
  # GE correlation
  GEC = {
    # Simulate E independently 
    E2 <- rbinom(Nsim, 1, 0.44)
    # Pull in the *same* global PRS you already computed
    PRS_raw <- PRS  
    # Shift PGS for E=1 by delta based on target correlation (rho)
    rho_target <- 0.5
    pE        <- mean(E2)
    sigma_prs <- sd(PRS_raw)
    delta     <- rho_target * sigma_prs / sqrt(pE * (1 - pE))
    
    # Induce a group‐mean shift for E=1
    # Choose delta so that mean(PRS | E=1) - mean(PRS | E=0) = desired offset
    PRS_shifted <- PRS_raw + delta * (E2 - mean(E2))
    
    #  Rescale back to exactly the original variance of PGS
    PRS_gec <- PRS_shifted * (sigma_prs / sd(PRS_shifted))
    # GE Model 
    y_gec <- PRS_gec + 0.5 * E2 + rnorm(Nsim)
 
    list(
      y   = y_gec,
      E   = E2,
      PRS = PRS_gec
    )
  }
  
  
  
)
#########################################################################################################################
# Single iteration 

# Compute metric of liabilit threshold disease risk modeling
results <- lapply(names(scenarios), function(name) {
  scen        <- scenarios[[name]]
  tmp <- if (name == "GEC") {
    generate_PP(scen$y, scen$E, custom_PRS = scen$PRS)
  } else {
    generate_PP(scen$y, scen$E)
  }
  
  PP_obs      <- tmp$PP
  overallprev <- mean(tmp$D)  
  print(overallprev)
  
  PP_add <- add_comp(PP_obs, overallprev)
  PP_add$ExpectedPrev_nostoch <- PP_add$ExpectedPrev_nostoch*100
  PP_int <- int_comp(PP_obs, overallprev)
  PP_int$ExpectedPrev_nostoch <- PP_int$ExpectedPrev_nostoch*100
  
  # Delta observed: deviations between H and L exposures at extremes of PGS
  r1   <- mean(PP_obs$Prev[PP_obs$PGS>98 & PP_obs$GROUP=="1"])
  l1   <- mean(PP_obs$Prev[PP_obs$PGS<3  & PP_obs$GROUP=="1"])
  r0   <- mean(PP_obs$Prev[PP_obs$PGS>98 & PP_obs$GROUP=="0"])
  l0   <- mean(PP_obs$Prev[PP_obs$PGS<3  & PP_obs$GROUP=="0"])
  deltaObs <- (r1-r0)-(l1-l0)
  
  # Delta null distribution
  sigma_ui <- summary(lm(ui ~ meanPRS + ENV, data = PP_add))$sigma
  PP_null_iter <- lapply(1:50, function(i) {
    est_ui <- PP_add$Estimated_ui + rnorm(nrow(PP_add), 0, sigma_ui)
    data.frame(
      PGS       = PP_add$PGS,
      GROUP     = PP_add$GROUP,
      Prev_iter = 100 * (1 - pnorm(qnorm(1 - overallprev),
                                   mean = est_ui,
                                   sd   = 1))
    )
  }) %>% bind_rows(.id="iter")
  deltas_null <- PP_null_iter %>%
    group_by(iter) %>%
    summarise(
      r1n = mean(Prev_iter[PGS>98 & GROUP=="1"]),
      l1n = mean(Prev_iter[PGS<3  & GROUP=="1"]),
      r0n = mean(Prev_iter[PGS>98 & GROUP=="0"]),
      l0n = mean(Prev_iter[PGS<3  & GROUP=="0"])
    ) %>%
    mutate(delta = (r1n - r0n) - (l1n - l0n))
  deltaNull <- mean(deltas_null$delta)
  sdNull    <- sd(deltas_null$delta)
  departure <- (deltaObs - deltaNull) / sdNull
  
  # Adjusted R2 additive model
  r2a    <- 1 - var(PP_add$Prev - PP_add$ExpectedPrev_nostoch)/var(PP_add$Prev)
  n <- nrow(PP_add)
  p <- 2 
  adjR2a <- 1 - ((1 - r2a) * (n - 1) / (n - p - 1))
  
  # Adjusted R2 interaction model
  r2i    <- 1 - var(PP_int$Prev - PP_int$ExpectedPrev_nostoch)/var(PP_int$Prev)
  n <- nrow(PP_int)
  p <- 3 
  adjR2i <- 1 - ((1 - r2i) * (n - 1) / (n - p - 1))
  
  data.frame(
    Scenario       = name,
    deltaObs       = round(deltaObs,1),
    departure      = round(departure,2),
    Additive_R2    = round(adjR2a,3),
    Interaction_R2 = round(adjR2i,3)
  )
}) %>% bind_rows()

print(results)

# get null dataframe
null_out    <- generate_PP(scenarios$Null$y, scenarios$Null$E)
overallprev <- mean(null_out$D)
PP_null_add <- add_comp(null_out$PP, overallprev)

# null layer from addtive model to be overlaid on observed curves
null_layer <- geom_smooth(
  data     = PP_null_add,
  aes(x = PGS, y = ExpectedPrev_nostoch * 100, group = GROUP),
  method   = "lm",
  formula  = y ~ poly(x, 3),
  se       = FALSE,
  color    = "grey30",
  linetype = "dashed",
  size     = 1
) 

plot_titles <- c(
  Null = "Null: Additive",
  Alt  = "Alt: Interaction",
  Het  = "Heteroscedasticity",
  GEC  = "G–E Correlation"
)

# Generate prevalence vs percentile PGS relationships across all scenarios
pps <- lapply(names(scenarios), function(name) {
  scen <- scenarios[[name]]
  if (name == "GEC") {
    generate_PP(scen$y, scen$E, custom_PRS = scen$PRS)$PP
  } else {
    generate_PP(scen$y, scen$E)$PP
  }
})
names(pps) <- names(scenarios)


# Step 2: Get max prevalence across all plots
max_prev <- max(sapply(pps, function(pp) max(pp$Prev, na.rm = TRUE)))
y_limit <- c(-4, ceiling(max_prev + 6))


plots <- lapply(names(scenarios), function(name) {
  scen <- scenarios[[name]]
  
  if (name == "GEC") {
    out <- generate_PP(scen$y, scen$E, custom_PRS = scen$PRS)
  } else {
    out <- generate_PP(scen$y, scen$E)
  }
  
  PP_obs <- out$PP
  overallprev <- mean(out$D)
  PP_add <- add_comp(PP_obs, overallprev)
  PP_add$ExpectedPrev_nostoch <- PP_add$ExpectedPrev_nostoch * 100
  
  
  ggplot(PP_obs, aes(x = PGS, y = Prev, color = GROUP)) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE) +
    geom_point(size = 0.8) +
    null_layer +
    labs(
      title = plot_titles[name],
      x     = "Percentile PGS",
      y     = "Prevalence (%)"
    ) +
    theme_minimal() +
    scale_y_continuous(limits = y_limit)
})


# Combine plots
sf <- wrap_plots(plotlist = plots, ncol = 2)
sf


#########################################################################################################################

# 100 iterations

nIter <- 100
all_results <- list()

# Get null distributions
null_scen <- scenarios$Null
null_deltas <- replicate(100, {
  tmp <- if ("PRS" %in% names(null_scen)) {
    generate_PP(null_scen$y, null_scen$E, custom_PRS = null_scen$PRS)
  } else {
    generate_PP(null_scen$y, null_scen$E)
  }
  PP_obs <- tmp$PP
  overallprev <- mean(tmp$D)
  PP_add <- add_comp(PP_obs, overallprev)
  
  r1 <- mean(PP_obs$Prev[PP_obs$PGS > 98 & PP_obs$GROUP == "1"])
  l1 <- mean(PP_obs$Prev[PP_obs$PGS < 3  & PP_obs$GROUP == "1"])
  r0 <- mean(PP_obs$Prev[PP_obs$PGS > 98 & PP_obs$GROUP == "0"])
  l0 <- mean(PP_obs$Prev[PP_obs$PGS < 3  & PP_obs$GROUP == "0"])
  (r1 - r0) - (l1 - l0)
})

deltaNull <- mean(null_deltas)
sdNull <- sd(null_deltas)


for (i in 1:nIter) {
  message("Running iteration ", i)
  iter_results <- lapply(names(scenarios), function(name) {
    scen        <- scenarios[[name]]
    tmp <- if (name == "GEC") {
      generate_PP(scen$y, scen$E, custom_PRS = scen$PRS)
    } else {
      generate_PP(scen$y, scen$E)
    }
    
    PP_obs      <- tmp$PP
    overallprev <- mean(tmp$D)
    
    PP_add <- add_comp(PP_obs, overallprev)
    PP_add$ExpectedPrev_nostoch <- PP_add$ExpectedPrev_nostoch * 100
    PP_int <- int_comp(PP_obs, overallprev)
    PP_int$ExpectedPrev_nostoch <- PP_int$ExpectedPrev_nostoch * 100
    
    # deltaObs
    r1 <- mean(PP_obs$Prev[PP_obs$PGS > 98 & PP_obs$GROUP == "1"])
    l1 <- mean(PP_obs$Prev[PP_obs$PGS < 3  & PP_obs$GROUP == "1"])
    r0 <- mean(PP_obs$Prev[PP_obs$PGS > 98 & PP_obs$GROUP == "0"])
    l0 <- mean(PP_obs$Prev[PP_obs$PGS < 3  & PP_obs$GROUP == "0"])
    deltaObs <- (r1 - r0) - (l1 - l0)
    
    # Standardized departure
    departure_sdu <- (deltaObs - deltaNull) / sdNull
    
    # R²
    r2a <- 1 - var(PP_add$Prev - PP_add$ExpectedPrev_nostoch) / var(PP_add$Prev)
    r2i <- 1 - var(PP_int$Prev - PP_int$ExpectedPrev_nostoch) / var(PP_int$Prev)
    
    n_a <- nrow(PP_add); n_i <- nrow(PP_int)

    #Adjusted R2
    adjR2a <- 1 - ((1 - r2a) * (n_a - 1) / (n_a - 3))
    adjR2i <- 1 - ((1 - r2i) * (n_i - 1) / (n_i - 4))
    
    data.frame(
      Iter = i,
      Scenario = name,
      departure = round(departure_sdu, 2),
      Additive_R2 = round(adjR2a, 3),
      Interaction_R2 = round(adjR2i, 3)
    )
  })
  all_results[[i]] <- bind_rows(iter_results)
}

final_results <- bind_rows(all_results)


###########################################################################################
# Ploting

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Set scenario levels and labels
final_results$Scenario <- factor(final_results$Scenario,
                                 levels = c("Null", "GEC", "Het", "Alt"),
                                 labels = c("Null (Additive)", "GEcorrelation",
                                            "Heteroscedasticity", "Alt (Interaction)"))

# Departure boxplot
p1 <- ggplot(final_results, aes(x = Scenario, y = departure, fill = Scenario)) +
  geom_boxplot(fill = "#DBD3D3") +
  labs(title = "Departure from Additive Null", y = "Departure (sdu)", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        legend.position = "none")

# Additive vs interaction R2
r2_long <- final_results %>%
  pivot_longer(cols = c(Additive_R2, Interaction_R2),
               names_to = "Model", values_to = "R2")

r2_long$Scenario <- factor(r2_long$Scenario,
                           levels = c("Null (Additive)", "GEcorrelation",
                                      "Heteroscedasticity", "Alt (Interaction)"))

p2 <- ggplot(r2_long, aes(x = Scenario, y = R2, fill = Model)) +
  geom_boxplot(position = position_dodge(), color = "black") +
  scale_fill_manual(values = c("Additive_R2" = "#E5E1DA", "Interaction_R2" = "#89A8B2"),
                    labels = c("Additive R²", "Interaction R²")) +
  labs(title = "Additive and Interaction R²", y = "Adjusted R²", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        legend.position = "right")



library(patchwork)
combined_plot <- (p1 + p2) +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Niter=100, Nsamples=100000, h2=0.2, Ncausal=100, Prevalence=12% ",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

#########################################################################################################################
