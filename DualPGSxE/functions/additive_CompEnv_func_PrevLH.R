# additve curves (10 iterations) for composite environment

# inputs - Prevperc PP, overallprev
add_comp <- function(PP, overallprev){
    # overall threshold t from overall population prevalence
    t_overall = qnorm(1-overallprev) #inverse cdf to get t 
    
    PP$ENV1[PP$ENVCODE == "00"] <- 0
    PP$ENV2[PP$ENVCODE == "00"] <- 0
    
    PP$ENV1[PP$ENVCODE == "01"] <- 0
    PP$ENV2[PP$ENVCODE == "01"] <- 1
    
    PP$ENV1[PP$ENVCODE == "10"] <- 1
    PP$ENV2[PP$ENVCODE == "10"] <- 0
    
    PP$ENV1[PP$ENVCODE == "11"] <- 1
    PP$ENV2[PP$ENVCODE == "11"] <- 1

    # Compute threshold t at each percentile PGS based on observed prevalence at each PGS.
    PP$Prev <- PP$Prev/100
    PP$Prev[PP$Prev == 0] <- 0.001
    PP$Prev[PP$Prev == 1] <- 0.99
    PP$t_PGS <- qnorm(1-PP$Prev)
    PP$overallP <- overallprev
    PP$t <- t_overall

    # Compute underlying mean ui of the liability model thresholds  
    PP$ui <- PP$t - PP$t_PGS

    # Fit linear regression model, ui = t-invcdf(1-Pi) = a meanPRSi + b Env(0/1) + c
    model <- lm(data = PP, ui ~ meanPRS + ENV1 + ENV2)
    modfit <- summary(model)
    modfit

    prscoef = 2
    env1coef = 3
    env2coef = 4
    intcpt = 1

    # get regression coefficients
    # mean PRS
    a1= model$coefficients[2]
    # env1 
    a2 = model$coefficients[3]
    # env2
    a3 = model$coefficients[4]
    # intercept
    a4 = model$coefficients[1]

    PP$meanPRScoeff <- a1
    PP$meanPRSpval <- modfit$coefficients[prscoef,4]

    PP$env1coeff <- a2
    PP$env1pval <- modfit$coefficients[env1coef,4]
    
    PP$env2coeff <- a3
    PP$env2pval <- modfit$coefficients[env2coef,4]

    PP$intcpt <- a4
    PP$intcptpval <- modfit$coefficients[intcpt,4]

    # estimated ui from regression coeff
    PP$Estimated_ui <- a1*PP$meanPRS + a2*PP$ENV1 +  a3*PP$ENV2  + a4

    PP$ExpectedPrev_nostoch <- (1 - pnorm(PP$t, mean = PP$Estimated_ui, sd = 1))
    
    PP$AdditiveR2 <- modfit$adj.r.squared
    # sd = residual error from linear regression model
    modsigma <- (modfit$sigma)

    #PP$env <- paste0(PP$ENV1, PP$ENV2)

    PP_full <- data.frame()
    if(nrow(PP[PP$Prev == 0,]) == 0){
      df <- data.frame()
      # 10 iterations to add stochasticity to underlying estimated ui
      # Computed expected prevalence from estimated ui:  1-cdf(N(ui,1), toverall)
      for(it in 1:50){
        exprev <- c()
        for(estui in PP$Estimated_ui){
          estui_se <- rnorm(1, estui, modsigma)
          # expected prevalence from estimated ui
          expecprev_sub <- (1 - pnorm(PP$t[1], mean = estui_se, sd = 1))
          exprev <- c(exprev, expecprev_sub)
        }
        PP$ExpectedPrev <- exprev
        PP$iteration <- it
        df <- rbind(df, PP) 
      }
      PP_full <- rbind(PP_full, df)
    }
    return(PP_full)
}

# save coeffcients and their p-values also