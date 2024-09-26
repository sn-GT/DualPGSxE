
int_comp <- function(PP, overallprev){
    # overall threshold t from overall population prevalence
    t_overall = qnorm(1-overallprev) #inverse cdf to get t 

    # Env (1/0) based on high vs low group of environmental exposure
    #PP$ENV1[PP$GROUP == paste0(e1[2],"-", e2[2], "\n(00)")] <- 0
    #PP$ENV2[PP$GROUP == paste0(e1[2],"-", e2[2], "\n(00)")] <- 0

    #PP$ENV1[PP$GROUP == paste0(e1[2],"-", e2[1], "\n(01)")] <- 0
    #PP$ENV2[PP$GROUP == paste0(e1[2],"-", e2[1], "\n(01)")] <- 1

    #PP$ENV1[PP$GROUP == paste0(e1[1],"-", e2[2], "\n(10)")] <- 1
    #PP$ENV2[PP$GROUP == paste0(e1[1],"-", e2[2], "\n(10)")] <- 0

    #PP$ENV1[PP$GROUP == paste0(e1[1],"-", e2[1], "\n(11)")] <- 1
    #PP$ENV2[PP$GROUP == paste0(e1[1],"-", e2[1], "\n(11)")] <- 1
    #######

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
    #PP <- PP[PP$Prev !=0,] ############ change this how to deal with 0 prev cases)
    model <- lm(data = PP, ui ~ meanPRS + ENV1 + ENV2 + meanPRS*ENV1 + meanPRS*ENV2 + meanPRS*ENV1*ENV2)
    model
    X = model.matrix(model)
    #drop = which(colnames(X) == 'ENV1:ENV2')
    X1 = X[,-7]
    newfit = lm(PP$ui ~ X1-1)

    #modfit <- summary(model)
    modfit <- summary(newfit)
    model <- newfit
    # get regression coefficients
    # mean PRS
    a1= model$coefficients[2]
    # env1
    a2 = model$coefficients[3]
    # env2
    a3 = model$coefficients[4]
    # interaction env1 
    a4 = model$coefficients[5]
    # interaction env2 
    a5 = model$coefficients[6]
    # interaction all
    a7 = model$coefficients[7]
    # intercept
    a8 = model$coefficients[1]

    PP$meanPRScoeff <- a1
    PP$meanPRSpval <- modfit$coefficients[2,4]

    PP$env1coeff <- a2
    PP$env1pval <- modfit$coefficients[3,4]
    
    PP$env2coeff <- a3
    PP$env2pval <- modfit$coefficients[4,4]

    PP$env1intcoeff <- a4
    PP$env1intpval <- modfit$coefficients[5,4]

    PP$env2intcoeff <- a5
    PP$env2intpval <- modfit$coefficients[6,4]

    PP$intallcoeff <- a7
    PP$intallpval <- modfit$coefficients[7,4]

    PP$intcpt <- a8
    PP$intcptpval <- modfit$coefficients[1,4]


    # estimated ui from regression coeff
    PP$Estimated_ui <- a1*PP$meanPRS + a2*PP$ENV1 + a3*PP$ENV2 + 
      (a4*PP$meanPRS*PP$ENV1) +(a5*PP$meanPRS*PP$ENV2) + 
      #(a6*PP$ENV1*PP$ENV2) + 
      (a7*PP$meanPRS*PP$ENV1*PP$ENV2) + a8

    PP$ExpectedPrev_nostoch <- (1 - pnorm(PP$t, mean = PP$Estimated_ui, sd = 1))
    
    PP$InteractionR2 <- modfit$adj.r.squared
    
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
