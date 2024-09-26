
delta_comp <- function(envcomb, PP_full){    
    delf <- data.frame()
    for(it in 1:50){
        deldf <- data.frame()
        count=1
        for(e in envcomb){
            print(e)
            PP_summ <- PP_full[PP_full$iteration == it,]
            obs_h <- PP_summ[PP_summ$ENVCODE == e[2],]
            obs_l <- PP_summ[PP_summ$ENVCODE == e[1],]

            rightdev <- mean(obs_h$Prev[obs_h$PGS %in% c(96, 97, 98, 99, 100)]) - mean(obs_l$Prev[obs_l$PGS %in% c(96, 97, 98, 99, 100)])
            leftdev <- mean(obs_h$Prev[obs_h$PGS %in% c(1, 2,3,4,5)]) - mean(obs_l$Prev[obs_l$PGS %in% c(1, 2,3,4,5)])
            delobs5 <- abs(rightdev) - abs(leftdev)

            rightdev <- mean(obs_h$Prev[obs_h$PGS %in% c(99, 100)]) - mean(obs_l$Prev[obs_l$PGS %in% c(99, 100)])
            leftdev <- mean(obs_h$Prev[obs_h$PGS %in% c(1,2)]) - mean(obs_l$Prev[obs_l$PGS %in% c(1, 2)])
            delobs2 <- abs(rightdev) - abs(leftdev)
            rightdev2 <- rightdev
            leftdev2 <- leftdev

            exp_h <- obs_h
            exp_l <- obs_l

            rightdev <- mean(exp_h$ExpectedPrev[exp_h$PGS %in% c(96, 97, 98, 99, 100)]) - mean(exp_l$ExpectedPrev[exp_l$PGS %in% c(96, 97, 98, 99, 100)])
            leftdev <- mean(exp_h$ExpectedPrev[exp_h$PGS %in% c(1, 2,3,4,5)]) - mean(exp_l$ExpectedPrev[exp_l$PGS %in% c(1, 2,3,4,5)])
            delexp5 <- abs(rightdev) - abs(leftdev)

            rightdev <- mean(exp_h$ExpectedPrev[exp_h$PGS %in% c(99, 100)]) - mean(exp_l$ExpectedPrev[exp_l$PGS %in% c(99, 100)])
            leftdev <- mean(exp_h$ExpectedPrev[exp_h$PGS %in% c(1,2)]) - mean(exp_l$ExpectedPrev[exp_l$PGS %in% c(1, 2)])
            delexp2 <- abs(rightdev) - abs(leftdev)
            rightdevExp2 <- rightdev
            leftdevExp2 <- leftdev


            rightdev_expstoc <- mean(exp_h$ExpectedPrev_nostoch[exp_h$PGS %in% c(99, 100)]) - mean(exp_l$ExpectedPrev_nostoch[exp_l$PGS %in% c(99, 100)])
            leftdev_expstoc <- mean(exp_h$ExpectedPrev_nostoch[exp_h$PGS %in% c(1,2)]) - mean(exp_l$ExpectedPrev_nostoch[exp_l$PGS %in% c(1, 2)])
            delexp2_expstoc <- abs(rightdev_expstoc) - abs(leftdev_expstoc)
            #rightdevExp2_expstoc <- rightdev_expstoc
            #leftdevExp2_expstoc <- leftdev_expstoc

            deldx <- data.frame("Delta" = paste0("Delta", count),
                                "ENV" = paste0(e[1],"-", e[2]),
                                "ENVCODE1" = e[1],"ENVCODE2" = e[2], 
                                "rightdevObs2"= rightdev2,"leftdevObs2" = leftdev2,
                                "rightdevExp2"= rightdevExp2,"leftdevExp2" = leftdevExp2,
                                "rightdevExp2_nostoc"= rightdev_expstoc,"leftdevExp2_nostoc" = leftdev_expstoc,
                              "deltaObs2" = delobs2, "deltaObs5" = delobs5, 
                               "deltaExp2_nostoc" = delexp2_expstoc, 
                              "deltaExp2" = delexp2, "deltaExp5" = delexp5, "iteration" = it, 
                              "ENVfield1" = PP_summ$field1[1], 
                              "ENVfield2" = PP_summ$field2[1], 
                              "ENVdesc" = PP_summ$Envdesc[1])
            deldf <- rbind(deldf, deldx)
            count=count+1
          
        }
    delf <- rbind(delf, deldf)
    }
    return(delf)
}