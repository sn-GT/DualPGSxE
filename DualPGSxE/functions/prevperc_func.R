# Prevalence vs percentile function

# Input is dataframe m with the following headers - 
# FID, GROUP (env groups), invGRS (PRS score), BCA (case-control status)
prevperc <- function(m){	
	PP <- data.frame()
    l <- names(table(m$GROUP))
    for(t in l){
      br_grs <- m
      table(br_grs$BCA)
      n=100
      df <- data.frame()
      for(i in 0:(n-1)){
        print(i)
        if(i == 0){
          dx <- br_grs[br_grs$GROUP == t,]
          dx1 <- dx[which(dx$invGRS >= quantile(dx$invGRS, i/n) & dx$invGRS <= quantile(dx$invGRS, (i+1)/n)),]
        }
        else {
          dx <- br_grs[br_grs$GROUP == t,]
          dx1 <- dx[which(dx$invGRS > quantile(dx$invGRS, i/n) & dx$invGRS <= quantile(dx$invGRS, (i+1)/n)),]
        }
        if(nrow(dx1) !=0){
        #prev = data.frame(table(dx1$BCA)/nrow(dx1))[2,2] * 100
        prev = (table(dx1$BCA)["1"]/nrow(dx1)) * 100
        if(is.na(prev)){prev = 0}
        df2 = data.frame("PGS" = i+1, "Prev" = prev , "meanPRS" = mean(dx1$invGRS) ,"n" = nrow(dx1), "ENVCODE"=dx1$ENVCODE[1]) 
        df <- rbind(df, df2)
        }
      }
      df$GROUP <- paste0(t)
      df$scaled <- scale(df$Prev, scale = F)
      df$zeroes_perGroup <- nrow(df[df$Prev == 0,])
      PP <- rbind(PP, df)
    }
    return(PP)
}