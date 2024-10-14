require(fitdistrplus)
require(pracma)
require(stringr)
require(berryFunctions)
require(xtable)
path <- str_split(rstudioapi::getSourceEditorContext()$path,"15_percent")[[1]][1]
setwd(path)
set.seed(1)

#############################
###  F U N C T I O N S ######
#############################

weibull_shift <-function(pars,mit,d_org){
  
  r_new <- rweibull(1000000,shape = pars ,scale = mit/gamma(1+ 1/pars))
  

  Xs_org <- d_org$x
  Ys_org <- d_org$y

  if(any(is.na(r_new))==F){
    d_new   <- density(r_new,n = length(Xs_org),from=min(Xs_org),to=max(max(Xs_org)))
    res     <- trapz(x=Xs_org,y=abs(Ys_org - d_new$y)) 
  }else{res <- 1000}
  return(res)
}

optim_classes<-function(pars,opti = T,classes,shares,m){
  
  s <- rweibull(1000000,shape= pars[1], scale = pars[2])
  
  diffs <- NULL
  perc <- NULL
  
  shares_modeled <- NULL
  
  # Calculate Difference for later evaluation
  for(i in 1:(nrow(classes)-1)){
    diffs <- c(diffs,((length(which(s >= classes[i,"from"] & s<= classes[i,"to"])) / (length(s))*100) - shares[i]))
  }
  diffs <- c(diffs,((length(which(s >= classes[nrow(classes),"from"])) / (length(s))*100) - shares[nrow(classes)]))
  diffs <- c(diffs,(mean(s) - m))
  
  # Calculate Difference in percentage of cumulative shares 
  for( i in 1: (nrow(classes)-1)){
    shares_modeled <- c(shares_modeled,length(which(s >= classes[i,"from"] & s<= classes[i,"to"])) / (length(s))*100)
  }
  shares_modeled <- c(shares_modeled,length(which(s >= classes[nrow(classes),"from"])) / (length(s))*100)
  
  diffs_perc <- 100 * (cumsum(shares) - cumsum(shares_modeled)) / cumsum(shares)
  diffs_perc <- c(diffs_perc,100 * (mean(s) - m)/ mean(s))

  
  # scaling
  scaled_diff_perc <-  c(diffs_perc,rep(diffs_perc[length(diffs_perc)],nrow(classes)-1))
  
  # sum of squares
  squareddiffs <- scaled_diff_perc^2
  
  if(opti == T)res <- sum(squareddiffs)
  if(opti == F){
    res <- list()
    
    res[["eval"]] <-  diffs
    res[["eval_perc"]] <- diffs_perc
    res[["eval_cum"]] <-  cumsum(shares_modeled)
    res[["val"]] <- s
  }
  all_s[[count]] <<- s
  count <<- count + 1
  return(res)
}



#######################################
###  R E S I D E N T I A L  ###########
#######################################
len <- 10

#### Residential ####

#load data
trend      <- read.csv2("./data/heizenergieverbrauch_trend.csv")
distro     <- read.csv2("./data/heizenergieverbrauch_verteilung.csv")

res_values <- NULL
head(distro)

for(i in 1:nrow(distro)){
  res_values <- c(res_values,rep(distro[i,1],distro[i,2]*100000))
}

#UBA

res_energy_2020    <- trend[which(trend$Jahr==2020),"barlabel.gesamt"] 

res_2045_zi_40     <- trend[which(trend$Jahr==2008),"barlabel.gesamt"] * 0.6
res_2045_zi_55     <- trend[which(trend$Jahr==2008),"barlabel.gesamt"] * 0.45
res_2045_zi_70     <- trend[which(trend$Jahr==2008),"barlabel.gesamt"] * 0.35

res_2045_agora_y45  <- trend[which(trend$Jahr==2018),"barlabel.gesamt"] * 0.68
res_2045_agora_y50  <- trend[which(trend$Jahr==2011),"barlabel.gesamt"] * 0.56


lin_interpol_m_zi_40     <- (res_energy_2020 - res_2045_zi_40) / (2050-2020)
lin_interpol_m_zi_55     <- (res_energy_2020 - res_2045_zi_55) / (2050-2020)
lin_interpol_m_zi_70     <- (res_energy_2020 - res_2045_zi_70) / (2050-2020)
lin_interpol_m_agora_y45 <- (res_energy_2020 - res_2045_agora_y45) / (2045-2020)
lin_interpol_m_agora_y50 <- (res_energy_2020 - res_2045_agora_y50) / (2050-2020)


res_2030_zi_40 <- res_energy_2020 -(lin_interpol_m_zi_40* (2030-2020))
res_2040_zi_40 <- res_energy_2020 -(lin_interpol_m_zi_40* (2040-2020))

res_2030_zi_55 <- res_energy_2020 -(lin_interpol_m_zi_55* (2030-2020))
res_2040_zi_55 <- res_energy_2020 -(lin_interpol_m_zi_55* (2040-2020))

res_2030_zi_70 <- res_energy_2020 -(lin_interpol_m_zi_70* (2030-2020))
res_2040_zi_70 <- res_energy_2020 -(lin_interpol_m_zi_70* (2040-2020))

res_2030_agora_y45 <- res_energy_2020 -(lin_interpol_m_agora_y45* (2030-2020))
res_2040_agora_y45 <- res_energy_2020 -(lin_interpol_m_agora_y45* (2040-2020))

res_2030_agora_y50 <- res_energy_2020 -(lin_interpol_m_agora_y50* (2030-2020))
res_2040_agora_y50 <- res_energy_2020 -(lin_interpol_m_agora_y50* (2040-2020))


y <- trend$barlabel.gesamt
x <- 1:length(y)


#fit distribution

classes_res <-  data.frame("from"= c(0,30,50,75,100,130,160,200,250),
                           "to"= c(29,49,74,99,129,159,199,249,NA))

share_res   <- c(2.5,2.32,7.97,11.84,17.31,12.91,14.04,12.67,18.44)                     
mean_res    <- trend[nrow(trend),"barlabel.gesamt"]

fit_org_res_raw <- fitdist(rep(classes_res$from+1,share_res*100),"weibull")

q_test        <- numeric(len)
val_test_list <- list()
par_opt_list  <- list()
goals         <- numeric(len)

for(i in 1:len){
  
  count <<- 1 
  all_s <<- list()
  
  print(i)
  optimized <- optim(par=fit_org_res_raw$estimate,
                     fn=optim_classes,
                     classes = classes_res,
                     shares=share_res,
                     m=mean_res,
                     control=list(trace=1,maxit = 200),
                     method="SANN")
  
  result<- optim_classes(pars=optimized$par,
                         classes = classes_res,
                         shares=share_res,
                         m=mean_res,
                         opti=F)
  
  fit_test <- fitdistr(result$val,"weibull")
  
  val_test_list[[i]] <- result$val
  par_opt_list[[i]]  <- optimized$par
  goals[i]           <- optimized$value
  q_test[i]          <- round(qweibull(0.15, shape = fit_test$estimate["shape"],scale = fit_test$estimate["scale"]))
}  

chosen_idx <- which(goals == min(goals))
d_org_res <- density(val_test_list[[chosen_idx]])

fit_org_res <- fitdistr(val_test_list[[chosen_idx]],"weibull")
q_org_res   <- round(qweibull(0.15, shape = fit_org_res$estimate["shape"],scale = fit_org_res$estimate["scale"]))


#update fit for future


p_res_2030_zi_40 <- optim(par=fit_org_res$estimate["shape"],
                          fn=weibull_shift,
                          mit=res_2030_zi_40,
                          d_org = d_org_res,
                          control=list(trace=1),
                          method="Brent",
                          lower=0,
                          upper=100)$par
p_res_2030_zi_40 <- c(p_res_2030_zi_40,res_2030_zi_40 / gamma(1+ 1/p_res_2030_zi_40) )

p_res_2040_zi_40 <- optim(par=fit_org_res$estimate["shape"],
                          fn=weibull_shift,
                          mit=res_2040_zi_40,
                          d_org = d_org_res,
                          control=list(trace=1),
                          method="Brent",
                          lower=0,
                          upper=100)$par
p_res_2040_zi_40 <- c(p_res_2040_zi_40,res_2040_zi_40 / gamma(1+ 1/p_res_2040_zi_40) )

p_res_2030_zi_70 <- optim(par=fit_org_res$estimate["shape"],
                          fn=weibull_shift,
                          mit=res_2030_zi_70,
                          d_org = d_org_res,
                          control=list(trace=1),
                          method="Brent",
                          lower=0,
                          upper=100)$par
p_res_2030_zi_70 <- c(p_res_2030_zi_70,res_2030_zi_70 / gamma(1+ 1/p_res_2030_zi_70) )


p_res_2040_zi_70 <- optim(par=fit_org_res$estimate["shape"],
                          fn=weibull_shift,
                          mit=res_2040_zi_70,
                          d_org = d_org_res,
                          control=list(trace=1),
                          method="Brent",
                          lower=0,
                          upper=100)$par
p_res_2040_zi_70 <- c(p_res_2040_zi_70,res_2040_zi_70 / gamma(1+ 1/p_res_2040_zi_70) )


q_res_2030_zi_40 <- qweibull(0.15, shape= p_res_2030_zi_40[1] ,scale = p_res_2030_zi_40[2])
q_res_2040_zi_40 <- qweibull(0.15, shape= p_res_2040_zi_40[1] ,scale = p_res_2040_zi_40[2])
q_res_2030_zi_70 <- qweibull(0.15, shape= p_res_2030_zi_70[1] ,scale = p_res_2030_zi_70[2])
q_res_2040_zi_70 <- qweibull(0.15, shape= p_res_2040_zi_70[1] ,scale = p_res_2040_zi_70[2])

q_df_res <-round(data.frame("UBA_zi_40"=c(q_res_2030_zi_40,q_res_2040_zi_40),"UBA_zi_70"=c(q_res_2030_zi_70,q_res_2040_zi_70)))


# Time dependency 
q_vec_zi_40_df <- as.data.frame(matrix(ncol=25,nrow=10))
q_vec_zi_70_df <- as.data.frame(matrix(ncol=25,nrow=10))


for(k in 1:len){

  
q_vec_zi_40 <- q_org_res
q_vec_zi_70 <- q_org_res

for(i in 1:24){
  print(i)
  trend_zi_40 <- res_energy_2020 -(lin_interpol_m_zi_40* ((2020+i)-2020))
  trend_zi_70 <- res_energy_2020 -(lin_interpol_m_zi_70* ((2020+i)-2020))
  
  
  res_zi_40 <- optim(par=fit_org_res$estimate["shape"],
                     fn=weibull_shift,
                     mit=trend_zi_40,
                     d_org = d_org_res,
                     control=list(trace=1),
                     method="Brent",
                     lower=0,
                     upper=100)$par
  
  res_zi_70 <- optim(par=fit_org_res$estimate["shape"],
                     fn=weibull_shift,
                     mit=trend_zi_70,
                     d_org = d_org_res,
                     control=list(trace=1),
                     method="Brent",
                     lower=0,
                     upper=100)$par

  p_zi_40 <- c(res_zi_40,trend_zi_40 / gamma(1+ 1/res_zi_40))
  p_zi_70 <- c(res_zi_70,trend_zi_70 / gamma(1+ 1/res_zi_70))
  
  q_zi_40 <- qweibull(0.15, shape= p_zi_40[1] ,scale = p_zi_40[2])
  q_zi_70 <- qweibull(0.15, shape= p_zi_70[1] ,scale = p_zi_70[2])
  
  q_vec_zi_40 <- c(q_vec_zi_40,q_zi_40)
  q_vec_zi_70 <- c(q_vec_zi_70,q_zi_70)
  
}

q_vec_zi_40_df[k,] <- q_vec_zi_40
q_vec_zi_70_df[k,] <- q_vec_zi_70

}

q_vec_zi_40    <- apply(q_vec_zi_40_df,2,FUN=mean,na.rm=T,simplify=T)
q_vec_zi_40_sd <- apply(q_vec_zi_40_df,2,FUN=sd,na.rm=T,simplify=T)

q_vec_zi_70    <- apply(q_vec_zi_70_df,2,FUN=mean,na.rm=T,simplify=T)
q_vec_zi_70_sd <- apply(q_vec_zi_70_df,2,FUN=sd,na.rm=T,simplify=T)



################################################
### N O N  -  R E S I D E N T I A L  ###########
################################################

#### Non-Residential ####
nonres_values      <- read.csv2("./data/nwg.csv", header=FALSE)[,1]

# UBA
nonres_energy_2020       <- mean(nonres_values)
nonres_2050_zi_25    <- 0.75*nonres_energy_2020 #126
nonres_2050_zi_45    <- 0.55*nonres_energy_2020 #96

lin_interpol_m_zi_25 <- (nonres_energy_2020 - nonres_2050_zi_25) / (2050-2020)
lin_interpol_m_zi_45 <- (nonres_energy_2020 - nonres_2050_zi_45) / (2050-2020)

nonres_2030_zi_25 <- nonres_energy_2020 -(lin_interpol_m_zi_25* (2030-2020))
nonres_2040_zi_25 <- nonres_energy_2020 -(lin_interpol_m_zi_25* (2040-2020))

nonres_2030_zi_45 <- nonres_energy_2020 -(lin_interpol_m_zi_45* (2030-2020))
nonres_2040_zi_45 <- nonres_energy_2020 -(lin_interpol_m_zi_45* (2040-2020))


#fit distribution

fit_org_nonres <- fitdist(nonres_values,"weibull")
d_org_nonres   <- density(rweibull(100000,shape = fit_org_nonres$estimate["shape"] ,scale = fit_org_nonres$estimate["scale"]),n = 1000,from=0,to=400)
q_org_nonres   <- round(qweibull(0.15, shape=fit_org_nonres$estimate["shape"],scale = fit_org_nonres$estimate["scale"]),0)

q_vec_zi_25_df <- as.data.frame(matrix(ncol=25,nrow=10))
q_vec_zi_45_df <- as.data.frame(matrix(ncol=25,nrow=10))


for(k in 1:len){

q_vec_zi_25 <- q_org_nonres
q_vec_zi_45 <- q_org_nonres


for(i in 1:24){
  print(i)
  trend_zi_25 <- nonres_energy_2020 -(lin_interpol_m_zi_25* ((2020+i)-2020))
  trend_zi_45 <- nonres_energy_2020 -(lin_interpol_m_zi_45* ((2020+i)-2020))
  
  
  nonres_zi_25 <- optim(par=fit_org_nonres$estimate["shape"],
                     fn=weibull_shift,
                     mit=trend_zi_25,
                     d_org = d_org_nonres,
                     control=list(trace=1),
                     method="Brent",
                     lower=0,
                     upper=100)$par
  
  nonres_zi_45 <- optim(par=fit_org_nonres$estimate["shape"],
                     fn=weibull_shift,
                     mit=trend_zi_45,
                     d_org = d_org_nonres,
                     control=list(trace=1),
                     method="Brent",
                     lower=0,
                     upper=100)$par
  
  p_zi_25 <- c(nonres_zi_25,trend_zi_25 / gamma(1+ 1/nonres_zi_25))
  p_zi_45 <- c(nonres_zi_45,trend_zi_45 / gamma(1+ 1/nonres_zi_45))
  
  q_zi_25 <- qweibull(0.15, shape= p_zi_25[1] ,scale = p_zi_25[2])
  q_zi_45 <- qweibull(0.15, shape= p_zi_45[1] ,scale = p_zi_45[2])
  
  q_vec_zi_25 <- c(q_vec_zi_25,q_zi_25)
  q_vec_zi_45 <- c(q_vec_zi_45,q_zi_45)
  
}

q_vec_zi_25_df[k,] <- q_vec_zi_25
q_vec_zi_45_df[k,] <- q_vec_zi_45
}

q_vec_zi_25    <- apply(q_vec_zi_25_df,2,FUN=mean,na.rm=T,simplify=T)
q_vec_zi_25_sd <- apply(q_vec_zi_25_df,2,FUN=sd,na.rm=T,simplify=T)

q_vec_zi_45    <- apply(q_vec_zi_45_df,2,FUN=mean,na.rm=T,simplify=T)
q_vec_zi_45_sd <- apply(q_vec_zi_45_df,2,FUN=sd,na.rm=T,simplify=T)

######################
### P L O T S   ######
######################

### Example Area Plot

png(paste(path,"tex/figures/fig3.png",sep=""),width=700,height=500)

mit = 150 
r_ex1 <- rgamma(1000000,shape = 4 ,rate = 4/mit)

mit = 140
r_ex2 <- rgamma(1000000,shape = 0.5 ,rate = 0.5/mit)
r_ex3 <- rgamma(1000000,shape = 4 ,rate = 4/mit)

d_ex1 <- density(r_ex1,from=0,to=400)
d_ex2 <- density(r_ex2,from=0,to=400)
d_ex3 <- density(r_ex3,from=0,to=400)

col_ex_2 <- rgb(0,0,1,alpha=0.5)
col_ex_3 <- rgb(1,0,0,alpha=0.5)

plot(d_ex1,col="black",type="l",xlim=c(0,400),ylim=c(0,0.013),
     main = "Areas between probability density functions",
     xlab="Energy Efficiency",
     lwd=6)
polygon(c(d_ex1$x,rev(d_ex1$x)),c(d_ex1$y,rev(d_ex2$y)),col=col_ex_2)
polygon(c(d_ex1$x,rev(d_ex1$x)),c(d_ex1$y,rev(d_ex3$y)),col=col_ex_3)
lines(d_ex1,col="black",lwd=6)
lines(d_ex2,col="blue",lwd=2)
lines(d_ex3,col="red",lwd=2)
legend("topright",c("Gamma 1 (mean: 150)","Gamma 2 (mean: 140)","Gamma 3 (mean: 140)","Area between 1 & 2","Area between 1 & 3"),col=c("black","blue","red",col_ex_2,col_ex_3),lty = c(1,1,1,NA,NA),pch=c(NA,NA,NA,15,15),lwd=c(6,2,2,NA,NA))
dev.off()

## Res 1 (Histogram)

col1 <- "grey"
col2 <- rgb(1,0,0,0.5)
col3 <- rgb(190/255,190/255,190/255,0.5)

png("./figures/fig6.png",width=700,height=500)
par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2,mar=c(5,5,6,3))
hist(res_values,prob=T,breaks=seq(from =0, to = max(res_values),by=10),ylim=c(0,0.008),col=col3,border=col3,main="",#"Histogram and fitted distributions \n (german residential buildings)",
     xlab="Final energy demand [kWh/(m^2*a)]")


lines(density(val_test_list[[chosen_idx]]),col=col2,lwd=3)
abline(v=q_test[chosen_idx],col="blue",lty=2,lwd=3)
axis(3,at=c(classes_res$from,400),padj=-20)
mtext(c("A+","A","B","C","D","E","F","G","H"),
      side=3,at=classes_res$from + (c(classes_res$to[-9],400) - classes_res$from)/2,line=0,cex=1.5)
legend("topright",c("weibull based on Statista","15% quantile"),
       col=c(col2,col2),lty=c(1,2),
       lwd=c(2,2),cex=1.5)
legend("topleft","co2online",fill=col1,cex=1.5)
dev.off()

## Res 2 (Fit classes)

result<- optim_classes(pars=fit_org_res$estimate,
                       classes = classes_res,
                       shares=share_res,
                       m=mean_res,
                       opti=F)

png("./figures/fig5_top.png",width=700,height=500)
par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2,mar=c(5,5,6,3))
layout(matrix(c(1,1,1,2), nrow = 1, ncol = 4, byrow = TRUE))
idx <- 1:(length(result$eval)-1)
barplot(share_res[idx] + result$eval[idx],col=col2
        ,main="Difference Energy Efficiency Classes",ylim=c(0,25),ylab="Share [%]",
        xlab="Energy Efficiency Class",width=1,
        names.arg = c("A+","A","B","C","D","E","F","G","H")[idx])
barplot(share_res[idx],col=col3,add=T,width=1)
#text(2,20,paste("Difference Mean:",abs(round(result$eval[length(result$eval)])),sep=" "),cex=2)
legend("topright",c("Fitted","Statista"),fill=c(col2,col1),cex=1.5)

idx <- length(result$eval)
barplot(mean_res + result$eval[idx],col=col2
        ,main="Difference Average",ylim=c(0,200),
        xlab="Average",width=1,ylab="Final energy demand [kWh/(m^2*a)]")
barplot(mean_res,col=col3,add=T,width=1)
legend("topright",c("Fitted","co2online"),fill=c(col2,col1),cex=1.5)
dev.off()

## Res 3 (Fit cumulated shares)

png("./figures/fig5_bot.png",width=700,height=500)
par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2,mar=c(5,5,6,3))
plot(cumsum(share_res),type="l",ylab="Cumulative Share [%]",xaxt="n",col=col1,lwd=2,xlab="Energy Efficiency Class")
points(cumsum(share_res),col=col1)

axis(1,at=c(1:9),c("A+","A","B","C","D","E","F","G","H"))
lines(result$eval_cum,col=col2,lwd=2)
points(result$eval_cum,col=col2)
abline(h=15,col="blue",lty=2,lwd=2)
legend("topleft",c("Fitted","Statista","15% quantile"),col=c(col2,col1,"red"),cex=1.5,lty=c(1,1,2))
dev.off()


### Trends 
png("./figures/fig4.png",width=700,height=500)
par(cex.axis=2, cex.lab=1.5, cex.main=2, cex.sub=2,mar=c(5,5,4,3))

plot(x+2001,y,main="Trends of average final energy demand \n according to scenario assumptions (interpolations)",
     xlim=c(2002,2045),ylim=c(0,180),xlab="Year",ylab="Average final energy demand [kWh/(m^2*a)]",lwd=2)
lines(c(x+2001,2020:2045),c(y,cumsum(c(y[length(y)],rep(-lin_interpol_m_zi_40,2045-2020)))),
      col="chartreuse1",
      lwd=2)
lines(c(x+2001,2020:2045),c(y,cumsum(c(y[length(y)],rep(-lin_interpol_m_zi_55,2045-2020)))),
      col="chartreuse3",
      lwd=2)
lines(c(x+2001,2020:2045),c(y,cumsum(c(y[length(y)],rep(-lin_interpol_m_zi_70,2045-2020)))),
      col="chartreuse4",
      lwd=2)
lines(c(x+2001,2020:2045),c(y,cumsum(c(y[length(y)],rep(-lin_interpol_m_agora_y45,2045-2020)))),
      col="blue",
      lwd=2)
lines(c(x+2001,2020:2045),c(y,cumsum(c(y[length(y)],rep(-lin_interpol_m_agora_y50,2045-2020)))),
      col="blue4",
      lwd=2)

lines(x+2001,y,col="black",lwd=2)

points(2030,res_2030_zi_40,pch=16,cex=2)
points(2040,res_2040_zi_40,pch=16,cex=2)

points(2030,res_2030_zi_70,pch=16,cex=2)
points(2040,res_2040_zi_70,pch=16,cex=2)

legend("bottomleft",c("UBA Target 40","UBA Target 55","UBA Target 70", "Agora year 2045", "Agora year 2050"),
       col=c("chartreuse1","chartreuse3","chartreuse4","blue","blue4"),
       lty=c(1,1,1,1,1),lwd=c(2,2,2,2,2),cex=2)
dev.off()


### shift
col_30 <- "purple"
col_40 <- "gold"

png("./figures/fig7_left.png",width=700,height=500)
par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2,mar=c(5,5,4,3))
plot(d_org_res, main="UBA target 40",
     ylim=c(0,0.01),xlab="Final energy demand [kWh/(m^2*a)]",xlim=c(0,400),lwd=2)
lines(density(rweibull(100000, p_res_2030_zi_40[1],p_res_2030_zi_40[2])),col=col_30,lwd=2)
lines(density(rweibull(100000, p_res_2040_zi_40[1],p_res_2040_zi_40[2])),col=col_40,lwd=2)
abline(v=q_org_res,lwd=2,lty=2)
abline(v=q_res_2030_zi_40,lwd=2,lty=2,col=col_30)
abline(v=q_res_2040_zi_40,lwd=2,lty=2,col=col_40)
legend("topright", c("2021","2030 ","2040","15% quantile"),lwd=c(2,2,2,2),lty=c(1,1,1,2),col=c("black",col_30,col_40,"grey"),cex=1.5)
dev.off()

png("./figures/fig7_right.png",width=700,height=500)
par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2,mar=c(5,5,4,3))
plot(d_org_res, main="UBA target 70",
     ylim=c(0,0.01),xlab="Final energy demand [kWh/(m^2*a)]",xlim=c(0,400),lwd=2)
lines(density(rweibull(100000, p_res_2030_zi_70[1],p_res_2030_zi_70[2])),col=col_30,lwd=2)
lines(density(rweibull(100000, p_res_2040_zi_70[1],p_res_2040_zi_70[2])),col=col_40,lwd=2)
abline(v=q_org_res,lwd=2,lty=2)
abline(v=q_res_2030_zi_70,lwd=2,lty=2,col=col_30)
abline(v=q_res_2040_zi_70,lwd=2,lty=2,col=col_40)
legend("topright", c("2021","2030 ","2040","15% quantile"),lwd=c(2,2,2,2),lty=c(1,1,1,2),col=c("black",col_30,col_40,"grey"),cex=1.5)
dev.off()


#### evolution threshold res 

png("./figures/fig8.png",width=700,height=500)
par(cex.axis=2, cex.lab=1.5, cex.main=1.5, cex.sub=2,mar=c(5,5,6,3))

col_trend <- rgb(col2rgb("grey")[1]/255,col2rgb("grey")[2]/255,col2rgb("grey")[3]/255,0.5) 
col_Aplus <- rgb(col2rgb("steelblue3")[1]/255,col2rgb("steelblue3")[2]/255,col2rgb("steelblue3")[3]/255,0.5) 
col_A     <- rgb(col2rgb("green4")[1]/255,col2rgb("green4")[2]/255,col2rgb("green4")[3]/255,0.5) 
col_B     <- rgb(col2rgb("green2")[1]/255,col2rgb("green2")[2]/255,col2rgb("green2")[3]/255,0.5) 
col_C     <- rgb(col2rgb("olivedrab3")[1]/255,col2rgb("olivedrab3")[2]/255,col2rgb("olivedrab3")[3]/255,0.5) 
col_D     <- rgb(col2rgb("yellow2")[1]/255,col2rgb("yellow2")[2]/255,col2rgb("yellow2")[3]/255,0.5) 

plot(NA,ylim=c(0,100),xaxt="n",xlim=c(0,25),xlab="Year",ylab="Final energy demand [kwh/(m^2*a)]",
     main = "Evolution of the 15% quantiles under different scenarios (residential)")

polygon(x=c(-5,-5,40,40),y=c(-5,30,30,-5),col=col_Aplus,border = F)  #A+
polygon(x=c(-5,-5,40,40),y=c(30,50,50,30),col=col_A,border = F)  #A+
polygon(x=c(-5,-5,40,40),y=c(50,75,75,50),col=col_B,border = F)  #A+
polygon(x=c(-5,-5,40,40),y=c(75,100,100,75),col=col_C,border = F)  #A+
polygon(x=c(-5,-5,40,40),y=c(100,150,150,100),col=col_D,border = F)  #A+

lines(q_vec_zi_40,lty=3,lwd=2)
lines(q_vec_zi_70,lty=2,lwd=2)

textField(x=0,y=25.5,"A+",cex=2)
textField(x=-0.4,y=45.5,"A",cex=2)
textField(x=-0.4,y=70.5,"B",cex=2)
textField(x=-0.4,y=95.5,"C",cex=2)
textField(x=-0.4,y=125.5,"D",cex=2)
axis(1,at=c(1,5,10,15,20,25),labels = c(2020,2025,2030,2035,2040,2045))
legend("bottomleft",c("UBA Target 40","UBA Target 70"),lty=c(3,2),lwd=c(2,2))
dev.off()

#### evolution threshold nonres 

png("./figures/fig9.png",width=700,height=500)
par(cex.axis=2, cex.lab=1.5, cex.main=1.5, cex.sub=2,mar=c(5,5,6,3))

col_E     <- rgb(col2rgb("orange")[1]/255,col2rgb("orange")[2]/255,col2rgb("orange")[3]/255,0.5) 
col_F     <- rgb(col2rgb("coral1")[1]/255,col2rgb("coral1")[2]/255,col2rgb("coral1")[3]/255,0.5) 
col_G     <- rgb(col2rgb("red")[1]/255,col2rgb("red")[2]/255,col2rgb("red")[3]/255,0.5) 

plot(NA,ylim=c(0,150),xaxt="n",xlim=c(0,25),xlab="Year",ylab="Final energy demand [kwh/(m^2*a)]",
     main = "Evolution of the 15% quantiles under different scenarios (non-residential)")

polygon(x=c(-5,-5,40,40),y=c(-5,30,30,-5),col=col_Aplus,border = F)  #A+
polygon(x=c(-5,-5,40,40),y=c(30,50,50,30),col=col_A,border = F)  #A+
polygon(x=c(-5,-5,40,40),y=c(50,75,75,50),col=col_B,border = F)  #A+
polygon(x=c(-5,-5,40,40),y=c(75,100,100,75),col=col_C,border = F)  #A+
polygon(x=c(-5,-5,40,40),y=c(100,130,130,100),col=col_D,border = F)  #A+
polygon(x=c(-5,-5,40,40),y=c(130,160,160,130),col=col_E,border = F)  #A+
polygon(x=c(-5,-5,40,40),y=c(160,200,200,160),col=col_F,border = F)  #A+
polygon(x=c(-5,-5,40,40),y=c(200,250,250,200),col=col_G,border = F)  #A+

lines(q_vec_zi_25,lty=3,lwd=2)
lines(q_vec_zi_45,lty=2,lwd=2)

textField(x=0,y=23.5,"A+",cex=2)
textField(x=-0.26,y=43.3,"A",cex=2)
textField(x=-0.26,y=68.5,"B",cex=2)
textField(x=-0.26,y=93.5,"C",cex=2)
textField(x=-0.26,y=123.5,"D",cex=2)
textField(x=-0.26,y=145,"E",cex=2)
textField(x=-0.26,y=199,"F",cex=2)
textField(x=-0.26,y=245,"G",cex=2)

axis(1,at=c(1,5,10,15,20,25,30),labels = c(2020,2025,2030,2035,2040,2045,2050))
legend("bottomleft",c("UBA Target 25","UBA Target 45"),lty=c(3,2),lwd=c(2,2))
dev.off()

round(q_res_2030_zi_40)
round(q_res_2040_zi_40)

round(q_res_2030_zi_70)
round(q_res_2040_zi_70)

################################
### L A T E X   T A B L E S ####
################################

tab <- data.frame("Year"=2021:2045,
           "UBA Target 40"=q_vec_zi_40,
           "UBA Target 70"=q_vec_zi_70,
           "UBA Target 25"=q_vec_zi_25,
           "UBA Target 45"=q_vec_zi_45)

print(xtable(tab), include.rownames=FALSE)
#save.image("./data/15perc.Rdata")
