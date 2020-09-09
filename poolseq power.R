setwd("/Users/jannawilloughby/GDrive/gray bats - alabama/poolseq_quant//")
library(scales)

popsize = 200
nSNPs   = 100

#initiate population
pop = matrix(data=NA, nrow=popsize, ncol=(nSNPs+1))
pop[ ,1] = seq(1, popsize, 1) #individual ids

#function to create pools of genotypes in HWE
geno.pool = function(popsize){
  p = sample(seq(from = 0, to = 1, by = 0.01), 1)
  pool = c(rep(0, round((popsize*p*p), 0)),                                                #homozygous p*p
           rep(1, round((popsize*(1-p)*(1-p)), 0)),                                        #homozygous 1-p^2
           rep(2, popsize-(round((popsize*p*p), 0)+round((popsize*(1-p)*(1-p)), 0)))       #heterozygotes 2pq
  )
  genos = sample(pool, popsize, replace=T)
  return(genos)
}

#pop[,2:(nSNPs+1)] = apply(pop[,2:(nSNPs+1)], 2, geno.pool, popsize=popsize)

#get genotypes; since above won't work I'll loop it but I'm not happy about it
for(c in 2:(nSNPs+1)){
  pop[,c] =geno.pool(popsize)
}



rep(geno.pool, popsize)

pop = apply(pop, 1, geno.pool, popsize=popsize)

pop = lapply(pop, geno.pool, popsize=popsize)







for(n in 1:nSNPs){
  
  
 

  gtype = sample(pool, popsize, replace=FALSE)
  for(i in 1:popsize){
    if(gtype[i]==0){                            #homozygous (0,0)
      genos[i,1]   = 0
      genos[i,2] = 0
      next
    }else if(gtype[i]==1){                      #heterozygous (0,1)
      genos[i,1]   = 0
      genos[i,2] = 1
    }else{                                      #homozygous (1,1)
      genos[i,1]   = 1
      genos[i,2] = 1
    }
  }
  pop[,(n+1)]  = apply(genos, 1, sum)
  pool = genos = NULL
}  



  
#allele 1 positions
positions = seq(1, (nSNPs*2), 2)
  
  # randomly sample either position 1 or 2 (add 0 or 1) to starting position
  fallele  <- positions + sample(0:1, nloci * noff, replace = TRUE)
  fallele2 <- fg[fallele]
  fallele3 <- matrix(fallele2, nrow = noff, ncol = nloci, byrow = TRUE)
  
  mallele  <- positions + sample(0:1, nloci * noff, replace = TRUE)
  mallele2 <- mg[mallele]
  mallele3 <- matrix(mallele2, nrow = noff, ncol = nloci, byrow = TRUE)
  
  offspringG[, positions]     <- fallele3
  offspringG[, positions + 1] <- mallele3
  
  offspring    = cbind(offspring, offspringG)
  gstartID     = gstartID + nrow(offspring)
  
  
  
}




#output dataframe
OUT = NULL

#iterate over number of collected samples (marks and recaptures)
for(c in 1:length(collects)){
  collect = collects[c]
  recollect = recollects[c]
  Nst = NULL
  #repeat each marking effort 100 times
  for(i in 1:1000){
    #create initial population
    pop1 = data.frame(uid = seq(1, popsize, 1), capture = rep(0, popsize), recapture = rep(0, popsize))
    #initial marking
    pop1$capture[(sample(x=c(1:popsize), size=collect, replace=F))] = 1
    #recapture
    pop1$recapture[(sample(x=c(1:popsize), size=recollect, replace=F))] = 1
    #estimate pop1 size and add to list
    Nst = c(Nst, ((sum(pop1$capture)+1)*(collect+1))/((sum(pop1$recapture[pop1$capture==1])+1)))
  }
  #record mean/sd for pop1 size estimates
  writeout = c(collect, recollect, mean(Nst), sd(Nst), quantile(Nst, probs=0.975), quantile(Nst, probs=0.025))
  #add these values to others
  OUT = rbind(OUT, writeout)
}
colnames(OUT) = c("collect", "recollect", "NstM", "NstSD", "NstUL", "NstLL")
rownames(OUT) = seq(1,nrow(OUT), 1)
o = as.data.frame(OUT)
o$NstSE = 1.96*o$NstSD/sqrt(i)
o$totalgenos = o$collect + o$recollect








#plot data nicely
off = 100
plot(-10000, -10000, xlim=c(0, 100000), ylim=c(0,3000000), xlab="total number of samples genotyped per year", ylab="estimated population size")
segments(x0=0, x1=100000, y0=popsize, y1=popsize, lty=2, col="grey50")
lines(x=o$totalgenos, y=o$NstM, lty=1, col=alpha("firebrick3", 0.5), lwd=2)
points(x=o$totalgenos, y=o$NstM, pch=19, col=alpha("firebrick3", 0.5), cex=1.5)
segments(x0=c(o$totalgenos - off), x1=c(o$totalgenos + off), y0=o$NstLL, y1=o$NstLL, col=alpha("firebrick3", 0.8))
segments(x0=c(o$totalgenos - off), x1=c(o$totalgenos + off), y0=o$NstUL, y1=o$NstUL, col=alpha("firebrick3", 0.8))
segments(x0=o$totalgenos, x1=o$totalgenos, y0=o$NstLL, y1=o$NstUL, col=alpha("firebrick3", 0.8))

write.table(o, "markrecapture_samp/simulationoutput.csv", row.names=F, col.names=T, sep=",")
