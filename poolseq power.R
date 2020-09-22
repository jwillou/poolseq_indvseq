setwd("/Users/jannawilloughby/GDrive/gray bats - alabama/poolseq_quant/")
library(scales)

popsize = 100                                    #total number of seqs in pool seq
seqsize = 30                                     #total number of seqs with individual genotype data 
pop.Ts  = c(20000, 100000, 500000, 1000000)      #true population size
nSNPs   = 10000                                  #number of SNPs

rarecutoff = 0.10                                #rare genotypes cutoff value
repsperpop.T = 100                               #number of replicates per true population size (pop.Ts)

####simulate sampling####
#set up output files
sink(file="watchrareprop.txt")
OUT = NULL

for(p in 1:length(pop.Ts)){
  pop.T = pop.Ts[p]
  
  for(r in 1:repsperpop.T){
    #initiate pop seq matrix
    pop = matrix(data=NA, nrow=popsize, ncol=(nSNPs+1))
    pop[ ,1] = seq(1, popsize, 1) #individual ids
    
    #initiate indv seq matrix
    ind = matrix(data=NA, nrow=seqsize, ncol=(nSNPs+1))
    ind[ ,1] = seq(1, seqsize, 1) #individual ids
    
    #function to create pools of genotypes in HWE
    geno.pool = function(popsize, seqsize){
      p = sample(seq(from = 0, to = rarecutoff, by = 0.001), 1)
      q = 1 - p
      lgpop = sample(c(2,1,0), pop.T, replace=T, prob=c((p*p), (2*p*q), (q*q)))
      pop.genos = sample(lgpop, popsize, replace=F)
      ind.genos = sample(lgpop, seqsize, replace=F)
      
      #genos = sample(pool, popsize, replace=F)
      return(list(pop.genos, ind.genos, p, q))
    }
    
    #get genotypes; since apply won't work I'll loop it but I'm not happy about it
    pfreq = qfreq = NULL
    for(c in 2:(nSNPs+1)){
      all.output = geno.pool(popsize, seqsize)
      pop[,c] = unlist(all.output[[1]])
      ind[,c] = unlist(all.output[[2]])
      pfreq = c(pfreq, unlist(all.output[[3]]))
      qfreq = c(qfreq, unlist(all.output[[4]]))
    }
    
    #estimate allele freqs for "pool seq" data
    estimateAF = apply(pop[,2:ncol(pop)], 2, sum, na.rm=T)/(popsize*2)
    #plot(pfreq, estimateAF, xlim=c(0,1), ylim=c(0,1))
    #segments(0,0,1,1)
    
    #what is the combined probability of genotypes for individuals, given pop freqs?
    freqraregenos = freqrarealleles = NULL
    for(i in 1:seqsize){
      ii = cbind(ind[i,2:ncol(ind)], estimateAF, rep(NA, ncol(ind)-1))
      #colnames(ii) = c("genotype", "AF", "geno.prob")
      ii[ii[,1]==2, 3] =      ii[ii[,1]==2, 2]  *      ii[ii[,1]==2, 2]
      ii[ii[,1]==1, 3] =  2 * ii[ii[,1]==1, 2]  * (1 - ii[ii[,1]==1, 2])
      ii[ii[,1]==0, 3] = (1 - ii[ii[,1]==0, 2]) * (1 - ii[ii[,1]==0, 2])
      freqraregenos = c(freqraregenos, nrow(ii[ii[,3]<rarecutoff,,drop=F])/nrow(ii))
      freqrarealleles = c(freqrarealleles, nrow(ii[ii[,3]<rarecutoff & ii[,1]>=1,,drop=F])/nrow(ii))
    }
    
    #how often do we find rare alleles and rare genotypes in the genotyped samples?
    freqfindalleles = freqfindgenos = NULL
    for(pp in 1:length(estimateAF)){
      if(estimateAF[pp]<rarecutoff | 1-estimateAF[pp]<rarecutoff){
        if(estimateAF[pp]<rarecutoff){
          freqfindgenos   = c(freqfindgenos, nrow(ind[ind[,(pp+1)]==2,,drop=F])/nrow(ind))
          freqfindalleles = c(freqfindalleles, nrow(ind[ind[,(pp+1)]>=1,,drop=F])/nrow(ind))
          next
        }
        if(estimateAF[pp]<rarecutoff){
          freqfindgenos   = c(freqfindgenos, nrow(ind[ind[,(pp+1)]==2,,drop=F])/nrow(ind))
          freqfindalleles = c(freqfindalleles, nrow(ind[ind[,(pp+1)]>=1,,drop=F])/nrow(ind))
          next
        }
      }
    }
    OUT = rbind(OUT, c(r, pop.T, popsize, seqsize, rarecutoff, mean(freqraregenos), sd(freqraregenos), mean(freqrarealleles), sd(freqrarealleles), mean(freqfindgenos), sd(freqfindgenos), mean(freqfindalleles), sd(freqfindalleles)))
    print(r)
  }
}
colnames(OUT) = c("replicate", "truepopsize", "afpopsize", "genopopsize", "rarefreqcutoff", "raregeno_mean", "raregeno_sd", "rarealleles_mean", "rarealleles_sd", "findgeno_mean", "findgeno_sd", "findalleles_mean", "findalleles_sd")
write.table(OUT, "rareprop.csv", col.names=T, row.names=F, sep=",")
sink()


####analyze output####
#sims = read.table("rareprop.csv", header=T, sep=",")
sims = as.data.frame(OUT)

#rare genotypes
for(pops in pop.Ts){
  hist(sims$raregeno_mean[sims$truepopsize==as.character(pops)], main=paste(pops, mean(sims$raregeno_mean[sims$truepopsize==as.character(pops)]), sep=" "), xlim=c(0.025, 0.05))
}

#rare genotypes
for(pops in pop.Ts){
  hist(sims$findgeno_mean[sims$truepopsize==as.character(pops)], main=paste(pops, mean(sims$findgeno_mean[sims$truepopsize==as.character(pops)]), sep=" "), xlim=c(0.0025, 0.005))
}

#exploratory plots
plot(sims$truepopsize, sims$raregeno_mean)
plot(sims$truepopsize, sims$raregeno_sd)
plot(sims$truepopsize, sims$rarealleles_mean)
plot(sims$truepopsize, sims$rarealleles_sd)
plot(sims$truepopsize, sims$findgeno_mean)
plot(sims$truepopsize, sims$findgeno_sd)
plot(sims$truepopsize, sims$findalleles_mean)
plot(sims$truepopsize, sims$findalleles_sd)



