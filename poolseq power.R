setwd("/Users/jannawilloughby/GDrive/gray bats - alabama/poolseq_quant//")
library(scales)

popsize = 2000 #total number of seqs in pool seq
seqsize = 50   #total number of seqs with individual genotype data 
nSNPs   = 100  #number of SNPs

#initiate pop seq matrix
pop = matrix(data=NA, nrow=popsize, ncol=(nSNPs+1))
pop[ ,1] = seq(1, popsize, 1) #individual ids

#initiate indv seq matrix
ind = matrix(data=NA, nrow=seqsize, ncol=(nSNPs+1))
ind[ ,1] = seq(1, seqsize, 1) #individual ids

#function to create pools of genotypes in HWE
geno.pool = function(popsize, seqsize){
  totpopsize = popsize + seqsize
  p = sample(seq(from = 0, to = 1, by = 0.001), 1)
  genos = c(rep(2, round((totpopsize*p*p), 0)),                                                #homozygous p*p
            rep(1, round((totpopsize*(p)*(1-p)), 0)),                                          #homozygous 1-p^2
            rep(0, totpopsize-(round((totpopsize*p*p), 0)+round((totpopsize*(1-p)*(1-p)), 0)))       #heterozygotes 2pq
  )
  #genos = sample(pool, popsize, replace=F)
  return(list(genos, p))
}

#pop[,2:(nSNPs+1)] = apply(pop[,2:(nSNPs+1)], 2, geno.pool, popsize=popsize)

#get genotypes; since above won't work I'll loop it but I'm not happy about it
pfreq = NULL
for(c in 2:(nSNPs+1)){
  all.output = unlist(geno.pool(popsize, seqsize))
  pop[,c] = all.output[1:popsize]
  ind[,c] = all.output[(1+popsize):(1+popsize+seqsize)]
  pfreq = c(pfreq, all.output[length(all.output)])
}

#estimate allele freqs for "pool seq" data
estimateAF = apply(pop[,2:ncol(pop)], 2, sum, na.rm=T)/(popsize*2)
plot(pfreq, estimateAF, xlim=c(0,1), ylim=c(0,1))
segments(0,0,1,1)

  


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
