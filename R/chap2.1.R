# GA for getting argmax(f(x)) 
# http://www.obitko.com/tutorials/genetic-algorithms/japanese/index.php

require(dplyr)
require(magrittr)
require(foreach)

# individual --------------------------------------------------------------

individual <- function(N){
  sample(c(0,1), N, replace=TRUE)
}
# example
(ind1 <- individual(64))

fitnessCurve <- function(x){
  return(
    2*sin(29*x-1) + 5*sin(1- 8*x) + 3*sin(1- 70*x)
  )
}
# example


X <- seq(0,1,0.001)
Y <- fitnessCurve(X)
plot(x=X, y=Y, type="l")

bin2dec <- function(chrom, norm = TRUE){
  stopifnot(all(chrom %in% c(0,1)) == TRUE)
  dec <- sum(2^(1:length(chrom) -1) * chrom)
  if(norm == TRUE){
    dec <- dec / sum(2^(1:length(chrom) -1))
  }
  
  return(dec)
}
# example
(ind1 <- individual(64))
fitnessCurve( bin2dec(ind1))


crossover <- function(p1, p2, show.pos = FALSE){
  len_chrom <- length(p1)
  stopifnot(len_chrom == length(p2))
  
  # position at crossover
  at <- sample(len_chrom - 1 , 1)
  
  child <- c(p1[1:at], p2[(at+1):len_chrom])
  if(show.pos){
    child <- list(at = at, p1 = p1, p2 = p2, chrom = child)
  }
  return(child)
}
# example
(chr1 <- rep(1, 10))
(chr2 <- rep(0, 10))
crossover(chr1, chr2, show.pos = TRUE)
crossover(chr1, chr2, show.pos = FALSE)


mutation <- function(chrom, mutate.prob){
  stopifnot(all(chrom %in% c(0,1)) == TRUE,
            !missing(mutate.prob))
  
  # position at mutation
  at <- which(runif(length(chrom)) < mutate.prob)
  chrom[at] <- 1 - chrom[at]
  return(chrom)
}
# example
chr3 <- rep(1, 100)
mutation(chr3, 0.05)

# population --------------------------------------------------------------
# 
# GEN_MAX     = 15  # number of generation
# POP_SIZE    = 10  # population size
# CHROM_SIZE  = 64   # size of problem
# N_ELITE     = 5   # number of elite individual for next chromration
# N_TOURNAMEST= 5   # number of individual for tournament selection
# MUTATE_PROB = 0.05 # mutation rate
# 
# population <- tibble(
# chrom = foreach(i = 1:POP_SIZE) %do% individual(CHROM_SIZE))
# 
# population$pheno <- sapply(population$chrom, bin2dec)
# population$fits  <- sapply(population$pheno, fitnessCurve)
# 
# population %<>% arrange(desc(fits))
# 
# # Preserve elite individuals
# nextInd.elite <- population %>% head(N_ELITE)
# nextInd.elite %>% head %>% print()
# 
# # Generate new children
# # tournament selection -> crossover
# nextChrom <- foreach(i = (N_ELITE+1):NROW(population)) %do% {
#   parents <- population %>% 
#     sample_n(N_TOURNAMEST) %>% 
#     arrange(desc(fits)) %>% 
#     head(2)
#   crossover(parents$chrom[1] %>% unlist,
#             parents$chrom[2] %>% unlist)
# }
# 
# # mutation
# nextInd.gen <- tibble(
#   chrom = lapply(nextChrom, mutation, mutate.prob = MUTATE_PROB))
# 
# nextInd.gen$pheno <- sapply(nextInd.gen$chrom, bin2dec)
# nextInd.gen$fits  <- sapply(nextInd.gen$pheno, fitnessCurve)
# 
# nextGen.all <- rbind(nextInd.elite, nextInd.gen) %>% 
#   arrange(desc(fits))
# 
# nextGen.all

# popuration function -----------------------------------------------------

initPopulation <- function(pop.size, chrom.size){
  stopifnot(!missing(pop.size), !missing(chrom.size))
  
  population <- tibble(
    chrom = foreach(i = 1:pop.size) %do% individual(chrom.size))
  
  population$pheno <- sapply(population$chrom, bin2dec)
  population$fits  <- sapply(population$pheno, fitnessCurve)
  
  return(population %>% arrange(desc(fits)))
}


alternate <- function(population, elite.size, tournament.size, mutate.prob = 0.01){
  stopifnot(!missing(population), 
            NROW(population) > elite.size )
  
  # Preserve elite individuals
  nextInd.elite <- population %>%
    arrange(desc(fits)) %>% 
    head(elite.size)
  
  # Generate new children
  # tournament selection -> crossover
  nextChrom <- foreach(i = (elite.size+1):NROW(population)) %do% {
    parents <- population %>% 
      sample_n(tournament.size) %>% 
      arrange(desc(fits)) %>% 
      head(2)
    crossover(parents$chrom[1] %>% unlist,
              parents$chrom[2] %>% unlist)
  }
  
  # mutation
  nextInd.gen <- tibble(
    chrom = lapply(nextChrom, mutation, mutate.prob = mutate.prob))
  
  nextInd.gen$pheno <- sapply(nextInd.gen$chrom, bin2dec)
  nextInd.gen$fits  <- sapply(nextInd.gen$pheno, fitnessCurve)
  
  nextGen.all <- rbind(nextInd.elite, nextInd.gen) %>% arrange(desc(fits))
  return(nextGen.all)
}


# exec --------------------------------------------------------------------
set.seed(6)
start_time <- Sys.time()

GEN_MAX     = 20  # number of generation
POP_SIZE    = 100  # population size
CHROM_SIZE  = 16   # size of problem
N_ELITE     = 5   # number of elite individual for next chromration
N_TOURNAMEST= 5   # number of individual for tournament selection
MUTATE_PROB = 0.05 # mutation rate

generation <- list(NULL)
pop <- initPopulation(pop.size = POP_SIZE, chrom.size = CHROM_SIZE)

for(i in 1:GEN_MAX){
  generation[[i]] <- pop
  pop <- alternate(pop, 
                   tournament.size = N_TOURNAMEST,
                   elite.size = N_ELITE, 
                   mutate.prob = MUTATE_PROB)
}

top1 <- NULL
for(i in 1:length(generation)){
  this <- generation[[i]]
  top1 <- rbind(top1, data.frame(gen=i, fits=this$fits[1]))
}

top1$diff <- max(Y) - top1$fits
Sys.time() - start_time


# plot & print --------------------------------------------------------------
library("animation")
X <- seq(0,1,1e-4); Y=fitnessCurve(X)

saveGIF({
  
  for(g in 1:length(generation)){
    par(mfrow = c(1,2))
    
    top5 <- generation[[g]]
    plot(x=X, y=Y, type="l",
         main = sprintf("generation = %i", g))
    points(x=X[which.max(Y)], max(Y), col="red", pch = 4)
    points(x=top5$pheno, y=top5$fits, col=6, pch=1)
    points(x=top5$pheno[1], y=top5$fits[1], col=6, pch=16)

    plot(fits~gen, top1, type="b",
         main = sprintf("generation = %i (diff = %f)",
                        g, top1$diff[g]))
    abline(h=max(Y), col=3)
    points(x=g, y=top1$fits[g], pch=16,
           col=ifelse(top1$diff[g]>0,"green","red"))
    # top1;max(Y)
    
    par(mfrow = c(1,1))
    
  }
}, interval = 1.0, movie.name = "stepGA.gif", ani.width=960, ani.height=480)

print(top1)



