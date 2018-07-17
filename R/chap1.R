require(dplyr)
require(magrittr)
require(foreach)

GEN_MAX     = 1000 # number of generation
POP_SIZE    = 1000 # population size
ELITE       = 10   # number of elite individual for next chromration
MUTATE_PROB = 0.01 # mutation rate
N           = 64   # size of problem (= number of element)


# individual --------------------------------------------------------------

individual <- function(N){
  stopifnot(N > 1)
  
  chrom <- sample(c(0,1), N, replace=TRUE)
  return(chrom)  
}
# example
(ind1 <- individual(64))


fitness <- function(chrom){
  stopifnot(all(chrom %in% c(0,1)) == TRUE)
  return(
    sum((chrom *2 -1) * seq_along(chrom)^0.5) %>% abs
  )
}
# example
(ind1 <- individual(64))
fitness(ind1)


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
c3 <- rep(1, 100)
mutation(c3, 0.05)

# population --------------------------------------------------------------

initPopulation <- function(pop.size, chrom.size){
  stopifnot(!missing(pop.size), !missing(chrom.size))
  
  population <- tibble(
    chrom = foreach(i = 1:pop.size) %do% individual(chrom.size))
  
  population$fits <- sapply(population$chrom, fitness)
  
  return(population %>% arrange(fits))
}

# 
# GEN_MAX     = 100  # number of generation
# POP_SIZE    = 100  # population size
# ELITE       = 10   # number of elite individual for next chromration
# MUTATE_PROB = 0.01 # mutation rate
# N           = 64   # size of problem (= number of element)
# 
# population <- initPopulation(pop.size = POP_SIZE, chrom.size = N)
# population
# 
# pop.size = POP_SIZE
# # Preserve elite individuals
# nextInd.pres <- population %>% head(ELITE)
# 
# # Generate children via: selection -> crossover -> mutation
# # roulette source based on rank for selection
# x <- 1:NROW(population)
# denom <- foreach(i = x, .combine = c) %do% rep(x[i], rev(x)[i])
# 
# # selection -> crossover
# nextChrom <- foreach(i = (ELITE+1):pop.size) %do% {
#   p1 <- sample(denom, 1)
#   p2 <- sample(denom, 1)
#   while(p1 == p2){
#     p2 <- sample(denom, 1)
#     # print(p2)  
#   }
#   crossover(population$chrom[p1] %>% unlist,
#             population$chrom[p2] %>% unlist)
# }
# 
# # mutation
# nextInd.gen <- tibble(
#   chrom = lapply(nextChrom, mutation, mutate.prob = MUTATE_PROB))
# nextInd.gen$fits <- sapply(nextInd.gen$chrom, fitness)
# 
# population <- bind_rows(nextInd.pres, nextInd.gen) %>% 
#   arrange(fits) 
# 
# population %>% head()

# popuration function -----------------------------------------------------

alternate <- function(population, elite.size, mutate.prob = 0.01){
  stopifnot(!missing(population), NROW(population) > elite.size )
  
  # Preserve elite individuals
  nextInd.pres <- population %>%
    arrange(fits) %>% 
    head(elite.size)
  # nextInd.pres %>% head %>% print()
  
  # Generate children via: selection -> crossover -> mutation
  # roulette source based on rank for selection
  x <- 1:NROW(population)
  denom <- foreach(i = x, .combine = c) %do% rep(x[i], rev(x)[i])
  
  # selection -> crossover
  nextChrom <- foreach(i = (elite.size+1):NROW(population)) %do% {
    p1 <- sample(denom, 1)
    p2 <- sample(denom, 1)
    while(p1 == p2){
      p2 <- sample(denom, 1)
      # print(p2)  
    }
    crossover(population$chrom[p1] %>% unlist,
              population$chrom[p2] %>% unlist)
  }
  
  # mutation
  nextInd.gen <- tibble(
    chrom = lapply(nextChrom, mutation, mutate.prob = mutate.prob))
  nextInd.gen$fits <- sapply(nextInd.gen$chrom, fitness)
  
  return( rbind(nextInd.pres, nextInd.gen) )
}

(pop <- initPopulation(pop.size = 150, chrom.size = 500))
topTeam <- pop %>% head(1)

for(i in 2:50){
  nextpop <- alternate(pop, elite.size = 10, mutate.prob = 0.05)
  topTeam[i, ] <- pop %>% head(1)

  pop <- nextpop
}

plot(x=1:length(topTeam$fits), y=topTeam$fits,
     type="b", xlab = "generation", ylab = "fitness")

