require(dplyr)
require(foreach)

GEN_MAX     = 1000 # number of chromration
POP_SIZE    = 1000 # population size
ELITE       = 1    # number of elite individual for next chromration
MUTATE_PROB = 0.01 # mutation rate
N           = 64   # size of problem (= number of element)


# individual --------------------------------------------------------------

individual <- function(N){
  stopifnot(N > 1)
  
  chrom <- sample(c(-1,1), N, replace=TRUE)
  return(chrom)  
}
# example
(ind1 <- individual(64))


fitness <- function(chrom){
  stopifnot(all(chrom %in% c(-1,1)) == TRUE)
  return(
    sum(chrom * seq_along(chrom)^0.5) %>% abs
  )
}
# example
(ind1 <- individual(64))
fitness(ind1)


crossover <- function(p1, p2){
  len_chrom <- length(p1)
  stopifnot(len_chrom == length(p2))
  
  # position at crossover
  at <- sample(len_chrom, 1)
  
  child <- list(at = at,
                chrom = c(p1[1:at], p2[(at+1):len_chrom]))
  return(child)
}
# example
(chr1 <- rep(1, 10))
(chr2 <- rep(0, 10))
crossover(chr1, chr2)


mutation <- function(chrom, mutate.prob){
  stopifnot(all(chrom %in% c(-1,1)) == TRUE,
            !missing(mutate.prob))
  
  # position at mutation
  at <- which(runif(length(chrom)) < mutate.prob)
  
  chrom[at] <- -chrom[at]
  return(chrom)
}

c3 <- rep(1, 100)
mutation(c3, 0.05)

# population --------------------------------------------------------------

(population <- foreach(i = 1:10) %do% individual(20))

fits <- sapply(population, fitness)

population2 <- foreach(i = order(fits)) %do% population[[i]]

sort(fits)
lapply(population2, fitness) %>% unlist

# 順位に基づくランキング選択で親個体を2つ選択する
# 戻り値: 選択した親個体の添え字
selectParent <- function(population){
  x <- seq_along(population)
  y <- rev(x)
  
  denom <- foreach(i = x, .combine = c) %do% rep(x[i], y[i])
  return(sample(denom, 1))
}

p1 <- selectParent(population)
p2 <- selectParent(population)
while(p1 == p2){
  p2 <- sample(denom, 1)
  print(p2)  
}
c(p1, p2)

sum(1:5)
