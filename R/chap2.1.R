# get argmax(f(x)) 
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
    4*sin(25*x-1) + 6*sin(1- 10*x)
  )
}

X <- seq(0,1,0.001); Y=fitnessCurve(X)
plot(x=X, y=Y, type="l")
points(x=X[which.max(Y)], max(Y), col="red")

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


# fitness <- function(chrom){
#   stopifnot(all(chrom %in% c(0,1)) == TRUE)
#   return(
#     sum((chrom *2 -1) * seq_along(chrom)^0.5) %>% abs
#   )
# }
# # example
# (ind1 <- individual(64))
# fitness(ind1)


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
# 
# population <- initPopulation(pop.size = POP_SIZE, chrom.size = N)
# population
# 
# pop.size = POP_SIZE
# # Preserve elite individuals
# nextInd.elite <- population %>% head(ELITE)
# nextInd.elite %>% head %>% print()
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
# population <- bind_rows(nextInd.elite, nextInd.gen) %>% 
#   arrange(fits) 
# 
# population %>% head()

# popuration function -----------------------------------------------------

initPopulation <- function(pop.size, chrom.size){
  stopifnot(!missing(pop.size), !missing(chrom.size))
  
  population <- tibble(
    chrom = foreach(i = 1:pop.size) %do% individual(chrom.size))
  
  population$pheno <- sapply(population$chrom, bin2dec)
  population$fits  <- sapply(population$pheno, fitnessCurve)
  
  return(population %>% arrange(desc(fits)))
}


alternate <- function(population, elite.size, mutate.prob = 0.01){
  stopifnot(!missing(population), 
            NROW(population) > elite.size )
  
  # Preserve elite individuals
  nextInd.elite <- population %>%
    arrange(desc(fits)) %>% 
    head(elite.size)
  
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
    }
    crossover(population$chrom[p1] %>% unlist,
              population$chrom[p2] %>% unlist)
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


GEN_MAX     = 30  # number of generation
POP_SIZE    = 20  # population size
CHROM_SIZE  = 64   # size of problem
N_ELITE     = 5   # number of elite individual for next chromration
MUTATE_PROB = 0.05 # mutation rate

generation <- list(NULL)
pop <- initPopulation(pop.size = POP_SIZE, chrom.size = CHROM_SIZE)

for(i in 1:GEN_MAX){
  generation[[i]] <- pop
  pop <- alternate(pop, elite.size = N_ELITE, mutate.prob = MUTATE_PROB)
}


X <- seq(0,1,1e-5); Y=fitnessCurve(X)
plot(x=X, y=Y, type="l")
points(x=X[which.max(Y)], max(Y), col="red")

bin2dec(rep(1,64), norm=F)

# generation %>% str(2)

X <- seq(0,1,1e-4); Y=fitnessCurve(X)


library("animation")

g=2

par(mfrow = c(1,2))
top5 <- generation[[g]]
plot(x=X, y=Y, type="l")
points(x=X[which.max(Y)], max(Y), col="red", pch = 4)

points(x=top5$pheno, y=top5$fits, col=6, pch=1)
points(x=top5$pheno[1], y=top5$fits[1], col=6, pch=16)



top1 <- NULL
for(i in 1:length(generation)){
  this <- generation[[i]]
  top1 <- rbind(top1, data.frame(gen=i, fits=this$fits[1]))
}

plot(fits~gen, top1, type="b",ylim=c(7.5,10))
abline(h=max(Y), lty = 2, col=3)
points(x=g, y=top1$fits[g], col=2, pch=16)
# top1;max(Y)

par(mfrow = c(1,1))
Sys.time() - start_time


saveGIF({
  #ŒJ‚è•Ô‚µ‰ñ”‚ðÝ’è‚·‚é 
  for (i in 1:5){
    #ˆ—“à—e‚ð‹Lq‚·‚é
    plot(runif(10), ylim = c(0,1))
  }
}, interval = 1.0, movie.name = "TEST.gif")


