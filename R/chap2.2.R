# GA for Travelling salesman problem (TSP) 
# https://qiita.com/yama1223xxx/items/b609086de919e3af34d3
# http://www.obitko.com/tutorials/genetic-algorithms/japanese/index.php

rm(list = ls())

require(dplyr)
require(magrittr)
require(foreach)
require(pforeach)

# location of cities and travel routes --------------------------------------
set.seed(7)
setCities <- function(n.cities, mesh = seq(0,1,1e-4)){
  cities <- data.frame(id = 1:n.cities, 
                       x  = sample(mesh, n.cities), 
                       y  = sample(mesh, n.cities))
  cities[1, ] <- data.frame(1,0,0)
  invisible(cities)
}
(cities <- setCities(10))

# start from id=1
travel <- data.frame(id = c(1, sample(2:NROW(cities)))) %>% 
  left_join(cities, by="id") %>% 
  slice(c(1:NROW(.), 1))
travel

## calculate total trip. 
total.distance <- foreach(i=1:NROW(cities), .combine = sum) %do%({
  from <- travel[i, ]
  to   <-  travel[i+1, ]
  
  distance <- ((from$x - to$x)^2 + (from$y - to$y)^2)^0.5
  return(distance)
})
total.distance

# plot travelling orders
plot(y~x, travel[-1,], cex=1.5, 
     main = sprintf("total trip = %f", total.distance))
points(y~x, data=travel[1,], pch=16, col="red", cex=1.5)
for(i in 1:NROW(cities)){
  
  from <- travel[i, ]
  to <-  travel[i+1, ]
  arrows(x0 = from$x, y0 = from$y, x1 = to$x, y1=to$y, col="blue", 
         length = 0.1, angle = 20)
  
}

# individual --------------------------------------------------------------

individual <- function(.cities){
  c(1, sample(2:NROW(.cities)))
}
# example
x1 <- individual(cities)
paste0(x1, collapse = "-")


# get total distance
fitness <- function(trip, cities){
  stopifnot(length(trip) == NROW(cities))
  
  travel <- data.frame(id = trip) %>% 
    left_join(cities, by="id") %>% 
    slice(c(1:NROW(.), 1))
  
  total.distance <- foreach(i=1:NROW(cities), .combine = sum) %do%({
    from <- travel[i, ]
    to   <- travel[i+1, ]
    
    distance <- ((from$x - to$x)^2 + (from$y - to$y)^2)^0.5
    return(distance)
  })
  
  return(total.distance)
}
# example
x1 <- individual(cities)
fitness(x1, cities)


# order crossover
crossover <- function(p1, p2, show.pos = FALSE){
  len_chrom <- length(p1)
  stopifnot(len_chrom >= 3,
            len_chrom == length(p2))
  
  # position at crossover
  at <- sample(2:(length(p1)-1),2)
  at <- sort(at)
  
  child <- rep(NA, length(p1))
  child[1] <- 1
  child[at[1]:at[2]] <- p1[at[1]:at[2]]
  chid.inherit <- child
  
  p2.rot <- p2[c((at[2]+1):length(p2),1:(at[2]))]
  child.omit <- setdiff(p2, child)
  child[which(is.na(child))] <- child.omit
  
  if(show.pos){
    child <- list(at = at, p1 = p1, p2 = p2, 
                  p1.inherit = chid.inherit, p2.rotated = p2.rot,
                  p2.omitted = child.omit, new.chrom = child)
  }
  return(child)
}
# example
chr1 <- individual(cities)
chr2 <- individual(cities)
crossover(chr1, chr2)
crossover(chr1, chr2, show.pos = TRUE)


# mutation (inversion)
mutation <- function(chrom, mutate.prob =1.0, show.pos=FALSE){
  stopifnot(length(chrom)>=3)
  if(runif(1) > mutate.prob){
    return(chrom)
  }
  
  at <- sample(2:length(chrom), 2)
  at <- sort(at)
  chrom[at] <- chrom[rev(at)]
  if(show.pos == TRUE){
    print(sprintf("inversion between %i <=> %i", at[1], at[2]))
  }
  return(chrom)
}
# example
mutation(1:10, show.pos = TRUE)
mutation(1:3, show.pos = TRUE) # must be [1] 1 3 2


# population --------------------------------------------------------------
printChroms <- function(l){
  sapply(l, paste0, collapse = "-") %>% tibble()
}
# 
# POP_SIZE    = 10  # population size
# N_ELITE     = 2   # number of elite individual for next chromration
# MUTATE_PROB = 0.5 # mutation rate
# 
# print(factorial(NROW(cities))) > POP_SIZE
# 
# chrom <- NULL
# while(length(chrom) < POP_SIZE){
#   new.ind <- individual(cities)
#   # print(length(chrom))
#   germ <- paste0(new.ind, collapse = "-")
#   if(! germ %in% chrom){
#     chrom <- c(chrom, list(new.ind))
#   }
# }
# 
# population <- tibble(chrom = chrom)
# printChroms(population$chrom)
# 
# population$fits  <- sapply(population$chrom, fitness, cities)
# population %<>% arrange(fits)
# 
# # Preserve elite individuals
# nextInd.elite <- population %>% head(N_ELITE)
# nextInd.elite %>% head %>% print()
# 
# # Generate children via: selection -> crossover -> mutation
# # 
# pf <- population$fits
# roulette.pie <-  (max(pf) -pf) / (max(pf) - min(pf))
# roulette.pie <- roulette.pie / sum(roulette.pie)
# 
# roulette.pie %>% tibble()
# pie(roulette.pie, clockwise = TRUE,
#     main = "roulette based on fitness")
# 
# roulette.thr <- cumsum(roulette.pie)
# roulette.thr %>% tibble()
# 
# 
# which(roulette.thr > print(x <- runif(1)))[1]
# which(roulette.thr > 0.999)[1]
# #
# 
# # selection -> crossover
# chrom <- NULL
# while(length(chrom) < NROW(population)- NROW(nextInd.elite)){
#   p1 <- which(roulette.thr > runif(1))[1]
#   p2 <- which(roulette.thr > runif(1))[1]
#   while(p1 == p2){
#     p2 <- which(roulette.thr > runif(1))[1]
#   }
#   
#   new.ind <- crossover(population$chrom[p1] %>% unlist,
#                        population$chrom[p2] %>% unlist)
#   
#   germ <- paste0(new.ind, collapse = "-")
#   if(! germ %in% chrom){
#     chrom <- c(chrom, list(new.ind))
#   }
# }
# printChroms(chrom)
# 
# 
# # mutation
# nextInd.gen <- tibble(
#   chrom = lapply(chrom, mutation, mutate.prob = MUTATE_PROB))
# 
# all.equal(chrom, nextInd.gen$chrom)
# cbind(chrom = printChroms(chrom), nextInd = printChroms(nextInd.gen$chrom))
# 
# population <- bind_rows(nextInd.elite, nextInd.gen)
# population$fits  <- sapply(population$chrom, fitness, cities)
# population %<>% arrange(fits)
# 
# population
# 

# popuration function -----------------------------------------------------

initPopulation <- function(pop.size, .cities){
  stopifnot(factorial(NROW(.cities)) > pop.size)
  
  chrom <- NULL
  while(length(chrom) < pop.size){
    new.ind <- individual(.cities)
    # print(length(chrom))
    germ <- paste0(new.ind, collapse = "-")
    if(! germ %in% chrom){
      chrom <- c(chrom, list(new.ind))
    }
  }
  population <- tibble(chrom = chrom)
  population$fits  <- sapply(population$chrom, fitness, .cities)
  
  return(population %>% arrange(fits))
}
# example
# pop <- initPopulation(20, cities)
# printChroms(pop$chrom)

alternate <- function(population,  .cities, elite.size, mutate.prob = 0.3){
  stopifnot(!missing(population), !missing(.cities), 
            NROW(population) > elite.size )
  
  # Preserve elite individuals
  nextInd.elite <- population %>% 
    arrange(fits) %>% 
    head(elite.size)

  # Generate children via: selection -> crossover -> mutation
  # roulette source based on fitness for selection
  pf <- population$fits
  roulette.pie <-  (max(pf) -pf) / (max(pf) - min(pf)) 
  roulette.pie <- roulette.pie / sum(roulette.pie)
  roulette.thr <- cumsum(roulette.pie)
  
  # selection -> crossover
  chrom <- NULL
  while(length(chrom) < NROW(population) - NROW(nextInd.elite)){
    p1 <- which(roulette.thr > runif(1))[1]
    p2 <- which(roulette.thr > runif(1))[1]
    while(p1 == p2){
      p2 <- which(roulette.thr > runif(1))[1]
    }
    
    new.ind <- crossover(population$chrom[p1] %>% unlist,
                         population$chrom[p2] %>% unlist)
    
    germ <- paste0(new.ind, collapse = "-")
    if(! germ %in% chrom){
      chrom <- c(chrom, list(new.ind))
    }
  }

  # mutation
  nextInd.gen <- tibble(
    chrom = lapply(chrom, mutation, mutate.prob = mutate.prob))

  population <- bind_rows(nextInd.elite, nextInd.gen)
  population$fits  <- sapply(population$chrom, fitness, cities)
  
  return(population %>% arrange(fits))
}
# example
# (pop <- initPopulation(30, cities))
# (pop.new <- alternate(pop, cities, elite.size = 3, mutate.prob = 0.8))


plot.trip <- function(trip, .cities){
  stopifnot(!missing(trip),
            !missing(.cities), 
            NROW(.cities) == length(trip))
  travel <- data.frame(id = trip) %>% 
    left_join(.cities, by="id") %>% 
    slice(c(1:NROW(.), 1))
  # print(trip)
  # print(travel)
  
  total.distance <- fitness(trip, .cities)
  # print(total.distance)
  
  # plot travelling orders
  plot(y~x, travel[-1,], cex=1.5, 
       main = sprintf("total trip = %f", total.distance))
  points(y~x, data=travel[1, ], pch=16, col="red", cex=1.5)
  for(i in 1:NROW(.cities)){
    from <- travel[i, ]
    to   <- travel[i+1, ]
    arrows(x0 = from$x, y0 = from$y, x1 = to$x, y1=to$y, col="blue", 
           length = 0.1, angle = 20)
  }
}
# example
# sample.trip <- individual(cities)
# fitness(sample.trip, cities)
# plot.trip(sample.trip, cities)

# exec --------------------------------------------------------------------
set.seed(7)
N_CITIES    = 25  # number of cities to travel
GEN_MAX     = 600  # number of generation
POP_SIZE    = 50  # population size
N_ELITE     = 5   # number of elite individual for next chromration
MUTATE_PROB = 1.0 # mutation rate

(cities <- setCities(N_CITIES))
sample.trip <- individual(cities)
plot.trip(sample.trip, cities)


start_time <- Sys.time()
generation <- list(NULL)

(pop <- initPopulation(POP_SIZE, cities))
for(i in 1:GEN_MAX){
  print(i)
  generation[[i]] <- pop
  pop <- alternate(pop,  cities, 
                   elite.size  = N_ELITE, 
                   mutate.prob = MUTATE_PROB)
}

generation[[i]]$chrom[1]
top1 <- NULL
for(i in 1:length(generation)){
  this <- generation[[i]]
  top1 <- rbind(top1, tibble(gen=i, 
                             chrom = this$chrom[1], 
                             fits  = this$fits[1]))
}
top1

Sys.time() - start_time


# plot & print --------------------------------------------------------------
library("animation")
saveGIF({
  
  for(g in 1:length(generation)){
    par(mfrow = c(1,2))
    
    this.trip <- unlist(top1$chrom[g])
    plot.trip(this.trip, cities)
    
    plot(fits~gen, top1, type="b",
         main = sprintf("generation = %i (fits = %f)", g, top1$fits[g]))
    points(x=g, y=top1$fits[g], pch=16, col="red", cex=1.5)
    # top1;max(Y)
    
    par(mfrow = c(1,1))
  }
}, interval = 0.4, movie.name = "stepGA.gif", ani.width=640, ani.height=320)


print(top1)


