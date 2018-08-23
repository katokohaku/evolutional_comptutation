---
title: "Genetic Algorithm for Travelling salesman problem (TSP) with R"
author: "Satoshi Kato"
date: "2018/08/23"
output:
  html_document:
    keep_md: yes
    fig_caption: yes
    pandoc_args:
    - --from
    - markdown+autolink_bare_uris+tex_math_single_backslash-implicit_figures
    toc: yes
    toc_depth: 3
  word_document:
    toc: yes
    toc_depth: 3
  pdf_document:
    toc: yes
    toc_depth: 3
editor_options: 
  chunk_output_type: inline
---



# problem setting

## Travelling salesman problem (TSP)

![example2](./output/stepGA_20cities.gif)

["Given a list of cities and the distances between each pair of cities, what is the shortest possible route that visits each city and returns to the origin city?" ](https://en.wikipedia.org/wiki/Travelling_salesman_problem)

### location of cities and travel routes


```r
set.seed(7)
setCities <- function(n.cities, mesh = seq(0,1,1e-4)){
  cities <- data.frame(id = 1:n.cities, 
                       x  = sample(mesh, n.cities), 
                       y  = sample(mesh, n.cities))
  cities[1, ] <- data.frame(1,0,0)
  invisible(cities)
}
(cities <- setCities(10))
```

```
##    id      x      y
## 1   1 0.0000 0.0000
## 2   2 0.3977 0.2314
## 3   3 0.1156 0.7727
## 4   4 0.0697 0.0962
## 5   5 0.2436 0.4533
## 6   6 0.7916 0.0846
## 7   7 0.3398 0.5603
## 8   8 0.9714 0.0086
## 9   9 0.1657 0.9850
## 10 10 0.4587 0.3163
```

### represent travel route


```r
# start from id=1
travel <- data.frame(id = c(1, sample(2:NROW(cities)))) %>% 
  left_join(cities, by="id")

travel <- rbind(travel, travel[1, ])
```

### calculate total trip (to be minizimzed)


```r
total.distance <- foreach(i=1:NROW(cities), .combine = sum) %do%({
  from <- travel[i, ]
  to   <-  travel[i+1, ]
  
  distance <- ((from$x - to$x)^2 + (from$y - to$y)^2)^0.5
  return(distance)
})
total.distance
```

```
## [1] 6.012823
```

### sample view: route (travelling orders)


```r
plot(y~x, travel[-1,], cex=1.5, 
     main = sprintf("total trip = %f", total.distance))
points(y~x, data=travel[1,], pch=16, col="red", cex=1.5)
for(i in 1:NROW(cities)){
  
  from <- travel[i, ]
  to <-  travel[i+1, ]
  arrows(x0 = from$x, y0 = from$y, x1 = to$x, y1=to$y, col="blue", 
         length = 0.1, angle = 20)
  
}
```

![](chap2.2_files/figure-html/unnamed-chunk-4-1.png)<!-- -->



# individual

## encodeing

### genetic representation


```r
individual <- function(.cities){
  c(1, sample(2:NROW(.cities)))
}
# example
trip <- individual(cities)
paste0(trip, collapse = "-")
```

```
## [1] "1-5-7-3-8-2-9-4-6-10"
```

### phenotype & fitness


```r
# get total distance
fitness <- function(trip, cities){
  stopifnot(length(trip) == NROW(cities))
  
  travel <- data.frame(id = trip) %>% 
    left_join(cities, by="id")
  
  travel <- rbind(travel, travel[1, ])
  
  total.distance <- foreach(i=1:NROW(cities), .combine = sum) %do%({
    from <- travel[i, ]
    to   <- travel[i+1, ]
    
    distance <- ((from$x - to$x)^2 + (from$y - to$y)^2)^0.5
    return(distance)
  })
  
  return(total.distance)
}
# example
(trip <- individual(cities))
```

```
##  [1]  1  9  5  7  4 10  3  8  2  6
```

```r
fitness(trip, cities)
```

```
## [1] 6.214204
```

## genetic operator

### crossover 

In this case, order-crossover for trip route is applied.



```r
crossover <- function(p1, p2, show.pos = FALSE){
  len_chrom <- length(p1)
  stopifnot(len_chrom >= 3,
            len_chrom == length(p2))
  
  # position at crossover
  at <- sample(2:(length(p1)-1),2)
  at <- sort(at)
  
  child <- rep(NA, length(p1))
  child[at[1]:at[2]] <- p1[at[1]:at[2]]
  chid.inherit <- child
  
  p2.rot <- p2[c((at[2]+1):length(p2),1:(at[2]))]
  p2.omit <- setdiff(p2.rot, child)
  child[which(is.na(child))] <- p2.omit
  
  result <- child
  pos <- which(child == 1)
  if(pos > 1){
    result <- child[c(pos:length(child), 1:(pos-1))]  
  }
  
  if(show.pos){
    result <- list(p1 = p1, p2 = p2, at = at, 
                  p1.inherit = chid.inherit, 
                  p2.rotated = p2.rot, p2.omitted = p2.omit, 
                  new.chrom = child, new.chrom.sorted = result)
  }
  return(result)
}
# example
chr1 <- individual(cities)
chr2 <- individual(cities)
crossover(chr1, chr2)
```

```
##  [1]  1  3 10  4  5  7  6  9  2  8
```

```r
crossover(chr1, chr2, show.pos = TRUE)
```

```
## $p1
##  [1]  1  9 10  8  5  7  6  3  4  2
## 
## $p2
##  [1]  1  3 10  4  9  2  8  5  7  6
## 
## $at
## [1] 4 5
## 
## $p1.inherit
##  [1] NA NA NA  8  5 NA NA NA NA NA
## 
## $p2.rotated
##  [1]  2  8  5  7  6  1  3 10  4  9
## 
## $p2.omitted
## [1]  2  7  6  1  3 10  4  9
## 
## $new.chrom
##  [1]  2  7  6  8  5  1  3 10  4  9
## 
## $new.chrom.sorted
##  [1]  1  3 10  4  9  2  7  6  8  5
```

### mutation

In this case, inversion (node exchange) is applied as mutation.


```r
# mutation (inversion)
mutation <- function(chrom, mutate.prob =0.05, show.pos = FALSE){
  stopifnot(length(chrom)>3)
  
  for(i in 2:length(chrom)){
    if(runif(1) < mutate.prob){
      base <- setdiff(2:length(chrom), i)
      at <- c(i, sample(base, 1))
      if(show.pos == TRUE){
        cat(sprintf("%i <=> %i \n", at[1], at[2]))
      }
      chrom[at] <- chrom[rev(at)]
    }
  }
  return(chrom)
}
# example
mutation(1:10, mutate.prob = 0.3, show.pos = TRUE)
```

```
## 3 <=> 7 
## 6 <=> 7
```

```
##  [1]  1  2  7  4  5  3  6  8  9 10
```


# population

## simple procedure

init 1st generation -> 2nd generation.


```r
printChroms <- function(l){
  sapply(l, paste0, collapse = "-") %>% tibble()
}
```


```r
POP_SIZE    = 10  # population size
N_ELITE     = 2   # number of elite individual for next chromration
MUTATE_PROB = 0.05 # mutation rate

print(factorial(NROW(cities))) > POP_SIZE
```

```
## [1] 3628800
```

```
## [1] TRUE
```

```r
chrom <- NULL
while(length(chrom) < POP_SIZE){
  new.ind <- individual(cities)
  # print(length(chrom))
  germ <- paste0(new.ind, collapse = "-")
  if(! germ %in% chrom){
    chrom <- c(chrom, list(new.ind))
  }
}

population <- tibble(chrom = chrom)
printChroms(population$chrom)
```

```
## # A tibble: 10 x 1
##    .                   
##    <chr>               
##  1 1-9-10-3-6-8-2-5-7-4
##  2 1-5-4-2-7-8-9-3-10-6
##  3 1-9-3-4-8-7-10-5-2-6
##  4 1-2-8-3-6-9-10-4-7-5
##  5 1-7-8-3-6-4-9-5-2-10
##  6 1-4-8-3-5-2-7-10-9-6
##  7 1-9-7-8-6-5-10-4-3-2
##  8 1-7-4-10-3-6-2-9-5-8
##  9 1-3-10-4-6-8-5-7-9-2
## 10 1-9-6-10-5-4-3-8-2-7
```

```r
population$fits  <- sapply(population$chrom, fitness, cities)
population %<>% arrange(fits)

# Preserve elite individuals
nextInd.elite <- population %>% head(N_ELITE)
nextInd.elite %>% head %>% print()
```

```
## # A tibble: 2 x 2
##   chrom       fits
##   <list>     <dbl>
## 1 <dbl [10]>  5.14
## 2 <dbl [10]>  5.42
```

```r
# Generate children via: selection -> crossover -> mutation
#
pf <- population$fits
roulette.pie <-  (max(pf) -pf) / (max(pf) - min(pf))
roulette.pie <- roulette.pie / sum(roulette.pie)

roulette.pie %>% tibble()
```

```
## # A tibble: 10 x 1
##          .
##      <dbl>
##  1 0.221  
##  2 0.183  
##  3 0.158  
##  4 0.151  
##  5 0.145  
##  6 0.101  
##  7 0.0224 
##  8 0.0123 
##  9 0.00750
## 10 0
```

```r
pie(roulette.pie, clockwise = TRUE,
    main = "roulette based on fitness")
```

![](chap2.2_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

```r
roulette.thr <- cumsum(roulette.pie)
roulette.thr %>% tibble()
```

```
## # A tibble: 10 x 1
##        .
##    <dbl>
##  1 0.221
##  2 0.404
##  3 0.561
##  4 0.712
##  5 0.857
##  6 0.958
##  7 0.980
##  8 0.992
##  9 1.000
## 10 1.000
```

```r
which(roulette.thr > print(x <- runif(1)))[1]
```

```
## [1] 0.3776336
```

```
## [1] 2
```

```r
which(roulette.thr > 0.999)[1]
```

```
## [1] 9
```

```r
#

# selection -> crossover
chrom <- NULL
while(length(chrom) < NROW(population)- NROW(nextInd.elite)){
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
printChroms(chrom)
```

```
## # A tibble: 8 x 1
##   .                   
##   <chr>               
## 1 1-3-10-2-7-8-9-4-6-5
## 2 1-10-7-8-9-3-6-2-5-4
## 3 1-4-8-3-5-2-6-10-7-9
## 4 1-5-2-7-4-8-9-3-10-6
## 5 1-5-4-2-7-8-10-3-9-6
## 6 1-3-10-5-4-2-7-8-9-6
## 7 1-8-7-6-5-10-4-3-9-2
## 8 1-9-3-4-5-7-8-10-2-6
```

```r
# mutation
nextInd.gen <- tibble(
  chrom = lapply(chrom, mutation, mutate.prob = MUTATE_PROB))

all.equal(chrom, nextInd.gen$chrom)
```

```
## [1] "Component 1: Mean relative difference: 0.1176471"
## [2] "Component 6: Mean relative difference: 0.4615385"
```

```r
cbind(chrom = printChroms(chrom), nextInd = printChroms(nextInd.gen$chrom))
```

```
##                      .                    .
## 1 1-3-10-2-7-8-9-4-6-5 1-3-10-2-7-9-8-4-6-5
## 2 1-10-7-8-9-3-6-2-5-4 1-10-7-8-9-3-6-2-5-4
## 3 1-4-8-3-5-2-6-10-7-9 1-4-8-3-5-2-6-10-7-9
## 4 1-5-2-7-4-8-9-3-10-6 1-5-2-7-4-8-9-3-10-6
## 5 1-5-4-2-7-8-10-3-9-6 1-5-4-2-7-8-10-3-9-6
## 6 1-3-10-5-4-2-7-8-9-6 1-3-10-8-4-2-7-5-9-6
## 7 1-8-7-6-5-10-4-3-9-2 1-8-7-6-5-10-4-3-9-2
## 8 1-9-3-4-5-7-8-10-2-6 1-9-3-4-5-7-8-10-2-6
```

```r
population <- bind_rows(nextInd.elite, nextInd.gen)
population$fits  <- sapply(population$chrom, fitness, cities)
population %<>% arrange(fits)

population
```

```
## # A tibble: 10 x 2
##    chrom       fits
##    <list>     <dbl>
##  1 <dbl [10]>  5.14
##  2 <dbl [10]>  5.19
##  3 <dbl [10]>  5.32
##  4 <dbl [10]>  5.34
##  5 <dbl [10]>  5.42
##  6 <dbl [10]>  5.72
##  7 <dbl [10]>  5.82
##  8 <dbl [10]>  5.97
##  9 <dbl [10]>  6.12
## 10 <dbl [10]>  6.32
```

## functionise

### initialize 1st generation  


```r
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
pop <- initPopulation(20, cities)
printChroms(pop$chrom)
```

```
## # A tibble: 20 x 1
##    .                   
##    <chr>               
##  1 1-10-5-2-6-8-9-3-7-4
##  2 1-4-2-8-7-6-10-5-3-9
##  3 1-2-3-9-4-10-7-5-8-6
##  4 1-5-3-4-8-6-7-9-2-10
##  5 1-5-8-4-2-6-7-3-9-10
##  6 1-10-3-5-4-7-2-8-6-9
##  7 1-10-7-6-8-3-4-2-5-9
##  8 1-8-6-4-9-2-5-3-10-7
##  9 1-9-2-7-3-6-8-4-10-5
## 10 1-3-5-2-8-4-9-7-10-6
## 11 1-4-9-6-2-7-3-10-5-8
## 12 1-5-7-2-3-8-10-6-9-4
## 13 1-4-2-8-9-10-7-5-6-3
## 14 1-6-3-9-2-5-8-4-10-7
## 15 1-6-3-5-7-8-4-10-2-9
## 16 1-7-3-6-9-4-10-2-8-5
## 17 1-8-7-2-10-4-3-5-6-9
## 18 1-3-10-4-9-5-8-2-7-6
## 19 1-9-10-7-8-4-3-6-5-2
## 20 1-2-5-6-3-8-7-10-4-9
```

### alternate to next generation



```r
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
  population$fits  <- sapply(population$chrom, fitness, .cities)
  
  return(population %>% arrange(fits))
}
# example
(pop <- initPopulation(30, cities))
```

```
## # A tibble: 30 x 2
##    chrom       fits
##    <list>     <dbl>
##  1 <dbl [10]>  4.65
##  2 <dbl [10]>  4.77
##  3 <dbl [10]>  5.12
##  4 <dbl [10]>  5.16
##  5 <dbl [10]>  5.17
##  6 <dbl [10]>  5.17
##  7 <dbl [10]>  5.20
##  8 <dbl [10]>  5.38
##  9 <dbl [10]>  5.48
## 10 <dbl [10]>  5.50
## # ... with 20 more rows
```

```r
(pop.new <- alternate(pop, cities, elite.size = 3, mutate.prob = 0.8))
```

```
## # A tibble: 30 x 2
##    chrom       fits
##    <list>     <dbl>
##  1 <dbl [10]>  4.43
##  2 <dbl [10]>  4.54
##  3 <dbl [10]>  4.65
##  4 <dbl [10]>  4.72
##  5 <dbl [10]>  4.74
##  6 <dbl [10]>  4.76
##  7 <dbl [10]>  4.77
##  8 <dbl [10]>  5.11
##  9 <dbl [10]>  5.12
## 10 <dbl [10]>  5.24
## # ... with 20 more rows
```

### plot trip route


```r
plot.trip <- function(trip, .cities){
  stopifnot(!missing(trip),
            !missing(.cities), 
            NROW(.cities) == length(trip))
  travel <- data.frame(id = trip) %>% 
    left_join(.cities, by="id")
  travel <- rbind(travel, travel[1, ])
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
  invisible(list(trip = paste0(trip, collapse = "-"), 
                 travel = travel, total.trip = total.distance))
}
# example
sample.trip <- individual(cities)
fitness(sample.trip, cities)
```

```
## [1] 4.634169
```

```r
this.trip <- plot.trip(sample.trip, cities)
```

![](chap2.2_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

```r
print(this.trip)
```

```
## $trip
## [1] "1-4-3-7-9-2-5-6-8-10"
## 
## $travel
##    id      x      y
## 1   1 0.0000 0.0000
## 2   4 0.0697 0.0962
## 3   3 0.1156 0.7727
## 4   7 0.3398 0.5603
## 5   9 0.1657 0.9850
## 6   2 0.3977 0.2314
## 7   5 0.2436 0.4533
## 8   6 0.7916 0.0846
## 9   8 0.9714 0.0086
## 10 10 0.4587 0.3163
## 11  1 0.0000 0.0000
## 
## $total.trip
## [1] 4.634169
```


# example 1
## exec

```r
set.seed(7)
N_CITIES    = 10  # number of cities to travel
GEN_MAX     = 15  # number of generation
POP_SIZE    = 50  # population size
N_ELITE     = 3   # number of elite individual for next chromration
MUTATE_PROB = 0.1# mutation rate

factorial(N_CITIES)
```

```
## [1] 3628800
```


```r
start_time <- Sys.time()

(cities <- setCities(N_CITIES))
sample.trip <- individual(cities)
plot.trip(sample.trip, cities)
```

![](chap2.2_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

```r
generation <- list(NULL)

(pop <- initPopulation(POP_SIZE, cities))
for(i in 1:GEN_MAX){
  print(i)
  generation[[i]] <- pop
  pop <- alternate(pop,  cities, 
                   elite.size  = N_ELITE, 
                   mutate.prob = MUTATE_PROB)
}
Sys.time() - start_time
```


## eval


```r
top1 <- NULL
for(i in 1:length(generation)){
  this <- generation[[i]]
  top1 <- rbind(top1, tibble(gen=i, 
                             chrom = this$chrom[1], 
                             fits  = this$fits[1]))
}
top1
```

```
## # A tibble: 15 x 3
##      gen chrom       fits
##    <int> <list>     <dbl>
##  1     1 <dbl [10]>  4.71
##  2     2 <dbl [10]>  4.44
##  3     3 <dbl [10]>  4.17
##  4     4 <dbl [10]>  4.12
##  5     5 <dbl [10]>  3.72
##  6     6 <dbl [10]>  3.72
##  7     7 <dbl [10]>  3.72
##  8     8 <dbl [10]>  3.72
##  9     9 <dbl [10]>  3.72
## 10    10 <dbl [10]>  3.56
## 11    11 <dbl [10]>  3.56
## 12    12 <dbl [10]>  3.56
## 13    13 <dbl [10]>  3.56
## 14    14 <dbl [10]>  3.56
## 15    15 <dbl [10]>  3.56
```

## plot animation


```r
library(animation)
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
}, movie.name = "./output/stepGA_10cities-50steps.gif", 
interval = 0.7, ani.width=960, ani.height=480)
```

![example1](./output/stepGA_10cities-50steps.gif)


# example 2 (large ver.)
## exec

```r
set.seed(1)
N_CITIES    = 20  # number of cities to travel
GEN_MAX     = 350  # number of generation
POP_SIZE    = 150  # population size
N_ELITE     = 5   # number of elite individual for next chromration
MUTATE_PROB = 0.05 # mutation rate

factorial(N_CITIES)
```

```
## [1] 2.432902e+18
```


```r
start_time <- Sys.time()

(cities.L <- setCities(N_CITIES))

generation.L <- list(NULL)
(pop <- initPopulation(POP_SIZE, cities.L))
for(i in 1:GEN_MAX){
  print(i)
  generation.L[[i]] <- pop
  pop <- alternate(pop,  cities.L, 
                   elite.size  = N_ELITE, 
                   mutate.prob = MUTATE_PROB)
}
Sys.time() - start_time
```

## eval & plot animation


```r
top1.L <- NULL
for(i in 1:length(generation.L)){
  this <- generation.L[[i]]
  top1.L <- rbind(top1.L, tibble(gen=i, 
                             chrom = this$chrom[1], 
                             fits  = this$fits[1]))
}

pos <- which(!duplicated(top1.L$fits))
pos <- c(pos, rep(tail(pos,1), 5))
saveGIF({
  
  for(g in pos){
    par(mfrow = c(1,2))
    
    this.trip <- unlist(top1.L$chrom[g])
    plot.trip(this.trip, cities.L)
    
    plot(fits~gen, top1.L, type="b",
         main = sprintf("generation = %i (fits = %f)", g, top1.L$fits[g]))
    points(x=g, y=top1.L$fits[g], pch=16, col="red", cex=1.5)
    # top1;max(Y)
    
    par(mfrow = c(1,1))
  }
}, movie.name = "./output/stepGA_20cities.gif", 
interval = 0.8, ani.width=960, ani.height=480)
```

![example2](./output/stepGA_20cities.gif)



# example 3 (more large ver.)
## exec

```r
set.seed(1)
N_CITIES    = 30  # number of cities to travel
GEN_MAX     = 1000  # number of generation
POP_SIZE    = 150  # population size
N_ELITE     = 5   # number of elite individual for next chromration
MUTATE_PROB = 0.05 # mutation rate

factorial(N_CITIES)
```

```
## [1] 2.652529e+32
```


```r
start_time <- Sys.time()

(cities.LL <- setCities(N_CITIES))

generation.LL <- list(NULL)
(pop <- initPopulation(POP_SIZE, cities.LL))
for(i in 1:GEN_MAX){
  print(i)
  generation.LL[[i]] <- pop
  pop <- alternate(pop,  cities.LL, 
                   elite.size  = N_ELITE, 
                   mutate.prob = MUTATE_PROB)
}
Sys.time() - start_time

top1.LL <- NULL
for(i in 1:length(generation.LL)){
  this <- generation.LL[[i]]
  top1.LL <- rbind(top1.LL, tibble(gen=i, 
                                  chrom = this$chrom[1], 
                                  fits  = this$fits[1]))
}

this.trip <- unlist(top1.LL$chrom[g])
plot.trip(this.trip, cities.LL)
```

![](chap2.2_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

## eval & plot animation


```r
pos <- which(!duplicated(top1.LL$fits))
pos <- c(pos, rep(tail(pos,1), 5))
saveGIF({

  for(g in pos){
    par(mfrow = c(1,2))

    this.trip <- unlist(top1.LL$chrom[g])
    plot.trip(this.trip, cities.LL)

    plot(fits~gen, top1.LL, type="b",
         main = sprintf("generation = %i (fits = %f)", g, top1.LL$fits[g]))
    points(x=g, y=top1.LL$fits[g], pch=16, col="red", cex=1.5)
    # top1;max(Y)

    par(mfrow = c(1,1))
  }
}, movie.name = "./output/stepGA_30cities.gif",
interval = 0.8, ani.width=960, ani.height=480)
```

![example2](./output/stepGA_30cities.gif)

