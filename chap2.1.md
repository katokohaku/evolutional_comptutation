---
title: "Genetic Algorithm for getting max f(x) and argmax(x, f) with R"
author: "Satoshi Kato"
date: "2018/08/19"
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

## maximize the fitness function

![View](./output/stepGA_maximCurve.gif)


```r
fitnessCurve <- function(x){
  return(
    2 * sin(29 * x - 1) + 5 * sin(1 - 8 * x) + 3 * sin(1 - 70 * x)
  )
}
# example
X <- seq(0,1,0.001)
Y <- fitnessCurve(X)
plot(x=X, y=Y, type="l")
```

![](chap2.1_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

# individual




## encoding 

### genetic representation


```r
individual <- function(N){
  sample(c(0,1), N, replace=TRUE)
}
# example
(ind1 <- individual(20))
```

```
##  [1] 0 0 0 1 1 0 1 0 0 0 1 0 0 1 1 1 1 1 0 0
```

### decode to phenotype (-> fitness evaluation)

*phenotype* is input to fitness function.


```r
bin2dec <- function(chrom, norm = TRUE){
  stopifnot(all(chrom %in% c(0,1)) == TRUE)
  dec <- sum(2^(1:length(chrom) -1) * chrom)
  if(norm == TRUE){
    dec <- dec / sum(2^(1:length(chrom) -1))
  }
  
  return(dec)
}
# example
(ind1 <- individual(20))
```

```
##  [1] 0 0 1 0 1 0 0 0 0 0 1 0 0 0 0 1 0 0 0 1
```

```r
(pheno <- bin2dec(ind1))
```

```
## [1] 0.5322461
```

```r
fitnessCurve(pheno)
```

```
## [1] 5.467544
```

## genetic operator

### crossover 

In this case, 2-point-crossover is applied.


```r
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
```

```
##  [1] 1 1 1 1 1 1 1 1 1 1
```

```r
(chr2 <- rep(0, 10))
```

```
##  [1] 0 0 0 0 0 0 0 0 0 0
```

```r
crossover(chr1, chr2, show.pos = TRUE)
```

```
## $at
## [1] 7
## 
## $p1
##  [1] 1 1 1 1 1 1 1 1 1 1
## 
## $p2
##  [1] 0 0 0 0 0 0 0 0 0 0
## 
## $chrom
##  [1] 1 1 1 1 1 1 1 0 0 0
```

```r
crossover(chr1, chr2, show.pos = FALSE)
```

```
##  [1] 1 1 1 1 1 1 1 1 0 0
```

### mutation

In this case, single-point-mutation (bit flipping) is applied.


```r
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
```

```
##   [1] 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 0 1 1
##  [36] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1
##  [71] 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1
```



# population

## simple procedure

init 1st generation -> 2nd generation.


```r
GEN_MAX     = 15  # number of generation
POP_SIZE    = 10  # population size
CHROM_SIZE  = 64   # size of problem
N_ELITE     = 5   # number of elite individual for next chromration
N_TOURNAMEST= 5   # number of individual for tournament selection
MUTATE_PROB = 0.05 # mutation rate

population <- tibble(
chrom = foreach(i = 1:POP_SIZE) %do% individual(CHROM_SIZE))

population$pheno <- sapply(population$chrom, bin2dec)
population$fits  <- sapply(population$pheno, fitnessCurve)

population %<>% arrange(desc(fits))

# Preserve elite individuals
nextInd.elite <- population %>% head(N_ELITE)
nextInd.elite %>% head %>% print()
```

```
## # A tibble: 5 x 3
##   chrom       pheno  fits
##   <list>      <dbl> <dbl>
## 1 <dbl [64]> 0.0784 6.66 
## 2 <dbl [64]> 0.0643 4.91 
## 3 <dbl [64]> 0.640  2.71 
## 4 <dbl [64]> 0.907  1.55 
## 5 <dbl [64]> 0.656  0.589
```

```r
# Generate new children
# tournament selection -> crossover
nextChrom <- foreach(i = (N_ELITE+1):NROW(population)) %do% {
  parents <- population %>%
    sample_n(N_TOURNAMEST) %>%
    arrange(desc(fits)) %>%
    head(2)
  crossover(parents$chrom[1] %>% unlist,
            parents$chrom[2] %>% unlist)
}

# mutation
nextInd.gen <- tibble(
  chrom = lapply(nextChrom, mutation, mutate.prob = MUTATE_PROB))

nextInd.gen$pheno <- sapply(nextInd.gen$chrom, bin2dec)
nextInd.gen$fits  <- sapply(nextInd.gen$pheno, fitnessCurve)

nextGen.all <- rbind(nextInd.elite, nextInd.gen) %>%
  arrange(desc(fits))

nextGen.all
```

```
## # A tibble: 10 x 3
##    chrom       pheno   fits
##    <list>      <dbl>  <dbl>
##  1 <dbl [64]> 0.0784  6.66 
##  2 <dbl [64]> 0.0643  4.91 
##  3 <dbl [64]> 0.627   4.54 
##  4 <dbl [64]> 0.640   2.71 
##  5 <dbl [64]> 0.907   1.55 
##  6 <dbl [64]> 0.907   1.55 
##  7 <dbl [64]> 0.907   1.55 
##  8 <dbl [64]> 0.656   0.589
##  9 <dbl [64]> 0.656   0.588
## 10 <dbl [64]> 0.375  -5.85
```


## functionize

### initialize 1st generation  


```r
initPopulation <- function(pop.size, chrom.size){
  stopifnot(!missing(pop.size), !missing(chrom.size))
  
  population <- tibble(
    chrom = foreach(i = 1:pop.size) %do% individual(chrom.size))
  
  population$pheno <- sapply(population$chrom, bin2dec)
  population$fits  <- sapply(population$pheno, fitnessCurve)
  
  return(population %>% arrange(desc(fits)))
}
# example
(pop <- initPopulation(pop.size = 20, chrom.size = 32))
```

```
## # A tibble: 20 x 3
##    chrom       pheno   fits
##    <list>      <dbl>  <dbl>
##  1 <dbl [32]> 0.702   8.42 
##  2 <dbl [32]> 0.798   6.59 
##  3 <dbl [32]> 0.811   4.74 
##  4 <dbl [32]> 0.767   3.97 
##  5 <dbl [32]> 0.751   3.76 
##  6 <dbl [32]> 0.898   2.78 
##  7 <dbl [32]> 0.877   1.96 
##  8 <dbl [32]> 0.596   1.47 
##  9 <dbl [32]> 0.648   1.47 
## 10 <dbl [32]> 0.0258  0.904
## 11 <dbl [32]> 0.959   0.581
## 12 <dbl [32]> 0.174  -0.514
## 13 <dbl [32]> 0.345  -1.33 
## 14 <dbl [32]> 0.453  -1.34 
## 15 <dbl [32]> 0.352  -1.44 
## 16 <dbl [32]> 0.934  -2.31 
## 17 <dbl [32]> 0.487  -3.19 
## 18 <dbl [32]> 0.291  -4.53 
## 19 <dbl [32]> 0.293  -4.81 
## 20 <dbl [32]> 0.305  -5.94
```
### alternate to next generation


```r
alternate <- function(population, elite.size = 1, tournament.size = 4, mutate.prob = 0.01){
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
# example
(pop <- initPopulation(pop.size = 20, chrom.size = 32))
```

```
## # A tibble: 20 x 3
##    chrom         pheno   fits
##    <list>        <dbl>  <dbl>
##  1 <dbl [32]> 0.791     6.67 
##  2 <dbl [32]> 0.737     5.90 
##  3 <dbl [32]> 0.535     5.37 
##  4 <dbl [32]> 0.526     5.15 
##  5 <dbl [32]> 0.000124  5.04 
##  6 <dbl [32]> 0.540     5.01 
##  7 <dbl [32]> 0.748     4.09 
##  8 <dbl [32]> 0.681     3.19 
##  9 <dbl [32]> 0.644     2.02 
## 10 <dbl [32]> 0.959     0.540
## 11 <dbl [32]> 0.591     0.501
## 12 <dbl [32]> 0.174    -0.590
## 13 <dbl [32]> 0.441    -1.31 
## 14 <dbl [32]> 0.460    -2.00 
## 15 <dbl [32]> 0.465    -2.64 
## 16 <dbl [32]> 0.247    -2.70 
## 17 <dbl [32]> 0.332    -2.74 
## 18 <dbl [32]> 0.324    -4.15 
## 19 <dbl [32]> 0.198    -5.63 
## 20 <dbl [32]> 0.386    -8.04
```

```r
(alternate(pop))
```

```
## # A tibble: 20 x 3
##    chrom         pheno   fits
##    <list>        <dbl>  <dbl>
##  1 <dbl [32]> 0.729     7.59 
##  2 <dbl [32]> 0.732     6.92 
##  3 <dbl [32]> 0.791     6.67 
##  4 <dbl [32]> 0.737     5.90 
##  5 <dbl [32]> 0.535     5.37 
##  6 <dbl [32]> 0.526     5.15 
##  7 <dbl [32]> 0.000226  5.03 
##  8 <dbl [32]> 0.540     5.01 
##  9 <dbl [32]> 0.541     4.78 
## 10 <dbl [32]> 0.748     4.09 
## 11 <dbl [32]> 0.681     3.22 
## 12 <dbl [32]> 0.681     3.19 
## 13 <dbl [32]> 0.681     3.19 
## 14 <dbl [32]> 0.644     2.17 
## 15 <dbl [32]> 0.555     2.10 
## 16 <dbl [32]> 0.650     1.15 
## 17 <dbl [32]> 0.665     0.545
## 18 <dbl [32]> 0.441    -1.31 
## 19 <dbl [32]> 0.458    -1.79 
## 20 <dbl [32]> 0.930    -2.18
```

# sample

getting max f(x) and argmax(x, f)

## exec


```r
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
Sys.time() - start_time
```

```
## Time difference of 2.366189 secs
```


## eval

```r
top1 <- NULL
for(i in 1:length(generation)){
  this <- generation[[i]]
  top1 <- rbind(top1, data.frame(gen=i, fits=this$fits[1]))
}

top1$diff <- max(Y) - top1$fits
print(top1)
```

```
##    gen     fits          diff
## 1    1 9.402023  0.0009377269
## 2    2 9.402023  0.0009377269
## 3    3 9.402023  0.0009377269
## 4    4 9.402023  0.0009377269
## 5    5 9.402023  0.0009377269
## 6    6 9.402399  0.0005608239
## 7    7 9.402783  0.0001767947
## 8    8 9.402783  0.0001767947
## 9    9 9.402783  0.0001767947
## 10  10 9.402783  0.0001767947
## 11  11 9.403708 -0.0007482118
## 12  12 9.403708 -0.0007482118
## 13  13 9.403708 -0.0007482118
## 14  14 9.403835 -0.0008743498
## 15  15 9.403835 -0.0008743498
## 16  16 9.403835 -0.0008743498
## 17  17 9.403835 -0.0008743498
## 18  18 9.403835 -0.0008743498
## 19  19 9.403835 -0.0008743498
## 20  20 9.403835 -0.0008743498
```


## plot animation


```r
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
}, interval = 1.0, movie.name = "./output/stepGA_maximCurve.gif", ani.width=960, ani.height=480)
```

![View](./output/stepGA_maximCurve.gif)










