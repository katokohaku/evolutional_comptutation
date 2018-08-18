---
title: "進化計算アルゴリズム入門(R版)"
author: "Satoshi Kato"
date: "2018/08/15"
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

## 目的

GA for getting argmax(f(x)) 
[http://www.obitko.com/tutorials/genetic-algorithms/japanese/index.php]

# individual






```r
individual <- function(N){
  sample(c(0,1), N, replace=TRUE)
}
# example
(ind1 <- individual(64))
```

```
##  [1] 0 1 0 0 0 0 0 1 0 1 1 1 1 0 1 0 0 0 1 1 0 1 1 1 1 0 1 1 1 1 1 1 1 1 0
## [36] 1 0 0 1 0 1 1 1 0 0 0 0 1 0 0 1 1 1 0 0 1 0 0 0 0 1 0 1 0
```



```r
fitnessCurve <- function(x){
  return(
    2*sin(29*x-1) + 5*sin(1- 8*x) + 3*sin(1- 70*x)
  )
}
# example
X <- seq(0,1,0.001)
Y <- fitnessCurve(X)
plot(x=X, y=Y, type="l")
```

![](chap2.1_files/figure-html/unnamed-chunk-2-1.png)<!-- -->



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
(ind1 <- individual(64))
```

```
##  [1] 0 0 1 0 1 0 0 0 1 0 1 0 1 1 0 0 1 1 1 1 0 0 1 1 1 0 1 0 0 1 0 0 1 0 0
## [36] 0 0 0 0 1 0 1 0 0 0 1 1 1 0 0 1 1 1 1 0 1 0 1 0 0 1 0 1 1
```

```r
fitnessCurve( bin2dec(ind1))
```

```
## [1] 1.443962
```



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
## [1] 5
## 
## $p1
##  [1] 1 1 1 1 1 1 1 1 1 1
## 
## $p2
##  [1] 0 0 0 0 0 0 0 0 0 0
## 
## $chrom
##  [1] 1 1 1 1 1 0 0 0 0 0
```

```r
crossover(chr1, chr2, show.pos = FALSE)
```

```
##  [1] 1 1 1 1 1 1 1 1 0 0
```



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
##   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
##  [36] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1
##  [71] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
```



# population


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
##   chrom      pheno  fits
##   <list>     <dbl> <dbl>
## 1 <dbl [64]> 0.780 5.56 
## 2 <dbl [64]> 0.775 4.96 
## 3 <dbl [64]> 0.813 4.29 
## 4 <dbl [64]> 0.979 1.96 
## 5 <dbl [64]> 0.653 0.827
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
##    chrom      pheno  fits
##    <list>     <dbl> <dbl>
##  1 <dbl [64]> 0.780 5.56 
##  2 <dbl [64]> 0.529 5.40 
##  3 <dbl [64]> 0.775 4.96 
##  4 <dbl [64]> 0.775 4.96 
##  5 <dbl [64]> 0.813 4.29 
##  6 <dbl [64]> 0.979 1.96 
##  7 <dbl [64]> 0.980 1.90 
##  8 <dbl [64]> 0.653 0.827
##  9 <dbl [64]> 0.655 0.691
## 10 <dbl [64]> 0.563 0.588
```


# popuration function


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
##    chrom      pheno    fits
##    <list>     <dbl>   <dbl>
##  1 <dbl [32]> 0.535  5.41  
##  2 <dbl [32]> 0.620  4.74  
##  3 <dbl [32]> 0.611  4.13  
##  4 <dbl [32]> 0.680  3.04  
##  5 <dbl [32]> 0.553  2.70  
##  6 <dbl [32]> 0.881  2.52  
##  7 <dbl [32]> 0.822  1.65  
##  8 <dbl [32]> 0.655  0.669 
##  9 <dbl [32]> 0.913  0.413 
## 10 <dbl [32]> 0.168 -0.0906
## 11 <dbl [32]> 0.169 -0.135 
## 12 <dbl [32]> 0.499 -0.831 
## 13 <dbl [32]> 0.119 -1.07  
## 14 <dbl [32]> 0.445 -1.12  
## 15 <dbl [32]> 0.345 -1.34  
## 16 <dbl [32]> 0.127 -2.19  
## 17 <dbl [32]> 0.463 -2.42  
## 18 <dbl [32]> 0.296 -5.32  
## 19 <dbl [32]> 0.235 -5.45  
## 20 <dbl [32]> 0.226 -7.26
```



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
##    chrom       pheno   fits
##    <list>      <dbl>  <dbl>
##  1 <dbl [32]> 0.798   6.60 
##  2 <dbl [32]> 0.0910  5.71 
##  3 <dbl [32]> 0.0602  4.05 
##  4 <dbl [32]> 0.758   3.46 
##  5 <dbl [32]> 0.512   2.44 
##  6 <dbl [32]> 0.959   0.598
##  7 <dbl [32]> 0.502  -0.220
##  8 <dbl [32]> 0.116  -0.371
##  9 <dbl [32]> 0.573  -0.646
## 10 <dbl [32]> 0.120  -1.20 
## 11 <dbl [32]> 0.345  -1.34 
## 12 <dbl [32]> 0.353  -1.49 
## 13 <dbl [32]> 0.182  -1.72 
## 14 <dbl [32]> 0.249  -2.29 
## 15 <dbl [32]> 0.332  -2.85 
## 16 <dbl [32]> 0.370  -4.50 
## 17 <dbl [32]> 0.419  -5.37 
## 18 <dbl [32]> 0.300  -5.65 
## 19 <dbl [32]> 0.414  -6.52 
## 20 <dbl [32]> 0.222  -7.72
```

```r
(alternate(pop))
```

```
## # A tibble: 20 x 3
##    chrom       pheno   fits
##    <list>      <dbl>  <dbl>
##  1 <dbl [32]> 0.798   6.60 
##  2 <dbl [32]> 0.0597  3.93 
##  3 <dbl [32]> 0.758   3.46 
##  4 <dbl [32]> 0.758   3.46 
##  5 <dbl [32]> 0.0531  2.51 
##  6 <dbl [32]> 0.512   2.44 
##  7 <dbl [32]> 0.512   2.44 
##  8 <dbl [32]> 0.967   1.67 
##  9 <dbl [32]> 0.0479  1.53 
## 10 <dbl [32]> 0.959   0.598
## 11 <dbl [32]> 0.959   0.596
## 12 <dbl [32]> 0.502  -0.188
## 13 <dbl [32]> 0.115  -0.212
## 14 <dbl [32]> 0.501  -0.382
## 15 <dbl [32]> 0.573  -0.646
## 16 <dbl [32]> 0.573  -0.646
## 17 <dbl [32]> 0.574  -0.691
## 18 <dbl [32]> 0.120  -1.19 
## 19 <dbl [32]> 0.120  -1.20 
## 20 <dbl [32]> 0.353  -1.49
```

# exec


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
## Time difference of 2.031412 secs
```


# eval

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


# animation


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

```
## Warning in system(cmd, intern = intern, wait = wait | intern,
## show.output.on.console = wait, : 命令 'C:\WINDOWS\system32\cmd.exe /c
## convert --version' の実行は状態 4 を持ちました
```

```
## I cannot find ImageMagick with convert = "convert"
```

```
## but I can find it from the Registry Hive: C:\ImageMagick-7.0.8-Q16
```

```
## Executing: 
## "C:\ImageMagick-7.0.8-Q16\convert.exe -loop 0 -delay 100
##     Rplot1.png Rplot2.png Rplot3.png Rplot4.png Rplot5.png
##     Rplot6.png Rplot7.png Rplot8.png Rplot9.png Rplot10.png
##     Rplot11.png Rplot12.png Rplot13.png Rplot14.png Rplot15.png
##     Rplot16.png Rplot17.png Rplot18.png Rplot19.png Rplot20.png
##     "stepGA_maximCurve.gif""
```

```
## Output at: ./output/stepGA_maximCurve.gif
```

```
## [1] TRUE
```

![様子](./output/stepGA_maximCurve.gif)










