# GA for Travelling salesman problem (TSP) 
# https://qiita.com/yama1223xxx/items/b609086de919e3af34d3
# http://www.obitko.com/tutorials/genetic-algorithms/japanese/index.php

set.seed(7)
rm(list = ls())

require(dplyr)
require(magrittr)
require(foreach)

# location of cities and travel routes --------------------------------------
mes <- seq(0,1,1e-4)
cities <- data.frame(id = 1:10, x =sample(mes, 10), y = sample(mes, 10))

travel <- data.frame(id = sample(NROW(cities))) %>% 
  left_join(cities, by="id") %>% 
  slice(c(1:NROW(.), 1))
travel

## calculate total trip. 
total.distance <- foreach(i=1:NROW(cities), .combine = sum) %do%({
  from <- travel[i, ]
  to <-  travel[i+1, ]
  
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

# example
x1 <- sample(NROW(cities))
x1


# order crossover
crossover <- function(p1, p2, show.pos = FALSE){
  len_chrom <- length(p1)
  stopifnot(len_chrom == length(p2))
  
  # position at crossover
  at <- sample(2:(length(p1)-1),2)
  at <- sort(at)
  
  child <- rep(NA, 10)
  child[at[1]:at[2]] <- p1[at[1]:at[2]]
  chid.inherit <- child
  
  p2.rot <- p2[c((at[2]+1):length(x2),1:(at[2]))]
  child.omit <- setdiff(p2, child)
  child[which(is.na(child))] <- child.omit
  
  if(show.pos){
    child <- list(at = at, p1 = p1, p2 = p2, 
                  p1.inherit = chid.inherit, p2.rotated = p2.rot,
                  p2.omitted = child.omit, new.chrom = child)
  }
  return(child)
}

chr1 <- sample(10)
chr2 <- sample(10)

crossover(chr1, chr2)
crossover(chr1, chr2, show.pos = TRUE)


# mutation (inversion)


x3 <- 1:10
at <- sample(2:(length(p1)-1),2)
at <- sort(at)
x3[pos] <- x3[rev(pos)]
x3
