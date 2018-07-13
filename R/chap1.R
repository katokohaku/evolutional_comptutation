require(dplyr)

GEN_MAX     = 1000 # number of generation
POP_SIZE    = 1000 # population size
ELITE       = 1    # number of elite individual for next generation
MUTATE_PROB = 0.01 # mutation rate
N           = 64   # size of problem (= number of element)


fitness <- function(gene){
  stopifnot(all(gene %in% c(-1,1)) == TRUE)
  return(
    (gene * seq_along(gene)^0.5) %>% sum %>% abs
  )
}

individual_init <- function(N){
  stopifnot(N > 1)

  gene <- sample(c(-1,1), N, replace=TRUE)
  return(gene)  
}

g1 <- rep(1,10)
g2 <- rep(0,10)

crossover <- function(p1, p2){
  len_gene <- length(p1)
  stopifnot(len_gene == length(p2))
  
  pos <- sample(len_gene, 1)
  child <- list(pos = pos,
                gene1 = c(p1[1:pos], p2[(pos+1):len_gene]),
                gene2 = c(p2[1:pos], p1[(pos+1):len_gene]))
  return(child)
}
g1;g2
crossover(g1,g2, show.pos = TRUE)


mutation <- function(gene, mutate.prob){
  stopifnot(all(gene %in% c(-1,1)) == TRUE)
  
  
}