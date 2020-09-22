#' export from DDD package, returns newick instead of tree
#' @param L ltable
#' @param dropextinct drop extinct
#' @return newick
#' @rawNamespace useDynLib(enviDiv)
#' @rawNamespace importFrom(RcppParallel, RcppParallelLibs)
#' @export
ltable_to_phylo <- function(L, dropextinct = T) {
  L <- L[order(abs(L[, 3])), 1:4]
  age <- L[1, 1]
  L[, 1] <- age - L[, 1]
  L[1, 1] <- -1
  notmin1 <- which(L[, 4] != -1)
  L[notmin1, 4] <- age - L[notmin1, 4]
  if (dropextinct == T) {
    sall <- which(L[, 4] == -1)
    tend <- age
  } else {
    sall <- which(L[, 4] >= -1)
    tend <- (L[, 4] == -1) * age + (L[, 4] > -1) * L[, 4]
  }
  L <- L[, -4]
  linlist <- cbind(data.frame(L[sall, ]), paste("t", abs(L[sall,
                                                           3]), sep = ""), tend)
  linlist[, 4] <- as.character(linlist[, 4])
  names(linlist) <- 1:5
  done <- 0
  while (done == 0) {
    j <- which.max(linlist[, 1])
    parent <- linlist[j, 2]
    parentj <- which(parent == linlist[, 3])
    parentinlist <- length(parentj)

    if (parentinlist == 1) {
      spec1 <- paste(linlist[parentj, 4], ":",
                     linlist[parentj, 5] - linlist[j, 1], sep = "")
      spec2 <- paste(linlist[j, 4], ":",
                     linlist[j, 5] - linlist[j, 1], sep = "")

      cat("spec1: ", spec1, "\n")
      cat("spec2: ", spec2, "\n")
      linlist[parentj, 4] <- paste("(", spec1, ",", spec2,
                                   ")", sep = "")
      linlist[parentj, 5] <- linlist[j, 1]
      linlist <- linlist[-j, ]
    }
    else {
      linlist[j, 1:3] <- L[which(L[, 3] == parent), 1:3]
    }
    if (nrow(linlist) == 1) {
      done <- 1
    }
  }
  linlist[4] <- paste(linlist[4], ":", linlist[5], ";", sep = "")
  return(linlist[1, 4])
}


#' internal ltab to ltable to phylo
#' @param ltab ltab
#' @param crown_age crown age
#' @param dropextinct drop extinct
#' @return newick
#' @export
ltab_to_ltable_to_phylo <- function(ltab,
                                    crown_age,
                                    dropextinct = T) {
  local_l_table <- ltab
  local_l_table[, 1] <- crown_age - local_l_table[, 1]
  local_l_table <-  local_l_table[order(abs(local_l_table[, 3])), 1:4]

  local_l_table[1, 2] <- 0
  local_l_table[which(local_l_table[, 1] < 0), 1] <- 0

  a <- subset(local_l_table, local_l_table[, 1] == crown_age)
  connected <- FALSE
  if (a[2, 3] == a[1, 2]) connected <- TRUE
  if (a[1, 3] == a[2, 2]) connected <- TRUE

  if (connected == FALSE) {
    parent_id <- local_l_table[1, 3]
    local_l_table[which(local_l_table[, 2] == -1), 2] <- parent_id
  }

  newick_string <- ltable_to_phylo(local_l_table, dropextinct)
  return(newick_string)
}


#' @keywords internal
split_into_blocks <- function(m,
                              block_size,
                              nb = ceiling(m / block_size))  {
       if (nb > m)
           nb <- m
       int <- m/nb
       upper <- round(1:nb * int)
       lower <- c(1, upper[-nb] + 1)
       size <- c(upper[1], diff(upper))
       cbind(lower, upper, size)
}
