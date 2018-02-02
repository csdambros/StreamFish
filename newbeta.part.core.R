##### Functional diversity of stream fishes in Brazil - additional functions
##### Created by CSDambros in Feb 1st 2018

#' Recreate the functional.betapart.core function from the betapart package
#' The function implements the prefix argument, necessary to run in parallel
#' When the function runs, it generates an external file. The implementation allows using different names for these files for each run when running in parallel

#' @param prefix a name to be attached at the beginning of the vertex files

functional.betapart.core <-
  function (x,
            traits,
            multi = TRUE,
            warning.time = TRUE,
            return.details = FALSE,
            prefix = "")
  {
    require(betapart)
    
    if (!is.matrix(x)) {
      x <- as.matrix(x)
    }
    if (!is.numeric(x))
      stop("The data in 'x' is not numeric.", call. = TRUE)
    xvals <- unique(as.vector(x))
    if (any(!is.element(xvals, c(0, 1))))
      stop(
        "The 'x' table contains values other than 0 and 1: data should be presence/absence.",
        call. = TRUE
      )
    if (!is.numeric(traits))
      stop("The data in 'traits' is not numeric.", call. = TRUE)
    if (any(is.na(traits)))
      stop("NA are not allowed in 'traits'", call. = TRUE)
    if (ncol(x) != nrow(traits))
      stop("Number of species in 'x' and 'traits' must be identical",
           call. = TRUE)
    D <- ncol(traits)
    Si <- apply(x, 1, sum)
    if (any(Si <= D))
      stop(paste(
        "'community ",
        row.names(x)[which(Si <= D)],
        " must contain at least ",
        D + 1,
        " species",
        sep = ""
      ))
    N <- nrow(x)
    if (N < 2)
      stop("Computing dissimilairty requires at least 2 communities",
           call. = TRUE)
    nb.step <- 2
    if (multi == T)
      nb.step <- N
    step.fbc <-
      as.data.frame(matrix("", nb.step, 1, dimnames = list(c(
        "           FRi",
        paste("intersection", 2:nb.step, sep = "_")
      ), c("iteration"))))
    step.fbc[, 1] <- as.character(step.fbc[, 1])
    step.fbc[1, 1] <- paste("0/", N, sep = "")
    for (k in 2:nb.step)
      step.fbc[k, 1] <- paste("0/", choose(N,
                                           k), sep = "")
    FRi <- rep(NA, N)
    names(FRi) <- row.names(x)
    coord_vert_i <- list()
    for (i in 1:N) {
      tr_i <- traits[which(x[i,] == 1),]
      vert0 <-
        convhulln(tr_i, paste0("Fx TO", " '", paste0(prefix, "vert.txt"), "'"))
      vert1 <- scan(paste0(prefix, "vert.txt"), quiet = T)
      verti <- (vert1 + 1)[-1]
      coord_vert_i[[i]] <- tr_i[verti,]
      FRi[i] <- convhulln(tr_i[verti,], "FA")$vol
      step.fbc["           FRi", 1] <- paste(i, "/", N, sep = "")
      step.fbc[, 1] <- as.character(step.fbc[, 1])
      write.table(
        step.fbc,
        file = paste0(prefix, "step.fbc.txt"),
        row.names = T,
        col.names = F,
        sep = "\t"
      )
    }
    sumFRi <- sum(FRi)
    intersect <- function(set1, set2) {
      set1rep <- d2q(cbind(0, cbind(1, set1)))
      set2rep <- d2q(cbind(0, cbind(1, set2)))
      polytope1 <- redundant(set1rep, representation = "V")$output
      polytope2 <- redundant(set2rep, representation = "V")$output
      H_chset1 <- scdd(polytope1, representation = "V")$output
      H_chset2 <- scdd(polytope2, representation = "V")$output
      H_inter <- rbind(H_chset1, H_chset2)
      V_inter <- scdd(H_inter, representation = "H")$output
      vert_1n2 <- q2d(V_inter[,-c(1, 2)])
      coord_vert_inter <- rep(NA, ncol(set1))
      vol_inter <- 0
      if (is.matrix(vert_1n2) == T)
        if (nrow(vert_1n2) > ncol(vert_1n2)) {
          coord_vert_inter <- vert_1n2
          vol_inter <- convhulln(vert_1n2, "FA")$vol
        }
      res <-
        list(coord_vert_inter = coord_vert_inter, vol_inter = vol_inter)
      return(res)
    }
    comb2 <- combn(1:N, 2, simplify = T)
    vol_inter2_mat <- matrix(0, N, N, dimnames = list(row.names(x),
                                                      row.names(x)))
    vol_inter2 <- rep(0, ncol(comb2))
    coord_vert_inter2 <- list()
    for (k in 1:ncol(comb2)) {
      i <- comb2[1, k]
      j <- comb2[2, k]
      seti <- traits[which(x[i,] == 1),]
      setj <- traits[which(x[j,] == 1),]
      interij <- intersect(seti, setj)
      vol_inter2_mat[j, i] <- interij$vol_inter
      vol_inter2[k] <- interij$vol_inter
      coord_vert_inter2[[k]] <- interij$coord_vert_inter
      step.fbc["intersection_2", 1] <- paste(k, "/", ncol(comb2),
                                             sep = "")
      write.table(
        step.fbc,
        file = paste0(prefix, "step.fbc.txt"),
        row.names = T,
        col.names = F,
        sep = "\t"
      )
    }
    matNN <-
      matrix(0, N, N, dimnames = list(row.names(x), row.names(x)))
    shared <- matNN
    not.shared <- matNN
    for (i in 1:(N - 1))
      for (j in (i + 1):N) {
        shared[j, i] <- vol_inter2_mat[j, i]
        not.shared[i, j] <- FRi[i] - vol_inter2_mat[j, i]
        not.shared[j, i] <- FRi[j] - vol_inter2_mat[j, i]
      }
    sum.not.shared <- not.shared + t(not.shared)
    max.not.shared <- pmax(not.shared, t(not.shared))
    min.not.shared <- pmin(not.shared, t(not.shared))
    comb_inter <- list()
    comb_inter[[1]] <- comb2
    coord_vert_inter <- list()
    coord_vert_inter[[1]] <- coord_vert_inter2
    vol_inter <- list()
    vol_inter[[1]] <- vol_inter2
    FRt <- NA
    a <- NA
    if (N > 2 & multi == T) {
      if (warning.time == T & N > 10)
        stop(
          paste(
            "Computing mulitple functional dissimilarity on more than 10 communities may take a long time. \n    \t\t\t\t\t\t\t\t\tSet 'multi' or 'warning.time' to FALSE"
          )
        )
      if (warning.time == T & D > 4)
        stop(
          paste(
            "Computing mulitple functional dissimilarity in a",
            D,
            "-dimensions functional space may take a long time. \n    \t\t\t\t\t\t\t\t\tSet 'multi' or 'warning.time' to FALSE"
          )
        )
      for (z in 3:N) {
        comb_z <- combn(1:N, z, simplify = T)
        vol_inter_z <- rep(0, ncol(comb_z))
        coord_vert_inter_z <- list()
        for (k in 1:ncol(comb_z)) {
          seti <- coord_vert_inter[[z - 2]][[which(apply(comb_inter[[z -
                                                                       2]], 2, identical, comb_z[1:(z - 1), k]) ==
                                                     T)]]
          setj <-
            coord_vert_inter[[z - 2]][[which(apply(comb_inter[[z -
                                                                 2]], 2, identical, comb_z[2:z, k]) == T)]]
          coord_vert_inter_z[[k]] <- rep(NA, D)
          if (is.na(sum(seti) + sum(setj)) == F) {
            interij <- intersect(seti, setj)
            vol_inter_z[k] <- interij$vol_inter
            coord_vert_inter_z[[k]] <- interij$coord_vert_inter
          }
          step.fbc[paste("intersection", z, sep = "_"),
                   1] <- paste(k, "/", ncol(comb_z), sep = "")
          write.table(
            step.fbc,
            file = paste0(prefix, "step.fbc.txt"),
            row.names = T,
            col.names = F,
            sep = "\t"
          )
        }
        comb_inter[[z - 1]] <- comb_z
        coord_vert_inter[[z - 1]] <- coord_vert_inter_z
        vol_inter[[z - 1]] <- vol_inter_z
      }
      sumvol_sign <- rep(NA, N - 1)
      for (k in 2:N) {
        sumvol_sign[k - 1] <- (-1) ^ (k - 1) * sum(vol_inter[[k -
                                                                1]])
      }
      FRt <- sumFRi + sum(sumvol_sign)
      a <- sumFRi - FRt
    }
    details <- NA
    if (return.details == T) {
      CH <- list(FRi = FRi, coord_vertices = coord_vert_i)
      intersections <-
        list(
          combinations = comb_inter,
          volumes = vol_inter,
          coord_vertices = coord_vert_inter
        )
      details <- list(CH = CH, intersections = intersections)
    }
    functional.computations <- list(
      sumFRi = sumFRi,
      FRt = FRt,
      a = a,
      shared = shared,
      not.shared = not.shared,
      sum.not.shared = sum.not.shared,
      max.not.shared = max.not.shared,
      min.not.shared = min.not.shared,
      details = details
    )
    class(functional.computations) <- "functional.betapart"
    return(functional.computations)
  }


#' Recreate the functional.beta.pair function from the betapart package
#' The function implements the prefix argument, necessary to run in parallel
#' When the function runs, it generates an external file. The implementation allows using different names for these files for each run when running in parallel
#' The parameter is passed to the function functional.betapart.core

#' @param prefix a name to be attached at the beginning of the vertex files

functional.beta.pair <-
  function (x,
            traits,
            index.family = "sorensen",
            prefix = "")
  {
    require(betapart)
    index.family <- match.arg(index.family, c("jaccard", "sorensen"))
    fbc <- x
    if (!inherits(x, "functional.betapart")) {
      fbc <- functional.betapart.core(
        x,
        traits,
        multi = FALSE,
        warning.time = FALSE,
        return.details = FALSE,
        prefix = prefix
      )
    }
    switch(index.family, sorensen = {
      funct.beta.sim <- fbc$min.not.shared / (fbc$min.not.shared +
                                                fbc$shared)
      funct.beta.sne <-
        ((fbc$max.not.shared - fbc$min.not.shared) / ((2 *
                                                         fbc$shared) + fbc$sum.not.shared)) * (fbc$shared /
                                                                                                 (fbc$min.not.shared +
                                                                                                    fbc$shared))
      funct.beta.sor <- fbc$sum.not.shared / (2 * fbc$shared +
                                                fbc$sum.not.shared)
      functional.pairwise <-
        list(
          funct.beta.sim = as.dist(funct.beta.sim),
          funct.beta.sne = as.dist(funct.beta.sne),
          funct.beta.sor = as.dist(funct.beta.sor)
        )
    }, jaccard = {
      funct.beta.jtu <-
        (2 * fbc$min.not.shared) / ((2 * fbc$min.not.shared) +
                                      fbc$shared)
      funct.beta.jne <-
        ((fbc$max.not.shared - fbc$min.not.shared) / (fbc$shared +
                                                        fbc$sum.not.shared)) * (fbc$shared /
                                                                                  ((2 * fbc$min.not.shared) +
                                                                                     fbc$shared))
      funct.beta.jac <-
        fbc$sum.not.shared / (fbc$shared + fbc$sum.not.shared)
      functional.pairwise <-
        list(
          funct.beta.jtu = as.dist(funct.beta.jtu),
          funct.beta.jne = as.dist(funct.beta.jne),
          funct.beta.jac = as.dist(funct.beta.jac)
        )
    })
    return(functional.pairwise)
  }

# Include the functional.betapart.core into the betapart environment (required)
library(betapart)
environment(functional.betapart.core) <- environment(betapart.core)

#' Null model
#' The function randomizes species traits while maintaining the changes in species taxonomic composition and linkage between traits intact (fixed algorithm)

#' @param i a single number to be used in a for or lapply loop


functional.beta.pair.null <-
  function(i,
           x,
           traits,
           index.family = "jaccard",
           prefix = "",gower=TRUE,nreps=99) {

    tryCatch({

      
      
      return(func.pair.random)
    }, error=function(e){})
  }

###

functional.beta.pair.random <-
  function(x,
           traits,
           index.family = "jaccard",
           prefix = "",gower=TRUE) {
    require(FD)
    require(betapart)
    
      tryCatch({
    traits.random <- traits[sample(1:nrow(traits)), ]
    rownames(traits.random) <- rownames(traits)
    
    if(gower==TRUE){
    gower.dist <- gowdis(traits.random)
    traits.random <- cmdscale(gower.dist, 2)
    }
    
    func.pair.random <-
      functional.beta.pair(
        x = x,
        traits = traits.random,
        index.family = index.family,
        prefix = paste0(prefix, i)
      )
    
    return(func.pair.random)
      }, error=function(e){})
  }



######### Original (working)


functional.beta.pair.random <-
  function(i,
           x,
           traits,
           index.family = "jaccard",
           prefix = "",gower=TRUE) {
    require(FD)
    require(betapart)
    
    tryCatch({
      traits.random <- traits[sample(1:nrow(traits)), ]
      rownames(traits.random) <- rownames(traits)
      
      if(gower==TRUE){
        gower.dist <- gowdis(traits.random)
        traits.random <- cmdscale(gower.dist, 2)
      }
      
      func.pair.random <-
        functional.beta.pair(
          x = x,
          traits = traits.random,
          index.family = index.family,
          prefix = paste0(prefix, i)
        )
      
      return(func.pair.random)
    }, error=function(e){})
  }


