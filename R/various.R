
#' Generating all possible haplotypes for a set of sample and their corresponding AS
#'
#' This function allows you to compute all possible haplotypes for an individual for a set of fSNPs and each rSNP
#' @param x data table with GT and AS for fSNPs and rSNPs
#' @param y name of samples to make haps
#' @export
#' @return list
#' haps()

haps <- function(x, y) {
  haps.samples <- list()
  for (i in seq_along(y)) {
    cols <- grep(y[i], names(x), value = T)
    # get # het rows in fSNP
    # split GT and AS columns into haplotype a and b
    z <- lapply(cols, function(i) {
      tmp <- x[SNP == "fSNP", i, with = F]
      tmp[, (i) := sub("\\|", ",", get(i))]
      hap.a <- as.numeric(unlist(strsplit(unlist(tmp), ","))[c(TRUE, FALSE)])
      hap.b <- as.numeric(unlist(strsplit(unlist(tmp), ","))[c(FALSE, TRUE)])
      DT <- data.table::data.table(hap.a, hap.b)
      names(DT) <- paste0(i, c(".hap.a", ".hap.b"))
      return(DT)
    })
    z <- data.table::data.table(do.call("cbind", z))
    het.rows <- which(x[SNP == "fSNP", get(cols[1])] == "0|1" | x[SNP == "fSNP", get(cols[1])] == "1|0")
    if (length(het.rows) <= 1) {
      haps.samples[[i]] <- z
    } else {
      haps.samples[[i]] <- haps_sub(z, het.rows)
    }
  }
  names(haps.samples) <- y
  return(haps.samples)
}


#' Generating all possible haplotypes for a sample: sub function for haps: based on Chris F code and on binomial theorem: (a+b)^n= Sumk=0,n nchoosek a^k.b^(n-k)
#'
#' This function allows you to compute all possible haplotypes for an individual for a set of SNPs
#' @param z data table with GT for 1 individual
#' @param het.rows positions in z that correspond to het SNPs
#' @keywords all haplotype
#' @return list of lists. Inner list haplotypes per individual, outer list each individual
#' haps_sub()

haps_sub <- function(z, het.rows) {
  l <- length(het.rows)
  haps <- list()
  haps[[1]] <- z
  index <- 2
  for (k in 1:l) {
    combs <- choose(l, k)
    tmp2 <- iterpc(l, k, c(1:l))
    for (j in 1:combs) {
      tmp <- copy(z)
      rws <- t(getnext(tmp2))
      tmp[het.rows[rws], c(1, 3) := z[het.rows[rws], c(2, 4), with = F]]
      tmp[het.rows[rws], c(2, 4) := z[het.rows[rws], c(1, 3), with = F]]
      haps[[index]] <- tmp
      index <- index + 1
    }
  }
  return(haps)
}


#' Sum haplotypes to get genotypes
#'
#' This function allows you to sum a pair of haplotypes (strings) and return a string with the sum
#' @param a vector with haplotypes 1
#' @param b vector with haplotypes
#' @keywords sum haplotype pairs
#' @export
#' @return vector  with the sum of the haplotype pair (genotype)
#' add.geno()

add.geno <- function(a, b) {
  an <- lapply(strsplit(a, ""), as.numeric)
  bn <- lapply(strsplit(b, ""), as.numeric)
  ## lapply(mapply("+",an,bn,SIMPLIFY=FALSE),paste,collapse="")
  u <- unlist(lapply(mapply("+", an, bn, SIMPLIFY = FALSE), paste, collapse = ""))
  return(u)
}


#' Calculate probability of population haplotype pairs and return sparse matrix
#'
#' This function allows you to  calculate frequency of haplotype pairs in population, returns a sparse matrix
#' @param h matrix with population haplotypes
#' @keywords simulation probability haplotype pairs
#' @export
#' @return sparse matrix with probabilities of haplotype pairs and genotype
#' p.hap.pair.s()

p.hap.pair.s <- function(h) {
  DT <- addstrings::p.hap.pair.dt(h)
  ## create matrix of frequency of haplotype by genotype
  M2 <- Matrix::sparseMatrix(i = 1:nrow(DT), j = match(DT$geno, unique(DT$geno)), x = DT$freq, dims = c(nrow(DT), length(unique(DT$geno))), dimnames = list(DT[, haps], unique(DT$geno)))
  return(M2)
}


#' Format a named list of stan output for baseqtl
#'
#' This function allows you to extract the parameter information from a list of stan summaries (output from stan.many.sim)
#' @param x list of summaries
#' @param y parameter to extract from summary, NULL if already extracted
#' @param rtag optional argument,whether snps were grouped using tag function
#' @param model, character vector indicating which model was run: full or neg.only
#' @param nhets, vector with the number of hets  for each rsnp
#' @param ASE.hets, vector with the number of hets with sufficient ASE counts for each rsnp
#' @param gene gene id of gene under study
#' @param EAF, data table with snp and eaf for tag snps, output from snp.eaf (real.data.R), defaults to NULL
#' @param info, named vector with info type score for tag snp, names snp_id, defaults to NULL
#' @param nfsnps, number of fsnps in analysis, defaults to NULL
#' @param min.pval data table with cols fsnp_id, OR (odds ratio), pvalue, gene_id, output from porp.het, for each gene selects the min(pvalue) and add it to output
#' @param probs extremes for the posterior probability mass, defaults to 2.5 and 97.5
#' @keywords stan multiple snps
#' @export
#' @return DT with formatted data
#' stan.bt()

stan.bt <- function(x, y = "bj", rtag = NULL, model = "NB-ASE", nhets = NA, ASE.het = NA, gene, EAF = NULL, info = NULL, nfsnps = NULL, min.pval = NULL, probs = NULL) {
  ind.error <- sapply(x, class) == "try-error"
  if (all(ind.error)) {
    stop("No runs succeeded")
  }
  if (any(ind.error)) {
    x_rep <- x[[which(!ind.error)[[1]]]]
    x_rep[] <- NA
    x[ind.error] <- replicate(
      sum(ind.error), x_rep,
      simplify = FALSE
    )
  }
  if (!is.null(y)) {
    l <- lapply(x, function(i) i[y, ])
  } else {
    l <- x
  }
  DT <- data.table::data.table(do.call(rbind, l))
  ## convert to log2 all cols except n_eff and Rhat
  cols.ex <- intersect(colnames(DT), c("n_eff", "Rhat", "post.prop.neg"))
  ## new version allowing extra probs cols
  DT2 <- DT[, lapply(.SD, function(i) i / log(2)), .SDcols = names(DT)[!names(DT) %in% cols.ex]]
  DT[, names(DT)[!names(DT) %in% cols.ex] := DT2]

  ## add col for whether 95-99%CI contains the null (0)

  if (is.null(probs)) {
    ex.prob <- c(0.025, 0.975)
  } else {
    ## take the extremes of probs
    ex.prob <- probs[c(1, length(probs))]
  }

  exp1 <- paste0("`", ex.prob * 100, "%", "`")
  exp <- paste0("sign(", exp1[1], ") == sign(", exp1[2], ")")

  null <- paste("null", 100 * diff(ex.prob), sep = ".")
  DT[, eval(null) := "yes"][eval(parse(text = exp)), eval(null) := "no"]

  ## Add significance column
  DT[, Signif := "no"][get(null) == "no", Signif := "yes"]

  ## Add distance column
  DT[get(null) == "no" & eval(parse(text = exp1[1])) > 0, d := eval(parse(text = exp1[1]))]
  DT[get(null) == "no" & eval(parse(text = exp1[1])) < 0, d := -eval(parse(text = exp1[2]))]
  DT[get(null) == "yes", d := abs(eval(parse(text = exp1[1])) - eval(parse(text = exp1[2])))]



  cols.ex <- c(null, "Signif", cols.ex)
  ## add log2 to relevant cols
  data.table::setnames(DT, names(DT)[!names(DT) %in% cols.ex], paste0("log2_aFC_", names(DT)[!names(DT) %in% cols.ex]))

  DT[, tag := names(x)]
  DT[, Gene_id := gene]
  data.table::setcolorder(DT, c("Gene_id", "tag", grep("log2", names(DT), value = T), cols.ex))
  DT[, model := model]
  DT[, nhets := nhets]
  DT[, ASE.hets := ASE.het]

  ## order EAF as DT$tag
  if (!is.null(EAF)) {
    EAF <- EAF[order(match(snp, DT$tag)), ]
    DT[, tag.EAF := EAF$eaf]
  }

  if (!is.null(info)) {
    info <- info[DT$tag]
    DT[, info := info]
  }
  if (!is.null(nfsnps)) {
    DT[, n.fsnps := nfsnps]
  }
  if (!is.null(min.pval)) {
    DT[, min.p.fsnp := min(min.pval$pvalue)]
  }
  if (is.null(rtag)) {
    data.table::setnames(DT, grep("tag", names(DT), value = T), gsub("tag", "SNP", grep("tag", names(DT), value = T)))
    DT[, SNP := names(x)]
  }
  data.table::setorderv(DT, c(null, "log2_aFC_d"), order = c(1, -1))

  ## if post.prop.neg in DT (proportion of the posterior that is <0) then calculate PEP (proportion of posterior of opposite sign to bj)

  if ("post.prop.neg" %in% names(DT)) {
    DT[, PEP := post.prop.neg][log2_aFC_mean < 0, PEP := 1 - post.prop.neg]
    DT[, post.prop.neg := NULL]
  }


  return(DT)
}


#' Format a named list of stan output for baseqtl2T.nogt
#'
#' This function allows you to extract the parameter information from a list of stan summaries (output from stan.many.sim)
#' @param x list of summaries
#' @param rtag optional argument,whether snps were grouped using tag function, defaults to yes
#' @param gene gene id of gene under study
#' @param EAF, data table with snp and eaf for tag snps, output from snp.eaf (real.data.R), defaults to NULL
#' @param info, named vector with info type score for tag snp, names snp_id, defaults to NULL
#' @param nfsnps, number of fsnps in analysis, defaults to NULL
#' @param min.pval data table with cols fsnp_id, OR (odds ratio), pvalue, gene_id, output from porp.het, for each gene selects the min(pvalue) and add it to output
#' @param probs extremes for the posterior probability mass, defaults to 2.5 and 97.5
#' @keywords stan multiple snps 2 tissues
#' @export
#' @return DT with formated data
#' stan.2T()

stan.2T <- function(x, rtag = NULL, gene, EAF = NULL, info = NULL, nfsnps = NULL, min.pval = NULL, probs = NULL) {
  ind.error <- sapply(x, class) == "try-error"
  if (all(ind.error)) {
    stop("No runs succeeded")
  }
  if (any(ind.error)) {
    x_rep <- x[[which(!ind.error)[[1]]]]
    x_rep[] <- NA
    x[ind.error] <- replicate(
      sum(ind.error), x_rep,
      simplify = FALSE
    )
  }

  DT <- rbindlist(x)
  ## convert to log2
  cols <- unlist(lapply(c("mean", "sd", "%"), function(i) grep(i, names(DT), value = T)))
  DT2 <- DT[, lapply(.SD, function(i) i / log(2)), .SDcols = cols]
  DT[, (cols) := DT2]
  if (!is.null(rtag)) data.table::setnames(DT, "rSNP", "tag")
  data.table::setnames(DT, cols, paste0("log2_aFC_", cols))
  DT[, Gene_id := gene]
  data.table::setcolorder(DT, names(DT)[c(ncol(DT), 1:(ncol(DT) - 1))])
  ## order EAF as DT$tag
  if (!is.null(EAF)) {
    EAF <- EAF[order(match(snp, DT$tag)), ]
    DT[, tag.EAF := EAF$eaf]
  }

  if (!is.null(info)) {
    if (is.list(info)) {
      info <- lapply(info, function(i) i[DT$tag])
      DT[, (paste("info", names(info), sep = "_")) := info]
    }
    if (is.numeric(info)) {
      DT[, info := info]
    }
  }

  if (!is.null(nfsnps)) {
    DT[, n.fsnps := nfsnps]
  }
  if (!is.null(min.pval)) {
    DT[, min.p.fsnp := min(min.pval$pvalue)]
  }

  return(DT)
}

#' Change names of bp, bn to more general case
#'
#' Make data tble summary more general by replacing bp and bn for bt1 and bt2 when dealing with treatments/conditions
#' @param a suffix for old names
#' @param b suffix for new names
#' @param x data table to change names
#' @return data table with new names

relab <- function(a, b, x) {
  old <- grep(paste0(a, "$"), names(x), value = T)
  new <- gsub(paste0(a, "$"), b, old)
  return(data.table::setnames(x, old, new))
}
