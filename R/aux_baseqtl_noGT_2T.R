## Prepares inputs for running BaseQTL withut genotypes with individuals from two conditions (healthy/diseased)

#' convenient wrap for testing function per skin
#'
#' @param l list to test, list has dta tables and I am testing for rows>0 in both elements
#' @param txt character with text for error message
#' @export
#' @return character with error message
#' list.err()

list.err <- function(l, txt) {
  s <- sapply(l, function(i) nrow(i) == 0)
  if (any(s)) {
    return(paste(txt, paste(names(s[s == FALSE]), collapse = " and ")))
  } else {
    return(NULL)
  }
}


#' Third aux function to prepare inputs to run baseqtl without rsnp GT: prepare stan inputs wrap
#'
#' This function allows you to prepare stan inputs final stages
#' @param gene gene id for the gene to run
#' @param ai list with estimates for allelic imbalance per skin, defaults to NULL
#' @param case, list with data table with ASE counts per skin
#' @param rp.f list with reference panel haps for fSNPs per skin
#' @param rp.r ref panel haps for rSNPs
#' @param f.ase list with data table with genotypes for fSNPs per skin
#' @param counts.g list with data table with total counts for gene per skin
#' @param covariates per skin
#' @param min.ase minimun number of ASE counts for an individual in order to be included, defaults to 5
#' @param min.ase.n minimun number individuals with the minimun of ASE counts, defaults to 5.
#' @param info numeric cut-off for var(E(G))/var(G), var(E(G)) is the expected variance for input and var(G) for reference panel, similar to info score in impute2, defaults to 0.3. rsnps with lower info wont be run by stan.
#' @param snps.ex data table with rsnps to exclude from running
#' @param prefix character with prefix to add for saving files, defaults to NULL
#' @param out path to save outputs, default to current directory
#' @param prior named list: mean= vector with the mean of Gaussians, sd= vector with Gaussians sd for eQTL effect prior, mix=vector with mixing proportions. Defaults to NULL, mixture of 3 components with mean (0,0,0) ;sd  c( 0.0436992,  0.3492696, 0.4920049); and mixing proportions  c(0.955, 2*0.015, 0.015).
#' @param mc.cores The number of parallel cores to use.
#' @return list with stan input
#' @rdname aux-in
#' @export
aux.in <- function(gene, ai = NULL, case, rp.f, rp.r, f.ase, counts.g, covariates,
                   min.ase = 5, min.ase.n = 5, info = 0.3, snps.ex, prefix = NULL, out = ".",
                   prior = NULL, mc.cores = getOption("mc.cores", 1)) {
  if (!is.null(ai)) {
    stan.f <- mapply(function(a, b, c, d) fsnp.prep2(a, b, c, min.ase, min.ase.n, d),
      a = rp.f,
      b = f.ase,
      c = case,
      d = ai,
      SIMPLIFY = FALSE
    )
  } else {
    stan.f <- mapply(function(a, b, c) fsnp.prep2(a, b, c, min.ase, min.ase.n, ai),
      a = rp.f,
      b = f.ase,
      c = case,
      SIMPLIFY = FALSE
    )
  }

  if (!is.null(prior)) {
    if (!is.list(prior)) stop("prior argument must be a list")
    if (any(!names(prior) %in% c("mean", "sd", "mix"))) stop("prior argument must be a named 'mean' and 'sd'")
    if (length(unique(sapply(prior, length))) != 1) stop("mean and sd for prior argument must have the same length")
  } else {
    ## use default prior
    prior <- c(0, 0, 0, sqrt(2) * 0.0309, sqrt(0.0309**2 + 0.3479**2), sqrt(2) * 0.3479, 0.955, 2 * 0.015, 0.015)
    k <- length(prior) / 3 ## number of gaussians
    s <- seq(1, length(prior), k)
    l <- lapply(1:3, function(i) as.numeric(prior[s[i]:(s[i] + k - 1)]))
    names(l) <- c("mean", "sd", "mix")
    prior <- l
  }

  if (any(sapply(stan.f, is.character))) {
    return(stan.f[[sapply(stan.f, is.character)]])
  }

  stan.noGT <- lapply(
    seq_along(stan.f),
    function(j) {
      parallel::mclapply(
        1:nrow(rp.r),
        function(i) {
          stan.trecase.rna.noGT.eff2(
            counts.g[[j]],
            rp.1r = rp.r[i, , drop = FALSE],
            rp.f[[j]],
            stan.f[[j]]
          )
        },
        mc.cores = mc.cores
      )
    }
  )

  stan.noGT <- lapply(stan.noGT, setNames, rownames(rp.r))

  ## remove rsnps when genotypes of fsnps are not compatible with reference panel
  w <- lapply(stan.noGT, function(i) sapply(i, function(j) is.character(j) | is.null(j)))

  if (all(sapply(w, all))) {
    return("No inputs made for rSNPs")
  }

  if (any(lapply(w, sum) > 0)) {
    ## I need to test the same rSNPs in both datasets, so I select the rSNPs that can be run in both methods and remove the rest
    keep <- Reduce(intersect, lapply(w, function(i) names(i)[!i]))
    snps.ex <- rbind(snps.ex, data.table::data.table(treatment = "Either/Both", id = rownames(rp.r)[!rownames(rp.r) %in% keep], reason = "Failed stan input"))
    stan.noGT <- lapply(stan.noGT, function(i) i[keep])

    ## update refrence panel
    rp.r <- rp.r[rownames(rp.r) %in% keep, , drop = F]
  }

  ## restrict rsnps to those with info over cut-off. Pool individual for same cis-SNP
  ## I need to pool stan.noGT...$NB$p.g per snp and then use info.cut, make nested list with appropiate names
  pg.rsnp <- lapply(rownames(rp.r), function(i) list(NB = list(p.g = c(stan.noGT[[1]][[i]]$NB$p.g, stan.noGT[[2]][[i]]$NB$p.g))))
  names(pg.rsnp) <- rownames(rp.r)
  info.ok <- info.cut(pg.rsnp, rp.r, info)

  if (length(info.ok) < nrow(rp.r)) {
    snps.ex <- rbind(snps.ex, data.table::data.table(treatment = "Both", id = rownames(rp.r)[!rownames(rp.r) %in% names(info.ok)], reason = "Below info cut-off"))
    stan.noGT <- lapply(stan.noGT, function(i) i[names(info.ok)])
  }
  if (nrow(snps.ex)) {
    ## add gene and re-order
    snps.ex[, Gene_id := gene]
    data.table::setcolorder(snps.ex, c("Gene_id", "id", "treatment", "reason"))
    if (!is.null(prefix)) {
      write.table(snps.ex, paste0(out, "/", prefix, ".2T.noGT.excluded.snps.txt"), row.names = FALSE)
    } else {
      write.table(snps.ex, paste0(out, "/", gene, ".2T.noGT.excluded.snps.txt"), row.names = FALSE)
    }
  }

  if (length(info.ok) >= 1) { ## rsnps to run

    stan.noGT2 <- lapply(seq_along(stan.noGT), function(j) {
      if (!is.matrix(covariates[[j]])) {
        lapply(stan.noGT[[j]], function(i) in.neg.beta.noGT.eff2(i, covar = covariates[[j]][names(i$NB$counts), , drop = F]))
      } else {
        lapply(stan.noGT[[j]], function(i) in.neg.beta.noGT.eff2(i))
      }
    })
    names(stan.noGT2) <- names(f.ase)

    ## format for 2T stan input
    add <- c("N", "G", "A", "L") ## variables to add
    same <- "K" ## no change
    conc <- c("Y", "sNB", "gNB", "pNB", "gase", "m", "n", "pH", "s") ## to concatenate
    if (!is.null(AI_estimate)) {
      conc <- c("Y", "sNB", "gNB", "pNB", "gase", "m", "n", "pH", "ai0", "sdai0", "s", "h2g")
    }

    bind <- "cov" ## to bind

    stan2T <- lapply(names(info.ok), function(j) {
      l <- list()
      counter <- 1
      ## extract same snps per tissue
      l.sub <- lapply(stan.noGT2, "[[", j)
      for (a in add) {
        l[[counter]] <- Reduce("+", lapply(l.sub, function(i) i[[a]]))
        counter <- counter + 1
      }

      for (s in same) {
        l[[counter]] <- Reduce(mean, lapply(l.sub, function(i) i[[s]]))
        counter <- counter + 1
      }
      for (x in conc) {
        l[[counter]] <- Reduce(c, lapply(l.sub, function(i) i[[x]]))
        counter <- counter + 1
      }
      for (b in bind) {
        l[[counter]] <- Reduce(rbind, lapply(l.sub, function(i) i[[b]]))
        counter <- counter + 1
      }
      ## ASEi: for each element of l.sub the counter of ASE inds (col2) starts in 1. This counter is used to get "m" counts. I need to make it sequential to get correct "m" counts as "m" counts for all individuals have been concatenated in a vector. I need to add to each list the number of individuals with ASE counts (A) from previous list, 0 to the first list.
      ## make a sequence starting with 0 and get number of individuals with ASE in each list ("A" element)

      seq <- c(0, sapply(l.sub, "[[", "A"))

      ## get cumsum, last element wont be used
      seq <- cumsum(seq)
      ## add number of ind with ASE from previous list to current list only to individuals with ASE counts (1 in col 1), add 0 to the first matrix.
      for (x in length(l.sub)) {
        l.sub[[x]][["ASEi"]][, 2][l.sub[[x]][["ASEi"]][, 1] == 1] <- l.sub[[x]][["ASEi"]][, 2][l.sub[[x]][["ASEi"]][, 1] == 1] + seq[x]
      }
      ## cbind and add to l
      l[[counter]] <- Reduce(rbind, lapply(l.sub, "[[", "ASEi"))
      counter <- counter + 1


      names(l) <- c(add, same, conc, bind, "ASEi")

      ## add prior
      l <- add.prior(prior, l)

      ## add new data, psoriasis indicators

      l[["I"]] <- rep(c(1, 0), sapply(l.sub, function(i) i[["N"]]))

      l[["IA"]] <- rep(c(1, 0), sapply(l.sub, function(i) i[["A"]]))

      return(l)
    })

    names(stan2T) <- names(info.ok)

    ## get number of fsnps used in each model
    nfsnps <- Reduce(function(a, b) paste(a, b, sep = ","), lapply(stan.f, function(i) length(unlist(i$f.comb))))

    if (!is.null(ai)) {
      return(list(stan.in = stan2T, nfsnps = nfsnps, info.ok = info.ok, min_AI = unique(min(rbindlist(ai)$AI_post))))
    } else {
      return(list(stan.in = stan2T, nfsnps = nfsnps, info.ok = info.ok))
    }
  } else {
    return("None of the snps  met the conditions to be run by model")
  }
}




#' Prepare inputs for baseqtl with 2 treatments in different samples  with unknown rsnp GT and missing values for GT fsnps with optional reference panel bias correction
#'
#' This function allows you to get inputs for baseqtl for one gene and multiple pre-selected snps.
#' @inheritParams baseqtl2T.nogt
#' @export
#' @return list with 1)c.ase and 2)stan.noGT object

baseqtl2Tnogt.in <- function(gene, chr, snps = 5 * 10^5, counts.f, covariates = 1, additional_cov = NULL, e.snps, u.esnps = NULL, gene.coord, vcf, sample.file = NULL, le.file, h.file, population = c("EUR", "AFR", "AMR", "EAS", "SAS", "ALL"), maf = 0.05, min.ase = 5, min.ase.snp = 5, min.ase.n = 5, tag.threshold = .9, info = 0.3, out = ".", prefix = NULL, prior = NULL, ex.fsnp = 0.01, prob = NULL, AI_estimate = NULL, pretotalReads = 100, treatments, fishjoin = "joint") {

  ## check inputs and extract inputs for gene

  if (is.null(additional_cov)) {
    ingene <- mapply(aux.in1,
      counts.f = counts.f,
      covariates = covariates,
      vcf = vcf,
      MoreArgs = list(gene,
        chr,
        additional_cov = NULL,
        snps = snps,
        e.snps = e.snps,
        u.esnps = u.esnps,
        gene.coord = gene.coord,
        sample.file = sample.file,
        le.file = le.file,
        h.file = h.file,
        population = population,
        maf = maf,
        nhets = NULL,
        min.ase = min.ase,
        min.ase.het = NULL,
        min.ase.snp = min.ase.snp,
        min.ase.n = min.ase.n,
        tag.threshold = tag.threshold,
        info = info,
        out = out,
        model = NULL,
        prob = prob,
        AI_estimate = AI_estimate,
        pretotalReads = pretotalReads
      ),
      SIMPLIFY = F
    )
  } else {
    ingene <- mapply(aux.in1,
      counts.f = counts.f,
      covariates = covariates,
      additional_cov = additional_cov,
      vcf = vcf,
      MoreArgs = list(gene,
        chr,
        snps = snps,
        e.snps = e.snps,
        u.esnps = u.esnps,
        gene.coord = gene.coord,
        sample.file = sample.file,
        le.file = le.file,
        h.file = h.file,
        population = population,
        maf = maf,
        nhets = NULL,
        min.ase = min.ase,
        min.ase.het = NULL,
        min.ase.snp = min.ase.snp,
        min.ase.n = min.ase.n,
        tag.threshold = tag.threshold,
        info = info,
        out = out,
        model = NULL,
        prob = prob,
        AI_estimate = AI_estimate,
        pretotalReads = pretotalReads
      ),
      SIMPLIFY = F
    )
  }


  names(ingene) <- treatments

  ## get inputs for counts, covariates and probs to use later
  counts.g <- lapply(ingene, function(i) i$counts)
  covariates <- lapply(ingene, function(i) i$covariates)
  probs <- lapply(ingene, function(i) i$probs)
  population <- unique(sapply(ingene, function(i) i$population))

  if (length(population) != 1) stop("Only one population type can be selected")


  ## create dt to collect snps excluded from analysis
  snps.ex <- data.table::data.table(id = character(), reason = character())

  ## get fSNPs (feature snps or exonic snps)
  fsnps <- tryCatch(
    {
      data.table::fread(cmd = paste("grep", gene, e.snps))
    },
    error = function(e) {
      paste("No entry for gene", gene, "in", e.snps)
    }
  )



  if (isTRUE(grep("No entry", fsnps) == 1) | nrow(fsnps) == 0) { ## no fsnps
    stop(paste("No entry for gene", gene, "in", e.snps))
  }

  fsnps[, id := paste(V2, V4, V5, sep = ":")]
  ## exclude fsnps if required

  if (!is.null(ex.fsnp)) {
    fsnps <- fsnps[!id %in% ex.fsnp, ]
  }

  if (nrow(fsnps) == 0) stop(paste("No valid fSNPs after excluding ", ex.fsnp))

  ## get fsnps and extract GT and ASE, remove non-informative snps (missing or hom in all samples)
  fsnps.window <- c(min(fsnps$V2), max(fsnps$V2))
  ## get by treatments
  gt.as <- lapply(vcf, function(i) vcf_w(i, chr, st = fsnps.window[1], end = fsnps.window[2], exclude = "yes"))
  names(gt.as) <- treatments

  if (any(sapply(gt.as, is.character))) stop(gt.as[sapply(gt.as, is.character)])

  snps.ex <- rbindlist(lapply(gt.as, function(i) i$excluded[id %in% fsnps$id, ]), idcol = "treatment")

  gt.as <- lapply(gt.as, function(i) i$keep[id %in% fsnps$id])
  s <- list.err(gt.as, "No exonic snps if vcf for ")
  if (!is.null(s)) stop(s)

  gcoord <- data.table::fread(gene.coord)
  ## Extract haps for fsnps and rsnp from reference panel ########
  if (is.numeric(snps)) {
    if ("percentage_gc_content" %in% names(gcoord)) { ## newer version of gcoord from psoriasis snakefile rule geno_info using start, end, chrom, longest transcript length and GC percentage
      cis_window <- tryCatch(
        {
          gcoord[gene_id == gene & chrom == chr, .(start, end)] + c(-snps, snps)
        },
        error = function(e) {
          paste("Gene ", gene, "and chromosome", chr, "are incompatibles in gene.coord input")
        }
      )
    } else { ## old version fp gcoord
      cis_window <- tryCatch(
        {
          cl_coord(file = gene.coord, chr, gene = gene, cw = snps)
        },
        error = function(e) {
          paste("Gene ", gene, "and chromosome", chr, "are incompatibles in gene.coord input")
        }
      )
    }

    if (is.character(cis_window)) stop(cis_window)
    if (is.list(cis_window)) cis_window <- unlist(cis_window)
    rp <- haps.range(file1 = le.file, file2 = h.file, cw = cis_window, population, maf)
    if (!nrow(rp)) stop("No snps extracted from the reference panel")
    rp.r <- rp
  } else {
    pos <- as.numeric(sapply(strsplit(snps, split = ":"), `[[`, 1))
    w <- which(!is.na(pos))
    if (!length(w)) stop(cat("Invalid format for snps ", snps[w]))
    ## get gene start and end, ciswindow=0
    if ("percentage_gc_content" %in% names(gcoord)) { ## newer version of gcoord from psoriasis snakefile rule geno_info using start, end, chrom, longest transcript length and GC percentage
      st_end <- tryCatch(
        {
          gcoord[gene_id == gene & chrom == chr, .(start, end)] + rep(0, 2)
        },
        error = function(e) {
          paste("Gene ", gene, "and chromosome", chr, "are incompatibles in gene.coord input")
        }
      )
    } else {
      st_end <- tryCatch(
        {
          cl_coord(gene.coord, chr, gene = gene, cw = 0)
        },
        error = function(e) {
          paste("Gene ", gene, "and chromosome", chr, "are incompatibles in gene.coord input")
        }
      )
    }

    if (is.character(st_end)) stop(st_end)
    if (is.list(st_end)) st_end <- unlist(st_end)
    ## construct cis-window with snps, making sure to include the whole gene
    m <- min(pos) < st_end[1]
    M <- max(pos) > st_end[2]
    cis_window <- ifelse(m, min(pos), st_end[1])
    cis_window <- c(cis_window, ifelse(M, max(pos), st_end[2]))
    rp <- haps.range(file1 = le.file, file2 = h.file, cis_window, population, maf)
    if (!nrow(rp)) stop("No snps extracted from the reference panel")
    w <- which(!snps %in% rownames(rp))
    if (length(w) == length(snps)) stop(cat("None of the snps ", snps, "are in the reference panel"))
    if (length(w)) {
      rp.r <- rp[snps[-w], , drop = FALSE] ## only rsnps requested and in ref panel
      snps.ex <- rbind(snps.ex, data.table::data.table(id = snps[w], reason = rep("snp not in reference panel", length(w))))
    } else {
      rp.r <- rp[snps, , drop = FALSE]
    }
  }
  w <- which(!fsnps$id %in% rownames(rp))
  if (length(w) == length(fsnps$id)) stop("None of the fsnps are in the reference panel")
  if (length(w)) {
    snps.ex <- rbind(snps.ex, data.table::data.table(treatment = "Both", id = fsnps$id[w], reason = rep("exonic snp is not in reference panel with maf filtering", length(w))))
  }

  ## only use fsnps in ref panel (filtered by maf)
  f.ase <- lapply(gt.as, function(i) i[id %in% rownames(rp), ])

  s <- list.err(gt.as, "None of the fsnps are in the reference panel for ")
  if (!is.null(s)) stop(s)

  rp.f <- rp[rownames(rp) %in% unique(unlist(lapply(f.ase, function(i) i$id))), , drop = FALSE] ## exonic snps to work with


  ## Compute Fisher test, either jointly or by treatments

  if (!is.null(ex.fsnp)) {
    if (is.numeric(ex.fsnp)) {
      if (!is.null(fishjoin)) { # join
        f.tot <- Reduce(function(a, b) merge(a, b, by = c("CHROM", "POS", "REF", "ALT", "ID", "id"), all = T, sort = F), f.ase)
        if (population == "ALL") {
          het.all <- suppressWarnings(prop_het(f.tot, rp.f, gene))
        } else {
          sam <- data.table::fread(sample.file)
          ## get rows in sam  for population
          rows.panel <- sam[GROUP == population, which = T]
          ## samples in hap file are in columns and each sample corresponds to 2 cols
          ## I need to get for each sample in which col of hap file it is.
          cols <- sort(c(rows.panel * 2, rows.panel * 2 - 1))

          ## only use reference panel for population
          het.all <- prop_het(f.tot, rp.f[, cols, drop = F], gene)
        }
        het.f <- het.all[pvalue > ex.fsnp, ]
        if (!nrow(het.f)) {
          stop("None of the fsnps are above Fisher cut-off p-value")
        } else {
          f.ase <- lapply(f.ase, function(i) i[id %in% het.f$fsnp, ])
          rp.f <- rp.f[rownames(rp.f) %in% het.f$fsnp, , drop = FALSE]
        }
      } else { # ind skins

        ## if pop!=ALL and sample.file !=NULL, select population from reference panel to check frequency of genotypes
        if (!is.null(sample.file) & population != "ALL") {
          sam <- data.table::fread(sample.file)
          ## get rows in sam  for population
          rows.panel <- sam[GROUP == population, which = T]
          ## samples in hap file are in columns and each sample corresponds to 2 cols
          ## I need to get for each sample in which col of hap file it is.
          cols <- sort(c(rows.panel * 2, rows.panel * 2 - 1))

          ## only use reference panel for population
          het.all <- lapply(f.ase, function(i) prop_het(i, rp.f[rownames(rp.f) %in% i$id, cols, drop = F], gene))
        } else {
          het.all <- lapply(f.ase, function(i) prop_het(i, rp.f[rownames(rp.f) %in% i$id, , drop = F], gene))
        }

        het.f <- lapply(het.all, function(i) i[pvalue > ex.fsnp, ])
        s <- list.err(het.f, "None of the fsnps are above Fisher cut-off for p-value")
        if (!is.null(s)) {
          stop(s)
        } else {
          f.ase <- mapply(function(i, j) i[id %in% j$fsnp, ],
            i = f.ase,
            j = het.f,
            SIMPLIFY = FALSE
          )
          het.f <- rbindlist(het.f)
        }
        rp.f <- rp.f[rownames(rp.f) %in% unique(unlist(lapply(f.ase, function(i) i$id))), , drop = FALSE]
      }
    } else {
      f.ase <- lapply(f.ase, function(i) i[!id %in% ex.fsnp, ])
      rp.f <- rp.f[rownames(rp.f) %in% unique(unlist(lapply(f.ase, function(i) i$id))), , drop = FALSE]
    }
  } else {
    het.f <- NULL
  }
  ## save het.all

  if (exists("het.all")) {
    if (!data.table::is.data.table(het.all)) het.all <- rbindlist(het.all, idcol = "treatment")
    if (!is.null(prefix)) {
      write.table(het.all, paste0(out, "/", prefix, ".2T.fsnps.het.fisher.test.txt"), row.names = FALSE)
    } else {
      write.table(het.all, paste0(out, "/", gene, ".2T.fsnps.het.fisher.test.txt"), row.names = FALSE)
    }
  }



  ## get rSNPs
  if (tag.threshold != "no") {
    ## Group rsnps by r2, use tags function from GUESSFM
    if (nrow(rp.r) == 1) stop("Only one regulatory snp to test, please set tag.threshold='no' \n Cannot cluster one snp only")
    ## exclude snps with SD=0
    x <- t(rp.r)
    SD.x <- apply(x, 2, sd)
    w <- SD.x == 0
    if (all(w)) stop("All SNPs have zero standard deviation, please correct maf cut-off and rerun")
    if (any(w)) {
      snps.ex <- rbind(snps.ex, data.table::data.table(treatment = "Both", id = colnames(x)[w], reason = rep("Snp with zero standard deviation", sum(w))))
      x <- x[, !w, drop = FALSE]
    }

    rtag <- tag.noGT(X = x, tag.threshold = tag.threshold)
    ## save rtag as data.table
    dt <- data.table::data.table(Gene_id = gene, tag = tags(rtag), SNP = rtag@.Data)
    if (!is.null(prefix)) {
      write.table(dt, paste0(out, "/", prefix, ".2T.noGT.tags.lookup.txt"), row.names = FALSE)
    } else {
      write.table(dt, paste0(out, "/", gene, ".2T.noGT.eqtl.tags.lookup.txt"), row.names = FALSE)
    }

    ## restrict rsnp to tag snps
    rp.r <- rp.r[unique(tags(rtag)), ]
  }
  r.tag <- switch(is.numeric(tag.threshold),
    "yes"
  ) ## to output results after stan, when tag.threshold is char, returns NULL

  ## get ASE counts per individual for fSNPs in reference panel

  fsnps <- fsnps[id %in% rownames(rp.f), ]

  c.ase <- lapply(f.ase, function(i) tot.ase_counts(x = i))
  if (any(sapply(c.ase, is.character))) stop(c.ase[[sapply(c.ase, is.character)]])

  c.ase <- lapply(c.ase, function(i) filt.fsnp(i, ase = 1, min.ase.snp = 0, n = min.ase.n, rem = NULL))
  if (any(sapply(c.ase, is.character))) stop(c.ase[[sapply(c.ase, is.character)]])


  ## Only use fSNPs with AI_estimates for SNPs with sufficient reads

  if (!is.null(AI_estimate)) {
    ai <- data.table::fread(AI_estimate)
    ai <- ai[CHROM == chr & Total_pre >= pretotalReads & Keep == "yes" & AI_post != 0, ]
    ai[, id := paste(POS, REF, ALT, sep = ":")]

    fsnps <- fsnps[id %in% ai$id, ]

    if (nrow(fsnps) == 0) {
      stop("No fSNPs with allelic imbalance estimates")
    }
  } else {
    ai <- NULL
  }


  if (nrow(fsnps)) { ### option for NB, go on for ASE

    ## select SNPs in reference panel and order c.ase as in f.ase
    c.ase <- mapply(function(a, b) a[, unlist(lapply(b$id, function(i) grep(i, colnames(a), value = TRUE)))],
      a = c.ase,
      b = f.ase,
      SIMPLIFY = FALSE
    )

    if (!any(sapply(c.ase, is.character))) { ## try ASE, second check

      ## look at ufsnps, the ones to use for ASE counts
      if (!is.null(u.esnps)) {
        c.ase <- lapply(c.ase, function(i) aux.in2(gene = gene, u.esnps, case = i, ai = ai))
      }


      if (!any(sapply(c.ase, is.character))) { ## go on with ASE, third check

        ## Check if there are sufficient hets with enough ASE counts
        c.ase <- lapply(c.ase, function(i) filt.fsnp(i, ase = min.ase, min.ase.snp, n = min.ase.n, rem = NULL))


        if (!any(sapply(c.ase, is.character))) { ## go for NB only, fourth check

          ##############################  Stan inputs ###############################

          ## fsnps, pre-compute what I need because applies to every snp

          message("Preparing stan inputs")
          ## ai <- if(is.null(AI_estimate)) NULL else ai[id %in% f.ase$id,.(id,  NREF_post, NALT_post, Total_post, AI_post)]
          if (!is.null(ai)) {
            ## make all inputs with the same fSNPS, by treatments
            ai <- lapply(f.ase, function(i) ai[id %in% i$id, .(id, NREF_post, NALT_post, Total_post, AI_post)])
            c.ase <- mapply(function(a, b) a[, unlist(lapply(b$id, function(i) grep(i, colnames(a), value = T)))],
              a = c.ase,
              b = ai,
              SIMPLIFY = F
            )

            rp.f <- lapply(ai, function(i) rp.f[i$id, , drop = F])
            f.ase <- mapply(function(a, b) a[id %in% b$id, ],
              a = f.ase,
              b = ai,
              SIMPLIFY = F
            )
          }

          message(paste("Effective number of exonic SNPs:", names(f.ase), sapply(f.ase, nrow)))

          inp <- aux.in(gene, ai, case = c.ase, rp.f, rp.r, f.ase, counts.g, covariates, min.ase, min.ase.n, info = info, snps.ex, prefix, out, prior)

          if (is.character(inp)) stop(inp)
          return(list(inp = inp, model = "NB-ASE", het.f = het.f, probs = probs, r.tag = r.tag, population = population))
        }
      }
    }
  } else {
    return("No ASE for runing model")
  }
}
