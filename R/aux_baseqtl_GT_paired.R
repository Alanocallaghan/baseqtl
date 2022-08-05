## Prepares inputs for running BaseQTL with paired samples

#' Aux for formatting 2 tretmenat inputs
#' @param l list with stan input for 2 treatments for 1 snps
#' @param dt data table to use for formatting
#' @param n name of elements on l to join
#' @export
#' @return data table with inputs for stan
#' treat.form()

treat.form <- function(l, dt, n) {
  tmp <- Reduce(
    function(a, b) merge(a, b, by = c("samples", "index"), all = T),
    mapply(function(x, y) {
      tmp2 <- data.table::data.table(samples = rep(dt$samples, dt[[x]]))
      tmp2[, eval(n) := y[[n]]]
      tmp2[, index := 1:.N, by = samples]
      return(tmp2)
    },
    x = names(dt)[names(dt) != "samples"],
    y = l,
    SIMPLIFY = FALSE
    )
  )
  ## change NA to 0
  tmp[is.na(tmp)] <- 0
  tmp[, index := NULL]
  return(tmp)
}

#' Aux for returning genotype info
#' @param dt data table with inputs
#' @export
#' @return vector to return to stan
#' g.for()

g.for <- function(dt) {
  n <- names(dt)[!names(dt) %in% "samples"]
  dt[, g := get(n[1])][get(n[1]) == 0, g := get(n[2])]
  return(dt$g)
}



#' Make inputs for running  baseqtl with known rsnp GT and paired design, optional refbias correction
#'
#' This function allows you to run baseqtl for one gene and multiple pre-selected snps. When there is no enough information to ASE counts or rSNP is not in the reference panel, the function will run bayesian negative binomial model only by default.
#' @export
#' @return data.table with summary of gene-snp associations. Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of excluded rsnps.
#' @inheritParams baseqtl.gt
btrecase.gt.paired.in <- function(gene, chr, snps = 5 * 10^5, counts.f, covariates = 1, additional_cov = NULL, e.snps, u.esnps = NULL, gene.coord, vcf, le.file, h.file, population = c("EUR", "AFR", "AMR", "EAS", "SAS", "ALL"), nhets = 5, min.ase = 5, min.ase.het = 5, tag.threshold = .9, out = ".", prefix = NULL, model = c("both", "NB-ASE", "NB"),
                                  prob = NULL, prior = NULL, ex.fsnp = NULL, AI_estimate = NULL, pretotalReads = 100) {


  ## check inputs and extract inputs for gene
  model <- tryCatch(match.arg(model), error = function(e) {
    e
  })

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


  ingene <- lapply(1:2, function(i) {
    aux.in1(gene,
      chr,
      snps,
      counts.f[i],
      covariates,
      additional_cov,
      e.snps,
      u.esnps,
      gene.coord,
      vcf[i],
      sample.file = NULL, ## only for noGT
      le.file,
      h.file,
      population = population,
      maf = NULL,
      nhets = nhets,
      min.ase = min.ase,
      min.ase.het = min.ase.het,
      min.ase.snp = NULL,
      min.ase.n = NULL,
      tag.threshold = tag.threshold,
      out = out,
      model = model,
      prob = prob,
      AI_estimate = AI_estimate,
      pretotalReads = pretotalReads
    )
  })

  counts.g <- lapply(ingene, "[[", "counts")
  covariates <- ingene[[1]]$covariates
  probs <- ingene[[1]]$probs
  population <- ingene[[1]]$population

  rsnps.ex <- data.table::data.table(id = character(), reason = character())
  gcoord <- data.table::fread(gene.coord)

  ## get gt:ase counts
  gt.rs <- lapply(vcf, function(i) aux.snps(gene, chr, snps, gcoord, gene.coord, i, rsnps.ex))
  rs <- gt.rs[[1]]$rs ## same in both
  gt.as <- lapply(gt.rs, "[[", "gt.as")
  rsnps.ex <- gt.rs[[1]]$rsnps.ex
  cis_window <- gt.rs[[1]]$cis

  ## QC vcf, make sure GT is the same in both files
  ## same names
  if (!Reduce(identical, lapply(gt.as, names))) stop("Names in vcf are not the same")
  ## check GT
  gt.cols <- names(gt.as[[1]])[-grep("_AS", names(gt.as[[1]]))]
  if (!Reduce(identical, lapply(gt.as, function(i) i[, gt.cols, with = F]))) stop("SNPs and/or genotypes are not the same in both vcf files")

  ## recode to 0,1,-1,2 scale
  rec.rs <- rec_mytrecase_rSNPs(x = rs$POS, y = rs)

  ## tagging
  if (tag.threshold != "no") {
    ## remove rsnps with var0, not possible to tag
    cols <- grep("_GT$", names(rec.rs), value = T)
    v <- apply(rec.rs[, cols, with = F], 1, function(i) var(abs(i)))
    w <- which(v == 0)
    if (length(w) == nrow(rec.rs)) stop("all SNPs have zero variance, cannot tag")
    pre <- ifelse(is.null(prefix), paste0(gene, ".paired.GT"), prefix)
    id.tags <- help.tags(gene, tag.threshold, rs, pre, rec.rs, out)

    if (length(w)) {
      rsnps.ex <- rbind(rsnps.ex, data.table::data.table(id = rec.rs[w, id], reason = rep("Snp with zero variance", length(w))))
    }
    rec.rs <- rec.rs[id %in% id.tags, ]
  }
  r.tag <- switch(is.numeric(tag.threshold),
    "yes"
  ) ## to output results after stan, when tag.threshold is char, returns NULL
  ## help with hets
  GT.aux <- rec.rs[, grep("_GT", names(rec.rs)), with = F] ## to make easier calculation of correlations, etc

  ## counts number of hets per rsnp
  rec.rs[, nhet := apply(GT.aux, 1, function(i) sum(abs(i) == 1))]
  ## remove snps with less than min hets
  w <- rec.rs[nhet < nhets, which = TRUE]
  rsnps.ex <- rbind(rsnps.ex, data.table::data.table(id = rec.rs[w, id], reason = rep(paste("rsnp with less than", nhets, "het ind."), length(w))))
  rec.rs <- rec.rs[!w, ]
  if (nrow(rec.rs) == 0) stop(cat("No rsnp with at least", nhets, "het ind."))

  if (model == "NB-ASE" | model == "both") {

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
      cat(fsnps, "\n No fSNPS, analysis will be done with total gene counts only")
      rsnps.ex <- rbind(rsnps.ex, data.table::data.table(id = rec.rs$id, reason = paste("No entry for gene", gene, "in", e.snps)))
    } else {
      fsnps[, id := paste(V2, V4, V5, sep = ":")]
      fsnps <- fsnps[id %in% gt.as[[1]]$id] # same in both treatments
      check1 <- nrow(fsnps)

      ## exclude fsnps if required
      if (!is.null(ex.fsnp)) {
        fsnps <- fsnps[!id %in% ex.fsnp, ]
      }
      check2 <- nrow(fsnps)

      ## Only use fSNPs with AI_estimates for SNPs with sufficient reads

      if (!is.null(AI_estimate)) {
        ai <- data.table::fread(AI_estimate)
        ai <- ai[CHROM == chr & Total_pre >= pretotalReads & Keep == "yes" & AI_post != 0, ]
        ai[, id := paste(POS, REF, ALT, sep = ":")]

        fsnps <- fsnps[id %in% ai$id, ]
      } else {
        ai <- NULL
      }


      if (nrow(fsnps)) { ## stay

        if (is.character(snps)) { ## analysis on pre-specified set of rsnps
          gt.as <- lapply(gt.as, function(i) i[id %in% c(rec.rs$id, fsnps$id), ]) ## snps and fsnps only
          ## make sure gt.as (here for fsnps) doesn't have missing values or unphased GT
          gt.qc <- lapply(gt.as, function(i) vcf.gt.qc(i))
          s <- sapply(gt.qc, function(i) sum(i[c(2, 4)]) != 0)
          if (any(s)) stop(cat("Invalid GT field for some exonic snps and samples \n", paste(names(gt.qc[[1]]), collapse = ","), "\n", gt.qc[s], "\n", "to remove snps or samples with wrong GT use vcf.gt.qc with appropiate arguments \n", "if all exonic snps are removed the analysis will be done with total counts only"))
        }

        ## get info from reference panel
        ## matrix for reference panel haps
        rp <- haps.range(file1 = le.file, file2 = h.file, cis_window, maf = 0)
        ## keep record of rsnps not in reference panel
        w <- which(!rec.rs$id %in% rownames(rp))
        rsnps.ex <- rbind(rsnps.ex, data.table::data.table(id = rec.rs$id[w], reason = rep("No rsnp in reference panel", length(w))))
        rec.rs2 <- rec.rs[!w, ] ## for full model, only rsnps in ref panel
        w <- which(rownames(rp) %in% rec.rs2$id)
        if (length(w) != 0) { ## proceed with full model, neg only for rsnps not in ref panel
          w <- which(rownames(rp) %in% gt.as[[1]]$id)
          rp <- rp[w, , drop = FALSE]
          ## make sure to select fsnps from reference panel
          f.ase <- lapply(gt.as, function(i) i[id %in% fsnps$id, ][id %in% rownames(rp), ])
          if (nrow(f.ase[[1]]) == 0) {
            print("No fsnps ref panel")
            rsnps.ex <- rbind(rsnps.ex, data.table::data.table(id = rec.rs2$id, reason = "No fsnps in reference panel"))
          } else {

            ## get counts per hap for fsnp
            counts.ase <- lapply(f.ase, tot.ase_counts)

            if (!is.null(u.esnps)) {
              counts.ase <- lapply(counts.ase, function(i) aux.in2(gene, u.esnps, i, ex.fsnp, ai))
            }

            if (!all(sapply(counts.ase, function(i) is.character(i)))) { ## enough fsnps for ASE

              ## check if sufficient indiviuals with ASE counts to proceed
              s <- sapply(counts.ase, function(i) is.character(filt.fsnp(i, min.ase.snp = 0)))
              if (!all(s)) { ## Enough ind with ase counts

                print(paste("Effective number of fSNPs:", nrow(f.ase[[1]])))

                ## select rsnps in ref panel for full model
                rs.full <- rs[id %in% rec.rs2$id, which(names(rs) %in% names(rec.rs2)), with = F]

                ###################  run stan full model #######################
                ###### prepare stan inputs
                ## reference panel fsnps and rsnps:
                rp.f <- rp[f.ase[[1]]$id, , drop = FALSE]
                rp.r <- rp[rs.full$id, , drop = FALSE]

                ## get haplotype counts from fsnps

                ## add AI_estimate per sample to stan inputs
                if (!is.null(ai)) {
                  ai <- ai[id %in% f.ase[[1]]$id, .(id, NREF_post, NALT_post, Total_post, AI_post)]
                  ## get min AI_post
                  min_AI <- min(ai$AI_post)
                }


                print("Preparing stan inputs")

                ## order fsnps in couts.ase as in f.ase

                counts.ase <- mapply(function(a, b) a[, unlist(lapply(b$id, function(i) grep(i, colnames(a), value = TRUE)))],
                  a = counts.ase,
                  b = f.ase,
                  SIMPLIFY = F
                )

                stan.f <- mapply(function(a, b) stan.fsnp.noGT.eff(rp.f, a, b, NB = "no", min.ase, min.ase.het, ai = ai),
                  a = f.ase,
                  b = counts.ase,
                  SIMPLIFY = F
                )


                if (!any(sapply(stan.f, is.character))) { ## proceed with full model

                  counts <- lapply(counts.g, unlist) ## to avoid repeating in mclapply

                  stan.in1 <- mapply(function(a, b, c) parallel::mclapply(rs.full$id, function(i) stan.trecase.eff2(a, rp.1r = rp.r[i, , drop = FALSE], rp.f, b, rs.hap = rs.full[id == i, ], rec.rsnp = rec.rs2[id == i, ], c, min.ase, min.ase.het)),
                    a = counts,
                    b = f.ase,
                    c = stan.f,
                    SIMPLIFY = F
                  )

                  stan.in1 <- lapply(stan.in1, setNames, rs.full$id)


                  ## remove rsnps with insufficient ase, min.ase.hets (haplotypes with fsnps not compatible with reference panel)
                  w <- lapply(stan.in1, function(i) sapply(i, is.character))
                  if (any(lapply(w, sum) > 0)) {

                    ## need to discard any rSNP that failed in any treatment
                    failed <- Reduce(union, lapply(w, function(a) rs.full$id[a]))

                    rsnps.ex <- rbind(rsnps.ex, data.table::data.table(id = failed, reason = "Not enough het individuals with sufficient ase counts or Hap not in reference panel"))

                    stan.in1 <- lapply(stan.in1, function(i) i[!rs.full$id %in% failed])
                  }
                  if (length(stan.in1[[1]]) >= 1) { ## same length both elements

                    stan.in2 <- lapply(stan.in1, function(j) {
                      lapply(j, function(i) {
                        l <- in.neg.beta.prob.eff2(i, covar = covariates)
                        ## add prior
                        l <- add.prior(prior, l)
                        return(l)
                      })
                    })


                    ## need to combine inputs from 2 treatments
                    ## need to identify overlapping/distinct samples
                    all <- Reduce(union, lapply(stan.f, function(i) names(i$m)))
                    ## make indicator to link individuals with ASE counts
                    ASEI <- do.call(data.table, lapply(stan.f, function(i) as.numeric(all %in% names(i$m))))
                    ## add sample names to ASEI
                    ASEI[, samples := all]

                    ## Combine inputs for each rsnp between the 2 treatments
                    stan2t <- lapply(names(stan.in2[[1]]), function(i) {
                      ## extract same snp per treatment
                      l.sub <- lapply(stan.in2, "[[", i)
                      ## ASE side
                      ## extract s and make data table to keep track of same individuals
                      s.dt <- treat.form(l.sub, ASEI, "s")
                      ## same with gase and pH
                      gase.dt <- treat.form(l.sub, ASEI, "gase")
                      pH.dt <- treat.form(l.sub, s.dt, "pH")
                      ## prepare list with elements to return as inputs
                      l <- mapply(function(l, dt, n) {
                        tmp <- treat.form(l, dt, n)
                        tmp[, samples := NULL]
                      },
                      dt = list(ASEI, s.dt),
                      n = c("m", "n"),
                      MoreArgs = list(l = l.sub),
                      SIMPLIFY = FALSE
                      )

                      names(l) <- c("m", "n")
                      if (!is.null(ai)) {
                        for (x in c("ai0", "sdai0")) {
                          l[[x]] <- treat.form(l.sub, s.dt, x)[, samples := NULL]
                        }
                      }

                      ## to return gase, pH and s
                      for (x in c("gase", "pH", "s")) {
                        l[[x]] <- g.for(get(paste0(x, ".dt")))
                      }
                      l[["A"]] <- length(all) ## union of ind with ASE counts across treatments
                      l[["L"]] <- length(l$pH)

                      ## NB side and same elements in l.sub
                      l[["Y"]] <- Reduce(cbind, lapply(l.sub, function(i) i$Y))

                      same <- c("N", "K", "g", "cov", "k", "aveP", "sdP", "mixP")
                      for (x in same) {
                        l[[x]] <- l.sub[[1]][[x]]
                      }
                      return(l)
                    })

                    names(stan2t) <- names(stan.in1[[1]])

                    ## count number of hets with sufficient ase counts to input in final report
                    ASE.hets <- Reduce(
                      function(a, b) paste(a, b, sep = ","),
                      lapply(stan.in1, function(j) sapply(j, function(i) nrow(i$gm[abs(g.ase) == 1, ])))
                    )
                    ## get eaf for tags run in model
                    eaf.t <- snp.eaf(le.file, names(stan2t), population)

                    trecase.in <- list(stanIn = stan2t, ASE.hets = ASE.hets, eaf.t = eaf.t, probs = probs, r.tag = r.tag, nhets = rec.rs2[id %in% names(stan2t), nhet], nfsnps = nchar(unlist(strsplit(colnames(stan.f[[1]]$n[[1]]), ","))[1]))

                    if (!is.null(AI_estimate)) trecase.in[["minAI"]] <- min_AI
                  } ## closing from failed stan.in1
                } else { ## all rsnps will fail because stan.f failed
                  rsnps.ex <- rbind(rsnps.ex, data.table::data.table(id = rs.full$id, reason = stan.f))
                } ## closing from failed stan.f
              } else { ## not Enough inds with ase counts

                w <- which(!rec.rs2$id %in% rsnps.ex$id)
                if (length(w)) rsnps.ex <- rbind(rsnps.ex, data.table::data.table(id = rec.rs2$id[w], reason = "No enough inds with ASE counts"))
              }
            } else { ##  no unique fSNPs, enough fsnps for ASE

              rsnps.ex <- rbind(rsnps.ex, data.table::data.table(id = rec.rs2$id, reason = "Not unique fsnps in gene"))
            } ## closing no unique fSNPs
          } ## closing from no fsnps in ref panel
        } ## closing from rsnps in ref panel
      } else { ## closing from no AI estimates
        ## go to NB if model ==both, otherwise finish
        if (!length(check1)) {
          cat("No fSNPs with genotypes in vcf file")
        } else if (!length(check2)) {
          cat("All fSNPs were excluded after ex.fsnp argument")
        } else {
          cat("No fSNPs with allelic imbalance estimates")
        }
        if (model == "both") {
          cat("\n analysis will be done with total gene counts only")
          rsnps.ex <- rbind(rsnps.ex, data.table::data.table(id = rec.rs$id, reason = "No fSNPs with GT or AI estimates"))
        }
      }
    } ## closing from fsnps in gene
    if (model == "NB-ASE") {
      if (nrow(rsnps.ex)) {
        if (!is.null(prefix)) {
          write.table(rsnps.ex, paste0(out, "/", prefix, ".paired.GT.excluded.rsnps.txt"), row.names = FALSE)
        } else {
          write.table(rsnps.ex, paste0(out, "/", gene, ".paired.GT.excluded.rsnps.txt"), row.names = FALSE)
        }
      }
      if (exists("stan2t")) {
        return(list(nbase = trecase.in))
      } else {
        return("No rsnps were run with NB-ASE")
      }
    }
  } else { ## when model==NB I add all rec.rs snps to prepare input

    rsnps.ex <- rbind(rsnps.ex, data.table::data.table(id = rec.rs$id, reason = "requested model NB"))
  }

  ## neg.binom only
  ## prepare input and run
  if (nrow(rsnps.ex) > 0) {
    ex <- c("Missing or homo GT all samples", "zero variance", "rsnp with less than")
    w <- unlist(sapply(ex, function(i) grep(i, rsnps.ex$reason)))
    id.keep <- unique(rsnps.ex[!w, id])

    if (length(id.keep)) {
      in.neg <- lapply(counts.g, function(j) parallel::mclapply(id.keep, function(i) input.neg.only.bj(DT1 = j, DT2 = rec.rs[id == i, ], covar = covariates)))

      in.neg <- lapply(in.neg, setNames, id.keep)

      neg2t <- lapply(id.keep, function(i) {
        ## extract  same snps
        l.sub <- lapply(in.neg, "[[", i)
        ## NB side and same elements in l.sub
        l <- list()
        l[["Y"]] <- Reduce(cbind, lapply(l.sub, function(i) i$Y))

        same <- c("N", "K", "g", "cov")

        for (x in same) {
          l[[x]] <- l.sub[[1]][[x]]
        }

        ## add prior
        l <- add.prior(prior, l)
        return(l)
      })
      names(neg2t) <- id.keep

      ## get eaf for tags run in model
      eaf.t.neg <- snp.eaf(le.file, names(neg2t), population)
      ## if tag not in ref panel (possible with neg binom side)
      if (length(id.keep) > nrow(eaf.t.neg)) {
        eaf.t.neg <- rbind(eaf.t.neg, data.table::data.table(snp = id.keep[which(!id.keep %in% eaf.t.neg$snp)], eaf = NA))
        ## sort as in id.keep
        eaf.t.neg <- eaf.t.neg[order(match(snp, id.keep))]
      }
      if (!is.null(prefix)) {
        write.table(rsnps.ex, paste0(out, "/", prefix, ".paired.GT.excluded.rsnps.txt"), row.names = FALSE)
      } else {
        write.table(rsnps.ex, paste0(out, "/", gene, ".paired.GT.excluded.rsnps.txt"), row.names = FALSE)
      }
      neg.in <- list(neg = neg2t, eaf.t = eaf.t.neg, r.tag = r.tag, probs = probs, nhets = rec.rs[id %in% id.keep, nhet])
      if (exists("stan2t")) {
        return(list(nbase = trecase.in, neg = neg.in))
      } else {
        return(list(neg = neg.in))
      }
    }
  } else {
    if (model == "NB") {
      return("No rsnps can be run with NB")
    }

    if (exists("stan2t")) {
      if (!is.null(prefix)) {
        write.table(rsnps.ex, paste0(out, "/", prefix, ".paired.GT.excluded.rsnps.txt"), row.names = FALSE)
      } else {
        write.table(rsnps.ex, paste0(out, "/", gene, "paired.GT.excluded.rsnps.txt"), row.names = FALSE)
      }
      return(list(nbase = trecase.in))
    } else {
      return("No snps can be run with NB-ASE or NB")
    }
  }
  if (!exists("stan2t")) {
    if (nrow(rsnps.ex)) {
      if (!is.null(prefix)) {
        write.table(rsnps.ex, paste0(out, "/", prefix, ".paired.GT.excluded.rsnps.txt"), row.names = FALSE)
      } else {
        write.table(rsnps.ex, paste0(out, "/", gene, ".paired.GT.excluded.rsnps.txt"), row.names = FALSE)
      }
    }
    return("No snps were run")
  }
}
