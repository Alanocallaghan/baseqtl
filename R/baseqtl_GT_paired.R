
#' Run BaseQTL with known rsnp GT, paired treatments optional refbias correction
#'
#' This function allows you to run BaseQTL with paired design for one gene and multiple pre-selected snps. When there is no enough information to ASE counts or rSNP is not in the reference panel, the function will run bayesian negative binomial model only.
#' @inheritParams baseqtl.gt
#' @examples
#' ## Retrieve input files for running baseqtl.gt (paired)
#' counts.f <- rep(
#'     system.file(
#'         "extdata/input", "counts.txt",
#'         package = "baseqtl", mustWork = TRUE
#'     ),
#'     2
#' )
#' covariates <- system.file(
#'     "extdata/input", "lbsize_gc.rds",
#'     package = "baseqtl", mustWork = TRUE
#' )
#' e.snps <- system.file(
#'     "extdata/input", "chr22.fSNPS.ENSG00000159958.txt",
#'     package = "baseqtl", mustWork = TRUE
#' )
#' u.snps <- system.file(
#'     "extdata/input", "chr22.unique.fSNPS.ENSG00000159958.txt",
#'     package = "baseqtl", mustWork = TRUE
#' )
#' gene.coord <- system.file(
#'     "extdata/input", "ENSG00000159958_data.txt",
#'     package = "baseqtl", mustWork = TRUE
#' )
#' vcf <- rep(
#'     system.file(
#'         "extdata/input", "chr22GT.86GEU.vcf.gz",
#'         package = "baseqtl", mustWork = TRUE
#'     ),
#'     2
#' )
#' le.file <- system.file(
#'     "extdata/input", "1000GP_Phase3_subset_chr22.legend.gz",
#'     package = "baseqtl", mustWork = TRUE
#' )
#' h.file <- system.file(
#'     "extdata/input", "1000GP_Phase3_subset_chr22.hap.gz",
#'     package = "baseqtl", mustWork = TRUE
#' )
#' AI_estimate <- system.file(
#'     "extdata/input", "AI_estimate.GT.txt",
#'     package = "baseqtl", mustWork = TRUE
#' )
#'
#' ## Choose your output directory
#' out <- tempdir()
#'
#' ## To minimise file sizes and computational load SNPs are within a 10^4 cis-window
#' ## baseqtl is computational intense, it is recommended to be run with several cores
#' ## cores are automatically detected by R
#'
#' \dontrun{
#' ## Run baseqtl.gt.paired:
#' baseqtl.gt.paired(
#'   gene = "ENSG00000159958",
#'   chr = 22,
#'   snps = 10^3,
#'   counts.f = counts.f,
#'   covariates = covariates,
#'   e.snps = e.snps,
#'   u.esnps = u.snps,
#'   gene.coord = gene.coord,
#'   vcf = vcf,
#'   le.file = le.file,
#'   h.file = h.file,
#'   out = out,
#'   AI_estimate = AI_estimate
#' )
#' }
#' @return Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of excluded rsnps.
#' @export
baseqtl.gt.paired <- function(gene, chr, snps = 5 * 10^5, counts.f, covariates = 1, additional_cov = NULL,
                              e.snps, u.esnps = NULL, gene.coord, vcf, le.file, h.file,
                              population = c("EUR", "AFR", "AMR", "EAS", "SAS", "ALL"), nhets = 5,
                              min.ase = 5, min.ase.het = 5, tag.threshold = .9, out = ".", prefix = NULL,
                              model = c("both", "NB-ASE", "NB"), stan.model = NULL, stan.negonly = NULL,
                              prob = NULL, prior = NULL, ex.fsnp = NULL, AI_estimate = NULL,
                              pretotalReads = 100, mc.cores = getOption("mc.cores", 1),
                              inference.method = c("sampling", "vb")) {

  if (!is.null(screen.method)) {
    screen.method <- match.arg(screen.method, choices = "vb")
    call <- match.call()
    call[["inference.method"]] <- call[["screen.method"]]
    call[["screen.method"]] <- NULL
    call[["prob"]] <- screen.prob
    screen.results <- eval(call, parent.frame())
    sig <- screen.results[[sprintf("null.%s", screen.prob)]] == "no"
    # assocs.test <- screen.results[sig, c("Gene_id", "tag")]
    # snps.list <- screen.results[sig, c("Gene_id", "tag")]
    snps.list <- screen.results[["tag"]][sig]
  } else {
    screen.results <- NULL
  }


  inference.method <- match.arg(inference.method)
  ## check for valid stan models
  if (is.null(stan.model)) {
    ## check if ref panelbias correction
    if (is.null(AI_estimate)) {
      stan.model <- stanmodels$paired_GT_nb_ase
    } else {
      stan.model <- stanmodels$paired_GT_nb_ase_refbias
    }
  } else { ## model provided by user
    if (!inherits(stan.model, "stanmodel")) {
      stop("Model provided by user in argument stan.model is not a stanmodel object")
    }
  }

  if (is.null(stan.negonly)) {
    stan.negonly <- stanmodels$paired_GT_nb
  } else {
    if (!inherits(stan.negonly, "stanmodel")) {
      stop("Model provided by user for argument stan.negonly is not a stanmodel object")
    }
  }

  base.in <- btrecase.gt.paired.in(gene,
    chr,
    snps,
    counts.f,
    covariates,
    additional_cov = additional_cov,
    e.snps,
    u.esnps,
    gene.coord,
    vcf,
    le.file,
    h.file,
    population,
    nhets,
    min.ase,
    min.ase.het,
    tag.threshold,
    out,
    prefix,
    model,
    prob,
    prior,
    ex.fsnp,
    AI_estimate,
    pretotalReads
  )

  if (is.character(base.in)) stop(base.in)
  if (!length(base.in)) stop("No SNPs in cis-window!")

  param <- c("ba", "bd", "bp", "bn")
  ## get inputs
  if (any(names(base.in) == "nbase")) { ## proceed with full model
    stan.in2 <- base.in$nbase$stanIn
    ASE.hets <- base.in$nbase$ASE.hets
    eaf.t <- base.in$nbase$eaf.t
    probs <- base.in$nbase$probs
    nhets <- base.in$nbase$nhets
    nfsnps <- base.in$nbase$nfsnps
    r.tag <- base.in$nbase$r.tag
    minAI <- base.in$nbase$minAI

    message("Running NB-ASE model")
    stan.full <- parallel::mclapply(
      1:length(stan.in2),
      function(i) {
        s <- run.stan(stan.model, data = stan.in2[[i]], pars = param, probs = probs, inference.method = inference.method)
        dt.wide <- stan.paired.for(s, names(stan.in2)[i])
        return(dt.wide)
      },
      mc.cores = mc.cores
    )

    full.sum <- stan.2T(
      x = stan.full, rtag = r.tag, gene = gene, EAF = eaf.t,
      nfsnps = nfsnps,
      probs = probs
    )
    full.sum[, model := "NB-ASE"][, nhets := nhets][, ASE.hets := ASE.hets]

    ## relabel "bp", "bn" to "bt1", "bt2"
    m <- mapply(relab,
      a = c("bp", "bn"),
      b = c("bt1", "bt2"),
      MoreArgs = list(full.sum),
      SIMPLIFY = F
    )

    ## add min AI_post
    if (!is.null(AI_estimate)) full.sum[, min_AI := minAI]
  }

  if (any(names(base.in) == "neg")) {
    in.neg <- base.in$neg$neg
    eaf.t <- base.in$neg$eaf.t
    r.tag <- base.in$neg$r.tag
    probs <- base.in$neg$probs
    nhets <- base.in$neg$nhets

    message("Running NB model")

    stan.neg <- parallel::mclapply(
      1:length(in.neg),
      function(i) {
        s <- run.stan(stan.negonly, data = in.neg[[i]], pars = param, probs = probs, inference.method = inference.method)
        dt.wide <- stan.paired.for(s, names(in.neg)[i])
        return(dt.wide)
      },
      mc.cores = mc.cores
    )

    neg.sum <- stan.2T(x = stan.neg, rtag = r.tag, gene = gene, EAF = eaf.t, nfsnps = "NA", probs = probs)
    neg.sum[, model := "NB"][, nhets := nhets]
    n <- mapply(relab,
      a = c("bp", "bn"),
      b = c("bt1", "bt2"),
      MoreArgs = list(neg.sum),
      SIMPLIFY = F
    )
  }





  if ((exists("full.sum") & exists("neg.sum")) | exists("neg.sum")) {
    if (exists("full.sum")) {
      neg.sum <- rbind(full.sum, neg.sum, fill = TRUE)
    }
    if (!is.null(screen.results)) {
        neg.sum <- rbind(neg.sum, screen.results, fill = TRUE)
        neg.sum$test.type <- c(rep("full", nrow(neg.sum) - nrow(screen.results)), rep("screen", nrow(screen.results)))
    }


    if (!is.null(prefix)) {
      write.table(neg.sum, paste0(out, "/", prefix, ".paired.GT.stan.summary.txt"), row.names = FALSE)
    } else {
      write.table(neg.sum, paste0(out, "/", gene, ".paired.GT.stan.summary.txt"), row.names = FALSE)
    }
  }


  if (exists("full.sum") & !exists("neg.sum")) {
    if (!is.null(screen.results)) {
        full.sum <- rbind(full.sum, screen.results, fill = TRUE)
        full.sum$test.type <- c(rep("full", nrow(full.sum) - nrow(screen.results)), rep("screen", nrow(screen.results)))
    }
    if (!is.null(prefix)) {
      write.table(full.sum, paste0(out, "/", prefix, ".paired.GT.stan.summary.txt"), row.names = FALSE)
    } else {
      write.table(full.sum, paste0(out, "/", gene, ".paired.GT.stan.summary.txt"), row.names = FALSE)
    }
  }
  full.sum
}
