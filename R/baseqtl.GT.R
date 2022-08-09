#' Run BaseQTL with known rsnp GT, optional refbias correction
#'
#' This function allows you to run BaseQTL for one gene and multiple pre-selected snps. When there is no enough information to ASE counts or rSNP is not in the reference panel, the function will run bayesian negative binomial model only.
#' @param prefix optional prefix for saving tables, if NULL gene_id will be used
#' @param model  whether to run NB-ASE (full model negative binomial and allele specific counts),NB (negative binomial only) or both (NB-ASE and NB for those associations with no ASE information)
#' @param stan.model compiled stanmodel object with stan model, defaults NULL to use  built-in NB-ASE model. When AI_estimate argument (below) is provided the model corrects for reference panel bias, otherwise it doesn't.
#' @param stan.negonly compiled stanmodel object with neg only side, deafults to NULL to use built-in model.
#' @param prior named list: mean= vector with the mean of Gaussians, sd= vector with Gaussians sd for eQTL effect prior, mix=vector with mixing proportions. Defaults to NULL, mixture of 2 components with mean (0,0); sd  c( 0.0309, 0.3479); and mixing proportions  c(0.97359164, 0.02640836).
#' @param ex.fsnp character vector with pos:ref:alt for fsnps to exclude,  defaults to NULL which corresponds to no exclusions
#' @param inference.method The inference method to use with Stan. Options are \code{c("sampling", "vb", "optimizing")}
#' @param snps.list A list of SNPs to restrict the number of tests performed. Useful for screening a specific subset of SNPs.
#' @param screen.method The method to used to perform screening with approximate inference.
#' @param screen.prob The probability threshold to use for screening with approximate inference.
#' @param mc.cores The number of parallel cores to use when performing tests.
#' 
#' @inheritParams aux.in1
#' @examples
#' ## example inputs saved in package directory
#' counts.f <- system.file("extdata/input", "counts.txt", package = "baseqtl", mustWork = TRUE)
#' covariates <- system.file("extdata/input", "lbsize_gc.rds", package = "baseqtl", mustWork = TRUE)
#' e.snps <- system.file("extdata/input", "chr22.fSNPS.ENSG00000159958.txt", package = "baseqtl", mustWork = TRUE)
#' u.snps <- system.file("extdata/input", "chr22.unique.fSNPS.ENSG00000159958.txt", package = "baseqtl", mustWork = TRUE)
#' gene.coord <- system.file("extdata/input", "ENSG00000159958_data.txt", package = "baseqtl", mustWork = TRUE)
#' vcf <- system.file("extdata/input", "chr22GT.86GEU.vcf.gz", package = "baseqtl", mustWork = TRUE)
#' le.file <- system.file("extdata/input", "1000GP_Phase3_subset_chr22.legend.gz", package = "baseqtl", mustWork = TRUE)
#' h.file <- system.file("extdata/input", "1000GP_Phase3_subset_chr22.hap.gz", package = "baseqtl", mustWork = TRUE)
#' AI_estimate <- system.file("extdata/input", "AI_estimate.GT.txt", package = "baseqtl", mustWork = TRUE)
#' 
#' ## output dir
#' out <- tempdir()
#' 
#' baseqtl.gt(
#'     gene = "ENSG00000159958",
#'     chr = 22,
#'     snps = 10^4,
#'     counts.f = counts.f,
#'     covariates = covariates,
#'     e.snps = e.snps,
#'     u.esnps = u.snps,
#'     gene.coord = gene.coord,
#'     vcf = vcf,
#'     le.file = le.file,
#'     h.file = h.file,
#'     out = out,
#'     inference.method = "vb",
#'     AI_estimate = AI_estimate
#' )
#' 
#' @return Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of excluded rsnps from model.
#' @export

baseqtl.gt <- function(gene, chr, snps = 5 * 10^5, counts.f, covariates = 1, additional_cov = NULL,
                       e.snps, u.esnps = NULL, gene.coord, vcf, le.file, h.file,
                       population = c("EUR", "AFR", "AMR", "EAS", "SAS", "ALL"), nhets = 5,
                       min.ase = 5, min.ase.het = 5, tag.threshold = .9, out = ".", prefix = NULL,
                       model = c("both", "NB-ASE", "NB"), stan.model = NULL, stan.negonly = NULL,
                       prob = NULL, prior = NULL, ex.fsnp = NULL, AI_estimate = NULL,
                       pretotalReads = 100, inference.method = "sampling",
                       snps.list = NULL,
                       screen.method = NULL, screen.prob = 0.5,
                       # backend = c("rstan", "cmdstanr"),
                       mc.cores = getOption("mc.cores", 1)) {

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
  
  ## check for valid stan models
  if (is.null(stan.model)) {
    ## check if ref panelbias correction
    if (is.null(AI_estimate)) {
      stan.model <- stanmodels$GT_nb_ase
    } else {
      stan.model <- stanmodels$GT_nb_ase_refbias
    }
  } else { ## model provided by user
    if (!inherits(stan.model, "stanmodel")) {
      stop("Model provided by user in argument stan.model is not a stanmodel object")
    }
  }

  if (is.null(stan.negonly)) {
    stan.negonly <- stanmodels$GT_nb
  } else {
    if (!inherits(stan.negonly, "stanmodel")) {
      stop("Model provided by user for argument stan.negonly is not a stanmodel object")
    }
  }

  ## prepare inputs
  base.in <- baseqtl.gt.in(
    gene = gene,
    chr = chr,
    snps = snps,
    counts.f = counts.f,
    covariates = covariates,
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

  if (!is.null(snps.list) & length(snps.list)) {
    nbase.snps <- names(base.in$nbase$stanIn)
    nb.snps <- names(base.in$nbase$neg)
    if (!all(snps.list) %in% c(nbase.snps, nb.snps)) {
      stop("snp.list must be a list of SNPs that could be tested.")
    }
    base.in$nbase$stanIn <- base.in$nbase$stanIn[intersect(snps.list, nbase.snps)]
    base.in$nbase$neg <- base.in$nbase$neg[intersect(snps.list, nb.snps)]
  }
  ## get inputs
  if (any(names(base.in) == "nbase")) { ## proceed with full model
    stan.in2 <- base.in$nbase$stanIn
    ASE.hets <- base.in$nbase$ASE.hets
    eaf.t <- base.in$nbase$eaf.t
    probs <- base.in$nbase$probs
    nhets <- base.in$nbase$nhets
    nfsnps <- base.in$nbase$nfsnps
    r.tag <- base.in$nbase$r.tag
    

    message("Running NB_ASE model")
    stan.full <- parallel::mclapply(stan.in2,
      function(i) {
        s <- run.stan(stan.model, data = i, pars = "bj", probs = probs, inference.method = inference.method)
        return(s)
      }, mc.cores = mc.cores
    )
    names(stan.full) <- names(stan.in2)
    full.sum <- stan.bt(x = stan.full, y = NULL, rtag = r.tag, model = "NB-ASE", nhets = nhets, ASE.het = ASE.hets, gene = gene, EAF = eaf.t, nfsnps = nfsnps, probs = probs)

    ## add min AI_post
    if (!is.null(AI_estimate)) full.sum[, min_AI := base.in$nbase$minAI]
  }

  if (any(names(base.in) == "neg")) {
    in.neg <- base.in$neg$neg
    eaf.t <- base.in$neg$eaf.t
    r.tag <- base.in$neg$r.tag
    probs <- base.in$neg$probs
    nhets <- base.in$neg$nhets

    message("Running NB model")
    ## get full posterior
    stan.neg <- parallel::mclapply(in.neg,
      function(i) {
        s <- run.stan(stan.negonly, data = i, pars = "bj", probs = probs, inference.method = inference.method)
        return(s)
      },
      mc.cores = mc.cores
    )
    names(stan.neg) <- names(in.neg)

    neg.sum <- stan.bt(x = stan.neg, y = NULL, rtag = r.tag, model = "NB", nhets = nhets, gene = gene, EAF = eaf.t, nfsnps = "NA", probs = probs)
  }

  if ((exists("full.sum") & exists("neg.sum")) | exists("neg.sum")) {
    if (exists("full.sum")) {
      neg.sum <- rbind(full.sum, neg.sum, screen.results, fill = TRUE)
    }
    if (!is.null(screen.results)) {
      neg.sum$test.type <- c(rep("full", nrow(neg.sum) - nrow(screen.results)), rep("screen", nrow(screen.results)))
    }

    if (!is.null(prefix)) {
      write.table(neg.sum, paste0(out, "/", prefix, ".GT.stan.summary.txt"), row.names = FALSE)
    } else {
      write.table(neg.sum, paste0(out, "/", gene, ".GT.stan.summary.txt"), row.names = FALSE)
    }
  }
  if (exists("full.sum") & !exists("neg.sum")) {
    if (!is.null(screen.results)) {
      full.sum <- rbind(full.sum, screen.results, fill = TRUE)
      full.sum$test.type <- c(rep("full", nrow(full.sum) - nrow(screen.results)), rep("screen", nrow(screen.results)))
    }

    if (!is.null(prefix)) {
      write.table(full.sum, paste0(out, "/", prefix, ".GT.stan.summary.txt"), row.names = FALSE)
    } else {
      write.table(full.sum, paste0(out, "/", gene, ".GT.stan.summary.txt"), row.names = FALSE)
    }
  }
  full.sum
}
