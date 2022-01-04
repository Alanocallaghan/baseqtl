
#' Run BaseQTL with known rsnp GT, optional refbias correction
#'
#' This function allows you to run BaseQTL for one gene and multiple pre-selected snps. When there is no enough information to ASE counts or rSNP is not in the reference panel, the function will run bayesian negative binomial model only.
#' @param gene gene id for the gene to run
#' @param chr chromosome where the gene is, example chr=22
#' @param snps either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window
#' @param counts.f path to file with filtered counts: rows genes, first col gene_id followed by samples, prepared in inputs.R
#' @param covariates path to matrix of covariates prepared in inputs.R, if no covariates, covariates =1, default
#' @param additional_cov full name to file with first column sample names and additional columns gene independent covariates, defaults to NULL
#' @param e.snps path to file listing exonic snps for the chromosome where the gene is, prepared in input.R
#' @param u.esnps whether to use unique exonic snps per gene, defaults to NULL when it is not necessary if strand info is known
#' @param gene.coord path to file listing gene coordinates and exons, prepared in input.R
#' @param vcf path to vcf file with ASE and GT for the chromosome where the gene is
#' @param le.file path to gz legend file (legend.gz) for the chromosome of interest for the reference panel (snp description)
#' @param h.file path to gz haplotpe file for the chromosome of interest for the reference panel (haplotypes for all samples in reference panel)
#' @param population ethnicity to get EAF for rsnp: AFR AMR EAS EUR SAS ALL, defaults to EUR
#' @param nhets minimun number of het individuals in order to run the minumn model (NB only), defaults to NULL
#' @param min.ase minimun number of ASE counts for an individual in order to be included, defaults to 5
#' @param min.ase.het minimun number of het individuals with the minimun of ASE counts in order to run the ASE side of the model, defaults to NULL
#' @param min.ase.n minimun number individuals with the minimun of ASE counts, defaults to NULL
#' @param tag.threshold numeric with r2 threshold (0-1) for grouping snps to reduce the number of running tests, to disable use "no"
#' @param out path to save outputs, default to current directory
#' @param prefix optional prefix for saving tables, if NULL gene_id will be used
#' @param model  whether to run NB-ASE (full model negative binomial and allele specific counts),NB (negative binomial only) or both (NB-ASE and NB for those associations with no ASE information)
#' @param stan.model compiled stanmodel object with stan model, defaults NULL to use  built-in NB-ASE model. When AI_estimate argument (below) is provided the model corrects for reference panel bias, otherwise it doesn't.
#' @param stan.negonly compiled stanmodel object with neg only side, deafults to NULL to use built-in model.
#' @param prob  number p∈(0,1) indicating the desired probability mass to include in the intervals, defaults to 0.99 -0.95 quantiles
#' @param prior named list: mean= vector with the mean of Gaussians, sd= vector with Gaussians sd for eQTL effect prior, mix=vector with mixing proportions. Defaults to NULL, mixture of 2 components with mean (0,0); sd  c( 0.0309, 0.3479); and mixing proportions  c(0.97359164, 0.02640836).
#' @param ex.fsnp, character vector with pos:ref:alt for fsnps to exclude,  defaults to NULL which corresponds to no exclusions
#' @param AI_estimate full name to data table with AI estimates for reference panel bias for fSNPs, defaults to NULL with no correction
#' @param pretotalReads numeric indicating a cut-off for total initial reads to consider AI estimates, defaults to 100
#' @export
#' @return Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of excluded rsnps from model.

baseqtl.gt <- function(gene, chr, snps = 5 * 10^5, counts.f, covariates = 1, additional_cov = NULL, e.snps, u.esnps = NULL, gene.coord, vcf, le.file, h.file, population = c("EUR", "AFR", "AMR", "EAS", "SAS", "ALL"), nhets = 5, min.ase = 5, min.ase.het = 5, tag.threshold = .9, out = ".", prefix = NULL, model = c("both", "NB-ASE", "NB"), stan.model = NULL,
                       stan.negonly = NULL,
                       prob = NULL, prior = NULL, ex.fsnp = NULL, AI_estimate = NULL, pretotalReads = 100, inference.method="sampling",
                       mc.cores=getOption("mc.cores", 1)) {

  ## check for valid stan models
  if (is.null(stan.model)) {
    ## check if ref panelbias correction
    if (is.null(AI_estimate)) {
      stan.model <- stanmodels$GT_nb_ase
    } else {
      stan.model <- stanmodels$GT_nb_ase_refbias
    }
  } else { ## model provided by user
    if (class(stan.model) != "stanmodel") {
      stop("Model provided by user in argument stan.model is not a stanmodel object")
    }
  }

  if (is.null(stan.negonly)) {
    stan.negonly <- stanmodels$GT_nb
  } else {
    if (class(stan.negonly) != "stanmodel") {
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

    stan.full <- parallel::mclapply(stan.in2, function(i) {
      s <- run.stan(stan.model, data = i, pars = "bj", probs = probs, method = inference.method)
      return(s)
    }, mc.cores=mc.cores)
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
    stan.neg <- parallel::mclapply(in.neg, function(i) {
      s <- run.stan(stan.negonly, data = i, pars = "bj", probs = probs, method = inference.method)
      return(s)
    }, mc.cores=mc.cores)
    names(stan.neg) <- names(in.neg)

    neg.sum <- stan.bt(x = stan.neg, y = NULL, rtag = r.tag, model = "NB", nhets = nhets, gene = gene, EAF = eaf.t, nfsnps = "NA", probs = probs)
  }

  if ((exists("full.sum") & exists("neg.sum")) | exists("neg.sum")) {
    if (exists("full.sum")) {
      neg.sum <- rbind(full.sum, neg.sum, fill = TRUE)
    }

    if (!is.null(prefix)) {
      write.table(neg.sum, paste0(out, "/", prefix, ".GT.stan.summary.txt"), row.names = FALSE)
    } else {
      write.table(neg.sum, paste0(out, "/", gene, ".GT.stan.summary.txt"), row.names = FALSE)
    }
  }



  if (exists("full.sum") & !exists("neg.sum")) {
    if (!is.null(prefix)) {
      write.table(full.sum, paste0(out, "/", prefix, ".GT.stan.summary.txt"), row.names = FALSE)
    } else {
      write.table(full.sum, paste0(out, "/", gene, ".GT.stan.summary.txt"), row.names = FALSE)
    }
  }
  full.sum
}
