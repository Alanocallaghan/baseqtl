
#' Run BaseQTL with known rsnp GT, paired treatments optional refbias correction
#'
#' This function allows you to run BaseQTL with paired design for one gene and multiple pre-selected snps. When there is no enough information to ASE counts or rSNP is not in the reference panel, the function will run bayesian negative binomial model only.
#' @param gene gene id for the gene to run
#' @param chr chromosome where the gene is, example chr=22
#' @param snps either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window
#' @param counts.f path to files with filtered counts for each treament: rows genes, first col gene_id followed by samples, prepared in inputs.R, for each treatment
#' @param covariates path to matrix of covariates prepared in inputs.R, if no covariates, covariates =1, default
#' @param additional_cov full name to file with first column sample names and additional columns gene independent covariates, defaults to NULL
#' @param e.snps path to file listing exonic snps for the chromosome where the gene is, prepared in input.R
#' @param u.esnps whether to use unique exonic snps per gene, defaults to NULL when it is not necessary if strand info is known
#' @param gene.coord path to file listing gene coordinates and exons, prepared in input.R
#' @param vcf path to vcf files with ASE and GT for the chromosome where the gene is, for each treat,ment, same order as counts.f
#' @param le.file path to gz legend file (legend.gz) for the chromosome of interest for the reference panel (snp description)
#' @param h.file path to gz haplotpe file for the chromosome of interest for the reference panel (haplotypes for all samples in reference panel)
#' @param population ethnicity to get EAF for rsnp: AFR AMR EAS EUR SAS ALL, defaults to EUR
#' @param nhets minimun number of het individuals in order to run the minumn model (NB only), defaults to NULL
#' @param min.ase minimun number of ASE counts for an individual in order to be included, defaults to 5
#' @param min.ase.het minimun number of het individuals with the minimun of ASE counts in order to run the ASE side of the model, defaults to NULL
#' @param min.ase.n minimun number individuals with the minimun of ASE counts, defaults to NULL
#' @param tag.threshold numeric with r2 threshold (0-1) for grouping snps to reduce the number of running tests, to disable use "no"
#' @param out path to save outputs, default to current directory
#' @param prefix optional prefix for saving tables, if NULL gene_id.eqtl will be used
#' @param model whether to run NB-ASE (full model negative binomial and allele specific counts),NB (negative binomial only) or both (NB-ASE and NB for those associations with no ASE information)
#' @param stan.model compiled stanmodel object with stan model, defaults NULL which uses  built-in NB-ASE model. When AI_estimate is provided the model corrects for reference panel bias, otherwise id doesn't.
#' @param stan.negonly compiled stanmodel object with neg only side, defaults to NULL which is built-in model.
#' @param prob  number p∈(0,1) indicating the desired probability mass to include in the intervals, defaults to 0.99, 0.95 quantiles
#' @param prior named list: mean= vector with the mean of Gaussians, sd= vector with Gaussians sd for eQTL effect prior, mix=vector with mixing proportions. Defaults to NULL, mixture of 3 components with mean (0,0,0) ;sd  c( 0.0436992,  0.3492696, 0.4920049); and mixing proportions  c(0.955, 2*0.015, 0.015).
#' @param ex.fsnp, if character: vector with pos:ref:alt for fsnps to exclude,  defaults to NULL
#' @param AI_estimate full name to data table with AI estimates for reference panel bias for fSNPs, defaults to NULL
#' @param pretotalReads numeric indicating a cut-off for total initial reads to consider AI estimates, defaults to 100
#' @export
#' @return Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of excluded rsnps.

baseqtl.gt.paired <- function(gene, chr, snps = 5 * 10^5, counts.f, covariates = 1, additional_cov = NULL,
                              e.snps, u.esnps = NULL, gene.coord, vcf, le.file, h.file,
                              population = c("EUR", "AFR", "AMR", "EAS", "SAS", "ALL"), nhets = 5,
                              min.ase = 5, min.ase.het = 5, tag.threshold = .9, out = ".", prefix = NULL,
                              model = c("both", "NB-ASE", "NB"), stan.model = NULL, stan.negonly = NULL,
                              prob = NULL, prior = NULL, ex.fsnp = NULL, AI_estimate = NULL,
                              pretotalReads = 100, mc.cores = getOption("mc.cores", 1L),
                              inference.method = c("sampling", "vb")) {

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

    print("Running NB model")

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

    if (!is.null(prefix)) {
      write.table(neg.sum, paste0(out, "/", prefix, ".paired.GT.stan.summary.txt"), row.names = FALSE)
    } else {
      write.table(neg.sum, paste0(out, "/", gene, ".paired.GT.stan.summary.txt"), row.names = FALSE)
    }
  }


  if (exists("full.sum") & !exists("neg.sum")) {
    if (!is.null(prefix)) {
      write.table(full.sum, paste0(out, "/", prefix, ".paired.GT.stan.summary.txt"), row.names = FALSE)
    } else {
      write.table(full.sum, paste0(out, "/", gene, ".paired.GT.stan.summary.txt"), row.names = FALSE)
    }
  }
  full.sum
}
