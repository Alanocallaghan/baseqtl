#' Run baseqtl for 2 treatments  in one model from different smples
#' This function allows you to run baseqtl for one gene and multiple pre-selected snps.
#' @param treatments character vector with the 2 treatments (tissues, diseases) to study, sam order as counts.f, covariates and vcf
#' @param fishjoin whether to run Fisher test jointly in all samples or by treatment, defaults to jointly but  for QC purposes, to use the same fSNPs in ech treatment as the ones used  when the treatments were run using independent models select  NULL
#' @inheritParams baseqtl.nogt
#' @export
#' @return Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of rsnps that couldn't be run.

baseqtl2T.nogt <- function(gene, chr, snps = 5 * 10^5, counts.f, covariates = 1, additional_cov = NULL,
                           e.snps, u.esnps = NULL, gene.coord, vcf, sample.file = NULL, le.file,
                           h.file, population = c("EUR", "AFR", "AMR", "EAS", "SAS", "ALL"),
                           maf = 0.05, min.ase = 5, min.ase.snp = 5, min.ase.n = 5, tag.threshold = .9,
                           info = 0.3, out = ".", prefix = NULL, stan.model = NULL, prob = NULL,
                           prior = NULL, ex.fsnp = 0.01, AI_estimate = NULL, pretotalReads = 100,
                           treatments, fishjoin = "joint",
                           inference.method = c("sampling", "vb"),
                           mc.cores = getOption("mc.cores", parallel::detectCores())) {
  if (is.null(stan.model)) {
    ## check if ref panelbias correction
    if (is.null(AI_estimate)) {
      stan.model <- stanmodels$noGT2T_nb_ase
    } else {
      stan.model <- stanmodels$noGT2T_nb_ase_refbias
    }
  } else { ## model provided by user
    if (!inherits(stan.model, "stanmodel")) {
      stop("Model provided by user is not a stanmodel object")
    }
  }

  skin <- treatments ## compatibility with previous definitions

  inputs <- baseqtl2Tnogt.in(
    gene,
    chr,
    snps,
    counts.f,
    covariates,
    additional_cov,
    e.snps,
    u.esnps,
    gene.coord,
    vcf,
    sample.file,
    le.file,
    h.file,
    population,
    maf,
    min.ase,
    min.ase.snp,
    min.ase.n,
    tag.threshold,
    info,
    out,
    prefix,
    prior,
    ex.fsnp,
    prob,
    AI_estimate,
    pretotalReads,
    skin,
    fishjoin
  )

  if (is.character(inputs)) stop(inputs)
  ## rename inputs
  probs <- Reduce(intersect, inputs$probs)
  population <- inputs$population
  r.tag <- inputs$r.tag
  het.f <- inputs$het.f
  nfsnps <- inputs$inp$nfsnps
  info.ok <- inputs$inp$info.ok
  stan2T <- inputs$inp$stan.in

  if (!is.null(AI_estimate)) min_AI <- inputs$inp$min_AI


  cat("Running stan for ", length(info.ok), " rSNPS")
  stan.full <- parallel::mclapply(
    1:length(stan2T),
    function(i) {
      s <- run.stan(stan.model, data = stan2T[[i]], pars = c("ba", "bd", "bp", "bn"), probs = probs, inference.method = inference.method)
      dt.wide <- stan.paired.for(s, names(stan2T)[i])
      return(dt.wide)
    },
    mc.cores = mc.cores
  )

  ## get maf for tags run in model
  maf.t <- snp.eaf(le.file, names(stan2T), population)

  full.sum <- stan.2T(
    x = stan.full, rtag = r.tag, gene = gene, EAF = maf.t, info = info.ok,
    nfsnps = nfsnps,
    min.pval = het.f,
    probs = probs
  )

  ## relabel "bp", "bn" to "bt1", "bt2"
  m <- mapply(relab,
    a = c("bp", "bn"),
    b = c("bt1", "bt2"),
    MoreArgs = list(full.sum),
    SIMPLIFY = F
  )

  full.sum[, model := inputs$model]
  if (!is.null(AI_estimate)) full.sum[, min_AI := min_AI]

  if (!is.null(prefix)) {
    write.table(dt, paste0(out, "/", prefix, ".2T.noGT.summary.txt"), row.names = FALSE)
  } else {
    write.table(full.sum, paste0(out, "/", gene, ".2T.noGT.summary.txt"), row.names = FALSE)
  }
  full.sum
}
