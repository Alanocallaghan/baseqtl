rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#' Run baseqtl for 2 treatments  in one model from different smples
#' This function allows you to run baseqtl for one gene and multiple pre-selected snps.
#' @param gene gene id for the gene to run
#' @param chr chromosome where the gene is, example chr=22
#' @param snps either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window
#' @param counts.f path to files with filtered counts: rows genes, first col gene_id followed by samples, prepared in inputs.R
#' @param covariates path to matrix of covariates prepared in inputs.R, if no covariates, covariates =1, default. Same  order of treatments as counts.f
#' @param additional_cov full name to file with first column sample names and additional columns gene independent covariates, defaults to NULL
#' @param e.snps path to file listing exonic snps for the chromosome where the gene is, prepared in input.R
#' @param u.esnps whether to use unique exonic snps per gene, defaults to NULL when it is not necessary if strand info is known
#' @param gene.coord path to file listing gene coordinates and exons, prepared in input.R
#' @param vcf path to vcf files with ASE and GT for the exonic snps for the chromosome where the gene is, same order of treatment as counts.f and covariates
#'  @param sample.file sample file for the reference panel (sample description), to be used if ex.fsnp test is required and population is not the whole reference panel
#' @param le.file path to gz legend file (legend.gz) for the chromosome of interest for the reference panel (snp description)
#' @param h.file path to gz haplotpe file for the chromosome of interest for the reference panel (haplotypes for all samples in reference panel)
#' @param population ethnicity to set a cut-off for maf: AFR AMR EAS EUR SAS ALL, defaults to EUR
#' @param maf cut-off for maf, defaults to 0.05
#' @param min.ase minimun number of ASE counts for an individual in order to be included, defaults to 5
#' @param min.ase.snp minum number of ASE counts for a single snp to be considered, for a particular individual
#' @param min.ase.n minimun number individuals with the minimun of ASE counts, defaults to 5.
#' @param tag.threshold numeric with r2 threshold (0-1) for grouping snps to reduce the number of running tests, to disable use "no"
#' @param info numeric cut-off for var(E(G))/var(G), var(E(G)) is the expected variance for input and var(G) for reference panel, similar to info score in impute2, defaults to 0.3. rsnps with lower info wont be run by stan.
#' @param out path to save outputs, default to current directory
#' @param prefix optional prefix for saving files, if NULL gene_id will be used
#' @param model compiled stanmodel object with stan model, defaults NULL which uses  built-in NB-ASE model. When AI_estimate is provided the model corrects for reference panel bias, otherwise id doesn't.
#' @param prob  number p∈(0,1) indicating the desired probability mass to include in the intervals, defaults to 0.95 and 0.99 quantiles
#' @param prior named list: mean= vector with the mean of Gaussians, sd= vector with Gaussians sd for eQTL effect prior, mix=vector with mixing proportions. Defaults to NULL, mixture of 3 components with mean (0,0,0) ;sd  c( 0.0436992,  0.3492696, 0.4920049); and mixing proportions  c(0.955, 2*0.015, 0.015).
#' @param ex.fsnp if character: vector with pos:ref:alt for fsnps to exclude, if numeric  p-value cut-off for fSNPs to exclude based on fisher test for suspected genotype error based on comparing the proportion of hets in the sample and reference panel,  defaults to 0.01. If NULL no fSNP will be excluded.
#' @param AI_estimate full name to data table with AI estimates for reference panel bias for fSNPs, defaults to NULL
#' @param pretotalReads numeric indicating a cut-off for total initial reads to consider AI estimates, defaults to 100
#' @param treatments character vector with the 2 treatments (tissues, diseases) to study, sam order as counts.f, covariates and vcf
#' @param fishjoin whether to run Fisher test jointly in all samples or by treatment, defaults to jointly but  for QC purposes, to use the same fSNPs in ech treatment as the ones used  when the treatments were run using independent models select  NULL
#' @export
#' @return Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of rsnps that couldn't be run.

baseqtl2T.nogt <- function(gene, chr, snps = 5 * 10^5, counts.f, covariates = 1, additional_cov = NULL, e.snps, u.esnps = NULL, gene.coord, vcf, sample.file = NULL, le.file, h.file, population = c("EUR", "AFR", "AMR", "EAS", "SAS", "ALL"), maf = 0.05, min.ase = 5, min.ase.snp = 5, min.ase.n = 5, tag.threshold = .9, info = 0.3, out = ".", prefix = NULL, model = NULL, prob = NULL, prior = NULL, ex.fsnp = 0.01, AI_estimate = NULL, pretotalReads = 100, treatments, fishjoin = "joint", inference.method="sampling") {
  if (is.null(model)) {
    ## check if ref panelbias correction
    if (is.null(AI_estimate)) {
      model <- stanmodels$noGT2T_nb_ase
    } else {
      model <- stanmodels$noGT2T_nb_ase_refbias
    }
  } else { ## model provided by user
    if (class(model) != "stanmodel") {
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
      s <- run.stan(model, data = stan2T[[i]], pars = c("ba", "bd", "bp", "bn"), probs = probs, method = inference.method)
      dt.wide <- stan.paired.for(s, names(stan2T)[i])
      return(dt.wide)
    }
  )

  ## get maf for tags run in model
  maf.t <- snp.eaf(le.file, names(stan2T), population)

  full.sum <- stan.2T(
    x = stan.full, rtag = r.tag, gene = gene, EAF = maf.t, info = info.ok,
    nfsnps = nfsnps,
    min.pval = het.f,
    probs = probs
  )

  full.sum[, model := inputs$model]
  if (!is.null(AI_estimate)) full.sum[, min_AI := min_AI]

  if (!is.null(prefix)) {
    write.table(dt, paste0(out, "/", prefix, ".2T.noGT.summary.txt"), row.names = FALSE)
  } else {
    write.table(full.sum, paste0(out, "/", gene, ".2T.noGT.summary.txt"), row.names = FALSE)
  }
}
