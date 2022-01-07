#' Run baseqtl with unknown rsnp GT and missing values for GT fsnps with optional  reference panel bias correction
#'
#' This function allows you to run baseqtl for one gene and multiple pre-selected snps.
#' @param gene gene id for the gene to run
#' @param chr chromosome where the gene is, example chr=22
#' @param snps either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window
#' @param counts.f path to file with filtered counts: rows genes, first col gene_id followed by samples, prepared in inputs.R
#' @param covariates path to matrix of covariates prepared in inputs.R. If using gc correction (each gene different value), the matrix has rownames= genes and cols=samples plus extra columns if other covariates are added. If only using lib size or gene independent covariates, rows are samples and columns are covariates. If no covariates, covariates =1, default
#' @param additional_cov full name to file with first column sample names and additional columns gene independent covariates, defaults to NULL
#' @param e.snps path to file listing exonic snps for the chromosome where the gene is, prepared in input.R
#' @param u.esnps whether to use unique exonic snps per gene, defaults to NULL when it is not necessary if strand info is known
#' @param gene.coord path to file listing gene coordinates and exons, prepared in input.R
#' @param vcf path to vcf file with ASE and GT for the exonic snps for the chromosome where the gene is
#' @param sample.file sample file for the reference panel (sample description), to be used if ex.fsnp test is required and population is not the whole reference panel
#' @param le.file path to gz legend file (legend.gz) for the chromosome of interest for the reference panel (snp description)
#' @param h.file path to gz haplotpe file for the chromosome of interest for the reference panel (haplotypes for all samples in reference panel)
#' @param population ethnicity to set a cut-off for maf: AFR AMR EAS EUR SAS ALL, and for testing fSNPs to exclude (ex.fsnp below) if applicable, defaults to EUR
#' @param maf cut-off for maf, defaults to 0.05
#' @param min.ase minimun number of ASE counts for an individual in order to be included, defaults to 5
#' @param min.ase.snp minum number of ASE counts for a single snp to be considered, for a particular individual
#' @param min.ase.n minimun number individuals with the minimun of ASE counts, defaults to 5.
#' @param tag.threshold numeric with r2 threshold (0-1) for grouping snps to reduce the number of running tests, to disable use "no"
#' @param info numeric cut-off for var(E(G))/var(G), var(E(G)) is the expected variance for input and var(G) for reference panel, similar to info score in impute2, defaults to 0.3. rsnps with lower info wont be run by stan.
#' @param out path to save outputs, default to current directory
#' @param prefix optional prefix for saving files, if NULL gene_id.eqtl will be used
#' @param model compiled stanmodel object with stan model, defaults NULL which uses  built-in NB-ASE model. When AI_estimate is provided the model corrects for reference panel bias, otherwise id doesn't.
#' @param model.negonly compiled stanmodel object with neg only side, defaults to NULL which is built-in model.
#' @param prob  number p∈(0,1) indicating the desired probability mass to include in the intervals, defaults to 0.99 and 0.95 quantiles
#' @param prior named list: mean= vector with the mean of Gaussians, sd= vector with Gaussians sd for eQTL effect prior, mix=vector with mixing proportions. Defaults to NULL, mixture of 2 components with mean (0,0) ;sd  c( 0.0309, 0.3479); and mixing proportions  c(0.97359164, 0.02640836).
#' @param ex.fsnp if character: vector with pos:ref:alt for fsnps to exclude, if numeric  p-value cut-off for fSNPs to exclude based on fisher test for suspected genotype error based on comparing the proportion of hets in the sample and reference panel,  defaults to 0.01. If NULL no fSNP will be excluded.
#' @param AI_estimate full name to data table with AI estimates for reference panel bias for fSNPs, if NULL no reference panel bias correction is applied
#' @param pretotalReads numeric indicating a cut-off for total initial reads to consider AI estimates, defaults to 100
#' @param save_input whether to save input to stan model for QC purposes, defaults to FALSE to save disk space. Object ending with "noGT.stan.input.rds" is a named list with each element the inputs for a cis-SNP. For each cis-SNP there is a list of 2 elements: "NB" and "ase". "NB" is a list with elements "counts" and "p.g". "Counts" is a data.table with columns sample names and one row corresponding to the gene, values total read counts. "p.g" is a named list with each element a sample. For each sample there is an array with names genotypes (0,1,2) and values the genotype probabilities. For the "ase" list they are for elements: "m" numeric vector with  total ASE counts per sample. "g" list with each element a sample and for each sample the genoptype of the cis SNP coded as 0,1,2 and -1, with -1 indicating that the alternative allele is in haplotype 1. "p" has the same structure as "g" and indicates the probability for each genotype. "n"  is similar to "g" and "p" but contains the mapped reads to haplotype 2. The file ending with "noGT.fsnps.counts.rds is a matrix with rows samples and columns fSNPS. When a fSNPs ends with ".n" correspond to the counts matching the alternative allele and ".m" indicates the total counts matching the SNP.
#' @export
#' @return data.table with summary of gene-snp associations. Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of excluded rsnps.

baseqtl.nogt <- function(gene, chr, snps = 5 * 10^5, counts.f, covariates = 1, additional_cov = NULL,
                         e.snps, u.esnps = NULL, gene.coord, vcf, sample.file = NULL, le.file,
                         h.file, population = c("EUR", "AFR", "AMR", "EAS", "SAS", "ALL"),
                         maf = 0.05, min.ase = 5, min.ase.snp = 5, min.ase.n = 5, tag.threshold = .9,
                         info = 0.3, out = ".", prefix = NULL, model = NULL, model.negonly = NULL,
                         prob = NULL, prior = NULL, ex.fsnp = 0.01, AI_estimate = NULL,
                         pretotalReads = 100, save_input = FALSE,
                         mc.cores = getOption("mc.cores", parallel::detectCores()), inference.method = "sampling") {
  ## check stan models
  if (is.null(model)) {
    ## check if ref panelbias correction
    if (is.null(AI_estimate)) {
      model <- stanmodels$noGT_nb_ase
    } else {
      model <- stanmodels$noGT_nb_ase_refbias
    }
  } else { ## model provided by user
    if (!inherits(model, "stanmodel")) {
      stop("Model provided by user for argument model is not a stanmodel object")
    }
  }

  if (is.null(model.negonly)) {
    model.negonly <- stanmodels$noGT_nb
  } else {
    if (!inherits(model.negonly, "stanmodel")) {
      stop("Model provided by user for argument model.negonly is not a stanmodel object")
    }
  }

  if (!is.null(prior)) {
    if (!is.list(prior)) stop("prior argument must be a list")
    if (any(!names(prior) %in% c("mean", "sd", "mix"))) stop("prior argument must be a named 'mean' and 'sd'")
    if (length(unique(sapply(prior, length))) != 1) stop("mean and sd for prior argument must have the same length")
  } else {
    ## use default prior
    prior <- c(0, 0, 0.0309, 0.3479, 0.97359164, 0.02640836)
    k <- length(prior) / 3 ## number of gaussians
    s <- seq(1, length(prior), k)
    l <- lapply(1:3, function(i) as.numeric(prior[s[i]:(s[i] + k - 1)]))
    names(l) <- c("mean", "sd", "mix")
    prior <- l
  }


  inputs <- baseqtl.nogt.in(
    gene = gene,
    chr = chr,
    snps = snps,
    counts.f = counts.f,
    covariates = covariates,
    additional_cov = additional_cov,
    e.snps = e.snps,
    u.esnps = u.esnps,
    gene.coord = gene.coord,
    vcf = vcf,
    sample.file = sample.file,
    le.file = le.file,
    h.file = h.file,
    population = population,
    maf = maf,
    min.ase = min.ase,
    min.ase.snp = min.ase.snp,
    min.ase.n = min.ase.n,
    tag.threshold = tag.threshold,
    info = info,
    out = out,
    model = NULL, ## compatibility with GT functions
    prefix = prefix,
    prob = prob,
    ex.fsnp = ex.fsnp,
    AI_estimate = AI_estimate,
    pretotalReads = pretotalReads,
    save_input
  )



  if (is.character(inputs$inp)) stop(inputs$inp)

  stan.noGT2 <- inputs$inp$stan.noGT2
  if (!is.null(ex.fsnp) & is.numeric(ex.fsnp)) {
    het.f <- inputs$inp$het.f
    het.fall <- inputs$inp$het.all
  } else {
    het.f <- NULL
  }

  nfsnps <- inputs$inp$nfsnps
  info.ok <- inputs$inp$info.ok
  probs <- inputs$probs
  r.tag <- inputs$r.tag


  stan.noGT2 <- lapply(stan.noGT2, function(l) {
    ## add prior
    l <- add.prior(prior, l)
    return(l)
  })

  if (inputs$model == "NB-ASE") {
    model2run <- model
  } else {
    model2run <- model.negonly
  }

  model <- inputs$model # "NB-ASE" or 'NB'

  message("Running stan for ", model)

  stan.full <- parallel::mclapply(
    1:length(stan.noGT2),
    function(i) {
      s <- run.stan(model2run, data = stan.noGT2[[i]], pars = "bj", probs = probs, inference.method = inference.method)
      return(s)
    },
    mc.cores = mc.cores
  )
  names(stan.full) <- names(stan.noGT2)

  ## get maf for tags run in model
  maf.t <- snp.eaf(le.file, names(stan.noGT2), population)
  if (model == "NB") {
    nfsnps <- NA
  }
  ## stan.bt creates PEP from post.prop.neg
  full.sum <- stan.bt(x = stan.full, y = "bj", rtag = r.tag, model = model, gene = gene, EAF = maf.t, info = info.ok, nfsnps = nfsnps, min.pval = het.f, probs = probs)
  ## add min AI_post
  if (!is.null(AI_estimate) & model == "NB-ASE") {
    full.sum[, min_AI := inputs$inp$min_AI]
  }

  if (!is.null(prefix)) {
    write.table(full.sum, paste0(out, "/", prefix, ".noGT.stan.summary.txt"), row.names = FALSE)
    if (exists("het.fall")) {
      write.table(het.fall, paste0(out, "/", prefix, ".fsnps.het.fisher.test.txt"), row.names = FALSE)
    }
  } else {
    write.table(full.sum, paste0(out, "/", gene, ".noGT.stan.summary.txt"), row.names = FALSE)
    if (exists("het.fall")) {
      write.table(het.fall, paste0(out, "/", gene, ".fsnps.het.fisher.test.txt"), row.names = FALSE)
    }
  }

  return(full.sum)
}
