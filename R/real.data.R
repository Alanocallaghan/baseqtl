
#' Open vcf in R via bcftools query
#'
#' This function allows you to extract part of a vcf file and read it as data table
#' @param cl1 command line bcftools query to extract body of vcf
#' @param cl2 command line bcftools query to extract header of vcf
#' @param sep separation fields based on bcftools query, defaults =" "
#' @keywords vcf fread
#' @export
#' @return data table
#' vcf_cl()

vcf_cl <- function(cl1, cl2, sep = " ") {
  tmp <- data.table::fread(cmd = cl1, header = F, sep = sep, colClasses = c(V4 = "character", V5 = "character"))
  names(tmp) <- gsub("^.*]", "", names(data.table::fread(cmd = cl2)))[-1]
  return(tmp)
}


#' command line to extract any field via bcftools query, allow to select samples or region if required
#'
#' This function allows you to write the command to extract any field for a region of a vcf, allows sample selection
#' @param vcf full path to vcf or bcf, compressed or uncompressed
#' @param chr chromosome to extract
#' @param st start position to extract
#' @param end end position to extract
#' @param samples character vector with samples to extract, defaults=NULL
#' @param f.arg character vector with -f argument for bcftools query, defaults to GT and ASE
#' @param part whether to extract body or header of vcf
#' @keywords command line bcftools query GT ASE
#' @export
#' @return character vector with bcftools command
#' cl_bcfq()

cl_bcfq <- function(vcf, chr = NULL, st = NULL, end = NULL, samples = NULL, f.arg = NULL, part = c("body", "header")) {
  if (!is.null(samples)) { ## select samples first
    samp <- paste(" -s ", paste0(samples, collapse = ","))
  }
  if (!is.null(chr) & !is.null(st) & !is.null(end)) {
    reg <- paste0(" -r ", chr, ":", st, "-", end)
  }


  if (part == "body") {
    if (is.null(samples) & !exists("reg") & is.null(f.arg)) {
      x <- paste('bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', vcf)
    }
    if (!is.null(samples) & !exists("reg") & is.null(f.arg)) {
      x <- paste(' bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', samp, vcf)
    }
    if (is.null(samples) & exists("reg") & is.null(f.arg)) {
      x <- paste('bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', reg, vcf)
    }
    if (!is.null(samples) & exists("reg") & is.null(f.arg)) {
      x <- paste(' bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', reg, samp, vcf)
    }
    if (is.null(samples) & !exists("reg") & !is.null(f.arg)) {
      x <- paste("bcftools query -f ", f.arg, vcf)
    }

    if (!is.null(samples) & !exists("reg") & !is.null(f.arg)) {
      x <- paste(" bcftools query -f ", f.arg, samp, vcf)
    }

    if (!is.null(samples) & exists("reg") & !is.null(f.arg)) {
      x <- paste(" bcftools query -f ", f.arg, reg, samp, vcf)
    }

    if (is.null(samples) & exists("reg") & !is.null(f.arg)) {
      x <- paste(" bcftools query -f ", f.arg, reg, vcf)
    }
  } else {
    if (is.null(samples) & is.null(f.arg)) {
      x <- paste('bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', vcf, " -H  | head -1 ")
    }
    if (!is.null(samples) & is.null(f.arg)) {
      x <- paste0(' bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n"', samp, vcf, " -H | head -1")
    }
    if (is.null(samples) & !is.null(f.arg)) {
      x <- paste("bcftools query -f  ", f.arg, vcf, " -H  | head -1 ")
    }
    if (!is.null(samples) & !is.null(f.arg)) {
      x <- paste("bcftools query -f", f.arg, samp, vcf, " -H | head -1")
    }
  }

  return(x)
}




#' Extract GT and ASE for a gene
#'
#' This function allows you to extract GT and ASE for a gene, wraper for cl_bcfq and vcf_cl plus some formatting, removes samples with missing data or unphased. also removes snps if homo or missing for all samples. Adds id col to ease matching with legend/hap format
#' @param vcf full path to vcf
#' @param chr chromosome to extract, null for whole vcf
#' @param st start position to extract, null for whole vcf
#' @param end end position to extract, null for whole vcf
#' @param samples character vector with samples to extract, defaults=NULL
#' @param f.arg character vector with -f argument for bcftools query, defaults to GT and ASE
#' @param qc use function for qc purpose only
#' @param exclude whether to return a list with snps excluded from vcf (all homo or missing)
#' @export
#' @return data table with GT and ASE, unless !is.null(excluded), returns list with first element data table with GT and ASE and second element a data table with exluded snps.
#' vcf_w()

vcf_w <- function(vcf, chr = NULL, st = NULL, end = NULL, samples = NULL, f.arg = NULL, qc = NULL, exclude = NULL) {
  body <- cl_bcfq(vcf, chr, st, end, samples, f.arg, part = "body")

  header <- cl_bcfq(vcf, chr, st, end, samples, f.arg, part = "header")


  ## open GT and ASE for the selected gene

  gt.ase <- tryCatch(
    {
      vcf_cl(body, header, sep = " ")
    },
    error = function(e) {
      paste("Region not found for chrom and positions", chr, st, end, sep = ":")
    }
  )
  if (is.character(gt.ase)) {
    return(gt.ase)
  } else {
    if (!is.null(qc)) {
      return(gt.ase)
    } else {

      ## recode names gt.ase to make it compatible with Cincinatti files and functions
      names(gt.ase) <- gsub(":", "_", names(gt.ase))

      ## replace unphased data with "." missing value

      for (col in grep("_GT", names(gt.ase))) {
        data.table::set(gt.ase, i = grep("[0-1]/[0-1]", gt.ase[[col]]), j = col, value = ".")
      }

      ## exclude  snps if homo or missing for all samples

      ex <- apply(gt.ase[, grep("_GT", names(gt.ase)), with = F], 1, function(i) {
        u <- unique(i)
        setequal(u, c("0|0", ".")) | setequal(u, c("1|1", ".")) | setequal(u, "0|0") | setequal(u, ".") | setequal(u, "1|1")
      })


      ## this col will help to match snps with legend/hap reference panel
      gt.ase[, id := paste(POS, REF, ALT, sep = ":")]

      if (is.null(exclude)) {
        ## select relevant snps
        gt.ase <- gt.ase[which(ex == FALSE), ]
        return(gt.ase)
      } else {
        excl <- gt.ase[which(ex == TRUE), ]
        excl[, reason := "Missing or homo GT all samples"]
        excl <- excl[, .(id, reason)]
        l <- list(keep = gt.ase[which(ex == FALSE), ], excluded = excl)
        return(l)
      }
    }
  }
}


#' Check if gt.ase, output from vcf_w has phased GT in GT field, can also save a new vcf with excluding wrong GT by snp or by sample, as required.
#'
#' This function allows you to test wether a vcf_w returned object is correctly formatted in GT field with the option to list snps or samples with wrong GT in format to be excluded from
#' @param gt.ase object returned from vcf_w
#' @param exclude optional argument,removes entries with wrong GT format by snps or by sample, options "snps" or "samples"
#' @param vcf.path path and file name of original vcf, argument used for preparing new vcf with wrong GT entries removed
#' @param path optional, path to save new vcf excluding wrongly formatted GT, the default corresponds to the working directory
#' @param vcf.out optional, prefix for new vcf with wrong GT entries removed
#' @keywords vcf gt qc
#' @export
#' @return named vector when the only argument used is gt.ase. The vector gives the total number of snps, number of snps with wrong GT format, total number of samples and number of samples with wrong GT format. When all arguments are used it saves and indexes a new vcf excluding wrongly GT entries by snps or samples in format vcf.gz. In this mode the function returns a DT with the chr:pos:ref:alt of snps excluded or the name of the samples excluded.
#' vcf.gt.qc()

vcf.gt.qc <- function(gt.ase, exclude = c("snps", "samples"), vcf.path, path = ".", vcf.out = "chr22.GTqc") {
  gt.col <- grep("GT$", names(gt.ase))
  ok <- c("0|1", "0|0", "1|1", "1|0")
  by.snp <- apply(gt.ase[, gt.col, with = FALSE], 1, function(i) sum(!i %in% ok) != 0)
  by.sample <- apply(gt.ase[, gt.col, with = FALSE], 2, function(i) sum(!i %in% ok) != 0)
  report <- c(nrow(gt.ase), sum(by.snp), length(gt.col), sum(by.sample))
  names(report) <- c("total snps", "snps with wrong GT", "total samples", "samples with wrong GT")
  if (!(missing(exclude) & missing(vcf.path) & missing(vcf.out))) {
    na.ex <- pmatch(exclude, c("snps", "samples"))
    if (is.na(na.ex) | missing(vcf.path) | missing(vcf.out)) {
      stop("invalid 'exclude','vcf.path' or 'vcf.out' argument")
    } else {
      out <- paste0(path, "/", vcf.out, ".vcf.gz")
      if (exclude == "snps") {
        DT <- gt.ase[which(by.snp), .(CHROM, POS, REF, ALT)]
        DT[, ex := paste0(CHROM, ":", POS)]
        DT[, snps.excluded := paste0(ex, ":", REF, ":", ALT)]
        del <- paste0("^", paste0(DT$ex, collapse = ","))
        bcf.f <- paste("bcftools view -Oz -t", del, vcf.path, "-o", out)
        report <- data.table::data.table(snps.excluded = DT$snps.excluded)
      } else {
        del <- names(by.sample)[by.sample]
        del <- gsub(".GT$", "", del)
        del2 <- paste0("^", paste0(del, collapse = ","))
        bcf.f <- paste("bcftools view -Oz -s", del2, vcf.path, "-o", out)
        report <- data.table::data.table(samples.excluded = del)
      }
      system(bcf.f)
      bcf.i <- paste("bcftools index -t", out)
      system(bcf.i)
    }
  }

  return(report)
}



#' get start and end of cis-window per gene to make a bcftools query from a vcf file
#'
#' This function allows you to write the command to extract GT and ASE for a region of a vcf
#' @param file full path to file with gene coordinates, as prepared in inputs.R
#' @param chr chromosome to extract
#' @param gene gene id
#' @param cw length of cis-window, defaults to 5*10^5
#' @keywords command line bcftools query cis window
#' @export
#' @return vector with start and end of cis-window for the selected gene

cl_coord <- function(file, chr, gene, cw = 500000) {

  ## cat(x) check if command looks ok then run with system, copy cat(x) output to shell and check if it works

  g.st <- paste0("awk '$2 ==", chr, "' ", file, " | grep ", gene, " | cut -d ' ' -f4 | sed 's/,.*//' ")

  g.end <- paste0("awk '$2 ==", chr, "' ", file, " | grep ", gene, " | cut -d ' ' -f5 | sed 's/.*,//' ")



  window.st <- as.numeric(system(g.st, intern = TRUE)) - cw

  window.end <- as.numeric(system(g.end, intern = TRUE)) + cw

  v <- c(window.st, window.end)
  names(v) <- c("start", "end")
  return(v)
}


#' Extract haps from hap file for the whole cis-window (range)
#'
#' This function allows you extract haplotypes for a range of snps
#' @param file1 full path to file to legend.gz file
#' @param file2 full path to hap.gz file
#' @param cw vector start and end position of snps to extract within the range
#' @param population ethnicity to set a cut-off for maf: AFR AMR EAS EUR SAS ALL
#' @param maf cut-off for maf
#' @keywords haplotypes reference panel range snps
#' @export
#' @return matrix with haplotypes as in hap format (each col one hap), rownames are snp id

haps.range <- function(file1, file2, cw, population = "EUR", maf = 0.05) {
  eth <- 6:11 ## field number in legend file for ethniticy
  names(eth) <- c("AFR", "AMR", "EAS", "EUR", "SAS", "ALL")

  ## get snp info for range including line number, from legend file
  snp.i <- system(paste0("gunzip -c ", file1, " | awk '{if ($2 >= ", cw[1], "&& $2 <= ", cw[2], ") {print NR \" \" $2\":\"$3 \":\" $4 \" \" $", unname(eth[names(eth) == population]), "} }'"), intern = TRUE) ## second field in legend file is POS,then ref then alt allele

  if (length(snp.i) == 0) {
    return("no snps in reference panel")
  }

  ## format snp.i as DT

  DT <- data.table::as.data.table(lapply(1:3, function(i) unlist(lapply(strsplit(snp.i, split = " "), `[[`, i))))
  names(DT) <- c("line", "snp", "maf")
  DT[, line := as.numeric(line) - 1] ## first line in legend file is headings, need to substract 1 to match hap.gz file
  DT[, maf := as.numeric(maf)]
  haps <- paste0("gunzip -c ", file2, " | sed -n '", DT$line[1], ",", DT$line[nrow(DT)], "p' ")

  rf <- data.table::fread(cmd = haps, header = F) ## referene panel for snps in hap format
  ## remove snps below maf cut-off

  keep <- which(DT$maf >= maf & DT$maf < (1 - maf))
  rf <- rf[keep, ]
  mat <- as.matrix(rf)
  rownames(mat) <- DT$snp[keep]
  return(mat)
}




#' Extracts snp eaf from legend file for a set of snps
#'
#' This function allows you extract snp info from legend file for a range of snps
#' @param file1 full path to file to legend.gz file
#' @param snps character vector with id of snps to extract information, id=POS:REF:ALT
#' @param population ethnicity to set a cut-off for maf: AFR AMR EAS EUR SAS ALL
#' @keywords snp eaf reference panel
#' @export
#' @return DT with snp id and EAF from legend file for selected population

snp.eaf <- function(file1, snps, population = "EUR") {
  eth <- 6:11 ## field number in legend file for ethniticy
  names(eth) <- c("AFR", "AMR", "EAS", "EUR", "SAS", "ALL")
  ## get first and last POS within snps
  tmp <- as.numeric(gsub(":.*", "", snps))
  cw <- c(min(tmp), max(tmp))

  ## get snp info for range including line number, from legend file
  snp.i <- system(paste0("gunzip -c ", file1, " | awk '{if ($2 >= ", cw[1], "&& $2 <= ", cw[2], ") {print NR \" \" $2\":\"$3 \":\" $4 \" \" $", unname(eth[names(eth) == population]), "} }'"), intern = TRUE) ## second field in legend file is POS,then ref then alt allele

  if (length(snp.i) == 0) {
    return("no snps in reference panel")
  }

  ## format snp.i as DT

  DT <- data.table::as.data.table(lapply(1:3, function(i) unlist(lapply(strsplit(snp.i, split = " "), `[[`, i))))
  names(DT) <- c("line", "snp", "eaf")
  DT[, eaf := as.numeric(eaf)]
  DT[, line := NULL]
  DT <- DT[snp %in% snps, ]
  ## DT has unique entries, tmp may have duplicated snps
  DT2 <- merge(data.table::data.table(snp = snps), DT, by = "snp", sort = F)
  return(DT2)
}



#' Total gene and ASE counts, per fsnp, per individual
#'
#' Get total and AS counts per snp per individual
#' @param x DT with ASE and GT created from reading vcf
#' @param y data table with total counts for samples
#' @param z data table with each row the genotype for 1 rsnp coded as 0,1,-1 or 2, output from a rec_mytrecase_rSNPs
#' @keywords counts gene ASE
#' @export
#' @return matrix if z=NULL or list of  data tables, each data table corresponds to each rsnp, cols are total counts (y), GT 0,1,-1,2 for the rSNP, ase counts per fsnps across samples
#' tot.ase_counts

tot.ase_counts <- function(x, y = NULL, z = NULL) {
  as <- grep("_AS", names(x), value = T)
  tmp <- x[, as, with = F]
  ## missing values in phaser for AS are ".", happens when dealing with rna GT, convert to 0,0, as they wont affect counts
  tmp[tmp == "."] <- "0,0"

  ## get counts for hap2 (n)
  tmp2 <- sapply(1:ncol(tmp), function(i) as.numeric(unlist(lapply(strsplit(tmp[[i]], ","), `[[`, 2))))
  ## get counts for hap1+2 (m)
  tmp12 <- sapply(1:ncol(tmp), function(i) as.numeric(unlist(lapply(strsplit(tmp[[i]], ","), function(j) sum(as.numeric(j))))))
  if (!is.matrix(tmp2)) {
    tmp3 <- matrix(data = c(tmp2, tmp12), nrow = 2, byrow = T)
    rownames(tmp3) <- paste0(x$id, c(".n", ".m"))
    colnames(tmp3) <- gsub("_AS", "", as)
    tmp3 <- t(tmp3)
  } else {
    rownames(tmp2) <- paste0(x$id, ".n")
    rownames(tmp12) <- paste0(x$id, ".m")
    colnames(tmp2) <- colnames(tmp12) <- gsub("_AS", "", as)
    tmp3 <- t(rbind(tmp2, tmp12))
    tmp3 <- tmp3[, sort(colnames(tmp3), decreasing = T)]
  }
  if (is.null(z)) {
    return(tmp3)
  } else {

    ## add gene counts to ase
    y2 <- y[, which(names(y) %in% rownames(tmp3)), with = FALSE]
    tmp4 <- cbind(t(y2), tmp3)
    ## add GT of rsnps
    l <- lapply(1:nrow(z), function(i) {
      tmp5 <- cbind(t(z[i, grep("_GT", names(z), value = T), with = F]), tmp4)
      colnames(tmp5)[1:2] <- c("rsnp", "y")
      rownames(tmp5) <- gsub("_GT", "", rownames(tmp4))
      tmp5 <- data.table::data.table(tmp5, keep.rownames = T)
      return(tmp5)
    })
    return(l)
  }
}


#' Recode GT from (0,1,-1,2) scale  to GUESSFM scale 0 M, 1 hom ref, 2 het, 3 hom alt
#'
#' This function allows you to recode GT for input in tag function from GUESSFM
#' @param DT data table GT coded in trecase scale, rows SNPS, cols samples plus additionals
#' @keywords recode GUESSFM
#' @export
#' @return matrix with rows samples and cols SNPS

rec.guess <- function(DT) {
  M <- t(as.matrix(DT[, grep("_GT", names(DT), value = T), with = F]))
  ## recode
  M[M == 2] <- 3
  M[abs(M) == 1] <- 2
  M[M == 0] <- 1
  M[is.na(M)] <- 0
  M[M == "."] <- 0
  colnames(M) <- DT$id
  rownames(M) <- grep("_GT", names(DT), value = T)
  M <- apply(M, 2, as.numeric)
  return(M)
}


#' function for filtering input for stan, also option to remove fsnps with 0 m counts in all samples
#'
#' This function allows you to test whether a rsnp has enough het inds with a certain number of ASE counts
#' @param geno.exp data table for an element of output list tot.ase_counts
#' @param ase ase cut-off default is 5 counts, trecase default 5
#' @param n minimun number of individuals het for rsnp with counts >=ase, trease default 5
#' @param rem whether to exclude fsnps with 0 m counts in all samples, defaults to NULL
#' @keywords filtering stan input
#' @export
#' @return data table with sample name, GT rsnp, total counts, n and m counts for each fSNP, same data table outputted as an element of tot.ase_counts
#'
#' filt.rsnp()

filt.rsnp <- function(geno.exp, ase = 5, n = 5, rem = NULL) {
  n.col <- grep("\\.n", names(geno.exp), value = T)
  m.col <- grep("\\.m", names(geno.exp), value = T)
  m.counts <- rowSums(geno.exp[, m.col, with = F])
  # select ASE input when total ase counts are above threshold and
  A <- which(m.counts >= ase)
  if (nrow(geno.exp[A, ][abs(rsnp) == 1, ]) < n) {
    return("Not enough individuals with ASE counts")
  }
  if (is.null(rem)) {
    return(geno.exp)
  }
  ## remove
  tmp <- geno.exp[, m.col, with = FALSE]
  rem <- names(tmp)[colSums(geno.exp[, m.col, with = F]) == 0]
  s <- sapply(rem, function(i) sub("\\.m", "", i))
  r <- sapply(s, function(i) grep(i, names(geno.exp)))
  return(geno.exp[, r := NULL])
}


#' function for filtering fsnps, allows to select filtering if n=0 or n=m (likely GT error), also option to remove fsnps with less than cut-off total ASE counts in all samples, INPUT matrix
#'
#' This function allows you to test whether a fsnp has enough het inds with a certain number of ASE counts
#' @param c.ase matrix output of tot.ase_counts
#' @param ase total ase cut-off default is 5 counts, trecase default 5
#' @param min.ase.snp per snp ase cut-off
#' @param n minimun number of individuals het for fsnp with counts >=ase, trecase default 5
#' @param gt.err optional parameter, if "yes" it sets to m=0 when corresponding n=0 or n=m (likely GT error, homo called het), defaults to NULL
#' @param rem whether to exclude fsnps with less than min.ase.snps ASE counts in all samples, defaults to "yes", otherwise set it to null
#' @keywords filtering matrix m counts
#' @export
#' @return matrix with rownames sample name, cols: n and m counts for each fSNP, option to remove snps with ASE counts below cut-off in all individuals. Also ASE count below threshold (ase.min.snp) are replaced with 0, so will be excluded in subsequent analysis
#'
#' filt.fsnp()

filt.fsnp <- function(c.ase, ase = 5, min.ase.snp = 5, n = 5, gt.err = NULL, rem = "yes") {
  n.col <- grep("\\.n", colnames(c.ase), value = T)
  m.col <- grep("\\.m", colnames(c.ase), value = T)
  ## remove all entries below cut-off, both for n and m cols
  c.ase[c.ase < min.ase.snp] <- 0 ## entries with less than cut-off converted to 0
  ## remove counts when n=0 or n=m
  if (!is.null(gt.err)) {
    for (x in 1:length(n.col)) {
      w <- c.ase[, n.col[x]] == 0 | c.ase[, n.col[x]] == c.ase[, m.col[x]]
      c.ase[w, m.col[x]] <- 0
    }
  }


  m.counts <- rowSums(c.ase[, m.col, drop = F])
  ## select ASE input when total ase counts are above threshold, only inlcudes snps with counts > min.ase.snps
  A <- which(m.counts >= ase)
  if (nrow(c.ase[A, , drop = FALSE]) < n) {
    return("Not enough individuals with sufficient ASE counts per exonic snp")
  }
  ## remove
  if (is.null(rem)) {
    return(c.ase)
  }
  tmp <- c.ase[, m.col, drop = FALSE]
  keep <- colnames(tmp)[colSums(tmp) != 0] ##
  s <- sapply(keep, function(i) sub("\\.m", "", i))
  k <- sapply(s, function(i) grep(i, colnames(c.ase), value = T))
  return(c.ase[, k])
}

##############################################################
############ Functions for working with unknown rSNP GT ######
##############################################################

#' Get correlation of haps from reference panel, to input in tag from GUESSFM, from Chris
#'
#' calculates correlation by cols
#' @param x matrix (example haplotypes for the snps, cols snps, rows samples)
#' @keywords correlation haplotypes reference panel
#' @export
#' @return matrix of correlations
#' cor2
cor2 <- function(x) {
  SD.x <- apply(x, 2, sd)
  if (any(SD.x == 0)) stop("For some snps SD=0, remove and re-run")
  tmp <- 1 / (NROW(x) - 1) * crossprod(scale(x, TRUE, TRUE))
  return(tmp)
}

##' Derive tag SNPs for a SnpMatrix object using heirarchical clustering
##'
##' Uses complete linkage and the \code{\link{hclust}} function to define clusters,
##' then cuts the tree at 1-tag.threshold. Based on Chris Wallace's function tags but allowing for haplotype input.
##' @param X matrix of haplotypes, each rows is an observed hap, cols snps
##' @param tag.threshold threshold to cut tree, default=0.99
##' @param quiet if FALSE (default), show progress messages
##' @param method method used for heirarchical clustering.  See hclust for options.
##' @return character vector, names are \code{snps}, values are the tag for each SNP
##' @export
##' tag.noGT

tag.noGT <- function(X, tag.threshold = 0.99, quiet = FALSE, method = "single") {
  r2 <- (cor2(X))^2
  tmp <- 1 - r2
  D <- stats::as.dist(1 - r2)
  hc <- stats::hclust(D, method = method)
  clusters <- stats::cutree(hc, h = 1 - tag.threshold)

  snps.use <- names(clusters)[!duplicated(clusters)]
  groups <- split(names(clusters), clusters)

  ## now process each group, picking best tag
  n <- sapply(groups, length)
  names(groups)[n == 1] <- unlist(groups[n == 1])
  for (i in which(n > 1)) {
    g <- groups[[i]]
    a <- apply(r2[g, g], 1, mean)
    names(groups)[i] <- g[which.max(a)]
  }
  groups <- new("groups", groups, tags = names(groups))

  ## check
  r2 <- r2[GUESSFM::tags(groups), GUESSFM::tags(groups)]
  diag(r2) <- 0
  if (!quiet) {
    message("max r2 is now ", max(r2))
  }
  return(methods::as(groups, "tags"))
}

group.tags <- function(tags, keep) {
  groups <- tags[names(tags) %in% keep]
  groups <- split(names(groups), groups)
}


#####################################################################################

##' Calculate variance of expected genotype from a list with each element an individual and for each individual a vector of probabilities for GT=0,1 or 2.
##'
##' @title var.eg
##' @param l list, each element is for one individual, for each ind a vector with p(G)
##' @return character vector, names are \code{snps}, values are the tag for each SNP
##' @export
##' var.eg

var.eg <- function(l) {
  var_eg <- c()
  for (i in names(l)) {
    tmp <- sapply(l[[i]]$NB$p.g, function(k) sum(as.numeric(names(k)) * k))
    ## correct var, in r is divided n-1
    var_eg <- c(var_eg, var(tmp) * (length(tmp) - 1) / length(tmp))
  }
  names(var_eg) <- names(l)
  return(var_eg)
}



##' Calculates variance of genotype for a group of snps based on reference panel
##'
##' @param file1 path to legend file
##' @param file2 path to hap file
##' @param x vector with snps coded as pos:ref:alt
##' @return named vector with var(G)
##' @export
##' var.rp

var.rp <- function(file1, file2, x) {
  cw <- as.numeric(gsub(":.*", "", x))
  cw <- c(min(cw), max(cw))
  tmp <- haps.range(file1, file2, cw)
  tmp <- tmp[x, ]
  ## get GT from hap
  tmp2 <- sapply(seq(1, ncol(tmp), 2), function(i) rowSums(tmp[, i:(i + 1)]))
  ## get var per snp
  var.snp <- apply(tmp2, 1, var)
  return(var.snp)
}


#' wrap functions to calculate var(E(G)), var(G), r2=var(E(G))/var(G) and abs(`log2(aFC)_mean.ngt`-`log2(aFC)_mean.gt`) into a data.frame
#'
#' This function allows you to to calculate var(E(G)), var(G), r2=var(E(G))/var(G) and abs(`log2(aFC)_mean.ngt`-`log2(aFC)_mean.gt`) into a data.frame
#' @param genes, character vector of genes to select inputs from, defaults to NULL to use all genes run by model
#' @param path character vector with path to files with stan.input
#' @param pattern character vector with pattern to match
#' @param noGT data table with Gene_id and tags to extract gene/snps pairs
#' @param le.file path to legend file for specific chromosome
#' @param hap.file path to haplotype file for specific chromosome
#' @keywords stan inputvar(E(G))
#' @export
#' @return input data table noGT with new cols: var(E(G)), var(G) and r2
#' var.e()

var.e <- function(genes = NULL, path, pattern, noGT, le.file, hap.file) {
  ## stan input
  if (is.null(genes)) {
    tmp <- lapply(list.files(path, pattern = pattern, full.names = T), readRDS)
    genes <- names(tmp) <- gsub(pattern, "", gsub(".*ENS", "ENS", list.files(path, pattern = pattern, full.names = T)))
  } else {
    tmp <- lapply(genes, function(i) readRDS(list.files(path = path, pattern = paste0(i, pattern), full.names = T)))
    names(tmp) <- genes
  }
  ## Extract snps per gene from noGT
  snps <- lapply(genes, function(i) which(names(tmp[[i]]) %in% noGT[Gene_id == i, tag.ngt]))
  names(snps) <- genes

  ngt.tgs <- mapply(function(x, y) x[y], tmp, snps, SIMPLIFY = FALSE)

  ## Calculate Var(E(G)) for each gene-snp pair
  var.e1 <- lapply(ngt.tgs, var.eg)

  ## transform to DT
  var.e1 <- rbindlist(lapply(var.e1, function(i) data.table::data.table(snp = names(i), var.exp.g = i)), idcol = "Gene_id")

  ## Calculate var(G)
  var.g <- var.rp(file1 = le.file, file2 = hap.file, x = unique(var.e1$snp))

  ## make new cols in noGT

  DT <- merge(noGT, data.table::data.table(snp = names(var.g), var.rp = var.g), by.x = "tag.ngt", by.y = "snp", all.x = TRUE)

  DT <- merge(DT, var.e1, by.x = c("Gene_id", "tag.ngt"), by.y = c("Gene_id", "snp"), all.x = T)

  DT[, r2 := var.exp.g / var.rp][, abs.dif.log2.aFCg.ngt := abs(`log2(aFC)_mean.ngt` - `log2(aFC)_mean.gt`)]


  return(unique(DT))
}

#' wrap function to calculate r2=var(E(G))/var(G) to use for info input for baseqtl noGT
#'
#' This function allows you to to calculate select snps with r2=var(E(G))/var(G) above specific threshold
#' @param stan.noGT input list for one gene
#' @param rp.r matrix with reference panel haplotype info, cols individuals, rows snps
#' @param info, cut-off to remove snps, when r2<info, remove
#' @export
#' @return named vector with r2 and names snp_id, r2 above threshold
#' info.cut()

info.cut <- function(stan.noGT, rp.r, info) {

  ## Calculate Var(E(G)) for each snp


  var.e1 <- sapply(stan.noGT, function(i) {
    pg <- sapply(i$NB$p.g, function(k) sum(as.numeric(names(k)) * k))
    var <- mean(pg^2) - (mean(pg))^2
    return(var)
  })


  ## Calculate var(G)

  ## get GT from reference panel haps
  tmp2 <- sapply(seq(1, ncol(rp.r), 2), function(i) Matrix::rowSums(rp.r[, i:(i + 1), drop = FALSE]))
  if (!is.matrix(tmp2)) tmp2 <- matrix(tmp2, nrow = 1, dimnames = list(unique(names(tmp2)), NULL))
  ## get var per snp, var in R uses n-1
  var.g <- apply(tmp2, 1, function(i) var(i) * (ncol(tmp2) - 1) / ncol(tmp2))

  ## r2
  r2 <- var.e1 / var.g

  r2 <- r2[r2 >= info]

  return(r2)
}


#' Test of proportions for fsnps
#'
#' This function tests equality of the proportion of hets for fsnps in 2 populations, sample and reference panel
#' @param f.ase data table with fsnps and GT for samples (f.ase), cols for GT must end with "_GT"
#' @param rp.f with rownames fsnps and col haplotypes of fsnps
#' @param gene character with gene_id
#' @keywords test hets proportion
#' @export
#' @return matrix with rownames fsnps and colnames OR (odds ratio) and pvalue for Fisher's exact test

prop_het <- function(f.ase, rp.f, gene) {
  ## get sample GT, count hets and totals
  sam.GT <- f.ase[, grep("_GT$", names(f.ase), value = T)]
  sam.rec <- rec_mytrecase_rSNPs(f.ase$POS, f.ase)
  sam.GT <- sam.rec[, grep("_GT$", names(f.ase), value = T), with = F]
  sam.hets <- apply(sam.GT, 1, function(i) sum(abs(i) == 1, na.rm = T))
  sam.all <- apply(sam.GT, 1, function(i) sum(!is.na(i)))
  names(sam.hets) <- names(sam.all) <- f.ase$id
  sam.nohets <- sam.all - sam.hets

  ## get GT from reference panel haps, count hets and total=nrow(tmp)
  tmp <- sapply(seq(1, ncol(rp.f), 2), function(i) Matrix::rowSums(rp.f[, i:(i + 1), drop = F]))
  if (!is.matrix(tmp)) tmp <- matrix(tmp, nrow = 1, dimnames = list(unique(names(tmp)), NULL))

  ref.hets <- apply(tmp, 1, function(i) sum(abs(i) == 1))

  ## run fisher's exact test for each fsnp, takes hets and no-hets to work out odds ratio
  mat <- rbind(sam.hets[names(ref.hets)], sam.nohets[names(ref.hets)], ref.hets, ncol(tmp) - ref.hets)
  fish <- apply(mat, 2, function(i) {
    f <- fisher.test(x = matrix(i, nrow = 2), alternative = "two.sided")
    v <- c(f$estimate, f$p.value)
    names(v) <- c("OR", "pvalue")
    return(v)
  })

  fish <- data.table::data.table(t(fish), keep.rownames = T)
  fish[, gene_id := gene]
  data.table::setnames(fish, "rn", "fsnp")
  return(fish)
}
