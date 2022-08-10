#' @importFrom basilisk BasiliskEnvironment
bcf_env <- BasiliskEnvironment("bcf_env", pkgname="baseqtl",
    packages=c("bcftools==1.14.18", htslib="1.14.8"),
    channels = "bioconda")
