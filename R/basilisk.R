#' @importFrom basilisk BasiliskEnvironment basiliskStart basiliskStop basiliskRun
bcf_env <- BasiliskEnvironment("bcf_env", pkgname="baseqtl",
    packages=c("bcftools==1.15.1", "htslib==1.15.1"),
    channels = "bioconda")
