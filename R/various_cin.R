
#' Format txt converted vcf files: one file per chr
#'
#' This function allows you to get good format from vcf text.
#' @param file file name with full path to vcf
#' @param head path and file with header, defaults to prefix.header.txt
#' @keywords vcf format
#' @export
#' @return named list with info from vcf files
#' name()
name <- function(file, head=NULL){
    
    temp <- data.table::fread(file, sep='\t')
    if(is.null(head)){
        temp2 <- data.table::fread(paste0(gsub(".txt","",file),".header.txt"))
    } else {
        temp2 <- data.table::fread(head, sep="\t")
    }
    
    names(temp)<- gsub("^.*]","",names(temp2))
    names(temp)<-gsub(":","_", names(temp))

    return(temp)
   
        
}


#' Get haps of fsnp plus rSNP, improved for more consistent output
#' 
#' extract haps from vcf file for fsnps and 1 rsnp
#' @param x data table of fsnps, rows fsnps and columns phased GT, plus any additional
#' @param w sample names, default to NULL as can be deduced from x
#' @param y rsnp id, id="POS:REF:ALT"
#' @param z data table or rsnps, rows fsnps and columns phased GT, plus any additional
#' @keywords sample haps 
#' @export
#' @return list with 2 matrices, first for hap1: row samples and each col GT for fsnps and last rsnp
#' hap_sam2

hap_sam2 <- function(x,w=NULL,y=NULL,z=NULL){
    
    gt <- grep("_GT",names(x), value=T)
    tmp <- x[,gt,with=F]
    ##get rsnp
    if(!is.null(y) & !is.null(z)){
        z[,id:=paste(POS,REF,ALT, sep=":")]
        rsnp <- z[id==y,gt, with=F]
        ## add rsnp to fsnps
        tmp <- rbind(tmp,rsnp)
    }
    ## get hap1 
    tmp1 <- sapply(1:ncol(tmp), function(i) as.numeric(unlist(lapply(strsplit(tmp[[i]], "|"), `[[`, 1))))
    ## get hap2 
    tmp2 <- sapply(1:ncol(tmp), function(i) as.numeric(unlist(lapply(strsplit(tmp[[i]], "|"), `[[`, 3))))

    if(ncol(tmp) >1 & !is.matrix(tmp1)){
        l <- lapply(list(tmp1,tmp2), function(i) matrix(i,ncol=1))
        names(l) <- c("hap1","hap2")
    } else {
        ## transpose 
        l <- list(hap1=t(tmp1),hap2=t(tmp2))
    }
    if(is.null(w)){
        w <- gsub("_GT","",gt)
    }
    ## name samples
    l <- lapply(l,function(i) {row.names(i)=w; return(i)})
    
    return(l)
}       

#' Selects SNPs and recode GT to my implementation of trecase scale: 0 AA, 1 AB, 2 BB BA -1 
#' 
#' recodes fixed genotypes to my implementation of trecase for selected SNPs, samples with "_GT" suffix.
#' @param x vector with SNP positions for a chr
#' @param y data.table with phased GP (0|1) for each sample, for the chr selected in x input
#' @param z vector with order fo samples, defaults to NULL
#' @param recode whether not to recode, defaults to yes
#' @keywords recode GT rSNP m=trecase
#' @export
#' @return data table with recoded GT per sample
#' rec_mytrecase_rSNPs

rec_mytrecase_rSNPs <- function(x,y, z=NULL, recode="yes"){
    tmp <- y[POS %in% x,]
    if(recode=="yes") {
        tmp2 <- tmp[,grep("_GT", names(tmp), value=T), with=F]
                                       
        for(i in seq_along(names(tmp2))) { #recode
                                        
            tmp2[get(names(tmp2)[i])=="0|0",names(tmp2)[i]:="0"]
            tmp2[get(names(tmp2)[i])=="0|1",names(tmp2)[i]:="1"]
            tmp2[get(names(tmp2)[i])=="1|0",names(tmp2)[i]:="-1"]
            tmp2[get(names(tmp2)[i])=="1|1",names(tmp2)[i]:="2"]
            tmp2[get(names(tmp2)[i])==".",names(tmp2)[i]:=NA]
            tmp2[get(names(tmp2)[i])=="./.",names(tmp2)[i]:=NA]
            
            tmp2[,names(tmp2)[i]:=as.numeric(get(names(tmp2)[i]))]
        }
        tmp[, grep("_GT", names(tmp), value=T) := tmp2]
    }
    if(!is.null(z)) {
        data.table::setcolorder(tmp,c(grep("_GT", names(tmp), invert=T, value=T), paste0(z,"_GT")))
    }
    return(tmp)
}       

