Compute per-barcode expression functional score
================
Tyler Starr
6/9/2022

This notebook reads in per-barcode counts from `count_variants.ipynb`
for expression Sort-seq experiments, computes functional scores for RBD
expression levels, and does some basic QC on variant expression
functional scores.

``` r
#list of packages to install/load
packages = c("yaml","data.table","tidyverse","fitdistrplus","gridExtra")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages],
                   lib=c(paste("/uufs/chpc.utah.edu/common/home/",Sys.getenv("USER"),"/RLibs/",Sys.getenv("R_VERSION"),sep="")),
                   repos=c("http://cran.us.r-project.org"))
}
#load packages
invisible(lapply(packages, library, character.only=T))

knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#read in config file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$expression_sortseq_dir)){
  dir.create(file.path(config$expression_sortseq_dir))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Rocky Linux 8.10 (Green Obsidian)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /uufs/chpc.utah.edu/sys/spack/v019/linux-rocky8-x86_64/gcc-13.1.0/openblas-0.3.21-wswkguokqgrxenycjlx4cxqijupp47mw/lib/libopenblasp-r0.3.21.so;  LAPACK version 3.9.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: America/Denver
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3      fitdistrplus_1.2-4 survival_3.5-8     MASS_7.3-60.2     
    ##  [5] lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4       
    ##  [9] purrr_1.0.2        readr_2.1.5        tidyr_1.3.1        tibble_3.2.1      
    ## [13] ggplot2_4.0.0      tidyverse_2.0.0    data.table_1.15.4  yaml_2.3.8        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Matrix_1.7-0       gtable_0.3.6       compiler_4.4.0     tidyselect_1.2.1  
    ##  [5] splines_4.4.0      scales_1.4.0       fastmap_1.1.1      lattice_0.22-6    
    ##  [9] R6_2.5.1           generics_0.1.3     knitr_1.46         pillar_1.9.0      
    ## [13] RColorBrewer_1.1-3 tzdb_0.4.0         rlang_1.1.3        utf8_1.2.4        
    ## [17] stringi_1.8.3      xfun_0.43          S7_0.2.0           timechange_0.3.0  
    ## [21] cli_3.6.2          withr_3.0.0        magrittr_2.0.3     digest_0.6.35     
    ## [25] grid_4.4.0         hms_1.1.3          lifecycle_1.0.4    vctrs_0.6.5       
    ## [29] evaluate_0.23      glue_1.7.0         farver_2.1.1       fansi_1.0.6       
    ## [33] rmarkdown_2.26     tools_4.4.0        pkgconfig_2.0.3    htmltools_0.5.8.1

## Setup

First, we will read in metadata on our sort samples, the table giving
number of reads of each barcode in each of the sort bins, and the
barcode-variant lookup tables, and merge these tables together.

``` r
#read dataframe with list of barcode runs
barcode_runs <- read.csv(file=config$barcode_runs,stringsAsFactors=F); barcode_runs <- subset(barcode_runs, select=-c(R1))

#eliminate rows from barcode_runs that are not from an expression sort-seq experiment
barcode_runs <- barcode_runs[barcode_runs$sample_type == "SortSeq",]

#read file giving count of each barcode in each sort partition
counts <- data.table(fread(config$variant_counts_file,stringsAsFactors=F))

#eliminate rows from counts that are not part of an expression sort-seq bin
counts <- subset(counts, sample %in% barcode_runs[barcode_runs$sample_type=="SortSeq","sample"])

#read in barcode-variant lookup tables
dt_BA32 <- data.table(read.csv(file=config$codon_variant_table_file_BA32,stringsAsFactors=F))
#dt_BA32 <- data.table(read.csv(file=config$codon_variant_table_file_BA32,stringsAsFactors=F))

#merge, eliminate barcodes duplicated within a library
dt <- rbind(dt_BA32); setkey(dt,barcode,library)
duplicates <- dt[duplicated(dt,by=c("barcode","library")),.(library,barcode)] #the data.table duplciates function annoyingly only flags the first of each duplicate so doesn't intrinsically allow removal of both of the entries of the duplicate. So, flag what are duplciates, and then remove
dt[,duplicate:=FALSE]
for(i in seq_len(nrow(duplicates))){
  dt[library==duplicates[i,library] & barcode==duplicates[i,barcode],duplicate:=TRUE]
}
#modified above from for (i in 1:nrow(duplicates)) bc I was getting a 0 error, I think bc I have just 1 lib so no need to merge and remove duplicates?
dt <- dt[duplicate==FALSE,]; dt[,duplicate:=NULL]

dt <- merge(counts, dt, by=c("library","barcode")); rm(dt_BA32);rm(dt_BA32);rm(counts); rm(duplicates)
```

    ## Warning in rm(dt_BA32): object 'dt_BA32' not found

Convert from Illumina read counts to estimates of the number of cells
that were sorted into a bin, and add some other useful information to
our data frame.

``` r
#for each bin, normalize the read counts to the observed ratio of cell recovery among bins
for(i in 1:nrow(barcode_runs)){
  lib <- as.character(barcode_runs$library[i])
  bin <- as.character(barcode_runs$sample[i])
  ratio <- sum(dt[library==lib & sample==bin,"count"])/barcode_runs$number_cells[i]
  if(ratio<1){ #if there are fewer reads from a sortseq bin than cells sorted
    dt[library==lib & sample==bin, count.norm := as.numeric(count)] #don't normalize cell counts, make count.norm the same as count
    print(paste("reads < cells for",lib,bin,", un-normalized (ratio",ratio,")")) #print to console to inform of undersampled bins
  }else{
    dt[library==lib & sample==bin, count.norm := as.numeric(count/ratio)] #normalize read counts by the average read:cell ratio, report in new "count.norm" column
    print(paste("read:cell ratio for",lib,bin,"is",ratio))
  }
}
```

    ## [1] "read:cell ratio for pool1 SortSeq_bin1 is 1.60388095238095"
    ## [1] "read:cell ratio for pool1 SortSeq_bin2 is 1.65507102376056"
    ## [1] "read:cell ratio for pool1 SortSeq_bin3 is 1.86083788226244"
    ## [1] "read:cell ratio for pool1 SortSeq_bin4 is 1.95212055089341"
    ## [1] "read:cell ratio for pool2 SortSeq_bin1 is 1.62286825396825"
    ## [1] "read:cell ratio for pool2 SortSeq_bin2 is 1.58901452033396"
    ## [1] "read:cell ratio for pool2 SortSeq_bin3 is 2.06185890170913"
    ## [1] "read:cell ratio for pool2 SortSeq_bin4 is 2.01035380124468"

``` r
#annotate each barcode as to whether it's a homolog variant, SARS-CoV-2 wildtype, synonymous muts only, stop, nonsynonymous, >1 nonsynonymous mutations
dt[,variant_class:=as.character(NA)]
dt[n_codon_substitutions==0, variant_class := "wildtype"]
dt[n_codon_substitutions > 0 & n_aa_substitutions==0, variant_class := "synonymous"]
dt[n_aa_substitutions>0 & grepl("*",aa_substitutions,fixed=T), variant_class := "stop"]
dt[n_aa_substitutions == 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := "1 nonsynonymous"]
dt[n_aa_substitutions > 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := ">1 nonsynonymous"]

#cast the data frame into wide format
dt <- dcast(dt, library + barcode + target + variant_class + aa_substitutions + n_aa_substitutions ~ sample, value.var="count.norm")

#add total count corresponding to count across the four bins for each barcode
dt[,expr_count := sum(SortSeq_bin1,SortSeq_bin2,SortSeq_bin3,SortSeq_bin4),by=c("library","barcode")]

#add indicator if count>1 in >1 bin
dt[,total_bins_w_count := sum(.(SortSeq_bin1,SortSeq_bin2,SortSeq_bin3,SortSeq_bin4)>0),by=c("library","barcode")]
```

## Calculating mean fluorescence

Next, for each barcode, calculate its mean fluorescence as an indicator
of RBD expression level. We will use a maximum likelihood approach to
determine the mean and standard deviation of fluorescence for a barcode,
given its distribution of cell counts across sort bins, and the known
fluorescence boundaries of those sort bins from the sorting log. The
package `fitdistcens` enables this ML estimation for these type of
*censored* observations, where we know we observed a cell within some
fluorescence interval but do not know the exact fluorescence value
attributed to that observation. The counts are multiplied by 20 so that
there is not a large rounding effect when they are rounded to integers.

``` r
#define function to calculate ML meanF
calc.MLmean <- function(b1,b2,b3,b4,min.b1,min.b2,min.b3,min.b4,max.b4,min.count=1){ #b1-4 gives observed cell counts in bins 1-4; remaining arguments give fluorescence boundaries of the respective bins; min.count gives minimum number of total observations needed across bins in order to calculate meanF (default 1)
  data <- data.frame(left=c(rep(min.b1,round(b1)),rep(min.b2,round(b2)),rep(min.b3,round(b3)),rep(min.b4,round(b4))),
                     right=c(rep(min.b2,round(b1)),rep(min.b3,round(b2)),rep(min.b4,round(b3)),rep(max.b4,round(b4)))) #define data input in format required for fitdistcens
  if(nrow(unique(data))>1 & nrow(data)>min.count){ #only fits if above user-specified min.count, and if the data satisfies the fitdistcens requirement that cells are observed in at least two of the censored partitions to enable ML estimation of identifiable parameters
    fit <- fitdistcens(data,"norm")
    return(list(as.numeric(summary(fit)$estimate["mean"]),as.numeric(summary(fit)$estimate["sd"])))
  } else {
    return(list(as.numeric(NA),as.numeric(NA)))
  }
}

#fit ML mean and sd fluorescence for each barcode, and calculate total cell count as the sum across the four bins. Multiply cell counts by a factor of 20 to minimize rounding errors since fitdistcens requires rounding to integer inputs
invisible(dt[library=="pool1",c("expression","ML_sdF") := tryCatch(calc.MLmean(b1=SortSeq_bin1*20,b2=SortSeq_bin2*20,
                                                                      b3=SortSeq_bin3*20,b4=SortSeq_bin4*20,
                                                                      min.b1=log(-642),min.b2=log(131),min.b3=log(558),
                                                                      min.b4=log(2195),max.b4=log(16954)),
                                                          error=function(e){return(list(as.numeric(NA),as.numeric(NA)))}),by=c("library","barcode")])
invisible(dt[library=="pool2",c("expression","ML_sdF") := tryCatch(calc.MLmean(b1=SortSeq_bin1*20,b2=SortSeq_bin2*20,
                                                                      b3=SortSeq_bin3*20,b4=SortSeq_bin4*20,
                                                                      min.b1=log(-701),min.b2=log(151),min.b3=log(659),
                                                                      min.b4=log(2485),max.b4=log(21659)),
                                                          error=function(e){return(list(as.numeric(NA),as.numeric(NA)))}),by=c("library","barcode")])

#save temp data file for downstream troubleshooting since the ML meanF took >1hr to calculate -- don't use these for final anlaysis though for reproducibility!
save(dt,file=paste(config$expression_sortseq_dir,"/dt.temp.Rda",sep=""))
```

## Basic plotting and QC

Let’s look at the distibution of expression scores by variant class for
each library.

``` r
#histogram of mean F, separated by class
hist(dt[variant_class %in% (c("1 nonsynonymous",">1 nonsynonymous")),expression],col="gray40",main="",breaks=50,xlab="ML mean fluorescence (a.u.)")
hist(dt[variant_class %in% (c("synonymous","wildtype")),expression],col="#92278F",add=T,breaks=50)
hist(dt[variant_class %in% (c("stop")),expression],col="#BE1E2D",add=T,breaks=50)
```

<img src="compute_expression_meanF_files/figure-gfm/unfiltered_expression_distribution-1.png" style="display: block; margin: auto;" />

Next let’s look at the distributon of cell counts across the four bins
for each barcode, for those for which estimates were generated on the
bottom.

``` r
#histograms
par(mfrow=c(2,2))
hist(log10(dt[library=="pool1",expr_count]+0.1),xlab="cell count (log10, plus 0.1 pseudocount)",main="pool1, all bc",col="gray50")
hist(log10(dt[library=="pool2",expr_count]+0.1),xlab="cell count (log10, plus 0.1 pseudocount)",main="pool2, all bc",col="gray50")
hist(log10(dt[library=="pool1" & !is.na(expression),expr_count]+0.1),xlab="cell count (log10, plus 0.1 pseudocount)",main="pool1, determined",col="gray50")
hist(log10(dt[library=="pool2" & !is.na(expression),expr_count]+0.1),xlab="cell count (log10, plus 0.1 pseudocount)",main="pool2, determined",col="gray50")
```

<img src="compute_expression_meanF_files/figure-gfm/cell_count_coverage-1.png" style="display: block; margin: auto;" />

Filter out expression measurements determined from \<10 estimated cells

``` r
min_count <- 10
dt[expr_count<min_count, c("expression","ML_sdF","expr_count") := NA]
```

Do as violin plots, faceted by each target. In next notebook, we’ll
evaluate count depth and possibly apply further filtering to remove
low-count expression estimates

``` r
p1 <- ggplot(dt[!is.na(expression),],aes(x=variant_class,y=expression))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("expression")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~target,nrow=4)

grid.arrange(p1,ncol=1)
```

<img src="compute_expression_meanF_files/figure-gfm/expression_distribution_vioplot-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$expression_sortseq_dir,"/violin-plot_meanF-by-target.pdf",sep="")))
```

We have generated expression measurements for 10.7% of the barcodes in
our libraries.

## Data Output

Finally, let’s output our measurements for downstream analyses.

``` r
dt[,.(library,barcode,target,variant_class,aa_substitutions,n_aa_substitutions,
     expr_count,expression)] %>%
  mutate_if(is.numeric, round, digits=6) %>%
  write.csv(file=config$expression_sortseq_file, row.names=F)
```
