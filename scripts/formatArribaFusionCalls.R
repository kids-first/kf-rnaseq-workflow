suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("optparse"))


option_list <- list(
  make_option(c("-f", "--fusionfile"),type="character",
              help="Fusion calls from Arriba"),
  make_option(c("-t", "--tumorid"), type="character",
              help="KF tumor id"),
  make_option(c("-o","--outputfile"),type="character",
              help="Formatted fusion calls from [STARfusion | Arriba] (.TSV)")
)

# Get command line options, if help option encountered print help and exit,
opt <- parse_args(OptionParser(option_list=option_list))
inputfile <- opt$fusionfile
tumorid <- opt$tumorid
outputfile <- opt$outputfile

df<-read_tsv(inputfile,col_types = readr::cols(breakpoint1 = readr::col_character(),breakpoint2 = readr::col_character()))
df$`gene1--gene2`<-paste(df$`#gene1`,df$gene2,sep="--")


# write to output file
write.table(df,outputfile,sep="\t",quote=FALSE,row.names = FALSE)
