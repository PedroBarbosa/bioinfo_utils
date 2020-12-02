library(snapcount)
library(argparse)

parser <- ArgumentParser(description='Query snaptron database based on a list of regions in bed format')
parser$add_argument(dest='bedfile', help='Inut file')
args <- parser$parse_args()
INFILE <<- args$bedfile
OUTFILE <<- sub('\\.bed$', '.rds', INFILE) 

bed_to_granges <- function(file){
   df <- read.table(file,
                    header=F,
                    stringsAsFactors=F)
 
   if(length(df) > 6){
      df <- df[,-c(7:length(df))]
   }
 
   if(length(df)<3){
      stop("File has less than 3 columns")
   }
 
   header <- c('chr','start','end','id','score','strand')
   names(df) <- header[1:length(names(df))]
 
   if('strand' %in% colnames(df)){
      df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
   }
 
   library("GenomicRanges")
 
   if(length(df)==3){
      gr <- with(df, GRanges(chr, IRanges(start, end)))
   } else if (length(df)==4){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
   } else if (length(df)==5){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
   } else if (length(df)==6){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
   }
   return(gr)
}

granges <- bed_to_granges(INFILE)

q_sra <- QueryBuilder("srav2", regions=granges)

#Returns junctions whose start and end coordinates are within the boundaries of the regions requested
q_sra <- set_coordinate_modifier(q_sra, Coordinates$Within)

#Returns junction with a minimum number of samples containing it
q_sra <- set_row_filters(q_sra, samples_count >= 20)

#Query junctions and remove intervals with not junctions found
jx_q_sra <- query_jx(q_sra, return_rse = TRUE, split_by_region = TRUE)
jx_q_sra <- jx_q_sra[lengths(jx_q_sra) != 0]

#Save output_object
saveRDS(jx_q_sra, file = OUTFILE)

