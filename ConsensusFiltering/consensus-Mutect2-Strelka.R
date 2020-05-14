#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 5) {
  stop("5 arguments must be provided, Mutect file, Strelka file and base file name, plus molecular_id and the ref_molecular_id.n", call.=FALSE)
}

Mutectfile <- args[1]
Strelkafile <- args[2]
baseName <- args[3]
molecular_id <- args[4]
ref_molecular_id <- args[5]

## Load Libraries
library(dplyr)


print("Pull down and process GATK variants")
G <- read.delim(file = Mutectfile, row.names = NULL, header=FALSE,
                  stringsAsFactors = FALSE, skip = 1)
annHead <- read.delim(file = Mutectfile, row.names = NULL,
                      nrows = 1, header=FALSE, stringsAsFactors = FALSE)
annHead <- annHead[is.na(annHead) == F]
extraHead <- c("MISC",	"DP.ANN",	"CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL.GATK",	"FILTER.GATK",	"INFO.GATK",	"FORMAT.GATK",	"NORMAL.GATK",	"TUMOR.GATK")
colnames(G) <- c(annHead, extraHead)
#Remove uninformative columns
G <- dplyr::select(G, -c(MISC, QUAL.GATK, ID, Otherinfo, DP.ANN))

G$AD.Ref.TUMOR.GATK <- as.integer(sub(",.*$", "", sub("^[^:]+:", "", G$TUMOR.GATK))); head( G$AD.Ref.TUMOR.GATK)
G$AD.Alt.TUMOR.GATK <- as.integer(sub(":.*$", "", sub("^[^:]+:[^,]+,", "", G$TUMOR.GATK)));head( G$AD.Alt.TUMOR.GATK)
G$VAF.TUMOR.GATK <- as.numeric(sub(":.*$", "", sub("^[^:]+:[^:]+:", "", G$TUMOR.GATK))); head(G$VAF.TUMOR.GATK)
G$DP.TUMOR.GATK <- as.integer(sub(":.*$", "", sub("^[^:]+:[^:]+:[^:]+:", "", G$TUMOR.GATK)));head(G$DP.TUMOR.GATK)
####
G$AD.Ref.NORMAL.GATK <- as.integer(sub(",.*$", "", sub("^[^:]+:", "", G$NORMAL.GATK))); head(G$AD.Ref.NORMAL.GATK)
G$AD.Alt.NORMAL.GATK <- as.integer(sub(":.*$", "", sub("^[^:]+:[^,]+,", "", G$NORMAL.GATK))); head(G$AD.Alt.NORMAL.GATK)
G$VAF.NORMAL.GATK <- as.numeric(sub(":.*$", "", sub("^[^:]+:[^:]+:", "", G$NORMAL.GATK))); head(G$VAF.NORMAL.GATK)
G$DP.NORMAL.GATK <- as.integer(sub(":.*$", "", sub("^[^:]+:[^:]+:[^:]+:", "", G$NORMAL.GATK))); head(G$DP.NORMAL.GATK)
### Leaves NAs for multiallelic calls :(
G <- Filter(function(y)!all(y == "."), G);
G[G == ""] <- NA
G <- G %>% dplyr::mutate_if(is.factor, as.character)
G$VariantID <- paste(G$CHROM, G$Gene.refGene, G$POS, G$REF, G$ALT, sep = "-")
G <- G[duplicated(G$VariantID)==F,] # take the first instance of any multiallelic calls. Sigh.


print("Pull down and process Strelka variants")
S <- read.delim(file = Strelkafile, row.names=NULL, header = FALSE,
                skip = 1)
annHead <- read.delim(file = Strelkafile, row.names = NULL,
                      nrows = 1, header=FALSE, stringsAsFactors = FALSE)
annHead <- annHead[is.na(annHead) == F]
extraHead <- c("MISC", "DP.ANN", "CHROM", "POS", "END", "REF", "ALT", "QUAL.STR", "FILTER.STR", "INFO.STR", "FORMAT.STR", "NORMAL.STR", "TUMOR.STR")
colnames(S) <- c(annHead, extraHead)

#Remove uninformative columns
S <- dplyr::select(S, -c(MISC, QUAL.STR, Otherinfo, DP.ANN))
S$FORMAT.STR <- sub("GT:", "", S$FORMAT.STR);
S$NORMAL.STR <- sub("./.:", "", S$NORMAL.STR);
S$TUMOR.STR <- sub("./.:", "", S$TUMOR.STR);
S$DP.NORMAL.STR <- as.integer(sub(":.*$", "", S$NORMAL.STR));
S$DP.TUMOR.STR <- as.integer(sub(":.*$", "", S$TUMOR.STR));

S$VariantID <- paste(S$CHR, S$Gene.refGene, S$POS, S$REF, S$ALT, sep = "-")
S <- Filter(function(y)!all(y == "."), S);
S[S == ""] <- NA
S <- S %>% dplyr::mutate_if(is.factor, as.character)

## Deal with flubbed annotation issues where the same variant is now annotated two different ways!  ARGH!
commonCols <- c("VariantID","CHROM", "POS", "REF", "ALT", "Chr", "Start", "End", "Ref", "Alt")
GATKCols <- colnames(G)[!colnames(G) %in% colnames(S)]
STRCols <- colnames(S)[!colnames(S) %in% colnames(G)]

GMerge <- G %>% select(c(commonCols, GATKCols))
SMerge <- S %>% select(c(commonCols, STRCols))
variants <- full_join(GMerge, SMerge)
variants$Type <- ifelse(nchar(variants$REF) == nchar(variants$ALT), "SNV", "INDEL");

reannotate <- left_join(variants, G)
reannotate <- left_join(reannotate, S)
reannotate$molecular_id <- molecular_id
reannotate$ref_molecular_id <- ref_molecular_id



write.table(reannotate, file = paste0(baseName, ".consensus.tsv"), sep = "\t",
            row.names = F, quote = FALSE)


