### RScript of the commands described in the tutorial found here: https://astrobiomike.github.io/amplicon/dada2_workflow_ex

##### PROCESSING #####

library(dada2)
packageVersion("dada2") 
path="~/Genomics/Testing_data/dada2_amplicon_ex_workflow"
setwd(path)

list.files() # make sure what we think is here is actually here

## first we're setting a few variables we're going to use ##

# one holding the file names of all the forward reads
forward_reads  <- sort(list.files(path, pattern="_R1_L100.fq.gz", full.names = TRUE))
# and one with the reverse
reverse_reads <- sort(list.files(path, pattern="_R2_L100.fq.gz", full.names = TRUE))

# one with all sample names, by scanning our "samples" file we made earlier
samples <- sapply(strsplit(basename(forward_reads), "_"), `[`, 1)

# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads <-  file.path(path, "filtered",paste0(samples, "_sub_R1_filtered.fq.gz"))
filtered_reverse_reads <-  file.path(path, "filtered",paste0(samples, "_sub_R2_filtered.fq.gz"))


## Quality trimming/filtering
plotQualityProfile(forward_reads[1:3])

plotQualityProfile(reverse_reads[1:3])

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                              rm.phix=TRUE, minLen=175, truncLen=c(250,200))



class(filtered_out) # matrix
dim(filtered_out) # 20 2

filtered_out


## don't worry if the numbers vary a little, this might happen due to different versions being used 
## from when this was initially put together

plotQualityProfile(filtered_reverse_reads[1:3])

## Generating an error model of our data
################################################################# saved object ###############################################

#err_forward_reads <- learnErrors(filtered_forward_reads)

load("err_forward_reads.RData")

################################################################# saved object ###############################################

#err_reverse_reads <- learnErrors(filtered_reverse_reads)

load("err_reverse_reads.RData")




#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance” equal to the number of reads with that unique sequence.
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#Dereplication in the DADA2 pipeline has one crucial addition from other pipelines: DADA2 retains a summary of the quality information associated with each unique sequence.
#The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent sample inference step, significantly increasing DADA2’s accuracy.


derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow

derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derep_reverse) <- samples

## Inferring ASVs
dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo")
# dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread=TRUE) # problem running this way if on Binder
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo")
# dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread=TRUE) # problem running this way if on Binder


## Merging forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, trimOverhang=TRUE, minOverlap=70)

# this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
class(merged_amplicons) # list
length(merged_amplicons) # 20 elements in this list, one for each of our samples
names(merged_amplicons) # the names() function gives us the name of each element of the list 

class(merged_amplicons$B1) # each element of the list is a dataframe that can be accessed and manipulated like any ordinary dataframe

names(merged_amplicons$B1) # the names() function on a dataframe gives you the column names

#########################################################################saved Object #######################################
## Generating a count table
seqtab <- makeSequenceTable(merged_amplicons)

#load("seqtab.RData")

class(seqtab) # matrix
dim(seqtab) # 20 2521

  ## don't worry if the numbers vary a little, this might happen due to different versions being used 
  ## from when this was initially put together

## Chimera identification

#########################################################################saved Object #######################################

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T) # Identified 17 bimeras out of 2521 input sequences.

load("seqtab.nochim.RData")

# though we only lost 17 sequences, we don't know if they held a lot in terms of abundance, this is one quick way to look at that
sum(seqtab.nochim)/sum(seqtab) # 0.9931372 # in this case we barely lost any in terms of abundance

  ## don't worry if the numbers vary a little, this might happen due to different versions being used 
  ## from when this was initially put together


## Overview of counts throughout
# set a little function
getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

summary_tab

  ## don't worry if the numbers vary a little, this might happen due to different versions being used 
  ## from when this was initially put together

write.table(summary_tab, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)


## Assigning taxonomy
## skipping this codeblock for time, and it will not run in the binder environment
## downloading DECIPHER-formatted SILVA v138 reference
 download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile="SILVA_SSU_r138_2019.RData")

## loading reference taxonomy object
# load("SILVA_SSU_r138_2019.RData")

## loading DECIPHER
# library(DECIPHER)
# packageVersion("DECIPHER") # v2.6.0 when this was initially put together, though might be different in the binder or conda installation, that's ok!

## creating DNAStringSet object of our ASVs
# dna <- DNAStringSet(getSequences(seqtab.nochim))

## and classifying
# tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)

# loading output taxonomy object
load("tax-info.RData")


## Extracting the standard goods from DADA2

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

## Loading libraries
  # don't worry if versions are different from what's listed here, shown are are just what was used when this was initially put together

library(phyloseq) ; packageVersion("phyloseq") 

taxa_data<- read.delim("ASVs_taxonomy.tsv")
asv_data<- read.delim("ASVs_counts.tsv")
metadata <-read.delim("sample_info.tsv")

asv_data <- asv_data %>%tibble::column_to_rownames("X") 

taxa_data <- taxa_data %>% tibble::column_to_rownames("X")

metadata <- metadata %>% tibble::column_to_rownames("X") 

asv_data <- as.matrix(asv_data)
taxa_data <- as.matrix(taxa_data)

OTU = otu_table(asv_data, taxa_are_rows = TRUE)
TAX = tax_table(taxa_data)
samples = sample_data(metadata)
phyloseq_object <- phyloseq(OTU, TAX, samples)

phyloseq_object
rank_names(phyloseq_object)
sample_variables(phyloseq_object)


taxa_bar_phylum <- plot_bar(phyloseq_object,  fill="phylum")
ggsave(filename = "taxa_bar_phylum.png", plot = taxa_bar_phylum,width = 20 , height = 15 ,units = "in" ,dpi = 300)

taxa_bar_type <- plot_bar(phyloseq_object, fill="phylum") + facet_wrap(~type, scales= "free_x", nrow=1)
ggsave(filename = "taxa_bar_type.png", plot = taxa_bar_type,width = 20 , height = 15 ,units = "in" ,dpi = 300)
