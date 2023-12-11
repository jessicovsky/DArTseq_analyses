# gl_filtering.R: steps on filtering DArTseq raw data for analyses.

# ---------------------------------------------------------------------------#
#                                                                            #
#                               DArTseq analyses                             #
#                      By Jessica FR Coelho & Julia T Verba                  #
#                             jessicovsky@gmail.com                          #
#                                   Nov 2023                                 #
#                                                                            #
# ---------------------------------------------------------------------------#
# Check working directory
getwd()

# Get error messages in English
Sys.setenv(lang = "en_US")

# Load libraries
Sys.which("make") #"C:\\rtools40\\usr\\bin\\make.exe"
library(dartR)
library(devtools)
library(terra)

# -------------------------------------------------------------------------- #
#                        Import & Visualize data                             #
# -------------------------------------------------------------------------- #
# Open csv file and convert it to genlight (gl)
dart <- gl.read.dart(filename = "dartseqdata_genomebased/Report_DHare21-6652_SNP_2_bylat.csv",
                     ind.metafile = "dartseqdata_genomebased/SampleFile-DHare21-6652-my_metadata.csv")
dart
# 91 genotypes,  98,313 binary SNPs, 30.35% missing data

# Producing a heatmap with the raw data
glPlot(dart,
       col = c("DeepSkyBlue",
               "DeepPink1",
               "Gold"),
       legend = T)

# Ordering SNPs by callrate then plot heatmap again
dart2 <- dart[ ,order(dart@other$loc.metrics$CallRate, decreasing = TRUE)] 
dart2

# -------------------------------------------------------------------------- #
#                   Drop putative sex-linked SNPs                            #
# -------------------------------------------------------------------------- #
# Import and remove 'outlier' alleles - identified in outlier_removal.R
dropalleles <- read.csv(file = "alleles2remove.csv",
                        header = TRUE)
View(dropalleles) #166 alleles

# Now remove chrs after . of the alleles name to match gl pattern
rname_alleles <- gsub("\\..*", "", dropalleles$rmv_alleles)
rmv_alleles <- cbind(rname_alleles, dropalleles[ ,3:5])
View(rmv_alleles)

# Remove alleles
gl_dropalleles <- gl.drop.loc(dart,
                              loc.list = rmv_alleles$rname_alleles,
                              verbose = 3)
#Summary of recoded dataset
#Original No. of loci: 98313 
#No. of loci deleted: 166 
#No. of loci retained: 98147 

# -------------------------------------------------------------------------- #
#                             Subseting dataset                              #
# -------------------------------------------------------------------------- #
# Harengula thrissina outgroup:
mex <- gl.keep.pop(dart, pop.list = c("MEX"))
mex  #98.19 % missing data -> remove

# Harengula sp.:
bra <- gl.keep.pop(dart, pop.list = c("FNO", "CE", "RN", "PB",
                                      "PE", "AL", "BA", "ABR",
                                      "ES", "RJ", "SP", "SC"))
bra  #89 genotypes,  98,147 SNPs, 28.86% missing data

# --------------------------------------------------------------------------- #
#                     Filtering dataset I: linked-SNPs                        #
#               Nucleotide diversity (pi, Fis) and Tajima's D                 #
#                  CR/ind (t=0.6) > Reproducibility (t=.99)                   #
# --------------------------------------------------------------------------- #
# Call rate per individual
gl0 <- gl.filter.callrate(bra,
                          method = "ind",
                          t = 0.6)
gl0 #86 genotypes, 98,147 binary SNPs, 28.11 % MD
# Removed CE08[CE], CE11[CE], SP02[SP]

# Reproducibility
gl1 <- gl.filter.reproducibility(gl0, t = 0.99)
gl1 #86 genotypes, 88,639 binary SNPs, 29.5 % MD

# Call rate per locus: dataset linked-SNPs-0MD
gl1D <- gl.filter.callrate(gl1, threshold = 1)
gl1D #86 genotypes,	9,369 binary SNPs, 0 % MD

# Export dataset to run DNAsp and calculate Tajima's D
gl2fasta(gl1D,
         method = 1,
         outfile = "DNAsp_again/tajD_DNAsp1.fasta")

# --------------------------------------------------------------------------- #
#                Filtering dataset II: unlinked-best-SNPs                     #
#               Population structure (PCoA, STRUCTURE, Fst)                   #
# CR/ind (t=0.6) > Reproducibility (t=.99) > Monomorphs > Secondaries (best)  #
# --------------------------------------------------------------------------- #
# Monomorphs
gl2 <- gl.filter.monomorphs(gl1, v = 5)
gl2 # 86 genotypes, 85,102 binary SNPs,	30.1 % MD

# Secondaries
gl3 <- gl.filter.secondaries(gl2, method = "best")
gl3 #86 genotypes, 59,992 binary SNPs, 34 % MD

# --------------------------------------------------------------------------- #
#                  Filtering dataset III: unlinked-random-SNPs                #
#                             Demography (dadi)                               #
# CR/ind (t=0.6) > Reproducibility (t=.99) > Monomorphs > Secondaries (random)#
# --------------------------------------------------------------------------- #
# Secondaries
gl4 <- gl.filter.secondaries(gl2, method = "random")
gl4 #86 genotypes, 59,992 binary SNPs, 34.01 % MD

# --------------------------------------------------------------------------- #
#                            Other useful commands                            #
# --------------------------------------------------------------------------- #
gl.report.secondaries(dart) #Average n of SNPs per locus 
gl.report.heterozygosity(gl1, method = "ind") #Heterozygosity

# --------------------------------------------------------------------------- #
#                            jessicovsky@gmail.com                            #
# --------------------------------------------------------------------------- #
# References:
# https://github.com/JMNeves/mugil_dart/blob/main/scripts/dartR_JessikaNeves33.r
# This script was developed under the dartR version 1.1.6 available at: 
# https://cran.r-project.org/src/contrib/Archive/dartR/
