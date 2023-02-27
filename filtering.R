# gl_filtering.R appliers filtering steps on DArTseq raw data for analyses.

# ---------------------------------------------------------------------------#
#                                                                            #
#                               DArTseq analyses                             #
#                             By Jessica FR Coelho                           #
#                             jessicovsky@gmail.com                          #
#                                   Mar 2022                                 #
#                                                                            #
# ---------------------------------------------------------------------------#
# Get working directory
getwd()

# Get messages in English
Sys.setenv(lang = "en_US")

# Load the package then libraries
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
# 91 genotypes,  98,313 binary SNPs, 30% missing data

# Producing a heatmap with the raw data
glPlot(dart,
       col = c("DeepSkyBlue",
               "DeepPink1",
               "Gold"),
       legend = T)

# Ordering SNPs by callrate then plot heatmap
dart2 <- dart[ ,order(dart@other$loc.metrics$CallRate, decreasing = TRUE)] 
dart2

# -------------------------------------------------------------------------- #
#                             Subseting data                                 #
# -------------------------------------------------------------------------- #
# Harengula thrissina outlier:
mex <- gl.keep.pop(dart, pop.list = c("MEX"))
mex  #98.19 % missing data -> remove

# Harengula sp.:
bra <- gl.keep.pop(dart, pop.list = c("FNO", "CE", "RN", "PB",
                                      "PE", "AL", "BA", "ABR",
                                      "ES", "RJ", "SP", "SC"))
bra  #89 genotypes, 98,313 binary SNPs, 28.82 % missing data

# --------------------------------------------------------------------------- #
#                             Filtering data                                  #
# CR/ind (t=0.6) > Reproducibility (t=.99) > Monomorphs > Unlinked > Paralogs #
# --------------------------------------------------------------------------- #
# Call rate per individual
gl0 <- gl.filter.callrate(bra,
                          method = "ind",
                          t = 0.6)
gl0 #86 genotypes, 98,313 binary SNPs, 28.06 % MD
# Removed CE08[CE], CE11[CE], SP02[SP]

# Reproducibility
gl1 <- gl.filter.reproducibility(gl0, t = 0.99)
gl1 #86 genotypes, 88,747 binary SNPs, 29.47 % MD

# Monomorphs
gl2 <- gl.filter.monomorphs(gl1, v = 5)
gl2 #86 genotypes, 85,210 binary SNPs, 30.07 % MD

# Secondaries
gl3 <- gl.filter.secondaries(gl2, method = "best")
gl3 #86 genotypes, 60,071 binary SNPs, 33.96 % MD

# Paralogs
gl4 <- gl.filter.hamming(gl3, threshold = 0.2, pb = T, v = 5)
gl4 <- gl3p #Both the same; I just renamed to keep pattern
gl4 #86 genotypes, 29,342 binary SNPs, 21.04 % MD

# After removing 'outlier' alleles
gl3p_dropalleles #86 genotypes, 29,176 binary SNPs, 21.14 % MD

# --------------------------------------------------------------------------- #
#                    Average n of SNPs per locus                              #
# --------------------------------------------------------------------------- #
gl.report.secondaries(dart2)

# --------------------------------------------------------------------------- #
#                                   PCoA                                      #
# --------------------------------------------------------------------------- #
pc_split <- gl.pcoa(gl4, nfactors = 5)
gl.pcoa.plot(pc_split,
             gl4,
             xaxis = 1,
             yaxis = 2)

# Remove alleles driving localities to split in outlier_removal.R
# --------------------------------------------------------------------------- #
#                            jessicovsky@gmail.com                            #
# --------------------------------------------------------------------------- #
# References:
# https://github.com/JMNeves/mugil_dart/blob/main/scripts/dartR_JessikaNeves33.r
# This script was developed under the dartR version 1.1.6 available at: 
# https://cran.r-project.org/src/contrib/Archive/dartR/