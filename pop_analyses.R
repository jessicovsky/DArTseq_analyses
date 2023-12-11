# pop_analyses.R 
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

# Get messages in English
Sys.setenv(lang = "en_US")

# Load libraries
Sys.which("make") #"C:\\rtools40\\usr\\bin\\make.exe"
library(dartR)
library(devtools)
library(terra)
library(adegenet)
library(pegas)
library(hierfstat)
library(poppr)

# Get datasets to run analyses:
# unliked-best-SNPs: PCA, Admixture, Fst, IBD
# linked-SNPs: pi-SNP, Fis 
# linked-SNPs-0MD: Tajima's D (export gl1D and run in DNAsp)
# unlinked-random-SNPs: demography (export and run dadi)

# --------------------------------------------------------------------------- #
#                                   PCA                                       #
# --------------------------------------------------------------------------- #
pc <- gl.pcoa(gl3, nfactors = 5)
gl.pcoa.plot(pc,
             gl3,
             xaxis = 1,
             yaxis = 2)

# -------------------------------------------------------------------------- #
#                                 ADMIXTURE                                  #
# -------------------------------------------------------------------------- #
gl3_admix <- gl2plink(gl3,
                      plink_path = getwd(),
                      outfile = "gl3_admix_plink",
                      outpath = tempdir())

# -------------------------------------------------------------------------- #
#                                   Fst                                      #
# -------------------------------------------------------------------------- #
# Reassign all coastal individuals to a single population
bra_fst <- gl.merge.pop(gl3, old = c("CE", "RN", "PB",
                                     "PE", "AL", "BA", "ABR",
                                     "ES", "RJ", "SP", "SC"),
                        new = "COAST")
gl3_fst <- gl.fst.pop(bra_fst)
gl3_fst 

# -------------------------------------------------------------------------- #
#                    IBD - Isolation by Distance                             #
# -------------------------------------------------------------------------- #
# Isolation by Distance along the Brazilian coast
bra_coast <- gl.keep.pop(gl3, pop.list = c("CE", "RN", "PB",
                                           "PE", "AL", "BA", "ABR",
                                           "ES", "RJ", "SP", "SC"))
bra_coast_ibd <- gl.ibd(bra_coast,
                        plot_theme = theme_classic(geom_smooth(method="auto",
                                                               color="orange")))

# -------------------------------------------------------------------------- #
#              Nucleotide diversity (pi-SNP): run in Pixy                    #
# -------------------------------------------------------------------------- #
# Specify the positions from within dartR object, convert to vcf, then run pixy
snp_pos <- gl1@other$loc.metrics$SnpPosition
snp_chr <- gl1@other$loc.metrics$CloneID
piSNP_pixy <- gl2vcf(gl1
                     plink_path = getwd(),
                     snp_pos,
                     snp_chr,
                     outpath = getwd())

# --------------------------------------------------------------------------- #
#                       Inbreeding coefficient: Fis                           #
# --------------------------------------------------------------------------- #
# Convert genlight object to genind
gl1_gi <- gl2gi(gl1)

# Calculate Fis per individual/population
gl1_stat = basic.stats(gl1_gi,
                       diploid = TRUE,
                       digits = 4)
apply(gl1_stat$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)
apply(gl1_stat$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)
gl1_stat

# Calculate Fis confidence intervals
boot.ppfis(dat = gl1_gi,
           nboot = 100)

# --------------------------------------------------------------------------- #
#                                     SFS                                     #
# --------------------------------------------------------------------------- #
# Convert genlight to vcf file using radiator 
library(remotes)
remotes::install_github("thierrygosselin/radiator")
library(radiator)
genomic_converter(gl3p_dropalleles, output = 'vcf')

# Go to radiator folder and read vcf
library(vcfR)
vcf <- read.vcfR('radiator/radiator_data_20220405@1045.vcf', verbose = T);vcf
#86 samples, 1 CHROMs, 29,176 variants, 21.14% MD

pega <- vcfR2DNAbin(vcf,
                    extract.indels = T,
                    consensus = T,
                    extract.haps = F,
                    unphased_as_NA = F,
                    asterisk_as_del = F,
                    verbose = T)
# After extracting indels, 29176 variants remain

# SFS plot
library(pegas)
sfs <- site.spectrum(pega, folded = T)

# S3 method for class 'spectrum'
par(mfrow = c(1,1))
plot(sfs, col = "red")

# --------------------------------------------------------------------------- #
#                            jessicovsky@gmail.com                            #
# --------------------------------------------------------------------------- #
# This script was developed under the dartR version 1.1.6
# References:
# https://github.com/JMNeves/mugil_dart/blob/main/scripts/dartR_JessikaNeves33.r
# https://cran.r-project.org/src/contrib/Archive/dartR/
# Tutorial IBD
# https://github.com/green-striped-gecko/dartRworkshop/blob/master/landscape_genetics_lecture.pdf
