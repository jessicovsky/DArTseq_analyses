# pop_analyses.R

# ---------------------------------------------------------------------------#
#                                                                            #
#                               DArTseq analyses                             #
#                             By Jessica FR Coelho                           #
#                             jessicovsky@gmail.com                          #
#                                   Apr 2022                                 #
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

# Import image from Julia's computer
load("pop_analyses-workspace.RData")

# Get dataset to run analyses
# Filtered dataset with putative sex-linked alleles removed
gl3p_dropalleles
gl4_dropalleles 

# --------------------------------------------------------------------------- #
#                                 Structure                                   #
# --------------------------------------------------------------------------- #
# Takes note of these:
nLoc(gl3p_dropalleles) 
nInd(gl3p_dropalleles) 

# Set up pop and name lines, STRUCTURE takes numeric pop-IDs
pop_lvls <- levels(gl3p_dropalleles@pop)
Pp <- gl3p_dropalleles@pop
POP = as.numeric(factor(Pp, levels = pop_lvls))
ind.names <- gl3p_dropalleles$ind.names

# Set parameters; if unsure run: ?gl2structure
gl2structure(gl3p_dropalleles,
             indNames = ind.names,
             addcolumns = c(POP),
             ploidy = 2, 
             exportMarkerNames = TRUE,
             outfile = "Hsp_gl3p_dropalleles.str")
# Now open Structure and load the above file

# --------------------------------------------------------------------------- #
#                                     SFS                                     #
# --------------------------------------------------------------------------- #
# Convert genlight to vcf file using radiator 
library(remotes)
remotes::install_github("thierrygosselin/radiator")
library(radiator)
genomic_converter(gl3p_dropalleles, output = 'vcf')

# Go to radiator folder and read vcf
install.packages("vcfR")
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
#                           Summary statisics                                 #
#                      Private alleles, Ho, Hs, Fis                           #
# --------------------------------------------------------------------------- #
library(hierfstat)
library(poppr)

# Private alleles in one population compared with a second population
# for all pops pairwise. Also reports a count of fixed allelic differences
# and the mean absolute allele freq differences between pairs of populations.
gl.report.pa(gl3p_dropalleles)
gl.basic.stats(gl3p_dropalleles)

# Heterozygosity by population, then by individual, then by batches
gl.report.heterozygosity(gl3p_dropalleles)
gl.report.heterozygosity(gl3p_dropalleles, method = "ind")

# For Harengula sp. all sites
hsp_gi <- gl2gi(gl3p_dropalleles) #Convert genlight object to genind object
basic.stats(hsp_gi,
            diploid = TRUE,
            digits = 4)           #calculates basic statistics
boot.ppfis(dat = hsp_gi,
           nboot = 100)           #calculates Fis confidence intervals

# All sites at the Brazilian coast - removed only FNO
coast_gi <- popsub(hsp_gi,
                   sublist = c("CE","RN","PB","PE","AL","BA","ABR",
                               "ES","RJ","SP","SC"),
                   exclude = NULL)
coast_gi #76 individuals; 29,176 loci; 57,438 alleles
basic.stats(coast_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = coast_gi, nboot = 100)

coast_stat <- gl.drop.pop(gl3p_dropalleles, recalc = TRUE, pop.list = c("FNO"))
gl.report.heterozygosity(coast_stat)

# Separate populations (per site)
levels(pop(gl3p_dropalleles))

# FNO
FNO_stat <- gl.keep.pop(gl3p_dropalleles, recalc = TRUE, pop.list = c("FNO"))
FNO_stat@ind.names
gl.report.heterozygosity(FNO_stat)

FNO_gi <- popsub(hsp_gi, sublist = "FNO", exclude = NULL) #Convert to genind
FNO_gi #10 individuals; 28,412 loci; 36,992 alleles
basic.stats(FNO_gi, diploid = TRUE, digits = 4)           #Basic statistics
boot.ppfis(dat = FNO_gi, nboot = 100)                     #Fis conf intervals

# CE
CE_gi <- popsub(hsp_gi, sublist = "CE", exclude = NULL) #Convert to genind
CE_gi #8 individuals; 28,923 loci; 39,656 alleles
basic.stats(CE_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = CE_gi, nboot = 100)

# RN
RN_gi <- popsub(hsp_gi, sublist = "RN", exclude = NULL)
RN_gi #10 individuals; 29,110 loci; 41,831 alleles
basic.stats(RN_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = RN_gi, nboot = 100)

# PB
PB_gi <- popsub(hsp_gi, sublist = "PB", exclude = NULL)
PB_gi #6 individuals; 28,824 loci; 38,229 alleles
basic.stats(PB_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = PB_gi, nboot = 100)

# PE
PE_gi <- popsub(hsp_gi, sublist = "PE", exclude = NULL)
PE_gi #9 individuals; 29,013 loci; 40,648 alleles
basic.stats(PE_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = PE_gi, nboot = 100)

# AL
AL_gi <- popsub(hsp_gi, sublist = "AL", exclude = NULL)
AL_gi #3 individuals; 27,944 loci; 34,256 alleles
basic.stats(AL_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = AL_gi, nboot = 100)

# PE + AL merged - because AL n=3
PEAL_gi <- popsub(hsp_gi, sublist = c("PE", "AL"), exclude = NULL)
PEAL_gi #12 individuals; 61,006 loci; 87,295 alleles
basic.stats(PEAL_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = PEAL_gi, nboot = 100)

# BA
BA_gi <- popsub(hsp_gi, sublist = "BA", exclude = NULL)
BA_gi #5 individuals; 28,595 loci; 37,007 alleles
basic.stats(BA_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = BA_gi, nboot = 100)

# ABR
ABR_gi <- popsub(hsp_gi, sublist = "ABR", exclude = NULL)
ABR_gi #7 individuals; 28,917 loci; 39,172 alleles
basic.stats(ABR_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = ABR_gi, nboot = 100)

# ES
ES_gi <- popsub(hsp_gi, sublist = "ES", exclude = NULL)
ES_gi #8 individuals; 28,963 loci; 40,001 alleles									
basic.stats(ES_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = ES_gi, nboot = 100)

# RJ
RJ_gi <- popsub(hsp_gi, sublist = "RJ", exclude = NULL)
RJ_gi #9 individuals; 29,035 loci; 40,859 alleles									
basic.stats(RJ_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = RJ_gi, nboot = 100)

# SP
SP_gi <- popsub(hsp_gi, sublist = "SP", exclude = NULL)
SP_gi #8 individuals; 28,998 loci; 40,492 alleles									
basic.stats(SP_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = SP_gi, nboot = 100)

# SC
SC_gi <- popsub(hsp_gi, sublist = "SC", exclude = NULL)
SC_gi #3 individuals; 27,482 loci; 33,349 alleles
basic.stats(SC_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = SC_gi, nboot = 100)

# SP + SC merged, SC n=3
SPSC_gi <- popsub(hsp_gi, sublist = c("SP", "SC"), exclude = NULL)
SPSC_gi #11 individuals; 60,893 loci; 86,327 alleles
basic.stats(SPSC_gi, diploid = TRUE, digits = 4)
boot.ppfis(dat = SPSC_gi, nboot = 100)

# -------------------------------------------------------------------------- #
#                                   Fst                                      #
# -------------------------------------------------------------------------- #
library(StAMPP)
bra_fst <- stamppFst(gl3p_dropalleles,
                     nboots = 100,
                     percent = 95,
                     nclusters = 1)
# Results saved in results/pop_metrics.xlsx [all.basic-stats]

# Reassign individuals to new populations: merge PE+AL and SP+SC
# then recalculate Fst
pop(gl3p_dropalleles)

# Use gl4 as backup to not mess up gl3p_dropalleles
gl4 <- gl3p_dropalleles
levels(pop(gl4))

gl4 <- gl.merge.pop(gl4,
                    old = c("SP", "SC"),
                    new = "SP_SC")
gl4 <- gl.merge.pop(gl4,
                    old = c("PE", "AL"),
                    new = "PE_AL")
levels(pop(gl4))

# Calculate Fst with the new reasigned populations
bra2_fst <- stamppFst(gl4,
                      nboots = 100,
                      percent = 95,
                      nclusters = 1)

# Re-calculate heterozygosity
# Each site at the coast + island
gl.report.heterozygosity(gl4)

# H at the coast
bra2 <- gl.merge.pop(gl4,
                     old = c("ABR","PE_AL","BA","CE","ES","PB",
                             "RJ","RN","SP_SC"),
                     new = "bra2")

levels(pop(bra2))
gl.report.heterozygosity(bra2)

# --------------------------------------------------------------------------- #
#                                   DNAsp                                     #
#                            pi, theta and Tajima's D                         #
# --------------------------------------------------------------------------- #
library(stringr)
gl3p_dropalleles@ind.names
gl2fasta(gl3p_dropalleles,
         method = 3,
         outfile = "/result/pi-theta-tajD_DNAsp3.fasta")
# Method 3 - heterozygous positions are replaced by the standard ambiguity codes. 
# The resultant SNP bases are concatenated across loci to generate a single 
# combined sequence to be used in subsequent MP phylogenetic analyses.

gl2fasta(gl3p_dropalleles,
         method = 1,
         outfile = "pi-theta-tajD_DNAsp1.fasta")

# Open this file at DNAsp as "Unfold a FASTA file" and then run the analyses.

# -------------------------------------------------------------------------- #
#                       F stats - Fit, Fst, Fis                              #
# -------------------------------------------------------------------------- #
library(adegenet)
library(pegas)
# https://cran.r-project.org/web/packages/pegas/vignettes/ReadingFiles.pdf

install.packages("hierfstat")
library(hierfstat)
fit_test <- beta.dosage(gl3p_dropalleles,Mb=TRUE)
fit_test2 <- betas(gl3p_dropalleles,
                   nboot=0,
                   lim=c(0.025,0.975),
                   diploid=TRUE,
                   betaijT=TRUE)
gl3p_drop_df <- genind2df(hsp_gi)

# -------------------------------------------------------------------------- #
#                    IBD - Isolation by Distance                             #
# -------------------------------------------------------------------------- #
# Isolation by Distance - all sites
gl.ibd(gl4)

# I want to re-order (latitude) it then re-import to run IBD
write.csv2(gl4_dropalleles, file = "gl4_dropalleles.csv")
bra_coast

# Isolation by Distance - coast (removing FNO)
bra_coast <- gl.drop.pop(gl4_dropalleles, pop.list = "FNO")
coast_ibd <- gl.ibd(bra_coast)

# Re-order by latitude
coast_ordered
class(coast_ordered)

coast_ibd2 <- gl.ibd(coast_ordered)

# --------------------------------------------------------------------------- #
# Matrix of genetic distance
fst_genmat <- read.csv2("fst_genmat.csv", header = TRUE)
View(fst_genmat)

# -------------------------------------------------------------------- #
#                             Circuitscape                             #
# -------------------------------------------------------------------- #
# Prepare occ and Circuitscape present model at Models_LGM folder

# Just to visualize Circuitscape output:
circuit_run_noFNONA <- raster("Circuitscape/cscape_noFNONA_cum_curmap.asc")
plot(circuit_run_noFNONA)

# After running Circuitscape:
# 1. Open folder with results and look for .OUT file
# 2. Save the matrix one as .txt file

# Load distance matrix - file .OUT from Circuitscape
# This cscape run excluded FNO and kept NA values on the raster model
mt_resist_noFNONA <- read.table("Circuitscape/cscape_noFNONA_resistancestxt.txt",
                                header = T, row.names = 1)
mt_resist_dist_noFNONA <- dist(mt_resist_noFNONA)

# Ordered locations manually by latitude
dgen_ord <- read.csv("Circuitscape/dgen_ordered2.csv",
                     header = T, row.names = 1)
dgen_ord_dist <- dist(dgen_ord)

# Mantel test to check if variables are correlated
mymantel_noFNONA <- mantel.randtest(dgen_ord_dist,
                                    mt_resist_dist_noFNONA,
                                    nrepet = 100000);mymantel_noFNONA

# --------------------------------------------------------------------------- #
#                            jessicovsky@gmail.com                            #
# --------------------------------------------------------------------------- #
# This script was developed under the dartR version 1.1.6

# References:
# https://github.com/JMNeves/mugil_dart/blob/main/scripts/dartR_JessikaNeves33.r
# https://cran.r-project.org/src/contrib/Archive/dartR/

# dartR
# file:///C:/Users/Jessica/Desktop/dartRseq_tutorial.pdf

# Tutorial IBD
# https://github.com/green-striped-gecko/dartRworkshop/blob/master/landscape_genetics_lecture.pdf

# Treemix
# https://bitbucket.org/nygcresearch/treemix/wiki/Home
# https://speciationgenomics.github.io/Treemix/
