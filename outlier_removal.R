# outlier_removal: remove putatively sex-linked SNPs. 
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
#                    Visualize duplicated pattern in PCA                     #
# -------------------------------------------------------------------------- #
gl.pcoa.plot(pc,
             dart,
             xaxis = 1,
             yaxis = 2)

# -------------------------------------------------------------------------- #
#                  Cleaning outlier alleles from dataset                     #
# -------------------------------------------------------------------------- #
# Import PCA with duplicated pattern:
pc <- readRDS("pc_duplicated_gl3p.Rdata")

# Visualize
par(mfrow=c(1,1))
gl.pcoa.plot(pc,
             dart,
             xaxis = 1,
             yaxis = 2)
pc$scores

# Identify the potentially problematic alleles
loadings <- as.data.frame(pc$loadings)
View(loadings)
hist(loadings$Axis2)

# Get the first row (index) as a list of alleles to remove from gl object
all_alleles <- cbind(alleles = rownames(loadings), loadings)
rownames(all_alleles) <- 1:nrow(all_alleles)
View(all_alleles)

# --------------------------------------------------------------------------- #
#                             Thresholds in PC2                               #
#                                module 0.035                                 #
# --------------------------------------------------------------------------- #
# Check and investigate alleles to remove
# Thresholds below and above which to remove alleles
t <- as.numeric(-0.035)
t2 <- as.numeric(0.035)

# Dataframe of alleles to be kept
loadings_keep <- subset(loadings[loadings$Axis2 > t & loadings$Axis2 < t2, ])
View(loadings_keep)
hist(loadings_keep$Axis2)

# Merge the alleles to exclude into a single dataframe
excluded_axis2_down <- subset(loadings[loadings$Axis2 < t, ])
excluded_axis2_up   <- subset(loadings[loadings$Axis2 > t2, ])
nrow(excluded_axis2_down) #134 alleles to remove
nrow(excluded_axis2_up) #32 alleles to remove
alleles2remove <- rbind(excluded_axis2_down, excluded_axis2_up)

# Rename first row
alleles2remove <- cbind(rmv_alleles = rownames(alleles2remove),
                        alleles2remove)
rownames(alleles2remove) <- 1:nrow(alleles2remove)
nrow(alleles2remove) #166 alleles to remove
View(alleles2remove)

# Distribution of alleles to remove according threshold PC 2 = 0.035 module
hist(alleles2remove$Axis2)
write.csv(alleles2remove,
          file = "alleles2remove.csv")

# Now remove chrs after . of the alleles name to match gl pattern
rname_alleles2 <- gsub("\\..*", "", alleles2remove$rmv_alleles)
rmv_alleles <- cbind(rname_alleles2, alleles2remove[ ,2:5])
View(rmv_alleles)

# Remove alleles
gl_new <- gl.drop.loc(dart,
                      loc.list = rmv_alleles$rname_alleles2,
                      verbose = 3)
#Original No. of loci: 29342 
#No. of loci deleted: 166 
#No. of loci retained: 29176 
#No. of individuals: 86 
#No. of populations:  12 

# Do a PCoA with the new genlight object to visualize new pattern
pc_rmvalleles <- gl.pcoa(gl_new, nfactors = 5)
gl.pcoa.plot(pc_rmvalleles,
             gl_new,
             xaxis = 1,
             yaxis = 2)
gl.pcoa.plot(pc_rmvalleles,
             gl_new,
             xaxis = 3,
             yaxis = 4)
gl.pcoa.plot(pc_rmvalleles,
             gl_new,
             xaxis = 1,
             yaxis = 3)

# Heatmap of the 166 removed alleles
gl_new_htmap <- gl.keep.loc(gl_new,
                            loc.list = rmv_alleles$rname_alleles2,
                            verbose = 3)
glPlot(gl_new_htmap,
       col = c("DeepSkyBlue",
               "DeepPink1",
               "Gold"),
       legend = T)
saveRDS(gl_new_htmap, file = "gl_new_htmap.Rdata")

# --------------------------------------------------------------------------- #
#                            jessicovsky@gmail.com                            #
# --------------------------------------------------------------------------- #
#References:
#https://github.com/JMNeves/mugil_dart/blob/main/scripts/dartR_JessikaNeves33.r
#This script was developed under the dartR version 1.1.6 available at: 
#https://cran.r-project.org/src/contrib/Archive/dartR/
