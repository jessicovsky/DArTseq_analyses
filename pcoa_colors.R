# ---------------------------------------------------------------------------#
#                                                                            #
#                               DArTseq analyses                             #
#                               PCoA color plots                             #
#                               Jessica FR Coelho                            #
#                             jessicovsky@gmail.com                          #
#                                   Feb 2022                                 #
#                                                                            #
# ---------------------------------------------------------------------------#
# Run before PCoA analyses
# Load libraries
library(dartR)
library(ggplot2)
library(adegenet)
library(directlabels)

gl.pcoa.plot <- function(glPca,
                         data, 
                         scale = FALSE,
                         ellipse = FALSE,
                         p = 0.95,
                         labels = "pop",
                         hadjust = 1.5, 
                         vadjust = 1,
                         xaxis = 1,
                         yaxis = 2) {
  if(class(glPca)!="glPca" | class(data)!="genlight") {
    cat("Fatal Error: glPca and genlight objects required for glPca and data parameters respectively!\n");
    stop()
  }
  
  if (labels == "smart") {
    hadjust <- 0; vadjust <- 0
  }
  
  # Create a dataframe to hold the required scores
  m <- cbind(glPca$scores[ ,xaxis], glPca$scores[ ,yaxis])
  df <- data.frame(m)
  
  # Convert the eigenvalues to percentages
  s <- sum(glPca$eig)
  e <- round(glPca$eig*100/s,1)
  
  # Labels for the axes
  xlab <- paste("PCoA Axis", xaxis, "(",e[xaxis],"%)")
  ylab <- paste("PCoA Axis", yaxis, "(",e[yaxis],"%)")
  
  if (labels == "none" | labels == FALSE) {
    cat("Plotting points with no labels\n")
  } else {
    pop <- factor(pop(data))
    df <- cbind(df, pop)
    colnames(df) <- c("PCoAx", "PCoAy", "pop")
    
    # Lat order: CE > RN > PB > PE > AL > BA > ES > RJ > SP > SC
    # Islands: FNO > ABR
    # Plot
    p <- ggplot(df, aes(x = PCoAx, y = PCoAy, colour = pop)) +
      geom_point(shape = 19,
                 size = 4,
                 colour = c("#c51b7d", #ABR
                            "#fee08b", #AL
                            "#e6f598", #BA
                            "#9e0142", #CE
                            "#abdda4", #ES
                            "#4d4d4d", #FNO
                            "#f46d43", #PB
                            "#fdae61", #PE
                            "#66c2a5", #RJ
                            "#d53e4f", #RN
                            "#5e4fa2", #SC
                            "#3288bd")[pop], aes(colour = pop)) + #SP
      #geom_dl(aes(label = pop), method = "first.points") +
      #ggtitle(paste("PCoA Plot")) +
      theme(axis.title = element_text(face = "bold.italic",
                                      size = "10",
                                      color = "black"),
            axis.text.x = element_text(face = "bold",
                                       angle = 0,
                                       vjust = 0.5,
                                       size = 14),
            axis.text.y = element_text(face = "bold",
                                       angle = 0,
                                       vjust = 0.5,
                                       size = 14),
            legend.title = element_text(colour = "black",
                                        size = 10,
                                        face = "bold"),
            legend.text = element_text(colour = "black",
                                       size = 10,
                                       face = "bold")) +
      labs(x = xlab, y = ylab) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      theme(legend.position = "left") +
      # Scale the axes in proportion to % explained, if requested
      if(scale == TRUE) {
        p <- p + coord_fixed(ratio = e[yaxis]/e[xaxis])
      }
    # Add ellipses if requested
    if(ellipse == TRUE) {
      p <- p + stat_ellipse(aes(colour = pop), type = "norm", level = 0.95)
    }
  }
  return (p)
}

# --------------------------------------------------------------------------- #
#                            jessicovsky@gmail.com                            #
# --------------------------------------------------------------------------- #
# Original script: https://github.com/JMNeves/mugil_dart/blob/main/scripts/
# > PCoA_color_JessikaNeves.R#L1, 
# developed under the dartR version 1.1.6 available at 
# https://cran.r-project.org/src/contrib/Archive/dartR/