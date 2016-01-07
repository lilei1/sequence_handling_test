# seqqs produces three output files
# 1) nucleotide table
# 2) read length by position
# 3) a Phred quality matrix by nucleotide position

#   To take command line arguments
args <- commandArgs(TRUE)
#   This creates a vector of character strings for arguments
#   we will just take two arguments here, the stats directory
#   and the sample name
statsdir <- args[1]
samplename <- args[2]
#   This is the 'YlOrRd' sequential 9-level palette from RColorBrewer
heatmap_colors <- c("#fff7ec", 
    "#fee8c8", 
    "#fdd49e", 
    "#fdbb84", 
    "#fc8d59", 
    "#ef6548", 
    "#d7301f", 
    "#b30000", 
    "#7f0000")

#####
#   Define some plotting functions to generate our plots
#####
#       Base Composition Plot
#           Take two positional arguments:
#               1) Dataframe of nucleotide counts by read position
#               2) Whether the reads are raw or trimmed
#           Other arguments (with defaults)
#               1) margin to apply on upper bound for plot
#               2) margin to apply on lower bound of plot
BaseCompositionPlot <- function(nuc, trim, upper_margin = 0.02, lower_margin = 0.02)
{
    #   Set up some plot parameters
    #       Character scaling, for size of points/text
    cex <- 0.3
    #       Plot character. 19 is a filled point
    pch <- 19
    #       Line type. l (lowercase L) is a solid line with no breaks
    ltype <- "l"
    #       Line width
    lwd <- 0.5
    #       The labels
    xlab <- "Position in Read (bp)"
    ylab <- "Nucleotide Proportion"
    #   Since this is the first plot in the grid, we have to label which column
    #   contains trimmed and which contains raw data. This is a bit kludge-y, but
    #   it works for now...
    main <- paste(trim, "Base Composition", sep="\n")
    #   Calculate the total number of bases at each position
    total_bases <- apply(nuc, 1, sum)
    #   Calculate the proportion of each nucleotide, ignoring N
    nuc.A <- nuc$A / total_bases
    nuc.C <- nuc$C / total_bases
    nuc.G <- nuc$G / total_bases
    nuc.T <- nuc$T / total_bases
    #   Set the upper and lower bounds for the plot
    upper_bound <- max(nuc.A, nuc.C, nuc.G, nuc.T) + upper_margin
    lower_bound <- min(nuc.A, nuc.C, nuc.G, nuc.T) - lower_margin
    #   build the plot!
    #       Start with A
    plot(nuc.A, 
        col="green", 
        ylim=c(lower_bound, upper_bound), 
        cex=cex, 
        xlab=xlab, 
        ylab=ylab, 
        main=main, 
        pch=pch)
    #   Then add points for the other nucleotides
    points(nuc.C, col="blue", cex=cex, pch=pch)
    points(nuc.G, col="black", cex=cex, pch=pch)
    points(nuc.T, col="red", cex=cex, pch=pch)
    #   Add lines connecting them
    lines(nuc.A, col="green", type="l", lwd=lwd)
    lines(nuc.C, col='blue', type="l", lwd=lwd)
    lines(nuc.G, col='black', type="l", lwd=lwd)
    lines(nuc.T, col='red', type="l", lwd=lwd)
    #   Put a line at 0.25 for the expected proportion in a random sequence
    #   note that this isn't the case for exome capture or repetitive sequences
    abline(h=0.25, col="orange")
    #   We need a legend, too
    legend(90, 
        upper_bound, 
        c("A", "C", "G", "T"), 
        pch=pch, 
        col=c("green", "blue", "black", "red"), 
        cex=0.7)
}

#       Length distribution plot
#           Takes two positional arguments:
#               1) The data frame containing numbers of reads at various lengths
#           Takes optional argument:
#               1) axis label step for plotting
LengthDistributionPlot <- function(len, step=5)
{
    #   Set up some plotting parameters
    xlab <- "Read Length"
    ylab <- "log(Read Count)"
    main <- "Read Length Distribution"
    #       A vector for spacing between bars. We set this to 0
    spaces <- rep(0, length(len))
    #   Get the log of the counts of reads
    counts <- log(len$count)
    #   If a class has 0 reads, then it causes some problems, as plot() doesn't seem to 
    #   handle -Inf values well. We set these to 0, even though it's not mathematically correct
    counts[counts == -Inf] <- 0
    #   Generate a barplot
    #   the return value for barplot() is the positions of the bars, so we save those
    at <- barplot(counts, 
        xlab=xlab, 
        ylab=ylab, 
        main=main, 
        col="blue", 
        space=spaces, 
        border="white")
    #   Generate the spacing and labels for the barplot
    #       From 1 to the length of the reads
    at.labs <- seq_along(at)
    #       Take the labels that occur every [STEP] numbers, and add the last class
    at.labs <- c(at.labs[seq(1, length(at.labs), step)], length(at.labs))
    #       And generate the positions of the labels
    at.pos <- c(at[seq(1, length(at), step)], at[length(at)])
    #       Throw on the axis
    axis(1, at=at.pos, labels=as.character(at.labs))
}

#       Heatmap of base qualities
#           Takes one argument:
#               1) The quality matrix
BaseQualHeatmap <- function(qual, colors)
{
    #   set up some plotting parameters
    xlab <- "Read Position (bp)"
    ylab <- "Base Quality"
    main <- "Base Quality Heatmap"
    #   We cast it to a matrix so that we can use image() to build a heatmap
    qual.mat <- as.matrix(qual)
    #   Then we remove columns (qualities) that have 0 bases in them, since they do not inform us
    qual.mat <- qual.mat[,apply(qual.mat, 2, sum) > 0]
    #   Save the remaining column names, as they are our quality scores
    qscores <- colnames(qual.mat)
    #   Then we start building the image
    #       This first call builds the base of the heatmap, with axis labels
    #       and title string
    image(qual.mat,
        col=colors,
        xlab=xlab,
        ylab=ylab,
        main=main,
        axes=F)
    #   Then we add our own axes
    #       x-axis. number of rows in our matrix corresponds to the read length
    axis(1, at=seq(0, 1, length.out=nrow(qual.mat)), labels=as.character(1:nrow(qual.mat)))
    #       y-axis. We use the number of quality scores and the character vector
    #       of the scores themselves to build the axis
    axis(2, at=seq(0, 1, length.out=length(qscores)), labels=qscores)
    #   Put a box on it
    box()
}

#####define a empty function!
#EmptyPlot <- function(x, y, upper_margin = 0.02, lower_margin = 0.02){
#    plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10))
#}
 



#####
#   Forward plots
#####
outputfile <- paste(statsdir, "/plots/", samplename, "_Forward_SeqqsPlots.pdf", sep="")
pdf(file=outputfile, width=8.8, 11)

#   Set it up so that we have multiple plots on one graphics device
par(mfrow=c(3, 2))
#   They fill top-to-bottom, left-to-right
#####
#   Base composition plots
#####
#       Forward before trimming
raw.forward.nucl.name <- paste(statsdir, "/raw_", samplename, "_R1_nucl.txt", sep="")
nucl <-try( read.table(raw.forward.nucl.name,sep="\t", header=TRUE))
if(class(nucl)!='try-error'){
    BaseCompositionPlot(nucl, "Raw")
}else{
    plot(1,type='n',axes=FALSE,ann=FALSE)
}
#       Forward after trimming
Trimmed.forward.nucl.name <- paste(statsdir, "/trimmed_", samplename, "_R1_nucl.txt", sep="")
nucl <- try(read.table(Trimmed.forward.nucl.name,sep="\t",header=TRUE))
if(class(nucl)!='try-error'){
BaseCompositionPlot(nucl, "Trimmed")
}else{
    plot(1,type='n',axes=FALSE,ann=FALSE)
}
#####
#   Length distribution
#####
#       Before trimming
raw.forward.len.name <- paste(statsdir, "/raw_", samplename, "_R1_len.txt", sep="")
len <- try(read.table(raw.forward.len.name,sep="\t", header=TRUE))
if(class(len)!='try-error'){
LengthDistributionPlot(len)
}else{
    plot(1,type='n',axes=FALSE,ann=FALSE)
}

#       After trimming
Trimmed.forward.len.name <- paste(statsdir, "/trimmed_", samplename, "_R1_len.txt", sep="")
len <- try(read.table(Trimmed.forward.len.name,sep="\t", header=TRUE))
if(class(len)!='try-error'){
   LengthDistributionPlot(len)
}else{
    plot(1,type='n',axes=FALSE,ann=FALSE)
}

#####
#   Quality heatmaps
#####
#       Before trimming
raw.forward.qual.name <- paste(statsdir, "/raw_", samplename, "_R1_qual.txt_adj", sep="")
qual <- try(read.table(raw.forward.qual.name,sep="\t", header=TRUE))
if(class(qual)!='try-error'){
BaseQualHeatmap(qual, heatmap_colors)
}else{
    plot(1,type='n',axes=FALSE,ann=FALSE)
}

#       After trimming
Trimmed.forward.qual.name <- paste(statsdir, "/trimmed_", samplename, "_R1_qual.txt_adj", sep="")
qual <- try(read.table(Trimmed.forward.qual.name,sep="\t",header=TRUE))
if(class(qual)!='try-error'){
BaseQualHeatmap(qual, heatmap_colors)
}else{
    plot(1,type='n',axes=FALSE,ann=FALSE)
}

dev.off()


#####
#   Reverse plots
#####
outputfile <- paste(statsdir, "/plots/", samplename, "_Reverse_SeqqsPlots.pdf", sep="")
pdf(file=outputfile, width=8.8, height=11)
#   Set it up so that we have multiple plots on one graphics device
par(mfrow=c(3, 2))

#####
#   Base composition plots
#####
#       Forward before trimming
raw.forward.nucl.name <- paste(statsdir, "/raw_", samplename, "_R2_nucl.txt", sep="")
nucl <- try(read.table(raw.forward.nucl.name,sep="\t", header=TRUE))
if(class(nucl)!='try-error'){
BaseCompositionPlot(nucl, "Raw")
}else{
    plot(1,type='n',axes=FALSE,ann=FALSE)
}
#       Forward after trimming
Trimmed.forward.nucl.name <- paste(statsdir, "/trimmed_", samplename, "_R2_nucl.txt", sep="")
nucl <- try(read.table(Trimmed.forward.nucl.name,sep="\t", header=TRUE))
if(class(nucl)!='try-error'){
BaseCompositionPlot(nucl, "Trimmed")
}else{
    plot(1,type='n',axes=FALSE,ann=FALSE)
}

#####
#   Length distribution
#####
#       Before trimming
raw.forward.len.name <- paste(statsdir, "/raw_", samplename, "_R2_len.txt", sep="")
len <- try(read.table(raw.forward.len.name,sep="\t", header=TRUE))
if(class(len)!='try-error'){
LengthDistributionPlot(len)
}else{
    plot(1,type='n',axes=FALSE,ann=FALSE)
}
#       After trimming
Trimmed.forward.len.name <- paste(statsdir, "/trimmed_", samplename, "_R2_len.txt", sep="")
len <- try(read.table(Trimmed.forward.len.name,sep="\t", header=TRUE))
if(class(len)!='try-error'){
LengthDistributionPlot(len)
}else{
    plot(1,type='n',axes=FALSE,ann=FALSE)
}
#####
#   Quality heatmaps
#####
#       Before trimming
raw.forward.qual.name <- paste(statsdir, "/raw_", samplename, "_R2_qual.txt_adj", sep="")
qual <- try(read.table(raw.forward.qual.name,sep="\t", header=TRUE))
if(class(qual)!='try-error'){
BaseQualHeatmap(qual, heatmap_colors)
}else{
    plot(1,type='n',axes=FALSE,ann=FALSE)
}
#       After trimming
Trimmed.forward.qual.name <- paste(statsdir, "/trimmed_", samplename, "_R2_qual.txt_adj", sep="")
qual <- try(read.table(Trimmed.forward.qual.name,sep="\t",header=TRUE))
if(class(qual)!='try-error'){
BaseQualHeatmap(qual, heatmap_colors)
}else{
    plot(1,type='n',axes=FALSE,ann=FALSE)
}
dev.off()
