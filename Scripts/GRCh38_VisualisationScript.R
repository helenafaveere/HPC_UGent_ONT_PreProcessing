lib_path <- "/user/gent/484/vsc48405/gONT/Tools/Rlibs"
.libPaths(c(lib_path))

# Load
library(BiocManager)
library(Biobase)
library(QDNAseq)
library(Rsamtools)
remotes::install_github("asntech/QDNAseq.hg38@main")
library(QDNAseq.hg38)
library(ggplot2)
library(plotly)
library(htmlwidgets)
options(bitmapType = "cairo")

# Loading bins
bins <- getBinAnnotations(binSize = 100, genome = "hg38")

args <- commandArgs(trailingOnly = TRUE)
workdir <- args[1]
barcode <- args[2]

# Construct bam_dir dynamically
bam_dir <- file.path(workdir, "basecalling")
output_dir <- file.path(workdir, "QDNAseq")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# List of all BAM files in the directory
bam_files <- list.files(bam_dir, pattern = paste0("sorted_simplex_mapped_barcode", barcode, "\\.bam$"), full.names = TRUE)

# Looping through each BAM file
for (bam_file in bam_files) {
  # Extracting the sample name from BAM file name
  sample_name <- gsub("sorted_simplex_mapped_|\\.bam$", "", gsub(".bam$", "", basename(bam_file)))
  sample_output_dir <- file.path(output_dir, paste0("barcode", barcode))
  dir.create(sample_output_dir, showWarnings = FALSE, recursive = TRUE)

  # Reading the BAM file to create a read count object
  readCounts <- binReadCounts(bins, bamfiles = bam_file)

  # Initial visualisation
  pdf(file.path(sample_output_dir, paste0(sample_name, "_Initial_read_count.pdf")), width = 50, height = 25)
  par(cex = 3)
  plot(readCounts)
  dev.off()

  # Filtering
  pdf(file.path(sample_output_dir, paste0(sample_name, "_highlighted_filters.pdf")), width = 50, height = 25)
  par(cex = 3)
  plot(readCounts, logTransform = TRUE, ylim = c(-2, 20))
  highlightFilters(readCounts, logTransform = TRUE, residual = TRUE, blacklist = TRUE)
  dev.off()

  readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE, chromosomes = "NA")
  
  # Visualisation after filtering
  pdf(file.path(sample_output_dir, paste0(sample_name, "_filtered_logtransform.pdf")), width = 50, height = 25)
  par(cex = 3)
  plot(readCountsFiltered, logtransform = TRUE, ylim = c(-2, 20))
  dev.off()

  # Visualizing the distribution of log2 ratios
  copynumber <- assayDataElement(readCountsFiltered, "counts")
  copynumber <- as.numeric(copynumber)
  log_ratios <- log2(copynumber)
  
  pdf(file.path(sample_output_dir, paste0(sample_name, "_log2_ratios_distribution.pdf")), width = 50, height = 25)
  par(cex = 3)
  plot(log_ratios, ylim = c(-1, 20))
  dev.off()

  # GC and mapability correction
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  pdf(file.path(sample_output_dir, paste0(sample_name, "_noiseplot_estimated_GC_map.pdf")), width = 50, height = 25)
  par(cex = 3)
  noisePlot(readCountsFiltered)
  dev.off()

  readCountsCorrected <- correctBins(readCountsFiltered)
  readCountsSmoothed <- smoothOutlierBins(readCountsCorrected)
  pdf(file.path(sample_output_dir, paste0(sample_name, "_Smoothing.pdf")), width = 50, height = 25)
  par(cex = 3)
  plot(readCountsSmoothed)
  dev.off()

  # Segmentation + normalisation of the segments
  segmentedReadCounts <- segmentBins(readCountsSmoothed)
  normalizedCounts <- normalizeSegmentedBins(segmentedReadCounts)

  # Calling the aberrations and visualising them
  copyNumbersCalled <- callBins(normalizedCounts)
  pdf(file.path(sample_output_dir, paste0(sample_name, "_QDNAseq_CNV_final.pdf")), width = 50, height = 25)
  par(cex = 3)
  plot(copyNumbersCalled)
  dev.off()

  ## Ggplot visualisation
  # Extracting data from the QDNAseq object
  binData <- data.frame(
  chromosome = as.character(fData(copyNumbersCalled)$chromosome),
  position = (fData(copyNumbersCalled)$start + fData(copyNumbersCalled)$end) / 2,
  log2ratio = log2(as.numeric(assayDataElement(copyNumbersCalled, "copynumber"))),  
  segmented = log2(as.numeric(assayDataElement(copyNumbersCalled, "segmented")))  
)

  
  # Removing na values
  binData <- na.omit(binData)
  
  # Including the sex chromosomes + ordering
  binData$chromosome[binData$chromosome %in% c("23", "chrX")] <- "X"
  binData$chromosome[binData$chromosome %in% c("24", "chrY")] <- "Y"
  present_chromosomes <- unique(binData$chromosome)
  chromosome_order <- c(as.character(1:22), "X", "Y")
  chromosome_order <- chromosome_order[chromosome_order %in% present_chromosomes]
  binData$chromosome <- factor(binData$chromosome, levels = chromosome_order)
  
  # Setting the chromosome boundaries for the X-axis
  chr_boundaries <- aggregate(position ~ chromosome, data = binData, max)
  chr_boundaries <- chr_boundaries[order(match(chr_boundaries$chromosome, chromosome_order)), ]
  chr_offset <- c(0, cumsum(chr_boundaries$position)[-nrow(chr_boundaries)])
  names(chr_offset) <- chromosome_order
  binData$cumulative_position <- binData$position + chr_offset[binData$chromosome]
  chr_boundaries$cumulative_position <- chr_boundaries$position + chr_offset[chr_boundaries$chromosome]
  chr_boundaries$midpoint <- (c(0, head(chr_boundaries$cumulative_position, -1)) + chr_boundaries$cumulative_position) / 2
  

  # Threshold for gain/loss
  threshold <- 0.5 
  binData$color_group <- cut(binData$segmented, breaks = c(-Inf, -threshold, threshold, Inf), labels = c("Loss", "Neutral", "Gain"))
  binData$cbs_color_group <- cut(binData$segmented, breaks = c(-Inf, -threshold, threshold, Inf), labels = c("CBS Loss", "CBS Neutral", "CBS Gain"))
  
  # Calculating the number of reads
  total_reads <- sum(assayData(readCountsFiltered)$counts, na.rm = TRUE)

  # Extract the chromosome names from the featureData
  chromosomes <- fData(readCountsFiltered)$chromosome

  # Extract read counts
  sex_counts <- assayData(readCountsFiltered)$counts

  # Compute total reads for X and Y chromosomes
  X_reads <- sum(sex_counts[chromosomes == "X"], na.rm = TRUE)
  Y_reads <- sum(sex_counts[chromosomes == "Y"], na.rm = TRUE)

  # Calculate the Y/X ratio
  Y_over_X_ratio <- Y_reads / X_reads

  # Assign sex based on the threshold
  sex_assignment <- ifelse(Y_over_X_ratio > 0.01, "male", "female")
  
  
  # Outlier marking
  outliers <- binData[binData$log2ratio > 5 | binData$log2ratio < -5, ]
    
  outliers$adjusted_y <- outliers$log2ratio # Start by keeping values unchanged
  outliers$adjusted_y[outliers$log2ratio > 5] <- 5 # Cap values above 5
  outliers$adjusted_y[outliers$log2ratio < -5] <- -5 # Cap values below -5

  outliers$indicator <- "Above/Below Limit"
  outliers$color_group <- "Above/Below Limit" # Assign to the same legend category

  
  # Plotting
  static_plot <- ggplot(binData, aes(x = cumulative_position, y = log2ratio, color = color_group)) +
    geom_point(size = 0.5, alpha = 0.8) +
    geom_point(aes(y = segmented, color = cbs_color_group), size = 0.7) +
    geom_hline(yintercept = c(-threshold, threshold), linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_vline(data = chr_boundaries, aes(xintercept = cumulative_position), linetype = "dashed", color = "black", linewidth = 0.3) +
    scale_color_manual(
      values = c("Gain" = "#0000CD", "Neutral" = "gray", "Loss" = "#CD5C5C", "CBS Gain" = "#00008B", "CBS Neutral" = "#505050",  "CBS Loss" = "#8B0000", "Above/Below Limit" = "black"),
      name = "Copy Number State"
    ) +
    scale_x_continuous(
      breaks = chr_boundaries$midpoint,
      labels = chr_boundaries$chromosome
    ) +
    theme_minimal() +
    labs(
      title = paste("CNV Profile for Sample", sample_name, "(", format(total_reads, big.mark = ","), "reads,", sex_assignment,")"),
      x = "Chromosome",
      y = "Log2 Copy Number Ratio"
    ) +
    theme(
      axis.text.x = element_text(size = 12, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
    
    static_plot <- static_plot + 
        geom_point(data = outliers, aes(x = cumulative_position, y = adjusted_y, color = color_group), shape = 17, size = 3)


  png(file.path(sample_output_dir, paste0(sample_name, "_ggplot.png")), width = 4000, height = 2000, res = 250)
  print(static_plot + ylim(-1,0.5))
  dev.off()

  # Creating the interactive plot and saving as html
  interactive_plot <- ggplotly(static_plot)
  saveWidget(interactive_plot, file.path(sample_output_dir, paste0(sample_name, ".html")))
}