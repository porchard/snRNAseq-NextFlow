suppressPackageStartupMessages(library("DropletUtils"))
library("optparse")

option_list <- list(
    make_option(c("--donor"), type = "character", help = "[Required] Library ID."),
    make_option(c("--barcodeList"), type = "character", help = "Absolute path containng matrix.txt, features.tsv and barcodes.tsv files."),
    make_option(c("--cbMetrics"), type = "character", help = "CellBender metrics file, often in the form *.cellbender_FPR_*_metrics.csv."),
    make_option(c("--lowerForKnee"), type = "integer", help = "[Required] Lower cut-off, used for EmptyDrops."),
    make_option(c("--outKnee"), type = "character", help = "File to save knee and inflection points."),
    make_option(c("--outPass"), type = "character", help = "File to save significant cells.")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

donor <- opts$donor
inputDir <- opts$barcodeList
cbMetrics <- opts$cbMetrics
outKnee <- opts$outKnee
outPass <- opts$outPass
lowerForKnee <- opts$lowerForKnee
lowerForKnee <- as.numeric(lowerForKnee)

knee_inflection_rank <- function(m, lower=50, fit.bounds=NULL, exclude.from=50, df=20) {
    #code adapted from https://github.com/MarioniLab/DropletUtils/blob/devel/R/barcodeRanks.R
    #m is a SummarizedExperiment containing such a matrix. (same as in EmptyDrops)
    #lower is a numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.

    totals <- unname(colSums(counts(m)))
    o <- order(totals, decreasing=TRUE)

    stuff <- rle(totals[o])
    run.rank <- cumsum(stuff$lengths) - (stuff$lengths-1)/2 # Get mid-rank of each run.
    run.totals <- stuff$values

    keep <- run.totals > lower
    if (sum(keep)<3) {
        stop("insufficient unique points for computing knee/inflection points")
    }
    y <- log10(run.totals[keep])
    x <- log10(run.rank[keep])

    # Numerical differentiation to identify bounds for spline fitting.
    edge.out <- find_curve_bounds(x=x, y=y, exclude.from=exclude.from)
    left.edge <- edge.out["left"]
    right.edge <- edge.out["right"]

    # As an aside: taking the right edge to get the total for the inflection point.
    # We use the numerical derivative as the spline is optimized for the knee.
    inflection <- 10^(y[right.edge])
    inflection_rank <- 10^(x[right.edge])

    # We restrict curve fitting to this region, thereby simplifying the shape of the curve.
    # This allows us to get a decent fit with low df for stable differentiation.
    if (is.null(fit.bounds)) {
        new.keep <- left.edge:right.edge
    } else {
        new.keep <- y > log10(fit.bounds[1]) & y < log10(fit.bounds[2])
    }

    # Smoothing to avoid error multiplication upon differentiation.
    # Minimizing the signed curvature and returning the total for the knee point.
    fitted.vals <- rep(NA_real_, length(keep))

    if (length(new.keep) >= 4) {
        fit <- smooth.spline(x[new.keep], y[new.keep], df=df)
        fitted.vals[keep][new.keep] <- 10^fitted(fit)

        d1 <- predict(fit, deriv=1)$y
        d2 <- predict(fit, deriv=2)$y
        curvature <- d2/(1 + d1^2)^1.5
        knee <- 10^(y[new.keep][which.min(curvature)]) #which.min(curvature) is knee
	knee_rank <- 10^(x[new.keep][which.min(curvature)])
    } else {
        # Sane fallback upon overly aggressive filtering by 'exclude.from', 'lower'.
        knee <- 10^(y[new.keep[1]])
	knee_rank <- 10^(x[new.keep[1]]) 
    }

    res <- list(inflection_rank, knee_rank)
    names(res) <- c("inflection_rank", "knee_rank")

    return(res)
}

endCliff <- function(m, inflection_rank, lowerForKnee=100, lowerRankForEndCliff=1000, fit.bounds=NULL, exclude.from=50, df=20) {
    #lowerForKnee required to be the same as `lower` used in barcodeRanks(). lower is a numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.
    #lowerRankForEndCliff is max *rank* of the barcodes to include, and anything with higher *rank* is assumed to be empty droplets.
    ##Should use *rank* of barcodes around plateau points that can be obtained from CellBender. Note that lowerRankForEndCliff is *barcode rank*, which is different from lowerForKnee.
    
    #m is a SummarizedExperiment containing such a matrix. (same as in EmptyDrops)
    
    totals <- unname(colSums(counts(m)))
    o <- order(totals, decreasing=TRUE)

    stuff <- rle(totals[o])
    run.rank <- cumsum(stuff$lengths) - (stuff$lengths-1)/2 # Get mid-rank of each run.
    run.totals <- stuff$values

    lower <- lowerRankForEndCliff
    if (totals[lower] < 100) { #if cellbender predicts plateau points too low, force to find end_cliff above point with UMIs=100, assuming anything with UMIs < 100 are empty barcodes - to avoid discreteness
           plateau <- 100
           keep <- run.totals > plateau
    } else if (totals[lower] > totals[inflection_rank]) {
	    plateau <- 100
	    keep <- run.totals > plateau
    } else
    {
	    plateau <- totals[lower]
	    keep <- run.totals > plateau
    }

    if (sum(keep)<3) {
        stop("insufficient unique points for computing end cliff points")
    }
    y <- log10(run.totals[keep])
    x <- log10(run.rank[keep])

    # Numerical differentiation to identify bounds for spline fitting.
    edge.out <- find_curve_bounds(x=x, y=y, exclude.from=inflection_rank) #exclude.from=inflection_rank to use points from inflection to end
    left.edge <- edge.out["left"]
    right.edge <- edge.out["right"]

    new.keep <- (right.edge):length(which(keep == TRUE))  #if using points from knee to end, then use (which.min(curvature)) instead of right.edge. right.edge means using inflection to end
    if (length(new.keep) >= 4) {
	    if (length(unique(x[new.keep])) < df) {
		    df <- length(unique(x[new.keep]))
	    }

        fitted.vals <- rep(NA_real_, length(keep))
        fit <- smooth.spline(x[new.keep], y[new.keep], df=df)
        fitted.vals[keep][new.keep] <- 10^fitted(fit)

        d1 <- predict(fit, deriv=1)$y
        d2 <- predict(fit, deriv=2)$y
        curvature <- d2/(1 + d1^2)^1.5 #fit new curve from inflection to end

        end_cliff <- 10^(y[new.keep][which.max(curvature)])
        end_cliff_rank <- 10^(x[new.keep][which.max(curvature)])
    } else {
        end_cliff <- 10^(y[new.keep][length(new.keep)])
        end_cliff_rank <- 10^(x[new.keep][length(new.keep)])
    }
    
    res <- list(end_cliff, end_cliff_rank, plateau)
    names(res) <- c("end_cliff", "end_cliff_rank", "plateau")

    return(res)
}

find_curve_bounds <- function(x, y, exclude.from) { #taken from EmptyDrops
    d1n <- diff(y)/diff(x)

    skip <- min(length(d1n) - 1, sum(x <= log10(exclude.from)))
    d1n <- tail(d1n, length(d1n) - skip)

    right.edge <- which.min(d1n)
    left.edge <- which.max(d1n[seq_len(right.edge)])

    c(left=left.edge, right=right.edge) + skip
}

set.seed(1234)
sce <- read10xCounts(inputDir)
br.out <- barcodeRanks(sce, lower = lowerForKnee)
tmp <- as.data.frame(metadata(br.out))

ranks <- knee_inflection_rank(sce, lower = lowerForKnee) #get ranks of the knee and inflection points
ranks <- as.data.frame(ranks)
tmp <- cbind(tmp, ranks)

metrics <- read.csv(cbMetrics, header = F)
colnames(metrics) <- c("metrics", "stats")
lowerForEndCliff <- (metrics[metrics$metrics == "found_cells", "stats"] + metrics[metrics$metrics == "found_empties", "stats"])

end_cliff <- endCliff(sce, inflection_rank = tmp$inflection_rank, lowerForKnee=lowerForKnee, lowerRankForEndCliff=lowerForEndCliff)
end_cliff <- as.data.frame(end_cliff)
tmp <- cbind(tmp, end_cliff)
write.table(tmp, outKnee, row.names = F, sep = "\t", quote = F)

e.out <- emptyDrops(sce, lower = lowerForKnee)
tmp <- as.data.frame(e.out)
tmp$barcode <- sce$Barcode
tmp <- tmp[!is.na(tmp$FDR) & tmp$FDR <= 0.005,] #maybe we should provide an argument to choose FDR threshold? In practice, I plot the adj. p-val distribution to choose a FDR threshold
write.table(tmp, outPass, col.names = T, row.names = F, sep = "\t", quote = F)

