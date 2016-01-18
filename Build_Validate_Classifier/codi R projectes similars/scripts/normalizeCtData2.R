normalizeCtData2<-function (q, norm = "deltaCt", deltaCt.genes = NULL, scale.rank.samples, 
          rank.type = "pseudo.median", Ct.max = 35, geo.mean.ref, verbose = TRUE) 
{
  data <- exprs(q)
  data.norm <- data
  method <- match.arg(norm, c("quantile", "scale.rankinvariant", 
                              "norm.rankinvariant", "deltaCt", "geometric.mean"))
  if (method %in% c("scale.rankinvariant", "norm.rankinvariant")) {
    Ct.index <- data > Ct.max
    data.Ctmax <- data
    data.Ctmax[Ct.index] <- NA
    if (rank.type == "pseudo.median") {
      ref.data <- apply(data.Ctmax, 1, median, na.rm = TRUE)
    }
    else if (rank.type == "pseudo.mean") {
      ref.data <- apply(data.Ctmax, 1, mean, na.rm = TRUE)
    }
    na.index <- is.na(ref.data)
    ref.data[na.index] <- 30
    data.rankinvar <- apply(data, 2, normalize.invariantset, 
                            ref = ref.data)
  }
  switch(method, quantile = {
    data.norm <- normalizeQuantiles(data)
  }, scale.rankinvariant = {
    ri.genes <- sapply(data.rankinvar, "[[", "i.set")
    ri.genes[Ct.index] <- FALSE
    ri.genes[na.index, ] <- FALSE
    ri.count <- rowSums(ri.genes)
    if (missing(scale.rank.samples)) scale.rank.samples <- ncol(data) - 
      1
    ri.index <- ri.count >= scale.rank.samples
    if (sum(ri.index) == 0) stop(paste("No rank invariant genes were found across", 
                                       scale.rank.samples, "samples"))
    ri.mean <- colMeans(data[ri.index, , drop = FALSE])
    ri.scale <- ri.mean/ri.mean[1]
    data.norm <- t(t(data) * ri.scale)
    if (verbose) {
      cat(c("Scaling Ct values\n\tUsing rank invariant genes:", 
            paste(featureNames(q)[ri.index], collapse = " "), 
            "\n"))
      cat(c("\tScaling factors:", format(ri.scale, digits = 3), 
            "\n"))
    }
  }, norm.rankinvariant = {
    if (verbose) cat("Normalizing Ct values\n\tUsing rank invariant genes:\n")
    for (i in 1:ncol(data)) {
      ri.sub <- data.rankinvar[[i]]
      ri.genes <- ri.sub[["i.set"]]
      ri.genes[Ct.index[, i]] <- FALSE
      ri.genes[na.index] <- FALSE
      if (sum(ri.genes) == 0) {
        warning(paste("\tNo rank invariant genes were found for sample ", 
                      sampleNames(q)[i], "; sample not normalized\n", 
                      sep = ""))
        next
      }
      if (verbose) cat(paste("\t", sampleNames(q)[i], ": ", 
                             sum(ri.genes), " rank invariant genes\n", sep = ""))
      data.norm[, i] <- as.numeric(approx(ri.sub$n.curve$y, 
                                          ri.sub$n.curve$x, xout = data[, i], rule = 2)$y)
    }
  }, deltaCt = {
    if (is.null(deltaCt.genes)) deltaCt.genes <- unique(featureNames(q)[featureType(q) == 
                                                                          "Endogenous Control"])
    c.index <- featureNames(q) %in% deltaCt.genes
    if (verbose) {
      cat(c("Calculating deltaCt values\n\tUsing control gene(s):", 
            paste(deltaCt.genes, collapse = " "), "\n"))
    }
    for (c in 1:ncol(data)) {
      refCt <- mean(data[c.index, c], na.rm = TRUE)
      refsd <- sd(data[c.index, c], na.rm = TRUE)
      data.norm[, c] <- data[, c] - refCt
      if (verbose) cat(paste("\tCard ", c, ":\tMean=", 
                             format(refCt, dig = 4), "\tStdev=", format(refsd, 
                                                                        dig = 3), "\n", sep = ""))
    }
  }, geometric.mean = {
    geo.mean <- apply(data, 2, function(x) {
      xx <- log2(subset(x, x < Ct.max))
      2^mean(xx)
    })
    if (missing(geo.mean.ref)) geo.mean.ref <- 1
    geo.scale <- geo.mean/geo.mean[geo.mean.ref]
    data.norm <- t(t(data) / geo.scale)
    if (verbose) {
      cat(c("Scaling Ct values\n\tUsing geometric mean within each sample\n"))
      cat(c("\tScaling factors:", format(geo.scale, digits = 3), 
            "\n"))
    }
  })
  exprs(q) <- data.norm
  if (nrow(getCtHistory(q)) == 0) 
    setCtHistory(q) <- data.frame(history = "Manually created qPCRset object.", 
                                  stringsAsFactors = FALSE)
  setCtHistory(q) <- rbind(getCtHistory(q), capture.output(match.call(normalizeCtData)))
  q
}