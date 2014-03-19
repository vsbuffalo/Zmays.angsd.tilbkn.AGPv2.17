# make an alias for zea mays transcriptDb
txdb <- TxDb.Zmays.Ensembl.AGPv2.17
txdb_filtered <- TxDb.Zmays.EnsemblFiltered.AGPv2.17

## Load ANGSD data for landraces, and append the TIL data (creating separate
## col names).
tmp_rdata_file <- "inst/extdata/angsd_tilbkn_pestpg.Rdata"
if (!file.exists(tmp_rdata_file)) {
	Zmays.angsd.tilbkn.AGPv2.17 <- readPestPG("inst/extdata/BKN.pestPG.gz")
	colnames(mcols(Zmays.angsd.tilbkn.AGPv2.17)) <- paste0(colnames(mcols(Zmays.angsd.tilbkn.AGPv2.17)), "_bkn")
	til <- readPestPG("inst/extdata/TIL.pestPG.gz", asGRanges=FALSE, asDataFrame=TRUE)
	colnames(til) <- paste0(colnames(til), "_til")
	mcols(Zmays.angsd.tilbkn.AGPv2.17) <- cbind(mcols(Zmays.angsd.tilbkn.AGPv2.17), til[, -c(1:4)])
	rm(til)

  # rename seqlevels for bkn
  new_names  <- local({
    nn <- names(seqlengths(txdb))
    names(nn)  <- nn
    names(nn)[names(nn) == "chloroplast"] <- "Pt"
    names(nn)[names(nn) == "mitochondrion"] <- "Mt"
    nn
  })

  # change seqlevels of bkn GRanges to match txdb (Mt=mitochondrion,
  # Pt=chloroplast), and add lengths.
  seqlevels(Zmays.angsd.tilbkn.AGPv2.17) <- unname(new_names[seqlevels(Zmays.angsd.tilbkn.AGPv2.17)])
  seqlengths(Zmays.angsd.tilbkn.AGPv2.17) <- seqlengths(txdb)[seqlevels(Zmays.angsd.tilbkn.AGPv2.17)]

  save(Zmays.angsd.tilbkn.AGPv2.17, file=tmp_rdata_file)
} else {
  load(tmp_rdata_file)
}

getNoncoding <- function(coding) {
	setdiff(as(seqinfo(coding), "GRanges"),
					reduce(coding, ignore.strand=TRUE))
}

coding <- genes(txdb)
noncoding <- getNoncoding(coding)
coding_filtered <- genes(txdb_filtered)
noncoding_filtered <- getNoncoding(coding_filtered)


Zmays.angsd.tilbkn.AGPv2.17$coding_wgs <- overlapWidths(Zmays.angsd.tilbkn.AGPv2.17, coding)
Zmays.angsd.tilbkn.AGPv2.17$noncoding_wgs <- overlapWidths(Zmays.angsd.tilbkn.AGPv2.17, noncoding)

Zmays.angsd.tilbkn.AGPv2.17$coding_fgs <- overlapWidths(Zmays.angsd.tilbkn.AGPv2.17, coding_filtered)
Zmays.angsd.tilbkn.AGPv2.17$noncoding_fgs <- overlapWidths(Zmays.angsd.tilbkn.AGPv2.17, noncoding_filtered)
rdata_file <- "inst/extdata/Zmays.angsd.tilbkn.AGPv2.17.Rdata"
save(Zmays.angsd.tilbkn.AGPv2.17, file=rdata_file)

