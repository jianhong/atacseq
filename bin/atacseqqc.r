#!/usr/bin/env Rscript

library(ATACseqQC)
library(ChIPpeakAnno)
library(MotifDb)
library(GenomicAlignments)
library(optparse)
library(GenomicFeatures)
library(rtracklayer)
library(Rsamtools)

option_list <- list(make_option(c("-n", "--name"), type="character", default=NULL, help="prefix name of output files", metavar="path"),
                    make_option(c("-p", "--peaks"), type="character", default=NULL, help="peak files", metavar="string"),
                    make_option(c("-b", "--bams"), type="character", default=NULL, help="bam files", metavar="string"),
                    make_option(c("-g", "--gtf"), type="character", default=NULL, help="filename of gtf file", metavar="path"),
                    make_option(c("-c", "--cores"), type="integer", default=1, help="Number of cores", metavar="integer"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$bams)){
  print_help(opt_parser)
  stop("Please provide bam file name", call.=FALSE)
}
if (is.null(opt$peaks)){
  print_help(opt_parser)
  stop("Please provide peak file name.", call.=FALSE)
}
if (is.null(opt$gtf)){
  print_help(opt_parser)
  stop("Please provide gtf file.", call.=FALSE)
}

txdb <- makeTxDbFromGFF(opt$gtf)
gtf <- import(opt$gtf)

bamfile.labels <- opt$name
if(is.null(bamfile.labels)) bamfile.labels <- gsub(".bam", "", basename(bamfile))

pf <- make.names(opt$name)
dir.create(pf)
## generate fragement size distribution
pdf(file.path(pf, paste0("fragmentSizeDistribution.", bamfile.labels, ".pdf")))
fragSize <- fragSizeDist(bamfile, bamfile.labels)
dev.off()
png(file.path(pf, paste0("fragmentSizeDistribution.", bamfile.labels, ".png")))
fragSize <- fragSizeDist(bamfile, bamfile.labels)
dev.off()

possibleTag <- combn(LETTERS, 2)
possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]),
                 paste0(possibleTag[2, ], possibleTag[1, ]))

bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                     param = ScanBamParam(tag=possibleTag, what = scanBamWhat()))[[1]]
bamTop100tag <- bamTop100$tag
tags <- names(bamTop100tag)[lengths(bamTop100tag)==100]
tags <- tags[tags!="PG"]

outPath <- file.path(pf, "splited")
dir.create(outPath)

which <- as(seqinfo(txdb), "GRanges")
seqlevelsStyle <- seqlevelsStyle(levels(bamTop100$rname))[1]
gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
shiftedBamfile <- file.path(outPath, "shifted.bam")
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
txs <- transcripts(txdb)
pt <- PTscore(gal1, txs)
pdf(file.path(pf, paste0("PT_score.", bamfile.labels, ".pdf")))
plot(pt$log2meanCoverage, pt$PT_score, 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")
dev.off()
png(file.path(pf, paste0("PT_score.", bamfile.labels, ".png")))
plot(pt$log2meanCoverage, pt$PT_score, 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")
dev.off()

nfr <- NFRscore(gal1, txs)
pdf(file.path(pf, paste0("NFRscore.", bamfile.labels, ".pdf")))
plot(nfr$log2meanCoverage, nfr$NFR_score, 
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5))
dev.off()
png(file.path(pf, paste0("NFRscore.", bamfile.labels, ".png")))
plot(nfr$log2meanCoverage, nfr$NFR_score, 
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5))
dev.off()

tsse <- TSSEscore(gal1, txs)
pdf(file.path(pf, paste0("TSSEscore.", bamfile.labels, ".pdf")))
plot(100*(-9:10-.5), tsse$values, type="b", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")
dev.off()
png(file.path(pf, paste0("TSSEscore.", bamfile.labels, ".png")))
plot(100*(-9:10-.5), tsse$values, type="b", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")
dev.off()

## split the reads into NucleosomeFree, mononucleosome, 
## dinucleosome and trinucleosome.
## and save the binned alignments into bam files.
objs <- splitGAlignmentsByCut(gal1, txs=txs, outPath = outPath)

bamfiles <- file.path(outPath,
                      c("NucleosomeFree.bam",
                        "mononucleosome.bam",
                        "dinucleosome.bam",
                        "trinucleosome.bam"))

TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)
(librarySize <- estLibSize(bamfiles))

NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", 
                                     "mononucleosome",
                                     "dinucleosome",
                                     "trinucleosome")], 
                          TSS=TSS,
                          librarySize=librarySize,
                          seqlev=seqlev,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)

## log2 transformed signals
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
#plot heatmap
pdf(file.path(pf, paste0("featureAlignedHeatmap.", bamfile.labels, ".pdf")))
featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE)
dev.off()

out <- featureAlignedDistribution(sigs, 
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l", 
                                  ylab="Averaged coverage")

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)
pdf(file.path(pf, paste0("featureAlignedTSScurve.", bamfile.labels, ".pdf")))
matplot(out, type="l", xaxt="n", 
        xlab="Position (bp)", 
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1, 
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()
png(file.path(pf, paste0("featureAlignedTSScurve.", bamfile.labels, ".png")))
matplot(out, type="l", xaxt="n", 
        xlab="Position (bp)", 
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1, 
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()


CTCF <- query(MotifDb, c("CTCF"))
CTCF <- as.list(CTCF)

pdf(file.path(pf, paste0("CTCF.footprint.", bamfile.labels, ".pdf")))
sigs <- factorFootprints(shiftedBamfile, pfm=CTCF[[1]], 
                         genome=genome,
                         min.score="90%", seqlev=seqlev,
                         upstream=100, downstream=100)
dev.off()
png(file.path(pf, paste0("CTCF.footprint.", bamfile.labels, ".png")))
sigs <- factorFootprints(shiftedBamfile, pfm=CTCF[[1]], 
                         genome=genome,
                         min.score="90%", seqlev=seqlev,
                         upstream=100, downstream=100)
dev.off()

pdf(file.path(pf, paste0("CTCF.vplot.", bamfile.labels, ".pdf")))
vp <- vPlot(shiftedBamfile, pfm=CTCF[[1]], 
            genome=genome, min.score="90%", seqlev=seqlev,
            upstream=200, downstream=200, 
            ylim=c(30, 250), bandwidth=c(2, 1))
dev.off()
pdf(file.path(pf, paste0("CTCF.vplot.distanceDyad.", bamfile.labels, ".pdf")))
distanceDyad(vp, pch=20, cex=.5)
dev.off()
png(file.path(pf, paste0("CTCF.vplot.distanceDyad.", bamfile.labels, ".png")))
distanceDyad(vp, pch=20, cex=.5)
dev.off()


