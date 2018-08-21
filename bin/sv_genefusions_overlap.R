#!/usr/bin/Rscript
# =========================================================
#
# Copyright (C) 2016-2018, Nuno A. Fonseca  (nuno dot fonseca at gmail dot com)
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# if not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================
argv <- commandArgs(trailingOnly = FALSE)
base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
base_dir <- paste(base_dir,"/../",sep="")
# Version 2
#####################################################
args <- commandArgs(trailingOnly=TRUE)

if ( length(args) !=5 ) {
  cat("ERROR: usage: sv_genefusions_overlap.R sv_bedpe genefusion.sum.tsv tolerance_distance force_orientation@{y,n} out_prefix\n")
  q(status=1)
}
#args
sv.bedpe.file <- args[1]
genefusions.sum.file <- args[2]
tolerance.distance.i <- args[3]
force.orientation.i <- args[4]
out.prefix <- args[5]

# Assumptions: number of fusions < 20k
if ( !file.exists(sv.bedpe.file) ) {
  cat("ERROR: file ",sv.bedpe.file," not found.\n")
  q(status=1)
}

if ( !file.exists(genefusions.sum.file) ) {
  cat("ERROR: file ",genefusions.sum.file," not found.\n")
  q(status=1)
}

tolerance.distance <- as.numeric(tolerance.distance.i)

if ( is.null(tolerance.distance) || tolerance.distance < 0 ) {
  cat("Invalid tolerance distance value:", tolerance.distance.i,"\n")
  q(status=1)
}

if ( force.orientation.i!="y" &&  force.orientation.i!="n" ) {
  cat("Invalid force orientation value:", force.orientation.i,"\n")
  q(status=1)
}

force.orientation <- FALSE
if (force.orientation.i=="y" ) { force.orientation <-  TRUE }

# Ready to go...
cat("SV:",sv.bedpe.file,"\n")
cat("Gene Fusions:",genefusions.sum.file,"\n")
cat("Tolerance:",tolerance.distance," bases\n")
cat("Same orientation:",force.orientation,"\n")
cat("Output prefix:",out.prefix,"\n")

#########################################################################
cat("Loading SV...\n")
sv.df <- read.table(sv.bedpe.file,sep="\t",quote="",check.names=FALSE, header=TRUE,as.is=T)

expected.cols <- c("chrom1","start1","end1","chrom2","start2","end2","sv_id","pe_support","strand1","strand2","svclass","svmethod")
cols.not.found <- expected.cols[!expected.cols %in% colnames(sv.df)]
if ( length(cols.not.found)>0 ) {
  cat("ERROR: missing cols (", cols.not.found,") in ",sv.bedpe.file,"\n")
  q(status=1)
}
cat("SVs found:",nrow(sv.df),"\n")
cat("SVs loaded\n")

cat("Loading gene fusions...\n")
genefusions.df <- read.table(genefusions.sum.file,sep="\t",quote="",check.names=FALSE, header=TRUE,as.is=TRUE)

# 
# The header names differ :(
# an extra step to make them uniform
# FusionCatcher: Fusion_gene_name        Gene_name1      Gene_name2      Gene1_ensembl_ID        Gene2_ensembl_ID        Strand  Chromosome1     Breakpoint1     Chromosome2     Breakpoint2     Frameshift      Fusion_junction_sequence        Splicing_pattern      Number_of_supporting_reads
                                        # FusionMap: FusionGene      KnownGene1      KnownGene2      GeneId1 GeneId2 Strand  Chromosome1     Breakpoint1     Chromosome2     Breakpoint2     FrameShift      FusionJunctionSequence  SplicePattern   Number of Supporting Reads      nfusion
expected.cols <- c("FusionGene","KnownGene1","KnownGene2","GeneId1","GeneId2","Strand","Chromosome1","Breakpoint1","Chromosome2","Breakpoint2","FrameShift","FusionJunctionSequence","SplicePattern","Number of Supporting Reads")
if ( "Fusion_gene_name" %in% colnames(genefusions.df) ) {
# rename the columns to have the same names as described in the PCAWG file format
    # no extra tests...be brave :)
    colnames(genefusions.df) <- expected.cols
    cat("FusionCatcher based file\n")
} else {
    cat("FusionMap based file\n")
}
cols.not.found <- expected.cols[!expected.cols %in% colnames(genefusions.df)]

if ( length(cols.not.found)>0 ) {
  cat("ERROR: missing cols (", cols.not.found,") in ",genefusions.sum.file,"\n")
  q(status=1)
}
cat("Gene fusions found:",nrow(genefusions.df),"\n")
cat("Gene fusions loaded\n")

################################
# return the minimum distance
match.dist <- function(sv.pos,sv.chr,fusion.pos,fusion.chr,tolerance) {

    #stopifnot(!is.na(sv.pos))
    #stopifnot(!is.na(fusion.pos))
    if (sv.chr!=fusion.chr) return(NULL)
    
    sv.pos<-as.numeric(sv.pos)

    res <- list(dist=9999999999999,pos=NA,chr=fusion.chr)
    # take the closest Breakpointpoint
    x <- strsplit(as.character(fusion.pos),split=":")
    #cat("match.dis:",sv.pos,":",fusion.pos,"\n")
    for ( xx in x[[1]]) {
        d <- abs(as.numeric(xx)-sv.pos)
        if ( d < res$dist ) {
            res$dist <- d
            res$pos=xx
        }
    }
    if ( is.na(res$pos) || res$dist>tolerance ) 
        return(NULL)
    #print(res)
    return(res)      
}

#
#
match.fusion2sv.dist <- function(fusion,sv,tolerance,same.order=FALSE) {
    
    m <- list()
    m$reversed <- FALSE
    m$m1 <- NULL
    m$m2 <- NULL
    m$dist1 <- NULL
    m$dist2 <- NULL
    m$dist <- NULL
    m$sv <- sv
    m$fusion <- fusion
    m$matchtype <- ""

    
    m$m1 <- match.dist(sv["start1"],sv["chrom1"],fusion["Breakpoint1"],fusion["Chromosome1"],tolerance)
    #
    if ( ! is.null(m$m1) )  {
        m$m2 <- match.dist(sv["start2"],sv["chrom2"],fusion["Breakpoint2"],fusion["Chromosome2"],tolerance)
        if ( is.null(m$m2) ) {
            m$m1 <- NULL;# retry again
        } else {
            m$dist1 <- m$m1$dist
        }
    }

    if ( is.null(m$m1) && ! same.order ) {
        # swap
        m$m1 <- match.dist(sv["start2"],sv["chrom2"],fusion["Breakpoint1"],fusion["Chromosome2"],tolerance)
        if ( is.null(m$m1) )  { return(NULL) }
        m$reversed <- TRUE
        m$dist1 <- m$m1$dist
        m$m2 <- match.dist(sv["start1"],sv["chrom1"],fusion["Breakpoint2"],fusion["Chromosome2"],tolerance)
    }
    if ( ! is.null(m$m2) ) {
        m$dist2 <- m$m2$dist
        # keep the closest breakpoint
        m$Breakpoint1 <- m$m1$pos
        m$Breakpoint2 <- m$m2$pos
        m$m1 <- NULL
        m$m2 <- NULL
        return(m)
    }
    return(NULL)
}

direct.match.fusion2sv <- function(gene.fusion,svs.df,tolerance,same.order=FALSE) {
    # SVs in the same chr has the fusion
    c11 <- svs.df$chrom1==as.character(gene.fusion["Chromosome1"])
    c22 <- svs.df$chrom2==as.character(gene.fusion["Chromosome2"])
    c12 <- svs.df$chrom1==as.character(gene.fusion["Chromosome2"])
    c21 <- svs.df$chrom2==as.character(gene.fusion["Chromosome1"])
    # same chr
    chr.sel <-(c11&c22) | (c12&c21)
    svs.sel <- svs.df[chr.sel,,drop=FALSE]

    # check each SV against the fusion
    best.match <- NULL
    if (nrow(svs.sel) > 0 ) {
        for ( sv.idx in seq(1,nrow(svs.sel)) ) {            
            match <- match.fusion2sv.dist(gene.fusion,svs.sel[sv.idx,],tolerance,same.order)
            if ( !is.null(match) ) {
                if ( is.null(best.match) ) {
                best.match <- match
                } else { # minimum distance
                    if ( best.match$dist1+best.match$dist2 > best.match$dist1+best.match$dist2 ) {
                        best.match <- match
                    }
                }
            }
        }
        
    }
    if (!is.null(best.match) ) {
        cat("Match!\n")
        #print(best.match)
        # set the gusion breakpoint
        gene.fusion["Breakpoint1"] <- best.match$Breakpoint1
        gene.fusion["Breakpoint2"] <- best.match$Breakpoint2
        v1 <- c(best.match$dist1,best.match$dist2,NA,"simple")
        #print(best.match)
        names(v1) <- c("SV_dist1","SV_dist2","SV_sv.sv.dist","SV_SVmatchtype")
        msv <- best.match$sv
        colnames(msv) <- paste("SV_",colnames(msv),sep="")
        sv.vals <- append(msv,append(rep(NA,ncol(svs.df)),v1))
        gene.fusion <- append(gene.fusion,sv.vals)
    } else {
        new.vals <- rep(NA,(ncol(svs.df)*2+4))
        names(new.vals) <- append(append(paste("SV_",colnames(svs.df),sep=""),paste("SV2_",colnames(svs.df),sep="")),c("SV_dist1","SV_dist2","SV_sv.sv.dist","SV_SVmatchtype"))
        gene.fusion <- append(gene.fusion,new.vals)
    }
    return(gene.fusion)
}
#   WIP
#    Fusion 
#    b1--b2
#    |   |
#SV1 b3  |              b4|
#SV2     b6             b5|
indirect.match.fusion2sv <- function(gene.fusion,svs.df,tolerance) {

    has.nearby.match <- function(svs.nearby.entry) {
        if (is.null(svs.nearby.entry$m11) &&
            is.null(svs.nearby.entry$m12) &&
            is.null(svs.nearby.entry$m21) &&
            is.null(svs.nearby.entry$m22)) {
            return(FALSE);
        }
        return(TRUE)
    }
    # for each nearby SV in one brk
    # for each nearby SV in the other brk
    # check if the two other SVs breakpoints are nearby
#tolerance <- 500000
    SV.nearby.fusions <- function(gene.fusion,svs.sel,tolerance) {
                                        # svs in the same
        nearby <- list()
        if (nrow(svs.sel) > 0 ) {
            for ( sv.idx in seq(1,nrow(svs.sel)) ) {
                # match b1->b1
                m11 <- match.dist(svs.sel[sv.idx,"start1"],svs.sel[sv.idx,"chrom1"],gene.fusion["Breakpoint1"],gene.fusion["Chromosome1"],tolerance)
                m21 <- match.dist(svs.sel[sv.idx,"start2"],svs.sel[sv.idx,"chrom2"],gene.fusion["Breakpoint1"],gene.fusion["Chromosome1"],tolerance)
                m12 <- match.dist(svs.sel[sv.idx,"start1"],svs.sel[sv.idx,"chrom1"],gene.fusion["Breakpoint2"],gene.fusion["Chromosome2"],tolerance)
                m22 <- match.dist(svs.sel[sv.idx,"start2"],svs.sel[sv.idx,"chrom2"],gene.fusion["Breakpoint2"],gene.fusion["Chromosome2"],tolerance)
                #m{FB}{SV_B}
                nearby[[sv.idx]] <- list(
                    m11=m11,
                    m12=m12,
                    m21=m21,
                    m22=m22
                )
            }
        }
        names(nearby) <- svs.sel$sv_id
        return(nearby)
    }

    # SVs in the same chr has the fusion
    c11 <- svs.df$chrom1==as.character(gene.fusion["Chromosome1"])
    c12 <- svs.df$chrom1==as.character(gene.fusion["Chromosome2"])
    c22 <- svs.df$chrom2==as.character(gene.fusion["Chromosome2"])
    c21 <- svs.df$chrom2==as.character(gene.fusion["Chromosome1"])

    #tolerance <- 1000000
    # same chr
    chr.sel <- (c11|c12)|(c22|c21)
    svs.sel <- svs.df[chr.sel,,drop=FALSE]
    rownames(svs.sel) <- svs.sel$sv_id
    svs.sel
    #tolerance
    # returns a list with all SVs nearby the SV breapoints
    svs.nearby <- SV.nearby.fusions(gene.fusion,svs.sel,tolerance)

    best.match <- NULL
    idxs <- seq(1,nrow(svs.sel))
    ##cat("Fusion:",gene.fusion$FusionGene,"\n")
    if (nrow(svs.sel)>=2) {
        for ( idx1 in idxs ) {
            sv1 <- names(svs.nearby)[idx1]
            svs.nearby1 <- svs.nearby[[sv1]]
            if (idx1==nrow(svs.sel)) break;
            if (!has.nearby.match(svs.nearby[[sv1]])) next;
            # 
            #cat("Matched SV1:",sv1,"\n")
            for ( idx2 in seq(idx1+1,nrow(svs.sel) )) {
                sv2 <- names(svs.nearby)[idx2]
                if (!has.nearby.match(svs.nearby[[sv2]])) next;
                svs.nearby2 <- svs.nearby[[sv2]]
                #cat("Matched:",sv1,"--",sv2,"\n")
                fusion.case <- NULL
                svm <- NULL
                #        Fusion               SV1           SV2
                #     b1        b2       sv1b1     sv1b2  sv2b1 sv2b2
                #                          a          b      c     d
                # A1  a          c                    d            b 
                # B1  c          a                    d            b 
                # A2  
                # B2
                # A3
                # B3
                # A4
                # B4
            # A) SV1.bp->Fusion.bp1 -> SV2.bp->Fusion.bp2
            # B) SV1.bp->Fusion.bp2 -> SV2.bp->Fusion.bp1
            # SV2.bp->Fusion.bp2
            # Fbp1--SV1bp1|SV1bp2 => Fbp2--SV2bp1|SV21bp2
            # case:
            # 1: SV1.b1 -> SV2.b2
            # 2: SV1.b2 -> SV2.b1
            # 3: SV1.b1 -> SV2.b1
            # 4: SV1.b2 -> SV2.b2
                # A 1
                if ( !is.null(svs.nearby[[sv1]]$m11) &&
                     !is.null(svs.nearby[[sv2]]$m22) &&
                     svs.sel[sv1,"chrom2"]==svs.sel[sv2,"chrom1"] ) {
                                        # the different cases
                    sv1.chrm <- svs.sel[sv1,"chrom2"]
                    sv1.pos <- svs.sel[sv1,"start2"]
                    sv2.chrm <- svs.sel[sv2,"chrom1"]
                    sv2.pos <- svs.sel[sv2,"start1"]
                    dist1 <- svs.nearby[[sv1]]$m11$dist
                    dist2 <- svs.nearby[[sv2]]$m22$dist
                    pos1 <- svs.nearby[[sv1]]$m11$pos                
                    pos2 <- svs.nearby[[sv2]]$m22$pos
                    sv.id1 <- sv1
                    sv.id2 <- sv2
                    fusion.case <- "A"
                    sv.case <- 1
                }
                # B 1
                if ( !is.null(svs.nearby[[sv2]]$m11) &&
                     !is.null(svs.nearby[[sv1]]$m22) &&
                     svs.sel[sv1,"chrom1"]==svs.sel[sv2,"chrom2"] ) {
                                        # the different cases
                    sv1.chrm <- svs.sel[sv1,"chrom1"]
                    sv1.pos <- svs.sel[sv1,"start1"]
                    sv2.chrm <- svs.sel[sv2,"chrom2"]
                    sv2.pos <- svs.sel[sv2,"start2"]
                    dist1 <- svs.nearby[[sv2]]$m11$dist
                    dist2 <- svs.nearby[[sv1]]$m22$dist
                    pos1 <- svs.nearby[[sv2]]$m11$pos                
                    pos2 <- svs.nearby[[sv1]]$m22$pos
                    sv.id1 <- sv2
                    sv.id2 <- sv1
                    fusion.case <- "B"
                    sv.case <- 1
                }
                # A 2
                if ( !is.null(svs.nearby[[sv1]]$m21) &&
                     !is.null(svs.nearby[[sv2]]$m12) &&
                     svs.sel[sv1,"chrom1"]==svs.sel[sv2,"chrom2"] ) {
                    # the different cases
                    sv1.chrm <- svs.sel[sv1,"chrom1"]
                    sv1.pos <- svs.sel[sv1,"start1"]
                    sv2.chrm <- svs.sel[sv2,"chrom2"]
                    sv2.pos <- svs.sel[sv2,"start2"]
                    dist1 <- svs.nearby[[sv1]]$m21$dist
                    dist2 <- svs.nearby[[sv2]]$m12$dist
                    pos1 <- svs.nearby[[sv1]]$m21$pos                
                    pos2 <- svs.nearby[[sv2]]$m12$pos
                    sv.id1 <- sv1
                    sv.id2 <- sv2
                    fusion.case <- "A"
                    sv.case <- 2
                }
                # B 2
                if ( !is.null(svs.nearby[[sv2]]$m21) &&
                     !is.null(svs.nearby[[sv1]]$m12) &&
                     svs.sel[sv1,"chrom2"]==svs.sel[sv2,"chrom1"] ) {
                    sv1.chrm <- svs.sel[sv1,"chrom2"]
                    sv1.pos <- svs.sel[sv1,"start2"]
                    sv2.chrm <- svs.sel[sv2,"chrom1"]
                    sv2.pos <- svs.sel[sv2,"start1"]
                    dist1 <- svs.nearby[[sv2]]$m21$dist
                    dist2 <- svs.nearby[[sv1]]$m12$dist
                    pos1 <- svs.nearby[[sv2]]$m21$pos                
                    pos2 <- svs.nearby[[sv1]]$m12$pos
                    sv.id1 <- sv2
                    sv.id2 <- sv1
                    fusion.case <- "B"
                    sv.case <- 2
                }

                # A 3
                if ( !is.null(svs.nearby[[sv1]]$m11) &&
                     !is.null(svs.nearby[[sv2]]$m12) &&
                     svs.sel[sv1,"chrom2"]==svs.sel[sv2,"chrom2"] ) {                                    
                    sv1.chrm <- svs.sel[sv1,"chrom2"]
                    sv1.pos <- svs.sel[sv1,"start2"]
                    sv2.chrm <- svs.sel[sv2,"chrom2"]
                    sv2.pos <- svs.sel[sv2,"start2"]
                    dist1 <- svs.nearby[[sv1]]$m11$dist
                    dist2 <- svs.nearby[[sv2]]$m12$dist
                    pos1 <- svs.nearby[[sv1]]$m11$pos                
                    pos2 <- svs.nearby[[sv2]]$m12$pos
                    sv.id1 <- sv1
                    sv.id2 <- sv2
                    fusion.case <- "A"
                    sv.case <- 3
                }
                # B 3
                if ( !is.null(svs.nearby[[sv2]]$m11) &&
                 !is.null(svs.nearby[[sv1]]$m12) &&
                 svs.sel[sv1,"chrom2"]==svs.sel[sv2,"chrom2"] ) {
                    sv1.chrm <- svs.sel[sv1,"chrom2"]
                    sv1.pos <- svs.sel[sv1,"start2"]
                    sv2.chrm <- svs.sel[sv2,"chrom2"]
                    sv2.pos <- svs.sel[sv2,"start2"]
                    dist1 <- svs.nearby[[sv2]]$m11$dist
                    dist2 <- svs.nearby[[sv1]]$m12$dist
                    pos1 <- svs.nearby[[sv2]]$m11$pos                
                    pos2 <- svs.nearby[[sv1]]$m12$pos
                    sv.id1 <- sv2
                    sv.id2 <- sv1
                    fusion.case <- "B"
                    sv.case <- 3

                }
                # A 4
                if ( !is.null(svs.nearby[[sv1]]$m21) &&
                     !is.null(svs.nearby[[sv2]]$m22) &&
                     svs.sel[sv1,"chrom1"]==svs.sel[sv2,"chrom1"] ) {                                    
                    sv1.chrm <- svs.sel[sv1,"chrom1"]
                    sv1.pos <- svs.sel[sv1,"start1"]
                    sv2.chrm <- svs.sel[sv2,"chrom1"]
                    sv2.pos <- svs.sel[sv2,"start1"]
                    dist1 <- svs.nearby[[sv1]]$m21$dist
                    dist2 <- svs.nearby[[sv2]]$m22$dist
                    pos1 <- svs.nearby[[sv1]]$m21$pos                
                    pos2 <- svs.nearby[[sv2]]$m22$pos
                    sv.id1 <- sv1
                    sv.id2 <- sv2
                    fusion.case <- "A"
                    sv.case <- 4

                }
                # B 4
                if ( !is.null(svs.nearby[[sv2]]$m21) &&
                     !is.null(svs.nearby[[sv1]]$m22) &&
                     svs.sel[sv1,"chrom1"]==svs.sel[sv2,"chrom1"] ) {
                    sv1.chrm <- svs.sel[sv1,"chrom1"]
                    sv1.pos <- svs.sel[sv1,"start1"]
                    sv2.chrm <- svs.sel[sv2,"chrom1"]
                    sv2.pos <- svs.sel[sv2,"start1"]
                    dist1 <- svs.nearby[[sv2]]$m21$dist
                    dist2 <- svs.nearby[[sv1]]$m22$dist
                    pos1 <- svs.nearby[[sv2]]$m21$pos                
                    pos2 <- svs.nearby[[sv1]]$m22$pos
                    sv.id1 <- sv2
                    sv.id2 <- sv1
                    fusion.case <- "B"
                    sv.case <- 4

                }
                
                if (!is.null(fusion.case)) {
                    svm <- match.dist(sv1.pos,sv1.chrm,sv2.pos,sv2.chrm,tolerance)
                }
                if ( !is.null(svm) ) {
                    # check if he SVs breakpoints are also nearby
                    # check SV type?
                    # check strand?
                    cat("Complex match")
                    cat(fusion.case,sv.case,"\n")  
                    mdist <- svm$dist+dist1+dist2
                    new.best.match <- list(
                        sv1=svs.sel[sv1,],
                        sv2=svs.sel[sv2,],
                        tdist=mdist,
                        sv_sv_dist=svm$dist,
                        sv_dist=dist1,
                        sv_dist2=dist2,
                        pos1=pos1,
                        pos2=pos2,
                        SVmatchtype=paste0("complex:",fusion.case,":",sv.case,":",svs.sel[sv1,"svclass"],":",svs.sel[sv2,"svclass"])
                    )
                    #print(best.match)
                    #print("------------------------------------------------")
                    #print(new.best.match)
                    if ( is.null(best.match) ) {
                        best.match <- new.best.match
                    } else {
                        if ( mdist<best.match$tdist ) {
                            cat("Alternative best complex match picked:",best.match$tdist,"->",mdist,"\n")                        
                            best.match <- new.best.match
                        }
                    }
                }
            }            
        }
    }

    # add a sv2_id info
    #"chrom1"     "start1"     "end1"       "chrom2"     "start2"    
    #"end2"       "sv_id"      "pe_support" "strand1"    "strand2"   
    # "svclass"    "svmethod"
    if (!is.null(best.match)  ) {
        cat("Complex Match!\n")
        #print(best.match)
        # set the gusion breakpoint
        gene.fusion["Breakpoint1"] <- best.match$pos1
        gene.fusion["Breakpoint2"] <- best.match$pos2
        v1 <- c(best.match$sv_dist,best.match$sv_dist2,best.match$sv_sv_dist,best.match$SVmatchtype)

        names(v1) <- c("SV_dist1","SV_dist2","SV_sv.sv.dist","SV_SVmatchtype")
        colnames(best.match$sv2) <- paste("SV2_",colnames(best.match$sv2),sep="")
        colnames(best.match$sv1) <- paste("SV_",colnames(best.match$sv1),sep="")
        sv.vals <- append(best.match$sv1,best.match$sv2)
        gene.fusion <- append(gene.fusion,append(sv.vals,v1))
    } else {
        new.vals <- rep(NA,(ncol(svs.df)*2+4))
        names(new.vals) <- append(append(paste("SV_",colnames(svs.df),sep=""),paste("SV2_",colnames(svs.df),sep="")),c("SV_dist1","SV_dist2","SV_sv.sv.dist","SV_SVmatchtype"))
        gene.fusion <- append(gene.fusion,new.vals)
    }

    return(gene.fusion)
}


match.fusion2sv <- function(genefusion,svs.df,tolerance=0,same.order=FALSE) {
    stopifnot ( !is.vector(genefusion) )
    m <- direct.match.fusion2sv(gene.fusion=genefusion,svs.df,tolerance,same.order)
    if ( is.na(m$SV_sv_id) ) {
        # indirect match        
        m1 <- indirect.match.fusion2sv(genefusion,svs.df,tolerance)
        if (!is.na(m1$SV_sv_id) ) {
            m <- m1
        }
    }
    return(m)
}

# normalize the column values (some chr. may include the prefix chr)
if (nrow(genefusions.df) > 0 ) {
    genefusions.df$Chromosome1 <- gsub("chr","",as.character(genefusions.df$Chromosome1),ignore.case = TRUE)
    genefusions.df$Chromosome2 <- gsub("chr","",as.character(genefusions.df$Chromosome2),ignore.case = TRUE)
}

######################################################
#idx <- 13
#sv.df
#tolerance.distance <- 2500000
#tolerance.distance <- 25000
res <- data.frame(matrix(nrow=nrow(genefusions.df),ncol=ncol(sv.df)*2+ncol(genefusions.df)+4))
colnames(res) <- append(append(colnames(genefusions.df),append(paste("SV_",colnames(sv.df),sep=""),paste("SV2_",colnames(sv.df),sep=""))),c("SV_dist1","SV_dist2","SV_sv.sv.dist","SV_SVmatchtype"))
#genefusions.df$nfusion
for ( idx in seq(1,nrow(genefusions.df)) ) {   
  fusionid <- as.character(genefusions.df[idx,"FusionGene"])
  cat("GeneFusion:",fusionid,"\n")
  res1 <- match.fusion2sv(genefusion=genefusions.df[idx,],svs.df=sv.df,tolerance=tolerance.distance,same.order=force.orientation)
  stopifnot(length(res1)==ncol(res))
  res[idx,] <- res1
}

## rename the types
x <- res[,"SV_SVmatchtype"]
x[is.na(x)] <- 0
x[x=="simple"] <- 1
x[grepl("complex",x)] <- 2
x[x==0] <- "none"
x[x==1] <- "direct"
x[x==2] <- "composite" ;# bridged

## TODO
res$SV_support_type <- x
# Save the results
ofile <- paste0(out.prefix,"_matched_sv.tsv")
write.table(res,file=ofile,
            sep="\t",row.names=FALSE,col.names=TRUE,
            quote=FALSE)
cat("Created ",ofile," ",sum(!is.na(res$SV_sv_id)),"\n")
q(status=0)







