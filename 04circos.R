library(RCircos)
rm(list=ls())

data(UCSC.HG19.Human.CytoBandIdeogram)

cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
chr.exclude <- NULL
tracks.inside <- 10
tracks.outside <- 1

RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside)  

RCircos.List.Plot.Parameters()

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()


data(RCircos.Gene.Label.Data)
data=RCircos.Gene.Label.Data
gene=read.table("16inter.txt",header=F,sep="\t",comment.char="",check.names=F)
intergene=intersect(data$Gene,gene$V1)

library(dplyr)
data <- data %>% filter(Gene %in% intergene)

side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(data, track.num, side)
name.col <- 4
track.num <- 2
RCircos.Gene.Name.Plot(data, name.col,track.num, side) 

data(RCircos.Heatmap.Data)

data.col <- 10
track.num <- 5
side <- "in"

RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col, track.num, side) 

data(RCircos.Scatter.Data) 

data.col <- 5 
track.num <- 6 
side <- "in"
by.fold <- 1 

RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col,track.num, side, by.fold)


data(RCircos.Line.Data) 
data.col <- 5 
track.num <- 7 
side <- "in" 

RCircos.Line.Data$chromosome = paste0("chr",RCircos.Line.Data$chromosome)

RCircos.Line.Plot(RCircos.Line.Data, data.col, track.num, side) 


data(RCircos.Histogram.Data) 
data.col <- 4 
track.num <- 8 
side <- "in"
RCircos.Histogram.Plot(RCircos.Histogram.Data, data.col, track.num, side)

data(RCircos.Tile.Data)
track.num <- 9
side <- "in"
RCircos.Tile.Plot(RCircos.Tile.Data, track.num, side)


data(RCircos.Link.Data) 
track.num <- 11
RCircos.Link.Plot(RCircos.Link.Data, track.num, TRUE)