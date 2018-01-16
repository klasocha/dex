load("dex-project/report/plot_data.RData") 

isoforms <- read.table(Sys.glob("dex-project/cufflinks/*/isoforms.fpkm_tracking")[1]) 
significant_isoforms = unlist(strsplit(labels(dend01_5cl),'_'))[2*(1:length(labels(dend01_5cl)))]
isoforms <- isoforms[match(significant_isoforms, isoforms$V1),]
transcript_strand <- read.table("dex-project/cufflinks/transcript-strand")

isoforms$V2 <- NULL
isoforms$V3 <- NULL
isoforms$V8 <- NULL
isoforms$V9 <- NULL
isoforms$V10 <- NULL
isoforms$V11 <- NULL
isoforms$V12 <- NULL
isoforms$V13 <- NULL

colnames(isoforms) <- c("tracking_id",	"gene_id",	"gene_short_name",	"tss_id",	"locus")


isoforms$strand <- transcript_strand$V1[match(isoforms$tracking_id, transcript_strand$V2)]

iso_chromosomes = c()
iso_starts = c()
iso_ends = c()

for (i in 1:nrow(isoforms)){
  iso_chromosomes = c(iso_chromosomes, strsplit(isoforms$locus[i], ":")[[1]][1])
  iso_starts = c(iso_starts, strsplit(strsplit(isoforms$locus[i], ":")[[1]][2], "-")[[1]][1])
  iso_ends =  c(iso_ends, strsplit(strsplit(isoforms$locus[i], ":")[[1]][2], "-")[[1]][2])
}

iso_starts =  as.numeric(iso_starts)
iso_ends =  as.numeric(iso_ends)

peaks1 = read.table(Sys.glob("dex-project/human-beds/*/hglft*.bed")[1]) 
peaks2 = read.table(Sys.glob("dex-project/human-beds/*/hglft*.bed")[2]) 
peaks3 = read.table(Sys.glob("dex-project/human-beds/*/hglft*.bed")[3]) 



for (i in 1:nrow(isoforms)){
  iso_start = iso_starts[i]
  if (isoforms$strand[i] == "-"){
    iso_start = iso_ends[i]
  }
  
  if (sum(peaks1$V1 == iso_chromosomes[i]) > 0){
    correct_chromosome = peaks1$V1 == iso_chromosomes[i]
    distances = abs((peaks1$V2[correct_chromosome] + peaks1$V3[correct_chromosome])/2 - iso_start)
    
    isoforms$peak1_chr[i] <- peaks1$V1[correct_chromosome][head(order(distances), n=1)]
    isoforms$peak1_start[i] <- peaks1$V2[correct_chromosome][head(order(distances), n=1)]
    isoforms$peak1_end[i] <- peaks1$V3[correct_chromosome][head(order(distances), n=1)]
    isoforms$peak1_dist[i] <- head(sort(distances), n=1)
  }
  if (sum(peaks2$V1 == iso_chromosomes[i]) > 0){
    correct_chromosome = peaks2$V1 == iso_chromosomes[i]
    distances = abs((peaks2$V2[correct_chromosome] + peaks2$V3[correct_chromosome])/2 - iso_start)
    
    isoforms$peak2_chr[i] <- peaks2$V1[correct_chromosome][head(order(distances), n=1)]
    isoforms$peak2_start[i] <- peaks2$V2[correct_chromosome][head(order(distances), n=1)]
    isoforms$peak2_end[i] <- peaks2$V3[correct_chromosome][head(order(distances), n=1)]
    isoforms$peak2_dist[i] <- head(sort(distances), n=1)
  }
  if (sum(peaks3$V1 == iso_chromosomes[i]) > 0){
    correct_chromosome = peaks3$V1 == iso_chromosomes[i]
    distances = abs((peaks3$V2[correct_chromosome] + peaks3$V3[correct_chromosome])/2 - iso_start)
    
    isoforms$peak3_chr[i] <- peaks3$V1[correct_chromosome][head(order(distances), n=1)]
    isoforms$peak3_start[i] <- peaks3$V2[correct_chromosome][head(order(distances), n=1)]
    isoforms$peak3_end[i] <- peaks3$V3[correct_chromosome][head(order(distances), n=1)]
    isoforms$peak3_dist[i] <- head(sort(distances), n=1)
  }
}


isoforms$cl5 <- as.numeric(cutree(hc01, 5)[labels(dend01_5cl)])
isoforms$cl3 <- as.numeric(cutree(hc01, 3)[labels(dend01_5cl)])
isoforms$cl2 <- as.numeric(cutree(hc01, 2)[labels(dend01_5cl)])


write.csv(isoforms, file = "isoforms_peak.csv",row.names=FALSE)


#Random stuff


random_isoforms <- read.table(Sys.glob("dex-project/cufflinks/*/isoforms.fpkm_tracking")[1]) 
random_sample <- sample(1:nrow(random_isoforms), nrow(isoforms), replace=F)
random_isoforms <- random_isoforms[random_sample,]

random_isoforms$V2 <- NULL
random_isoforms$V3 <- NULL
random_isoforms$V8 <- NULL
random_isoforms$V9 <- NULL
random_isoforms$V10 <- NULL
random_isoforms$V11 <- NULL
random_isoforms$V12 <- NULL
random_isoforms$V13 <- NULL

colnames(random_isoforms) <- c("tracking_id",	"gene_id",	"gene_short_name",	"tss_id",	"locus")


random_isoforms$strand <- transcript_strand$V1[match(random_isoforms$tracking_id, transcript_strand$V2)]

iso_chromosomes = c()
iso_starts = c()
iso_ends = c()

for (i in 1:nrow(random_isoforms)){
  iso_chromosomes = c(iso_chromosomes, strsplit(random_isoforms$locus[i], ":")[[1]][1])
  iso_starts = c(iso_starts, strsplit(strsplit(random_isoforms$locus[i], ":")[[1]][2], "-")[[1]][1])
  iso_ends =  c(iso_ends, strsplit(strsplit(random_isoforms$locus[i], ":")[[1]][2], "-")[[1]][2])
}

iso_starts =  as.numeric(iso_starts)
iso_ends =  as.numeric(iso_ends)


for (i in 1:nrow(random_isoforms)){
  iso_start = iso_starts[i]
  if (isoforms$strand[i] == "-"){
    iso_start = iso_ends[i]
  }
  
  if (sum(peaks1$V1 == iso_chromosomes[i]) > 0){
    correct_chromosome = peaks1$V1 == iso_chromosomes[i]
    distances = abs((peaks1$V2[correct_chromosome] + peaks1$V3[correct_chromosome])/2 - iso_start)
    
    random_isoforms$peak1_chr[i] <- peaks1$V1[correct_chromosome][head(order(distances), n=1)]
    random_isoforms$peak1_start[i] <- peaks1$V2[correct_chromosome][head(order(distances), n=1)]
    random_isoforms$peak1_end[i] <- peaks1$V3[correct_chromosome][head(order(distances), n=1)]
    random_isoforms$peak1_dist[i] <- head(sort(distances), n=1)
  }
  if (sum(peaks2$V1 == iso_chromosomes[i]) > 0){
    correct_chromosome = peaks2$V1 == iso_chromosomes[i]
    distances = abs((peaks2$V2[correct_chromosome] + peaks2$V3[correct_chromosome])/2 - iso_start)
    
    random_isoforms$peak2_chr[i] <- peaks2$V1[correct_chromosome][head(order(distances), n=1)]
    random_isoforms$peak2_start[i] <- peaks2$V2[correct_chromosome][head(order(distances), n=1)]
    random_isoforms$peak2_end[i] <- peaks2$V3[correct_chromosome][head(order(distances), n=1)]
    random_isoforms$peak2_dist[i] <- head(sort(distances), n=1)
  }
  if (sum(peaks3$V1 == iso_chromosomes[i]) > 0){
    correct_chromosome = peaks3$V1 == iso_chromosomes[i]
    distances = abs((peaks3$V2[correct_chromosome] + peaks3$V3[correct_chromosome])/2 - iso_start)
    
    random_isoforms$peak3_chr[i] <- peaks3$V1[correct_chromosome][head(order(distances), n=1)]
    random_isoforms$peak3_start[i] <- peaks3$V2[correct_chromosome][head(order(distances), n=1)]
    random_isoforms$peak3_end[i] <- peaks3$V3[correct_chromosome][head(order(distances), n=1)]
    random_isoforms$peak3_dist[i] <- head(sort(distances), n=1)
  }
}

random_isoforms$cl5 <- numeric(nrow(random_isoforms))
random_isoforms$cl3 <- numeric(nrow(random_isoforms))
random_isoforms$cl2 <- numeric(nrow(random_isoforms))


write.csv(random_isoforms, file = "random_isoforms_peak.csv",row.names=FALSE)