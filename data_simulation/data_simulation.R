# Run this first (before running hafoe) to generate the input files

library(seqinr)
library(stringr)

dir.create("data_simulation/plots", showWarnings = F)
fileName <- "input_files/AAV_all16_new.clustal_num"
format <- "clustal"

alignment <- seqinr::read.alignment(fileName, format = format, forceToLower = F)

set.seed(1)
reads_num <- 300
aav_num <- alignment$nb
fragment_size_min <- 100 
fragment_size_max <- 700 
min_orf <- 600
aln_len <- nchar(alignment$seq[1])


serotype_name <- c(
  "NC_002077.1",
  "NC_001401.2",  #old "ENA|J01901|J01901.1"
  "ENA|U48704|U48704.1" ,
  "NC_001829.1",
  "NC_006152.1",
  "ENA|AF028704|AF028704.1",
  "ENA|AF513851|AF513851.1",
  "ENA|AAN03857|AAN03857.1",
  "ENA|AAS99264|AAS99264.1",
  "ENA|AAT46337|AAT46337.1",
  "AY631966.1",
  "DQ813647.1",
  "EU285562.1",
  "AY242997.1", #AAVrh8
  "AY243015.1", #AAVrh10
  "AY243003.1") #AAVrh32

biased_aav_nums <- c(rep(2, 20), rep(14, 15), rep(6, 10), rep(7, 8), rep(9, 6), 
                     rep(8, 5), rep(13, 5), rep(1, 4), rep(5, 2), rep(3, 1), 
                     rep(15, 1))

seqs <- c()
seq_labels <- c()

i <- 1
while(i <= reads_num){
  #for (i in seq(1, reads_num)){
  cursor <- 1
  seq <- ""
  seq_label_list <- c()
  
  while (cursor < aln_len) {
    
    #choose random serotype 
    sero_idx <- sample(biased_aav_nums, 1, replace=T)  #1:aav_num uniform
    #choose random cut position
    cut <- sample((cursor + fragment_size_min) : (cursor + fragment_size_max), 1, replace=F) #cursor:aln_len
    #to reduce diversity make limitation on cut positions to be divisible by 10 
    while (cut %% 100 != 0){
      cut <- sample((cursor + fragment_size_min) : (cursor + fragment_size_max), 1, replace=F)
    }
    
    #if cut position is near the end, take the fragment until the end
    if (aln_len - cut <= fragment_size_min){
      cut <- aln_len
    }
    
    #append seq, to generate the chimeric sequence
    seq <- paste0(seq, substring(alignment$seq[which(alignment$nam == serotype_name[sero_idx])], cursor, cut))  #sero_idx -> 
    #append seq_label, to keep track of the chimeric sequence composition
    seq_label_list <- c(seq_label_list, rep(sero_idx, (cut - cursor + 1)))
    
    #update cursor
    cursor <- cut
    
  }
  
  #remove gap symbols from sequence and label list
  gap_positions <- unlist(gregexpr('-', seq))
  seq <- gsub("-", "", seq)
  
  seq_label_list <- seq_label_list[- gap_positions]
  
  seq_label <- paste(seq_label_list, collapse = " ")
  
  #check for the ORF length
  pos <- ORFik::findORFs(seq, startCodon = c("ATG"), minimumLength = min_orf)
  neg <- ORFik::findORFs(microseq::reverseComplement(seq), startCodon = c("ATG"), minimumLength = min_orf)
  if (length(pos) > 0 || length(neg) > 0){
    print(pos)
    print(neg)
    seqs <- c(seqs, seq)
    seq_labels <- c(seq_labels, seq_label)
    i <- i + 1
  } 
}  

df <- data.frame(chimeric_seq = seqs,
                 composition = seq_labels)
df['count'] <- sample(1:50, reads_num, replace = T)

save.image("./data_simulation/data_simulation.RData")


#CHIMERIC
#make chimeric library csv file
chimeric_library <- df[c('chimeric_seq', 'count', 'composition')]
chimeric_library <- as.data.frame(lapply(chimeric_library, rep, chimeric_library$count))
chimeric_library['count'] <- 1
index <- c(1:nrow(chimeric_library))
chimeric_library[, "index"] <- index
chimeric_library[, "X"] = paste0("AAV.", 100000 + as.numeric(chimeric_library$index))
chimeric_true_labels <- chimeric_library[c('X', 'composition', 'chimeric_seq')] #chimeric_seq?
colnames(chimeric_true_labels) <- c( 'X', 'Composition', 'Sequence')

write.csv(chimeric_true_labels, "input_files/Chimeric_lib_simulated_labels.csv", row.names = F)

chimeric_library <- chimeric_library[c('X', 'chimeric_seq','count')]
colnames(chimeric_library) <- c( 'X', 'Sequence', 'Count')

write.csv(chimeric_library, "input_files/Chimeric_lib_simulated.csv", row.names = F)


#ENRICHED
set.seed(1)
x0 <- rnorm(1000, mean = -1, sd = 0.5)
x <- x0[x0 >= -1]
x <- sample(x, 300)
x <- round(x, 2)

enriched_count <- floor((1+x) * df$count)
df['enriched_count'] <- enriched_count


write.fastq <- function(fastq_file, output_path){
  c = 1
  while (c <= nrow(fastq_file)){
    cat(paste0("@", fastq_file$Header[c]), 
        fastq_file$Sequence[c],
        "+",
        fastq_file$Quality[c],
        sep = "\n",
        file = output_path,
        append = TRUE)
    c = c + 1
  } 
}

#make enriched 1 library fastq file
enriched_library <- df[c('chimeric_seq', 'enriched_count')]
enriched_library <- as.data.frame(lapply(enriched_library, rep, enriched_library$enriched_count))
enriched_library['Quality'] <- strrep('~', nchar(enriched_library$chimeric_seq)) #highest quality score for PacBio
index <- c(1:nrow(enriched_library))
enriched_library[, "index"] <- index
enriched_library[, "Header"] = paste0("enriched.read", as.numeric(enriched_library$index))

enriched_library <- enriched_library[c('Header', 'chimeric_seq', 'Quality')]
colnames(enriched_library) <- c('Header', 'Sequence', 'Quality')

write.fastq(enriched_library, "input_files/Enriched_lib_simulated.fastq")
system("rm input_files/Enriched_lib_simulated.fastq.gz")
system("gzip input_files/Enriched_lib_simulated.fastq")


##PLOTS

#1. percent change distribution in enriched library with respect to chimeric library

x_df <- data.frame(x = x0)
p1 <- ggplot(x_df, aes(x=x0)) + 
  geom_histogram(aes(y=..density..), color="black", fill="grey", binwidth = 0.1) +
  geom_density() +
  annotate("rect", xmin = -3, xmax = -1, ymin = 0, ymax = 1, alpha = .8, fill = "white") +
  xlab("Fraction change") + 
  theme_minimal()

graphics.off()
pdf(file.path("data_simulation/plots/percent_change_distribution.pdf"), width=8, height=5)
print(p1)
dev.off()

#2. enriched library abundance distribution 

p2 <- ggplot(df, aes(x=enriched_count)) + 
  geom_histogram(color="black", fill="grey") +
  xlab("Abundance in enriched library") + 
  ylab("Count") + 
  theme_minimal()

graphics.off()
pdf(file.path("data_simulation/plots/enriched_library_abundance_distribution.pdf"), width=8, height=5)
print(p2)
dev.off()

#3. Serotype distribution barplot on the whole chimeric library data 

#get frequencies
s_all <- stringr::str_split(unlist(chimeric_true_labels[,'Composition'], 1), " ")

col_num_all <- max(unlist(lapply(s_all, length)))
for (i in seq_len(length(s_all))){
  if (length(s_all[[i]]) < col_num_all){
    s_all[[i]] <- c(s_all[[i]], rep("18", col_num_all - length(s_all[[i]]))) #gap
  }
}

matrix_all <- matrix(as.numeric(unlist(s_all)), ncol = col_num_all, byrow = TRUE)
rownames(matrix_all) <- chimeric_true_labels[, 'X']

frequencies = matrix(nrow = nrow(matrix_all), ncol = 0)
for (i in (0:18)){
  frequencies = cbind(frequencies, rowSums(matrix_all == i))
}

colnames(frequencies) <- c("no alignment","AAV1","AAV2","AAV3","AAV4","AAV5","AAV6","AAV7","AAV8","AAV9",
                           "AAV10","AAV11","AAV12","AAV13","AAVrh8","AAVrh10","AAVrh32", "multiple alignment", "gap")
frequencies <- subset(frequencies, select = -c(gap))

frequencies_final = colSums(frequencies)

#Abundance of AAV serotypes in the chimeric library
serotypes_freq <- as.data.frame(frequencies_final)
colnames(serotypes_freq) <- c("Freq")
serotypes_freq['Freq(%)'] <- round(serotypes_freq$Freq*100/sum(serotypes_freq$Freq), 2)
serotypes_freq <- serotypes_freq[order(serotypes_freq$`Freq(%)`, decreasing = T), ]
#print(serotypes_freq[-1])

#plot

col = c("#D3D3D3", "#A6CEE3", "#1F78B4", "#B2DF8A", "#555555", "#33A02C",
        "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
        "#FFFF99", "#B15928", "yellow", "#999999", "#a70000", "black", "white")
legend_name <- c("no alignment", "AAV1", "AAV2", "AAV3", "AAV4", "AAV5",
                 "AAV6", "AAV7", "AAV8", "AAV9", "AAV10", "AAV11", "AAV12", 
                 "AAV13", "AAVrh8", "AAVrh10", "AAVrh32", "multiple alignment", "gap")
col_df <- data.frame(col = col)
rownames(col_df) <- legend_name


plot.serotype.frequency <- function(serotypes_freq, col_df, library_name = ""){
  serotypes_freq['Name'] <- rownames(serotypes_freq)
  col_ordered <- col_df[serotypes_freq[order(serotypes_freq$`Freq(%)`, decreasing = T), "Name"],]
  
  p <- ggplot2::ggplot(data=serotypes_freq, 
                       aes(x=factor(Name, serotypes_freq[order(`Freq(%)`, decreasing = T), "Name"]), 
                           y=`Freq(%)`, 
                           order=`Freq(%)`)) +
    geom_bar(stat="identity", fill=col_ordered)+
    geom_text(aes(label=`Freq(%)`), vjust=-0.3, size=3.5) + 
    labs(#title=paste0("Distribution of AAV serotypes in ", library_name),
      x ="", y = "Frequency (%)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10)) 
  print(p)
}

graphics.off()
pdf(file.path("data_simulation/plots/serotype_distribution_chimeric_lib_true.pdf"), width=8, height=5)
plot.serotype.frequency(serotypes_freq, col_df = col_df, library_name = "chimeric library")
dev.off()

