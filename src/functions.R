
#plot sequence length distribution histogram
plot.histogram <- function(seq_data, library_name = "", variant_orf = "", additional_info = "", bin_size = 20){
  fig <- plotly::plot_ly(x = nchar(seq_data), type = "histogram", xbins = list(size=bin_size)) 
  fig <- fig %>% layout(title= list(text = paste0(library_name, " library ", variant_orf, " length distribution", additional_info))
                        #,xaxis = list(range=c(0,1000))
  )
  fig
}

plot.histogram_ <- function(x, library_name = "", variant_orf = "", additional_info = "", bin_size = 20){
  fig <- plotly::plot_ly(x = x, type = "histogram", xbins = list(size=bin_size)) 
  fig <- fig %>% layout(title= list(text = paste0(library_name, " library ", variant_orf, " length distribution", additional_info))
                        #,xaxis = list(range=c(0,1000))
  )
  fig
}

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

#filter by minimum ORF length and max sequence length (ORF length >= min_orf (600 aa) and sequence length <= max_seq(3000))
filter.by.orf <- function(file, file_path, min_orf = 600, max_seq = 3000, format = "fasta", library_name = ""){
  cat(paste0("File: ", file_path, "\n"))
  cat(paste0("Minimum ORF length: ", min_orf, "\n"))
  cat(paste0("Maximum sequence length: ", max_seq, "\n"))
  cat(paste0("Number of sequences before filtering: ", nrow(file), "\n"))
  
  dna <- readDNAStringSet(file_path, format = format)
  dna_rc <- Biostrings::reverseComplement(dna)
  
  #Positive strand
  #minimumLength -> number of amino acids not bases
  pos <- findORFs(dna, longestORF = TRUE, startCodon = c("ATG"), minimumLength = min_orf)
  #Get sequences 
  gr_pos <- unlist(pos, use.names = TRUE)
  gr_pos <- GRanges(seqnames = names(dna)[as.integer(names(gr_pos))], ranges(gr_pos), strand = "+")
  names(gr_pos) <- make.names(gr_pos@ranges@NAMES, unique=TRUE)
  gr_pos_df <- as.data.frame(gr_pos)
  #Give proper names:
  #Ex. m54278_211112_013512/4391649/ccs/ORF/1706-3910/2205/+   
  names(gr_pos) <- paste0(gr_pos_df$seqnames, "/", "ORF/", gr_pos_df$start, "-", gr_pos_df$end, "/", gr_pos_df$width, "/+")
  orf_seq_pos <- Biostrings::getSeq(dna, gr_pos)
  
  orf_info_df_pos <- as.data.frame(gr_pos)
  orf_seq_df_pos <- as.data.frame(orf_seq_pos)
  
  #Negative strand
  neg <- findORFs(dna_rc, longestORF = TRUE, startCodon = c("ATG"), minimumLength = min_orf)
  if (length(neg) > 0){
  #Get sequences
    gr_neg <- unlist(neg, use.names = TRUE)
    gr_neg <- GRanges(seqnames = names(dna_rc)[as.integer(names(gr_neg))], ranges(gr_neg), strand = "-")
    names(gr_neg) <- make.names(gr_neg@ranges@NAMES, unique=TRUE)
    gr_neg_df <- as.data.frame(gr_neg)
    #Give proper names: 
    names(gr_neg) <- paste0(gr_neg_df$seqnames, "/", "ORF/", gr_neg_df$start, "-", gr_neg_df$end, "/", gr_neg_df$width, "/-")
    orf_seq_neg <- getSeq(dna_rc, gr_neg)
    
    orf_info_df_neg <- as.data.frame(gr_neg)
    orf_seq_df_neg <- as.data.frame(orf_seq_neg)
  } else {
    orf_info_df_neg <- c()
    orf_seq_df_neg <- c()
  }
  
  #Combine pos & neg
  orf_info_df_all <- rbind(orf_info_df_pos, orf_info_df_neg)
  orf_seq_df_all <-  rbind(orf_seq_df_pos, orf_seq_df_neg)
  
  #Take max ORF (group by seqname, take the row with max orf length) (some variants have multiple orfs > 1.8k)
  group <- as.data.table(orf_info_df_all, keep.rownames = TRUE)
  orf_info_df_all_max <- group[group[, .I[width == max(width)], by=seqnames]$V1]
  #Check for duplicates
  #which(duplicated(orf_info_df_all$seqnames))
  
  #Remove duplicate rows by seqname leaving only first: leave only one ORF of same variant having same length
  orf_info_df_all_max <- orf_info_df_all_max[!duplicated(orf_info_df_all_max$seqnames),]
  
  orf_seq_df_all['name'] <- rownames(orf_seq_df_all)
  orf_seq_df_all_max <- orf_seq_df_all[rownames(orf_seq_df_all) %in% orf_info_df_all_max$rn,]
  
  #Contains orfs
  fasta_orf_seq <- orf_seq_df_all_max[c('name', 'x')]
  colnames(fasta_orf_seq) <- c("Header", "Sequence")
  
  #Contains original sequences
  fasta_orf <- dplyr::inner_join(file, orf_info_df_all_max, by = c("Header" = "seqnames")) 
  fasta_orf <- fasta_orf[c('rn', 'Sequence')]  
  colnames(fasta_orf) <- c('Header', 'Sequence')
  
  #Filter by sequence length (max_seq)
  fasta_orf_filtered <- fasta_orf[nchar(fasta_orf$Sequence) <= max_seq, ] 
  fasta_orf_seq_filtered <- fasta_orf_seq[fasta_orf_seq$Header %in% fasta_orf_filtered$Header,]
  
  #Store filtered orf sequences into fasta file
  writeFasta(fasta_orf_seq_filtered, file.path(output.dir, "files", paste0(library_name, "_ORF.fasta")))
  
  if (format == "fastq"){
    #Add sequence original name, orf start, end positions as columns
    orf_seq_df_all_max$name_original <- sub("/ORF.*", "", orf_seq_df_all_max$name)
    ranges <- str_split(orf_seq_df_all_max$name, pattern = "/", simplify = TRUE)[,3]  #[,5]
    start <- str_split(ranges, pattern = "-", simplify = TRUE)[,1]
    end <- str_split(ranges, pattern = "-", simplify = TRUE)[,2]
    orf_seq_df_all_max['start'] <- start
    orf_seq_df_all_max['end'] <- end
    
    #Contains orfs
    fastq_orf_seq <- inner_join(file, orf_seq_df_all_max, by = c("Header" = "name_original")) 
    fastq_orf_seq['quality_orf'] <- substring(fastq_orf_seq$Quality, fastq_orf_seq$start, fastq_orf_seq$end)
    fastq_orf_seq <- fastq_orf_seq[c('name', 'x', 'quality_orf')]  
    colnames(fastq_orf_seq) <- c('Header', 'Sequence', 'Quality')
    
    fastq_orf_seq_filtered <- fastq_orf_seq[fastq_orf_seq$Header %in% fasta_orf_filtered$Header,]
    
    # write.fastq(fastq_file = fastq_orf_seq_filtered, 
    #             output_path = file.path(output.dir, "files", paste0(library_name, "_ORF.fastq")))
    
    #quote error
    #writeFastq(fastq_orf_seq_filtered, file.path(output.dir, "files", paste0(library_name, "_ORF.fastq")))
  }
  
  #Store dataframe outputs in list (for plotting histogram)
  out <- list(fasta_orf_filtered, fasta_orf_seq_filtered)      
  
  cat(paste0("Number of sequences after filtering: ", nrow(fasta_orf_seq_filtered), "\n"))
  cat("ORF length summary statistics\n")
  print(summary(nchar(fasta_orf_seq_filtered$Sequence)))
  cat("\n")
  
  return(out)
}

choose.new.representatives <- function(sizes_path, members_path, chimeric_library, fasta_orf_seq_filtered_chim, fasta_path, csv_path){
  cat("Choosing a new representative ORF sequence for each cluster based on: \n")
  cat("\t1. greatest sequence abundance in the input chimeric library\n")
  cat("\t2. maximum ORF length\n")
  cat(paste0("Output path: ", fasta_path, "\n"))

  #use the generated members_ordered.csv and sizes.csv to separate variants by clusters
  sizes <- read.csv(sizes_path, header = FALSE)
  sizes['sum'] <- cumsum(sizes$V1)
  members <- read.csv(members_path, header = FALSE)

  #for each cluster choose a representative sequence (1.max count, 2.max ORF length)
  #for loop which works for each cluster
  for (i in seq_len(nrow(sizes))){
    #only for 1st cluster
    if (i == 1){
      #make dataframe containing only the variants which are in the 1st cluster (cluster 0)
      #example of a variant name: AAV.211507/ORF/754-2958/2205/-      I split by / to get   AAV.211507  ORF  754-2958  2205  -     this info in separate columns, to be able to use width/length of orfs
      temp <- as.data.frame(stringr::str_split(members[1:sizes$sum[1],], pattern = "/", simplify = TRUE))
      colnames(temp) <- c("name", "ORF", "pos", "width", "strand")


      #join to initial data by variant name to also take the count column
      temp <- dplyr::inner_join(temp, chimeric_library[c('X', 'Count')],  by = c("name" = "X"))
      #take only the row(s)/variant(s) having max count
      max_count <- temp[temp$Count == max(temp$Count),]
      #from them take only the row(s)/variant(s) having max orf length
      max_width <- max_count[max_count$width == max(max_count$width),]
      #make a dataframe rep containing cluster number, representative's name, its orf length and abundance count
      rep <- data.frame(Cluster = i-1,
                        Representative = max_width[1, 'name'], #[1,'name'] 1 means I take only 1st row just in case multiple having max count
                        Width = max_width[1,'width'],
                        Count = max_width[1, 'Count'])
      #make df of cluster members
      cluster_members <- data.frame(Cluster = i-1,
                                    Representative = max_width[1, 'name'])
      cluster_members$Members <- list(temp$name)
    }
    #for the rest of clusters
    else{
      #same thing for the rest of clusters (cluster 1, then cluster 2, etc.)
      temp <- as.data.frame(stringr::str_split(members[(sizes$sum[i-1]+1) : sizes$sum[i], ],
                                            pattern = "/", simplify = TRUE))
      colnames(temp) <- c("name", "ORF", "pos", "width", "strand")

      temp <- dplyr::inner_join(temp, chimeric_library[c('X', 'Count')],  by = c("name" = "X"))
      max_count <- temp[temp$Count == max(temp$Count),]
      max_width <- max_count[max_count$width == max(max_count$width),]
      #populates the dataframe rep created above for the 1st cluster
      rep <- rbind(rep,
                   c(i-1,
                     max_width[1, 'name'],
                     max_width[1, 'width'],
                     max_width[1, 'Count']))
      
      cluster_members <- rbind(cluster_members,
                               c(i-1,
                                 max_width[1, 'name'],
                                 NA))
      cluster_members$Members[i] <- list(temp$name)
      
    }
  }
  
  
  #from initial fasta file take only the representative variants
  fasta_orf_seq_filtered_names <- stringr::str_split(fasta_orf_seq_filtered_chim$Header, pattern = "/", simplify = TRUE)[,1]
  
  rep_fasta <- fasta_orf_seq_filtered_chim[fasta_orf_seq_filtered_names %in% rep$Representative,]
  
  #create a new fasta file containing representatives
  microseq::writeFasta(rep_fasta,
                       fasta_path)
  write.table(rep$Representative,
              csv_path,
              row.names = F,
              col.names = F,
              quote = F,
              sep = ",")
  return(cluster_members)
}

#change based on fastq file
make.fastq.files <- function(read_length, fasta.path, library_name = "", overlap = F, step_size = read_length){
  cat("Chopping each representative sequence into fragments of fixed size and storing into fastq files \n")
  cat(paste0("Fragment length: ", read_length, "\n"))
  cat(paste0("Use overlapping fragments (TRUE/FALSE): ", overlap, "\n"))
  if (overlap == F){
    step_size = read_length
  } 
  cat(paste0("Step size: ", step_size, "\n"))
  cat(paste0("Output path: ", file.path(output.dir, "files/variant_description", library_name, "fastq"), "\n"))
  if (step_size > read_length){
    cat("Warning: step size is greater than fragment length, some regions of the sequence will be skipped!")
  }

  data <- microseq::readFasta(fasta.path)  
  data$name <- matrix(unlist(stringr::str_split(data$Header, "/")), ncol = 5, byrow = TRUE)[,1]
  
    
  dir.create(file.path(output.dir, "files/variant_description"), showWarnings = F)
  dir.create(file.path(output.dir, "files/variant_description", library_name), showWarnings = F)
  dir.create(file.path(output.dir, "files/variant_description", library_name, "fastq"), showWarnings = F)
  splitter <- function(rl, ss, c){
    pos = 0     
    read_counter = 10001 #10^nchar(as.character(ceiling(max(nchar(data$Sequence))/ss)))+1
    
    # GENERAL CASE
    while(pos < nchar(data[c, "Sequence"])){  #fragments of length less than rl may be generated at the end of the sequence
      cat(paste0("@", data[c, "Header"], "_", read_counter),
          substring(data[c, "Sequence"], pos + 1, pos + as.numeric(rl)),
          "+",
          strrep("h", nchar(substring(data[c, "Sequence"], pos + 1, pos + as.numeric(rl)))),
          sep = "\n",
          file = file.path(output.dir, "files/variant_description", library_name, "fastq",
                           paste0(data[c,"name"], ".fastq")),
          append = TRUE)
      pos = pos + as.numeric(ss)
      read_counter = read_counter + 1
    }
  }
  
  c = 1
  while (c <= nrow(data)){
    splitter(read_length, step_size, c)
    c = c + 1
  }
}



#parents_df change general case from parents fasta or separate input paramenter from .sh
neighbor.joining <- function(parents_df, path_fasta, path_mapped, path_unmapped, read_length, library_name = "", overlap = F, step_size = read_length){
  start_time <- Sys.time()
  
  
  cat(paste0("Output path: ", file.path(output.dir, 
                                        "files/variant_description", 
                                        library_name, 
                                        paste0(library_name, "_variant_description.csv")), "\n"))
  
  if (overlap == F){
    step_size = read_length
  }
  
  #make directory for output file if not created yet
  dir.create(file.path(output.dir, "files/variant_description"), showWarnings = F)
  dir.create(file.path(output.dir, "files/variant_description", library_name), showWarnings = F)
  
  data <- microseq::readFasta(path_fasta)  
  
  #path_mapped - path of directory containing csv files generated from bam files(only mapped reads). 
  #In these csv files there are 2 columns(read name and serotype name)
  
  #list of file names in that directory
  files_mapped <- list.files(path = path_mapped, pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE)
  
  #path_unmapped - path of directory containing csv files generated from fastq files(only unmapped reads). 
  #In these csv files there is 1 column(read name)
  
  #list of file names in that directory
  files_unmapped <- list.files(path = path_unmapped, pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE)
  
  files_names <- sub('\\.csv$', '', list.files(path = path_mapped, pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE))
  #binds these two lists of files by column (resulting 2 column and each row corresponding to one variant). Needed to connect mapped and unmapped files for each variant
  file_pairs = cbind(files_mapped, files_unmapped, files_names)  

  isEmpty <- function(x) {
    return(length(x)==0)
  }
  
  
  #empty matrix which will be filled in the for loop
  m <- matrix(data=list(),
              nrow=nrow(file_pairs),
              ncol=2)
  read <- c()
  options <- list()
  
  #for loop iterating over all query variants' mapped and unmapped files simultaneously 
  for (j in seq_len(nrow(file_pairs))) {
    
    #there where some empty files in mapped and unmapped files, so I used if/else statements to ignore them
    #MAPPED
    if (file.size(file.path(path_mapped, file_pairs[[j, 1]])) > 0){
      
      #open csv file, do not use first row as header, use tab to seperate columns. read as matrix
      data_mapped <- as.matrix(read.csv(file.path(path_mapped, 
                                                  file_pairs[[j, 1]]), 
                                        header = F, 
                                        sep = "\t"))
      #specify column names of matrix
      colnames(data_mapped) <- c("variant_read", "serotype")
      
      #extract read numbers from the first column of data_mapped matrix (ex: AAV.100001.1003 --> 3)
      #read_num <- as.numeric(substring(data_mapped[,1], nchar(data_mapped[,1])-4+1, nchar(data_mapped[,1]))) - 1000
      read_num <- as.numeric(substring(data_mapped[,1], 
                                       nchar(data_mapped[,1]) - 5 + 1, 
                                       nchar(data_mapped[,1]))) - 10000
      
      #add this list of read numbers to data_mapped matrix as a column
      data_mapped <- cbind(data_mapped, read_num)
    } else {
      read_num <- c(-1)
    }
    
    #UNMAPPED
    #if file not empty
    if (file.size(file.path(path_unmapped, file_pairs[[j, 2]])) > 0){ 
      #same for unmapped files
      data_unmapped <- as.matrix(read.csv(file.path(path_unmapped, file_pairs[[j, 2]]), 
                                          header = F, 
                                          sep = "\t"))  
      
      colnames(data_unmapped) <- c("variant_read")
      
      #read_num_un <- as.numeric(substring(data_unmapped[,1], nchar(data_unmapped[,1])-4+1, nchar(data_unmapped[,1]))) - 1000
      read_num_un <- as.numeric(substring(data_unmapped[,1], 
                                          nchar(data_unmapped[,1])- 5 + 1, 
                                          nchar(data_unmapped[,1]))) - 10000
    } else {
      read_num_un <- c(-1)
    }
    
    # #empty vector to be filled in the for loop
    # v = list()
    # #for loop iterates over read numbers of corresponding variant 
    # for (i in 1: max(max(read_num), max(read_num_un))) {
    #   #if the number of reads is present in the read_num list (list for mapped reads)
    #   if (i %in% read_num) {
    #     #then the vector v is appended by the serotype number in corresponding position 
    #     s <- data_mapped[data_mapped[,"read_num"] == i, "serotype"]
    #     v[[i]] <- parents_df[serotype_name %in% s, "serotype_num"]
    #     
    #   }
    #   #if the number of reads is present in the read_num_un list(list for unmapped reads)
    #   if (i %in% read_num_un) {
    #     #then  the vector v is appended by 0 in corresponding position 
    #     v[[i]] = 0 
    #   } 
    # }
    

    #new
    variant_length <- nchar(data[j,]$Sequence) #!!! change, order may not be the same
    v = list()
    for (p in seq_len(variant_length)){
      serotypes_all <- c()
      for (z in 1: max(max(read_num), max(read_num_un))){
        if (p %in% seq(((z-1)*step_size)+1, ((z-1)*step_size)+read_length)){ #read start end positions in variant sequence
          
          #if the number of reads is present in the read_num list (list for mapped reads)
          if (z %in% read_num) {
            #then the vector v is appended by the serotype number in corresponding position 
            s <- data_mapped[data_mapped[,"read_num"] == z, "serotype"]
            s_num <- parents_df[serotype_name %in% s, "serotype_num"]
            
          }
          #if the number of reads is present in the read_num_un list(list for unmapped reads)
          if (z %in% read_num_un) {
            #then  the vector v is appended by 0 in corresponding position 
            s_num = 0 
          } 
          
          serotypes_all <- c(serotypes_all, s_num)
          serotypes_all <- unique(serotypes_all)
          
          #remove 0 if there are serotypes identified for that position
          if (length(serotypes_all) > 1){
            serotypes_all <- serotypes_all[serotypes_all != 0]
          }
          
          v[[p]] <- serotypes_all
        }
      }
    }
    
    # get unique consecutive values from v and lengths of them in a separate vector
    seen <- v[[1]]
    v_consecutive_unique <- v[1]
    count <- 1 
    consecutive_counts <- c()
    for (i in v[2:length(v)]){
      if (!identical(i, seen)){
        v_consecutive_unique <- append(v_consecutive_unique, list(i))
        consecutive_counts <- c(consecutive_counts, count)
        seen <- i
        count <- 1
      } else {
        count <- count + 1 
      }
    }
    consecutive_counts <- c(consecutive_counts, count)
    
    # this gives back the list of sequence length
    #rep(v_consecutive_unique, consecutive_counts)

    
    ###neighbor joining
    # v <- v_consecutive_unique
    # s <- 1
    # while(s > 0){
    #   
    #   for (i in 1:(length(v)-1)){
    #     #if there is intersect and they are not the same assign   the intersect to both neighbors
    #     if(!setequal(v[i], v[i+1])){
    #       if(!isEmpty(intersect(v[[i]], v[[i+1]]))){
    #         v[[i]] = intersect(v[[i]], v[[i+1]])
    #         v[[i+1]] = intersect(v[[i]], v[[i+1]])
    #       }
    #     }
    #   }
    #   
    #   s <- 0
    #   for (i in 1:(length(v)-1)){
    #     if( !setequal(v[i], v[i+1]) ) {
    #       s <- s + !isEmpty(intersect(v[[i]], v[[i+1]]))
    #     }
    #   }
    # }
    # 
    
    ###neighbor joining
    #new
    v <- v_consecutive_unique
    
    # check only the positions having multiple serotypes that need to be resolved
    multiple_idx <- which(lapply(v, length) > 1)
    
    s <- 1
    while(s > 0){
      for (i in multiple_idx){
        #if there is intersect and they are not the same assign the intersect to both neighbors
        #compare with the left neighbor
        if((i != 1) & (!setequal(v[i], v[i-1]))){
          if((!isEmpty(intersect(v[[i]], v[[i-1]])))){
            intersect_ <- intersect(v[[i]], v[[i-1]])
            v[[i]] <- intersect_
            v[[i-1]] <- intersect_
          }
        }
        #then compare with the right neighbor
        if((i != length(v)) & (!setequal(v[i], v[i+1]))){
          if(!isEmpty(intersect(v[[i]], v[[i+1]]))){
            intersect_ <- intersect(v[[i]], v[[i+1]])
            v[[i]] <- intersect_
            v[[i+1]] <- intersect_
          }
        }
      }
      
      # update multiple_idx to get rid of unnecessary positions
      multiple_idx <- which(lapply(v, length) > 1)
      
      s <- 0
      for (i in multiple_idx){
        if (s > 0) 
        {break}
        if ((i != 1) & (!setequal(v[i], v[i-1]))) {      # need to also check the left neighbor
          s <- s + !isEmpty(intersect(v[[i]], v[[i-1]])) 
        } 
        if((i != length(v)) & (!setequal(v[i], v[i+1]) )) { # then the right neighbor
          s <- s + !isEmpty(intersect(v[[i]], v[[i+1]]))
        } 
      }
    }
    
    
    for (i in 1:length(v)){
      if (length(v[[i]]) > 1){
        v[[i]] <- 17
      }
    }
    
    # make v of sequence length
    v <- rep(v, consecutive_counts)
    
    #after the vector v is made, the empty matrix mat defined above is appended. 1st column is the name of variant,
    #2nd column is a string of serotype numbers separated by spaces
    m[[j, 1]] <- file_pairs[[j, 3]]
    m[[j, 2]] <- paste(v, collapse = " ")
    
  }
  
  end_time <- Sys.time()
  # end_time - start_time
  

  s <- stringr::str_split(unlist(m[,2], 1), " ")
  # col_num <- ceiling(max(nchar(data$Sequence))/step_size)
  #new
  col_num <- max(nchar(data$Sequence))
  for (i in seq_len(length(s))){
    if (length(s[[i]]) < col_num){ 
      s[[i]] <- c(s[[i]], rep("18", col_num - length(s[[i]]))) #gap
    }
  }
  
  m_new <- matrix(as.numeric(unlist(s)), ncol = col_num, byrow = TRUE)
  rownames(m_new) <- m[,1]
  
  write.table(m_new, 
              file.path(output.dir, 
                        "files/variant_description", 
                        library_name, 
                        paste0(library_name, "_variant_description.csv")), 
              quote = F, 
              col.names = F)
}


plot.cluster.size.distribution <- function(cluster_members, size_thresh, library_name = ""){
  cluster_members["Size"] <- unlist(lapply(cluster_members$Members, length))
  data <- cluster_members[-3]
  data <- data[data$Size > size_thresh,]
  
  p <- ggplot2::ggplot(data=data, 
                      aes(x=factor(Representative, data[order(Size, decreasing = T), "Representative"]), 
                          y=Size, 
                          order=Size)) +
        geom_bar(stat="identity")+
        geom_text(aes(label=Size), vjust=-0.3, size=2) +
        labs(title=paste0("Cluster size distribution in ", library_name, "\n(Cluster size > ", size_thresh, ")"),
             x ="", y = "Cluster size") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) 
  print(p)
}



plot.variant.description <- function(matrix, col_df, library_name = ""){
  gplots::heatmap.2(matrix, 
                    dendrogram='none', 
                    Colv=FALSE, 
                    Rowv=FALSE, #TRUE 
                    trace="none", 
                    breaks = seq(-0.5, 18.5, 1), 
                    col = col_df$col, 
                    key = FALSE, 
                    cexRow=0.7)
  title(paste0("Variant description of ", library_name), 
        line = -2, 
        adj = 0.6)
  legend(x="bottomleft", 
         legend=rownames(col_df), 
         fill=col_df$col,  
         title = "AAV serotypes", 
         title.adj = 0.2, 
         inset=c(-.07, -.07), 
         xpd=TRUE,
         box.lwd = 0, 
         cex = 0.7)

}


plot.variant.description.conserved <- function(matrix,  col_df, identity_ranges, library_name = ""){
  gplots::heatmap.2(matrix, 
                    dendrogram='none', 
                    Colv=FALSE, 
                    Rowv=TRUE, 
                    trace="none", 
                    breaks = seq(-0.5, 18.5, 1),
                    col = col_df$col,
                    key = FALSE, 
                    cexRow=0.7,
                    add.expr = list(rect(xleft = identity_ranges$start_nt, 
                                         xright = identity_ranges$end_nt, 
                                         ybottom = par("usr")[3], ytop = par("usr")[4], 
                                         border = NA, 
                                         col = adjustcolor("blue", alpha = 0.2)),
                                    rect(xleft = aav2_vr_ranges$start_nt,   #VRs
                                         xright = aav2_vr_ranges$end_nt, 
                                         ybottom = par("usr")[3], ytop = par("usr")[4], 
                                         border = NA, #"red", 
                                         density = 20, 
                                         col = adjustcolor("red", alpha = 0.5))))
  
  title(paste0("Variant description of ", library_name), 
        line = -2, 
        adj = 0.6)
  legend(x="bottomleft", 
         legend=rownames(col_df), 
         fill=col_df$col,  
         title = "AAV serotypes", 
         title.adj = 0.2, 
         inset=c(-.07, -.07), 
         xpd=TRUE,
         box.lwd = 0, 
         cex = 0.7)
}



get.frequency.table <- function(chimeric_library, nj_matrix, cluster_members, output_path){
  #combined counts of chimeric lib representatives
  counts_chim <- data.frame(Representative = NA, 
                            Chimeric.Count = NA)
  
  for (rep in cluster_members$Representative){
    s <- sum(chimeric_library[chimeric_library$X %in% unlist(cluster_members$Members[cluster_members$Representative == rep]), c('Count')])
    counts_chim <- rbind(counts_chim, c(rep, s))
  }
  
  counts_chim <- na.omit(counts_chim)
  counts_chim$Chimeric.Count <- as.numeric(counts_chim$Chimeric.Count)
  
  counts_chim['Chimeric.Percentage'] <- round((100 * (counts_chim$Chimeric.Count / sum(counts_chim$Chimeric.Count))), 2)
  counts_chim <- counts_chim %>% arrange(desc(Chimeric.Count))
  
  write.csv(counts_chim, file.path(output_path, "chimeric_lib_rep_counts.csv"), row.names = FALSE, quote = F)
  
  #serotype frequencies per representative variant
  frequencies = matrix(nrow = nrow(nj_matrix), ncol = 0)
  for (i in (0:18)){
    frequencies = cbind(frequencies, rowSums(nj_matrix == i))
  }
  
  #serotype freq * representative combined count
  frequencies_with_count = matrix(nrow = 0, ncol = ncol(frequencies))
  for (i in rownames(frequencies)){
    frequencies_with_count = rbind(frequencies_with_count, 
                                   counts_chim[counts_chim$Representative == i, "Chimeric.Count"] * frequencies[i,])
  }
  
  frequencies_with_count_df = as.data.frame(frequencies_with_count)
  colnames(frequencies_with_count_df) <- c("no alignment","AAV1","AAV2","AAV3","AAV4","AAV5","AAV6","AAV7","AAV8","AAV9",
                                          "AAV10","AAV11","AAV12","AAV13","AAVrh8","AAVrh10","AAVrh32", "multiple alignment", "gap")
  rownames(frequencies_with_count_df) <- rownames(frequencies)
  
  frequencies_with_count_df <- subset(frequencies_with_count_df, select = -c(gap))
  
  write.csv(frequencies_with_count_df, file.path(output_path, "chimeric_lib_serotype_counts.csv"), quote = F)
  
  #for all representatives
  frequencies_final = colSums(frequencies_with_count_df)
  
  #summary frequency table
  cat("\nAbundance of AAV serotypes in the chimeric library\n")
  serotypes_freq <- as.data.frame(frequencies_final)
  colnames(serotypes_freq) <- c("Freq")
  serotypes_freq['Freq(%)'] <- round(serotypes_freq$Freq*100/sum(serotypes_freq$Freq), 2)
  serotypes_freq <- serotypes_freq[order(serotypes_freq$`Freq(%)`, decreasing = T), ]
  print(serotypes_freq[-1])
  
  return(serotypes_freq)
}




plot.serotype.frequency <- function(serotypes_freq, col_df, library_name = ""){
  serotypes_freq['Name'] <- rownames(serotypes_freq)
  col_ordered <- col_df[serotypes_freq[order(serotypes_freq$`Freq(%)`, decreasing = T), "Name"],]
  
  p <- ggplot2::ggplot(data=serotypes_freq, 
                        aes(x=factor(Name, serotypes_freq[order(`Freq(%)`, decreasing = T), "Name"]), 
                            y=`Freq(%)`, 
                            order=`Freq(%)`)) +
          geom_bar(stat="identity", fill=col_ordered)+
          geom_text(aes(label=`Freq(%)`), vjust=-0.3, size=3.5) + 
          labs(title=paste0("Distribution of AAV serotypes in ", library_name),
               x ="", y = "Frequency (%)") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  print(p)
}
  

plot.cluster.abundance <- function(file_path, size_thresh, library_name = ""){
  data = read.csv(file_path)
  data <- data[data$Chimeric.Count > size_thresh,]
  
  p <- ggplot2::ggplot(data=data, 
                        aes(x=factor(Representative, data[order(Chimeric.Count, decreasing = T), "Representative"]), 
                            y=Chimeric.Count, 
                            order=Chimeric.Count)) +
          geom_bar(stat="identity")+
          geom_text(aes(label=Chimeric.Count), vjust=-0.3, size=2) +
          labs(title=paste0("Cluster abundance in ", library_name, "\n(Combined count > ", size_thresh, ")"),
               x ="", y = "Combined count of variants") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) 
  
  print(p)
}

get.counts.table <- function(chim_file_path, output_path){
  counts_chim = read.csv(chim_file_path)
  
  sizes_en1 <- read.csv(file.path(output_path,"clstr_enriched1_lib/cluster_sizes.csv"), header = FALSE)
  sizes_en1$V1 <- sizes_en1$V1 - 1
  sizes_en1['sum'] <- cumsum(sizes_en1$V1)
  members <- read.csv(file.path(output_path,"clstr_enriched1_lib/members_ordered.csv"), header = FALSE)
  reps <- read.csv(file.path(output_path,"clstr_enriched1_lib/representatives.csv"), header = FALSE)
  
  
  rep_size_en1 <- data.frame(Representative = reps$V1, Size = sizes_en1$V1)
  rep_size_en1_ordered <- rep_size_en1[order(rep_size_en1$Size, decreasing = T),]
  row.names(rep_size_en1_ordered) <- NULL
  
  counts_en1 <- data.frame(Representative = stringr::str_split(rep_size_en1_ordered$Representative, pattern = "/", simplify = TRUE)[,1], 
                      Enriched1.Count = rep_size_en1_ordered$Size, 
                      Enriched1.Percentage = round((100 * (rep_size_en1_ordered$Size / sum(rep_size_en1_ordered$Size))), 2))
  counts_en1 <- counts_en1 %>% arrange(desc(Enriched1.Count))
  
  counts_all <- plyr::join_all(list(counts_chim, counts_en1), by='Representative', type='full')
  counts_all <- counts_all %>% arrange(desc(Chimeric.Count))
  
  #ratio/ enrichment rate
  counts_all['Ratio_En1Chim'] <- round(counts_all$Enriched1.Percentage / counts_all$Chimeric.Percentage, 2)
  counts_all$Ratio_En1Chim[which(!is.finite(counts_all$Ratio_En1Chim))] <- 0
  
  
  chimeric_orf <- microseq::readFasta(file.path(output_path,"Chimeric_ORF.fasta"))
  names <- as.data.frame(matrix(unlist(stringr::str_split(chimeric_orf$Header, pattern = "/")), 
                                ncol = 5, byrow = TRUE))[[1]]
  chimeric_orf['X'] <- names
  
  counts_all <- dplyr::left_join(counts_all, chimeric_orf[c("X", "Sequence")], by = c("Representative" = "X"))

  write.csv(counts_all, file.path(output_path, "counts.csv"), quote = F, row.names = F)
}

plot.top.reps.in.enrichedlib <- function(counts_file_path, topn_thresh, library_name = ""){
  counts <- read.csv(counts_file_path)
  counts <- counts %>% arrange(desc(Enriched1.Count))
  fig <- plotly::plot_ly(x = counts[1:30, "Enriched1.Count"] , y = counts[1:30, "Representative"], type = 'bar') 
  fig <- fig %>% layout(title= list(text = paste0("Top ", topn_thresh, " representative sequences after clustering \n", library_name)),  
                        xaxis = list(title = 'Number of reads in the cluster'),
                        yaxis = list(title = list(text ='Representative sequence'), 
                                     categoryorder = "array",
                                     categoryarray = ~counts$Representative)
  )
  print(fig)
}


plot.top.reps.in.enrichedlib.ggplot <- function(counts_file_path, topn_thresh, library_name = ""){
  counts <- read.csv(counts_file_path)
  counts <- counts %>% arrange(desc(Enriched1.Count))
  counts <- counts[1:30,]
  level_order <- counts[, 'Representative']
  
  p <- ggplot2::ggplot(data=counts, 
                       aes(x=counts[1:30, "Enriched1.Count"], 
                           y=factor(counts[1:30, "Representative"], levels = level_order), 
                           order=Enriched1.Count)) +
    geom_bar(stat="identity")+
    geom_text(aes(label=Enriched1.Count), hjust=-0.2, size=2.5) +
    labs(title=paste0("Top ", topn_thresh, " representative sequences after clustering \n", library_name), 
         x = 'Number of reads in the cluster', 
         y = 'Representative sequence') +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) #+
    #theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) 
  print(p)
}


filter.enriched.reduced.reps <- function(counts_file_path, output_path, enriched_chim_thresh, enriched_ratio_thresh, reduced_chim_thresh, reduced_ratio_thresh){
  counts <- read.csv(counts_file_path)
  
  enriched <- counts[counts$Chimeric.Percentage >= enriched_chim_thresh & counts$Ratio_En1Chim >= enriched_ratio_thresh,]
  enriched <- enriched[order(enriched$Ratio_En1Chim, decreasing = T),]
  rownames(enriched) <- NULL
  
  reduced <- counts[counts$Chimeric.Percentage >= reduced_chim_thresh & counts$Ratio_En1Chim < reduced_ratio_thresh,]
  #reduced <- na.omit(reduced) 
  reduced <- reduced[order(reduced$Ratio_En1Chim, decreasing = T),]
  rownames(reduced) <- NULL
  
  write.csv(enriched, file.path(output_path, "enriched_representatives.csv"), quote = F, row.names = F)
  write.csv(reduced, file.path(output_path, "reduced_representatives.csv"), quote = F, row.names = F)
  
  enriched_fa <- enriched[c("Representative", "Sequence")]
  # limit max 20 sequences
  enriched_fa <- enriched_fa[enriched_fa$Representative %in% enriched[1:min(nrow(enriched), 20),]$Representative,]
  colnames(enriched_fa) <- c("Header", "Sequence")
  microseq::writeFasta(enriched_fa, file.path(output.dir, "files/enriched_representatives.fasta"))
  
  reduced_fa <- reduced[c("Representative", "Sequence")]
  # limit max 20 sequences
  reduced_fa <- reduced_fa[reduced_fa$Representative %in% reduced[1:min(nrow(reduced), 20),]$Representative,]
  colnames(reduced_fa) <- c("Header", "Sequence")
  microseq::writeFasta(reduced_fa, file.path(output.dir, "files/reduced_representatives.fasta"))
  
  out <- list(enriched, reduced)
  return(out)
}

plot.enrichment.tiles <- function(enriched, reduced, type = "Enriched", axis_text_size = 6){
  if (type == "Enriched"){
    counts_filtered = enriched
    counts_filtered_other = reduced
  } else if (type == "Reduced"){
    counts_filtered = reduced
    counts_filtered_other = enriched
  } else {
    cat("\nPlease provide valid type: enriched or reduced.")
  }
  
  level_order <- counts_filtered[order(counts_filtered$Ratio_En1Chim, decreasing = F), 'Representative']
  mid = (round(min(min(counts_filtered$Chimeric.Percentage), min(counts_filtered_other$Chimeric.Percentage),
                   min(counts_filtered$Enriched1.Percentage), min(counts_filtered_other$Enriched1.Percentage)) +
                 max(max(counts_filtered$Chimeric.Percentage), max(counts_filtered_other$Chimeric.Percentage),
                     max(counts_filtered$Enriched1.Percentage), max(counts_filtered_other$Enriched1.Percentage)))+1)/2
  
  ggp1 <- ggplot(counts_filtered, aes(factor(Representative, level = level_order), "")) +    
    geom_tile(aes(fill = Chimeric.Percentage, height = 1)) +
    geom_text(aes(label = Chimeric.Percentage)) + 
    scale_fill_gradient2(low = "#1b98e0", mid = "white", high = "#EE4B2B", 
                         midpoint = mid,
                         name = "Chim[%]") +
    ylab("") +
    xlab("") +  
    ggtitle(paste0(type, " representative sequences after clustering \nRatio = Enriched 1/Chimeric")) + 
    theme_bw() +
    theme(axis.ticks.y = element_blank(), axis.text=element_text(size=axis_text_size), plot.title = element_text(hjust = 0.5)) +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 5.5))
  
  
  ggp2 <- ggplot(counts_filtered, aes(factor(Representative, level = level_order), "")) +    
    geom_tile(aes(fill = Enriched1.Percentage, height = 1)) +
    geom_text(aes(label = Enriched1.Percentage)) + 
    scale_fill_gradient2(low = "#1b98e0", mid = "white", high = "#EE4B2B", 
                         midpoint = mid,
                         name = "En1[%]") +
    ylab("") +
    xlab("") +  
    theme_bw() +
    theme(axis.ticks.y = element_blank(), axis.text=element_text(size=axis_text_size), plot.title = element_text(hjust = 0.5)) +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 5.5))
  
  
  ggp3 <- ggplot(counts_filtered, aes(factor(Representative, level = level_order), "")) +    
    geom_tile(aes(fill = Ratio_En1Chim, height = 1)) +
    geom_text(aes(label = Ratio_En1Chim)) + 
    scale_fill_gradient2(low = "#1b98e0", mid = "white", high = "#EE4B2B", 
                         midpoint = (min(min(counts_filtered$Ratio_En1Chim), min(counts_filtered_other$Ratio_En1Chim)) + 
                                       max(max(counts_filtered$Ratio_En1Chim), max(counts_filtered_other$Ratio_En1Chim)))/2, 
                         name = "Ratio") +
    ylab("") +
    xlab("") +  
    theme_bw() +
    theme(axis.ticks.y = element_blank(), axis.text=element_text(size=axis_text_size), plot.title = element_text(hjust = 0.5)) +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 5.5))
  
  

  print(cowplot::plot_grid(ggp1, ggp2, ggp3, align = "hv", ncol = 1))
  
}

# helper functions called in get.conserved.positions()
# function which finds the ranges given list of numbers
findRanges_helper <- function(run){
  rundiff <- c(1, diff(run))
  difflist <- split(run, cumsum(rundiff!=1))
  unlist(lapply(difflist, function(x){
    if(length(x) == 1) as.character(x) else paste0(x[1], "-", x[length(x)]) #length(x) == 1:2
  }), use.names=FALSE)
}

#duplicate single positions
dup_helper <- function(x) {
  if(length(x) == 1) {
    x = rep(x, 2)
  } else {
    x = x
  }
}

get.conserved.positions <- function(aln_file_path, format = "clustal"){
  alignment <- seqinr::read.alignment(aln_file_path, format = format, forceToLower = F)
  alignment$seq <- gsub("[\r\n\t]", "", alignment$seq)
  
  conn <- file(aln_file_path, open="r")
  linn <- readLines(conn)
  close(conn)
  
  #remove not needed lines
  linn <- linn[which(!(startsWith(linn, "CLUSTAL") | linn == ""))]
  
  #get identity pattern line
  identity_pattern = ""
  name_plus_6 <- max(nchar(alignment$nam)) + 6     #name with max length + 6 spaces
  for (i in seq(alignment$nb + 1, length(linn), alignment$nb + 1)){
    identity_pattern = paste0(identity_pattern, substring(linn[i], name_plus_6 + 1, nchar(linn[i])))
  }
  
  #find positions of * : . these symbols in identity_pattern string
  identity_symbol_idx  <- lapply(strsplit(identity_pattern, ''), 
                                 function(x) which(x == '*' | x == ":" | x == "."))
  identity_symbol_idx <- unlist(identity_symbol_idx)
  
  
  #get the ranges where there are identity symbols
  identity_ranges <- findRanges_helper(identity_symbol_idx)
  
  #remove single positions
  #identity_ranges = identity_ranges[unlist(lapply(identity_ranges, function(x) grepl("-", x)))]
  
  #split by - to separate start and end positions
  identity_ranges_splitted <- str_split(identity_ranges, pattern = "-")
  

  identity_ranges_splitted <- lapply(identity_ranges_splitted, dup_helper)
  
  #make dataframe of identity ranges
  identity_ranges_df = as.data.frame(matrix(unlist(identity_ranges_splitted), 
                                            ncol = 2, byrow = TRUE))
  colnames(identity_ranges_df) = c("start_nt", "end_nt")
  identity_ranges_df$start_nt = as.integer(identity_ranges_df$start_nt)
  identity_ranges_df$end_nt = as.integer(identity_ranges_df$end_nt)
  
  #filter to take only ranges with length >= 10
  #identity_ranges_df = identity_ranges_df[(identity_ranges_df$end_aa - identity_ranges_df$start_aa) >= 10,]
  
  #convert aa positions to nt
  # identity_ranges_df['start_nt'] <- identity_ranges_df$start_aa * 3 - 2
  # identity_ranges_df['end_nt'] <- identity_ranges_df$end_aa * 3 
  
  #for heatmap
  # identity_ranges_df['start_nt_'] <- identity_ranges_df$start_nt / 100 #read_length
  # identity_ranges_df['end_nt_'] <- identity_ranges_df$end_nt /100
  
  out = list(identity_ranges_df, alignment)
  return(out)
}


add.gap.info <- function(alignment, nj_matrix, step_size, nj_by_nt = F){
  aln_length <- nchar(alignment$seq[1])
  #split to characters
  aln_seq <- lapply(alignment$seq , function(y) unlist(strsplit(y, '')))
  aln_seq <- lapply(aln_seq, unlist)
  
  #gap positions in alignment sequences
  aln_gap_positions  <- lapply(aln_seq, 
                               function(x) which(x == '-'))
  
  aln_all_positions <- seq_len(aln_length) 
  
  aln_not_gap_positions <- lapply(aln_gap_positions, 
                                  function(x) setdiff(aln_all_positions, x))
  
  # adds gap positions for one variant
  fun <- function(gap_i, not_gap_i, name_i){
    l <- c()
    l[gap_i] <- 18
    
    if (nj_by_nt){
      for (i in seq_len(length(nj_matrix[name_i, ][nj_matrix[name_i, ] != 18]))) {
        ng_idx <- not_gap_i[i]
        l[ng_idx] <- nj_matrix[name_i, i]
      }
    } else {
      j = 0
      for (i in seq_len(length(nj_matrix[name_i, ][nj_matrix[name_i, ] != 18]))){ 
        ng_idx <- not_gap_i[(j + 1):min(length(not_gap_i), (j+step_size))]   #read_length for no overlap
        l[ng_idx] <- nj_matrix[name_i, i]
        j = j + step_size
      }
    }
    return(l)
  }
  
  list_all_nt <- mapply(fun, gap_i = aln_gap_positions, not_gap_i = aln_not_gap_positions, name_i = alignment$nam, SIMPLIFY = F)
  
  matrix_nt <- matrix(as.numeric(unlist(list_all_nt)), nrow = alignment$nb, byrow = TRUE)
  rownames(matrix_nt) <- alignment$nam
  
  return(matrix_nt)
}

#VP positions from https://www.ncbi.nlm.nih.gov/nuccore/110645916
plot.vp.positions <- function(){
  p <- ggplot() + 
        geom_segment(aes(x=1,xend=2208,y=2.5,yend=2.5), color = "blue", size=1.5) + #	2,203..4,410 Length:	2208 nt
        geom_segment(aes(x=412,xend=2208,y=2,yend=2), color = "green", size=1.5) + #2,614..4,410 Length:	1,797 nt
        geom_segment(aes(x=607,xend=2208,y=1.5,yend=1.5), color = "red", size=1.5) + #2,809..4,410 Length:	1,602 nt
        geom_hline(yintercept = 1, color = "white") + 
        geom_hline(yintercept = 3, color = "white") + 
        annotate("text", x = -61, y = 2.5, label = "VP1") +
        annotate("text", x = 352, y = 2, label = "VP2") +
        annotate("text", x = 547, y = 1.5, label = "VP3") +
        theme_classic() +
        xlab("") + ylab("") + 
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.line.y=element_blank()
        )
  print(p)
}


get.reps.nj.matrix.nt <- function(output_path, path_fasta, nj_matrix, nj_by_nt = F){
  data <- microseq::readFasta(path_fasta)
  ranges <- stringr::str_split(data$Header, pattern = "/", simplify = TRUE)[,3]  
  data['Start_orf'] <- stringr::str_split(ranges, pattern = "-", simplify = TRUE)[,1]
  data['End_orf'] <- stringr::str_split(ranges, pattern = "-", simplify = TRUE)[,2]
  data$Header <- sub("/ORF.*", "", data$Header)
  col_num <- max(nchar(data$Sequence))
  
  if(!nj_by_nt){
    fun <- function(name_i){
      positions <- seq_len(nchar(data[sub("/ORF.*", "", data$Header) == name_i,]$Sequence))
      l <- c()
      j = 0
      for (i in seq_len(length(nj_matrix[name_i, ][nj_matrix[name_i, ] != 18]))){
        idx <- positions[(j + 1):min(length(positions), (j+step_size))]
        l[idx] <- nj_matrix[name_i, i]
        j = j + step_size
      }
      l <- c(l, rep(18, col_num - length(l)))
      return(l)
    }
    
    list_all_nt <- lapply(rownames(nj_matrix), fun)
    reps_nj_mat_nt <- matrix(as.numeric(unlist(list_all_nt)), nrow = nrow(nj_matrix), byrow = TRUE)
    rownames(reps_nj_mat_nt) <- rownames(nj_matrix)
    
    seq_labels <- c()
    for (i in list_all_nt){
      i <- i[i != 18]
      seq_label <- paste(i, collapse = " ")
      seq_labels <- c(seq_labels, seq_label)
    }
  } else {
    seq_labels <- c()
    for (name_i in rownames(nj_matrix)){ 
      temp <- nj_matrix[name_i, ][nj_matrix[name_i, ] != 18]
      seq_label <- paste(temp, collapse = " ")
      seq_labels <- c(seq_labels, seq_label)
    }
  }
  
  chimeric_rep_predicted_labels <- data.frame(X = rownames(nj_matrix), 
                                              Composition = seq_labels)
  chimeric_rep_predicted_labels <- dplyr::inner_join(chimeric_rep_predicted_labels, 
                                                     data[c("Header", "Sequence", "Start_orf", "End_orf")], 
                                                     by = c("X" = "Header")) 
  
  write.csv(chimeric_rep_predicted_labels, file.path(output_path, "Chimeric_rep_predicted_labels.csv"), 
            row.names = F, quote = F)
  
}

