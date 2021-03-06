# Run this after running hafoe to get accuracy and some plots
# Data simulation plots, accuracy

output.dir <- "hafoe_out"

##ACCURACY

# Takes predicted and true labels and calculates the accuracy
get.accuracy <- function(predicted_labels, chimeric_true_labels){
  accuracies <- c()
  for (i in predicted_labels$X) {
    pred <- as.numeric(unlist(stringr::str_split(predicted_labels[predicted_labels$X == i,]$Composition, " ")))
    pred <- pred[pred != 18]
    pred_start_orf <- predicted_labels[predicted_labels$X == i,]$Start_orf
    pred_end_orf <- predicted_labels[predicted_labels$X == i,]$End_orf
    
    true <- as.numeric(unlist(stringr::str_split(chimeric_true_labels[chimeric_true_labels$X == i,]$Composition, " ")))
    #true has labels for the whole sequence, we need only for orf
    true <- true[pred_start_orf:pred_end_orf]  
    
    acc <- round(sum(pred == true)*100/length(true), 2)
    accuracies <- c(accuracies, acc)
  }
  avg_accuracy <- round(sum(accuracies)/length(accuracies), 3)
  return(avg_accuracy)
}


chimeric_true_labels <- read.csv("input_files/Chimeric_lib_simulated_labels.csv")
predicted_labels <- read.csv(file.path(output.dir, "files/Chimeric_rep_predicted_labels.csv"))

get.accuracy(predicted_labels, chimeric_true_labels)


#analyse errors 
get.error.stats <- function(predicted_labels, chimeric_true_labels){
  no_alns <- c()
  multi_alns <- c()
  wrong_alns <- c()
  for (i in predicted_labels$X) {
    pred <- as.numeric(unlist(stringr::str_split(predicted_labels[predicted_labels$X == i,]$Composition, " ")))
    pred <- pred[pred != 18]
    pred_start_orf <- predicted_labels[predicted_labels$X == i,]$Start_orf
    pred_end_orf <- predicted_labels[predicted_labels$X == i,]$End_orf
    
    true <- as.numeric(unlist(stringr::str_split(chimeric_true_labels[chimeric_true_labels$X == i,]$Composition, " ")))
    #true has labels for the whole sequence, we need only for orf
    true <- true[pred_start_orf:pred_end_orf]  
    
    no_aln <- sum(pred == 0)*100/length(true)
    multi_aln <- sum(pred == 17)*100/length(true)
    wrong_aln <- sum(((pred != true)&(pred != 0)&(pred != 17)))*100/length(true)
    
    no_alns <- c(no_alns, no_aln)
    multi_alns <- c(multi_alns, multi_aln)
    wrong_alns <- c(wrong_alns, wrong_aln)
  }
  
  print(paste0("No alignments: ", round(mean(no_alns),3), "%"))
  print(paste0("Multiple alignments: ", round(mean(multi_alns), 3), "%"))
  print(paste0("Wrong alignments: ", round(mean(wrong_alns), 3), "%"))
}

get.error.stats(predicted_labels, chimeric_true_labels)

##PLOTS

#1. Variant description on chimeric library representatives (chosen by hafoe) true

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

s_all <- stringr::str_split(unlist(chimeric_true_labels[,'Composition'], 1), " ")

col_num_all <- max(unlist(lapply(s_all, length)))
for (i in seq_len(length(s_all))){
  if (length(s_all[[i]]) < col_num_all){
    s_all[[i]] <- c(s_all[[i]], rep("18", col_num_all - length(s_all[[i]]))) #gap
  }
}

matrix_all <- matrix(as.numeric(unlist(s_all)), ncol = col_num_all, byrow = TRUE)
rownames(matrix_all) <- chimeric_true_labels[, 'X']
representatives <- read.csv(file.path(output.dir, "files/clstr_chimeric_lib/representatives.csv"), header = F)
matrix_rep <- matrix_all[rownames(matrix_all) %in% representatives$V1,]


for (name_i in representatives$V1){
  pred <- as.numeric(unlist(stringr::str_split(predicted_labels[predicted_labels$X == name_i,]$Composition, " ")))
  pred <- pred[pred != 18]
  pred_start_orf <- predicted_labels[predicted_labels$X == name_i,]$Start_orf
  pred_end_orf <- predicted_labels[predicted_labels$X == name_i,]$End_orf
  matrix_rep[name_i,][1:length(pred)] <- matrix_rep[name_i,][pred_start_orf:pred_end_orf]
  matrix_rep[name_i,][(length(pred)+1):length(matrix_rep[name_i,])] <- 18
}


graphics.off()
pdf(file.path("data_simulation/plots/variant_description_chimeric_rep_true.pdf"), width=8, height=5)
plot.variant.description(matrix_rep, col_df = col_df,
                         library_name = "generated library \n(representatives chosen by the program)\n")
dev.off()

