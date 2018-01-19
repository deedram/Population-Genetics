args = commandArgs(trailingOnly = TRUE)

fst <- function(path_in, path_out = path_in, input){
  fst_data <- read.table(path_in, header = TRUE, sep = ",")
  fst_data <- fst_data[,c("sample", "pos", "af")]
  
  # Split the data between the input sample and the other samples
  df_length <- nrow(fst_data)
  input_data <- fst_data[fst_data$sample == input,]
  other_samples <- fst_data[fst_data$sample != input,]

  # Modify fst_data to add af = 0 at every position not listed.
  unique_samples <- unique(unlist(other_samples$sample))
  for (i in 1:length(unique_samples)){
    data <- subset(other_samples, sample == unique_samples[i])
    data_adds <- input_data$pos[!input_data$pos%in%data$pos]
    new_data<- data.frame(sample = rep(data$sample[1], length(data_adds)), pos = data_adds, af = rep(0, length(data_adds)))
    other_samples <- rbind(other_samples, new_data)
    
    input_adds <- data$pos[!data$pos%in%input_data$pos]
    new_input <- data.frame(sample = rep(input_data$sample[1], length(input_adds)), pos = input_adds, af = rep(0, length(input_adds)))
    input_data <- rbind(input_data, new_input)
  }
  fst_data <- rbind(input_data, other_samples)

  # Compare the positions on the input data to the corresponding positions on each other sample
  ni = 1000
  nj = 1000

  datalist <- list()
  # d <- NULL
  for(i in 1:nrow(input_data)){
    for(j in 1:nrow(other_samples)){
      if (input_data$pos[i] == other_samples$pos[j]){
        if (!(input_data$af[i] == 0 & other_samples$af[j] == 0)){
          pi_hat = input_data$af[i]
          pj_hat = other_samples$af[j]
          p_hat = (pi_hat + pj_hat)/2
        
          alphai <- 2*pi_hat*(1-pi_hat)
          alphaj <- 2*pj_hat*(1-pj_hat)
        
          # Genetic variation within population
          b_s <- (ni*alphai + nj*alphaj)/(ni + nj - 1)
        
          # Genetic variation between populations
          a_s <- ((4*ni)*(pi_hat - p_hat)^2 + (4*nj)*(pj_hat - p_hat)^2 - b_s)/(2*((2*ni*nj)/(ni+nj)))
        
          # Fst for one site
          ab_s <- a_s + b_s
          fst <- a_s/ab_s
        }
        } else {next}
      datalist[[i]] <- list(as.character(input_data$sample[i]), as.character(other_samples$sample[j]), 
                              as.numeric(round(input_data$pos[i], 0)), as.numeric(round(fst, 5)))
    }
  }
  big_data <- as.data.frame(do.call(rbind, datalist))
  colnames(big_data) <- c("input", "sample", "pos", "fst")
  write.csv(big_data, out_path, colnames = TRUE)
}


