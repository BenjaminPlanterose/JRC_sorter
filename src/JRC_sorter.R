# Rscript --vanilla $MY_PATH/JRC_sorter.R 
# --target $target_file 
# --samples_dir $input_dir 
# --bin_length $bin_length 
# --output_dir $output_dir 
# --Lflank $Lflank 
# --nCores $nCores
# --epsilon $epsilon
# --tol $tol

# Read arguments
options = commandArgs(trailingOnly = TRUE)
index_names = startsWith(options, '--')
arguments = options[!index_names]
names(arguments) = options[index_names]
target_file = arguments['--target']
samples_dir = arguments['--samples_dir']
bin_length = as.numeric(arguments['--bin_length'])
nCores = as.numeric(arguments['--nCores'])
Lflank = as.numeric(arguments['--Lflank'])
epsilon = as.numeric(arguments['--epsilon'])
tol = as.numeric(arguments['--tol'])
res_folder = arguments['--output_dir']
res_dir = paste(getwd(), res_folder, sep = '/')


setwd(res_dir)
capture.output(
print(paste('target_file = ', target_file, '; samples_dir = ', samples_dir, '; bin_length = ', bin_length, 
            '; res_dir = ', res_dir, '; Lflank = ', Lflank, '; nCores = ', nCores, 
            '; epsilon = ', epsilon, '; tol = ', tol, sep = '')), file = 'params.txt')

# Load dependencies
library(data.table)
library(parallel)
library(MASS)

# Load functions
process_file <- function(file_i, START, END, breaks)
{
  #message(file_i)
  command = paste("awk ", "'", "$2>=", START, " && ", "$2<", END, "' ", file_i, sep = "")
  sub_data = fread(cmd = command, sep = "\t", header = F)
  sub_data = data.frame(start = sub_data$V2, M = sub_data$V4 + sub_data$V11, U = sub_data$V5 + sub_data$V12)
  sub_data$bins = cut(sub_data$start, breaks = breaks, labels = F)
  X = Reduce(f = cbind, list(aggregate(start ~ bins, data = sub_data, FUN = mean)[,-1],
                             aggregate(M ~ bins, data = sub_data, FUN = sum), 
                             aggregate(U ~ bins, data = sub_data, FUN = sum)[,-1]))
  colnames(X) = c('pos', 'bins', 'M', 'U')
  return(X)
}

rowMin <- function(mat)
{
  sapply(1:nrow(mat), function(x) which.min(mat[x,]))
}

JRC_sorter <- function(U, M, epsilon)
{
  total = U+M
  
  # Model 1
  p1 = sum(M)/sum(total)
  L1 = sum(dbinom(x = M, size = total, prob = p1, log = T))
  
  # Model 2
  D_mat = cbind((M-0)^2, (M-total)^2)
  cluster2 = rowMin(D_mat)-1
  indicator2 = cluster2 + (cluster2 == 0)*epsilon - (cluster2 == 1)*epsilon
  L2 = sum(dbinom(x = M, size = total, prob = indicator2, log = T))
  
  # Model 3
  D_mat = cbind((M-0.5*total)^2, (M-total)^2)
  cluster3 = (rowMin(D_mat))/2
  indicator3 = cluster3 + (cluster3 == 0)*epsilon - (cluster3 == 1)*epsilon
  L3 = sum(dbinom(x = M, size = total, prob = indicator3, log = T))
  
  # Model 4
  D_mat = cbind((M-0)^2, (M-0.5*total)^2)
  cluster4 = (rowMin(D_mat)-1)/2
  indicator4 = cluster4 + (cluster4 == 0)*epsilon - (cluster4 == 1)*epsilon
  L4 = sum(dbinom(x = M, size = total, prob = indicator4, log = T))
  
  # Model 5
  D_mat = Reduce(cbind, list((M-0)^2, (M-0.5*total)^2, (M-total)^2))
  cluster5 = (rowMin(D_mat)-1)/2
  #p5 = mean(cluster5)
  indicator5 = cluster5 + (cluster5 == 0)*epsilon - (cluster5 == 1)*epsilon
  L5 = sum(dbinom(x = M, size = total, prob = indicator5, log = T))
  
  L_vec = c(L1 = L1, L2 = L2, L3 = L3, L4 = L4, L5 = L5)
  L_vec
}

test_region <- function(files, START, END, bin_length, epsilon, tol)
{
  breaks = seq(START-1, END, bin_length)
  X_list = lapply(files, function(x) process_file(x, START, END, breaks))
  
  M = matrix(0, nrow = length(breaks), ncol = length(files))
  U = matrix(0, nrow = length(breaks), ncol = length(files))
  rownames(M) = rownames(U) = 1:length(breaks)
  colnames(M) = colnames(U) = sapply(strsplit(files, split = '_'), function(x) x[1])
  for(i in 1:length(X_list))
  {
    M[X_list[[i]]$bins, i] = X_list[[i]]$M
    U[X_list[[i]]$bins, i] = X_list[[i]]$U
  }
  
  # Likelihoods
  L_mat = sapply(1:nrow(M), function(i) JRC_sorter(M[i,], U[i,], epsilon = epsilon))
  colnames(L_mat) = 1:length(breaks)
  L_mat = t(na.omit(t(L_mat)))
  s_vec = rowSums(L_mat)
  out_category = names(which(s_vec >= max(s_vec)-tol))
  out_label = paste(out_category, collapse = '_')
  HWE_pval = NA
  if(out_label %in% c('L5', 'L2_L5', 'L3_L5', 'L4_L5'))
  {
    rows = colnames(L_mat)[names(sapply(1:ncol(L_mat), function(i) which.max(L_mat[,i]))) %in% out_category]
    HWE_pval = test_HWE(colSums(M[rows, ]), colSums(U[rows, ]))
  }
  return(list(out_label = out_label, HWE_pval = HWE_pval))
}

test_HWE <- function(M, U)
{
  # The Pearson goodness of fit test for HWE
  total = M+U
  D_mat = Reduce(cbind, list((M-0)^2, (M-0.5*total)^2, (M-total)^2))
  cluster5 = (rowMin(D_mat)-1)/2
  cluster5 = table(factor(cluster5, levels = c('0', '0.5', '1')))
  n = sum(cluster5)
  freq_obs = cluster5/sum(cluster5)
  p = unname(freq_obs['1'] + freq_obs['0.5']/2)
  freq_exp = c(`0` = (1-p)^2, `0.5` = 2*p*(1-p), `1` = p^2)
  print(freq_obs)
  print(freq_exp)
  chisq.test(x = cluster5, p = freq_exp, correct = T)$p.value
}

# Read target_file
regions = fread(input = target_file, header = F)$V1
Chr = sapply(strsplit(regions, split = ':'), function(x) x[1])
pos = sapply(strsplit(regions, split = ':'), function(x) x[2])
Start = as.integer(sapply(strsplit(pos, split = '-'), function(x) x[1])) - Lflank # Add flanks
End = as.integer(sapply(strsplit(pos, split = '-'), function(x) x[2])) + Lflank # Add flanks

# Prep names of sample files per CHR
setwd(samples_dir)
file.list = list.dirs(); file.list = file.list[file.list != '.']
file.list = lapply(Chr, function(x) paste(file.list, x, sep = '/'))

# Run JRC_sorter
RES = mclapply(1:length(regions), function(i) tryCatch(test_region(file.list[[i]], Start[i], End[i], 
                                                                        bin_length, epsilon, tol), error = function(x) NA), 
               mc.cores = nCores)

out_labels = unlist(sapply(RES, function(x) x[1]))
names(out_labels) = regions
pval_HWE = unlist(sapply(RES, function(x) x[2]))
names(pval_HWE) = regions

# Export results
setwd(res_dir)
write.table(x = data.frame(V1 = pval_HWE), file = 'p_values_HWE.txt', quote = F, col.names = F, sep = '\t')
write.table(x = data.frame(V1 = out_labels), file = 'out_labels.txt', quote = F, col.names = F, sep = '\t')