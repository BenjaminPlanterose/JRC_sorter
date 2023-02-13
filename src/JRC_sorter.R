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
bin_tol = as.numeric(arguments['--bin_tol'])
mod_tol = as.numeric(arguments['--mod_tol'])
res_folder = arguments['--output_dir']
av_count_bin_sample = as.numeric(arguments['--av_count_bin_sample'])
res_dir = paste(getwd(), res_folder, sep = '/')

setwd(res_dir)
capture.output(
print(paste('target_file = ', target_file, '; samples_dir = ', samples_dir, '; bin_length = ', bin_length, 
            '; res_dir = ', res_dir, '; Lflank = ', Lflank, '; nCores = ', nCores, 
            '; epsilon = ', epsilon, '; bin_tol = ', bin_tol,  '; mod_tol = ', mod_tol, 
            '; av_count_bin_sample = ', av_count_bin_sample, sep = '')), file = 'params.txt')

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
  
  if(nrow(sub_data) > 1)
  {
    sub_data = data.frame(start = sub_data$V2, M = sub_data$V4 + sub_data$V11, U = sub_data$V5 + sub_data$V12)
    sub_data$bins = cut(sub_data$start, breaks = breaks, labels = F)
    X = Reduce(f = cbind, list(aggregate(start ~ bins, data = sub_data, FUN = mean)[,-1],
                               aggregate(M ~ bins, data = sub_data, FUN = sum), 
                               aggregate(U ~ bins, data = sub_data, FUN = sum)[,-1]))
    colnames(X) = c('pos', 'bins', 'M', 'U')
    return(X)
  }
  else
  {
    return(NA)
  }
}

rowMin <- function(mat)
{
  sapply(1:nrow(mat), function(x) which.min(mat[x,]))
}

extract_UM <- function(files, START, END, bin_length)
{
  breaks = seq(START-1, END, bin_length)
  X_list = lapply(files, function(x) process_file(x, START, END, breaks))
  remove = which(sapply(X_list, function(x) mean(is.na(x))) == 1)
  X_list = X_list[!(1:length(X_list) %in% remove)]
  M = matrix(0, nrow = length(breaks), ncol = length(files))
  U = matrix(0, nrow = length(breaks), ncol = length(files))
  rownames(M) = rownames(U) = 1:length(breaks)
  colnames(M) = colnames(U) = sapply(strsplit(files, split = '_'), function(x) x[1])
  for(i in 1:length(X_list))
  {
    M[X_list[[i]]$bins, i] = X_list[[i]]$M
    U[X_list[[i]]$bins, i] = X_list[[i]]$U
  }
  # Remove empty rows
  remove = which(rowSums(M+U, na.rm = T) == 0)
  M = M[-remove,]; U = U[-remove,]
  
  return(list(M = M, U = U))
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

compute_likelihoods <- function(RES_i, epsilon, tol)
{
  M = RES_i[[1]]; U = RES_i[[2]]
  # Likelihoods
  L_mat = sapply(1:nrow(M), function(i) JRC_sorter(M[i,], U[i,], epsilon = epsilon))
  colnames(L_mat) = 1:ncol(L_mat)
  state_bin = colMax_tol(L_mat, tol)
  state_bin[grepl('L1', state_bin)] = 'L1'
  
  s_vec = numeric(length = 5); names(s_vec) = paste('L', 1:5, sep='')
  s_vec['L1'] = sum(L_mat[1,])
  
  for(i in paste('L', 2:5, sep=''))
  {
    where = grepl(i, state_bin)
    if(sum(where) > 0)
    {
      values = L_mat[1,]
      values[where] = L_mat[i, where]
      s_vec[i] = sum(values)
    }
    else
    {
      s_vec[i] = -Inf
    }
  }
  #s_vec = rowSums(L_mat)
  return(likelihoods = s_vec)
}

colMax_tol <- function(mat, tol)
{
  sapply(1:ncol(mat), function(i) paste(names(which(mat[,i] > max(mat[,i]) - tol*abs(max(mat[,i])))), collapse = '_'))
}

classify <- function(s_vec, tol)
{
  names(s_vec) = paste('L', 1:length(s_vec), sep = '')
  M_L = max(s_vec)
  out_category = names(which(s_vec >= M_L - tol*abs(M_L)))
  out_label = paste(out_category, collapse = '_')
  return(out_label)
}

test_HWE <- function(RES_i, epsilon, tol, method = 'complete', distance = 'euclidean')
{
  # Extract M and U
  M = RES_i$M
  U = RES_i$U
  
  # Re-compute likelihoods
  L_mat = sapply(1:nrow(M), function(i) JRC_sorter(M[i,], U[i,], epsilon = epsilon))
  colnames(L_mat) = 1:ncol(L_mat)
  state_bin = colMax_tol(L_mat, tol)
  state_bin[grepl('L1', state_bin)] = 'L1'
  
  # Extract contested states
  rows = which(grepl('L5', state_bin))
  M = M[rows, ]
  U = U[rows, ]
  # total = M+U
  #
  beta = na.omit(M/(M+U))
  mean_beta = colMeans(beta)
  dd = dist(t(beta), method = distance)
  hc = hclust(dd, method = method)
  assignations = cutree(hc, k = 3)
  summary = aggregate(mean_beta ~ assignations, mean, data = data.frame(mean_beta, assignations))
  summary$cluster = (rowMin(sapply(c(0, 0.5, 1), function(x) (summary$mean_beta - x)^2))-1)/2
  assignations = factor(assignations, levels = summary$assignations)
  levels(assignations) = summary$cluster
  assignations = as.character(assignations)
  
  # The Pearson goodness of fit test for HWE
  obs_freq = table(factor(assignations, levels = c('0', '0.5', '1')))
  n = sum(obs_freq)
  freq_obs = obs_freq/sum(obs_freq)
  p = unname(freq_obs['1'] + freq_obs['0.5']/2)
  freq_exp = c(`0` = (1-p)^2, `0.5` = 2*p*(1-p), `1` = p^2)
  #chisq.test(x = cluster5, p = freq_exp, correct = T)$p.value
  print(paste('n =', n, '; p =', round(p, 3)), sep = '')
  print(chisq.test(x = obs_freq, p = freq_exp, correct = T))
  print(n*freq_obs)
  
  # Colour assignations
  col = factor(assignations, levels = c('0', '0.5', '1'))
  levels(col) = scales::alpha(c('dodgerblue3', 'darkmagenta', 'brown3'), 0.5)
  
  return(as.character(col))
}


# Predict sex
setwd(samples_dir)
file.list = list.dirs(); file.list = file.list[file.list != '.']
file.list = paste(file.list, 'chrY', sep = '/')
sizes.Y = sapply(file.list, function(x) tryCatch(file.info(x)$size, error = function(x) 0))
file.list = list.dirs(); file.list = file.list[file.list != '.']
file.list = paste(file.list, 'chrX', sep = '/')
sizes.X = sapply(file.list, function(x) tryCatch(file.info(x)$size, error = function(x) 0))
set.seed(1); kmeans_res = kmeans(x = sizes.Y/(sizes.Y + sizes.X + 1), centers = 2)$cluster
size_tab = aggregate(sizes.Y ~ factor(kmeans_res, levels = c(1,2)), FUN = mean)
colnames(size_tab) = c('levels', 'size')
male_label = as.character(size_tab$levels)[which.max(size_tab$size)]
sex = kmeans_res
sex[kmeans_res == male_label] = 'M'
sex[kmeans_res != male_label] = 'F'

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

# Only male samples for chrY
chry_pos = Chr %in% c('chrY')
if(sum(chry_pos) > 0)
{
  file.list[chry_pos] = lapply(file.list[chry_pos], function(x) x[sex == 'M'])
}

# Extract and bin U/M counts for regions across samples.
RES = mclapply(1:length(regions), function(i) tryCatch(extract_UM(file.list[[i]], Start[i], End[i], 
                                                                        bin_length), error = function(x) NA), 
               mc.cores = nCores)
names(RES) = regions
setwd(res_dir)
saveRDS(object = RES, file = 'M_U.Rds')


# Compute Model likelihoods.
likelihoods = mclapply(1:length(RES), function(i) tryCatch(compute_likelihoods(RES[[i]], epsilon = epsilon, tol = bin_tol), 
                                                           error = function(x) NA), mc.cores = nCores)
NA_pos = which(sapply(likelihoods, length) != 5); length(NA_pos) # 106
likelihoods[NA_pos] = rep(list(rep(NA, 5)), length = length(NA_pos))
likelihoods = Reduce(rbind, likelihoods)
rownames(likelihoods) = regions
likelihoods = as.data.frame(likelihoods)
likelihoods$av_binfile = sapply(1:length(RES), function(i) tryCatch(mean(RES[[i]]$M + RES[[i]]$U), error = function(x) NA))
setwd(res_dir)
fwrite(data.table(likelihoods, keep.rownames = T), 'likelihoods.csv', nThread = nCores)

# Predict model that best fits the data.
likelihoods$predicted = sapply(1:nrow(likelihoods), function(i) classify(as.numeric(likelihoods[i, 1:5]), tol = mod_tol))
likelihoods$predicted[likelihoods$av_binfile < av_count_bin_sample] = ''
likelihoods$predicted[grepl('L1', likelihoods$predicted)] = 'L1'
setwd(res_dir)
fwrite(data.table(likelihoods, keep.rownames = T), 'predictions.csv', nThread = nCores)

