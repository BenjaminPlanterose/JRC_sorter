
## Load libraries
library(gplots)
library(data.table)
library(parallel)
library(MASS)

## Load functions
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

test_region0 <- function(RES_i, epsilon, tol)
{
  M = RES_i[[1]]; U = RES_i[[2]]
  # Likelihoods
  L_mat = sapply(1:nrow(M), function(i) JRC_sorter(M[i,], U[i,], epsilon = epsilon))
  colnames(L_mat) = 1:ncol(L_mat)
  s_vec = rowSums(L_mat)
  
  # M_L = max(s_vec)
  # out_category = names(which(s_vec >= M_L - tol*abs(M_L)))
  # out_label = paste(out_category, collapse = '_')
  # HWE_pval = NA
  # if(out_label %in% c('L5', 'L2_L5', 'L3_L5', 'L4_L5'))
  # {
  #   rows = which(names(sapply(1:ncol(L_mat), function(i) which.max(L_mat[,i]))) %in% out_category)
  #   if(length(rows) == 1) # R collapses matrices of 1 row into vectors!
  #   {
  #     HWE_pval = test_HWE(M[rows, ], U[rows, ])
  #   }
  #   else
  #   {
  #     HWE_pval = test_HWE(colSums(M[rows, ]), colSums(U[rows, ]))
  #   }
  # }
  #return(list(out_label = out_label, likelihoods = s_vec, HWE_pval = HWE_pval))
  return(likelihoods = s_vec)
}



colMax_tol <- function(mat, tol)
{
  sapply(1:ncol(mat), function(i) paste(names(which(mat[,i] > max(mat[,i]) - tol*abs(max(mat[,i])))), collapse = '_'))
}

test_region <- function(RES_i, epsilon, tol)
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

classify <- function(s_vec, tol)
{
  names(s_vec) = paste('L', 1:length(s_vec), sep = '')
  M_L = max(s_vec)
  out_category = names(which(s_vec >= M_L - tol*abs(M_L)))
  out_label = paste(out_category, collapse = '_')
  return(out_label)
}

Visualize_JRC <- function(RES, i, sex.col = NULL, row.dend = T)
{
  if(!is.null(sex.col))
  {
    if(ncol(RES[[i]]$M) < length(sex.col))
    {
      print(likelihoods$rn[i])
      sex.col = rep('deepskyblue2', ncol(RES[[i]]$M))
    }
  }
  breaks=seq(-0.05, 1.05, 0.05)
  my_palette <- colorRampPalette(c("ghostwhite", "darkcyan"))(n = length(breaks) - 1)
  if(row.dend)
  {
    heatmap.2(t(RES[[i]]$M/(RES[[i]]$M + RES[[i]]$U)), trace = 'none', Colv = 'none', breaks = breaks,
              na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
              offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
              RowSideColors = sex.col, main = i, cexCol = 0.8)
  }
  else
  {
    M = (RES[[i]]$M)[,order(sex.col)]
    U = (RES[[i]]$U)[,order(sex.col)]
    beta = M/(M+U)
    heatmap.2(t(beta), trace = 'none', Colv = 'none', Rowv = 'none', breaks = breaks,
              na.color = 'lightcoral', dendrogram = 'none', col = my_palette,
              offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
              RowSideColors = sex.col[order(sex.col)], main = i, cexCol = 0.8)
  }
  
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

# Create JRC_WB.txt file
setwd('/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/WB_pool/binokulars_output/wb_results/')
pvals = fread('p_values.txt')
pvals = pvals[!is.na(pvals$V2),]
thr = 0.05/nrow(pvals)
pvals = pvals[order(pvals$V2),]
sig = pvals$V1[pvals$V2 < thr]
length(sig) # 4215
setwd('/media/ultron/2tb_disk1/ben/JRC_project/JRC_sorter/')
write.table(sig, 'JRC_WB.txt', quote = F, row.names = F, col.names = F)

# Extract U/M from cord blood files; following bash script:
# START=$(date +%s.%N)
# /home/ultron/Git/JRC_sorter/src/JRC_sorter -t /media/ultron/2tb_disk1/ben/JRC_project/JRC_sorter/JRC_WB.txt -i /media/ultron/2tb_disk1/ben/JRC_project/bed_files/WGBS_cordblood -l 200 -c 4 -f 200
# END=$(date +%s.%N)
# DIFF=$(echo "$END - $START" | bc)
# echo $DIFF # 3.83292 days

# Read results
setwd('/media/ultron/2tb_disk1/ben/JRC_project/JRC_sorter/results/')
RES = readRDS('M_U.Rds')
Names = fread('JRC_WB.txt', header = F)$V1
names(RES) = Names

# Likelihoods for all JRC
likelihoods = mclapply(1:length(RES), function(i) tryCatch(test_region(RES[[i]], epsilon = 0.05, tol = 0.15), 
                                                           error = function(x) NA), mc.cores = 4)
NA_pos = which(sapply(likelihoods, length) != 5); length(NA_pos) # 106
likelihoods[NA_pos] = rep(list(rep(NA, 5)), length = length(NA_pos))
likelihoods = Reduce(rbind, likelihoods)
rownames(likelihoods) = Names
setwd('/media/ultron/2tb_disk1/ben/JRC_project/JRC_sorter/results/')
fwrite(data.table(likelihoods, keep.rownames = T), 'likelihoods.csv')

#################################### Analysis ####################################

# Get sex
setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/WGBS_cordblood/')
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
sex.col = factor(sex, levels = c('M', 'F')); levels(sex.col) = c('deepskyblue2', 'hotpink'); sex.col = as.character(sex.col)

# Read results
setwd('/media/ultron/2tb_disk1/ben/JRC_project/JRC_sorter/results/')
RES = readRDS('M_U.Rds')
Names = fread('JRC_WB.txt', header = F)$V1
names(RES) = Names
likelihoods = fread('likelihoods.csv')
setwd('/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/WB_pool/binokulars_output/wb_results/')
pvals = fread('p_values.txt')
tmp = pvals$V1; pvals = pvals$V2; names(pvals) = tmp; rm(tmp)
pvals = pvals[Names]
likelihoods$pvals = pvals + .Machine$double.xmin
likelihoods$av_binfile = sapply(1:length(RES), function(i) tryCatch(mean(RES[[i]]$M + RES[[i]]$U), error = function(x) NA))

# Check dimensions
DIM = lapply(1:length(RES), function(x) dim(RES[[x]][[1]]))
DIM = Reduce(rbind, DIM)
hist(DIM[,1]) # negative exponential distribution
hist(DIM[,2]) # Y-chr markers only males
nrow(DIM); length(RES) # 4109; 4215; 106 NAs

# Check top JRC example
# RB1 expressed from the maternal allele (imprinted), linked to a CpG island at intron 2 (https://pubmed.ncbi.nlm.nih.gov/20551090/)
i = 1; Names[i]; test_region(RES[[i]], epsilon = 0.05, tol = 0.15); likelihoods[i,]
heatmap.2(t(RES[[i]]$M/(RES[[i]]$M + RES[[i]]$U)), trace = 'n', Colv = 'n', na.color = 'blue')

# What happened at missing values?
NA_pos = which(is.na(likelihoods$L1))
RES[NA_pos[3]]; RES[NA_pos[6]] # Basically, zero coverage.
# This could be computational bias (bisulfite aligner)
# Or technical bias (special protocol to obtain more homogeneous coverage)

# Sex?
Visualize_JRC(RES, 'chr20:30891184-30893429', sex.col) # sex?
Visualize_JRC(RES, 'chr20:30895482-30896605', sex.col) # sex?


############################### Export imgs ###############################

setwd('/media/ultron/2tb_disk1/ben/JRC_project/figures/')

# NINJ2
tiff(filename = 'NINJ2_HD1.tiff', width = 9, height = 5.5, units = 'in', res = 300)
assignations = test_HWE(RES$`chr12:629400-633600`, epsilon = 0.05, tol = 0.4)
Visualize_JRC(RES, 'chr12:629400-633600', assignations, row.dend = F) # NINJ2-AS1 (35th result)
dev.off()

tiff(filename = 'NINJ2_HD2.tiff', width = 9, height = 9, units = 'in', res = 300)
assignations = test_HWE(RES$`chr12:629400-633600`, epsilon = 0.05, tol = 0.4)
Visualize_JRC(RES, 'chr12:629400-633600', assignations, row.dend = F) # NINJ2-AS1 (35th result)
dev.off()

# VTRNA2-1
tiff(filename = 'VTRNA2_HD1.tiff', width = 9, height = 5.5, units = 'in', res = 300)
assignations = test_HWE(RES$`chr5:136067400-136083800`, epsilon = 0.05, tol = 0.15)
Visualize_JRC(RES, 'chr5:136067400-136083800', assignations, row.dend = F)
dev.off()

tiff(filename = 'VTRNA2_HD2.tiff', width = 9, height = 9, units = 'in', res = 300)
assignations = test_HWE(RES$`chr5:136067400-136083800`, epsilon = 0.05, tol = 0.15)
Visualize_JRC(RES, 'chr5:136067400-136083800', assignations, row.dend = F)
dev.off()

# WDR27
tiff(filename = 'WDR27_HD1.tiff', width = 9, height = 5.5, units = 'in', res = 300)
i = 'chr6:169653600-169656200'
assignations = test_HWE(RES[[i]], epsilon = 0.2, tol = 0.3) # WDR27
Visualize_JRC(RES, 'chr6:169653600-169656200', assignations, row.dend = F)
dev.off()

tiff(filename = 'WDR27_HD2.tiff', width = 9, height = 9, units = 'in', res = 300)
i = 'chr6:169653600-169656200'
assignations = test_HWE(RES[[i]], epsilon = 0.2, tol = 0.3) # WDR27
Visualize_JRC(RES, 'chr6:169653600-169656200', assignations, row.dend = F)
dev.off()


############################### Manually Curating a set of gold standards ###############################

# Problem: maybe not verifiable in cord blood dataset.

verified_mQTL0 = c('chr20:62950800-62960600', 'chr1:205849000-205851600', 'chr16:87645200-87655200', 
                  'chr2:208110200-208115600', 'chr19:13763400-13766400', 'chr17:6654200-6655993', 
                  'chr6:2950000-2953600', 'chr10:132228600-132233800', 'chr9:122225800-122228600', 
                  'chr12:629400-633600')
###### mQTL
# chr20:62950800-62960600
assignations = test_HWE(RES$`chr20:62950800-62960600`, epsilon = 0.2, tol = 0.98)
Visualize_JRC(RES, 'chr20:62950800-62960600', assignations, row.dend = F) # SLC17A9
# [1] "n = 130 ; p = 0.427"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 0.061935, df = 2, p-value = 0.9695
# 0 0.5   1 
# 42  65  23 
likelihoods[likelihoods$rn == 'chr20:62950800-62960600',]

# chr1:205849000-205851600
assignations = test_HWE(RES$`chr1:205849000-205851600`, epsilon = 0.05, tol = 0.15)
Visualize_JRC(RES, 'chr1:205849000-205851600', assignations, row.dend = F) # PM20D1
# [1] "n = 130 ; p = 0.346"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 10.654, df = 2, p-value = 0.004859
# 0 0.5   1 
# 64  42  24 
likelihoods[likelihoods$rn == 'chr1:205849000-205851600',]

# chr17:6654200-6655993
assignations = test_HWE(RES$`chr17:6654200-6655993`, epsilon = 0.05, tol = 0.15)
Visualize_JRC(RES, 'chr17:6654200-6655993', assignations, row.dend = F) #MIR4520-1
# [1] "n = 130 ; p = 0.515"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 0.034738, df = 2, p-value = 0.9828
# 0 0.5   1 
# 30  66  34 
likelihoods[likelihoods$rn == 'chr17:6654200-6655993',]

# chr12:629400-633600
assignations = test_HWE(RES$`chr12:629400-633600`, epsilon = 0.05, tol = 0.4)
Visualize_JRC(RES, 'chr12:629400-633600', assignations, row.dend = F) # NINJ2-AS1 (35th result)
# [1] "n = 130 ; p = 0.473"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 3.2126, df = 2, p-value = 0.2006
# 0 0.5   1 
# 31  75  24
likelihoods[likelihoods$rn == 'chr12:629400-633600',]

# chr9:122225800-122228600 (below)
# assignations = test_HWE(RES$`chr9:122225800-122228600`, epsilon = 0.05, tol = 0.15)
# Visualize_JRC(RES, 'chr9:122225800-122228600', assignations, row.dend = F) # LHK6
# likelihoods$predicted[likelihoods$rn == 'chr9:122225800-122228600'] # L5


###### Maybe not mQTL but polymorphic
# chr2:208110200-208115600
assignations = test_HWE(RES$`chr2:208110200-208115600`, epsilon = 0.05, tol = 0.15)
Visualize_JRC(RES, 'chr2:208110200-208115600', assignations, row.dend = F) # NR_038437
# [1] "n = 130 ; p = 0.708"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 22.179, df = 2, p-value = 1.527e-05
# 0 0.5   1 
# 0  76  54 
likelihoods$predicted[likelihoods$rn == 'chr2:208110200-208115600'] # L5


# chr6:2950000-2953600
assignations = test_HWE(RES$`chr6:2950000-2953600`, epsilon = 0.05, tol = 0.15)
Visualize_JRC(RES, 'chr6:2950000-2953600', assignations, row.dend = F) # SERPINB6
# [1] "n = 130 ; p = 0.35"
 # Chi-squared test for given probabilities
 # data:  obs_freq
# X-squared = 37.692, df = 2, p-value = 6.535e-09
# 0 0.5   1 
# 39  91   0 
likelihoods$predicted[likelihoods$rn == 'chr6:2950000-2953600'] # L3_L5



## Unclear
# chr19:13763400-13766400
assignations = test_HWE(RES$`chr19:13763400-13766400`, epsilon = 0.05, tol = 0.1)
Visualize_JRC(RES, 'chr19:13763400-13766400', assignations, row.dend = F) # MRI1
likelihoods$predicted[likelihoods$rn == 'chr19:13763400-13766400'] # L3_L5

# chr16:87645200-87655200
assignations = test_HWE(RES$`chr16:87645200-87655200`, epsilon = 0.15, tol = 0.8)
Visualize_JRC(RES, 'chr16:87645200-87655200', assignations, row.dend = F) # JPH3
likelihoods$predicted[likelihoods$rn == 'chr16:87645200-87655200'] # L3_L5
# chr10:132228600-132233800
assignations = test_HWE(RES$`chr10:132228600-132233800`, epsilon = 0.05, tol = 0.15)
Visualize_JRC(RES, 'chr10:132228600-132233800', assignations, row.dend = F) # STK32C
likelihoods$predicted[likelihoods$rn == 'chr10:132228600-132233800'] # L1
#

###### Polymorphic imprinting
ver_polymorf_imprinted0 = c('chr10:133523800-133534200', 'chr19:857357-864400', 'chr2:113235979-113236184',
                           'chr5:136067400-136083800', 'chr6:284200-300400', 'chr7:27129600-27135000')

# chr10:133523800-133534200, CYP2E1
assignations = test_HWE(RES$`chr10:133523800-133534200`, epsilon = 0.1, tol = 0.15)
Visualize_JRC(RES, 'chr10:133523800-133534200', assignations, row.dend = F)
# [1] "n = 130 ; p = 0.212"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 9.3575, df = 2, p-value = 0.009291
# 0 0.5   1 
# 75  55   0 
likelihoods$predicted[likelihoods$rn == 'chr10:133523800-133534200'] # L3_L5


# chr5:136067400-136083800, VTRNA2-1
assignations = test_HWE(RES$`chr5:136067400-136083800`, epsilon = 0.05, tol = 0.15)
Visualize_JRC(RES, 'chr5:136067400-136083800', assignations, row.dend = F)
# [1] "n = 130 ; p = 0.392"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 54.179, df = 2, p-value = 1.719e-12
# 0 0.5   1 
# 28 102   0 

# chr6:284200-300400, DUSP22
assignations = test_HWE(RES$`chr6:284200-300400`, epsilon = 0.05, tol = 0.14)
Visualize_JRC(RES, 'chr6:284200-300400', assignations, row.dend = F)
# [1] "n = 130 ; p = 0.381"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 49.154, df = 2, p-value = 2.12e-11
# 0 0.5   1 
# 31  99   0 
likelihoods$predicted[likelihoods$rn == 'chr6:284200-300400'] # L3_L5

# chr7:27129600-27135000, HOXA4
assignations = test_HWE(RES$`chr7:27129600-27135000`, epsilon = 0.19, tol = 0.49) # Its an mQTL but only one fully-meth; inapproapiate for dendrogram-based
Visualize_JRC(RES, 'chr7:27129600-27135000', assignations, row.dend = F)

jkl = RES$`chr7:27129600-27135000`$M/(RES$`chr7:27129600-27135000`$M + RES$`chr7:27129600-27135000`$U)

plot(jkl['7',], jkl['6',], ylim = c(0,1), xlim = c(0,1))

# [1] "n = 130 ; p = 0.227"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 11.201, df = 2, p-value = 0.003696
# 0 0.5   1 
# 71  59   0 
likelihoods$predicted[likelihoods$rn == 'chr7:27129600-27135000'] # L3_L4_L5

### UNCLEAR
# chr19:857357-864400
assignations = test_HWE(RES$`chr19:857357-864400`, epsilon = 0.05, tol = 0.1)
Visualize_JRC(RES, 'chr19:857357-864400', assignations, row.dend = F)
# chr2:113235979-113236184, PAX8
assignations = test_HWE(RES$`chr2:113235979-113236184`, epsilon = 0.05, tol = 0.15)
Visualize_JRC(RES, 'chr2:113235979-113236184', assignations, row.dend = F)


ver_imprinted = c("chr1:3655000-3670800", "chr1:3694200-3721400", "chr1:68046200-68048600", "chr1:68049200-68052600",
                  "chr11:1997000-2005624", "chr11:2131000-2137600", "chr11:2145400-2163000", "chr11:2696800-2706200",
                  "chr11:2777200-2783800", "chr11:2790800-2793000", "chr11:2867400-2874000", "chr11:3147200-3164600",
                  "chr11:32425200-32427400", "chr11:32427600-32433800", "chr11:70140000-70151800", "chr11:131684600-131700400",
                  "chr13:48317600-48324000", "chr14:100722000-100731600", "chr14:100822800-100828600", "chr15:23551200-23568000",
                  "chr15:23646600-23655000", "chr15:23683000-23697938", "chr15:24837800-24849400", "chr15:24953600-24959600",
                  "chr15:25858000-25862800", "chr15:25863400-25864800", "chr15:93024600-93038800", "chr16:798400-809000",
                  "chr16:3430000-3433000", "chr18:80157000-80162200", "chr19:56836400-56849400", "chr20:37518400-37522800",
                  "chr20:43512600-43516400", "chr20:58831400-58865400", "chr20:58887600-58891000", "chr4:88696400-88699083", 
                  "chr6:3846800-3852800", "chr6:106507800-106510400", "chr6:144006000-144009813", 
                  "chr7:50779200-50785800", "chr7:94655600-94658262", "chr7:94658297-94659400",
                  "chr7:130482800-130495400", "chr8:865800-871200", "chr8:891200-895400", "chr8:943800-950800", 
                  "chr8:1128200-1145200", "chr8:1157800-1167600", "chr8:1371000-1375280", "chr8:1548200-1549600", "chr8:1669600-1671600",
                  "chr8:1701000-1702600", "chr8:140094600-140102200")
length(ver_imprinted) # 53

confirmed = logical(length = length(ver_imprinted))
for(i in 1:length(ver_imprinted))
{
  message(ver_imprinted[i])
  tryCatch(Visualize_JRC(RES, ver_imprinted[i], sex.col, row.dend = T), error = function(x) plot(rnorm(1), rnorm(1)))
  confirmed[i] = as.logical(as.numeric(readline()))
}
sum(confirmed) # 16
ver_imprinted[confirmed]

# L1

length(ver_imprinted) # 16

assignations = test_HWE(RES$`chr1:68046200-68048600`, epsilon = 0.05, tol = 0.01)

i = ver_imprinted[1] # chr1:68046200-68048600; DIRAS3
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[3] # chr11:2131000-2137600; IGF2
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[4] # chr11:2696800-2706200; KCNQ1
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[5] # chr13:48317600-48324000; RB1
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[6] # chr14:100822800-100828600; MEG3
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[7] # chr15:24953600-24959600; SNRPN
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[8] # chr19:56836400-56849400; MIMT1
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[9] # chr20:43512600-43516400; L3MBTL1
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[10] # chr20:58887600-58891000; GNAS
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[11] # chr4:88696400-88699083; HERC3
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[12] # chr6:3846800-3852800; FAM50B
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[13] # chr6:144006000-144009813; PLAGL1
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[14] # chr7:50779200-50785800; GRB10
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[15] # chr7:130482800-130495400; MEST
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

i = ver_imprinted[16] # chr8:140094600-140102200; TRAPPC9
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1


# FOXK1 (cell type)
i = 'chr7:4712000-4718200'
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L1

# HOXA4 (ageing)
i = 'chr7:27129600-27135000'
Visualize_JRC(RES, i, assignations, row.dend = F)
likelihoods$predicted[likelihoods$rn == i] # L3_L4_L5



############################### Setting parameters ###############################

# VERIFIED mQTL, polymorphic imprinting and imprinting
verified_mQTL = c('chr20:62950800-62960600', 'chr1:205849000-205851600', 'chr17:6654200-6655993',
                  'chr12:629400-633600')
ver_polymorf_imprinted = c('chr10:133523800-133534200', 'chr5:136067400-136083800', 'chr6:284200-300400', 
                           'chr7:27129600-27135000')
ver_imprinted = c("chr1:68046200-68048600",    "chr1:68049200-68052600",    "chr11:2131000-2137600",     "chr11:2696800-2706200",    
                  "chr13:48317600-48324000",   "chr14:100822800-100828600", "chr15:24953600-24959600",   "chr19:56836400-56849400",  
                  "chr20:43512600-43516400",   "chr20:58887600-58891000",   "chr4:88696400-88699083",    "chr6:3846800-3852800" ,    
                  "chr6:144006000-144009813", "chr7:50779200-50785800",    "chr7:130482800-130495400",  "chr8:140094600-140102200" )
df = data.frame(rn = c(verified_mQTL, ver_imprinted, ver_polymorf_imprinted))
df = cbind(df, (likelihoods[match(df$rn, likelihoods$rn),-1]))
df$groundtruth = c(rep('L5', length(verified_mQTL)), rep('L1', length(ver_imprinted)), rep('L2orL4', length(ver_polymorf_imprinted)))
df$predicted = sapply(1:nrow(df), function(i) classify(as.numeric(df[i, 2:6]), tol = 0.1))
df$predicted[df$av_binfile < 18] = ''
df$predicted[grepl('L1', df$predicted)] = 'L1'
df$predicted[df$predicted %in% c('L3_L5', 'L4_L5')] = c('L3/L4_L5')

table(df$groundtruth, df$predicted)
#        L1 L3_L4_L5 L3/L4_L5 L5
# L1     16        0        0  0
# L2orL4  0        1        3  0
# L5      0        0        0  4

# df$rn[df$groundtruth == 'L2orL4' & df$predicted == 'L3_L4_L5'] # chr7:27129600-27135000
# assignations = test_HWE(RES$`chr7:27129600-27135000`, epsilon = 0.05, tol = 0.3)
# Visualize_JRC(RES, 'chr7:27129600-27135000', assignations, row.dend = F) 
# could be low_freq mQTL. Hclust does not allow a cluster of 1 sample


############################### Classification ###############################

# There are too tolerances: i) bin_tolerance: to run JRC_seeker; ii) classification_tolerance
likelihoods$predicted = sapply(1:nrow(likelihoods), function(i) classify(as.numeric(likelihoods[i, 2:6]), tol = 0.1))
likelihoods$predicted[likelihoods$av_binfile < 18] = ''
likelihoods$predicted[grepl('L1', likelihoods$predicted)] = 'L1'

table(likelihoods$predicted)
# NA        L1 L2_L4_L5   L3_L4_L5    L3_L5    L4_L5     L5
# 720     3286        1       26       70       56       56 

RES = RES[order(likelihoods$av_binfile, decreasing = T)]
likelihoods = likelihoods[order(likelihoods$av_binfile, decreasing = T),]

setwd('/media/ultron/2tb_disk1/ben/JRC_project/JRC_sorter/results/')
fwrite(data.table(likelihoods, keep.rownames = T), 'predictions.csv')



############################### Highlighting examples ###############################

which(likelihoods$rn == 'chr11:32425200-32427400')
which(likelihoods$rn == 'chr11:32427600-32433800')
Visualize_JRC(RES, likelihoods$rn[903], sex.col)
Visualize_JRC(RES, likelihoods$rn[939], sex.col)


######## L1
L1 = likelihoods[likelihoods$predicted == 'L1',] # full M, full U, boundary U/M, imprinting
Visualize_JRC(RES, L1$rn[10], sex.col)
Visualize_JRC(RES, L1$rn[13], sex.col)
Visualize_JRC(RES, L1$rn[14], sex.col)
Visualize_JRC(RES, L1$rn[15], sex.col)
Visualize_JRC(RES, L1$rn[45], sex.col)

 
# chr1:77163734-77341988
# samtools view -b output_38.bam "chr1:77163734-77341988" > random_WB.bam
# samtools index random_WB.bam


######## L3_L5
L3_L5 = likelihoods[likelihoods$predicted == 'L3_L5',]
head(L3_L5)
Visualize_JRC(RES, 'chr13:50126000-50130000', sex.col) # DLEU1
Visualize_JRC(RES, 'chr5:1593854-1596000', sex.col) # SDHAP3
Visualize_JRC(RES, 'chr7:158956800-158959600', sex.col) # WDR60
Visualize_JRC(RES, 'chr1:247936400-247939000', sex.col) # OR2L13

L3_L5[21:40, ]
######## L4_L5
L4_L5 = likelihoods[likelihoods$predicted == 'L4_L5',]
L4_L5[18:28,]
Visualize_JRC(RES, 'chr12:187600-191600', sex.col) # SLC6A12 (maybe mQTL?)
Visualize_JRC(RES, 'chr8:648400-654200', sex.col) # ERICH1
Visualize_JRC(RES, 'chr6:169653600-169656200', sex.col) # WDR27
Visualize_JRC(RES, 'chr9:97574200-97577800', sex.col) # TMOD1




## CFD, SVOTL, MRI1, 
Visualize_JRC(RES, 'chr7:138661400-138668600', sex.col) # 
likelihoods[likelihoods$rn == 'chr7:138661400-138668600',]




L4_L5[21:40,]
## X-inact
Visualize_JRC(RES, 'chrX:40140000-40153000', sex.col) # 
Visualize_JRC(RES, 'chrX:49071400-49075200', sex.col) # 
Visualize_JRC(RES, 'chrX:153031800-153034000', sex.col) # 


######## L5
L5 = likelihoods[likelihoods$predicted == 'L5',]
L5[31:40,]
Visualize_JRC(RES, 'chr9:122225800-122228600', sex.col) # LHX6
Visualize_JRC(RES, 'chr14:92686200-92687800', sex.col) # RIN3
Visualize_JRC(RES, 'chr8:687000-689200', sex.col) # ERICH1
Visualize_JRC(RES, 'chr5:77849200-77852400', sex.col) # NR_105012
Visualize_JRC(RES, 'chr5:77843400-77849000', sex.col) # NR_105012
Visualize_JRC(RES, 'chr1:147074800-147080400', sex.col) # RNVU1-8
Visualize_JRC(RES, 'chr19:12765000-12767200', sex.col) # HOOK2
Visualize_JRC(RES, 'chr17:6654200-6655993', sex.col) # MIR4520-1
Visualize_JRC(RES, 'chr17:4708946-4710200', sex.col) # ARRB2
Visualize_JRC(RES, 'chr5:180313000-180317200', sex.col) # GFPT2
Visualize_JRC(RES, 'chr20:62950800-62960600', sex.col) # SLC17A9
Visualize_JRC(RES, 'chr22:50033400-50037400', sex.col) # TTLL8
Visualize_JRC(RES, 'chr20:29315821-29321546', sex.col) # IG (ASD????)
Visualize_JRC(RES, 'chr1:205849000-205851600', sex.col) # PM20D1
Visualize_JRC(RES, 'chr2:9738000-9740128', sex.col) # IG
Visualize_JRC(RES, 'chr12:629400-633600', sex.col) # NINJ2-AS1 (35th result)
Visualize_JRC(RES, 'chr17:5768800-5771823', sex.col) # LOC339166
Visualize_JRC(RES, 'chr2:208110200-208115600', sex.col) # 

L5[35,]

#
rest = likelihoods[likelihoods$predicted %in% c('L2_L4_L5', 'L3_L4_L5'),]
rest[1:10,]
Visualize_JRC(RES, 'chr8:143928600-143930400', sex.col) # PLEC
Visualize_JRC(RES, 'chr1:180951200-180955400', sex.col) # 
Visualize_JRC(RES, 'chr7:27129600-27135000', sex.col) # HOXA4
Visualize_JRC(RES, 'chr1:37994200-37998000', sex.col) # FHL3
Visualize_JRC(RES, 'chr15:67060000-67065600', sex.col) # 
Visualize_JRC(RES, 'chr22:49479400-49483400', sex.col) # 
#


############################### Breaking ties ###############################


########################  L5
i = 'chr9:122225800-122228600'; i %in% verified_mQTL # F
assignations = test_HWE(RES[[i]], epsilon = 0.2, tol = 0.9) # LHX6
Visualize_JRC(RES, i, assignations, row.dend = F)
# [1] "n = 130 ; p = 0.673"
 # Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 2.6781, df = 2, p-value = 0.2621
# 0 0.5   1 
# 18  49  63 
likelihoods[likelihoods$rn == i,]


i = 'chr14:92686200-92687800'; i %in% verified_mQTL # F
assignations = test_HWE(RES[[i]], epsilon = 0.2, tol = 0.9) # RIN3
Visualize_JRC(RES, i, assignations, row.dend = F)
# [1] "n = 130 ; p = 0.308"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 23.175, df = 2, p-value = 9.28e-06
# 0 0.5   1 
# 74  32  24 
likelihoods[likelihoods$rn == i,]

i = 'chr5:77843400-77849000'; i %in% verified_mQTL # F
assignations = test_HWE(RES[[i]], epsilon = 0.2, tol = 0.73) # NR_105012
Visualize_JRC(RES, i, assignations, row.dend = F)
# [1] "n = 130 ; p = 0.542"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 4.8503, df = 2, p-value = 0.08846
# 0 0.5   1 
# 21  77  32 
likelihoods[likelihoods$rn == i,]

i = 'chr1:147074800-147080400'; i %in% verified_mQTL # F
assignations = test_HWE(RES[[i]], epsilon = 0.2, tol = 0.73) # RNVU1-8
Visualize_JRC(RES, i, assignations, row.dend = F)
# [1] "n = 130 ; p = 0.335"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 0.030622, df = 2, p-value = 0.9848
# 0 0.5   1 
# 58  57  15 
likelihoods[likelihoods$rn == i,]

i = 'chr19:12765000-12767200'; i %in% verified_mQTL # F
assignations = test_HWE(RES[[i]], epsilon = 0.2, tol = 0.73) # HOOK2
Visualize_JRC(RES, i, assignations, row.dend = F)
# [1] "n = 130 ; p = 0.338"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 4.0029, df = 2, p-value = 0.1351
# 0 0.5   1 
# 62  48  20 
likelihoods[likelihoods$rn == i,]

i = 'chr17:6654200-6655993'; i %in% verified_mQTL # T; MIR4520-1

i = 'chr17:4708946-4710200'; i %in% verified_mQTL # F
assignations = test_HWE(RES[[i]], epsilon = 0.2, tol = 0.73) # ARRB2
Visualize_JRC(RES, i, assignations, row.dend = F)
# [1] "n = 130 ; p = 0.735"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 2.9909, df = 2, p-value = 0.2241
# 0 0.5   1 
# 13  43  74 
likelihoods[likelihoods$rn == i,]

i = 'chr20:62950800-62960600'; i %in% verified_mQTL # T

i = 'chr22:50033400-50037400'; i %in% verified_mQTL # F
assignations = test_HWE(RES[[i]], epsilon = 0.19, tol = 0.83) # TTLL8
Visualize_JRC(RES, i, assignations, row.dend = F)
# [1] "n = 130 ; p = 0.473"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 2.9793, df = 2, p-value = 0.2255
# 0 0.5   1 
# 41  55  34 
likelihoods[likelihoods$rn == i,]

i = 'chr20:29315821-29321546'; i %in% verified_mQTL # F
assignations = test_HWE(RES[[i]], epsilon = 0.19, tol = 0.9) # IG
Visualize_JRC(RES, i, assignations, row.dend = F)
# [1] "n = 130 ; p = 0.4"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 0.43269, df = 2, p-value = 0.8055
# 0 0.5   1 
# 45  66  19 
likelihoods[likelihoods$rn == i,]

i = 'chr1:205849000-205851600'; i %in% verified_mQTL # T
i = 'chr12:629400-633600'; i %in% verified_mQTL # T


i = 'chr17:5768800-5771823'; i %in% verified_mQTL # F
assignations = test_HWE(RES[[i]], epsilon = 0.2, tol = 0.4) # LOC339166
Visualize_JRC(RES, i, assignations, row.dend = F)
# [1] "n = 130 ; p = 0.504"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 3.0749, df = 2, p-value = 0.2149
# 0 0.5   1 
# 37  55  38 
likelihoods[likelihoods$rn == i,]

i = 'chr2:208110200-208115600'; i %in% verified_mQTL # F
i = 'chr2:208110200-208115600'; i %in% verified_mQTL0 # T
assignations = test_HWE(RES[[i]], epsilon = 0.2, tol = 0.9) # NR_038437
Visualize_JRC(RES, i, assignations, row.dend = F) # 
likelihoods[likelihoods$rn == i,]


########################  L3
Visualize_JRC(RES, 'chr13:50126000-50130000', sex.col) # DLEU1
Visualize_JRC(RES, 'chr5:1593854-1596000', sex.col) # SDHAP3
Visualize_JRC(RES, 'chr7:158956800-158959600', sex.col) # WDR60
Visualize_JRC(RES, 'chr1:247936400-247939000', sex.col) # OR2L13

i = 'chr13:50126000-50130000'
assignations = test_HWE(RES[[i]], epsilon = 0.1, tol = 0.55) # DLEU1
Visualize_JRC(RES, i, assignations, row.dend = F) # 
# [1] "n = 130 ; p = 0.135"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 3.1457, df = 2, p-value = 0.2075
# 0 0.5   1 
# 95  35   0 
likelihoods$predicted[likelihoods$rn == i] # L3_L5

i = 'chr5:1593854-1596000'; i %in% ver_polymorf_imprinted # F
assignations = test_HWE(RES[[i]], epsilon = 0.14, tol = 0.2) # SDHAP3
Visualize_JRC(RES, i, assignations, row.dend = F) # 
# [1] "n = 130 ; p = 0.1"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 1.6049, df = 2, p-value = 0.4482
# 0 0.5   1 
# 104  26   0 
likelihoods$predicted[likelihoods$rn == i] # L3_L5

i = 'chr7:158956800-158959600'; i %in% ver_polymorf_imprinted # F
assignations = test_HWE(RES[[i]], epsilon = 0.14, tol = 0.2) # WDR60
Visualize_JRC(RES, i, assignations, row.dend = F) # 
# [1] "n = 130 ; p = 0.246"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 13.861, df = 2, p-value = 0.0009776
# 0 0.5   1 
# 66  64   0 
likelihoods$predicted[likelihoods$rn == i] # L3_L5

i = 'chr1:247936400-247939000'; i %in% ver_polymorf_imprinted # F
assignations = test_HWE(RES[[i]], epsilon = 0.14, tol = 0.2) # OR2L13
Visualize_JRC(RES, i, assignations, row.dend = F) # 
# [1] "n = 130 ; p = 0.096"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 1.4713, df = 2, p-value = 0.4792
# 0 0.5   1 
# 105  25   0 
likelihoods$predicted[likelihoods$rn == i] # L3_L5




########################  L4

i = 'chr12:187600-191600'; i %in% ver_polymorf_imprinted # F
assignations = test_HWE(RES[[i]], epsilon = 0.05, tol = 0.3) # SLC6A12; MAYBE L5
Visualize_JRC(RES, i, assignations, row.dend = F) # 
# [1] "n = 130 ; p = 0.777"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 10.718, df = 2, p-value = 0.004707
# 0 0.5   1 
# 0  58  72 

i = 'chr6:169653600-169656200'; i %in% ver_polymorf_imprinted # F
assignations = test_HWE(RES[[i]], epsilon = 0.2, tol = 0.3) # WDR27
Visualize_JRC(RES, i, assignations, row.dend = F) 
# [1] "n = 130 ; p = 0.542"
# Chi-squared test for given probabilities
# data:  obs_freq
# X-squared = 92.597, df = 2, p-value < 2.2e-16
# 0 0.5   1 
# 0 119  11 







