library(gplots)

test_HWE <- function(RES_i, epsilon, tol, method = 'complete', distance = 'euclidean')
{
  
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
  
  colMax_tol <- function(mat, tol)
  {
    sapply(1:ncol(mat), function(i) paste(names(which(mat[,i] > max(mat[,i]) - tol*abs(max(mat[,i])))), collapse = '_'))
  }
  
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

Visualize_JRC <- function(RES, i, RowSideColors = NULL)
{
  breaks=seq(-0.05, 1.05, 0.05)
  my_palette <- colorRampPalette(c("ghostwhite", "darkcyan"))(n = length(breaks) - 1)
  if(is.null(RowSideColors))
  {
    heatmap.2(t(RES[[i]]$M/(RES[[i]]$M + RES[[i]]$U)), trace = 'none', Colv = 'none', breaks = breaks,
              na.color = 'lightcoral', dendrogram = 'row', col = my_palette,
              offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
              main = i, cexCol = 0.8)
  }
  else
  {
    M = (RES[[i]]$M)[,order(RowSideColors)]
    U = (RES[[i]]$U)[,order(RowSideColors)]
    beta = M/(M+U)
    heatmap.2(t(beta), trace = 'none', Colv = 'none', Rowv = 'none', breaks = breaks,
              na.color = 'lightcoral', dendrogram = 'none', col = my_palette,
              offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, key.title = '', key.xlab = 'Methylation',
              RowSideColors = RowSideColors[order(RowSideColors)], main = i, cexCol = 0.8)
  }
  
}

RES = readRDS("M_U.Rds") # File at /test_run/results/
assignations = test_HWE(RES_i = RES[[10]], epsilon = 0.05, tol = 0.15)
Visualize_JRC(RES = RES, i = 10, RowSideColors = assignations)
