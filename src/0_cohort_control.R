library(data.table)
library(GenomicRanges)

coordinates2datatable = function(coordinates)
{
  chr = sapply(strsplit(coordinates, split = ':'), function(x) x[1])
  pos = sapply(strsplit(coordinates, split = ':'), function(x) x[2])
  start = as.integer(sapply(strsplit(pos, split = '-'), function(x) x[1]))
  end = as.integer(sapply(strsplit(pos, split = '-'), function(x) x[2]))
  data.table(chr = chr, start = start, end = end, range = end - start)
}

datatable2grange = function(mat)
{
  GRanges(mat$chr, IRanges(start = mat$start, end = mat$end))
}

# Argument:
# Correct way is to stratify gender by adjusting as a covariate (allows better estimation of variance)
# Thus, we employ Table S8 (discovery) and Table S11 (replication).
# Overlap  is 1 (min 200 bp) or 2 (min 0 bp)
# Additionally, authors filtered by permutation p-value, not by q-value (thus, p-vals are uncorrected for multiple testing)
# Thus: not strong association exists between TD/ASD; thus, why we include ASD samples.


# Get WB_JRCs
setwd("/media/ultron/2tb_disk1/ben/JRC_project/binokulars_output/WB_pool/binokulars_output/wb_results/")
IM_WB = fread('im_regions.txt', header = F)$V1
WB = fread('p_values.txt'); dim(WB) # 662659      2
WB = WB[!is.na(WB$V2),]
BT_WB = 0.05/nrow(WB) # 7.936269e-08
sum(WB$V2 < BT_WB) # 4215
JRCs = WB[WB$V2 < BT_WB,]$V1
JRCs_grange = datatable2grange(coordinates2datatable(JRCs))
length(JRCs_grange) # 4215

# Get DMRs
# setwd('/media/ultron/2tb_disk1/ben/JRC_project/annotations/ewas_asd/')
# males = fread('males_ASD.txt')[,1:11]
# males = datatable2grange(males)
# females = fread('females_ASD.txt')[,1:11]
# females = datatable2grange(females)
# length(males) # 634
# length(females) # 1940
combined = fread('combined_ASD.txt')
combined = datatable2grange(combined)
length(combined) # 186
replication = fread('combined_ASD_replication.txt')[,1:11]
replication = datatable2grange(replication)
length(replication) # 3755

findOverlaps(query = combined, subject = replication, minoverlap = 200)

combined[39]
replication[800]


# 
# findOverlaps(query = males, subject = females)
# males[140]; females[403]
# 
# m_common = as.matrix(findOverlaps(query = males, subject = JRCs_grange, minoverlap = 200))
# nrow(m_common) # 26
# males[4]; JRCs_grange[1961]
# # <Rle>         <IRanges>  <Rle>
# #   [1]     chr1 18690978-18691361      *
# # <Rle>         <IRanges>  <Rle>
# #   [1]     chr1 18678800-18706800      *
# f_common = as.matrix(findOverlaps(query = females, subject = JRCs_grange, minoverlap = 200))
# nrow(f_common) # 136
# overlaps = JRCs[c(m_common[,2], f_common[,2])]
# length(overlaps) # 162
# length(overlaps)/length(JRCs) # 4 %
# 
# 
# combined_common = as.matrix(findOverlaps(query = combined, subject = JRCs_grange, minoverlap = 200))
# overlaps2 = JRCs[combined_common[,2]]
# 
# replicated = as.matrix(findOverlaps(query = combined, subject = replication, minoverlap = 200))
# combined[39]; replication[800]
# # NO ISSUE; only 1 was replicated.
# 
# 
# # Compare overlaps with likelihood obtained
# setwd('/media/ultron/2tb_disk1/ben/JRC_project/JRC_sorter/results/')
# pred = fread('predictions.csv')
# JRC_seeker_ASD = pred[pred$rn %in% overlaps,]
# 
# table(JRC_seeker_ASD$predicted)
# # NA      L1  L3_L4_L5    L3_L5    L4_L5       L5 
# # 21      116        3        7        4        9 
# 
# potential_confounders = JRC_seeker_ASD[!(JRC_seeker_ASD$predicted %in% c('', 'L1')),]
# # rn         L1         L2         L3         L4         L5         pvals av_binfile predicted
# # 1:    chr21:8257168-8257397  -3657.484       -Inf  -3147.849       -Inf  -3110.441 1.798715e-196  122.91282     L3_L5
# # 2: chr9:122225800-122228600 -23862.736       -Inf -23785.107 -17755.267 -10989.726  1.870069e-39  120.25096        L5
# # 3: chr2:130037103-130037369  -3328.740       -Inf       -Inf  -3023.924  -2607.654  2.498147e-19   82.33590        L5
# # 4: chr2:130037672-130037900  -3770.910       -Inf       -Inf  -3211.593  -2547.579  5.812826e-13   74.85128        L5
# # 5:    chr17:4708946-4710200  -6971.378       -Inf  -6878.772  -5203.713  -2778.367  1.332395e-24   58.17500        L5
# # 6:      chr12:187600-191600 -12162.590       -Inf -12026.824  -8781.941  -8580.331  1.432069e-11   57.52727     L4_L5
# # 7: chr1:180951200-180955400 -11006.657       -Inf -10089.666 -10907.323  -9989.689  2.404070e-09   55.86014  L3_L4_L5
# # 8:  chr20:30891184-30893429 -10419.540       -Inf  -8037.898 -10388.367  -7819.592  1.404403e-38   55.65148     L3_L5
# # 9:  chr20:30895482-30896605  -3229.352       -Inf  -2444.761       -Inf  -2403.502  3.758774e-12   52.71923     L3_L5
# # 10:  chr22:50054600-50062800 -14887.376       -Inf -13433.156 -14810.619 -13266.060  3.165078e-10   50.80281     L3_L5
# # 11:      chr11:781200-784600  -8538.588       -Inf  -7778.323  -8377.567  -7601.240  9.451335e-09   49.52186     L3_L5
# # 12:   chr1:37994200-37998000 -10586.387 -10580.716 -10553.542  -9697.575  -9597.022  1.095972e-32   45.80586  L3_L4_L5
# # 13:  chr20:29315821-29321546 -12915.683       -Inf -12369.034 -12523.463 -10989.070  2.285442e-43   42.67625        L5
# # 14:   chrX:40140000-40153000 -25805.584       -Inf -25725.319 -20106.932 -19804.660  4.992942e-22   36.97308     L4_L5
# # 15:       chr4:119400-126800 -10324.060       -Inf  -9105.979 -10175.045  -8854.959  3.011992e-08   28.61175     L3_L5
# # 16:      chr12:629400-633600  -9917.571       -Inf  -9077.569  -9500.898  -5905.614  5.289079e-26   28.06573        L5
# # 17:  chr12:95544226-95548400  -3781.478       -Inf  -3677.133  -3677.583  -3286.035  1.394384e-10   26.85348        L5
# # 18:  chr10:62923600-62925600  -2794.261       -Inf       -Inf  -2108.663  -1517.028  9.437445e-10   25.05077        L5
# # 19:   chr4:55156200-55161400  -6007.788       -Inf       -Inf  -5748.567  -5179.358  2.295455e-08   21.76189        L5
# # 20:  chr13:27440000-27443400  -3589.190       -Inf  -2877.688  -3498.951  -2666.586  6.742457e-16   21.37343     L3_L5
# # 21:  chr15:81114797-81121600  -6395.192       -Inf  -6302.860  -6038.829  -5706.295  9.035221e-16   20.52981     L4_L5
# # 22:  chr14:61184800-61192200  -6112.730  -6104.702  -5959.591  -5790.632  -5548.984  2.985667e-08   19.21034  L3_L4_L5
# # 23:  chr15:23551200-23568000 -11229.928 -11216.829 -11212.135 -10473.360 -10001.254  3.195169e-16   18.82412     L4_L5
# # rn         L1         L2         L3         L4         L5         pvals av_binfile predicted
# 
# JRC_seeker_ASD2 = pred[pred$rn %in% overlaps2,]
# potential_confounders2 = JRC_seeker_ASD2[!(JRC_seeker_ASD2$predicted %in% c('', 'L1')),]
# # rn         L1   L2         L3         L4         L5        pvals av_binfile predicted
# # 1: chr20:29315821-29321546 -12915.683 -Inf -12369.034 -12523.463 -10989.070 2.285442e-43   42.67625        L5 #### DELICATE
# # 2: chr14:52696200-52699200  -3241.300 -Inf  -3127.669  -3058.010  -2911.301 2.950886e-09   23.46462  L3_L4_L5
# # 3:  chr4:55156200-55161400  -6007.788 -Inf       -Inf  -5748.567  -5179.358 2.295455e-08   21.76189        L5
# 
# 
# # 
# # potential_confounders
# # setwd('/media/ultron/2tb_disk1/ben/JRC_project/annotations/ewas_asd/')
# # males = fread('males_ASD.txt')
# # females = fread('females_ASD.txt')
# # males = males[m_common[,1],]
# # females = females[f_common[,1],]
# # 
# # males$pval
# # females$pval
# 
# 
# 
# setwd('/media/ultron/2tb_disk1/ben/JRC_project/annotations/ewas_asd/')
# pheno = fread('phenotype_PRJNA590702.txt')
# head(pheno)
# 
# table(pheno$characteristics_ch1.2, pheno$characteristics_ch1.3, pheno$sample_set)
# # , ,  = Discovery
# # Sex: F Sex: M
# # diagnosis: ASD     15     33
# # diagnosis: TD      15     35
# # , ,  = Replication
# # Sex: F Sex: M
# # diagnosis: ASD      5     15
# # diagnosis: TD       1     11
# 

##########################################################################

library(data.table)
library(gplots)

# Sex control
setwd('/media/ultron/2tb_disk1/ben/JRC_project/annotations/ewas_asd/')
pheno = as.data.frame(fread('phenotype_PRJNA590702.txt'))
rownames(pheno) = pheno$`Sample Name`

# Predicted sex
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

head(file.list)
ID = sapply(strsplit(file.list, split = '_'), function(x) x[2])

head(pheno)

table(sex, pheno[ID,]$characteristics_ch1.3)
# sex Sex: F Sex: M
# F     36      0
# M      0     94

balloonplot(as.table(t(table(sex, pheno[ID,]$characteristics_ch1.3))), xlab = 'Predicted', ylab = 'Reported', main = '')

#########################################################################################################

setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/WGBS_cordblood/')
file.list = list.dirs(); file.list = file.list[file.list != '.']
chrs = paste('chr', c(1:22, 'X', 'Y', 'M'), sep = '')


cov_mat = matrix(nrow = length(file.list), ncol = length(chrs))
rownames(cov_mat) = file.list; colnames(cov_mat) = chrs
setwd('/media/ultron/2tb_disk1/ben/JRC_project/bed_files/WGBS_cordblood/')
for(i in 1:length(file.list))
{
  print(file.list[i])
  setwd(file.list[i])
  for(j in 1:length(chrs))
  {
    print(chrs[j])
    res = tryCatch(fread(chrs[j], nThread = 4), error = function(x) NA)
    cov_mat[i,j] = tryCatch(sum(res$V4+res$V5+res$V11+res$V12), error = function(x) NA)
  }
  setwd('..')
}

X = log10(cov_mat)
rownames(X) = sapply(strsplit(rownames(X), split = '_'), function(x) x[2])
sex.col = factor(sex, levels = c('M', 'F')); levels(sex.col) = c('deepskyblue2', 'hotpink'); sex.col = as.character(sex.col)
breaks = seq(3.5, 7.5, 0.25/2)
my_palette <- colorRampPalette(c("ghostwhite", "darkcyan"))(n = length(breaks) - 1)
heatmap.2(t(X)[,order(sex.col)], Rowv = 'n', Colv = 'n', trace = 'n', dendrogram = 'n', ColSideColors = sex.col[order(sex.col)], 
          key.title = '', key.xlab = 'log10(sum(U+M))', breaks = breaks, col = my_palette, cexCol = 0.5)


# study.col = factor(pheno$study, levels = c('MARBLES', 'EARLI')); levels(study.col) = c('green2', 'skyblue2'); study.col = as.character(study.col)
# heatmap.2(t(X)[,order(study.col)], Rowv = 'n', Colv = T, trace = 'n', dendrogram = 'col', ColSideColors = study.col[order(study.col)], 
#           key.title = '', key.xlab = 'log10(sum(U+M))')


# study.col = factor(pheno$characteristics_ch1.2, levels = c('diagnosis: TD', 'diagnosis: ASD')); levels(study.col) = c('green2', 'skyblue2'); study.col = as.character(study.col)
# heatmap.2(t(X)[,order(study.col)], Rowv = 'n', Colv = 'n', trace = 'n', dendrogram = 'n', ColSideColors = study.col[order(study.col)], 
#           key.title = '', key.xlab = 'log10(sum(U+M))')
# head(pheno)





