library(tidyverse)
library(cluster)

options(stringsAsFactors = F)

set.seed(2020)

knn <- as.numeric(snakemake@params[['knn']])
max_dist <- as.numeric(snakemake@params[['max_dist']])
k <- as.numeric(snakemake@params[['k']])

src_fls <- snakemake@input[['mixing']]
#src_fls <- Sys.glob("results/grid-search/larval-w1118-testes/overall-1/10-components/rep*/mixing.csv.gz")

usage_fls <- snakemake@input[['source']]
#usage_fls <- Sys.glob("results/grid-search/larval-w1118-testes/overall-1/10-components/rep*/source.csv.gz")

df <- src_fls %>%
  set_names(.,str_extract(.,"(?<=\\/rep-)\\d+(?=\\/mixing)")) %>%
  map_df(., read_csv,.id = "rep") %>%
  #dplyr::rename(index=X1) %>%
  gather(module,score,-rep,-index)

usage <- usage_fls %>%
  set_names(.,str_extract(.,"(?<=\\/rep-)\\d+(?=\\/source)")) %>%
  map_df(.,read_csv,.id='rep') %>%
  gather(module,usage,-X1,-rep)

aligner <- usage %>%
  group_by(rep,module) %>%
  summarize(aligner=-1*sign(median(usage))) %>%
  tidyr::unite(module, rep,module)

df <- pivot_wider(df,names_from = c(rep,module),values_from = score)

df <- gather(df,module,score,-index) %>%
  left_join(aligner, by='module') %>%
  mutate(score = score * aligner) %>%
  dplyr::select(-aligner) %>%
  spread(module,score)

mat <- t(column_to_rownames(df,'index'))

mat.dist <- as.matrix(dist(mat))

mat.dist.df <- mat.dist %>%
  as_tibble(rownames = 'comp1') %>%
  gather(comp2,dst,-comp1) %>%
  filter(comp1 != comp2)

mat.dist.df <- mat.dist.df %>%
  group_by(comp1) %>%
  top_n(knn,desc(dst)) %>%
  group_by(comp1) %>%
  summarize(outlier.score = mean(dst))

inliers <- filter(mat.dist.df ,outlier.score < max_dist)

mat2 <- mat[which(rownames(mat) %in% inliers$comp1),]

print(dim(mat))
print(mat.dist.df)
print(max_dist)
print(length(inliers))
print(inliers)
print(dim(mat2))

km <- kmeans(mat2, k, iter.max=500)

df2 <- gather(df, cluster,score,-index)

km_df <- km$cluster %>% enframe(name = 'cluster',value = 'cons')

sil <- cluster::silhouette(km$cluster,dmatrix=mat.dist[inliers$comp1,inliers$comp1])

sil_df <- as_tibble(sil[,1:3]) %>% mutate(component = names(km$cluster))

df2 <- df2 %>%
  left_join(km_df, by="cluster")

df2 <- df2 %>%
  filter(!is.na(cons)) %>%
  group_by(cons,index) %>%
  summarize(score=median(score))

# export
df2 %>% spread(cons,score) %>% write_csv(snakemake@output[['ica']])

# df2 %>% group_by(cons) %>%
#  mutate(l1norm = norm(as.matrix(unlist(score)),type = "1")) %>%
#  mutate(score = score/l1norm) %>%
#  dplyr::select(-l1norm) %>%
#  spread(cons,score) %>%
#	write_csv(snakemake@output[['components']])

usage2 <- unite(usage,cluster,rep,module) %>%
  left_join(km_df, by="cluster") %>%
  left_join(aligner, by=c(cluster='module')) %>%
  mutate(usage = usage*aligner) %>%
  group_by(cons,X1) %>%
  summarize(usage=median(usage))

left_join(mat.dist.df, km_df,by=c(comp1='cluster')) %>%
  dplyr::rename(cluster='comp1') %>%
  mutate(comps=k, rep = snakemake@wildcards[['ovr']]) %>%
  write_csv(snakemake@output[['dists']])

write_csv(usage2, snakemake@output[['usage']])

mutate(sil_df,comps=k, rep = snakemake@wildcards[['ovr']]) %>%
  write_csv(snakemake@output[['silhouette']])
