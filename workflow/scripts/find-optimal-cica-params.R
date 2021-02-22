library(tidyverse)
library(jsonlite)

metrics.grid <- read_csv(snakemake@input[[1]]) %>%
  mutate_at(vars('cov','pct_unique_term','pct_enr'), scales::rescale) %>%
  mutate(score = cov * pct_unique_term)

optimal <- metrics.grid %>%
  group_by(comps, qval) %>%
  summarize(score = mean(score)) %>%
  arrange(desc(score)) %>%
	head(1) %>%
	dplyr::select(qval,comps) %>%
  gather(metric,value) %>%
	deframe() %>%
	as.list()

write_json(optimal, snakemake@output[[1]])
