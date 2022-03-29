#Deeper test
library(OUwie)
library(tidyverse)
devtools::load_all()

#First simulate a tree and label it
tree <- sim.bdtree(n = 128, seed = 69)
tree$node.label <- c(rep.int(1, 64), rep.int(2, 63))
plot.phylo(tree, show.node.label = TRUE)

#Make dataframe
df_species <- tree$tip.label
df_regime <- c(rep.int(1, 64), rep.int(2, (64)))
df <- data.frame(df_species, df_regime) %>% rename(Genus_species = df_species, Reg = df_regime)

#Now sim
sim_dat <- OUwie.sim(phy = tree, data = df, alpha = c(1e-10,1e-10), sigma.sq = c(0.45,0.90), theta0 = 1.0, theta = c(0,0))

#Fit the data
sim_fit <- OUwie(phy = tree, data = sim_dat, model = "BMS")

sim_and_fit <- function(){
  tr <- sim.bdtree(n = 128)
  tr$node.label <- c(rep.int(1, 64), rep.int(2, 63))
  datasim <- OUwie.sim(phy = tr, data = df, alpha = c(1e-10,1e-10), sigma.sq = c(0.45,0.90), theta0 = 10.0, theta = c(0,0))
  datafit <- OUwie(phy = tr, data = datasim, model = "BMS")
  datafit$solution["sigma.sq",] <- c(0.45, 0.90)
  res <- arbutus(datafit)
  res
}

run <- replicate(1000, sim_and_fit())
saveRDS(run, file = "deepertest_result")

run_df <- run[1,] %>% map_df(function(x)x) 

run_df %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
  ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
ggsave("arbutus_plot.png")
