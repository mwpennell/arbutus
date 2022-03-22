#Deeper test
library(OUwie)
library(tidyverse)
devtools::load_all()

#First simulate a tree and label it
tree <- sim.bdtree(n = 128, seed = 69)
tree$node.label <- c(rep.int(1, 64), rep.int(2, (90-65)), rep.int(3, 38))
plot.phylo(tree, show.node.label = TRUE)

#Make dataframe
df_species <- tree$tip.label
df_regime <- c(rep.int(1, 64), rep.int(2, (90-65)), rep.int(3, 39))
df <- data.frame(df_species, df_regime) %>% rename(Genus_species = df_species, Reg = df_regime)

#Now sim
sim_dat <- OUwie.sim(phy = tree, data = df, alpha = c(1e-10,1e-10,1e-10), sigma.sq = c(0.45,0.90,0.23), theta0 = 1.0, theta = c(0,0,0))

#Fit the data
sim_fit <- OUwie(phy = tree, data = sim_dat, model = "BMS")

sim_and_fit <- function(){
  datasim <- OUwie.sim(phy = tree, data = df, alpha = c(1e-10,1e-10,1e-10), sigma.sq = c(0.45,0.90,0.23), theta0 = 1.0, theta = c(0,0,0))
  datafit <- OUwie(phy = tree, data = datasim, model = "BMS")
  res <- arbutus(datafit) %>% arbutus::pvalue_arbutus()
}

run <- replicate(1000, sim_and_fit)
saveRDS(run, file = "deepertest_result")