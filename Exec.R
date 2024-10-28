#Execute in order:
######
source("scripts/analyses/PCA.R") #Fig 1
source("scripts/analyses/RN_distrib_barplot.R") #Fig supp 1
######
#WGCNA network inference on the 3 sets of data:
source("scripts/networks/WGCNA_module_latitudinal.R")
source("scripts/networks/WGCNA_module_ancestors.R")
source("scripts/networks/WGCNA_module_latitudinal_and_ancestors.R")
######
#Network analyses Plots
source("scripts/networks/Module_sets_comparison.R")
source("scripts/networks/Modules_alluvial.R")
source("scripts/networks/Mod_GO.R")
source("scripts/networks/Evol_RN_mod.R")
######