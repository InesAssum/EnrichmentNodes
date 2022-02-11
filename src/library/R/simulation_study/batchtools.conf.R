# load cluster functions
cluster.functions <- makeClusterFunctionsSlurm(template="/home/icb/ines.assum/ext_tools/batchtools/template_slurm.tmpl",
                                               array.jobs=TRUE)
mail.start = "none"
mail.done = "none"
mail.error = "none"

# set some default resources
default.resources <- list(partition="cpu_p",
                          #exclude=" ",
                          memory="8G",
                          ncpus = 1,
                          measure.memory = TRUE,
                          walltime="48:00:00")
