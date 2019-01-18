test <- readRDS("/Users/ines/Documents/ICB/PhD/projects/KNIME_multiOMICs_enrichment/local_test/example_data/sim_data_84.RDS")
test2 <- test[["84"]][["0.1_0.1"]]
test2 <- test2[["mRNA"]]
saveRDS(test2,
        file="/Users/ines/Documents/ICB/PhD/projects/KNIME_multiOMICs_enrichment/local_test/example_data/mRNA.RDS")
