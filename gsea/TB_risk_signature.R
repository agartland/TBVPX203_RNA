if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("HGNChelper")
suppressPackageStartupMessages({
  library(TBSignatureProfiler)
  library(SummarizedExperiment)
})

library(TBSignatureProfiler)
library(HGNChelper)
data("TBsignatures")
names(TBsignatures)


tb_data<-SummarizedExperiment(assays=list(counts=counts))

tb_data <- mkAssay(tb_data, log = TRUE, counts_to_CPM = TRUE)

### Check to see that we now have 4 assays
assays(tb_data)


siglist_hivtb <- names(TBsignatures)
## Run the TBSignatureProfiler to score the signatures in the data
out <- capture.output(ssgsea_result <- runTBsigProfiler(input =tb_data,
                                                        useAssay = "log_counts_cpm",
                                                        signatures = TBsignatures,
                                                        algorithm = "ssGSEA",
                                                        combineSigAndAlgorithm = TRUE,
                                                        parallel.sz = 1))





## New colData entries from the Profiler
sigs <- c("Darboe_RISK_11", "Sweeney_OD_3", "Thompson_FAIL_13")
ssgsea_print_results <- as.data.frame(
  colData(ssgsea_result))[, c(sigs)]
ssgsea_print_results[, 1:3] <- round(ssgsea_print_results[, 1:3], 4)

DT::datatable(ssgsea_print_results)

write.csv(ssgsea_print_results,"Tb_risk_signature_results.csv")
