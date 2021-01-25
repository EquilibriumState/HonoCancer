# HonoCancer

Shiny app for basic analysis RNA-seq data derived from TCGA database
You can check it on https://www.shinyapps.io/admin/#/application/3317039

RudaTCGA.Rmd is a code which shows how data frames used in the app was created. Two of them:

- dataClin - has information about patients (where submitter_id, tumor_stage, prior_malignancy, vital_status and treatments_radiation_treatment_or_therapy are needed)
- cpm_fixed - cpm information of gene expression

This two data frames were loaded to shortdata.RData - it is example file, which is already loaded to app.

RData file input should look like shortdata.RData
