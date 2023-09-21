if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version='3.16')
BiocManager::install("crisprScore")

