# evolution_TuMV_in_A.thaliana_Epigenetic_Mutants {.tabset}

Scripts used to analyse the RNAseq data of the project

- From preprocessing the raw fastqs to call the variants with `LoFreq`: `preprocess_rawData.sh`

## Session Info {-}

R version 4.2.1 (2022-06-23)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.5.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] es_ES.UTF-8/es_ES.UTF-8/es_ES.UTF-8/C/es_ES.UTF-8/es_ES.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods  
[7] base     

other attached packages:
[1] apex_1.0.4            phangorn_2.10.0      
[3] ape_5.6-2             vcfR_1.13.0          
[5] RVenn_1.1.0           clusterProfiler_4.4.4
[7] biomaRt_2.52.0       

loaded via a namespace (and not attached):
  [1] estimability_1.4.1          rappdirs_0.3.3             
  [3] coda_0.19-4                 tidyr_1.2.1                
  [5] ggplot2_3.4.0               bit64_4.0.5                
  [7] knitr_1.41                  multcomp_1.4-20            
  [9] DelayedArray_0.22.0         data.table_1.14.6          
 [11] KEGGREST_1.36.3             RCurl_1.98-1.9             
 [13] doParallel_1.0.17           generics_0.1.3             
 [15] BiocGenerics_0.42.0         TH.data_1.1-1              
 [17] RSQLite_2.2.20              shadowtext_0.1.2           
 [19] correlation_0.8.3           chron_2.3-58               
 [21] bit_4.0.5                   enrichplot_1.16.2          
 [23] xml2_1.3.3                  httpuv_1.6.8               
 [25] ggsci_2.9                   venn_1.11                  
 [27] SummarizedExperiment_1.26.1 assertthat_0.2.1           
 [29] viridis_0.6.2               STRINGdb_2.8.4             
 [31] xfun_0.36                   hms_1.1.2                  
 [33] evaluate_0.20               promises_1.2.0.1           
 [35] fansi_1.0.3                 progress_1.2.2             
 [37] caTools_1.18.2              dbplyr_2.3.0               
 [39] htmlwidgets_1.6.1           igraph_1.3.5               
 [41] DBI_1.1.3                   geneplotter_1.74.0         
 [43] stats4_4.2.1                paletteer_1.5.0            
 [45] purrr_1.0.1                 hash_2.2.6.2               
 [47] ellipsis_0.3.2              dplyr_1.0.10               
 [49] ggpubr_0.5.0                backports_1.4.1            
 [51] insight_0.18.8              permute_0.9-7              
 [53] annotate_1.74.0             MatrixGenerics_1.8.1       
 [55] vctrs_0.5.1                 Biobase_2.56.0             
 [57] abind_1.4-5                 cachem_1.0.6               
 [59] withr_2.5.0                 ggforce_0.4.1              
 [61] emmeans_1.8.4-1             vegan_2.6-4                
 [63] treeio_1.20.2               prettyunits_1.1.1          
 [65] cluster_2.1.4               DOSE_3.22.1                
 [67] pacman_0.5.1                lazyeval_0.2.2             
 [69] crayon_1.5.2                genefilter_1.78.0          
 [71] pkgconfig_2.0.3             tweenr_2.0.2               
 [73] GenomeInfoDb_1.32.4         nlme_3.1-161               
 [75] statsExpressions_1.4.0      rlang_1.0.6                
 [77] lifecycle_1.0.3             sandwich_3.0-2             
 [79] downloader_0.4              seqinr_4.2-23              
 [81] filelock_1.0.2              BiocFileCache_2.4.0        
 [83] adegenet_2.1.9              polyclip_1.10-4            
 [85] matrixStats_0.63.0          datawizard_0.6.5           
 [87] Matrix_1.5-3                aplot_0.1.9                
 [89] carData_3.0-5               zoo_1.8-11                 
 [91] GlobalOptions_0.1.2         png_0.1-8                  
 [93] viridisLite_0.4.1           rjson_0.2.21               
 [95] parameters_0.20.1           bitops_1.0-7               
 [97] KernSmooth_2.23-20          Biostrings_2.64.1          
 [99] blob_1.2.3                  shape_1.4.6                
[101] stringr_1.5.0               qvalue_2.28.0              
[103] rstatix_0.7.1               gridGraphics_0.5-1         
[105] S4Vectors_0.34.0            ggsignif_0.6.4             
[107] scales_1.2.1                memoise_2.0.1              
[109] magrittr_2.0.3              plyr_1.8.8                 
[111] gplots_3.1.3                zlibbioc_1.42.0            
[113] compiler_4.2.1              scatterpie_0.1.8           
[115] RColorBrewer_1.1-3          plotrix_3.8-2              
[117] clue_0.3-63                 DESeq2_1.36.0              
[119] cli_3.6.0                   ade4_1.7-20                
[121] XVector_0.36.0              patchwork_1.1.2            
[123] MASS_7.3-58.1               mgcv_1.8-41                
[125] tidyselect_1.2.0            stringi_1.7.12             
[127] nVennR_0.2.3                yaml_2.3.6                 
[129] GOSemSim_2.22.0             locfit_1.5-9.7             
[131] ggrepel_0.9.2               grid_4.2.1                 
[133] fastmatch_1.1-3             tools_4.2.1                
[135] parallel_4.2.1              circlize_0.4.15            
[137] rstudioapi_0.14             foreach_1.5.2              
[139] gridExtra_2.3               farver_2.1.1               
[141] ggraph_2.1.0                digest_0.6.31              
[143] shiny_1.7.4                 proto_1.0.0                
[145] quadprog_1.5-8              Rcpp_1.0.9                 
[147] GenomicRanges_1.48.0        car_3.1-1                  
[149] broom_1.0.2                 performance_0.10.2         
[151] later_1.3.0                 httr_1.4.4                 
[153] AnnotationDbi_1.58.0        ComplexHeatmap_2.12.1      
[155] colorspace_2.0-3            XML_3.99-0.13              
[157] IRanges_2.30.1              splines_4.2.1              
[159] yulab.utils_0.0.6           rematch2_2.1.2             
[161] tidytree_0.4.2              graphlayouts_0.8.4         
[163] ggplotify_0.1.0             xtable_1.8-4               
[165] jsonlite_1.8.4              ggtree_3.4.4               
[167] tidygraph_1.2.2             UpSetR_1.4.0               
[169] zeallot_0.1.0               ggfun_0.0.9                
[171] R6_2.5.1                    gsubfn_0.7                 
[173] pillar_1.8.1                htmltools_0.5.4            
[175] mime_0.12                   glue_1.6.2                 
[177] fastmap_1.1.0               pinfsc50_1.2.0             
[179] DT_0.27                     BiocParallel_1.30.4        
[181] codetools_0.2-18            fgsea_1.22.0               
[183] mvtnorm_1.1-3               utf8_1.2.2                 
[185] lattice_0.20-45             tibble_3.1.8               
[187] sqldf_0.4-11                curl_5.0.0                 
[189] gtools_3.9.4                GO.db_3.15.0               
[191] survival_3.5-0              admisc_0.30                
[193] rmarkdown_2.19              munsell_0.5.0              
[195] DO.db_2.9                   GetoptLong_1.0.5           
[197] GenomeInfoDbData_1.2.8      iterators_1.0.14           
[199] ggstatsplot_0.10.0          reshape2_1.4.4             
[201] gtable_0.3.1                bayestestR_0.13.0 
