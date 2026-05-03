## mes_nums.csv
Data from temperature rearing experiment.

Column names:
- f0_box: box ID for F0
- 
- temp: temp at which f1 box reared
- date_collected: date at which f1 egg pads collected
- sex: f1 individual sex
- pronotum: f1 pronotum length
- rwing_length: f1 right wing length
- rwing_width: f1 right wing width (if collected - males only)
- lwing_length: f1 left wing length(if collected)
- lwing_width: f1 right wing width (if collected)
- rep: concatenated form of f0_box, temp, date_collected
- pron.w.ratio: pronotum:wing length ratio
- sw.yn: whether individual was recorded as expressing small-wing (Sw) or typical long-wing (Lw) phenotypes
- rearing_density: number of f1 in box prior to culling (to approx 100 remaining f1)
- days_to_hatch: days from egg collection until first hatchling observe (if recorded)

## analysis_plot.R
R script used to analyse mes_nums.csv and plot data.

### R session info
R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.3

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] GenomicRanges_1.56.2 GenomeInfoDb_1.40.1  IRanges_2.38.1       S4Vectors_0.42.1     BiocGenerics_0.50.0  ggridges_0.5.6      
 [7] ggbeeswarm_0.7.2     lubridate_1.9.4      forcats_1.0.0        purrr_1.0.4          readr_2.1.5          tibble_3.2.1        
[13] tidyverse_2.0.0      tidyr_1.3.1          patchwork_1.3.0      dplyr_1.1.4          stringr_1.5.1        ggplot2_4.0.2       

loaded via a namespace (and not attached):
 [1] jsonlite_1.8.9          gtable_0.3.6            compiler_4.4.2          tidyselect_1.2.1        systemfonts_1.3.1      
 [6] scales_1.4.0            textshaping_1.0.0       XVector_0.44.0          R6_2.5.1                labeling_0.4.3         
[11] generics_0.1.3          GenomeInfoDbData_1.2.12 pillar_1.10.1           RColorBrewer_1.1-3      tzdb_0.4.0             
[16] rlang_1.1.5             stringi_1.8.4           S7_0.2.1                timechange_0.3.0        cli_3.6.3              
[21] withr_3.0.2             magrittr_2.0.3          zlibbioc_1.50.0         grid_4.4.2              rstudioapi_0.17.1      
[26] hms_1.1.3               beeswarm_0.4.0          lifecycle_1.0.4         vipor_0.4.7             vctrs_0.6.5            
[31] glue_1.8.0              farver_2.1.2            ragg_1.3.3              httr_1.4.7              UCSC.utils_1.0.0       
[36] tools_4.4.2             pkgconfig_2.0.3    
