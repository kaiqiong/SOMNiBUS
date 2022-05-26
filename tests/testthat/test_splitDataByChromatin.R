context("testing splitData by chromatin statess")
source("helper.R", local = TRUE)

# set seed for reproducibility
set.seed(1234)

# create input
positions <- c(41224017:41224026,41224084:41224103,41234773:41234782,
               41243842:41243847,41248883:41248888)
my_data <- data.frame(Meth_Counts = sample.int(n = 20, size = length(positions), 
                                               replace = TRUE), 
                      Total_Counts = sample.int(n = 70, 
                                                size = length(positions), 
                                                replace = TRUE), 
                      Position = positions, ID = 1)
my_data <- my_data[order(my_data$Position),]

# require("annotatr")
# annotations = states_annotations(genome = "hg19", "hg19_Hmec-chromatin")
# chr17:41,196,312-41,277,381(GRCh37/hg19 by Entrez Gene)
# Size:81,070 basesOrientation:Minus strand
# r <- GRanges(seqnames = "chr17", ranges = IRanges(start = 41196312, end = 41277381))
# hits <- findOverlaps(query = annotations, subject = r)
# brca_marks <- sort(annotations[queryHits(hits)])
#
# All BRCA1 marks from hg19 genome
# > knitr::kable(as.data.frame(brca_marks))
#   |seqnames |    start|      end| width|strand |id                                      |tx_id |gene_id |symbol |type                               |
#   |:--------|--------:|--------:|-----:|:------|:---------------------------------------|:-----|:-------|:------|:----------------------------------|
#   |chr17    | 41181674| 41206474| 24801|*      |hg19_chromatin_Hmec-WeakTxn:38371       |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn        |
#   |chr17    | 41206474| 41209474|  3001|*      |hg19_chromatin_Hmec-TxnElongation:9428  |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation  |
#   |chr17    | 41209474| 41213674|  4201|*      |hg19_chromatin_Hmec-WeakTxn:38372       |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn        |
#   |chr17    | 41213674| 41214674|  1001|*      |hg19_chromatin_Hmec-TxnElongation:9429  |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation  |
#   |chr17    | 41214674| 41214874|   201|*      |hg19_chromatin_Hmec-WeakTxn:38373       |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn        |
#   |chr17    | 41214874| 41218874|  4001|*      |hg19_chromatin_Hmec-TxnElongation:9430  |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation  |
#   |chr17    | 41218874| 41220674|  1801|*      |hg19_chromatin_Hmec-WeakTxn:38374       |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn        |
#   |chr17    | 41220674| 41223874|  3201|*      |hg19_chromatin_Hmec-TxnElongation:9431  |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation  |
#   |chr17    | 41223874| 41224274|   401|*      |hg19_chromatin_Hmec-WeakTxn:38375       |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn        |
#   |chr17    | 41224274| 41228074|  3801|*      |hg19_chromatin_Hmec-TxnElongation:9432  |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation  |
#   |chr17    | 41228074| 41231274|  3201|*      |hg19_chromatin_Hmec-WeakTxn:38376       |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn        |
#   |chr17    | 41231274| 41232474|  1201|*      |hg19_chromatin_Hmec-TxnElongation:9433  |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation  |
#   |chr17    | 41232474| 41237274|  4801|*      |hg19_chromatin_Hmec-WeakTxn:38377       |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn        |
#   |chr17    | 41237274| 41248274| 11001|*      |hg19_chromatin_Hmec-TxnElongation:9434  |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation  |
#   |chr17    | 41248274| 41253874|  5601|*      |hg19_chromatin_Hmec-WeakTxn:38378       |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn        |
#   |chr17    | 41253874| 41258274|  4401|*      |hg19_chromatin_Hmec-TxnElongation:9435  |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation  |
#   |chr17    | 41258274| 41271274| 13001|*      |hg19_chromatin_Hmec-WeakTxn:38379       |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn        |
#   |chr17    | 41271274| 41274274|  3001|*      |hg19_chromatin_Hmec-TxnElongation:9436  |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation  |
#   |chr17    | 41274274| 41275674|  1401|*      |hg19_chromatin_Hmec-WeakTxn:38380       |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn        |
#   |chr17    | 41275674| 41275874|   201|*      |hg19_chromatin_Hmec-WeakEnhancer:83964  |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakEnhancer   |
#   |chr17    | 41275874| 41276074|   201|*      |hg19_chromatin_Hmec-WeakEnhancer:83965  |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakEnhancer   |
#   |chr17    | 41276074| 41276274|   201|*      |hg19_chromatin_Hmec-WeakPromoter:9953   |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakPromoter   |
#   |chr17    | 41276274| 41278474|  2201|*      |hg19_chromatin_Hmec-ActivePromoter:5485 |NA    |NA      |NA     |hg19_chromatin_Hmec-ActivePromoter |
#
# # Weak transcription BRCA1 marks from hg19 genome
#   |seqnames |    start|      end| width|strand |id                                |tx_id |gene_id |symbol |type                        |
#   |:--------|--------:|--------:|-----:|:------|:---------------------------------|:-----|:-------|:------|:---------------------------|
#   |chr17    | 41181674| 41206474| 24801|*      |hg19_chromatin_Hmec-WeakTxn:38371 |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn |
#   |chr17    | 41209474| 41213674|  4201|*      |hg19_chromatin_Hmec-WeakTxn:38372 |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn |
#   |chr17    | 41214674| 41214874|   201|*      |hg19_chromatin_Hmec-WeakTxn:38373 |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn |
#   |chr17    | 41218874| 41220674|  1801|*      |hg19_chromatin_Hmec-WeakTxn:38374 |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn |
#   |chr17    | 41223874| 41224274|   401|*      |hg19_chromatin_Hmec-WeakTxn:38375 |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn |
#   |chr17    | 41228074| 41231274|  3201|*      |hg19_chromatin_Hmec-WeakTxn:38376 |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn |
#   |chr17    | 41232474| 41237274|  4801|*      |hg19_chromatin_Hmec-WeakTxn:38377 |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn |
#   |chr17    | 41248274| 41253874|  5601|*      |hg19_chromatin_Hmec-WeakTxn:38378 |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn |
#   |chr17    | 41258274| 41271274| 13001|*      |hg19_chromatin_Hmec-WeakTxn:38379 |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn |
#   |chr17    | 41274274| 41275674|  1401|*      |hg19_chromatin_Hmec-WeakTxn:38380 |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakTxn |
#  
# # Weak Elongation BRCA1 marks from hg19 genome
#   |seqnames |    start|      end| width|strand |id                                     |tx_id |gene_id |symbol |type                              |
#   |:--------|--------:|--------:|-----:|:------|:--------------------------------------|:-----|:-------|:------|:---------------------------------|
#   |chr17    | 41206474| 41209474|  3001|*      |hg19_chromatin_Hmec-TxnElongation:9428 |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation |
#   |chr17    | 41213674| 41214674|  1001|*      |hg19_chromatin_Hmec-TxnElongation:9429 |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation |
#   |chr17    | 41214874| 41218874|  4001|*      |hg19_chromatin_Hmec-TxnElongation:9430 |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation |
#   |chr17    | 41220674| 41223874|  3201|*      |hg19_chromatin_Hmec-TxnElongation:9431 |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation |
#   |chr17    | 41224274| 41228074|  3801|*      |hg19_chromatin_Hmec-TxnElongation:9432 |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation |
#   |chr17    | 41231274| 41232474|  1201|*      |hg19_chromatin_Hmec-TxnElongation:9433 |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation |
#   |chr17    | 41237274| 41248274| 11001|*      |hg19_chromatin_Hmec-TxnElongation:9434 |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation |
#   |chr17    | 41253874| 41258274|  4401|*      |hg19_chromatin_Hmec-TxnElongation:9435 |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation |
#   |chr17    | 41271274| 41274274|  3001|*      |hg19_chromatin_Hmec-TxnElongation:9436 |NA    |NA      |NA     |hg19_chromatin_Hmec-TxnElongation |
#   
# # Weak Enhancer BRCA1 marks from hg19 genome  
#   |seqnames |    start|      end| width|strand |id                                     |tx_id |gene_id |symbol |type                             |
#   |:--------|--------:|--------:|-----:|:------|:--------------------------------------|:-----|:-------|:------|:--------------------------------|
#   |chr17    | 41275674| 41275874|   201|*      |hg19_chromatin_Hmec-WeakEnhancer:83964 |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakEnhancer |
#   |chr17    | 41275874| 41276074|   201|*      |hg19_chromatin_Hmec-WeakEnhancer:83965 |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakEnhancer |
#   
#   
# # Weak Promoter BRCA1 marks from hg19 genome
#   |seqnames |    start|      end| width|strand |id                                    |tx_id |gene_id |symbol |type                             |
#   |:--------|--------:|--------:|-----:|:------|:-------------------------------------|:-----|:-------|:------|:--------------------------------|
#   |chr17    | 41276074| 41276274|   201|*      |hg19_chromatin_Hmec-WeakPromoter:9953 |NA    |NA      |NA     |hg19_chromatin_Hmec-WeakPromoter |
#   
#   
# # Strong Promoter BRCA1 marks from hg19 genome
#   |seqnames |    start|      end| width|strand |id                                      |tx_id |gene_id |symbol |type                               |
#   |:--------|--------:|--------:|-----:|:------|:---------------------------------------|:-----|:-------|:------|:----------------------------------|
#   |chr17    | 41276274| 41278474|  2201|*      |hg19_chromatin_Hmec-ActivePromoter:5485 |NA    |NA      |NA     |hg19_chromatin_Hmec-ActivePromoter |  


# my_data
#   | Meth_Counts| Total_Counts| Position| ID| overlap
#   |-----------:|------------:|--------:|--:|
#   |          11|            4| 41224017|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           6|           43| 41224018|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |          19|           61| 41224019|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |          15|           15| 41224020|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |          14|           17| 41224021|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |          19|           43| 41224022|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |          19|           39| 41224023|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           7|           42| 41224024|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           9|           35| 41224025|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |          16|           58| 41224026|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           8|           18| 41224084|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           3|           11| 41224085|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |          16|           56| 41224086|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           2|            8| 41224087|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           5|           22| 41224088|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |          16|           42| 41224089|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           9|           63| 41224090|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           3|           42| 41224091|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           8|           20| 41224092|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           3|           33| 41224093|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           8|           61| 41224094|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |          19|           70| 41224095|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           4|           59| 41224096|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |          20|            3| 41224097|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |          19|           36| 41224098|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           7|           29| 41224099|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |          17|           65| 41224100|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           3|           30| 41224101|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           2|           57| 41224102|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |           5|           44| 41224103|  1| gap = -1: WeakTxn:38375; gap = 1000: TxnElongation:9431; WeakTxn:38375; TxnElongation:9432;
#   |          15|           61| 41234773|  1| gap = -1: WeakTxn:38377; gap = 1000; WeakTxn:38377;
#   |          15|           10| 41234774|  1| gap = -1: WeakTxn:38377; gap = 1000; WeakTxn:38377;
#   |          10|           21| 41234775|  1| gap = -1: WeakTxn:38377; gap = 1000; WeakTxn:38377;
#   |          17|           13| 41234776|  1| gap = -1: WeakTxn:38377; gap = 1000; WeakTxn:38377;
#   |           3|           29| 41234777|  1| gap = -1: WeakTxn:38377; gap = 1000; WeakTxn:38377;
#   |          16|           33| 41234778|  1| gap = -1: WeakTxn:38377; gap = 1000; WeakTxn:38377;
#   |           9|           45| 41234779|  1| gap = -1: WeakTxn:38377; gap = 1000; WeakTxn:38377;
#   |          20|           10| 41234780|  1| gap = -1: WeakTxn:38377; gap = 1000; WeakTxn:38377;
#   |          12|           31| 41234781|  1| gap = -1: WeakTxn:38377; gap = 1000; WeakTxn:38377;
#   |           5|           66| 41234782|  1| gap = -1: WeakTxn:38377; gap = 1000; WeakTxn:38377;
#   |          17|           53| 41243842|  1| gap = -1: TxnElongation:9434; gap = 1000; TxnElongation:9434;
#   |           3|           14| 41243843|  1| gap = -1: TxnElongation:9434; gap = 1000; TxnElongation:9434;
#   |          15|           38| 41243844|  1| gap = -1: TxnElongation:9434; gap = 1000; TxnElongation:9434;
#   |           5|           37| 41243845|  1| gap = -1: TxnElongation:9434; gap = 1000; TxnElongation:9434;
#   |           6|           32| 41243846|  1| gap = -1: TxnElongation:9434; gap = 1000; TxnElongation:9434;
#   |          18|           46| 41243847|  1| gap = -1: TxnElongation:9434; gap = 1000; TxnElongation:9434;
#   |          19|           60| 41248883|  1| gap = -1: WeakTxn:38378; gap = 1000; TxnElongation:9434;WeakTxn:38378;
#   |          18|           52| 41248884|  1| gap = -1: WeakTxn:38378; gap = 1000; TxnElongation:9434;WeakTxn:38378;
#   |          11|           33| 41248885|  1| gap = -1: WeakTxn:38378; gap = 1000; TxnElongation:9434;WeakTxn:38378;
#   |          16|           15| 41248886|  1| gap = -1: WeakTxn:38378; gap = 1000; TxnElongation:9434;WeakTxn:38378;
#   |           6|           34| 41248887|  1| gap = -1: WeakTxn:38378; gap = 1000; TxnElongation:9434;WeakTxn:38378;
#   |           7|           25| 41248888|  1| gap = -1: WeakTxn:38378; gap = 1000; TxnElongation:9434;WeakTxn:38378;


#--test incorrect values--# 
test_that(desc = "Test incorrect values",{
  expect_error(object = splitDataByChromatin(dat = my_data, chr = "CHR", gap = -1, 
                                             min.cpgs = 1, cell.line = "hmec", states = "all", verbose = FALSE), 
               regexp = "\\'chr\\' should be a character vector which length equals the number of rows in \\'dat\\'.")
  
  expect_error(object = splitDataByChromatin(dat = my_data, chr = NULL, gap = -1, 
                                             min.cpgs = 1, cell.line = "hmec", states = "all", verbose = FALSE), 
               regexp = "\\'chr\\' should be a character vector which length equals the number of rows in \\'dat\\'.")
  
  expect_error(object = splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                             min.cpgs = 1, cell.line = "hela", states = "all", verbose = FALSE), 
               regexp = "\\'hela\\' cell line is not supported\\. Choose one among these supported cell lines: \\'gm12878\\', \\'h1hesc\\', \\'hepg2\\', \\'hmec\\', \\'hsmm\\', \\'huvec\\', \\'k562\\', \\'nhek\\', \\'nhlf\\'.")
  
  expect_error(object = splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                             min.cpgs = 1, cell.line = "hmec", states = "ctcf", verbose = FALSE), 
               regexp = "None of the requested states are supported\\. Choose among these supported states: \\'all\\', \\'ActivePromoter\\', \\'Heterochrom\\', \\'Insulator\\', \\'PoisedPromoter\\', \\'RepetitiveCNV\\', \\'Repressed\\', \\'StrongEnhancer\\', \\'TxnElongation\\', \\'TxnTransition\\', \\'WeakEnhancer\\', \\'WeakPromoter\\', \\'WeakTxn\\'.")
  
  expect_warning(object = splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                             min.cpgs = 1, cell.line = "hmec", states = c("ctcf","ActivePromoter"), verbose = FALSE), 
               regexp = "The following unsupported states will be ignored: \\'ctcf\\'\\. Choose among these supported states: \\'all\\', \\'ActivePromoter\\', \\'Heterochrom\\', \\'Insulator\\', \\'PoisedPromoter\\', \\'RepetitiveCNV\\', \\'Repressed\\', \\'StrongEnhancer\\', \\'TxnElongation\\', \\'TxnTransition\\', \\'WeakEnhancer\\', \\'WeakPromoter\\', \\'WeakTxn\\'.")
  
  expect_warning(object = splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = "0", 
                                               min.cpgs = 1, cell.line = "hmec", states = "ActivePromoter", verbose = FALSE), 
                 regexp = "\\'gap\\' should be an integer >= -1\\. We will use the default value -1.")
  
  expect_warning(object = splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -2, 
                                               min.cpgs = 1, cell.line = "hmec", states = "ActivePromoter", verbose = FALSE), 
                 regexp = "\\'gap\\' should be an integer >= -1\\. We will use the default value -1.")
  
  expect_warning(object = splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                               min.cpgs = "-2", cell.line = "hmec", states = "ActivePromoter", verbose = FALSE), 
                 regexp = "\\'min.cpgs\\' should be a positive integer\\. We will use the default value 50.")
  
  expect_warning(object = splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                               min.cpgs = -2, cell.line = "hmec", states = "ActivePromoter", verbose = FALSE), 
                 regexp = "\\'min.cpgs\\' should be a positive integer\\. We will use the default value 50.")
  
  expect_warning(object = splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                               max.cpgs = "-2", cell.line = "hmec", states = "ActivePromoter", verbose = FALSE), 
                 regexp = "\\'max.cpgs\\' should be a positive integer\\. We will use the default value 2000.")
  
  expect_warning(object = splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                               max.cpgs = -2, cell.line = "hmec", states = "ActivePromoter"), 
                 regexp = "\\'max.cpgs\\' should be a positive integer\\. We will use the default value 2000.")
  
  expect_warning(object = splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                               min.cpgs = 50, max.cpgs = 2, cell.line = "hmec", states = "ActivePromoter", verbose = FALSE), 
                 regexp = "\\'max.cpgs\\' should be a higher than \\'min.cpgs\\'\\. We will use the default values: min.cpgs=50 and max.cpgs=2000.")
})


test_that("The length of results matches the expected length", {
  skip_if_run_bitbucket_pipeline()
  #-- min.cpgs = 1, gap = -1, cell.line = "hmec", states = "all" => 4 regions 
  # ([41224017-41224103],[41234773-41234782],[41243842-41243847],[41248883-41248888])--#
  mc_1_intersection_all <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                                cell.line = "hmec", states = "all", 
                                                gap = -1, min.cpgs = 1, verbose = FALSE)
  
  #-- min.cpgs = 1, gap = -1, cell.line = "hmec", states = "TxnElongation" => 1 region 
  # ([41243842-41243847])--#
  mc_1_intersection_TxnElongation <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                                          cell.line = "hmec", states = "TxnElongation", 
                                                          gap = -1, min.cpgs = 1, verbose = FALSE)
  
  
  #-- min.cpgs = 1, gap = 1000, cell.line = "hmec", states = "all" => 6 regions 
  # ([41224017-41224103] x 3,[41234773-41234782],[41243842-41248888], [41248883-41248888])--#
  mc_1_gap_1kb <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                       cell.line = "hmec", states = "all", 
                                       gap = 1000, min.cpgs = 1, verbose = FALSE)
  
  
  #-- min.cpgs = 1, gap = 1000, cell.line = "hmec", states = "TxnElongation" => 3 regions 
  # ([41224017-41224103] x 2,[41243842-41248888] x 1)--#
  mc_1_gap_1kb_TxnElongation <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                                     cell.line = "hmec", states = "TxnElongation", 
                                                     gap = 1000, min.cpgs = 1, verbose = FALSE)
  
  
  #-- min.cpgs = 1, gap = 1000, cell.line = "hmec", states = "WeakTxn" => 3 regions 
  # ([41224017-41224103],[41234773-41234782],[41243842-41248888])--#
  mc_1_gap_1kb_WeakTxn <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                               cell.line = "hmec", states = "WeakTxn", 
                                               gap = 1000, min.cpgs = 1, verbose = FALSE)
  
  #-- min.cpgs = 1, gap = 1000, cell.line = "hmec", states = "ActivePromoter" => 0 region
  mc_1_gap_1kb_ActivePromoter <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                                      cell.line = "hmec", states = "ActivePromoter", 
                                                      gap = 1000, min.cpgs = 1, verbose = FALSE)
  
  #-- min.cpgs = 1, gap = 1000, cell.line = "hmec", states = "TxnElongation" => 3 regions 
  # ([41224017-41224103] x 2,[41243842-41248888] x 1)--#
  #-- min.cpgs = 1, gap = 1000, cell.line = "hmec", states = "WeakTxn" => 3 regions 
  # ([41224017-41224103],[41234773-41234782],[41243842-41248888])--#
  #-- min.cpgs = 1, gap = 1000, cell.line = "hmec", states = "ActivePromoter" => 0 region
  # ===> 3 + 3 + 0 = 6
  mc_1_gap_1kb_Combined <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                                      cell.line = "hmec", states = c("WeakTxn","TxnElongation","ActivePromoter"), 
                                                      gap = 1000, min.cpgs = 1, verbose = FALSE)
  #-------#
  
  
  #-- min.cpgs = 10, gap = -1, cell.line = "hmec", states = "all" => 4 - 2 dropped regions = 2
  # ([41224017-41224103],[41234773-41234782])--#
  mc_10_intersection_all <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                                 cell.line = "hmec", states = "all", 
                                                 gap = -1, min.cpgs = 10, verbose = FALSE)
  
  #-- min.cpgs = 10, gap = -1, cell.line = "hmec", states = "TxnElongation" => 1 - 1 dropped region = 0
  mc_10_intersection_TxnElongation <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                                           cell.line = "hmec", states = "TxnElongation", 
                                                           gap = -1, min.cpgs = 10, verbose = FALSE)
  
  #-- min.cpgs = 10, gap = -1, cell.line = "hmec", states = "WeakTxn" => 3 - 1 dropped region = 2
  # ([41224017-41224103],[41234773-41234782])--#
  mc_10_intersection_WeakTxn <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                                     cell.line = "hmec", states = "WeakTxn", 
                                                     gap = -1, min.cpgs = 10, verbose = FALSE)
  
  #-- min.cpgs = 10, gap = 1000, cell.line = "hmec", states = "all" => 6 - 1 dropped regions  = 5
  # ([41224017-41224103] x 3,[41234773-41234782],[41243842-41248888])--#
  mc_10_gap_1kb <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                        cell.line = "hmec", states = "all", 
                                        gap = 1000, min.cpgs = 10, verbose = FALSE)
  
  
  #-- min.cpgs = 10, gap = 1000, cell.line = "hmec", states = "TxnElongation" => 3 regions 
  # ([41224017-41224103] x 2,[41243842-41248888] x 1)--#
  mc_10_gap_1kb_TxnElongation <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                                      cell.line = "hmec", states = "TxnElongation", 
                                                      gap = 1000, min.cpgs = 10, verbose = FALSE)
  
  
  #-- min.cpgs = 10, gap = 1000, cell.line = "hmec", states = "WeakTxn" => 3 - 1 dropped regions  = 2
  # ([41224017-41224103],[41234773-41234782])--#
  mc_10_gap_1kb_WeakTxn <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                                cell.line = "hmec", states = "WeakTxn", 
                                                gap = 1000, min.cpgs = 10, verbose = FALSE)
  
  ######
  
  #-- min.cpgs = 10, max.cpgs = 20, gap = -1, cell.line = "hmec", states = "all" => 4 - 3 dropped regions = 1
  # ([41234773-41234782])--#
  mc_10_max_20_intersection_all <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                                        cell.line = "hmec", states = "all", 
                                                        gap = -1, min.cpgs = 10, max.cpgs = 20, verbose = FALSE)
  
  
  #-- min.cpgs = 10, gap = 1000, cell.line = "hmec", states = "all" => 6 - 4 dropped regions  = 2
  # ([41234773-41234782],[41243842-41248888])--#
  mc_10_max_20_gap_1kb <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                               cell.line = "hmec", states = "all", 
                                               gap = 1000, min.cpgs = 10, max.cpgs = 20, verbose = FALSE)
  
  
  #-- min.cpgs = 5, max.cpgs = 10, gap = 1000, cell.line = "hmec", states = "TxnElongation" => 3 - 3 dropped regions  = 0
  mc_5_max_10_gap_1kb_TxnElongation <- splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                                            cell.line = "hmec", states = "TxnElongation", 
                                                            gap = 1000, min.cpgs = 5, max.cpgs = 10, verbose = FALSE)
  
  expect_true(length(mc_1_intersection_all) == 4)
  expect_true(length(mc_1_intersection_TxnElongation) == 1)
  expect_true(length(mc_1_gap_1kb) == 6)
  expect_true(length(mc_1_gap_1kb_TxnElongation) == 3)
  expect_true(length(mc_1_gap_1kb_WeakTxn) == 3)
  expect_true(length(mc_1_gap_1kb_ActivePromoter) == 0)
  expect_true(length(mc_1_gap_1kb_Combined) == 6)
  expect_true(length(mc_10_intersection_all) == 2)
  expect_true(length(mc_10_intersection_TxnElongation) == 0)
  expect_true(length(mc_10_intersection_WeakTxn) == 2)
  expect_true(length(mc_10_gap_1kb) == 5)
  expect_true(length(mc_10_gap_1kb_TxnElongation) == 3)
  expect_true(length(mc_10_gap_1kb_WeakTxn) == 2)
  expect_true(length(mc_10_max_20_intersection_all) == 1)
  expect_true(length(mc_10_max_20_gap_1kb) == 2)
  expect_true(length(mc_5_max_10_gap_1kb_TxnElongation) == 0)
})


test_that(desc = "Test region dropping",{
  skip_if_run_bitbucket_pipeline()
  # #--test region dropping--#
  m1 <- capture_messages(splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                              cell.line = "hmec", states = "TxnElongation", 
                                              gap = -1, min.cpgs = 10))
  m2 <- capture_messages(splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                              cell.line = "hmec", states = "WeakTxn", 
                                              gap = 1000, min.cpgs = 10))
  m3 <- as.character(capture_messages(splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                                           cell.line = "hmec", states = "all", 
                                                           gap = 1000, min.cpgs = 1, max.cpgs = 20)))
  # verbose = FALSE
  m4 <- capture_messages(splitDataByChromatin(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                              cell.line = "hmec", states = "TxnElongation", 
                                              gap = -1, min.cpgs = 10, verbose = FALSE))
  
  expect_match(object = m1[2], regexp = ".*17:41243842-41243847.*6.*10")
  expect_match(object = m2[2], regexp = ".*17:41248883-41248888.*6.*10")
  expect_match(object = m3[2], regexp = ".*17:41224017-41224103.*30.*20")
  expect_true(length(m4) == 0)
})
