context("testing splitData by genomic regions")
source("helper.R", local = TRUE)

# set seed for reproducibility
set.seed(1234)

# create input
positions <- c(43072000:43072009,43072067:43072086,43082756:43082765,
               43091825:43091830,43096866:43096871)
my_data <- data.frame(Meth_Counts = sample.int(n = 20, size = length(positions), 
                                               replace = TRUE), 
                      Total_Counts = sample.int(n = 70, 
                                                size = length(positions), 
                                                replace = TRUE), 
                      Position = positions, ID = 1)
my_data <- my_data[order(my_data$Position),]

# require("annotatr")
# annotations = build_annotations(genome = "hg38", 
# annotations = "hg38_genes_promoters")
# prom <- sort(annotations[!is.na(annotations$symbol) & 
# annotations$symbol == "BRCA1"])[1:5]
#
# BRCA1 promoters from hg38 genome
# GAP = -1 (intersection)
# seqnames         ranges strand |              id             tx_id     gene_id      symbol                 type
# <Rle>         <IRanges>  <Rle> |     <character>       <character> <character> <character>          <character>
# chr17 43071067-43072066      - | promoter:184718 ENST00000472490.1         672       BRCA1 hg38_genes_promoters
# chr17 43082761-43083760      - | promoter:184719 ENST00000621897.1         672       BRCA1 hg38_genes_promoters
# chr17 43091825-43092824      - | promoter:184721 ENST00000461574.1         672       BRCA1 hg38_genes_promoters
# chr17 43091852-43092851      - | promoter:184708 ENST00000644379.1         672       BRCA1 hg38_genes_promoters
# chr17 43095867-43096866      - | promoter:184724 ENST00000412061.3         672       BRCA1 hg38_genes_promoters
#
#
# GAP = 0 (adjacent)
# seqnames            ranges strand |              id             tx_id     gene_id      symbol                 type
# <Rle>         <IRanges>  <Rle> |     <character>       <character> <character> <character>          <character>
# chr17 43071066-43072067      - | promoter:184718 ENST00000472490.1         672       BRCA1 hg38_genes_promoters
# chr17 43082760-43083761      - | promoter:184719 ENST00000621897.1         672       BRCA1 hg38_genes_promoters
# chr17 43091824-43092825      - | promoter:184721 ENST00000461574.1         672       BRCA1 hg38_genes_promoters
# chr17 43091851-43092852      - | promoter:184708 ENST00000644379.1         672       BRCA1 hg38_genes_promoters
# chr17 43095866-43096867      - | promoter:184724 ENST00000412061.3         672       BRCA1 hg38_genes_promoters
#
#
# GAP = 5
# seqnames         ranges strand |              id             tx_id     gene_id      symbol                 type
# <Rle>         <IRanges>  <Rle> |     <character>       <character> <character> <character>          <character>
# chr17 43071062-43072071      - | promoter:184718 ENST00000472490.1         672       BRCA1 hg38_genes_promoters
# chr17 43082756-43083765      - | promoter:184719 ENST00000621897.1         672       BRCA1 hg38_genes_promoters
# chr17 43091820-43092829      - | promoter:184721 ENST00000461574.1         672       BRCA1 hg38_genes_promoters
# chr17 43091847-43092856      - | promoter:184708 ENST00000644379.1         672       BRCA1 hg38_genes_promoters
# chr17 43095862-43096871      - | promoter:184724 ENST00000412061.3         672       BRCA1 hg38_genes_promoters

# my_data
# | Meth_Counts| Total_Counts|Chrom | Position| ID| overlap
# |-----------:|------------:|:-----|--------:|--:|
# |          15|            8|17    | 43072000|  1|gap = -1 [promoter:184718]
# |           4|           21|17    | 43072001|  1|gap = -1 [promoter:184718]
# |          13|           34|17    | 43072002|  1|gap = -1 [promoter:184718]
# |          12|           21|17    | 43072003|  1|gap = -1 [promoter:184718]
# |          10|           64|17    | 43072004|  1|gap = -1 [promoter:184718]
# |          12|           39|17    | 43072005|  1|gap = -1 [promoter:184718]
# |          10|           39|17    | 43072006|  1|gap = -1 [promoter:184718]
# |          17|           19|17    | 43072007|  1|gap = -1 [promoter:184718]
# |          17|           45|17    | 43072008|  1|gap = -1 [promoter:184718]
# |           6|           30|17    | 43072009|  1|gap = -1 [promoter:184718]
# |          11|           66|17    | 43072067|  1|gap =  0 [promoter:184718]
# |           6|           61|17    | 43072068|  1|gap =  5 [promoter:184718]
# |           8|           22|17    | 43072069|  1|gap =  5 [promoter:184718]
# |          19|           46|17    | 43072070|  1|gap =  5 [promoter:184718]
# |          15|           65|17    | 43072071|  1|gap =  5 [promoter:184718]
# |           3|           25|17    | 43072072|  1|gap =  5 [promoter:184718]
# |           2|           27|17    | 43072073|  1|
# |          13|           11|17    | 43072074|  1|
# |          10|           69|17    | 43072075|  1|
# |          13|            2|17    | 43072076|  1|
# |           7|           27|17    | 43072077|  1|
# |          14|           44|17    | 43072078|  1|
# |          17|           17|17    | 43072079|  1|
# |          17|            1|17    | 43072080|  1|
# |          20|           70|17    | 43072081|  1|
# |           8|           28|17    | 43072082|  1|
# |           1|           20|17    | 43072083|  1|
# |           9|           68|17    | 43072084|  1|
# |          17|           65|17    | 43072085|  1|
# |          14|           14|17    | 43072086|  1|
#
# |          10|            3|17    | 43082756|  1|gap =  5 [promoter:184719]
# |          19|           68|17    | 43082757|  1|gap =  5 [promoter:184719]
# |          16|           37|17    | 43082758|  1|gap =  5 [promoter:184719]
# |           2|           17|17    | 43082759|  1|gap =  5 [promoter:184719]
# |          11|           40|17    | 43082760|  1|gap =  0 [promoter:184719]
# |          20|           25|17    | 43082761|  1|gap = -1 [promoter:184719]
# |          11|           43|17    | 43082762|  1|gap = -1 [promoter:184719]
# |          16|           65|17    | 43082763|  1|gap = -1 [promoter:184719]
# |          15|           31|17    | 43082764|  1|gap = -1 [promoter:184719]
# |           6|           18|17    | 43082765|  1|gap = -1 [promoter:184719]
#
# |          10|           24|17    | 43091825|  1|gap = -1 [promoter:184721] [3UTR:125803]
# |           7|           66|17    | 43091826|  1|gap = -1 [promoter:184721] [3UTR:125803]
# |           6|           29|17    | 43091827|  1|gap = -1 [promoter:184721] [3UTR:125803]
# |          17|           38|17    | 43091828|  1|gap = -1 [promoter:184721] [3UTR:125803]
# |           4|            7|17    | 43091829|  1|gap = -1 [promoter:184721] [3UTR:125803]
# |          15|            2|17    | 43091830|  1|gap = -1 [promoter:184721] [3UTR:125803]
#
# |          19|           60|17    | 43096866|  1|gap = -1 [promoter:184724]
# |          18|           52|17    | 43096867|  1|gap =  0 [promoter:184724]
# |          11|           33|17    | 43096868|  1|gap =  5 [promoter:184724]
# |          16|           15|17    | 43096869|  1|gap =  5 [promoter:184724]
# |           6|           34|17    | 43096870|  1|gap =  5 [promoter:184724]
# |           7|           25|17    | 43096871|  1|gap =  5 [promoter:184724]


#--test incorrect values--# 
test_that(desc = "Test incorrect values",{
  expect_error(object = splitDataByGene(dat = my_data, chr = "CHR", gap = -1, 
                                        types = "exon", min.cpgs = 1, organism = "human", build = "hg19", verbose = FALSE), 
               regexp = "\\'chr\\' should be a character vector which length equals the number of rows in \\'dat\\'.")
  expect_error(object = splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                        types = "exon", min.cpgs = 1, organism = "lupus", build = "hg19", verbose = FALSE), 
               regexp = "\\'lupus\\' organism is not supported\\. Choose one among these supported types: \\'human\\', \\'mouse\\', \\'rat\\', \\'fly\\'.")
  expect_error(object = splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                        types = "exon", min.cpgs = 1, organism = "human", build = "hg18", verbose = FALSE), 
               regexp = "\\'hg18\\' build is not supported\\. Choose one among these supported builds: \\'hg38\\', \\'hg19\\'.")
  expect_error(object = splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                        types = "exon", min.cpgs = 1, organism = "mouse", build = "mm10", verbose = FALSE), 
               regexp = "Please install the following package\\(s\\) in order to perform the requested partiniong: \\'TxDb.Mmusculus.UCSC.mm10.knownGene\\', \\'org.Mm.eg.db\\'.")
  expect_warning(object = splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = "0", 
                                          types = "exon", min.cpgs = 1, organism = "human", build = "hg38", verbose = FALSE), 
                 regexp = "\\'gap\\' should be an integer >= -1\\. We will use the default value -1.")
  expect_warning(object = splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -2, 
                                          types = "exon", min.cpgs = 1, organism = "human", build = "hg38", verbose = FALSE), 
                 regexp = "\\'gap\\' should be an integer >= -1\\. We will use the default value -1.")
  expect_warning(object = splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                          types = "exon", min.cpgs = "-2", organism = "human", build = "hg38", verbose = FALSE), 
                 regexp = "\\'min.cpgs\\' should be a positive integer\\. We will use the default value 50.")
  expect_warning(object = splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                          types = "exon", min.cpgs = -2, organism = "human", build = "hg38", verbose = FALSE), 
                 regexp = "\\'min.cpgs\\' should be a positive integer\\. We will use the default value 50.")
  expect_warning(object = splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                          types = "exon", max.cpgs = "-2", organism = "human", build = "hg38", verbose = FALSE), 
                 regexp = "\\'max.cpgs\\' should be a positive integer\\. We will use the default value 2000.")
  expect_warning(object = splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                          types = "exon", max.cpgs = -2, organism = "human", build = "hg38", verbose = FALSE), 
                 regexp = "\\'max.cpgs\\' should be a positive integer\\. We will use the default value 2000.")
  expect_warning(object = splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                          types = "exon", min.cpgs = 50, max.cpgs = 2, organism = "human", build = "hg38", verbose = FALSE), 
                 regexp = "\\'max.cpgs\\' should be a higher than \\'min.cpgs\\'\\. We will use the default values: min.cpgs=50 and max.cpgs=2000.")
  expect_warning(object = splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                        types = c("promoter", "coding"), min.cpgs = 2, organism = "human", build = "hg19", verbose = FALSE), 
               regexp = "The following unsupported types will be ignored: \\'coding\\'\\. Choose among these supported types: \\'all\\', \\'upstream\\', \\'threeprime\\', \\'fiveprime\\', \\'promoter\\', \\'exon\\', \\'intron\\'.")
  expect_error(object = splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), gap = -1, 
                                        types = "coding", min.cpgs = 2, organism = "human", build = "hg19", verbose = FALSE), 
               regexp = "None of the requested types are supported\\. Choose among these supported types: \\'all\\', \\'upstream\\', \\'threeprime\\', \\'fiveprime\\', \\'promoter\\', \\'exon\\', \\'intron\\'.")
})

test_that("The length of results matches the expected length", {
  skip_if_run_bitbucket_pipeline()
  
  #-- min.cpgs = 1, gap = -1, organism = "human", build = "hg38" => 4 regions 
  # ([43072000-43072009],[43082761-43082765],[43091825-43091830],[43096866])--#
  mc_1_intersection <- splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                       organism = "human", build = "hg38", 
                                       types = "promoter", gap = -1, min.cpgs = 1, verbose = FALSE)
  
  #-- min.cpgs = 1, gap = -1, organism = "human", build = "hg38" => 4 regions 
  # promoter: ([43072000-43072009],[43082761-43082765],[43091825-43091830],[43096866])--#
  # threeprime: [43091825-43091830]
  mc_1_combined <- splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                       organism = "human", build = "hg38", 
                                       types = c("promoter","threeprime"), gap = -1, min.cpgs = 1, verbose = FALSE)
  
  
  #-- min.cpgs = 1, gap = 0, organism = "human", build = "hg38" => 4 regions 
  # ([43072000-43072067],[43082760-43082765],[43091825-43091830], [43096866-43096867])--#
  mc_1_gap_adjacent <- splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                       organism = "human", build = "hg38", 
                                       types = "promoter", gap = 0, min.cpgs = 1, verbose = FALSE)
  
  #-- min.cpgs = 1, gap = 5, organism = "human", build = "hg38" => 4 regions 
  # ([43072000-43072072],[43082756-43082765],[43091825-43091830], [43096866-43096871])--#
  mc_1_gap_5 <- splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                organism = "human", build = "hg38", 
                                types = "promoter", gap = 5, min.cpgs = 1, verbose = FALSE)
  
  #-- min.cpgs = 5, gap = -1, organism = "human", build = "hg38" => 3 regions (4 - 1 drop)
  # ([43072000-43072009],[43082761-43082765],[43091825-43091830])--#
  mc_5_intersection <- splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                       organism = "human", build = "hg38", 
                                       types = "promoter", gap = -1, min.cpgs = 5, verbose = FALSE)
  
  #-- min.cpgs = 5, gap = 0, organism = "human", build = "hg38" => 3 regions (4 - 1 drop)
  # ([43072000-43072067],[43082760-43082765],[43091825-43091830])--#
  mc_5_gap_adjacent <- splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                       organism = "human", build = "hg38", 
                                       types = "promoter", gap = 0, min.cpgs = 5, verbose = FALSE)
  
  #-- min.cpgs = 1, gap = 5, organism = "human", build = "hg38" => 4 regions 
  # ([43072000-43072072],[43082756-43082765],[43091825-43091830], [43096866-43096871])--#
  mc_5_gap_5 <- splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                organism = "human", build = "hg38", 
                                types = "promoter", gap = 5, min.cpgs = 5, verbose = FALSE)
  
  #-- min.cpgs = 10, gap = -1, organism = "human", build = "hg38" => 1 region (4 - 3 drop)
  # ([43072000-43072009])--#
  mc_10_intersection <- splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                        organism = "human", build = "hg38", 
                                        types = "promoter", gap = -1, min.cpgs = 10, verbose = FALSE)
  
  #-- min.cpgs = 5, gap = 0, organism = "human", build = "hg38" => 1 region (4 - 3 drop)
  # ([43072000-43072067])--#
  mc_10_gap_adjacent <- splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                        organism = "human", build = "hg38", 
                                        types = "promoter", gap = 0, min.cpgs = 10, verbose = FALSE)
  
  #-- min.cpgs = 1, gap = 5, organism = "human", build = "hg38" => 2 regions (4 - 2 drop) 
  # ([43072000-43072072],[43082756-43082765])--#
  mc_10_gap_5 <- splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                 organism = "human", build = "hg38", 
                                 types = "promoter", gap = 5, min.cpgs = 10, verbose = FALSE)
  
  #-- min.cpgs = 5, max.cpgs = 10, gap = 5, organism = "human", build = "hg38" => 3 regions + 1 drop [43072000-43072072]
  # ([43082756-43082765],[43091825-43091830], [43096866-43096871])--#
  mc_5_max_10_gap_5 <- splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                       organism = "human", build = "hg38", 
                                       types = "promoter", gap = 5, 
                                       min.cpgs = 5, max.cpgs = 10, verbose = FALSE)
  
  #-- min.cpgs = 1, max.cpgs = 10, gap = -1, organism = "human", build = "hg38" => 4 regions
  mc_1_max_10_intersection <- splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                              organism = "human", build = "hg38", 
                                              types = "promoter", gap = -1, min.cpgs = 1,
                                              max.cpgs = 10, verbose = FALSE)
  
  #-- min.cpgs = 1, gap = 5, organism = "human", build = "hg38" => 1 regions (4 - 3 drop) 
  # ([43082756-43082765])--#
  mc_10_max_15_gap_5 <- splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                        organism = "human", build = "hg38", 
                                        types = "promoter", gap = 5, min.cpgs = 10, 
                                        max.cpgs = 15, verbose = FALSE)
  
  expect_true(length(mc_1_intersection) == 4)
  expect_true(length(mc_1_combined) == 5)
  expect_true(length(mc_1_gap_adjacent) == 4)
  expect_true(length(mc_1_gap_5) == 4)
  expect_true(length(mc_5_intersection) == 3)
  expect_true(length(mc_5_gap_adjacent) == 3)
  expect_true(length(mc_5_gap_5) == 4)
  expect_true(length(mc_10_intersection) == 1)
  expect_true(length(mc_10_gap_adjacent) == 1)
  expect_true(length(mc_10_gap_5) == 2)
  expect_true(length(mc_5_max_10_gap_5) == 3)
  expect_true(length(mc_1_max_10_intersection) == 4)
  expect_true(length(mc_10_max_15_gap_5) == 1)
})

test_that(desc = "Test region dropping",{
  skip_if_run_bitbucket_pipeline()
  
  # #--test region dropping--# 
  m1 <- capture_messages(splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                         organism = "human", build = "hg38", 
                                         types = "promoter", gap = -1, min.cpgs = 5))
  m2 <- capture_messages(splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                         organism = "human", build = "hg38", 
                                         types = "promoter", gap = 0, min.cpgs = 5))
  m3 <- capture_messages(splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                         organism = "human", build = "hg38", 
                                         types = "promoter", gap = 5, 
                                         min.cpgs = 5, max.cpgs = 10))
  m4 <- capture_messages(splitDataByGene(dat = my_data, chr = rep("chr17", nrow(my_data)), 
                                         organism = "human", build = "hg38", 
                                         types = "promoter", gap = 5, 
                                         min.cpgs = 5, max.cpgs = 10,
                                         verbose = FALSE))
  
  expect_match(object = m1[2], regexp = ".*17:43096866.*1.*5")
  expect_match(object = m2[2], regexp = ".*17:43096866-43096867.*2.*5")
  expect_match(object = m3[2], regexp = ".*17:43072000-43072072.*16.*10")
  expect_true(length(m4) == 0)
})
