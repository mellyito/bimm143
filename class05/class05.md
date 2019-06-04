Class 5
================
Melissa Ito
6/4/2019

Class 5 R graphics
------------------

Section 2A : line plot
======================

``` r
# can type read.table ("bimm then click "tab" and it will autofill
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
plot (weight, type ="o", 
      pch = 15, cex =1.5, lwd = 2, ylim = c(2,10),
      xlab = "Age (months)", ylab = "weight (kg)",
      main = "Baby Weight with Age")
```

![](class05_files/figure-markdown_github/unnamed-chunk-1-1.png)

``` r
# Barry's work
plot(weight$Age, weight$Weight)
```

![](class05_files/figure-markdown_github/unnamed-chunk-1-2.png)

Section 2B: barplot
===================

``` r
# sep = "\t" because a tab separates the two differnt entries
mouse <- read.table ("bimm143_05_rstats/feature_counts.txt", header = TRUE, sep ="\t")
par (mar = c (3.1, 11.1, 4.1,  2))
barplot (mouse$Count, horiz = TRUE,
         ylab ="", names.arg = mouse$Feature,
         main = "Number of features in mouse GRCm38 genome",
         las = 1, xlim = c(0, 80000))
```

![](class05_files/figure-markdown_github/unnamed-chunk-2-1.png)

Section 2C: Histogram
=====================

``` r
x <- c (rnorm(1000), rnorm(1000) +4)
hist(x, breaks = 100)
```

![](class05_files/figure-markdown_github/unnamed-chunk-3-1.png)

Section 3
=========

``` r
counts <- read.table ("bimm143_05_rstats/male_female_counts.txt", header = TRUE, sep = "\t")
# OR i can use read.delim function
counts <- read.delim("bimm143_05_rstats/male_female_counts.txt")
barplot (counts$Count, names.arg = counts$Sample, las = 2,
         col = rainbow (10))
```

![](class05_files/figure-markdown_github/unnamed-chunk-4-1.png)
