# class 12
Patrick Nguyen (ID: A17680785)

How many samples do we have?

``` r
exp <- read.table("rs8067378_ENSG00000172057.6.txt")
head(exp)
```

       sample geno      exp
    1 HG00367  A/G 28.96038
    2 NA20768  A/G 20.24449
    3 HG00361  A/A 31.32628
    4 HG00135  A/A 34.11169
    5 NA18870  G/G 18.25141
    6 NA11993  A/A 32.89721

``` r
nrow(exp)
```

    [1] 462

``` r
table(exp$geno)
```


    A/A A/G G/G 
    108 233 121 

``` r
library(ggplot2)
```

Boxplot

``` r
ggplot(exp) + aes(geno, exp, fill=geno) +
  geom_boxplot(notch=TRUE)
```

![](class12_files/figure-commonmark/unnamed-chunk-5-1.png)
