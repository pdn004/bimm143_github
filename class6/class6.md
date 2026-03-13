# class06
Patrick Nguyen (A17680785)

- [Background](#background)
- [A first function](#a-first-function)
- [A new cool function](#a-new-cool-function)

## Background

Functions are at the heart of using R.Everything we do involves calling
and using functions (from data input, analysis to results).

All functions in R have at least 3 things:

- A **name** the thing we use to call the function.
- One or more input **arguments** that are comma seperated
- The \***body**, liens of code between curly brackets {} that does the
  work of the function.

## A first function

Let’s write a silly wee function to add some numbers:

``` r
add <- function(x) {
  x + 1
}
```

Let’s try it out

``` r
add(100)
```

    [1] 101

Will this work?

``` r
add( c(100, 200, 300) )
```

    [1] 101 201 301

Modify to be more useful and add more than just 1

``` r
add <- function(x, y=1) {
  x + y
}
```

``` r
add(100, 10)
```

    [1] 110

Will this work?

``` r
add(100)
```

    [1] 101

``` r
plot(1:10, col="blue", typ="b")
```

![](class6_files/figure-commonmark/unnamed-chunk-7-1.png)

> **N.B.** Input arguments can be either **required** or **optional**.
> The later have a fall-back default that is specifed in the function
> code with an equals sign.

``` r
#add(x=100, y=200, z=300)
```

\##A second function

All functions in R look like this

    name <- function(arg) {
      body
    }

The ‘sample()’ function in R…

``` r
sample(1:10, size=4)
```

    [1] 2 4 5 3

> Q. Return 12 numbers picked randomly from the input 1:10

``` r
sample(1:10, size=12, replace=TRUE)
```

     [1]  4  6  2  4  5  3  4  1  7  8  9 10

> Q. Write the code that generates a random 12 nucleotide long DNA
> sequence?

``` r
sample( c("A","T","G","C"), size=12, replace=TRUE )
```

     [1] "C" "C" "G" "G" "A" "A" "A" "T" "T" "T" "T" "T"

> Q. Write a first version function called ‘generate_dna()’ that
> generates a user specified length ‘n’ random DNA sequence?

    name <- function(arg) {
      body
    }

``` r
generate_dna <- function(n=3) {
  sample( c("A","T","G","C"), size=n, replace=TRUE )
}
```

``` r
generate_dna(100)
```

      [1] "G" "G" "T" "T" "C" "T" "A" "G" "T" "C" "A" "A" "C" "A" "A" "G" "T" "A"
     [19] "C" "G" "T" "A" "T" "C" "C" "T" "T" "C" "C" "G" "G" "C" "T" "C" "T" "C"
     [37] "G" "T" "G" "G" "C" "A" "C" "C" "C" "C" "T" "T" "C" "G" "C" "C" "A" "G"
     [55] "G" "T" "C" "C" "T" "G" "C" "G" "G" "G" "A" "C" "A" "G" "T" "A" "G" "C"
     [73] "C" "C" "C" "G" "A" "T" "G" "A" "A" "T" "G" "C" "T" "T" "G" "G" "C" "C"
     [91] "C" "C" "C" "C" "A" "G" "G" "G" "A" "A"

> Q. Modify your function to return a FASTA like sequence so rather than
> \[1\] “C” “G” “C” “A” “A” “A” “C” “T” “A” “C” “C” “T” we want “GCAAT”

``` r
generate_dna <- function(n=3) {
  ans <- sample( c("A","T","G","C"), size=n, replace=TRUE )
  ans <- paste(ans, collapse ="")
  return(ans)
}
```

``` r
generate_dna(10)
```

    [1] "AACTGGGCGA"

An example

``` r
# Example pattern (not using your bases)
x <- c("H","E","L","L","O")

paste(x, collapse = "****")
```

    [1] "H****E****L****L****O"

``` r
# returns "HELLO"
```

> Q. Give the user an option to return FASTA format output sequence or
> standard multi-element vector format?

``` r
generate_dna <- function(n=3, fasta=TRUE) {
  ans <- sample( c("A","T","G","C"), size=n, replace=TRUE )
  
if(fasta) {
  ans <- paste(ans, collapse ="")
  cat("Hello...")
} else {
  cat("...is it me you are looking for...") 
  }
  
  return(ans)
  }
```

``` r
generate_dna(10)
```

    Hello...

    [1] "AGATCTGCGA"

``` r
generate_dna(10, fasta=FALSE)
```

    ...is it me you are looking for...

     [1] "A" "C" "C" "C" "T" "A" "C" "C" "T" "T"

## A new cool function

> Q.Write a function called ‘generate_protein()’ that generates a user
> specified length protein sequence in FASTA like format?

``` r
generate_protein <- function(n=3, fasta=TRUE) {
  aa <- c( "A","R","N","D",
           "C","Q","E","G",
           "H","I","L","K",
           "M","F","P","S",
           "T","W","Y","V")
  gen <- sample(aa,size=n, replace=TRUE)
  
  if(fasta) {
    prn <- paste(gen, collapse = "")
  }
  
  return(prn)
  
}
```

``` r
generate_protein(10)
```

    [1] "RGNENMSHLL"

> Q. Use your new ‘generate_protein()’ function to generate sequences
> between lengths 6 and 12 amino acids in length and check of any of
> these are unique in nature (i.e. found in the MR database at NCBI)?

``` r
generate_protein(6)
```

    [1] "IKYTPN"

``` r
generate_protein(7)
```

    [1] "FTITPDY"

``` r
generate_protein(8)
```

    [1] "TRKFADTW"

``` r
generate_protein(9)
```

    [1] "VQKEAKTHM"

``` r
generate_protein(10)
```

    [1] "HGCDASEKYV"

``` r
generate_protein(11)
```

    [1] "LTNNNTSPMGN"

``` r
generate_protein(12)
```

    [1] "FQPCDFPSDNVY"

Or we could do a ‘for()’ loop:

``` r
for(i in 6:12) {
  cat(">", i, sep="", "\n")
  cat( generate_protein(i), "\n" )
}
```

    >6
    NQLKQR 
    >7
    LAQAFCM 
    >8
    FCKIQTRL 
    >9
    MAVSADGPA 
    >10
    CWFYHTATFS 
    >11
    WLCGIPMKYVH 
    >12
    KIFHPIVTWQFP 
