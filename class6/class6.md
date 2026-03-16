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

    [1] 1 2 8 4

> Q. Return 12 numbers picked randomly from the input 1:10

``` r
sample(1:10, size=12, replace=TRUE)
```

     [1] 8 6 5 5 4 8 7 9 7 4 8 8

> Q. Write the code that generates a random 12 nucleotide long DNA
> sequence?

``` r
sample( c("A","T","G","C"), size=12, replace=TRUE )
```

     [1] "G" "C" "A" "C" "T" "C" "T" "G" "G" "A" "A" "T"

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

      [1] "C" "T" "A" "G" "C" "T" "G" "T" "C" "A" "C" "A" "G" "A" "A" "A" "A" "C"
     [19] "A" "G" "A" "G" "T" "G" "G" "T" "T" "G" "G" "T" "A" "T" "A" "C" "C" "T"
     [37] "T" "G" "C" "G" "T" "A" "C" "C" "C" "G" "T" "C" "T" "C" "C" "G" "A" "A"
     [55] "C" "T" "T" "C" "G" "A" "T" "G" "A" "A" "G" "C" "C" "C" "T" "T" "A" "G"
     [73] "G" "T" "G" "C" "C" "T" "T" "T" "C" "A" "T" "C" "A" "G" "C" "T" "T" "C"
     [91] "G" "G" "C" "C" "C" "G" "A" "G" "G" "T"

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

    [1] "TTAACGCGCG"

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

    [1] "TCGGCCGTGA"

``` r
generate_dna(10, fasta=FALSE)
```

    ...is it me you are looking for...

     [1] "G" "C" "T" "A" "A" "A" "C" "T" "C" "G"

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

    [1] "QHVRGPWPCR"

> Q. Use your new ‘generate_protein()’ function to generate sequences
> between lengths 6 and 12 amino acids in length and check of any of
> these are unique in nature (i.e. found in the MR database at NCBI)?

``` r
generate_protein(6)
```

    [1] "SKNPCE"

``` r
generate_protein(7)
```

    [1] "QAGARMD"

``` r
generate_protein(8)
```

    [1] "VDFYILQQ"

``` r
generate_protein(9)
```

    [1] "HLMNLCIEI"

``` r
generate_protein(10)
```

    [1] "WKVRLKTHDM"

``` r
generate_protein(11)
```

    [1] "DEDMIVNTECD"

``` r
generate_protein(12)
```

    [1] "HKLSNDLMAPTL"

Or we could do a ‘for()’ loop:

``` r
for(i in 6:12) {
  cat(">", i, sep="", "\n")
  cat( generate_protein(i), "\n" )
}
```

    >6
    VSFPWK 
    >7
    KFDCTDH 
    >8
    ACKNCCSA 
    >9
    CYKNEEEGL 
    >10
    KDQLQHWDHP 
    >11
    HVLGEDKLAAC 
    >12
    IKGGDVCWCKSV 
