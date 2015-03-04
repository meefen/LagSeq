# LagSeq

LagSeq is an R package for Lag-sequential Analysis.

## Installation

You can install `LagSeq` from `github` using the `devtools` package

```coffee
require(devtools)
install_github('meefen/LagSeq')
```

## Usage

Read help information of `LagSeq` and `LagSeq_Groups` first. Transform your sequential data into the required format. Basically, if you plan to compare two groups for potential transitional differences, organize your data as a data frame with columns mapping to `group`, `sequence`, and `codes`. 

```r
library(LagSeq)
load("lagseq_example_data.Rdata")
Lag_Seq_Groups(df, group=6, seq=1, codes=5)
```

## Resources

- [A blog post](http://meefen.github.io/blog/2014/04/17/temporailty-in-dialogues/)


## Disclaimer

I will not be responsible for any consequences of using this R package. Use LagSeq at your own risk.
