#My dir.
setwd("~/Dropbox/laboratorio/scripts/fujita/")

# Up main function
source("~/Dropbox/laboratorio/scripts/european_pop/load_functions.R")

#My packages
library(data.table)
library(splitstackshape)
library(ggplot2)
library(GGally)
library(knitr)
library(dplyr)
library(dtplyr)
library(ggthemes)
library(tidyverse)


pvalue <- fread("~/Dropbox/laboratorio/scripts/paper_1/data/pvalue.csv")

head(pvalue)
m <- ggplot(pvalue, aes(x=value))
m + geom_histogram(bins = 35, colour = "grey30", fill = "lightsteelblue4") + facet_grid(~ measure ~  boot) + opa()

