# CloneStrat
A package for scaled and allelic imbalance adjusted clonal deconvolution

## Installation instructions


```{r, eval=FALSE, echo=TRUE}
install.packages(c("mclust","fpc","sequenza","vcfR","bootcluster","devtools",
                   "factoextra","FactoMineR","RcppArmadillo","installr"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("copynumber")
```


If no error message has popped up so far we are now ready to install the package. Best practice would be to change your working directory to the directory where all the downloaded files are saved. Then just remove *~PATH* part from the code below. OR, another way would be to replace the *~PATH* in the following line with the directory address where the package zip file is saved.

*Note* Here on out it is assumed that the working directory has been changed to where the downloaded files are. Hence no additional path is defined.

```{r, eval=FALSE, echo=TRUE}
install.packages("~PATH/CloneStrat_0.1.1.tar.gz", repos = NULL, type = "source")
```

We can now load the package in global environment.

```{r}
require(CloneStrat)
```

## Preliminary usage

Let's take a look at the imaginary whole exome sequencing data built into the package.

```{r, eval=FALSE, echo=TRUE}
?test.dat
head(test.dat)

```

We will use this data to perform our first clonal deconvolution. A set of example user input is given below but there are lot of other choices in methods and inputs that can be provided.

```{r, eval=FALSE, echo=TRUE}
res.1 <- cluster.doc(test.dat, sample = 1, vaf = 2, 
                   optimization.method = 'GMM', clustering.method = 'hkm')

## example user input:
## What is the suspected chromosomal segmentation profile of the sample: 2 + 2
## How many clonal VAF clouds do you think are present: 2
## How many sub-clonal VAF clouds do you think are present: 2
## Would you like to see my suggestion instead? (yes/no): no
```

**Figure 1** shows how the distribution of variant allele frequencies look for the 8 simulated samples. Depending on the structure of the spread it is concievable that the allelic segmentation is a balanced 2 + 2

**Figure 2** shows the changes in *Bayesian Information criteria* (BIC) estimated based on expectation and variance of the clustering fit.

**Figure 3** shows the clustered samples here depict the distribution of clonal and sub-clonal variants.

<center>

![VAF distribution](https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.1.png)

![BIC changes](https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.2.png)

![Clonal deconvolution](https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.3.png)

</center>

## Scaling

This data is rather clean and the clusters are more or less obviously identifiable with a visual inspection. But real data tend to be much more noisy and the clonal clusters  are unassumingmore often than not. Let's use the data set from the neuroblastoma patient sample *ES* to see how messy real data can be.

```{r, eval=FALSE, echo=TRUE}
es<-read.table("ES_all_2+2.txt",header=T,stringsAsFactors=F)

## Let's check the VAF distribution of the samples
ggplot(es, aes(x=sample, y=mut, col=as.factor(sample))) + geom_point()
```

<center>

![VAF distribution for ES](https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.4.png)

</center>

What happend if we deconvolute this data assuming there are two clones and one subclone?

```{r, eval=FALSE, echo=TRUE}
res.2 <- cluster.doc(es, sample = 1, vaf = 2, 
                   optimization.method = 'GMM', clustering.method = 'hkm')

## example user input:
## What is the suspected chromosomal segmentation profile of the sample: 2 + 2
## How many clonal VAF clouds do you think are present: 2
## How many sub-clonal VAF clouds do you think are present: 1
```

<center>

![Clonal deconvolution for ES](https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.5.png)

</center>

As it is evident that the VAF distribution here is quite varied among samples, let's now use *Probabilistic Quotient Normalization* with the *Cancer Cell Fractions* (CCF)

```{r, eval=FALSE, echo=TRUE}
es.sc<-seqn.scale(es,2,3)

## Let's check the scaled VAF distribution of the samples
ggplot(es.sc, aes(x=sample, y=scaled.vaf, col=as.factor(sample))) + geom_point()
```

<center>

![VAF distribution for scaled ES](https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.6.png)

</center>

This step has pretty much rescaled the VAFs for the samples that had a relatively lower CCFs. Let's try the deconvolution now.

```{r, eval=FALSE, echo=TRUE}
res.2 <- cluster.doc(es.sc, sample = 1, vaf = 3, 
                   optimization.method = 'GMM', clustering.method = 'hkm')

## example user input:
## What is the suspected chromosomal segmentation profile of the sample: 2 + 2
## How many clonal VAF clouds do you think are present: 2
## How many sub-clonal VAF clouds do you think are present: 2
## Would you like to see my suggestion instead? (yes/no): no
```

<center>

![Clonal deconvolution for scaled ES](https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.7.png)

</center>

Can you plot the difference in the subclonal distributions for some of the samples now?

Now lets us assume we are not very satisfied as how things stand for sample 6 and sample 7. Instead of 2 clonal and 1 sub-clonal clouds we would rather see 

```{r, eval=FALSE, echo=TRUE}
res.3 <- cluster.doubt(res.2,1,3,c("sample_6","sample_7"),c(2,2,2,2))
```

<center>

![User rectified clonal deconvolution for scaled ES](https://github.com/Subhayan18/CloneStrat/blob/master/source/test.dat.8.png)

</center>

## Estimation of allelic segmentation

When the allelic make up is unavilable to the user, it can be estimated given the sequence reads from the constitutional DNA is also present. This can generally be obtained from a .vcf file before the variant calling.

```{r, eval=FALSE, echo=TRUE}
m16 <- vcfR::read.vcfR("case_m16_P1883_147_mutect.vcf")
CN.est.1 <- CopySeg(m16,tumor.sample="P1883_147",normal.sample="P1883_118",
                  DP="DP",AD="AD",file.name="sample")
## Here we have estimated allelic segmentation for the tumor sample 'P1883_147'
## Only first 30 rows are visible below
```

## Auxiliary functions

**Mutect2** is a popular variant caller. The package offers an extension that can convert a Mutect2 output to a data file compatible with the functions described.

```{r, eval=TRUE, echo=TRUE}
WES <- readxl::read_excel("WES.xlsx")  ##This step requires the package 'readxl'. If unavailable please install.
sample.name <- c("3407_18","6962_18","6963_18","6964_18","6965_18")
CS.dat <- mutect2.qc(WES,sample.name)
```

**A test** is provided to check the goodness of the cluster fit. This will indicate existence of outliers.
```{r, eval=FALSE, echo=TRUE}
CS.test<-T.goodness.test(es)$rej
## The variants that can only belong to one clone is shown
```
##Coming soon

Functions to plot copy number estimation and figure out allelic composition. I am also testing a new method to estimating copy numbers which gives user more freedom to tweak.

More features will be added gradually. If you have feature that you'd like to see incorporated in `CloneStrat`, please send 
a request!

---
title: "Your title"
output: 
  html_document:
    includes:
      after_body: footer.html
---

&nbsp;
<hr />
<p style="text-align: center;">A work by <a href="https://github.com/holtzy/">Yan Holtz</a></p>
<p style="text-align: center;"><span style="color: #808080;"><em>Yan.holtz.data@gmail.com</em></span></p>

<!-- Add icon library -->
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css">

<!-- Add font awesome icons -->
<p style="text-align: center;">
    <a href="https://twitter.com/r_graph_gallery?lang=en" class="fa fa-twitter"></a>
    <a href="https://www.linkedin.com/in/yan-holtz-2477534a/" class="fa fa-linkedin"></a>
    <a href="https://github.com/holtzy/" class="fa fa-github"></a>
</p>

&nbsp;
