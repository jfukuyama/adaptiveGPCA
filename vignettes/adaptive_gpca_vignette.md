# Adaptive gPCA vignette

Here we will describe how to use the adaptiveGPCA package. We'll
demonstrate how to use it with a dataset which is included in the
package which describes the results of a microbiome experiment.  The
data included is a list with three components: otutab is a matrix
containing filtered variance-stabilized OTU proportions from each
sample, Q is an inner product matrix for this data which comes from
the phylogenetic tree describing the relationships between the OTUs,
and sampledata gives some characteristics of the sample. We first load
in the package and the data.

```r
library(adaptiveGPCA)
library(ggplot2)
data(abt_data)
```

Then we can look at adaptive gPCA on this same data set using the
adaptivegPCA command.

```r
Qeig = eigen(abt_data$Q, symmetric = TRUE)
out.agpca = adaptivegPCA(abt_data$otutab, Qeig, k = 2)
```

Alternately, if we want to use the shiny interface to choose how much
of the tree constraints to use, we can first use the gpcaFullFamily
function to create a full set of ordinations and then use the
visualizeFullFamily function to visualize the biplots at each value of
the constraint. The visualizeFullFamily function is a "shiny gadget",
and so will open a browser window where you can visualize the data set
with various constraints. Clicking "done" in this window will give as
output an object of the same format as that given by the adaptivegPCA
function, the difference being that the value of r was chosen by you
instead of automatically. 

```r
out.ff = gpcaFullFamily(abt_data$otutab, Qeig, k = 2)
out.agpca = visualizeFullFamily(out.ff,
                    sampleData = abt_data$sampledata,
                    sample_mapping = aes(x = Axis1, y = Axis2, color = condition),
                    varData = abt_data$variabledata,
                    var_mapping = aes(x = Axis1, y = Axis2, color = Phylum))
```

In either case, we can plot the results. As desired, we get a nice
biplot representation where similar species are located in similar
positions. 

```r
ggplot(data.frame(out.agpca$U, abt_data$sampledata)) +
    geom_point(aes(x = Axis1, y = Axis2, color = type, shape = ind))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

```r
ggplot(data.frame(out.agpca$QV, abt_data$variabledata)) +
    geom_point(aes(x = Axis1, y = Axis2, color = Phylum))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png)

```r
out.agpca$r
```

```
## [1] 0.4624751
```