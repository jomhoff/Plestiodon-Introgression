# **_brevirostris_ Group Reticulation from UCE Data**

Ahoy.

This pipeline is the workflow used to estimate reticulation in the _Plestiodon brevirostris_ group using 3,282 Ultra-Conserved Elements (UCEs) from [Bryson Jr. et al. 2017](https://doi.org/10.1111/jbi.12989)

## **Infer Gene Trees with IQ-TREE**

The first step in this project is to take the aligned UCE data and use [IQ-TREE](https://github.com/iqtree/iqtree2) to infer unrooted gene trees for each of the 3,282 UCEs
```

```

## **Infer Species Tree with ASTRALIII**

Now that we have gene trees from IQtree2, we can use [ASTRALIII](https://github.com/smirarab/ASTRAL) to infer an unrooted species tree. ASTRALIII is a is a tool for estimating an unrooted species tree given a set of unrooted gene trees. ASTRAL is statistically consistent under the multi-species coalescent model, allowing for appropriate inference of Incomplete Lineage Sorting (ILS) and introgression. 
```

```

## **Visualize Gene Tree Discordance**


## **Estimate Introgression with QuIBL**
