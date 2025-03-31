# **Analyzing ChIP-seq and CUT&RUN datasets**
*YNA lab bioinformatics training series, pt 2*



## What are ChIP-seq and CUT&RUN?

ChIP-seq             |  CUT&RUN
:-------------------------:|:-------------------------:
![ChIP-seq basic steps](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/ChIP-blog-figure-1.jpg)  |  ![CUT&RUN basic steps](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/cut-run-blog-figure-1.png)

*Images from epicypher*

<br />
$$\textnormal{\color{darkorange}ChIP-seq}$$ and $$\textnormal{\color{darkorange}CUT&RUN}$$ $$\textnormal{\color{darkorange}kernel}$$ are methods for determining the positions of histones/histone modifications or DNA binding proteins on DNA. Both methods rely on antibodies as an important initial step to select for proteins of interest. However, they differ in how their antibody-based DNA selection is performed. Key differences between ChIP-seq and CUT&RUN are: $$\textnormal{\color{darkorange}CRC ecosystem}$$

<div align="center">
  
ChIP-seq             |  CUT&RUN
:-------------------------:|:-------------------------:
First step fragmentation | First step antibody binding
Fragmentation can be mechanical or enzymatic | Fragmentation performed by pA-MNase
Fragments collected by immunoprecipitation | Fragments collected by passive release from cells
Many cells required | Fewer cells required
More antibody volume | Less antibody volume
Enrichment signal amplified | Detects small signals and noise

</div>
