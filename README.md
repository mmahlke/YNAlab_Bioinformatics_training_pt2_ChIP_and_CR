# **Analyzing ChIP-seq and CUT&RUN datasets** 
*YNA lab bioinformatics training series, pt 2*



## What are ChIP-seq and CUT&RUN? 

ChIP-seq             |  CUT&RUN
:-------------------------:|:-------------------------:
![ChIP-seq basic steps](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/ChIP-blog-figure-1.jpg)  |  ![CUT&RUN basic steps](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/cut-run-blog-figure-1.png)

*Images from epicypher* 

<br />

$$\textnormal{\color{aqua}ChIP-seq}$$ and $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$ are methods for determining the positions of histones/histone modifications or DNA binding proteins on DNA. Both methods rely on antibodies as an important initial step to select for DNA of interest. However, they differ in how their antibody-based DNA selection is performed. Key differences between $$\textnormal{\color{aqua}ChIP-seq}$$ and $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$ are: 
<br />
<br />
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
<br />

In the end, both $$\textnormal{\color{aqua}ChIP-seq}$$ and $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$ DNA libraries will be sequenced and follow a similar path for analysis. Let's comapre the analysis steps: 
<br />
<br />

<div align="center">

ChIP-seq             |  CUT&RUN
:-------------------------:|:-------------------------:
Read trimming | Read trimming
Align with BWA | Align with bowtie2
Normalize to input | Normalize to scaling factor, spike-in control, or - control 
Create tracks with deeptools | Create tracks with deeptools
Call peaks with MACS2 | Call peaks with MACS2 and/or SEACR
Visualize tracks and peaks (IGV, UCSC) | Visualize tracks and peaks (IGV, UCSC)

</div>


## Analyzing ChIP-seq and CUT&RUN data
In our last training session, our final task was to submit a $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$ file for alignment. We will continue to work on that file today, but we will also grab some files to work with alongside it. 

+ a PD-NC4 ChIP-seq file
+ two additional PD-NC4 CUT&RUN files
+ a CUT&RUN -control file

Let's request a session on the cluster and grab the files with these commands:
```
srun -t 2:00:00 --cluster htc --partition htc --spus-per-task 16 --pty bash

cd /path/to/your/working/directory

cp <file_path> <destination>
cp <file_path> <destination>
```
And now we have them in our own working directory. These are five different files that have all been aligned to the same genome assembly (hg38p.14).

Great, we have these aligned files but we can't see anything. Let's generate something we can see. We can visualize read alignment using $$\textnormal{\color{gold}bigWig}$$ tracks. $$\textnormal{\color{gold}BigWig}$$ is a type of compressed, indexed, binary format used to efficiently store and visualize genome-wide signal data, like read coverage or signal intensity, in genome browsers. It's useful because it allows a browser to call to parts of the file at a time rather than loading the very large dataset across the entire genome. So let's generate some $$\textnormal{\color{gold}bigWig}$$ tracks for our files.

Before we start, it's important to think about how we will **normalize** our files. Normalizing files, calling $$\textnormal{\color{violet}peaks}$$, and creating $$\textnormal{\color{gold}bigWig}$$ tracks are intertwined.

## Normalizing files for analysis
**Normalizing** is a downsampling process that allows us to start with files that have equivalent sequencing coverage. Ideally, we would have the same ammount of coverage for every sample that we submit for sequencing, but in reality that is not the case. To create maningful visualizations ($$\textnormal{\color{gold}bigWig}$$ tracks) and robustly identify enrichment ($$\textnormal{\color{violet}peak}$$ calling with MACS2, SEACR) between samples, we need to start with samples that have equal sequencing coverage. Otherwise, a sample with higher coverage could appear to have more enrichment than another sample. Thus, we need to normalize our samples before we create $$\textnormal{\color{gold}bigWig}$$ or call $$\textnormal{\color{violet}peaks}$$. 

There are several ways to normalize depending on your data type and different types of normalization can be used to address different concerns. 

For $$\textnormal{\color{aqua}ChIP-seq}$$:
+ normalize to spike-in control
  + a spike-in control normalizes for differences in library preparation and sequencing outside of biological vairation between samples
+ normalize to input sample
  + an input sample is a control for background binding and tells us what part of a sample's enrichment is not due to randomness

For $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$:
+ normalize to spike-in control
  + again, normalize for differences in library preparation and sequencing
+ normalize to a - control sample (like an input)
  + again, controlling for background levels
  + without a - control, you can use thresholding to set a background level


**In general, the steps for normalizing are:**
1) Align your sequencing reads to the spike-in genome
2) Count total reads aligned to spike-in genome for each sample
    +  Optionally count total sequencing coverage*
4) Calculate scaling factors for each sample
5) Apply scaling factors to samples and control


Let's assess our samples and calculate normalization ratios.
```

