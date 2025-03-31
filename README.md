# **Analyzing ChIP-seq and CUT&RUN datasets**
*YNA lab bioinformatics training series, pt 2*



## What are ChIP-seq and CUT&RUN?

ChIP-seq             |  CUT&RUN
:-------------------------:|:-------------------------:
![ChIP-seq basic steps](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/ChIP-blog-figure-1.jpg)  |  ![CUT&RUN basic steps](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/cut-run-blog-figure-1.png)

*Images from epicypher*

<br />
ChIP-seq and CUT&RUN are methods for determining the positions of histones/histone modifications or DNA binding proteins on DNA. Both methods rely on antibodies as an important initial step to select for DNA of interest. However, they differ in how their antibody-based DNA selection is performed. Key differences between ChIP-seq and CUT&RUN are: 
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
In the end, both ChIP-seq and CUT&RUN DNA libraries will be sequenced and follow a similar path for analysis. Let's comapre the analysis steps:
<br />
<br />
<div align="center">

ChIP-seq             |  CUT&RUN
:-------------------------:|:-------------------------:
Read trimming | Read trimming
Align with BWA | Align with bowtie2
Normalize to input | Normalize to scaling factor, spike-in control, or - control 
Call peaks with MACS2 | Call peaks with MACS2 and/or SEACR
Create tracks with deeptools | Create tracks with deeptools
Visualize tracks and peaks (IGV, UCSC) | Visualize tracks and peaks (IGV, UCSC)

</div>



In our last training session, our final task was to submit a CUT&RUN file for alignment. We will continue to work on that file today, but we will also grab some files to work with alongside it. 

+ a PD-NC4 ChIP-seq file
+ another CUT&RUN file we can use to test normalization

Let's request a session on the cluster and grab the files with these commands:
```
srun -t 2:00:00 --cluster htc --partition htc --spus-per-task 16 --pty bash

cd /path/to/your/working/directory

cp <file_path> <destination>
cp <file_path> <destination>
```
And now we have them in our own working directory. These are four different files that have all been aligned to the same genome assembly (hg38p.14).

Great, we have these aligned files but we can't see anything. Let's generate something we can see. We can visualize read alignment using bigWig tracks. BigWig is a type of compressed, indexed, binary format used to efficiently store and visualize genome-wide signal data, like read coverage or signal intensity, in genome browsers. It's useful because it allows a browser to call to parts of the file at a time rather than loading the very large dataset across the entire genome. So let's generate some tracks for our files.

Before we start, it's important to think about how we will normalize our files. Normalizing files, calling peaks, and creating bigWig tracks are intertwined.

'Normalizing' is a downsampling process that allows us to start with files that have equivalent sequencing coverage. Ideally, we would have the same ammount of coverage for every sample that we submit for sequencing, but in reality that is not the case. To create maningful visualizations (bigWig tracks) and robustly identify enrichment (peak calling with MACS2, SEACR) between samples, we need to start with samples that have equal sequencing coverage. Otherwise, a sample with higher coverage could appear to have more enrichment than another sample. Thus, we need to normalize our samples before we create bigWigs or call peaks. 




The point of normalizing is so that we can compare samples accurately and the way that we visualize sample data is with tracks. We also use peak calling to identify statistically enriched reagions among our samples, and we need to ensure we are performing statistical tests on samples that are normalized to each other. So rather than being a standalone process, normalization is something that is performed when we call for peaks and when we generate tracks for visualization.  

What are we even normalizing though? Technically, we are normalizing the read counts or read coverage of each sample to each other. 
