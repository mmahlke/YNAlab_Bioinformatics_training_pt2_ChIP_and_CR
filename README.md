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

+ two additional PD-NC4 CUT&RUN files
  + each file has been converted to bam format, sorted, and indexed
  + each bam file comes with it's own index (.bai)   
+ a CUT&RUN -control file (created for this analysis)
+ Each file has also been aligned to the E. coli genome assembly (more on that below)

 
To view all of the steps taken for creating/processing these files, check [here](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/CUT.RUN_training_alignments.bash).

Let's request a session on the cluster and grab the files with these commands:
```
srun -t 2:00:00 --cluster htc --partition htc --spus-per-task 16 --pty bash

#Set the path to the directory you want to work in, use the actual path to your directory
cd /path/to/your/working/directory

#Copy a text file with all the file names we want to get to our current directory '.' is shorthand for current directory location
#The syntax here is 'copy' 'path to the file you want' 'destination to copy it to'
cp /ix1/yarbely/Data/list_of_files.txt .

## We are writing a loop (or set of instructions) with 'for' and 'do'
# We create a variable named FILE, and define the variable as $(look inside ./list_of_files.txt)
# Copy each file in list_of_files.txt to the distination directory

for FILE in $(cat ./list_of_files.txt)
do
    cp ${FILE} /destination/directory
done

```
And now we have them in our own working directory. These sample files have all been aligned to the human genome assembly (hg38p.14) and the E.coli genome assembly.

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
  + if you have individualized inputs for each sample, you can use them to normalize for IP efficiency 

For $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$:
+ normalize to spike-in control
  + again, normalize for differences in library preparation and sequencing
+ normalize to read coverage
  + control for IP efficiency or other factors before library prep that may effect final coverage
+ normalize to a - control sample (like an input)
  + again, controlling for background levels
  + without a - control, you can use thresholding to set a background level

What is a spike-in control? It's a small ammount of DNA form another species that is added to each sample before library preparation. Here, we are using E. coli DNA as the spike-in control for our CUT&RUN samples. These sample files have all been aligned to the human genome assembly (hg38p.14) and the E.coli genome assembly. 

**In general, the steps for normalizing are:**
1) Align your sequencing reads to the spike-in genome
2) Count total reads aligned to spike-in genome for each sample
    +  And total reads aligned to the target genome (hg38p.14)
4) Calculate scaling factors for each sample
5) Apply scaling factors to samples
6) Use scaling factors to generate bigWigs and to call peaks


Let's assess our samples and calculate normalization ratios. First, check where we are and load our modules. 
```
pwd

module load gcc/8.2.0
module load samtools/1.14
```

We can use ```samtools stats``` command to return statistics about our bam files, including the number of mapped reads. We will perform this for our hg38p.14 alignments and our E.coli alignments.

```
samtools stats PDNC4_test_Ec_sorted.bam > PDNC4_test_stats.txt
samtools stats PDNC4_test_Ec_sorted.bam > PDNC4_test_stats.txt
samtools stats PDNC4_test_Ec_sorted.bam > PDNC4_test_stats.txt
```

Now we can view the stats reports and create a table for our samples coverage. The stats report will look like this:


Let's extract the data we need and place that in a table

Sample            |  Mapped Reads | Normalization ratio | Scaling factor
:-------------------------:|:-------------------------:|:---:|:---:|
PDNC4_test | XXXX | 0.3 |
C4Y | XXXX | .5 |
E2_12m  | XXXX | 1 |

When calculating the ratios, always set the sample with the lowest coverage to 1 and make it's # of mapped reads the denominator for scaling all other samples. It's always better to scale down your existing data than to scale up, creating non-existent arbitrary data. 

Now we have the scaling factors and we can use them to create normalized bigWig tracks with ```deeptools```. Great documentation for this package can be found here.

Let's load our modules and run the command to make a bigWig track for each of our samples.
```
module purge
module load deeptools/3.3.0

bamCoverage -b BBB_cG1_sorted.bam -o BBB_cG1.bw --scaleFactor 1 -p max/2
bamCoverage -b BBB_cG1_sorted.bam -o BBB_cG1.bw --scaleFactor 1 -p max/2
bamCoverage -b BBB_cG1_sorted.bam -o BBB_cG1.bw --scaleFactor 1 -p max/2
```
Now download those bigWigs (.bw) to your computer. We will discuss the two ways we can view them. 

First, let's view the bigWigs with IGV (Integrated Genome Viewer). IGV is available as a software you can install or as a website you can visit. Let's visit the website together. 
https://igv.org/app/  

IGV is like UCSC Genome Browser lite. It serves the same general function but doesn't offer all the same features. Most importantly, you can set up a hub for viewing tracks indefinitely that you can share to others with UCSC genome browser. We will try it out next. 

Let's load all of our bigWig files onto IGV and take a look at Chr4. You can scroll around and look at the data across the entire genome. These CUT&RUN for CENP-A should show enrichment at the centromeres and at NeoCEN4 on Chr4. Here's a snapshot of the data at NeoCEN4:


We can compare CENP-A position and see how it changes, but how can we tell if the changes are significant? We do that by calling peaks. 


