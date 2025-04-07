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
cp /ix1/yarbely/Data/Training/list_of_files.txt .

## We are writing a loop (or set of instructions) with 'for' and 'do'
# We create a variable named FILE, and define the variable as $(look inside ./list_of_files.txt)
# Copy each file in list_of_files.txt to the distination directory

for FILE in $(cat ./list_of_files.txt)
do
    cp ${FILE} .
done

```
And now we have them in our own working directory. These sample files have all been aligned to the human genome assembly (hg38p.14) and the E.coli genome assembly.

Great, we have these aligned files but we can't see anything. Let's generate something we can see. We can visualize read alignment using $$\textnormal{\color{gold}bigWig}$$ tracks. $$\textnormal{\color{gold}BigWig}$$ is a type of compressed, indexed, binary format used to efficiently store and visualize genome-wide signal data, like read coverage or signal intensity, in genome browsers. It's useful because it allows a browser to access and display only parts of the file at a time rather than loading the very large dataset across the entire genome. So let's generate some $$\textnormal{\color{gold}bigWig}$$ tracks for our files.

Before we start, it's important to think about how we will **normalize** our files. Normalizing files, calling $$\textnormal{\color{violet}peaks}$$, and creating $$\textnormal{\color{gold}bigWig}$$ tracks are intertwined.

## Normalizing files for analysis
**Normalizing** is a downsampling process that allows us to start with files that have equivalent sequencing coverage. Ideally, we would have the same ammount of coverage for every sample that we submit for sequencing, but in reality that is not the case. To create maningful visualizations ($$\textnormal{\color{gold}bigWig}$$ tracks) and robustly identify enrichment ($$\textnormal{\color{violet}peak}$$ calling with MACS2, SEACR) between samples, we need to start with samples that have been adjusted to control for technical variability. Thus, we need to normalize our samples before we create $$\textnormal{\color{gold}bigWig}$$ or call $$\textnormal{\color{violet}peaks}$$. 

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

What is a spike-in control? It's a small ammount of DNA from another species that is added to each sample before library preparation. Here, we are using E. coli DNA as the spike-in control for our $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$ samples. These sample files have all been aligned to the human genome assembly (hg38p.14) and the E.coli genome assembly. 

**In general, the steps for normalizing are:**
1) Align your sequencing reads to the target genome
2) Align your sequencing reads to the spike-in genome
3) Count total reads aligned to the target and spike-in genomes
4) Calculate scaling factors for each sample
5) Use scaling factors to generate bigWigs and to call peaks


Let's assess our samples and calculate scaling factors. First, check where we are and load our modules. 
```
pwd

module load gcc/8.2.0
module load samtools/1.14
```

We can use ```samtools stats``` command to return statistics about our sam/bam files, including the number of mapped reads. We will perform this for our hg38p.14 alignments and our E.coli alignments.

```
samtools stats PDNC4_test.sam > PDNC4_test_stats.txt
samtools stats PDNC4_CA-HJ-LAP_cC4_Y.sam > PDNC4_CA-HJ-LAP_cC4_Y_stats.txt
samtools stats E2_12m_sc_D4.sam > E2_12m_sc_D4_stats.txt

samtools stats PDNC4_test_ecoli.sam > PDNC4_test_ecoli_stats.txt
samtools stats PDNC4_CA-HJ-LAP_cC4_Y_ecoli.sam > PDNC4_CA-HJ-LAP_cC4_Y_ecoli_stats.txt
samtools stats E2_12m_sc_D4_ecoli.sam > E2_12m_sc_D4_ecoli_stats.txt
```

Now we can view the stats reports and create a table for our samples coverage. The stats report will look like this:

![Samtools stats example output](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Stats_example.png)

Let's extract the data we need and place that in a table. For you own analysis, you should prepare separate normalization tables for any groups of samples you want to compare. 

Sample            |  Mapped Reads hg38p.14 | Mapped Reads E.coli | Normalization Ratio | Scaling Factor
:-------------------------:|:-------------------------:|:---:|:---:|:---:|
PDNC4_test | 16093650 | 34046 | 34046/16093650 = 0.00212 | 0.00097/0.00212 = 0.458
CA/HJ_LAP_C4Y | 23376346 | 56918 | 56918/23376346 = 0.00243 | 0.00097/0.00243 = 0.399
E2_12m_scD4  | 9628838 | 9342 | 9342/9628838 = 0.00097 | 1

When calculating the ratios, always set the sample with the lowest coverage to 1 and make it's normalization value the numerator when scaling all other samples. It's always better to scale down your existing data than to scale up because scaling up creates non-existent arbitrary data. 

When we look at the ratios and scaling factors here, what do we interpret? We can see:

+ E2_12m_scD4 encountered errors or struggled during library prep, because the ratio of E.coli spike-in to total reads is much lower than the other two samples
  + This is unlikely to reflect a biological difference, even though the sample also has less mapped reads to hg38p.14
  + Biological increase or decrease in target protein should not affect the number of E.coli reads present in the sample
+ While the ratio between E.coli and hg38p.14 reads is similar between PDNC4_test and CA/HJ_LAP_C4Y, CA/HJ_LAP_C4Y has more target protein reads
  + We expect this because it is overexpressing CENP-A
  + We see that the calculated ratio will downsample CA/HJ_LAP_C4Y the **most** to bring it to a level comparable to E2_12m_scD4
  + This will allow us to see if there is a real biological difference or just higher sequencing coverage      

Now we have the scaling factors and we can use them to create normalized $$\textnormal{\color{gold}bigWig}$$ tracks with ```deeptools```. Excellent documentation for this package can be found [here](https://deeptools.readthedocs.io/en/latest/content/list_of_tools.html).

Let's load our modules and run the command to make a $$\textnormal{\color{gold}bigWig}$$ track for each of our samples.
```
module purge
module load deeptools/3.3.0

bamCoverage -b PDNC4_test_sorted.bam -o PDNC4_test.bw --scaleFactor 0.458 -p max/2
bamCoverage -b PDNC4_CA-HJ-LAP_cC4_Y_sorted.bam -o CA-HJ-LAP_C4Y.bw --scaleFactor 0.399 -p max/2
bamCoverage -b E2_12m_sc_D4_sorted.bam -o E2_12m_sc_D4.bw --scaleFactor 1 -p max/2
```
Now download those $$\textnormal{\color{gold}bigWigs}$$ (.bw) to your computer.  

First, let's view the $$\textnormal{\color{gold}bigWigs}$$ with IGV (Integrated Genome Viewer). IGV is available as a software you can install or as a website you can visit. Let's visit the [website](https://igv.org/app/) together.  

IGV web app has support for some genomes like hg38. If you are using a unique genome, you can upload that reference to IGV before viewing your $$\textnormal{\color{gold}bigWig}$$ files. 

Let's load all of our $$\textnormal{\color{gold}bigWig}$$ files onto IGV. You can scroll around and look at the data across the entire genome. These $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$ for CENP-A should show enrichment at the centromeres and at NeoCEN4 on Chr4. Here's a snapshot of the data genome-wide:
![CENP-A genome-wide](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/IGV_snap.png)

Let's look at Chromosome 4 and let's set the scale of each track to be the same (0-60) so we can get a clearer view of the data. 

We can compare CENP-A position/abundance and see how it changes, but how can we tell if the changes are significant? We do that by calling $$\textnormal{\color{violet}peaks}$$. 

There are two different modules we will use to call $$\textnormal{\color{violet}peaks}$$. One is ```MACS2``` and the other is ```SEACR```. ```MACS2``` is a $$\textnormal{\color{violet}peak}$$ caller suitable for $$\textnormal{\color{aqua}ChIP-seq}$$ but ```SEACR``` is designed for calling $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$ $$\textnormal{\color{violet}peaks}$$.

When we call peaks with either ```MACS2``` or ```SEACR```, we can use a - control or set a user-defined threshold. For this analysis, we have a - control. 

Let's first call some peacks with ```MACS2```. Please see the documentation for MACS2 [here](https://pypi.org/project/MACS2/).
To call peaks, we need to prepare **normalized** bam files with the same scale factors we calculated above.
```
module purge
module load gcc/8.2.0
module load samtools/1.14

samtools view -h -b -s 0.458 PDNC4_test_sorted.bam > PDNC4_test_macs.bam
samtools view -h -b -s 0.399 PDNC4_CA-HJ-LAP_cC4_Y_sorted.bam > CA-HJ-LAP_cC4_Y_macs.bam

## For E2_12m_scD4, we set the scale factor to 1, so we do not need to downsample and prepare a new .bam file
```
Now we call peaks on the normalized bam files with ```MACS2```.
```
module load macs/2.2.7.1

macs2 callpeak -t PDNC4_test_macs.bam -c PDNC4_control_sorted.bam -n PDNC4_test_10-6 -f BAMPE -g 3.1e9 -q 0.0001
macs2 callpeak -t CA-HJ-LAP_cC4_Y_macs.bam -c -n PDNC4_CA-HJ-LAP_cG1_Y_sc_pdnc4_nl_10-6 -f BAMPE --nolambda -g 3.1e9 -q 0.0001
macs2 callpeak -t PDNC4_CA-HJ-LAP_cG1_Y_sc_pdnc4_macs.bam -c -n PDNC4_CA-HJ-LAP_cG1_Y_sc_pdnc4_nl_10-6 -f BAMPE --nolambda -g 3.1e9 -q 0.0001

```
Here:
+ -t is the input bam
+ -n is the name for the output prefix
+ -f is the format, here bam paired-end
+ -c is the control file
+ --nolambda inactivates dynamic lambda background detection (suggested for CUT&RUN)
+ -g is the genome size
+ -q is the q-value

The q-value represents the False Discovery Rate (FDR) and is an adjusted p-value, with a lower q-value indicating a more significant peak. A q-value of 0.05, for example, means that 5% of the called peaks are expected to be false positive. 

We have also turned off ```MACS2``` local dynamic lambda, which samples the background at each candidate peak. CUT&RUN can have a lot of background noise which can reduce peak calling with active dynamic lambda. See [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html) for more information about MACS2 lambda and sliding window peak calling algorithm.

To illustrate how peak calling parameters affect the output of called peaks, we are taking two approaches. Without a control, we are using a very stringent q-value threshold to call peaks ( 0.001% of called peaks are expected to be false positives). With a control, we are using a more permissive q-value threshold (

MACS2 will give three output files all named with the prefix you specified in the ```callpeak``` command:
+ _peaks.narrowPeak: BED6+4 format file which contains the peak locations together with peak summit, pvalue and qvalue
  + peaks.narrowPeak can be directly visualized in IGV 
+ _peaks.xls: a tabular file which contains information about called peaks including pileup and **fold enrichment**
  + the key information here is **fold enrichment** value. If you want to sort or filter your peaks based on a fold change cutoff, you can get that information here
+ _summits.bed: peak summits locations for every peak
  + the summit is where binding is highest/most likely to occur and can be used to find binding motifs

I highly recommend taking some time to familiarize yourself with different data format structures [here](https://genome.ucsc.edu/FAQ/FAQformat.html). 

Let's look at the content of our peak files. <br/>
PDNC4_test_10-6_peaks.narrowPeak ([BED6+4 format](https://genome.ucsc.edu/FAQ/FAQformat.html#format12)): 
1) chrom - chromosome name
2) chromStart - starting position in the chromosome 
3) chromEnd - ending position in the chromosome 
4) name - Name given to a region 
5) score - integer part of 9th (+3) column (-log10qvalue) multiplied by 10. int(-10*log10qvalue) Indicates how dark the peak will be displayed in the browser (0-1000)
6) strand - +/- to denote strand or orientation, or "." if no orientation is assigned
1) signalValue - Measurement of overall (average) enrichment for the region
2) pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
3) qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
4) peak - Peak summit; part of peak with highest count

PDNC4_test_10-6_peaks.xls: 
1) chromosome name
2) start position of peak
3) end position of peak
4) length of peak region
5) absolute peak summit position
6) pileup height at peak summit
7) -log10(pvalue) for the peak summit (e.g. pvalue =1e-10, then this value should be 10)
8) fold enrichment for this peak summit against random Poisson distribution with **local lambda**
9) -log10(qvalue) at peak summit

PDNC4_test_10-6_summits.bed: 
1) chromosome name
2) start position of summit
3) end position of summit
4) length of summit region
5) score - same as the narrowpeak value 5

Let's download our .narrowPeak files and upload them to IGV to look at the distribution of peaks. You should see that calling without a - control using MACS2 allows for very generous peak calling even with a stringent q-value. Calling with a - control provides more control on peak calling. 

Now let's see how ```SEACR``` performs when calling peaks on these same samples. See the documentation for ```SEACR``` [here](https://github.com/FredHutch/SEACR).

First, we need to prepare our files for peak calling with ```SEACR```.

```
module load bedtools/2.29.0

#First, we convert the bam file to a bed file

bedtools bamtobed -i PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_sorted.bam > PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr.bed

#Now we are selecting the information from the bed file that we want to keep, sorting it and sending the output to a new bed file

cut -f 1,2,3,5 PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr.bed | sort -k1,1 -k2,2n -k3,3n > PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr.clean.bed

#Next we are scaling the clean bed file using our calculated scaling factor and the size of the genome we aligned to, then sending the output to a bedgraph format
## Bedgraph format is the required file format for peak calling with SEACR 

bedtools genomecov -bg -scale 0.830699709 -i PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr.clean.bed -g /bgfs/yarbely/Data/Altemose_Oneill/PDNC4/PDNC4.hifiasm.v0.16.1_V2-2022NOV_singleLine.fasta > PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr_sc.bedgraph

#Now, repeat for the other two files
### Remember that the last file will not be scaled, so remove the -scale option

#
bedtools bamtobed -i PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_sorted.bam > PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr.bed
cut -f 1,2,3,5 PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr.bed | sort -k1,1 -k2,2n -k3,3n > PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr.clean.bed
bedtools genomecov -bg -scale 0.830699709 -i PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr.clean.bed -g /bgfs/yarbely/Data/Altemose_Oneill/PDNC4/PDNC4.hifiasm.v0.16.1_V2-2022NOV_singleLine.fasta > PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr_sc.bedgraph

#
bedtools bamtobed -i PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_sorted.bam > PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr.bed
cut -f 1,2,3,5 PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr.bed | sort -k1,1 -k2,2n -k3,3n > PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr.clean.bed
bedtools genomecov -bg -i PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr.clean.bed -g /bgfs/yarbely/Data/Altemose_Oneill/PDNC4/PDNC4.hifiasm.v0.16.1_V2-2022NOV_singleLine.fasta > PDNC4_CA-HJ-LAP_cG1_Y_pdnc4_seacr_sc.bedgraph


```
IGV is like UCSC Genome Browser lite. It serves the same general function but doesn't offer all the same features. Most importantly, you can set up a hub for viewing tracks indefinitely that you can share to others with UCSC genome browser. We will try it out next. 

