# **Analyzing ChIP-seq and CUT&RUN datasets** 
*YNA lab bioinformatics training series, pt 2*



## What are $$\textnormal{\color{aqua}ChIP-seq}$$ and $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$? 

ChIP-seq             |  CUT&RUN
:-------------------------:|:-------------------------:
![ChIP-seq basic steps](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/ChIP-blog-figure-1.jpg)  |  ![CUT&RUN basic steps](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/cut-run-blog-figure-1.png)

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
Normalize to scaling factor | Normalize to scaling factor 
Create tracks with deeptools | Create tracks with deeptools
Call peaks with MACS2 | Call peaks with MACS2 and/or SEACR
Visualize tracks and peaks (IGV, UCSC) | Visualize tracks and peaks (IGV, UCSC)

</div>


## Analyzing $$\textnormal{\color{aqua}ChIP-seq}$$ and $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$ data
In our last training session, our final task was to submit a $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$ file for alignment. We will continue to work on that file today, but we will also grab some files to work with alongside it. 

+ two additional PD-NC4 CUT&RUN files
  + each file has been converted to bam format, sorted, and indexed
  + each bam file comes with it's own index (.bai)   
+ a CUT&RUN - control file (created for this analysis)
+ Each file has also been aligned to the E. coli genome assembly (more on that below)

 
To view all of the steps taken for creating/processing these files, check [here](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/scripts/CUT.RUN_training_alignments.bash).

Before we start, let's talk a little about what we expect to see from our script that we ran in the last session. When we align .fasta/.fastq files to a genome assembly, the resulting file is a .sam format file. **Sequence Alignment/Map (SAM)** format is a TAB-delimited text format starting with an optional header denoted by '@'. If present, the header always precedes the alignment section that contains information about each read and how it aligns to the reference genome. Each alignment line has 11 mandatory fields for essential alignment information such as mapping position, and variable number of optional fields for flexible or aligner specific information.

The 11 mandatory fields are depicted below:

<div align="center">
 <img src="https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/sam_format.png" alt="sam format fields" style="width:75%; height:75%;">
</div>

**BAM** format files are a compressed, binary version of a **SAM** file. This makes them smaller and easier to handle, but also means we can't directly investigate their contents because they are not human-readable. However, BAM files (.bam) can be indexed (.bai) and used during downstream processing steps to make data handling more efficient. 

For today, I've already converted our **SAM** files to **BAM** format and generated indexes (.bai) for them. 

Let's request a session on the cluster and grab the files with these commands:
```ruby
srun -t 2:00:00 --cluster htc --partition htc --cpus-per-task 16 --pty bash

#Set the path to the directory you want to work in, use the actual path to your directory
cd /ix1/yarbely/<your_user>/training/CR_PDNC4

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
And now we have them in our own working directory. 

Great, we have these aligned files but we can't see anything. Let's generate something we can see. We can visualize read alignment using $$\textnormal{\color{gold}bigWig}$$ tracks. $$\textnormal{\color{gold}BigWig}$$ is a type of compressed, indexed, binary format used to efficiently store and visualize genome-wide signal data, like read coverage or signal intensity, in genome browsers. Genome browsers are interactive interfaces for genomic data. $$\textnormal{\color{gold}BigWig}$$ tracks useful because they allow a browser to access and display only parts of the file at a time rather than loading the very large dataset across the entire genome. So let's generate some $$\textnormal{\color{gold}bigWig}$$ tracks for our files.

Before we start, it's important to think about how we will **normalize** our files. Normalizing files, calling $$\textnormal{\color{violet}peaks}$$, and creating $$\textnormal{\color{gold}bigWig}$$ tracks are intertwined.

## Normalizing files for analysis
**Normalizing** is a downsampling process that allows us to start with files that have equivalent sequencing coverage. Ideally, we would have the same ammount of coverage for every sample after sequencing, but in reality that is not the case. To create maningful visualizations ($$\textnormal{\color{gold}bigWig}$$ tracks) and robustly identify enrichment ($$\textnormal{\color{violet}peak}$$ calling with MACS2, SEACR) between samples, we need to start with samples that have been adjusted to control for technical variability. Thus, we need to normalize our samples before we create $$\textnormal{\color{gold}bigWig}$$ or call $$\textnormal{\color{violet}peaks}$$. 

There are several ways to normalize depending on your data type and different types of normalization can be used to address different concerns. 

For $$\textnormal{\color{aqua}ChIP-seq}$$:
+ normalize to spike-in control
  + a spike-in control normalizes for differences in library preparation and sequencing outside of biological vairation between samples
+ normalize to input sample
  + an input sample is a control for background binding and tells us what part of a sample's enrichment is not due to randomness
  + input samples are generally prepared as a fraction of a sample so that the input and IP (immunoprecipitated) sample are made from the same exact starting material
  + if you have individualized inputs for each sample, you can use them to normalize for IP efficiency
+ normalize to read coverage
  + control for IP efficiency or other factors (like uneven cell count) before library prep that may effect final coverage

For $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$:
+ normalize to spike-in control
  + again, normalize for differences in library preparation and sequencing
+ normalize to a - control sample (like an input)
  + again, controlling for background levels
  + without a - control, you can use thresholding to set a background level
+ normalize to read coverage
  + again, control for IP efficiency or other factors 

What is a spike-in control? It's a small ammount of DNA from another species that is added to each sample before library preparation. Those DNA will also go through the library prep process and be amplified in the same PCR reaction as the target sample. Here, we are using E. coli DNA as the spike-in control for our $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$ samples. 

**In general, the steps for normalizing are:**
1) Align your sequencing reads to the target genome
2) Align your sequencing reads to the spike-in genome
3) Count total reads aligned to the target and spike-in genomes
4) Calculate scaling factors for each sample
5) Use scaling factors to generate $$\textnormal{\color{gold}bigWigs}$$ and to call $$\textnormal{\color{violet}peaks}$$


Let's assess our samples and calculate scaling factors. First, check where we are and load our modules. 
```bash
pwd

#We want to be in /ix1/yarbely/<your_username>/training/CR_PDNC4

module load gcc/8.2.0
module load samtools/1.14

#gcc is a compiler that translates code between human-readable and machine-readable formats. We need it to run samtools
```

We can use ```samtools stats``` command to return statistics about our sam/bam files, including the number of mapped reads. We will perform this for our hg38p.14 alignments and our E.coli alignments.

```bash
samtools stats PDNC4_test_sorted.bam > PDNC4_test_stats.txt
samtools stats PDNC4_CA-HJ-LAP_cC4_Y.bam > PDNC4_CA-HJ-LAP_cC4_Y_stats.txt
samtools stats E2_12m_sc_D4.bam > E2_12m_sc_D4_stats.txt

samtools stats PDNC4_test_ecoli.bam > PDNC4_test_ecoli_stats.txt
samtools stats PDNC4_CA-HJ-LAP_cC4_Y_ecoli.bam > PDNC4_CA-HJ-LAP_cC4_Y_ecoli_stats.txt
samtools stats E2_12m_sc_D4_ecoli.bam > E2_12m_sc_D4_ecoli_stats.txt
```

Now we can view the stats reports and create a table for our samples coverage. The stats report will look like this:

![Samtools stats example output](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/Stats_example.png)

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


## Preparing $$\textnormal{\color{gold}bigWig}$$ tracks and visualizing data

Now we have the scaling factors and we can use them to create normalized $$\textnormal{\color{gold}bigWig}$$ tracks with ```deeptools```. Excellent documentation for this package can be found [here](https://deeptools.readthedocs.io/en/latest/content/list_of_tools.html).

Let's load our modules and run the command to make a $$\textnormal{\color{gold}bigWig}$$ track for each of our samples.
```bash
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
![CENP-A genome-wide](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/IGV_snap.png)

Let's look at Chromosome 4 and let's set the scale of each track to be the same (0-60) so we can get a clearer view of the data. 

If you recall, we anticipated CA/HJ_LAP_C4Y to have more CENP-A than PDNC4_test. At CEN4, the position and abundance of CENP-A look similar between these two samples. E2_12m_scD4 actually looks to have more CENP-A reads piled up at CEN4, but remember that we scaled the other two samples down to match with scD4. Normalizing is not perfect and there are other considerations--E2_12m_scD4 is 12 months old and likely has some level of aneuploidy that can affect reads at CEN4; we know there is selective pressure to gain multiple copies of Chr4. 

If we look at NeoCEN4, we see that PDNC4_test has a strong 3 peaks CENP-A signal, with CENP-A/HJURP overexpression in CA/HJ_LAP_C4Y having a destabilizing effect on NeoCEN4. Interestingly, E2_12m_scD4 also has low CENP-A that is spreading from the 3 peaks position. 

We can compare CENP-A position/abundance and see how it changes, but how can we tell if the changes are significant? We do that by calling $$\textnormal{\color{violet}peaks}$$. 


## $$\textnormal{\color{violet}Peak}$$ calling

There are two different modules we will use to call $$\textnormal{\color{violet}peaks}$$. One is ```MACS2``` and the other is ```SEACR```. ```MACS2``` is a $$\textnormal{\color{violet}peak}$$ caller suitable for $$\textnormal{\color{aqua}ChIP-seq}$$ but ```SEACR``` is designed for calling $$\textnormal{\color{aqua}CUT}$$ & $$\textnormal{\color{aqua}RUN}$$ $$\textnormal{\color{violet}peaks}$$.

When we call peaks with either ```MACS2``` or ```SEACR```, we can use a - control or set a user-defined threshold. For this analysis, we have a - control. 

Let's first call some peacks with ```MACS2```. Please see the documentation for MACS2 [here](https://pypi.org/project/MACS2/).
To call peaks, we need to prepare **normalized** bam files with the same scale factors we calculated above.
```bash
module purge
module load gcc/8.2.0
module load samtools/1.14

samtools view -h -b -s 0.458 PDNC4_test_sorted.bam > PDNC4_test_macs.bam
samtools view -h -b -s 0.399 PDNC4_CA-HJ-LAP_cC4_Y_sorted.bam > CA-HJ-LAP_cC4_Y_macs.bam

## For E2_12m_scD4, we set the scale factor to 1, so we do not need to downsample and prepare a new .bam file
```
Now we call $$\textnormal{\color{violet}peaks}$$ on the normalized bam files with ```MACS2```.
```bash
module load macs/2.2.7.1

mkdir ./macs2_peaks

macs2 callpeak -t PDNC4_test_macs.bam -c Neg_control.bam -n PDNC4_test_ctrl_10-3 -f BAMPE --nolambda --outdir ./macs2_peaks -g 3.1e9 -q 0.001
macs2 callpeak -t PDNC4_test_macs.bam -n PDNC4_test_10-6 -f BAMPE --nolambda --outdir ./macs2_peaks -g 3.1e9 -q 0.000001

macs2 callpeak -t CA-HJ-LAP_cC4_Y_macs.bam -c Neg_control.bam -n CA-HJ-LAP_cC4_Y_ctrl_10-3 -f BAMPE --nolambda --outdir ./macs2_peaks -g 3.1e9 -q 0.001
macs2 callpeak -t CA-HJ-LAP_cC4_Y_macs.bam -n CA-HJ-LAP_cC4_Y_10-6 -f BAMPE --nolambda --outdir ./macs2_peaks -g 3.1e9 -q 0.000001

macs2 callpeak -t E2_12m_sc_D4_sorted.bam -c Neg_control.bam -n E2_12m_scD4_ctrl_10-3 -f BAMPE --nolambda --outdir ./macs2_peaks -g 3.1e9 -q 0.001
macs2 callpeak -t E2_12m_sc_D4_sorted.bam -n E2_12m_scD4_10-6 -f BAMPE --nolambda --outdir ./macs2_peaks -g 3.1e9 -q 0.000001

```
Here:
+ -t is the input bam
+ -n is the name for the output prefix
+ -f is the format, here bam paired-end
+ -c is the control file
+ --nolambda inactivates dynamic lambda background detection (suggested for CUT&RUN)
+ -g is the genome size
+ -q is the q-value

The q-value represents the False Discovery Rate (FDR) and is an adjusted p-value, with a lower q-value indicating a more significant $$\textnormal{\color{violet}peak}$$. A q-value of 0.05, for example, means that 5% of the called $$\textnormal{\color{violet}peaks}$$ are expected to be false positive. 

We have also turned off ```MACS2``` local dynamic lambda, which samples the background at each candidate $$\textnormal{\color{violet}peak}$$. CUT&RUN can have a lot of background noise which can reduce $$\textnormal{\color{violet}peak}$$ calling with active dynamic lambda. See [here](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html) for more information about MACS2 lambda and sliding window peak calling algorithm.

To illustrate how $$\textnormal{\color{violet}peak}$$ calling parameters affect the output of called $$\textnormal{\color{violet}peaks}$$, we are taking two approaches. Without a control, we are using a very stringent q-value threshold to call $$\textnormal{\color{violet}peaks}$$ ( 0.001% of called peaks are expected to be false positives). With a control, we are using a more permissive q-value threshold ( 0.1% of called peaks are expected to be false positives). 

MACS2 will give three output files all named with the prefix you specified in the ```callpeak``` command:
+ _peaks.narrowPeak: BED6+4 format file which contains the peak locations together with peak summit, pvalue and qvalue
  + peaks.narrowPeak can be directly visualized in IGV 
+ _peaks.xls: a tabular file which contains information about called peaks including pileup and **fold enrichment**
  + the key information here is **fold enrichment** value. If you want to sort or filter your peaks based on a fold change cutoff, you can get that information here
+ _summits.bed: peak summits locations for every peak
  + the summit is where binding is highest/most likely to occur and can be used to find binding motifs

I highly recommend taking some time to familiarize yourself with different data format structures [here](https://genome.ucsc.edu/FAQ/FAQformat.html). 

Let's look at the content of our $$\textnormal{\color{violet}peak}$$ files. Download the files and open them to see what the structure is. Below is a guide to what is inside each. 
<br/>

**PDNC4_test_10-6_peaks.narrowPeak ([BED6+4 format](https://genome.ucsc.edu/FAQ/FAQformat.html#format12)):**
1) chrom - chromosome name
2) chromStart - starting position in the chromosome 
3) chromEnd - ending position in the chromosome 
4) name - Name given to a region 
5) score - integer part of 9th column (-log10qvalue) multiplied by 10. int(-10*log10qvalue) Indicates how dark the peak will be displayed in the browser (0-1000)
6) strand - +/- to denote strand or orientation, or "." if no orientation is assigned
1) signalValue - Measurement of overall (average) enrichment for the region
2) pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
3) qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
4) peak - Peak summit; part of peak with highest count

**PDNC4_test_10-6_peaks.xls:**
1) chromosome name
2) start position of peak
3) end position of peak
4) length of peak region
5) absolute peak summit position
6) pileup height at peak summit
7) -log10(pvalue) for the peak summit (e.g. pvalue =1e-10, then this value should be 10)
8) fold enrichment for this peak summit against random Poisson distribution with **local lambda**
9) -log10(qvalue) at peak summit

**PDNC4_test_10-6_summits.bed:**
1) chromosome name
2) start position of summit
3) end position of summit
4) length of summit region
5) score - same as the narrowpeak value 5

Let's upload our .narrowPeak files to IGV to look at the distribution of $$\textnormal{\color{violet}peaks}$$. You should see that calling without a - control using MACS2 allows for very generous $$\textnormal{\color{violet}peak}$$ calling even with a stringent q-value. Calling with a - control provides more control on $$\textnormal{\color{violet}peak}$$ calling. 

You can also see that $$\textnormal{\color{violet}peak}$$ calling can perform differently based on the context. It behaves one way at CEN4 and another way at NeoCEN4. Keep that in mind. 

Now let's see how ```SEACR``` performs when calling $$\textnormal{\color{violet}peaks}$$ on these same samples. See the documentation for ```SEACR``` [here](https://github.com/FredHutch/SEACR).

First, we need to prepare our files for $$\textnormal{\color{violet}peak}$$ calling with ```SEACR```. Let's save some time by submitting the script found [here](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/scripts/seacr_peaks_prep.bash) and then go over the steps in the script outlined below. 

To submit the script, download the file and open it. You'll need to update the Bash instructions to include the path to your user name. In the beginning of the script, we also specify our working directory. Check and make sure you have updated that to show your path to your working directory. 

Next, upload the **updated** file to your folder at ```/ix1/yarbely/<your_user>/CR_PDNC4```. 

Then submit the script with the command 
```
sbatch seacr_peaks_prep.bash
```

Ok, let's examine what's included in the script.

```ruby
module load gcc/8.2.0
module load bedtools/2.29.0
module load samtools/1.14

cd /ix1/yarbely/<your_user>/CR_PDNC4

#Make a new directory called seacr_peaks in our current working directory
mkdir ./seacr_peaks

#First, we convert the bam file to a bed file

bedtools bamtobed -bedpe -i PDNC4_test_sorted.bam > PDNC4_test_seacr.bed

#Now we are selecting the information from the bed file that we want to keep, sorting it and sending the output to a new bed file

cut -f 1,2,3,5 PDNC4_test_seacr.bed | sort -k1,1 -k2,2n -k3,3n > PDNC4_test_seacr.clean.bed

#Next we are scaling the clean bed file using our calculated scaling factor and the size of the genome we aligned to, then sending the output to a bedgraph format
## Bedgraph format is the required file format for peak calling with SEACR 
## We need to prepare an index for our alignment file to use in the genomecov command
## We need to unzip the file to build the index

gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz

samtools faidx GCF_000001405.40_GRCh38.p14_genomic.fna.gz

#Now we can use that index to prepare the bedgraph with bedtools command genomecov

bedtools genomecov -bg -scale 0.458 -i PDNC4_test_seacr.clean.bed -g GCF_000001405.40_GRCh38.p14_genomic.fna.fai > PDNC4_test_seacr.bedgraph

#Now, repeat for the other two files
### Remember that the last file will not be scaled, so remove the -scale option

#CA-HJ-LAP_cC4
bedtools bamtobed -i PDNC4_CA-HJ-LAP_cC4_Y_sorted.bam > CA-HJ-LAP_cC4_seacr.bed
cut -f 1,2,3,5 CA-HJ-LAP_cC4_seacr.bed | sort -k1,1 -k2,2n -k3,3n > CA-HJ-LAP_cC4_seacr.clean.bed
bedtools genomecov -bg -scale 0.399 -i CA-HJ-LAP_cC4_seacr.clean.bed -g GCF_000001405.40_GRCh38.p14_genomic.fna.fai > CA-HJ-LAP_cC4_seacr.bedgraph

#E2_12m_scD4
bedtools bamtobed -i E2_12m_sc_D4_sorted.bam > E2_12m_sc_D4_seacr.bed
cut -f 1,2,3,5 E2_12m_sc_D4_seacr.bed | sort -k1,1 -k2,2n -k3,3n > E2_12m_sc_D4_seacr.clean.bed
bedtools genomecov -bg -i E2_12m_sc_D4_seacr.clean.bed -g GCF_000001405.40_GRCh38.p14_genomic.fna.fai > E2_12m_sc_D4_seacr.bedgraph

#We also need to prepare our - control 
bedtools bamtobed -i Neg_control.bam > Neg_control_seacr.bed
cut -f 1,2,3,5 Neg_control_seacr.bed | sort -k1,1 -k2,2n -k3,3n > Neg_control_seacr.clean.bed
bedtools genomecov -bg -i Neg_control.clean.bed -g GCF_000001405.40_GRCh38.p14_genomic.fna.fai > Neg_control_seacr.bedgraph

```
Now we are ready to call $$\textnormal{\color{violet}peaks}$$ with ```SEACR```. 
Let's again download a script [here](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/scripts/seacr_peak_calling.bash), update it and upload it to your current working directory.

Then submit the script with the command 
```
sbatch seacr_peak_calling.bash
```

Let's look inside the script to see what parameters we are calling peaks with:

```bash
module load seacr/1.3

#With a - control
SEACR_1.3.sh PDNC4_test_seacr.bedgraph Neg_control_seacr.bedgraph norm stringent PDNC4_test_ctrl
SEACR_1.3.sh PDNC4_test_seacr.bedgraph Neg_control_seacr.bedgraph norm relaxed PDNC4_test_ctrl

#Without a - control
SEACR_1.3.sh PDNC4_test_seacr.bedgraph 0.00001 non stringent PDNC4_test_0.00001
SEACR_1.3.sh PDNC4_test_seacr.bedgraph 0.00001 non relaxed PDNC4_test_0.00001
SEACR_1.3.sh PDNC4_test_seacr.bedgraph 0.01 non stringent PDNC4_test_0.01
SEACR_1.3.sh PDNC4_test_seacr.bedgraph 0.01 non relaxed PDNC4_test_0.01

##Repeat for other samples
```
 ```SEACR``` $$\textnormal{\color{violet}peaks}$$ are named either <output prefix>.stringent.bed OR <output prefix>.relaxed.bed depending on the parameters in the original command. ```SEACR``` $$\textnormal{\color{violet}peaks}$$ also have a simple .bed format:
1) Chromosome
2) Start coordinate
3) End coordinate
4) Total signal contained within denoted coordinates
5) Maximum bedgraph signal attained at any base pair within denoted coordinates
6) Region representing the farthest upstream and farthest downstream bases within the denoted coordinates that are represented by the maximum bedgraph signal

Now let's download our ```SEACR``` $$\textnormal{\color{violet}peaks}$$ and upload them to IGV. We can see that the $$\textnormal{\color{violet}peaks}$$ are quite different between ```SEACR``` and ```MACS2```. 

When working in IGV, you can save a snapshot of what you currently see in the browser window. Select **Save Image** and choose .svg format of .png format. Scalable vector graphics (.svg) is preferrable for publication because it retains its resolution better if you (or others) need to resize it. Portable network graphics (.png) is generally preferred for high resolution photos or high-res images that will not be resized or re-scaled. 

**Takehome:** $$\textnormal{\color{violet}peak}$$ calling is an art. You need to select the right $$\textnormal{\color{violet}peak}$$ calling strategy that works for your dataset and apply it the same way across your samples. You may need to test many different cutoffs, settings and approaches to find the right strategy. What another person did might only provide initial guidance for your final approach. 

If you are interested in more detailed approaches to peak calling and analysis, I'll be adding links to bash scripts for peak intersection/overlap and R scripts for filtering using fold change data and plotting peak data later on. Check back!! You can also apply R tools like DeSeq2 and others to analyze your peak data to help identify differential enrichment. 

## Creating track sessions with UCSC Genome Browser

The last thing we need to know how to do is to share this data interactively with others in a *pleasant* way. We can do that by setting up a track session on UCSC genome browser. 

IGV is like UCSC Genome Browser lite. It serves the same general function but doesn't offer all the same features. Most importantly, you can set up a session on UCSC genome browser for viewing tracks indefinitely that you can share to others. We will try it out next. 

You might be wondering at some point, Why the heck do I need to use this more complicated UCSC browser? I'll just use IGV, it's simple! 

**Here are some scenarios where you will need to use UCSC browser:**
+ share tracks to others without requiring others to manually load all your bigWigs and bed files and put them in an order that makes sense
+ share tracks to others while preserving colors, scaling, and other options you like for your data
  + Or for yourself! It's a lot of work to get all the settings right and then lose it when you close IGV 
+ create publicly available hubs when publishing so others can interactively look through your data


With the option to indefinitely have track sessions open and accessible comes a caveat--UCSC does not want to host your data. You can directly upload and save smaller data files like $$\textnormal{\color{violet}peak}$$ files in .bed or .narrowPeak format, but you cannot upload $$\textnormal{\color{gold}bigWig}$$ files. Instead, you need to find a remote host to store your $$\textnormal{\color{gold}bigWig}$$ files and then direct UCSC genome browser to access your files at the remote server. 

I use [Cyverse](https://cyverse.org/) to host my $$\textnormal{\color{gold}bigWig}$$ tracks for UCSC track sessions. It is free (up to 5 Gb) and I've registered with two emails to be able to host many tracks, then removed tracks that I no longer used. If you have a different free remote host you prefer, go for it! This was the first one I found and it works for me. A key feature is that the host must provide publicly-accessible links to your data so that UCSC genome browser can freely 'talk' to your track file at its location. 

First, let's upload our $$\textnormal{\color{gold}bigWigs}$$ to our Cyverse space. Navigate to the **Discovery Environment** and on the left option bar you will see an icon for your storage (a pancake stack!). Let's go here, set up a folder for this training, then upload our three $$\textnormal{\color{gold}bigWigs}$$ there. 

![Getting_Started_cyverse](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/Cyverse_upload.png)

Now let's go to the [UCSC Genome Browser](https://genome.ucsc.edu/). We will each build our own track session and save it. Let's first select an empty browser for our reference genome.  

**Select Genome browser**

![Go to UCSC and select Genome Browser](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/genome_browser.png) 

**Select the Genome Assembly you want to view from the Drop-down menu**
<br/>
We will select the genome we aligned to (hg38). You can see the initial release date and when they updated with the T2T patch (p14) highlighted in green.

![Select your assembly](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/GA_assembly_selection.png) 

**Now you should see the Genome Browser and annotation tracks**
<br/>
We can explore the browser a bit--take a look at the annotation tracks available below the browser.

![View the genome browser](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/GA_browser_view.png)



We want to add our tracks to this assembly. To do that, we select **My Data** and **Custom Tracks** to get to the **Manage Custom Tracks** page, then click on **Add Custom Tracks**.

![Navigate to Manage Custom tracks](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/GA_my_data.png)

![Add tracks](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/GA_custom_tracks.png)

![Add tracks](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/Add_browser_lines.png)

On the **Add Custom Tracks** page, we will enter a browser line that either:
+ specifies the information about and location of our $$\textnormal{\color{gold}bigWig}$$ track
+ or specifies the information about and is followed by the $$\textnormal{\color{violet}peak}$$ data included in our .bed or .narrowPeak file

To add a $$\textnormal{\color{gold}bigWig}$$ track, type a browser line like this into the 'Paste URLs or Data' box:
```
track type=bigWig name="PDNC4_test_CENP-A" bigDataUrl=<insert here a public link to your .bw file>
```
There are many options you can add to this browser line to control the way your track will appear in the browser. My most used option is adding ```color=<insert RGB color code>```. RGB color codes can be found using any online RGB color picker. 

Please look at the [Sessions User Guide](https://genome.ucsc.edu/goldenPath/help/hgSessionHelp.html) and the [Custom Tracks Guide](https://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html#CustomTracks) for more information on the available options. 

Let's use the example above to add our three $$\textnormal{\color{gold}bigWig}$$ tracks. In the space for **bigDataUrl**, we need to add the location of the file from Cyverse. We can get that link from Cyverse. 

![Getting public links](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/Cyverse_public_links.png)

Then click **Submit**. 

Repeat those steps using the links for the other two files, then navigate back to the active browser by selecting **Return to current position** to see the loaded bigWig tracks and to move around in the browser. 

Next, let's add $$\textnormal{\color{violet}peak}$$ data from a .narrowPeak file. Navigate back to the **Add Custom Tracks** page and enter a browser line like:
```
track type=narrowPeak name="PDNC4_test_peaks_macs2"
```
Press enter/return and paste the data from your .narrowPeak file below the browser line. Simply open the file, copy everything, then return to this page and paste it below the browser line. Then click **Submit**. 

You should now be able to see your tracks and this set of $$\textnormal{\color{violet}peaks}$$ in your browser window. 

**Tip:** NarrowPeak files contain the same but more data then a typical .bed file. If you want to see peaks with less information, you can convert your .narrowPeak or peaks.xls files to a simple .bed format and upload that.

Lastly, let's save this session so we can revisit it and edit it later. To save, select **My Data** then **My Sessions**. You will be prompted to log in at this point. We can log in to the YNA lab sessions management area. 
<br/>
<br/>
User: Arbelylab <br/>
Pass: 
<br/>
<br/>
![Session management](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/My_sessions.png)

Feel free to set up your own **Session Management** area if you'd like! 

You can add any session you want to this space. Always be aware of other user's sessions. Never save over an existing session unless it is your own. 

To save your session, scroll down on the **Session Management** page to where you see **Save Settings**. Type in a name for your active session and click **Submit**.

![Save your session](https://github.com/mmahlke/YNAlab_Bioinformatics_training_pt2_ChIP_and_CR/blob/main/Figures/Save_session.png)

Now you can access your session with your tracks added at any time. You will see your session listed on the **Arbely Lab Session Management** area. Remember that your tracks are linked to their location in Cyverse. If you move them around, you'll need to update your tracks browser line with a new location. 

Finally, you can save a snapshot of your UCSC browser window as well--just select **View** then select **PDF** to get options to download the current browser view and the chromosome ideogram with current location highlighted. PDF is a scalable vector format, so it's also good for publication. 

Time to celebrate!! :D



