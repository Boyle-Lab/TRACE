# TRACE
Transcription Factor Footprinting Using DNase I Hypersensitivity Data and DNA Sequence

# System Requirements

## Hardware Requirements

## Software Requirements


# Installation
Clone a copy of the TRACE repository:  
  
```
git clone https://github.com/Boyle-Lab/TRACE.git
```
To make sure that correct version of Python are used and all required packages are installed, we recommend using Conda and creating the environment from the environment.yml file: 
 
```
conda env create -f environment.yml
```
Build TRACE: 
  
```
make
```
 
# Usage information
To call TFBSs, TRACE requies a file of regions of interest, files of sequence infomation, read counts, and slopes at each position and a file containing intial model.    
    
### Generate required files for TRACE 
To generate data files in required format, you can use our python script dataProcessing.py.      
   
Required input:      
- `<peak_3.file>`: A file containing regions of interest. The 3 columns are chromosome number, start position and end position of regions of interest. To avoid potential errors in our main program, please make sure there are no repetitive regions.  
- `<bam.file>`: A bam file of aligned reads from DNase-seq or ATAC-seq.   
- `<genome.size>`: A genome file defining the length of each chromosome.   
- `<fasta.file>`: A sequence file in FASTA format.    
Output:   
- `<seq.file>`: A file containing sequence information of regions from <peak_3.file>, with required format. (see ./data/E2F1_seq.txt).   
- `<count.file>`: A file contains processed read counts at each position in regions from <peak_3.file>.   
- `<slope.file>`: A file contains processed deritives at each position in regions from <peak_3.file>.     
   
To generat data files using our python script:   
  
```
./src/dataProcessing.py <peak_3.file> <bam.file> <genome.size> <fasta.file> 
```
The default setting will use DNase-seq based protocol. To use ATAC-seq data instead, include ```--ATAC-seq``` argument and choose from pe (pair-end) and se (single-end). If you have preferred output directory and name, set argument `--prefix`.     
    
```
./src/dataProcessing.py <peak_3.file> <atac-seq.bam.file> <genome.size> <fasta.file> --ATAC-seq pe --prefix ./out/example
```
 
To build an intial TRACE model: 
 
```
./src/init_hmm.py <TF>
```
It will generate a starting model `<init.model.file>` for TF of your choice.  The default setting will generate a 10-motif model, to change the number of extra motifs, set argument `--motif-number`.       
  
```
./src/init_hmm.py <TF> --motif-number <N>
``` 
   
### Perform footprinting by TRACE
Besides  `<seq.file> <count.file> <slope.file> <init.model.file>`,  the main TRACE program also requires a file `<peak_3.file>` containing regions of interest. Please make sure they are the same regions that were used in data processing.
 
To perform footprinting:   
   
```
./esthmm <seq.file> <count.file> <slope.file> <init.model.file> --final-model <final.model.file> --peak-file <peak_3.file> 
```
   
`<seq.file> <count.file> <slope.file> <init.model.file>` are four required input files and they need to be in correct order. It will generate the final model `<final.model.file>`, and a output file that contains all binding sites predicton from provided regions. If `--peak-file` is not set, the program will only learn the model but will not generate binding sites predictions.    
   
If you want to apply TRACE like a motif-centric approach, set argument `--motif-file` and provide `<peak_7.file>`.   
- `<peak_7.file>`: A file containing regions of interest and motifs sites inside these regions. The first 3 columns and next 3 columns are chromosome number, start position and end position of regions of interest and motif sites, the last column is number of bases overlapping between motif sites and peaks, which can be easily obtained by bedtools intersect.   
TRACE will then generate a file containing all motif sites included in `<peak_7.file>` and their marginal posterior probabilities of being active binding sites and inactive binding sites.  
  
You can also set `--thread` and  `--max-inter` for max threads and iterations used (default is 40 and 200).   
   
```
./esthmm <seq.file> <count.file> <slope.file> <init.model.file> --final-model <final.model.file> --peak-file <peak_3.file> --motif-file <peak_7.file> --thread <N> --max-inter <N>
```
  
If you already have a trained TRACE model and only want to call binding sites based on an exsiting model, you can run decoding step directly: 
 
```
./viterbi <seq.file> <count.file> <slope.file> <final.model.file> --peak-file <peak_3.file>
```
 
## Demo
The data folder contains example data for DNase-seq on K562 cell and a initial model to start with for E2F1 binding sites prediction. For simplicity, we randomly selected 500 DNase-seq peaks in chr1.  
 
```
./esthmm ./data/E2F1_seq.txt ./data/E2F1_slope_2.txt ./data/E2F1_count.txt ./data/init_hmm.txt --final-model ./data/E2F1_hmm.txt --peak-file ./data/E2F1_peak_3.bed --motif-file ./data/E2F1_peak_7.bed
```
