# TRACE
Transcription Factor Footprinting Using DNase I Hypersensitivity Data and DNA Sequence

# Installation
Clone a copy of the TRACE repository:  
  
```
git clone https://github.com/Boyle-Lab/TRACE.git
```
To make sure that correct version of Python are used and all required packages are installed, we recommend using Conda and creating the environment from the environment.yml file: 
 
```
conda env create -f environment.yml
source activate TRACE_env
```
Our C program requires GNU Scientific Library (GSL). You can download here: [https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/). 
Build TRACE: 
  
```
make
```
 
# Usage information
To call TFBSs, TRACE requies a file of regions of interest, files of sequence infomation, read counts, and slopes at each position and a file containing intial model.    
    
### Generate required files for TRACE 
To generate required files in correct format, you can use our python script dataProcessing.py and init_hmm.py.      
  
To generat data files using our python script:    
```
./scripts/dataProcessing.py <peak_3.file> <bam.file> <genome.size> <fasta.file> 
```
Required input:       
- `<peak_3.file>`: A file containing regions of interest. The 3 columns are chromosome number, start position and end position of regions of interest. To avoid potential errors in our main program, please make sure there are no repetitive regions.  
- `<bam.file>`: A bam file of aligned reads from DNase-seq or ATAC-seq.   
- `<genome.size>`: A genome file defining the length of each chromosome.   
- `<fasta.file>`: A sequence file in FASTA format.     

Output:   
- `<seq.file>`: A file containing sequence information of regions from <peak_3.file>, with required format. (see ./data/E2F1_seq.txt).   
- `<count.file>`: A file contains processed read counts at each position in regions from <peak_3.file>.   
- `<slope.file>`: A file contains processed deritives at each position in regions from <peak_3.file>.     
   
The default setting will use DNase-seq based protocol. To use ATAC-seq data instead, include ```--ATAC-seq``` argument and choose from 'pe' (pair-end) and 'se' (single-end). If you have preferred output directory and name, set argument `--prefix`.  Otherwise all files will be saved in ./data.   
    
```
./scripts/dataProcessing.py <peak_3.file> <atac-seq.bam.file> <genome.size> <fasta.file> --ATAC-seq pe --prefix ./out/example
```
  
To build an initial TRACE model: 
  
```
./scripts/init_hmm.py <TF>
```
It will generate a starting model `<init.model.file>` for TF of your choice.  The default setting will generate a 10-motif model, to change the number of extra motifs, set argument `--motif-number`. All PWMs including root motifs are built-in. If your TF of interested is not in `./data/motif_cluster_info_2020.txt`, that means your TF of interested was not included the JASPAR cluster. You will need to use your own motif file in tranfac format as in `./data/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac.txt` and set `--motif-number` to 1.         
  
```
./scripts/init_hmm.py <TF> --motif-number <N> --motif-info ../data/motif_cluster_info_2020.txt --motif-transfec ./data/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac.txt
``` 
 
For original Boyle method which doesn't include motif information, there is an initial model available and universal for all TFs in `./data/Boyle_model.txt`.
  
### Perform footprinting by TRACE
Besides  `<seq.file> <count.file> <slope.file> <init.model.file>`,  the main TRACE program also requires a file `<peak_3.file>` containing regions of interest. Please make sure they are the same regions that were used in data processing.
 
To perform footprinting:   
   
```
./TRACE <seq.file> <count.file> <slope.file> <init.model.file> --final-model <final.model.file> --peak-file <peak_3.file> 
```
   
   `<seq.file> <count.file> <slope.file>` are three required input files and they need to be in correct order. Training step requires an initial model file `<init.model.file>` and will generate the final model `<final.model.file>`. If `--peak-file` is not set, the program will only learn the model but will not generate binding sites predictions. if `<peak_3.file>` is provided, it will generate an output file that contains all binding sites predicton from provided regions.       
 
If you want to apply TRACE like a motif-centric approach, set argument `--motif-file` and provide `<peak_7.file>`.   
- `<peak_7.file>`: A file containing regions of interest and motifs sites inside these regions. The first 3 columns and next 3 columns are chromosome number, start position and end position of regions of interest and motif site inside that region, the last column is number of bases overlapping between motif sites and peaks, which can be easily obtained by bedtools intersect.   
 
TRACE will then generate a file containing all motif sites included in `<peak_7.file>` and their marginal posterior probabilities of being active binding sites and inactive binding sites.  
  
You can also set `--thread` and  `--max-inter` for max threads and iterations used (default is 40 and 200).   
   
```
./TRACE <seq.file> <count.file> <slope.file> --initial-model <init.model.file> --final-model <final.model.file> --peak-file <peak_3.file> --motif-file <peak_7.file> --thread <N> --max-inter <N>
```
  
If you already have a trained TRACE model and only want to call binding sites based on an exsiting model, you can run decoding step directly by setting `--viterbi`: 
 
```
./TRACE --viterbi <seq.file> <count.file> <slope.file> --final-model <final.model.file> --peak-file <peak_3.file> --motif-file <peak_7.file> --thread <N> --max-inter <N>
```
 
## Demo
The data folder contains example data for DNase-seq on K562 cell and a initial model to start with for E2F1 binding sites prediction. For simplicity, we randomly selected 500 DNase-seq peaks in chr1.  
 
```
./env/init_hmm.py E2F1
./env/dataProcessing.py ./data/E2F1_peak_3.bed ./data/ENCFF826DJP.bam ./data/hg19.chrom.sizes --prefix ./data/E2F1_
./TRACE ./data/E2F1_seq.txt ./data/E2F1_slope_2.txt ./data/E2F1_count.txt --initial-model ./data/E2F1_init_model.txt --final-model ./data/E2F1_hmm.txt --peak-file ./data/E2F1_peak_3.bed --motif-file ./data/E2F1_peak_7.bed
```

## Interprete TRACE’s Output
Our demo shown above will generate two files: `E2F1_peak_7.bed_with_probs.txt` and `E2F1_hmm.txt_viterbi_results.txt`.   
  
`E2F1_peak_7.bed_with_probs.txt` contains all provided motif sites followed with states probability for all motifs included in the model as well as generic footprints. You can only use the first two scores (fourth and fifth colunm) which are probabilities of being actve binding sites or inactive binding sites for the first motif (your TF of interest). For assessment, we recommend using the value of fourth colunm minus fifth colunm.  
  
`E2F1_hmm.txt_viterbi_results.txt` contains all positions in the provided peak regions, with their assigned states and probabilities. The fourth colunm is the labeled states, 1-10 represent corresponding motifs in the model, so state 1 will be the sites that you want. State numbers that are are greater than the number of motifs are the peak states that you can ignore. The fifth and sixth colunms are the probabilities of being active or inactive binding sites.
