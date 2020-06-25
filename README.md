# TRACE
Transcription Factor Footprinting Using DNase I Hypersensitivity Data and DNA Sequence

Read the TRACE manuscript on [bioRxiv](https://www.biorxiv.org/content/10.1101/801001v1.full).

  
## Installation

Clone a copy of the TRACE repository:  
  
```bash
$ git clone https://github.com/Boyle-Lab/TRACE.git
```
To make sure that the correct version of Python is used and all required packages are installed, we recommend using Conda and creating the environment from the `environment.yml` file: 
 
```bash
$ conda env create -f environment.yml
$ source activate TRACE_env
```
Our C program requires the GNU Scientific Library (GSL). You can download it here: [https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/).  
Build TRACE:  
  
```bash
$ make
```

## Demo
We have provided a demo containing example data for DNase-seq in K562 cells and an initial model to start with for E2F1 binding sites prediction. For simplicity, we randomly selected 500 DNase-seq peaks in chr1.  
 
```bash
$ ./scripts/init_hmm.py E2F1
$ ./scripts/dataProcessing.py ./data/E2F1_peak_3.bed ./data/ENCFF826DJP.bam ./data/hg19.fa --prefix ./data/E2F1 --genome hg19
$ ./TRACE ./data/E2F1_seq.txt ./data/E2F1_slope_2.txt ./data/E2F1_count.txt --initial-model ./data/E2F1_init_model.txt --final-model ./data/E2F1_hmm.txt --peak-file ./data/E2F1_peak_3.bed --motif-file ./data/E2F1_peak_7.bed
```
(Note: `ENCFF826DJP.bam` and `hg19.fa` files were not provided.) 
  
## Usage information

To call TFBSs, TRACE requies a file of regions of interest, files of sequence infomation, read counts, and slopes at each position and a file containing the intial model.    
    
### Generate required files for TRACE 
To generate the required files in the correct format, you can use our python script `dataProcessing.py` and `init_hmm.py`.      

```bash
$ ./scripts/dataProcessing.py <peak_3.file> <bam.file> <fasta.file> 
```
Required input:       
- `<peak_3.file>`: A file containing regions of interest in [BED3](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. The 3 columns are chromosome number, start position and end position of regions of interest. To avoid potential errors in our main program, please make sure there are no repetitive regions.  
- `<bam.file>`: A BAM file of aligned reads from DNase-seq or ATAC-seq.   
- `<fasta.file>`: A genome sequence file in FASTA format. 
 
Output:   
- `<seq.file>`: A file containing sequence information of regions from `<peak_3.file>`, with required format. (see `./data/E2F1_seq.txt`).   
- `<count.file>`: A file contains processed read counts at each position in regions from `<peak_3.file>`.   
- `<slope.file>`: A file contains processed derivatives at each position in regions from `<peak_3.file>`.     
  
You can set argument `--genome` as `hg19` or `hg38`, default is `hg38`.   
The default setting will use the DNase-seq based protocol. To use ATAC-seq data instead, include the ```--ATAC-seq``` argument and choose from `pe` (pair-end) (recommended) and `se` (single-end). If you have a preferred output directory or prefix, set argument `--prefix`.  Otherwise all the files will be saved in `./data`.   
    
```bash
$ ./scripts/dataProcessing.py <peak_3.file> <atac-seq.bam.file> <fasta.file> --ATAC-seq pe --prefix ./out/example
```
  
### Build an initial TRACE model
  
```bash
$ ./scripts/init_hmm.py <TF>
```
It will generate a starting model `<init.model.file>` for the TF of your choice.  The default setting will generate a 10-motif model. To change the number of extra motifs, set argument `--motif-number`. All PWMs, including the root motifs, are built-in. If your TF of interest is not included in [`./data/motif_cluster_info_2020.txt`](data/motif_cluster_info_2020.txt), that means it was not part of the JASPAR cluster. You will need to use your own motif file, in tranfac format as in [`./data/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac.txt`](data/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac.txt), and set `--motif-number` to 1.         
  
```bash
$ ./scripts/init_hmm.py <TF> --motif-number <N> --motif-info ../data/motif_cluster_info_2020.txt --motif-transfec ./data/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac.txt
``` 
 
For the original Boyle method, which doesn't include motif information, there is an initial, universal model available for all TFs in `./data/Boyle_model.txt`.
  
### Perform footprinting with TRACE
In addition to  `<seq.file> <count.file> <slope.file> <init.model.file>`,  the main TRACE program also requires a file `<peak_3.file>` containing the regions of interest in [BED3](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. Please make sure they are the same regions that were used in the data processing.
   
```bash
$ ./TRACE <seq.file> <count.file> <slope.file> --initial-model <init.model.file> --final-model <final.model.file> --peak-file <peak_3.file> 
```
   
   `<seq.file> <count.file> <slope.file>` are three required input files and they need to be in correct order. The training step requires an initial model file `<init.model.file>` and will generate the final model in `<final.model.file>`. If `--peak-file` is not set, the program will only learn the model and will not generate binding site predictions. If `<peak_3.file>` is provided, it will generate an output file that contains all the binding site predictons from the provided regions.       
 
If you want to run TRACE like a motif-centric method, set the argument `--motif-file` and provide a `<peak_7.file>`.   
- [`<peak_7.file>`](data/E2F1_peak_7.bed): This is a file containing regions of interest and motif sites within these regions. The first 3 columns and next 3 columns are chromosome number, start position and end position of region of interest, and motif site inside that region, the last column is the number of bases overlapping between the motif site and peak, which can be easily obtained by BEDTools intersect.   
 
TRACE will then generate a file containing all motif sites included in `<peak_7.file>` and their marginal posterior probabilities of being active binding sites and inactive binding sites.  
  
You can also set `--thread` and  `--max-inter` for max threads and iterations used (default are 40 and 200).   
   
```bash
$ ./TRACE <seq.file> <count.file> <slope.file> --initial-model <init.model.file> --final-model <final.model.file> --peak-file <peak_3.file> --motif-file <peak_7.file> --thread <N> --max-inter <N>
```
  
If you already have a trained TRACE model and only want to call binding sites based on an existing model, you can run the decoding step directly by setting `--viterbi`. 
 
### Decoding
To perform the decoding step with the trained TRACE model, use the viterbi option. Once you have a trained model, this is the only step that needs to be performed on any new open chromatin data.
 
```bash
$ ./TRACE --viterbi <seq.file> <count.file> <slope.file> --final-model <final.model.file> --peak-file <peak_3.file> --motif-file <peak_7.file> --thread <N> --max-inter <N>
```


### Interpret TRACE Output
Our demo, shown above, will generate three files: `E2F1_peak_7.bed_with_probs.txt`,  `E2F1_hmm.txt_viterbi_results.txt` and a TRACE model file `./data/E2F1_hmm.txt`.   
  
- `E2F1_peak_7.bed_with_probs.txt` contains all the provided motif sites in `./data/E2F1_peak_7.bed` followed by columns of state probabilities for all motifs included in the model, as well as generic footprints. You can only use the first two scores (fourth and fifth column), which are the probabilities of being actve binding sites or inactive binding sites, for the first motif (your TF of interest). For the assessment, we recommend using the value of bound states minus unbound states.  
   
   - `E2F1_hmm.txt_viterbi_results.txt` contains all positions in the provided peak regions `./data/E2F1_peak_3.bed`, with their assigned states and probabilities. The fourth column is the labeled states. For a 10-motif model, 1-20 represent corresponding motifs 1-10 in the model, so state 1 and 2 will be the sites that you want. The first two states represent bound and unbound binding sites, depending on their parameters. Larger state numbers are the peak states that you can ignore. The fifth and sixth columns are the probabilities of being active or inactive binding sites.   
   
`./data/model_file/` includes models trained from K562 for TFs from the JASPAR CORE vertebrates non-redundant set . They can be applied to all the other cell lines.   
`./data/prediction/` includes predictions using the K562 models. Only binding site states were included here, and we labeled state 1 and 2 as state_bound or state_unbound. The fifth column is the likelihood ratio of being an active binding site.    
      
## WDL TRACE Workflow
This pipeline is designed to chain together all the required steps for TRACE in a workflow, written in Workflow Description Language ([WDL](https://github.com/openwdl/wdl)). With the required input parameters, this automated pipeline will generate binding site predictions and the corresponding TRACE model. If you have multiple TFs of interest, you can simply run the pipeline once and WDL will run each TF in parallel. Installation of the pipeline is easy as most dependencies are automatically installed.  
 
System requirements: 
* [Java 8](https://www.java.com/en/download/) or higher. 
* [Docker CE](https://docs.docker.com/install/) 
* Python 3.4.1 or higher. 

install Python packages: 
* [Caper](https://github.com/ENCODE-DCC/caper#installation).  
* [Croo](https://github.com/ENCODE-DCC/croo#installation) 
```bash
$ pip install caper
$ pip install croo
```
Download [`TRACE.wdl`](TRACE.wdl) and [`input.json`](input.json). Then add paprameters in json file.
```js
{
    "TRACE.skipTrain": false,
    "TRACE.THREAD": 40,
    "TRACE.ITER": 200,
    "TRACE.model_size": 10,
    "TRACE.genome": "hg19",
    "TRACE.seq_file": "./data/hg19.fa",
    "TRACE.bam_file": "./data/ENCFF826DJP.bam",
    "TRACE.bam_index_file" : "./data/ENCFF826DJP.bam.bai",
    "TRACE.peak_file": "./data/E2F1_peak_3.bed",
    "TRACE.peak_motif_file": "./data/E2F1_peak_7.bed",    
    "TRACE.model_file_list": [
        
    ],
    "TRACE.motif_list": [
        "E2F1"
    ]
}
```
Parameter|Default|Description
---------|-------|-----------
`TRACE.THREAD`| 40 | Number of threads.
`TRACE.ITER`| 200 | Number of interations in learning algorithm
`TRACE.model_size` | 10 | Number of motif in TRACE model
`TRACE.genome` | hg38 | Genome, `hg19` or `hg38`
`TRACE.seq_file` | N/A | Genome sequence file in FASTA format
`TRACE.bam_file` | N/A | DNase-seq or ATAC-seq bam file
`TRACE.bam_index_file` | N/A | Index file for bam file
`TRACE.peak_file` | N/A | File of open chromatin regions, format as [`<peak_3.file>`](data/E2F1_peak_3.bed) shown [above](#perform-footprinting-by-trace)
`TRACE.peak_motif_file` | N/A | File of open chromatin regions and motif sites within each peak, format as [`<peak_7.file>`](data/E2F1_peak_7.bed) shown [above](#perform-footprinting-by-trace)
`TRACE.prefix` | N/A | Index file for bam file
`TRACE.motif_list` | N/A | List of TFs that you want to predict binding sites for. Must be in [this list](data/motif_cluster_info_2020.txt), otherwise follow [Single run](#demo) 
`TRACE.skipTrain` | false | Set to `ture` if you want to skip training step and only run viterbi step with trained models
`TRACE.model_file_list` | N/A | List of final models for each TF in motif_list, must set skipTrain to true

Run WDL workflow using `input.json`, Cromwell, and Docker backend using Caper.
```bash
$ caper run TRACE.wdl -i input.json --docker
```  
   
## Computational run times
Running time and memory cost varies depending on the size of the training data and size of the model. Having a longer training set total length, and more motifs in the model will cost more computational time and memory. Here are a few examples:  
- training step: 
 
Size of training set (kilobases)|Number of states|Computational time|Memory
-------|-------|-----------|-----------
180.6 | 34 | 1min | 0.59G
180.6 | 316 | 9min | 3.5G
883.1 | 296 | 63min | 15.6 G
1308.6 | 316 | 90min | 20.2G

- viterbi step: 
 
Size of training set (kilobases)|Number of states|Computational time|Memory
-------|-------|-----------|-----------
180.6 | 34 | <1s | 0.5G
180.6 | 316 | 7s | 2.9G
883.1 | 296 | 29s | 12.8G
1308.6 | 316 | 39s | 16.8G
6507.2 | 326 | 215s | 120G
                                                          
- CPU: Intel Xeon E5-2696 v4 @ 3.7GHz
