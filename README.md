# TRACE
HMM for DNase footprinting and motif matching

# System Requirements

## Hardware Requirements

## Software Requirements


# Installation
Clone a copy of the TRACE repository

```
git clone --recurse-submodules https://github.com/Boyle-Lab/TRACE.git
```

Build TRACE

```
make all
```

# Usage information
To call TFBSs in desired regions, TRACE requies a file of these regions, files of sequence infomation, read counts, and slopes at each position of these regions and a file containing intial model.
To generate data files in required format, run the python script:
```
dataProcessing.py -i <peak_3.file> -o <seq.file> -b <bamfile> -s <genome.size> -f <fastafile> -r <count.file> -d <slope.file> -a <span> -t
```
To perform footprinting:

```
./esthmm -Q <sequence.file> -L <slope.file> -C <counts.file> -I <init.model.file> -O <final.model.file> -P <peak_6.file> -A <output1.file> -B <output2.file> -T <thread.num>
```

## Demo
The data folder contains example data for DNase-seq on K562 cell to predict binding sites of E2F1.  For simplicity, we randomly selected 200 DNase-seq peaks. the first 3 columns in peak.bed are ..

```
./esthmm -Q ./data/sequence.txt -L ./data/slope.txt -C ./data/counts.txt -I ./data/init_hmm.txt -O ./data/E2F1_hmm.txt -P ./data/peak.bed -A ./data/output1.bed -B ./data/output2.bed -T 20
```
