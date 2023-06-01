#!/usr/bin/env python3

##-----------------------------------------------------------------------------##
##  Processes DNase-seq or ATAC-seq data, generate data files with required    ##
##  format for TRACE                                                           ##
##                                                                             ##
##-----------------------------------------------------------------------------##

import os
from pysam import AlignmentFile
import argparse
import pybedtools
from scipy.signal import savgol_filter
import numpy as np
import pandas
from scipy.stats import scoreatpercentile
from math import exp
from sklearn.preprocessing import scale
from rpy2.robjects import FloatVector, IntVector, globalenv
from rpy2.robjects.packages import importr
stats = importr('stats')
base = importr('base')

#Internal
from biasCorrection import bias_correction_dnase, bias_correction_atac, \
  bias_correction_atac_pe

def seq2int(seq):
  d = {'A': '0', 'C': '1', 'G': '2', 'T': '3'}
  seqInt = ''.join(d[s] if s in d else s for s in seq.upper())
  return seqInt

# a generator to returns one header and sequence at a time
def read_fasta(fileObject):
  header = ''
  seq = ''
  # skip any useless leading information
  for line in fileObject:
    if line.startswith('>'):
      header = line.strip()
      break
  for line in fileObject:
    if line.startswith('>'):
      if 'N' in seq.upper():
        yield header, seq, False
      else:
        yield header, seq, True
      header = line.strip()
      seq = ''
    else:
      seq += line.strip()
  if header:
    if 'N' in seq.upper():
      yield header, seq, False
    else:
      yield header, seq, True


def get_fasta_info(inFile):
  sumCG = 0
  sumAll = 0
  sequence = ""
  pos = []
  pos.append(1)
  k = 1
  for header, seq, bool in read_fasta(inFile):
    # keep adding C, G bases counted to get total CG number
    sumCG += seq.upper().count("C")
    sumCG += seq.upper().count("G")
    sumAll += len(seq)
    sequence += seq2int(seq.upper())
    pos.append(sumAll + 1)
    k += 1
  return (sumCG / sumAll), sequence, pos

#sigmoid transformation#
def sigmoid(x):
  neg_x = [-i for i in x]
  return ((np.exp(x)-np.exp(neg_x))/(np.exp(x)+np.exp(neg_x))).tolist()

#loess function from R, get fitted value#
def loess_fromR(x, y, f, d = 2):
  x_vector = IntVector(x)
  y_vector = FloatVector(y)
  globalenv["x_vector"] = x_vector
  globalenv["y_vector"] = y_vector
  globalenv["f"] = f
  a = round(f, 2) if round(f, 2) > 0.0 else f
  model = stats.loess('y_vector~x_vector', span = a, degree = d)
  return model.rx2('fitted')

class Signal:
  def __init__(self, bamFile, bedFile, sizeFile):
    """
    Initializes Signal.
    """
    self.bam = AlignmentFile(bamFile, "rb")
    self.bed = pybedtools.BedTool(bedFile)
    self.size = pandas.read_csv(sizeFile, sep='\t', header=None)
    self.counts = None
    self.fastaFile = None

  #get sequence information file in the correct format for TRACE#
  def load_sequence(self, fastaFile, outputFile):
    self.bed = self.bed.sequence(fi=fastaFile)
    outFile = open(outputFile, 'w')
    GC, seq, pos = get_fasta_info(open(self.bed.seqfn))
    print("T=", len(seq), "GC: ", (1.0 - GC) / 2.0, "\t", GC / 2.0, "\t",
          GC / 2.0, "\t", (1.0 - GC) / 2.0, file = outFile)
    print("\n", file = outFile, end='')
    for i in range(len(seq)):
      print(seq[i], file = outFile)
    print("\n", file = outFile, end='')
    print("P= ", len(pos) - 1, file = outFile)
    for i in range(len(pos)):
      print(pos[i], file = outFile)
    outFile.close()

  #count number of 5' end cut at each position#
  def count_read(self, region, ext_l, ext_r, shift = 0):
    reads = self.bam.fetch(reference=region[0], start=int(region[1]) - ext_l,
                           end=int(region[2]) + ext_r + 1)
    tagList = []
    for read in reads:
      if (read.is_reverse):
        tagList.append(int(read.reference_end - 1 - shift))
      else:
        tagList.append(int(read.reference_start + shift))
    list = range(int(region[1]) - ext_l + 1, int(region[2]) + ext_r + 1)
    counts = []
    for i in range(len(list)):
      counts.append(tagList.count(list[i]))
    return counts

  #normalized by 10kb surrounding window#
  def within_norm(self, counts):
    mu = np.mean(counts[np.nonzero(counts)])
    return counts / mu

  #normalized by std and percentile#
  def between_norm(self, counts, perc, std):
    list = [np.sign(x) * (1.0 / (1.0 + (exp(-(np.sign(x) * x - perc) / std)))) if x != 0.0 else 0.0 for x in counts]
    return list

  def get_slope(self, counts, window=9, derivative=1):
    return savgol_filter(np.array(counts), window, 2, deriv=derivative).tolist()


  def bias_correction(self, counts, region, ext_l, ext_r, forward_shift = 0,
                      reverse_shift = 0):
    signal = bias_correction_dnase(self, counts, region[0],
                                   int(region[1]) - ext_l,
                                   int(region[2]) + ext_r, forward_shift,
                                   reverse_shift)
    return signal

  def bias_correction_atac(self, counts, region, ext_l, ext_r,
                           forward_shift = 0, reverse_shift = 0):
    signal_f, signal_r = bias_correction_atac(self, counts, region[0],
                                              int(region[1]) - ext_l,
                                              int(region[2]) + ext_r,
                                              forward_shift, reverse_shift)
    return (np.array(signal_f)+np.array(signal_r)).tolist()

  def bias_correction_atac_pe(self, counts, region, ext_l, ext_r,
                              forward_shift = 0, reverse_shift = 0):
    signal = bias_correction_atac_pe(self, counts, region[0],
                                     int(region[1]) - ext_l,
                                     int(region[2]) + ext_r, forward_shift,
                                     reverse_shift)
    return signal

  def get_signal(self, span, is_atac, shift):
    counts_raw=[]
    loessSignal = []
    bc_normed = []
    bc_loess=[]
    slope_2nd = []
    slope_1st = []
    for peak in self.bed:
      maximum = int(self.size[1][self.size[0] == peak[0]])
      # Size of flanking region not exceeding genome size
      ext_l = 5000 if int(peak[1]) > 5000 else int(peak[1])
      ext_r = 5000 if int(peak[2]) + 5000 < maximum else maximum - int(peak[2])
      ext_l_50 = 50 if int(peak[1]) > 50 else int(peak[1])
      ext_r_50 = 50 if int(peak[2]) + 50 < maximum else maximum - int(peak[2])
      # Count aligned reads at each position
      counts = self.count_read(peak, ext_l, ext_r, shift)
      counts_raw += counts[(ext_l + 1):(len(counts) - ext_r + 1)]
      mean = np.array(counts).mean()
      std = np.array(counts).std()
      counts = [min(x, mean + 10 * std) for x in counts]
      # Bias correction
      if is_atac == "se":
        counts_bc = self.bias_correction_atac(counts, peak, ext_l, ext_r)
      elif is_atac == "pe":
        counts_bc = self.bias_correction_atac_pe(counts, peak, ext_l, ext_r,
                                                 shift, -shift)
      else:
        counts_bc = self.bias_correction(counts, peak, ext_l, ext_r)
      # Normalize read counts without bias correction by mean of surrounding 10kb
      counts = (self.within_norm(np.array(counts)).tolist())[(ext_l - ext_l_50):(len(counts) - ext_r + ext_r_50)]
      # Smooth read counts without bias correction by local regression from R
      smoothed = loess_fromR(range(len(counts)), counts, 600 * span/len(counts))
      loessSignal += smoothed[(ext_l_50 + 1):(len(smoothed) - ext_r_50 + 1)]

      # Normalize and smooth bias corrected read counts
      counts = self.within_norm(np.array(counts_bc)).tolist()
      perc = scoreatpercentile(np.array(counts), 98)
      std = np.array(counts).std()
      normedSignals = self.between_norm(counts[(ext_l - ext_l_50):(len(counts) - ext_r + ext_r_50)], perc, std)
      bc_normed += normedSignals[(ext_l_50 + 1):(len(normedSignals) - ext_r_50 + 1)]
      ext_l = ext_l_50
      ext_r = ext_r_50
      smoothedSignal = loess_fromR(range(len(normedSignals)), normedSignals,
                                   600 * span / len(normedSignals))
      bc_loess += smoothedSignal[(ext_l_50 + 1):(len(smoothedSignal) - ext_r_50 + 1)]
      # Get the first and second derivatives
      slope_2nd += self.get_slope(smoothedSignal, derivative=2)[
                   (ext_l + 1):(len(smoothedSignal) - ext_r + 1)]
      slope_1st += self.get_slope(smoothedSignal, derivative=1)[
                   (ext_l + 1):(len(smoothedSignal) - ext_r + 1)]

    slope_2nd = sigmoid([-x * 100 for x in slope_2nd])
    return loessSignal, slope_2nd, slope_1st, bc_loess


def main():
  parser = argparse.ArgumentParser()
  # Optional parameters
  parser.add_argument("--atac-seq", type=str, dest = "is_atac", default=False,
                      choices=["pe", "se"],
                      help="If set, ATAC-seq based data processing will be used. "
                           "choose between two types: pair-end(pe) and single-end(se). "
                           "DEFAULT: False")
  parser.add_argument("--span", type=float, dest="span", default=0.05,
                      help='span number for loess, DEFAULT: 0.05')
  parser.add_argument("--prefix", type=str, dest = "prefix",
                      default="TRACE",
                      help="The prefix for results files. DEFAULT: TRACE")
  parser.add_argument("--shift", type=int, dest = "shift", default=0,
                      help="Number of bases for reads to be shifted")
  parser.add_argument("--genome", type=str, dest = "genome", default="hg38",
                      help="Specify genome. Currently hg19 and hg38 are available. default:hg38")
  # Required input
  parser.add_argument(dest = "input_files", metavar="peak_3.file bam.file",
                      type=str, nargs='*', help='BED files of interesting regions, fasta.file'
                                                'BAM file of reads, sequence file')
  args = parser.parse_args()
  if (len(args.input_files) < 3):
    print("Error: missing requird files")
    parser.print_help()
    exit(1)

  genome_size = os.path.dirname(__file__) + '/../data/' + args.genome + '.chrom.sizes'
  signal = Signal(args.input_files[1], args.input_files[0], genome_size)
  # Generate sequence file
  signal.load_sequence(args.input_files[2], args.prefix + "_seq.txt")
  signal.fastaFile = args.input_files[2]
  # Process DNase-seq or ATAC-seq data
  loessSignal, deriv2nd, deriv1st, bc_loess = signal.get_signal(args.span, args.is_atac, args.shift)

  # Generate count and slope files
  with open(args.prefix + "_count.txt", "w") as outFile:
    s = scale(loessSignal)
    for i in range(len(loessSignal)):
      print(s[i], file=outFile)
  with open(args.prefix +'_slope_2.txt', "w") as outFile:
    for i in range(len(deriv2nd)):
      print(deriv2nd[i], file=outFile)
  with open(args.prefix +'_slope_1.txt', "w") as outFile:
    for i in range(len(deriv1st)):
      print(deriv1st[i], file=outFile)
  return


if __name__ == "__main__":
  main()
