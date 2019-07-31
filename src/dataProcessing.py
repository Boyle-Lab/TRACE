import os
import sys
from pysam import AlignmentFile
import getopt
import pybedtools
from scipy.signal import savgol_filter
import numpy as np
import pandas
from scipy.stats import scoreatpercentile
from math import ceil, floor, exp
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
          GC / 2.0, "\t", (1.0 - GC) / 2.0, "\n", file = outFile)
    for i in range(len(seq)):
      print(seq[i], "\t", file = outFile)
    print("P= ", len(pos) - 1, "\n", file = outFile)
    for i in range(len(pos)):
      print(pos[i], "\n", file = outFile)
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
    #list = (np.sign(counts) * 1.0 / (1.0 + (np.exp(-(counts - perc) / std)))).tolist()
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
      ext_l = 5000 if int(peak[1]) > 5000 else int(peak[1])
      ext_r = 5000 if int(peak[2]) + 5000 < maximum else maximum - int(peak[2])
      ext_l_50 = 50 if int(peak[1]) > 50 else int(peak[1])
      ext_r_50 = 50 if int(peak[2]) + 50 < maximum else maximum - int(peak[2])
      counts = self.count_read(peak, ext_l, ext_r, shift)
      counts_raw += counts[(ext_l + 1):(len(counts) - ext_r + 1)]
      mean = np.array(counts).mean()
      std = np.array(counts).std()
      counts = [min(x, mean + 10 * std) for x in counts]
      if is_atac == 1:
        counts_bc = self.bias_correction_atac(counts, peak, ext_l, ext_r)
      elif is_atac == 2:
        counts_bc = self.bias_correction_atac_pe(counts, peak, ext_l, ext_r,
                                                 shift, -shift)
      else:
        counts_bc = self.bias_correction(counts, peak, ext_l, ext_r)

      counts = (self.within_norm(np.array(counts)).tolist())[(ext_l - ext_l_50):(len(counts) - ext_r + ext_r_50)]
      smoothed = loess_fromR(range(len(counts)), counts, 600 * span/len(counts))
      loessSignal += smoothed[(ext_l_50 + 1):(len(smoothed) - ext_r_50 + 1)]

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
      slope_2nd += self.get_slope(smoothedSignal, derivative=2)[
                   (ext_l + 1):(len(smoothedSignal) - ext_r + 1)]
      slope_1st += self.get_slope(smoothedSignal, derivative=1)[
                   (ext_l + 1):(len(smoothedSignal) - ext_r + 1)]

    slope_2nd = sigmoid([-x * 100 for x in slope_2nd])
    return loessSignal, slope_2nd, slope_1st, bc_loess


def main():
  try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:b:s:f:r:d:a:t:l:e:",
                ["ifile=", "ofile=", "bfile=", "sfile=", "ffile=", "rfile=",
                 "dfile=", "span=", "lFile=", "shift="])
  except getopt.GetoptError:
    print('dataProcessing.py -i <peak.file> -o <seq.file> -b <bamfile> '
          '-s <genome.size> -f <fastafile> -r <count.file> -d <slope.file> '
          '-a <span> -t')
    sys.exit(2)
  outputFile = False
  is_atac = False
  loessFile = False
  shift = 0
  for opt, arg in opts:
    if opt == '-h':
      print('dataProcessing.py -i <peak.file> -o <seq.file> -b <bamfile> '
            '-s <genome.size> -f <fastafile> -r <count.file> -d <slope.file> '
            '-a <span> -t')
      sys.exit()
    elif opt in ("-i", "--ifile"):
      bedFile = arg
    elif opt in ("-o", "--ofile"):
      outputFile = arg
    elif opt in ("-b", "--bfile"):
      bamFile = arg
    elif opt in ("-s", "--sfile"):
      sizeFile = arg
    elif opt in ("-f", "--ffile"):
      fastaFile = arg
    elif opt in ("-r", "--rfile"):
      countFile = arg
    elif opt in ("-d", "--dfile"):
      slopeFile = arg
    elif opt in ("-l", "--lFile"):
      loessFile = arg
    elif opt in ("-a", "--span"):
      span = float(arg)
    elif opt in ("-t", "--atac"):
      is_atac = int(arg)
    elif opt in ("-e", "--shift"):
      shift = int(arg)
  signal = Signal(bamFile, bedFile, sizeFile)
  signal.load_sequence(fastaFile, outputFile)
  signal.fastaFile = fastaFile
  loessSignal, deriv2nd, deriv1st, bc_loess = signal.get_signal(span, is_atac, shift)
  print(np.array(counts_raw).mean(), np.array(loessSignal).mean(),
        np.array(deriv2nd).mean(), np.array(bc_loess).mean())
  print(len(counts_raw), len(loessSignal), len(deriv2nd), len(deriv1st),
        len(bc_counts), len(bc_within_normed), len(bc_loess))
  if loessFile:
    s = scale(bc_loess)
    with open(loessFile, "w") as outFile:
      for i in range(len(bc_loess)):
        print(s[i], file=outFile)
  with open(countFile, "w") as outFile:
    s = scale(loessSignal)
    for i in range(len(loessSignal)):
      print(s[i], file=outFile)
  with open(slopeFile+'_2.txt', "w") as outFile:
    for i in range(len(deriv2nd)):
      print(deriv2nd[i], file=outFile)
  with open(slopeFile+'_1.txt', "w") as outFile:
    for i in range(len(deriv1st)):
      print(deriv1st[i], file=outFile)
  return

if __name__ == "__main__":
  main()
