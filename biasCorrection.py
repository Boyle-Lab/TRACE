##
##  The bias correction codes and tables are borrowed and modified from
##  RGT - Regulatory Genomics Toolbox from CostaLab, which is licensed GPLv3.
##


from pysam import AlignmentFile, Samfile, Fastafile
from math import log, ceil, floor, isnan

def load_table(table_file_name_F, table_file_name_R):
  """
  Creates a bias table from a tab separated file with a k-mer and bias estimate in each line.
  Keyword arguments:
  table_file_name -- Table file name.

  Return:
  bias_table_F, bias_table_R -- Bias tables.
  """
  bias_table_F = dict()
  table_file_F = open(table_file_name_F, "r")
  for line in table_file_F:
    ll = line.strip().split("\t")
    bias_table_F[ll[0]] = float(ll[1])
  table_file_F.close()
  bias_table_R = dict()
  table_file_R = open(table_file_name_R, "r")
  for line in table_file_R:
    ll = line.strip().split("\t")
    bias_table_R[ll[0]] = float(ll[1])
  table_file_R.close()
  return [bias_table_F, bias_table_R]

def bias_correction_dnase(signal_class, signal, chrName, start, end,
                          forward_shift, reverse_shift):
  table_file_name_F = '/home/nouyang/single_hit_bias_table_F.txt'
  table_file_name_R = '/home/nouyang/single_hit_bias_table_R.txt'
  bias_table = load_table(table_file_name_F, table_file_name_R)
  if not bias_table: return signal
  # Parameters
  window = 50
  defaultKmerValue = 1.0
  genome_file_name = signal_class.fastaFile
  # Initialization
  fastaFile = Fastafile(genome_file_name)
  fBiasDict = bias_table[0]
  rBiasDict = bias_table[1]
  k_nb = len(list(fBiasDict.keys())[0])
  p1 = start
  p2 = end
  p1_w = int(p1 - (window / 2))
  p2_w = int(p2 + (window / 2))
  p1_wk = p1_w - int(floor(k_nb / 2.))
  p2_wk = p2_w + int(ceil(k_nb / 2.))
  if p1 <= 0 or p1_w <= 0 or p1_wk <= 0: return signal

  # Raw counts
  nf = [0.0] * int(p2_w - p1_w)
  nr = [0.0] * int(p2_w - p1_w)
  for read in signal_class.bam.fetch(chrName, p1_w, p2_w):
    # check if the read is unmapped, according to issue #112
    if read.is_unmapped:
      continue

    if not read.is_reverse:
      cut_site = read.pos + forward_shift
      if p1_w <= cut_site < p2_w:
        nf[cut_site - p1_w] += 1.0
    else:
      cut_site = read.aend + reverse_shift - 1
      if p1_w <= cut_site < p2_w:
        nr[cut_site - p1_w] += 1.0

  # Smoothed counts
  Nf = []
  Nr = []
  fSum = sum(nf[:window])
  rSum = sum(nr[:window])
  fLast = nf[0]
  rLast = nr[0]
  for i in range(int(window / 2), int(len(nf) - (window / 2))):
    Nf.append(fSum)
    Nr.append(rSum)
    fSum -= fLast
    fSum += nf[i + int(window / 2)]
    fLast = nf[i - int(window / 2) + 1]
    rSum -= rLast
    rSum += nr[i + int(window / 2)]
    rLast = nr[i - int(window / 2) + 1]

  # Fetching sequence
  currStr = str(fastaFile.fetch(chrName, p1_wk, p2_wk - 1)).upper()
  currRevComp = revcomp(str(fastaFile.fetch(chrName, p1_wk + 1, p2_wk)).upper())

  # Iterating on sequence to create signal
  af = []
  ar = []
  for i in range(int(ceil(k_nb / 2.)), len(currStr) - int(floor(k_nb / 2)) + 1):
    fseq = currStr[i - int(floor(k_nb / 2.)):i + int(ceil(k_nb / 2.))]
    rseq = currRevComp[
           len(currStr) - int(ceil(k_nb / 2.)) - i:len(currStr) + int(
             floor(k_nb / 2.)) - i]
    try:
      af.append(fBiasDict[fseq])
    except Exception:
      af.append(defaultKmerValue)
    try:
      ar.append(rBiasDict[rseq])
    except Exception:
      ar.append(defaultKmerValue)

  # Calculating bias and writing to wig file
  fSum = sum(af[:window])
  rSum = sum(ar[:window])
  fLast = af[0]
  rLast = ar[0]
  bias_corrected_signal = []
  for i in range(int(window / 2), int(len(af) - (window / 2))):
    nhatf = Nf[i - int(window / 2)] * (af[i] / fSum)
    nhatr = Nr[i - int(window / 2)] * (ar[i] / rSum)
    zf = log(nf[i] + 1) - log(nhatf + 1)
    zr = log(nr[i] + 1) - log(nhatr + 1)
    bias_corrected_signal.append(zf + zr)
    fSum -= fLast
    fSum += af[i + int(window / 2)]
    fLast = af[i - int(window / 2) + 1]
    rSum -= rLast
    rSum += ar[i + int(window / 2)]
    rLast = ar[i - int(window / 2) + 1]

  # Termination
  fastaFile.close()
  return bias_corrected_signal


def bias_correction_atac(signal_class, signal, chrName, start, end,
                         forward_shift, reverse_shift):
  table_file_name_F = '/home/nouyang/atac_bias_table_F.txt'
  table_file_name_R = '/home/nouyang/atac_bias_table_R.txt'
  bias_table = load_table(table_file_name_F, table_file_name_R)
  # Parameters
  window = 50
  defaultKmerValue = 1.0
  genome_file_name = signal_class.fastaFile
  # Initialization
  fastaFile = Fastafile(genome_file_name)
  fBiasDict = bias_table[0]
  rBiasDict = bias_table[1]
  k_nb = len(list(fBiasDict.keys())[0])
  p1 = start
  p2 = end
  p1_w = int(p1 - (window / 2))
  p2_w = int(p2 + (window / 2))
  p1_wk = p1_w - int(floor(k_nb / 2.))
  p2_wk = p2_w + int(ceil(k_nb / 2.))

  if (p1 <= 0 or p1_w <= 0 or p2_wk <= 0):
    # Return raw counts
    nf = [0.0] * (p2 - p1)
    nr = [0.0] * (p2 - p1)
    for read in signal_class.bam.fetch(chrName, p1, p2):
      # check if the read is unmapped, according to issue #112
      if read.is_unmapped:
        continue

      if not read.is_reverse:
        cut_site = read.pos + forward_shift
        if p1 <= cut_site < p2:
          nf[cut_site - p1] += 1.0
      else:
        cut_site = read.aend + reverse_shift - 1
        if p1 <= cut_site < p2:
          nr[cut_site - p1] += 1.0

    return nf, nr

  # Raw counts
  nf = [0.0] * (p2_w - p1_w)
  nr = [0.0] * (p2_w - p1_w)
  for read in signal_class.bam.fetch(chrName, p1_w, p2_w):
    # check if the read is unmapped, according to issue #112
    if read.is_unmapped:
      continue

    if not read.is_reverse:
      cut_site = read.pos + forward_shift
      if p1_w <= cut_site < p2_w:
        nf[cut_site - p1_w] += 1.0
    else:
      cut_site = read.aend + reverse_shift - 1
      if p1_w <= cut_site < p2_w:
        nr[cut_site - p1_w] += 1.0

  # Smoothed counts
  Nf = []
  Nr = []
  fSum = sum(nf[:window])
  rSum = sum(nr[:window])
  fLast = nf[0]
  rLast = nr[0]
  for i in range(int(window / 2), len(nf) - int(window / 2)):
    Nf.append(fSum)
    Nr.append(rSum)
    fSum -= fLast
    fSum += nf[i + int(window / 2)]
    fLast = nf[i - int(window / 2) + 1]
    rSum -= rLast
    rSum += nr[i + int(window / 2)]
    rLast = nr[i - int(window / 2) + 1]

  # Fetching sequence
  currStr = str(fastaFile.fetch(chrName, p1_wk, p2_wk - 1)).upper()
  currRevComp = revcomp(str(fastaFile.fetch(chrName, p1_wk + 1,
                        p2_wk)).upper())

  # Iterating on sequence to create signal
  af = []
  ar = []
  for i in range(int(ceil(k_nb / 2.)), len(currStr) - int(floor(k_nb / 2)) + 1):
    fseq = currStr[i - int(floor(k_nb / 2.)):i + int(ceil(k_nb / 2.))]
    rseq = currRevComp[
           len(currStr) - int(ceil(k_nb / 2.)) - i:len(currStr) + int(
             floor(k_nb / 2.)) - i]
    try:
      af.append(fBiasDict[fseq])
    except Exception:
      af.append(defaultKmerValue)
    try:
      ar.append(rBiasDict[rseq])
    except Exception:
      ar.append(defaultKmerValue)

  # Calculating bias and writing to wig file
  fSum = sum(af[:window])
  rSum = sum(ar[:window])
  fLast = af[0]
  rLast = ar[0]
  bias_corrected_signal_forward = []
  bias_corrected_signal_reverse = []
  for i in range(int(window / 2), len(af) - int(window / 2)):
    nhatf = Nf[i - int(window / 2)] * (af[i] / fSum)
    nhatr = Nr[i - int(window / 2)] * (ar[i] / rSum)
    bias_corrected_signal_forward.append(nhatf)
    bias_corrected_signal_reverse.append(nhatr)
    fSum -= fLast
    fSum += af[i + int(window / 2)]
    fLast = af[i - int(window / 2) + 1]
    rSum -= rLast
    rSum += ar[i + int(window / 2)]
    rLast = ar[i - int(window / 2) + 1]

  # Termination
  fastaFile.close()
  return bias_corrected_signal_forward, bias_corrected_signal_reverse


def bias_correction_atac_pe(signal_class, signal, chrName, start, end,
                            forward_shift, reverse_shift):
  table_file_name_F = '/home/nouyang/atac_bias_table_F.txt'
  table_file_name_R = '/home/nouyang/atac_bias_table_R.txt'
  bias_table = load_table(table_file_name_F, table_file_name_R)
  # Parameters
  window = 50
  defaultKmerValue = 1.0
  genome_file_name = signal_class.fastaFile
  # Initialization
  fastaFile = Fastafile(genome_file_name)
  fBiasDict = bias_table[0]
  rBiasDict = bias_table[1]
  k_nb = len(list(fBiasDict.keys())[0])
  p1 = start
  p2 = end
  p1_w = p1 - int(window / 2)
  p2_w = p2 + int(window / 2)
  p1_wk = p1_w - int(floor(k_nb / 2.))
  p2_wk = p2_w + int(ceil(k_nb / 2.))
  if (p1 <= 0 or p1_w <= 0 or p2_wk <= 0):
    # Return raw counts
    signal = [0.0] * (p2 - p1)
    for read in signal_class.bam.fetch(chrName, p1, p2):
      # check if the read is unmapped, according to issue #112
      if read.is_unmapped:
        continue

      if not read.is_reverse:
        cut_site = read.pos + forward_shift
        if p1 <= cut_site < p2:
          signal[cut_site - p1] += 1.0
      else:
        cut_site = read.aend + reverse_shift - 1
        if p1 <= cut_site < p2:
          signal[cut_site - p1] += 1.0

    return signal

  # Raw counts
  nf = [0.0] * (p2_w - p1_w)
  nr = [0.0] * (p2_w - p1_w)
  for read in signal_class.bam.fetch(chrName, p1_w, p2_w):
    # check if the read is unmapped, according to issue #112
    if read.is_unmapped:
      continue

    if not read.is_reverse:
      cut_site = read.pos + forward_shift
      if p1_w <= cut_site < p2_w:
        nf[cut_site - p1_w] += 1.0
    else:
      cut_site = read.aend + reverse_shift - 1
      if p1_w <= cut_site < p2_w:
        nr[cut_site - p1_w] += 1.0

  # Smoothed counts
  Nf = []
  Nr = []
  fSum = sum(nf[:window])
  rSum = sum(nr[:window])
  fLast = nf[0]
  rLast = nr[0]
  for i in range(int(window / 2), len(nf) - int(window / 2)):
    Nf.append(fSum)
    Nr.append(rSum)
    fSum -= fLast
    fSum += nf[i + int(window / 2)]
    fLast = nf[i - int(window / 2) + 1]
    rSum -= rLast
    rSum += nr[i + int(window / 2)]
    rLast = nr[i - int(window / 2) + 1]

  # Fetching sequence
  currStr = str(fastaFile.fetch(chrName, p1_wk, p2_wk - 1)).upper()
  currRevComp = revcomp(str(fastaFile.fetch(chrName, p1_wk + 1,
                        p2_wk)).upper())

  # Iterating on sequence to create signal
  af = []
  ar = []
  for i in range(int(ceil(k_nb / 2.)), len(currStr) - int(floor(k_nb / 2)) + 1):
    fseq = currStr[i - int(floor(k_nb / 2.)):i + int(ceil(k_nb / 2.))]
    rseq = currRevComp[
           len(currStr) - int(ceil(k_nb / 2.)) - i:len(currStr) + int(
             floor(k_nb / 2.)) - i]
    try:
      af.append(fBiasDict[fseq])
    except Exception:
      af.append(defaultKmerValue)
    try:
      ar.append(rBiasDict[rseq])
    except Exception:
      ar.append(defaultKmerValue)

  # Calculating bias and writing to wig file
  fSum = sum(af[:window])
  rSum = sum(ar[:window])
  fLast = af[0]
  rLast = ar[0]
  bc_signal = []
  for i in range(int(window / 2), len(af) - int(window / 2)):
    nhatf = Nf[i - int(window / 2)] * (af[i] / fSum)
    nhatr = Nr[i - int(window / 2)] * (ar[i] / rSum)
    bc_signal.append(nhatf + nhatr)
    fSum -= fLast
    fSum += af[i + int(window / 2)]
    fLast = af[i - int(window / 2) + 1]
    rSum -= rLast
    rSum += ar[i + int(window / 2)]
    rLast = ar[i - int(window / 2) + 1]

  # Termination
  fastaFile.close()
  return bc_signal




def revcomp(s):
  """Revert complement string.
  *Keyword arguments:*
      - s -- String.
  """
  revDict = dict([("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"), ("N", "N")])
  return "".join([revDict[e] for e in s[::-1]])
