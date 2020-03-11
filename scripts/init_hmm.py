#!/usr/bin/env python

##-----------------------------------------------------------------------------##
##  generate initial model for TRACE                                           ##
##                                                                             ##
##-----------------------------------------------------------------------------##

import numpy as np
import pandas
import argparse
import re
import os
from Bio import motifs

def pfm2pwm(pfm):
  pwm =[]
  for i in range(len(pfm[0])):
    row = []
    sum = int(pfm[0][i]) + int(pfm[1][i]) + int(pfm[2][i]) + int(pfm[3][i]) + 1
    for j in range(4):
      row.append((int(pfm[j][i]) + 0.25) / sum)
    pwm.append(row)
  return pwm

def convertPWM(m):
  m.pseudocounts={'A': 0.25,'C': 0.25, 'G': 0.25, 'T':0.25}
  pwm=[]
  for i in range(len(m.pwm['A'])):
    pwm.append([m.pwm['A'][i], m.pwm['C'][i], m.pwm['G'][i], m.pwm['T'][i]])
  return pwm

def read_pfm(inFile):
  pfm = []
  #lAdjust = re.compile('[ATCG]\t\[')
  p = re.compile('\d+')
  for line in inFile:
    if line.startswith(">"):
      header = line.lstrip('>').split()
    else:
      row = p.findall(line)
      pfm.append(row)
  return header, pfm


def printPWMFile(outFile, header, pwm):
  print('DE\t', header[1], '\t', header[0], "\n", file = outFile, end = '')
  for i in range(len(pwm)):
    print(i, "\t", file = outFile, end = '')
    for j in range(4):
      print(pwm[i][j], "\t", file = outFile, end = '')
    print("X\n", file = outFile, end = '')
  return


def printPWM(outFile, pwm):
  print('PWM: n=', len(pwm), file = outFile)
  for i in range(len(pwm)):
    for j in range(4):
      print(pwm[i][j], "\t", file = outFile, end='')
    print("\n", file = outFile, end='')
  return


def printMatrix(outFile, a):
  print('A:\n', file = outFile, end = '')
  for i in range(len(a)):
    for j in a[i]:
      print(j, "\t", file = outFile, end = '')
    print("\n", file = outFile, end = '')
  return


def tMatrix_table(table):
  total = len(table)
  matrix = []
  for n in range(total):
    line = [0.0 for x in range(total)]
    matrix.append(line)
  for i in range(total):
    num = 1.0 / len(table[i])
    for j in table[i]:
      matrix[i][j] = num
  return matrix, 19

def build_transition_all_top(tfLengthList):
#two sets of peak states
  transition = []
  index = 0
  doubleLengthList = []
  for i in range(len(tfLengthList)):
    doubleLengthList.append(tfLengthList[i])
    doubleLengthList.append(tfLengthList[i])
  sumLen = sum(doubleLengthList)
  for i in range(len(doubleLengthList)):
    for j in range(doubleLengthList[i]):
      if j != doubleLengthList[i] - 1:
        transition.append([index + 1])
      else:
        if (i % 2) == 0:
          transition.append([sumLen])
        else:
          transition.append([sumLen+3])
      index += 1
  bg_1 = sumLen + 6 + 2
  fp_1 = sumLen + 6 + 0
  fp_2 = sumLen + 6 + 1
  index1 = []
  index2 = []
  index = 0
  for i in range(len(tfLengthList)):
    index1.append(index)
    index += tfLengthList[i]
    index2.append(index)
    index += tfLengthList[i]
  index1.append(fp_1)
  index2.append(fp_2)
  up1 = [sumLen, sumLen + 1, bg_1]  # up sumLen
  top1 = [sumLen + 1, sumLen + 4, sumLen + 2, bg_1] # top sumLen+1
  down1 = index1 + [sumLen + 2, bg_1]  # down sumLen+2
  up2 = [sumLen + 3, sumLen + 4, bg_1]  # up sumLen+3
  top2 = [sumLen + 4, sumLen + 1, sumLen + 5, bg_1] # top sumLen+4
  down2 = index2 + [sumLen + 5, bg_1]  # down sumLen+5
  transition.append(up1)
  transition.append(top1)
  transition.append(down1)
  transition.append(up2)
  transition.append(top2)
  transition.append(down2)
  fp1 = [fp_1, sumLen, bg_1]
  fp2 = [fp_2, sumLen + 3, bg_1]
  bg1 = [sumLen, sumLen + 1, sumLen + 2, sumLen + 3, sumLen + 4, sumLen + 5, fp_1, fp_2, bg_1, bg_1+1]
  transition.append(fp1)  # fp sumLen + 6
  transition.append(fp2)
  transition.append([bg_1])  # bg sumLen + 8
  transition.append(bg1)  # bg sumLen + 9
  return transition


def eMatrix_all_top(outFile, lenList):
  extraFP = 2
  extraState = 6 + 2 + 2
  print('pi:', file = outFile)
  print('0.0 ' * (sum(lenList) + extraState - 1), file = outFile, end = '')
  print('1.0', file = outFile)
  print('mu:', file = outFile)
  # initial parameters for mean of PWM scores for each motif
  for i in range(len(lenList)):
    j = 0
    while j < i:
      print('-2.5 ' * lenList[j], file = outFile, end='')
      j += 1
    print('2.0 ' * lenList[i], file = outFile, end = '')
    j = i + 1
    while j > i and j < len(lenList):
      print('-2.5 ' * lenList[j], file = outFile, end = '')
      j += 1
    print('-5.0 ' * extraState, file = outFile)
  # initial parameters for mean of cut counts and derivatives
  for j in range(len(lenList)):
    print('-0.25 ' * int(lenList[j] / 2), file = outFile, end='')
    print('-0.0 ' * int(lenList[j] / 2), file = outFile, end='')
  print('0.4 0.5 0.4 0.0 0.0 0.0 -0.1 0.0 0.0 0.0', file = outFile)
  for j in range(len(lenList)):
    print('-0.2 ' * int(lenList[j] / 2), file = outFile, end='')
    print('0.0 ' * int(lenList[j] / 2), file = outFile, end='')
  print('0.6 0.7 0.6 0.2 0.2 0.2 -0.1 -0.1 0.0 0.0', file = outFile)

  # initial parameters for covariance matrix
  print('sigma:', file = outFile)
  for i in range(len(lenList) + 2):
    for j in range(sum(lenList) + extraState):
      print('2.0 ', file = outFile, end='')
    print("\n", file = outFile, end='')
  print('rho:', file = outFile)
  for i in range(int((len(lenList) + 2)* ((len(lenList) + 2) - 1) / 2)):
    print('0.5 ' * (sum(lenList) + extraState), file=outFile)
  return


def main():
  parser = argparse.ArgumentParser()
  # Optional parameters
  parser.add_argument("--motif-number", type=int, dest="motif_num", default=10,
                      help='number of extra motifs in model, DEFAULT: 10')
  parser.add_argument("--cluster-transfec", type=str, dest="cluster_transfec",
                      default=os.path.join(os.path.dirname(__file__),
                      '../data/JASPAR_2020_matrix_clustering_vertebrates_cluster_root_motifs.txt'),
                      help='cluster roots motifs')
  parser.add_argument("--motif-transfec", type=str, dest="motif_transfec",
                      default=os.path.join(os.path.dirname(__file__),
                      '../data/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac.txt'),
                      help='motifs from JASPAR database')
  parser.add_argument("--motif-info", type=str, dest="motif_info",
                      default=os.path.join(os.path.dirname(__file__),
                      '../data/motif_cluster_info_2020.txt'),
                      help='motifs information with names and clusters they belong')
  parser.add_argument("--cluster-order", type=str, dest="cluster_order",
                      default=os.path.join(os.path.dirname(__file__),
                      '../data/cluster_order_2020.txt'),
                      help='order of cluster motifs')
  parser.add_argument("--prefix", type=str, dest = "prefix",
                      default=os.path.join(os.path.dirname(__file__), '../data/'),
                      help="The prefix for model file.")
  parser.add_argument("--file-path", type=str, dest = "file_path",
                      default=os.path.join(os.path.dirname(__file__), '../data/'),
                      help="The file path for all input file")
  # Required input
  parser.add_argument(dest="TF", metavar="transcription factor",
                      type=str, help='transcription factor of interest')

  args = parser.parse_args()

  if not os.path.isdir(args.file_path):
    print('file path invalid:' + args.file_path)
    exit(1)
  if not os.path.isfile(args.motif_info):
    print('file path invalid:' + args.motif_info)
    exit(1)
  if not os.path.isfile(args.cluster_order):
    print('file path invalid:' + args.cluster_order)
    exit(1)
  if not os.path.isfile(args.cluster_transfec):
    print('file path invalid:' + args.cluster_transfec)
    exit(1)
  if not os.path.isfile(args.motif_transfec):
    print('file path invalid:' + args.motif_transfec)
    exit(1)
  
  with open(args.cluster_transfec) as handle:
    root_transfac = motifs.parse(handle, "TRANSFAC", strict=False)
  with open(args.motif_transfec) as handle:
    motif_transfac = motifs.parse(handle, "TRANSFAC", strict=False)
  with open(args.motif_info, "r") as inFile:
    motif_info = pandas.read_table(inFile,header=None)
  with open(args.cluster_order, "r") as inFile:
    cluster_info = pandas.read_table(inFile, header=None)
  
  lenList = []
  pwmList = []
  motifList = []
  fileList = []# a list of pwm files #
  
  motif_pwm = next((item for item in motif_transfac if item['ID']==args.TF), None)
  print(motif_pwm['ID'], motif_pwm['AC'])
  if motif_pwm==None:
    print('Error: '+ args.TF +'not exists in built-in motif list')
    return
  fileList.append(motif_pwm)
  cluster = motif_info.iloc[np.where(motif_info[0] == args.TF)].iloc[0, 2]
  rank = np.where(cluster_info[0] == cluster)[0][0]
  cluster_info = cluster_info.drop([rank])[:(args.motif_num-1)]
  for motif in cluster_info[:(args.motif_num-1)][0]:
    m = next((item for item in root_transfac if item['AC']==motif), None)
    fileList.append(m)
    print(m['ID'])
    if m==None:
      print('Error: '+ motif +'not exists in list')
      return
  for motif in fileList:
    pwm = convertPWM(motif)
    pwmList.append(pwm)
    lenList.append(len(pwm))
    
  #build transition matrix
  transition_t = build_transition_all_top(lenList)
  matrix, sumLen = tMatrix_table(transition_t)
  #fileName = os.path.dirname(__file__) + '/../data/' + args.TF + '_init_model.txt'
  fileName = args.prefix + args.TF + '_init_model.txt'
  try:
    with open(fileName, "w") as outFile:
      print('M =', len(lenList), file = outFile) # M is number of motifs in model
      print('N =', len(transition_t), file = outFile) # N is number of hidden states in model
      print('P = 3', file = outFile) # number of states in a peak
                                   # surrounding FP is 3 in our model
      print('D = 1', file = outFile) # D > 0 means there are double sets of TFBSs,
                                   # including acttive and inavtive
      printMatrix(outFile, matrix)
      for pwm in pwmList:
        printPWM(outFile, pwm)
      doubleLengthList = []
      for i in range(len(lenList)):
        doubleLengthList.append(lenList[i] * 2)
      eMatrix_all_top(outFile, doubleLengthList)
  except IOError:
    print('file path invalid:' + fileName)
  return


if __name__ == "__main__":
    main()
