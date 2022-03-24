import sys
import numpy as np
import pandas as pd

gt_in=sys.argv[1]
gt_out=sys.argv[2]


def biallele(hmp):
  # keep only biallelic sites
  hmp['alleleN']=[len(x) for x in hmp.alleles]
  hmp=hmp.loc[hmp.alleleN==3].copy() # eg. len("A/T")==3
  hmp=hmp.drop(columns=['alleleN'])
  return hmp


def gt2num(hmp):
  # recode homo-ref to 0 and homo-alt to 2 and hetero to 1
  for col in range(11,len(hmp.columns)):
    hmp.iloc[:,col]=[0 if x==2*a[0] else 2 if x==2*a[2] else 1 for (x,a) in zip(hmp.iloc[:,col], hmp.alleles)]
  return hmp
 

def num_prep(num):
  # filter out irrelevant columns
  num=num.iloc[:,np.r_[0,11:len(num.columns)]]
  num=num.set_index('rs#')
  # sample names as rows, genotypic calls as columns
  num=num.T
  return num

# drop any SNP sites with nan
# gt=gt.dropna(axis='columns').copy()


def main():
  gt=pd.read_csv(gt_in,sep="\t")
  gt=biallele(gt)
  num=gt2num(gt)
  num=num_prep(num)
  num.to_csv(gt_out)
  

if __name__ == "__main__":
  main()
