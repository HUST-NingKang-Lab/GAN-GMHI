#!/usr/bin/env python3
import pandas as pd
import sys

path = 'tables'
test = pd.read_csv('%s/4347_final_relative_abundances.txt' % path, sep='\t',index_col=0).T
metadata = pd.read_csv('%s/Final_metadata_4347.csv' % path, index_col=0).T

test.insert(0,'batch',metadata['Author (year)'].tolist())
test.insert(0,'Type','')
test['Type']=pd.Series(test.index.values).str.split('_',expand=True)[0].tolist()

CD = pd.read_csv('demo/output/Crohns disease-BER.csv',index_col=1).drop(columns=['Unnamed: 0'])
CRC = pd.read_csv('demo/output/CRC-BER.csv',index_col=1).drop(columns=['Unnamed: 0'])
#H = pd.read_csv('demo/output/Healthy-BER.csv',index_col=1).drop(columns=['Unnamed: 0'])

final = pd.concat([test[test['Type']=='ACVD'],test[test['Type']=='advanced adenoma'],CRC,CD],ignore_index=False,axis=0,join='inner')
#final = pd.concat([final,H.iloc[:2635,:],test.iloc[2585:2586,:]],ignore_index=False,axis=0,join='inner')

OB = pd.read_csv('demo/output/Obesity-BER.csv',index_col=1).drop(columns=['Unnamed: 0'])
OW = pd.read_csv('demo/output/Overweight-BER.csv',index_col=1).drop(columns=['Unnamed: 0'])
T2D = pd.read_csv('demo/output/T2D-BER.csv',index_col=1).drop(columns=['Unnamed: 0'])
final = pd.concat([final,test[test['Type']=='IGT'],OB,OW,test[test['Type']=='Rheumatoid Arthritis'],test[test['Type']=='Symptomatic atherosclerosis'],T2D,test[test['Type']=='Ulcerative colitis'],test[test['Type']=='Underweight']],ignore_index=False,axis=0,join='inner')

final.to_csv(sys.argv[1])

