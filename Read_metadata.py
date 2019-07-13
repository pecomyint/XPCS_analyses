# -*- coding: utf-8 -*-
"""
Last modified date 08/11/2018

@author: Mark and Peco
"""
#read in information for an .imm file 

from io_imm import Immfile,rdmask,ldimgs
import numpy as np
import glob, os


"""
Edit the following path if necessary
"""
DD.DD['parentfolder'] = '/projectnb/ludwiggrp/Peco Myint/Mahsa_APS Nov 2017/ludwig201711/'


#Search the folder's name where the data exists 
files = os.listdir(DD.DD['parentfolder'])
for name in files:
    if name[:4] == DD.DD['Datalabel']:
        DD.DD['Datafolder'] = name


DD.DD['pdir'] = DD.DD['parentfolder'] + DD.DD['Datafolder'] + '/'


#Search imm, batchinfo and hdf
import glob, os
os.chdir(DD.DD['pdir'])
for file in glob.glob("*.batchinfo"):
    DD.DD['batchinfo'] = file 

#Read batchinfo's metada 
f=open( DD.DD['pdir'] + DD.DD['batchinfo'] , "r")
#print(f.read())
lines = [line.strip() for line in f]
for i in range(len(lines)):
    splittedinfo = lines[i].split("=")
    splittedinfo[0] = splittedinfo[0].strip()
    splittedinfo[1] = splittedinfo[1].strip()
    splittedinfo[1] = splittedinfo[1].strip('\"')
    try:
        splittedinfo[1] = np.int_(splittedinfo[1])
    except:
        pass
    if type(splittedinfo[1]) is not np.int64:
        try:
            splittedinfo[1] = np.float_(splittedinfo[1])
        except:
            pass
    #print(splittedinfo[0].strip())
    DD.DD[splittedinfo[0]] = splittedinfo[1]



DD.DD['filename'] = DD.DD['pdir']+ DD.DD['datafilename']
DD.DD['darkname'] = DD.DD['pdir']+ DD.DD['datafilename']
DD.DD['firstfile']=DD.DD['ndata0']
DD.DD['numberfiles']=DD.DD['ndataend'] 
DD.DD['firstdark']=1
DD.DD['numberdarks']=1

#use full .imm file (not compressed).
#Read in the dark
DD.FD = Immfile(DD.DD['darkname'],DD.DD['firstdark'],DD.DD['firstdark']+DD.DD['numberdarks']-1)
#read dark images by faking no processing
DD.THRESHOLD=0.0
DD.FD.DK=0.0 #boot strap up dark
imgs=ldimgs(DD.FD)

#now setup for data file
DD.FD = Immfile(DD.DD['filename'],DD.DD['firstfile'],DD.DD['firstfile']+DD.DD['numberfiles']-1)
DD.FD.DK = np.mean(imgs,axis=0)
imgs=[]

#DD.DD['detector'] = 'Lambda'
#DD.DD['xbar']= (96-57)/.02+1232.45-1
#DD.DD['ybar']= 1245.14  #offset by 1 from yorick
#DD.DD['dpix']=.055 #millimeters

DD.DD['ccdroi']=[1,DD.FD.rows,1,DD.FD.cols]
DD.DD['rows']=DD.FD.rows
DD.DD['cols']=DD.FD.cols
DD.mask=1+0*np.ones((DD.DD['rows'], DD.DD['cols'])) #rdmask('pstorage/mask1.bin',DD.DD,pos=188,type=np.int32)