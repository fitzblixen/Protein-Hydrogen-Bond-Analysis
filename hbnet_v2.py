#!/usr/bin/python
#
# Description
# ============================================
# python routine for calculating hydrogen bond 
# network in a protein structure
#
# V1 - 2012-04-30
# ---------------
# backbone hydrogen bonds
#
# V2 - 2012-06-28
# ---------------
# include side chain hbonds 
# energy term corrections
# and new output
#
# ============================================

import fnmatch
import string
import time
import sys
import os

from collections import defaultdict
from collections import deque

#self written files
from hbNetFunct import *
from hbNetClass import *

tic = time.clock()

#sidechain Definitions      res num  a1       a2     aa     d1      d2      dd   ha hd
sSlist = {'TRP': SideChain('TRP',0 ,'non' , 'non' , 'non', 'HE1' , 'non' , 'NE1', 0, 2),
	  'ASN': SideChain('ASN',0 ,'OD1' , 'non' , 'CG' , '1HD2', '2HD2', 'ND2', 2, 2),
	  'GLN': SideChain('GLN',0 ,'OE1' , 'non' , 'CD' , '1HE2', '2HD2', 'NE2', 2, 2),
	  'ARG': SideChain('ARG',0 ,'1HH2', '2HH2', 'NH2', '1HH1', '2HH1', 'NH1', 0, 2),
	  'HIS': SideChain('HIS',0 ,'ND1' , 'non' , 'CG' , 'HE2' , 'non' , 'NE2', 2, 2),
	  'LYS': SideChain('LYS',0 ,'non' , 'non' , 'HZ1', 'HZ2' , 'HZ3' , 'NZ' , 0, 3),
	  'SER': SideChain('SER',0 ,'OG'  , 'non' , 'CB' , 'HG'  , 'non' , 'OG' , 3, 3),
	  'THR': SideChain('THR',0 ,'OG1' , 'non' , 'CB' , 'HG1' , 'non' , 'OG1', 3, 3),
	  'TYR': SideChain('TYR',0 ,'OH'  , 'non' , 'CZ' , 'HH'  , 'non' , 'OH' , 2, 2),
	  'ASP': SideChain('ASP',0 ,'OD1' , 'OD2' , 'CG' , 'non' , 'non' , 'non', 2, 0),
	  'GLU': SideChain('GLU',0 ,'OE1' , 'OE2' , 'CD' , 'non' , 'non' , 'non', 2, 0),
	  'MET': SideChain('MET',0 ,'SD'  , 'non' , 'CG' , 'non' , 'non' , 'non', 3, 0),
	  'CYS': SideChain('CYS',0 ,'SG'  , 'non' , 'CB' , 'non' , 'non' , 'non', 3, 0)}

pep = BackBone('PEP',0,'O','C','H','N')

#HB criteria
length = [2.0,3.5]
angA_sp2 = [90,180]
angA_sp3 = [60,180]
angD_sp2 = [90,180]
angD_sp3 = [90,180]

# other variables
pairlist = deque() 
enelist  = deque()
filelist = []

for file in os.listdir('.'):
	if fnmatch.fnmatch(file, '*.pdb'):
		filelist.append(file)
#
#begin main loop over snapshots
#

fCnt = 1
fTot = len(filelist)
for filename in filelist:
	
	donArray,accArray = file_read(filename,fCnt,fTot,sSlist,pep)
	
	#or indx in range(0,len(accArray)-1):
		#rint indx,accArray[indx].res,accArray[indx].reN, accArray[indx].ahyb
	#
        # MAIN LOOP, search Array/Donor Matricies for Hbonds
        #

        for indx in range(0,len(donArray)-1):
		atmD  = donArray[indx].crd1
		atmDD = donArray[indx].crd2
		
		if donArray[indx].ahyb == 2:
			angA = angA_sp2
		elif donArray[indx].ahyb == 3:
			angA = angA_sp3

		for jndx in range(0,len(accArray)-1):
			atmA  = accArray[jndx].crd1 
			atmAA = accArray[jndx].crd2
	
			if accArray[jndx].ahyb == 2:
				angD = angD_sp2
			elif accArray[jndx].ahyb == 3:
				angD = angD_sp3
			if dist_check(atmDD,atmA,length):
				if geom_check(atmD,atmA,atmAA,angA) and \
				   geom_check(atmA,atmD,atmDD,angD):
					#if plan_check():
					hbEn = ene_cal(atmD,atmDD,atmA,atmAA)
					pairStr = (donArray[indx].res + ',' +  donArray[indx].reN + ',' + donArray[indx].nam1 + ',' +
						   accArray[jndx].res + ',' +  accArray[jndx].reN + ',' + accArray[jndx].nam1)
					pairlist.append(str(pairStr))
					enelist.append(hbEn)
	fCnt +=1

#
# OUTPUT
#

pairs = defaultdict(int)
energ = defaultdict(list)

iCnt = 0 
for item in pairlist:
	pairs[item] +=1
	energ[item].append(enelist[iCnt])
	iCnt +=1

# output results
fileout = open('HBnet.dat','w')

pCnt = 0
header1 = 'Res   Atom   Res   Atom   Lifetime   Energy \n '
header2 = '=========================================== \n'
fileout.write(header1)
fileout.write(header2)
for pair in pairs.keys():
	temp = pair.split(',')
	life = float(pairs[pair])/fTot
	ener = float(sum(energ[pair])/len(energ[pair]))
	fileout.write('%4s   %4s   %4s  %4s   %8.4f %8.4f \n'%(str(temp[1]),str(temp[2]),str(temp[4]),str(temp[5]),life,ener))
	pCnt += 1
fileout.close()
toc = time.clock()
print 'TIME RUN: ',toc-tic
