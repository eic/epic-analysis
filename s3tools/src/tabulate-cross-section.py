#!/usr/bin/env python3

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

# make a table of cross section values, for each path to xsec files;
# given a path, the average cross section of the files will be calculated

import os, re, pprint
from numpy import mean

# main directory tree root #############
xsecDir = "datarec/xsec"
########################################

# tree to hold cross sections and errors
# - structure: { path : { 'val':[list], 'err':[list] } }
# - one list element = one file's numbers, for each file in the path
# - later we will take list averages for each path, and build the table
xsecDict = {}

# traverse directories
for (path,dirs,files) in os.walk(xsecDir):

    # truncated path, to be written to table, used as a key in `xsecDict`
    pathTrun = re.sub(xsecDir,'',path)
    pathTrun = re.sub('^/','',pathTrun)

    # loop over xsec files, assume the last line is the most accurate
    # calculation of the cross section, and add the value and error
    # to `xsecDict` in the appropriate lists
    for xsecFile in files:
        if re.search(r'\.xsec$',xsecFile):

            if not pathTrun in xsecDict:
                print(f'\nPATH = {pathTrun}')
                xsecDict[pathTrun] = {'val':[],'err':[]}

            with open(path+'/'+xsecFile) as xf:
                print(f'  FILE = {xsecFile}')
                firstLine = ""
                for line in xf:
                  if firstLine=="": firstLine=line
                  pass # sets `line` to be the last line
                xsecVal = float(line.split()[3])
                xsecErr = float(line.split()[4])
                xsecDelta = xsecVal - float(firstLine.split()[3]) # difference between first and last line's xsec
                print(f'     xsec = {xsecVal} pb')
                print(f'      err = {xsecErr} pb')
                print(f' delta/xsec = {xsecDelta/xsecVal}')
                xsecDict[pathTrun]['val'].append(xsecVal)
                xsecDict[pathTrun]['err'].append(xsecErr)

print('\n')
pprint.pprint(xsecDict)
print('\n')

# compute average cross sections and tabulate
tableFileName = xsecDir+"/xsec.dat"
print(f'\nCROSS SECTION TABLE:  {tableFileName}\n'+'-'*60)
tableFileTmpName = tableFileName+".tmp"
tableFile = open(tableFileTmpName,'w+')
tableFile.write('#label cross_section_[pb] relative_uncertainty\n')
for path,nums in xsecDict.items():
  xsecVal = mean(nums['val'])
  xsecErr = mean(nums['err'])
  tableFile.write(f'{path} {xsecVal:.5g} {xsecErr/xsecVal:.3g}\n')
tableFile.close()
os.system(f'column -t {tableFileTmpName} | sort -n | tee {tableFileName}')
os.remove(tableFileTmpName)
