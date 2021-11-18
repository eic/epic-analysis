#!/usr/bin/env python3
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
                for line in xf: pass # set `line` to be the last line
                xsecVal = line.split()[3]
                xsecErr = line.split()[4]
                print(f'    xsec = {xsecVal} pb')
                print(f'     err = {xsecErr} pb')
                xsecDict[pathTrun]['val'].append(float(xsecVal))
                xsecDict[pathTrun]['err'].append(float(xsecErr))

print('\n')
pprint.pprint(xsecDict)
print('\n')

# compute average cross sections and tabulate
tableFileName = xsecDir+"/xsec.dat"
print(f'\nCROSS SECTION TABLE:  {tableFileName}\n'+'-'*60)
tableFileTmpName = tableFileName+".tmp"
tableFile = open(tableFileTmpName,'w+')
tableFile.write('#label cross_section_[pb] uncertainty_[pb]\n')
for path,nums in xsecDict.items():
  xsecVal = mean(nums['val'])
  xsecErr = mean(nums['err'])
  tableFile.write(f'{path} {xsecVal} {xsecErr}\n')
tableFile.close()
os.system(f'column -t {tableFileTmpName} | sort -n | tee {tableFileName}')
os.remove(tableFileTmpName)
