#!/usr/bin/env python
import os
import numpy as np
from astropy.table import Table,vstack 

def getInfo(fname):
    tractList   =   []
    fieldList   =   []
    while True:
        field   =   input('please input field name: \n')
        if not len(field)>0:
            break
        while True:
            tractMin=   input('please input tractMin: \n')
            tractMax=   input('please input tractMax: \n')
            if not len(tractMin)>0:
                break
            tractMin=   eval(tractMin)
            tractMax=   eval(tractMax)
            for i in range(tractMin,tractMax+1):
                tractList.append(i)
                fieldList.append(field)

    tableOut=   Table([tractList,fieldList],names=['tract','field'])
    tableOut.write(fname)
    return tableOut
def main():
    fname   =  '/work/xiangchong.li/work/S18ASimHSC/s18/fieldTractInfo.csv'
    if not os.path.exists(fname):
        t   =   getInfo(fname)
    else:
        t   =   Table.read(fname)
    to      =   Table.read('info.csv')
    to['field']=None
    for row in to:
        row['field']    =   t[t['tract']==row['tract']]['field'][0]
    print(to)
    to.write('info2.csv')
    return

if __name__=="__main__":
     main()
