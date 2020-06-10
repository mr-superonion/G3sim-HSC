#!/usr/bin/env python
import os
import argparse
import numpy as np

#HOWTOUSE: runS18Sim.py --exe rgcSimBatchRun.py --minG 2 --maxG 16 --ndiv 8 --nodes 5 --procs 20

def main(exe,ndiv,nodes,procs,minG,maxG):
    nt  =   maxG-minG
    per =   nt//ndiv+1
    if nodes==1:
        qq='mini'
    if nodes>1:
        qq='small'
    if nodes>9:
        qq='large'
    for i in range(ndiv):
        begin   =   i*per+minG
        end     =   min((i+1)*per,nt)+minG
        if end<= begin:
            continue
        command =   '%s $s19w --output ./ --batch-type=pbs --nodes %d --procs %d --time 100000000 --job s19G3-%d --queue %s  --minGroup %d --maxGroup %d --clobber-versions ' %(exe,nodes,procs,i,qq,begin,end)
        print(command)
        os.system(command)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", type=str, help="execute command")
    parser.add_argument("--ndiv", type=int, help="divide reference into ndiv parts")
    parser.add_argument("--nodes", type=int, help="number of nodes")
    parser.add_argument("--procs", type=int, help="number of procs")
    parser.add_argument("--minG", type=int, help="minimum of group number")
    parser.add_argument("--maxG", type=int, help="maximum of group number")
    args = parser.parse_args()
    print('parsering arguements')
    main(args.exe,args.ndiv,args.nodes,args.procs,args.minG,args.maxG)
