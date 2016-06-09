#!/usr/bin/env python
# -*- coding: utf-8 -*-

#================================================================#

import argparse

#================================================================#

def getArgs():
    parser = argparse.ArgumentParser(description="",version="1.0.0")
    parser.add_argument('-u',dest="up",type=argparse.FileType('r'),required=True,help='upstream csv file')
    parser.add_argument('-d',dest="down",type=argparse.FileType('r'),required=True,help='downstream csv file')
    
    args = parser.parse_args()

    return args

def splitcsv(csv):
    csvdict = {}
    
    for line in csv:
        elem = line.split(',')
        id = elem[0]
        pos = elem[1:-1]
        csvdict[id] = pos
    
    return csvdict

def main(args):
    up = splitcsv(args.up)
    down = splitcsv(args.down)
    
    for sample in sorted(up.keys()):
        
        upPos = up[sample]
        upPos.reverse()
        downPos = down[sample]
        profiles = [sample] + upPos[:-1] + downPos
        
        print ','.join(profiles)

if __name__ == '__main__':
    args = getArgs()
    main(args)