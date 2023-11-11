#!/usr/bin/env python3
import sys

dsets=[
    'DEN1',
    'ZAFR',
    'ZBRA',
    'ZILM',
    'ZSIN',
]

for ds in dsets:
    f=open("%s.chimera-flt.csv").readlines()