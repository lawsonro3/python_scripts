# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 14:19:50 2014

@author: mlawson
"""

import sys
import os
 
bufsize = 8192
 
lines = int(sys.argv[1])
fname = sys.argv[2]
fsize = os.stat(fname).st_size
 
iter = 0
with open(sys.argv[2]) as f:
    if bufsize > fsize:
        bufsize = fsize-1
    data = []
    while True:
        iter +=1
        f.seek(fsize-bufsize*iter)
        data.extend(f.readlines())
        if len(data) >= lines or f.tell() == 0:
            print(''.join(data[-lines:]))
            break