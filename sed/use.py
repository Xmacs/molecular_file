import sys  
import os
import linecache
import pandas as pd
import numpy as np
import glob
import math
import string


for i in ['200']:
    it=int(i)
    nt=it+9
    for j in ['5']:
        frame=4001
        #fp=open('%s-%s.lammpstrj' % (i,j),'r')
        filename=("%s-%s.lammpstrj" % (i,j))
        fo=open ('%s-%s.csv' % (i,j),'w')
        for loop in range (0,frame):
            atomheadline=linecache.getline(filename,nt*loop+10)
            atomendline=linecache.getline(filename,nt*loop+nt)
            t=""
            t=linecache.getline(filename,nt*loop + 2 )
            te=0
            te=int(t)
            atomheadx=""
            atomheady=""
            atomendx=""
            atomendy=""
            atomheadx=atomheadline.split()[2]
            atomheady=atomheadline.split()[3]
            atomendx=atomendline.split()[2]
            atomendy=atomendline.split()[3]
            headx=0
            heady=0
            endx=0
            endy=0
            headx=float(atomheadx)
            heady=float(atomheady)
            endx=float(atomendx)
            endy=float(atomendy)
            x1=headx-endx
            y1=heady-endy
            x2=math.pow(x1,2)
            y2=math.pow(y1,2)
            e=math.sqrt(x2+y2)
            #print("%s,%f" % (te,e),file=filenameoutput)
            fo.write("%d,%f\n" % (te,e))
            linecache.clearcache()
        fo.close

fi=open('ave.csv','w')
for k in ['200']:
    for q in ['5']:
        tlen=0
        num=0
        fh=open('%s-%s.csv' % (k,q),'r')
        fh.readline()
        ave=0
        for line in fh:
            data=line.split(",")
            tlen+=float(data[1])
            num+=1
        ave=tlen/float(num)
        fi.write("%s-%s,%f\n" %(k,q,ave))
        fh.close()
fi.close()


