#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sep 12 16:21:35 2019
Rev. 2023/02/15
@author: Alfonso Salinas
"""
import numpy as np
d2r= np.pi/180.0

def esf2car(cE):
    r,th,ph= cE
    x= r*np.sin(th)*np.cos(ph)
    y= r*np.sin(th)*np.sin(ph)
    z= r*np.cos(th)
    return np.array([x,y,z])

def esf2carV(fE,th,ph):
    fr,fth,fph= fE
    fC=np.zeros(3,dtype=float)
    fC[0]= np.sin(th)*np.cos(ph)*fr+np.cos(th)*np.cos(ph)*fth-np.sin(ph)*fph
    fC[1]= np.sin(th)*np.sin(ph)*fr+np.cos(th)*np.sin(ph)*fth+np.cos(ph)*fph
    fC[2]= np.cos(th)*fr-np.sin(th)*fth
    return fC

f1=open('input_Earth.txt','r')
f2=open('input_box.txt','w')

ninter=int(f1.readline())
print(ninter,file=f2)
dl=float(f1.readline())
print(dl,file=f2)
ra1=float(f1.readline())
print(ra1,file=f2)
ra2=float(f1.readline())
print(ra2,file=f2)
ncel=int(np.ceil(ra2/dl)*2)
print(ncel)
des=ncel/2*dl  
vdes=np.array([des,des,des]) 
sigmatip=int(f1.readline())
print(sigmatip,file=f2)
nali=int(f1.readline())
print(nali,file=f2)

alicoo=np.zeros((nali,3),dtype=float)
for i in range(nali):
    coo=np.array(list(map(float,(f1.readline()).split())))
    alicoo[i,0]=coo[0]+ra1
    alicoo[i,1:3]=coo[1:3]*d2r
    coocar=esf2car(alicoo[i,:])+vdes
    print(*coocar,file=f2)
for i in range(nali):
  alitipleido=list(map(float,(f1.readline()).split()))
  alitipesf=esf2carV(alitipleido,*alicoo[i,1:3])
  print(*alitipesf,file=f2)
nsal=int(f1.readline())
print(nsal,file=f2)
salcoo=np.zeros((nsal,3),dtype=float)
for i in range(nsal):
    coo=np.array(list(map(float,(f1.readline()).split())))
    salcoo[i,0]=coo[0]+ra1
    salcoo[i,1:3]=coo[1:3]*d2r
    coocar=esf2car(salcoo[i,:])+vdes
    print(*coocar,file=f2)
f1.close()
f2.close() 

