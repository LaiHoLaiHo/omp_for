#!/hetghome/hetgsoft/anaconda3/bin/python3.9
#sed -i 's/\r//g' RK45_w_tur_C_fin.py
#chmod 777 RK45_w_tur_C_fin.py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math as m
import math
from scipy.optimize import fsolve
from scipy import linalg
import csv
import pandas as pd
import sys
import numpy as np
import os
from scipy.integrate import RK45
from datetime import date
import time
from scipy.interpolate import interp1d
from os.path import exists
#####################################################
from ctypes import *
import ctypes
lib_file = "/hetghome/jordan/NeuOsc/Cfin_repro_2007_f1_right/cf_genmat_para_fin.so"
#lib_file = "/hetghome/jordan/NeuOsc/C_repro_2007_f1_right/cf_genmat_para.so"
cfuns = CDLL(lib_file)
#####################################################
cur_t = time.strftime("%H-%M-%S",time.localtime())
today = date.today()
date = today.strftime("%b-%d-%Y")
print(str(date)+"-"+str(cur_t))
cur_dat = str(date)+"-"+str(cur_t)
##################################################################
z0 = 0 #set for the initial z
t1=3 #final z
ts = 0.01 #maxstep
acu = 1e-5
test = False
step_num = 500
test_run = False
print("acu ", acu)
print("t1 ",t1)
#---------------------------------------------------
ar = 0.7
omg = 1
omga = -omg
Omg = 0
#--------------------
theta = m.pi/3
v = 1
#----------------------
vxp=m.sin(theta)*v
vxn=m.sin(theta)*(-v)
vz=m.cos(theta)*v
#-----------------------------------------------------------------
cd_ar = c_double(ar)
cd_omg = c_double(omg)
cd_omga = c_double(omga)
cd_Omg = c_double(Omg)
cd_vxp = c_double(vxp)
cd_vxn = c_double(vxn)
cd_vz = c_double(vz)
#-------------------------------
ms = cfuns.rms(c_int(-99))
dk = cfuns.rdk(c_int(-99))
hk = cfuns.rhk(c_int(-99))
nk = cfuns.rnk(c_int(-99))
vConsMu = cfuns.rConsMu(c_int(-99))
print("dk",dk)
if vConsMu == 0:
    ConsMu = False
else:
    ConsMu = True
##############################################################################################################
class  DoubleArrayType : 
    def  from_param ( self ,  param ): 
        typename  =  type ( param ) . __name__ 
        if  hasattr ( self ,  'from_'  +  typename ): 
            return  getattr ( self ,  'from_'  +  typename )( param ) 
        elif  isinstance ( param ,  ctypes. Array ): 
            return  param 
        else : 
            raise  TypeError ( "Can't convert %s "  %  typename )

    # Cast from array.array objects 
    def  from_array ( self ,  param ): 
        if  param . typecode  !=  'd' : 
            raise  TypeError ( 'must be an array of doubles' ) 
        ptr ,  _  =  param . buffer_info () 
        return  ctypes . cast ( ptr ,  ctypes . POINTER ( ctypes . c_double ))

    # Cast from lists/tuples 
    def  from_list ( self ,  param ): 
        val  =  (( ctypes . c_double ) * len ( param ))( * param ) 
        return  val

    from_tuple  =  from_list

    # Cast from a numpy array 
    #this is having some problems so there I have use list as input
    def  from_ndarray ( self ,  param ): 
        return  param . ctypes . data_as ( ctypes . POINTER ( ctypes . c_double ))

DoubleArray  =  DoubleArrayType () 
#--------------------------------------------------------------------------
class POINTDA(Structure):
    _fields_ = [('x', c_double*ms)]
cfuns.get_pointda.argtypes = (
    DoubleArray, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double
)
cfuns.get_pointda.restype = c_void_p
######################################################################################################
def tbif(t,y):
    #equation form d/dt(X+iY) = {A(X+iY)}/i
    z=t
    Q=y
    Qr = np.ndarray.tolist(Q.real)
    Qi = np.ndarray.tolist(Q.imag)
    #out double array real
    #-----------------------------------------------------------------------------
    outdar = POINTDA.from_address(
        cfuns.get_pointda(
            Qr,c_double(z), cd_ar, cd_omg, cd_omga, cd_Omg, cd_vxp, cd_vxn, cd_vz
        )
    )
    #out double array real numpy
    outdar_np0 = np.ctypeslib.as_array(outdar.x)
    outdar_np = outdar_np0 + np.zeros(ms)
    cfuns.free_pointda(byref(outdar))
    del outdar
    #----------------------------------------------------------------
    outdai = POINTDA.from_address(
        cfuns.get_pointda(
            Qi,c_double(z), cd_ar, cd_omg, cd_omga, cd_Omg, cd_vxp, cd_vxn, cd_vz
        )
    )
    outdai_np0 = np.ctypeslib.as_array(outdai.x)
    outdai_np = outdai_np0 + np.zeros(ms)
    cfuns.free_pointda(byref(outdai))
    del outdai
    #-----------------------------------------------------------------------------
    st = time.time()
    #print(outdar_np)
    #print(outdai_np)
    #print("outdar_np",outdar_np)
    #print("outdai_np",outdai_np)
    dqdz = outdar_np * (0-1j)+outdai_np
    #print(time.time()-st, "A i + B")
    return dqdz#,outdar,outdai
print("begin")
if test == True:
    st = time.time()
    for i in range(1):
        #print("pre sent")
        A = (tbif(0,np.array([1.0+0j]*ms)))
        #print("back sent")
    for i in A:
        print(i.real)
    #print((A),"1st A")
    #to_free(A[1],A[2])
    #print((A),"2nd A")
    print(time.time()-st, "should be fast")
    print(ms)

    sys.exit()
path = "/hetghome/jordan/NeuOsc/Cfin_repro_2007_f1_right/"
if ConsMu == True:
    sys,exit()
    os.mkdir(path + cur_dat+"_tur_evo_Cons")
if ConsMu == False:
    if test_run == True:
        trash = 1
    else:
        os.mkdir(path + cur_dat+"_tur_evo_dk_"+str(dk)+"_acu_"+str(int(m.log10(acu))))
        vec_path = path + cur_dat+"_tur_evo_dk_"+str(dk)+"_acu_"+str(int(m.log10(acu)))+"/evo_vec"
        os.mkdir(vec_path)
    ###############################################################
    print("start to solve")
    Sol = RK45((tbif),z0,np.array([(1.0+0j)/(ms**(0.5))]*ms),t1,max_step=ts, rtol=acu, atol=acu * (10**(-3)), vectorized=False, first_step=None)
    print('sol input-ed')
    Z = []
    Qz = []
    Count_z = 0
    pre_z = Sol.t
    st = time.time()    
    for j in range(step_num):
        if Count_z%1000 == 0:
            print("Sol.t ",Sol.t)
            if test_run == False:
                np.save(vec_path + "/vec_"+str(Sol.t),Sol.y)
                print("save a result")
            print("Sol.t ",Sol.t)
        #print("step size ", Sol.t-pre_z)
        #print("Sol.t ",Sol.t)
        pre_z = Sol.t
        Sol.step()
        #print(Count_z)
        Count_z +=1
        if Sol.status == 'finished':
            break
        #print("loop time", time.time()-st)
    print("tol time ",time.time()-st)
print("step_num ",step_num)



