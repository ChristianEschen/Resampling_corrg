#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 20:27:24 2017

@author: christian
"""
import numpy as np
import numpy.matlib
import re
import dicom
from operator import truediv
import copy
from argparse import ArgumentParser
import os.path
from dicom.UID import UID, generate_uid, pydicom_root_UID, InvalidUID
import time

def CSI_2_DCM(MRI_PATH,CSI_PATH,output):
    #args.MRI_PATH
    #args.CSI_PATH
    file=open(CSI_PATH, 'r')
    a=file.read()


    # Read image position patient

    IPP_loc=a.find('ImagePositionPatient')
    IPP=re.findall(r"[-+]?\d*\.\d+|\d+", a[IPP_loc+30:IPP_loc+300]) # 30 was nessesary
    IPP=map(float,IPP)
    IPP=IPP[0:3]
    print('readed image position patient', IPP)
    #time.sleep(9999999)
    

    # Read image orientation patient
    IOP_string='ImageOrientationPatient'
    IOP_loc=a.find(IOP_string)
    IOP=re.findall(r"[-+]?\d*\.\d+|\d+", a[IOP_loc:IOP_loc+300])
    IOP=map(float,IOP)
    IOP=IOP[0:3]
    print('readed image orientation patient',IOP)
    
    # Pixel spacing

    PS_loc=a.find('PixelSpacing')
    PS=re.findall(r"[-+]?\d*\.\d+|\d+", a[PS_loc:PS_loc+300])
    PS=map(float,PS)
    PS=PS[0:3]
    print('readed Pixel spacing', PS)


    # read matrix size

    matrix_size_1_loc =a.find('sSpecPara.lFinalMatrixSizePhase')
    matrix_dim=re.findall(r"[-+]?\d*\.\d+|\d+", a[matrix_size_1_loc:matrix_size_1_loc+300])
    matrix_dim=map(float,matrix_dim)
    matrix_dim=matrix_dim[0:2]
    print('The matrix dimensions of the CSI image is:',matrix_dim)
    if matrix_dim[0]!=matrix_dim[1]:
        print('warning matrix dimensions of CSI image does not agree')
        
    MheaderORG = dicom.read_file(MRI_PATH)
    CSIheaderORG = dicom.read_file(CSI_PATH)
    



    Mheader=copy.copy(MheaderORG)

    CSIheader=copy.copy(MheaderORG)
    # Find phase encoding direction

    phase_direction=MheaderORG.InplanePhaseEncodingDirection
    print('phase_direction',phase_direction)

    ############################
    if phase_direction== 'COL':
        FOV_loc2=a.find('sSpecPara.sVoI.dPhaseFOV')# 2
        FOV_loc1=a.find('sSpecPara.sVoI.dReadoutFOV') # 1
        print("Phase encoding direction is in column direction")
        print('FOV_loc2',FOV_loc2)
        print('FOV_loc1',FOV_loc1)
    else:
        ##################### changed to ::::
        FOV_loc1=a.find('sSpecPara.sVoI.dPhaseFOV') #1 
        FOV_loc2=a.find('sSpecPara.sVoI.dReadoutFOV') # 1
        #################
        print('FOV_loc2',FOV_loc2)
        print('FOV_loc1',FOV_loc1)
        print("Phase encoding direction is in row direction")
   
    #############################

    # Localize FOV
    FOV_loc3=a.find('sSpecPara.sVoI.dThickness')

    FOV1=re.findall(r"[-+]?\d*\.\d+|\d+", a[FOV_loc1:FOV_loc1+300])
    FOV1=map(float, FOV1)
    FOV1=FOV1[0]
    FOV2=re.findall(r"[-+]?\d*\.\d+|\d+", a[FOV_loc2:FOV_loc2+300])
    FOV2=map(float, FOV2)
    FOV2=FOV2[0]
    FOV3=re.findall(r"[-+]?\d*\.\d+|\d+", a[FOV_loc3:FOV_loc3+300])
    FOV3=map(float, FOV3)
    FOV3=FOV3[0]


     ############################ # check if patient is feet first prone
    if MheaderORG.PatientPosition== 'FFP':
        FOV=[FOV2,FOV1,FOV3]
    else:
        ##################### most likely case
        FOV=[FOV1,FOV2,FOV3]
   
    #############################

    #FOV=[FOV1,FOV2,FOV3]


    print('The field of view in CSI is:', FOV)



    # Add files from Mheader

    fields=len(Mheader)
    for iterator in range(fields):
        keynum=Mheader.keys()[iterator]
        mystring=Mheader[keynum]
        try:
            CSIheaderORG[keynum]
        except KeyError:  
            CSIheaderORG.add_new(mystring.tag, mystring.VR, mystring.value)


    

    # add specific pixel spacing and slice thickness 
    PixelSpacing=map(truediv,FOV[0:2],matrix_dim)
    PixelSpacing=map(float,PixelSpacing)
    CSIheaderORG.PixelSpacing=PixelSpacing
    CSIheaderORG.SliceThickness=int(FOV3)
    CSIheaderORG.SeriesDescription='CSI_dicom'

    ##### CALCULATE patient coordinates
    R=np.zeros((4,4))
    R[0:3,0] = MheaderORG.ImageOrientationPatient[0:3]
    R[0:3,1] = MheaderORG.ImageOrientationPatient[3:6]
    R[0:3,2] = np.cross(MheaderORG.ImageOrientationPatient[0:3],MheaderORG.ImageOrientationPatient[3:6])
    R[3,3] = 1
    #print('R',R,MheaderORG.ImageOrientationPatient)
    nCols = MheaderORG.Columns
    nRows= MheaderORG.Rows
    m=nRows*nCols
    i=np.matlib.repmat(np.array(range(0,nCols)),nRows,1)
    a=np.expand_dims(np.transpose(np.array(range(0,nRows))),1)
    #print('a',a,a.shape)
    j=np.matlib.repmat(a,1, nCols)

    M = R;
    M[0:3,0] = R[0:3,0] * MheaderORG.PixelSpacing[1];
    M[0:3,1] = R[0:3,1] * MheaderORG.PixelSpacing[0];
    M[0:3,3] = np.transpose(MheaderORG.ImagePositionPatient)

    i_new = np.reshape(np.transpose(i),(1,m))
    print('i_new',i_new.shape,i_new[0][509:514])
    #time.sleep(999999)
    j_new= np.reshape(np.transpose(j),(1,m))
    print('j_new',j_new.shape,j_new[0][0:10])
    zeross=np.zeros((1,m))
    oness=np.ones((1,m))
    imageCoordinates = np.concatenate((i_new,j_new,zeross,oness),axis=0)
    patientCoordinates = np.matmul(M,imageCoordinates)



    # Estimate center of mass
    center_of_mass = np.sum(patientCoordinates,1)/m
    center_of_mass = np.squeeze(center_of_mass[0:3])
    print('center_of_mass',center_of_mass)
    #print('([coords[0,(MheaderORG.Columns-1)/2+1]',coords[0,(MheaderORG.Columns-1)/2+1])
    
    
    ## Extract FOV i world coordiantes
    [Xq1,Yq1,Zq1]=np.meshgrid(np.linspace(center_of_mass[0]-np.divide(FOV[0],2),center_of_mass[0]+np.divide(FOV[0],2),matrix_dim[0]),
    np.linspace(center_of_mass[1]-np.divide(FOV[1],2),center_of_mass[1]+np.divide(FOV[1],2),matrix_dim[1]),np.linspace(center_of_mass[2]-np.divide(FOV[2],2),center_of_mass[2]+np.divide(FOV[2],2),matrix_dim[1]))

    #print('Xq1START AND END',Xq1[0][0][0],Xq1[-1][-1][-1])
    #print('Yq1START AND END',Yq1[0][0][0],Yq1[-1][-1][-1])
    #print('Yq1START AND END',Zq1[0][0][0],Zq1[-1][-1][-1])
    #print('IPP_MR',MheaderORG.ImagePositionPatient)
    #print('IPP',IPP)
    #print('IPP_cal_x',np.array(MheaderORG.PixelSpacing[0])/2+Xq1[0][0][0])
    #print('IPP_cal_y',np.array(MheaderORG.PixelSpacing[1])/2+Yq1[0][0][0])
    #print('IPP_cal_z',(Zq1[0][0][0]+Zq1[-1][-1][-1])/2)

    #print('matrix_dim[1])',int(matrix_dim[1]))
    i=np.matlib.repmat(np.array(range(0,int(matrix_dim[1]))),int(matrix_dim[0]),1)
    a=np.expand_dims(np.transpose(np.array(range(0,int(matrix_dim[0])))),1)
    j=np.matlib.repmat(a,1, matrix_dim[1])

    #print('MheaderORG.ImagePositionPatient',MheaderORG.ImagePositionPatient)
    #print('center_of_mass',list(center_of_mass))
    M = R;
    M[0:3,0] = R[0:3,0] * FOV[0]/matrix_dim[0];
    M[0:3,1] = R[0:3,1] * FOV[1]/matrix_dim[1];
    M[0:3,3] = np.transpose(list(center_of_mass))
    
    i_new = np.reshape(np.transpose(i),(1,int(matrix_dim[0]*matrix_dim[1])))
    j_new= np.reshape(np.transpose(j),(1,int(matrix_dim[0]*matrix_dim[1])))
    zeross=np.zeros((1,int(matrix_dim[0]*matrix_dim[1])))
    oness=np.ones((1,int(matrix_dim[0]*matrix_dim[1])))
    imageCoordinates= np.concatenate((i_new ,j_new ,oness,zeross),axis=0)
    patientCoordinates = np.matmul(M,imageCoordinates)
    #print('patientCoordinates',np.amax(patientCoordinates,1)+np.amax(patientCoordinates,1)/2)
    #time.sleep(9999)
 
 # OLD###################################################
    # Calculate center of mass

    #xc=np.divide(np.sum(Xq1[1,:,1]),matrix_dim[0])
    #yc=np.divide(np.sum(Yq1[:,1,1]),matrix_dim[1])
    #zc=np.divide(np.sum(Zq1[1,1,:]),matrix_dim[1])

    #xt=np.transpose(Xq1[1,:,1])-xc
    #yt=np.transpose(Yq1[:,1,1])-yc
    #zt=np.transpose(Zq1[1,1,:])-zc


    #rota=np.dot(R[0:3,0:3],np.row_stack((xt,yt,zt)))
    #rota2=np.array([np.matlib.repmat(xc,int(matrix_dim[0]),1),
    #                np.matlib.repmat(yc,int(matrix_dim[1]),1),
    #np.matlib.repmat(zc,int(matrix_dim[1]),1)])

    #rota2=np.squeeze(rota2, axis=(2,))#.shape
    #rota2=rota+rota2
    ##################################################
    IPP_x=np.array(MheaderORG.PixelSpacing[0])/2+Xq1[0][0][0]
    IPP_y = np.array(MheaderORG.PixelSpacing[1])/2+Yq1[0][0][0]
    IPP_z = (Zq1[0][0][0]+Zq1[-1][-1][-1])/2
    IPP_cal=np.array((IPP_x,IPP_y,IPP_z))
    result=np.allclose(IPP, IPP_cal)
    if result is None:
        print("########################################################################################")
        print("WARNING - No correspondence between calculated ImagePositioPatient and readed ImagePositionPaiten")
        print("########################################################################################")

    elif result:
        print("Correspondence calculated between ImagePositioPatient and readed ImagePositionPaitens")
        #time.sleep(9999)
    else:
        print("########################################################################################")
        print("WARNING - No correspondence between calculated ImagePositioPatient and readed ImagePositionPaiten")
        print("########################################################################################")

    print('IPP',IPP,IPP_cal)
    IPP_cal=IPP_cal.tolist()
    CSIheaderORG.ImagePositionPatient=IPP
    
    #time.sleep(9999999)  
    CSIheaderORG.Rows=int(matrix_dim[0])
    CSIheaderORG.Columns=int(matrix_dim[1])


    new_img=np.zeros((int(matrix_dim[0]),int(matrix_dim[1]))).astype('uint16')
    print('create CSI dicom image with zeros')


    CSIheaderORG[0x0028,0x0107].value = 4095
    CSIheaderORG[0x0028,0x0106].value = 4095
    CSIheaderORG[0x0028,0x0107].VR = 'US'
    CSIheaderORG[0x0028,0x0106].VR = 'US'
    CSIheaderORG[0x7fe0,0x0010].VR = 'OW'

    #print('Mheader',Mheader)
    #time.sleep(99999999)

    #  Generate UID

    n2=9


    UIDnr2=CSIheaderORG[0x0020,0x00e].value#[0:30

    groups = UIDnr2.split('.')
    UID_gen_prefix='.'.join(groups[:n2])+'.'#, '.'.join(groups[n:])


    uid_gen_series=generate_uid(prefix=UID_gen_prefix, truncate=True)
    print('Generated UID series:',uid_gen_series)
    CSIheaderORG[0x0020,0x00e].value=uid_gen_series
    CSIheaderORG[0x008,0x060].value='MR'
    CSIheaderORG[0x008,0x103e].value='CSI metabolite image'


    print('CSI header ORG',CSIheaderORG.ImagePositionPatient)
    print('IPP',IPP)
    
    def write_dicom(ds,pixel_array,filename):
        """
        INPUTS:
        pixel_array: 2D numpy ndarray.  If pixel_array is larger than 2D, errors.
        filename: string name for the output file.
        Image array are stored as uint16
        """

        ## These are the necessary imaging components of the FileDataset object.
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.PixelRepresentation = 0
        ds.HighBit = 15 #63
        ds.BitsStored = 16 #64
        ds.BitsAllocated = 16 #64
        #ds.SmallestImagePixelValue = '\\x00\\x00'
        #ds.LargestImagePixelValue = '\\xff\\xff'
        ds.Columns = pixel_array.shape[0]
        ds.Rows = pixel_array.shape[1]
        if pixel_array.dtype != np.uint16:
            pixel_array = pixel_array.astype(np.uint16)
    #    if pixel_array.dtype != np.float64:
    #        pixel_array = pixel_array.astype(np.float64)
        ds.PixelData = pixel_array.tostring()

        ds.save_as(filename)
        return
    write_dicom(CSIheaderORG, new_img, output)
    #CSIheader.save_as(args.outputFolder+args.output)

    print("CSI header saved as",output)

