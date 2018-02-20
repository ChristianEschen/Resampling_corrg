#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 22:24:17 2017

@author: christian
"""

import numpy as np
import dicom
#import nibabel as nib
from dicom.UID import UID, generate_uid, pydicom_root_UID, InvalidUID
from argparse import ArgumentParser
import os.path

def generate_dic(PETraw_path,dicom_path,pet_dicom):

#parser = ArgumentParser(description="ikjMatrix multiplication")
#parser.add_argument("-ResampledPETrawpath", dest="PETraw_path", required=True,
#    help="Resampled raw PET image", metavar="FILE")
#parser.add_argument("-CSIdicom_path", dest="dicom_path", required=True,
#    help="Pseudo dicom csi image", metavar="FILE")
#parser.add_argument("-PETdicom", dest="PETdcm", required=True,
#    help="pseudo valid dicom PET image", metavar="FILE")

#args = parser.parse_args()



    output=pet_dicom

    new_img = np.fromfile(PETraw_path, dtype='float64', sep="")
    new_img = new_img.reshape([int(np.sqrt(len(new_img))), int(np.sqrt(len(new_img)))])

    ds = dicom.read_file(csi_dicom_path)
    ds.pixel_array=new_img
    ds.PixelData=new_img.tostring()
    UIDnr=ds[0x0020,0x00d].value#[0:30]
    #
    temp=ds[0x028, 0x030].value
    temp2=[temp[1], temp[0]]
    ds[0x028, 0x030].value=temp2

    n=8
    groups = UIDnr.split('.')
    UID_gen_prefix='.'.join(groups[:n])+'.'#, '.'.join(groups[n:])


    uid_gen_Study=generate_uid(prefix=UID_gen_prefix, truncate=True)

    print('Generated UID studt:',uid_gen_Study)


    n2=9


    UIDnr2=ds[0x0020,0x00e].value#[0:30

    groups = UIDnr2.split('.')
    UID_gen_prefix='.'.join(groups[:n2])+'.'#, '.'.join(groups[n:])


    uid_gen_series=generate_uid(prefix=UID_gen_prefix, truncate=True)
    print('Generated UID series:',uid_gen_series)
    ds[0x0020,0x00e].value=uid_gen_series
    ds[0x008,0x060].value='PT'
    ds[0x008,0x103e].value='PET resampled in CSI space'

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
        
    write_dicom(ds, new_img, output)


