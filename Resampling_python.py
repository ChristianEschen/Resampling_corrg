#!/usr/bin/env python

import os
from argparse import ArgumentParser
from CSI2DICOM_func import CSI_2_DCM
from showresampling_func import showresampling
from generate_dicom_func import generate_dic

# Inputs
parser = ArgumentParser(description="ikjMatrix multiplication")
parser.add_argument("-var_file", dest="var_PATH", required=True,
    help="var input file path", metavar="FILE")


args = parser.parse_args()
# Iterate in list files
path=args.var_PATH
content_list=list()
with open(path) as f:
    for line in f:
    	currentline = line.split(",")
    	content = [x.strip() for x in currentline]
    	content_list.append(content)
flat_list = [item for sublist in content_list for item in sublist]
flat_list[:] = [item for item in flat_list if item != '']

#######################
##################################### create valid CSI.dcm  header ###########################
CSIdicom="CSI_pseudo_dicom.dcm"
CSI_2_DCM(flat_list[1],flat_list[2],CSIdicom)

print('DONE WITH CREATING CSI file')
########################################### Convert minc files ###########################################################
PET_file = "PET"
cmd = "dcm2mnc {0}  -dname '' -fname {1} . -usecoordinates -clobber".format(flat_list[0],PET_file)
os.system(cmd)
print('DONE converting PET to minc')

# Convert CSI header to MINC
CSI_file = "CSI"

cmd = "dcm2mnc {0}  -dname '' -fname {1} . -usecoordinates -clobber".format(CSIdicom,CSI_file)
os.system(cmd)

# Convert CSI header to MINC
MRI_file = "MRI"

cmd = "dcm2mnc {0}  -dname '' -fname {1} . -usecoordinates -clobber".format(flat_list[1],MRI_file)
os.system(cmd)

print('DONE converting CSI to minc')
############################################ Resampling ###################################################
# Resample PET to CSI
PET_resample_CSI="PET_resampled_CSI.mnc"
cmd="mincresample {0} -like {1} {2} -{3} -clobber".format(PET_file+".mnc",CSI_file+".mnc",PET_resample_CSI,flat_list[3])
os.system(cmd)

print('DONE Resampling')
PET_resample_CSI_raw="PET_resampled_CSI_raw.raw"
cmd="minctoraw {0} -nonormalize -double > {1}".format(PET_resample_CSI,PET_resample_CSI_raw)
os.system(cmd)

print('DONE converting to raw')
showresampling(PET_resample_CSI_raw)

print('done show Resampling')

pet_dicom="pseudo_PET_DICOM_res.dcm"
generate_dic(PET_resample_CSI_raw,CSIdicom,pet_dicom)

print('done generate pet_dicom')