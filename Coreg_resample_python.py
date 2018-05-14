#!/usr/bin/env python

import os
from argparse import ArgumentParser
from CSI2DICOM_func import CSI_2_DCM
#from showresampling_func import showresampling
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
CSI_2_DCM(flat_list[1],flat_list[6],CSIdicom)
print('DONE WITH CREATING CSI file')



####################### Convert minc files ############################
PET_file = "PET"
cmd = "dcm2mnc {0}  -dname '' -fname {1} . -usecoordinates -clobber".format(flat_list[0],PET_file)
os.system(cmd)
print('DONE converting PET to minc')

MRI_at_CSI_2D = "MRI_at_CSI_2D"
cmd = "dcm2mnc {0}  -dname '' -fname {1} . -usecoordinates -clobber".format(flat_list[1],MRI_at_CSI_2D)
os.system(cmd)
print('DONE converting MRI_at_CSI_2D to minc')

MRI_at_CSI_3D = "MRI_at_CSI_3D"
cmd = "dcm2mnc {0}  -dname '' -fname {1} . -usecoordinates -clobber".format(flat_list[2],MRI_at_CSI_3D)
os.system(cmd)
print('DONE converting MRI_at_CSI_3D to minc')

MRI_at_PET_file = "MRI_at_PET"
cmd = "dcm2mnc {0}  -dname '' -fname {1} . -usecoordinates -clobber".format(flat_list[3],MRI_at_PET_file)
os.system(cmd)
print('DONE converting MRI_at_PET to minc')

CSI_file = "CSI"
cmd = "dcm2mnc {0}  -dname '' -fname {1} . -usecoordinates -clobber".format(CSIdicom,CSI_file)
os.system(cmd)
print('DONE converting pseudo CSI dicom to minc')


##################################### Register##############################

matrix_trans='matrix_trans'
cmd= "minctracc {0} {1} {2} -{3} -identity -est_center -clobber".format(MRI_at_PET_file+'.mnc',MRI_at_CSI_3D+'.mnc',matrix_trans,flat_list[4])
os.system(cmd)
print('DONE registration')
# Resample PET to CSI
#time.sleep(99999)
PET_coreg_MRI_CSI_2D="PET_coreg_MRI_CSI_2D.mnc"
PET_coreg_MRI_CSI_3D="PET_coreg_MRI_CSI_3D.mnc"
MRI_PET_coreg_MRI_CSI_3D="MRI_PET_coreg_MRI_CSI_3D.mnc"
MRI_PET_coreg_MRI_CSI_2D="MRI_PET_coreg_MRI_CSI_2D.mnc"

##################################### Resampling ##############################
# Resample PET to CSI
cmd="mincresample {0} -transformation {1} -like {2} {3} -{4} -clobber".format(PET_file+'.mnc',matrix_trans+'.xfm',MRI_at_CSI_3D+'.mnc',PET_coreg_MRI_CSI_3D,flat_list[5])
os.system(cmd)
print('DONE resampling PET to CSI')

# Resample MRI images (for check)
cmd="mincresample {0} -transformation {1} -like {2} {3} -{4} -clobber".format(MRI_at_PET_file+'.mnc',matrix_trans+'.xfm',MRI_at_CSI_3D+'.mnc',MRI_PET_coreg_MRI_CSI_3D,flat_list[5])
os.system(cmd)
print('DONE resampling MRI 3D to CSI (for check)')

# Resample PET to slice
cmd="mincresample {0} -like {1} {2} -{3} -clobber".format(PET_coreg_MRI_CSI_3D,MRI_at_CSI_2D+".mnc",PET_coreg_MRI_CSI_2D,flat_list[5])
os.system(cmd)
print('DONE resampling PET to 2D MRI at CSI')

# Resample MRI to slice
cmd="mincresample {0} -like {1} {2} -{3} -clobber".format(MRI_PET_coreg_MRI_CSI_3D,MRI_at_CSI_2D+".mnc",MRI_PET_coreg_MRI_CSI_2D,flat_list[5])
os.system(cmd)
print('DONE resampling MRI to 2D MRI at CSI (for check)')


####################################### Resample PET to CSI space  ##############


# Resample PET to CSI
PET_resample_CSI="PET_resampled_CSI.mnc"
cmd="mincresample {0} -like {1} {2} -{3} -clobber".format(PET_coreg_MRI_CSI_2D,CSI_file+".mnc",PET_resample_CSI,flat_list[5])
os.system(cmd)

print('DONE Resampling')

# convert to raw
PET_resample_CSI_raw="PET_resampled_CSI_raw.raw"
cmd="minctoraw {0} -nonormalize -double > {1}".format(PET_resample_CSI,PET_resample_CSI_raw)
os.system(cmd)

#print('DONE converting to raw')
#showresampling(PET_resample_CSI_raw)

#print('DONE show Resampling')

pet_dicom="pseudo_PET_DICOM_res_.dcm"
generate_dic(PET_resample_CSI_raw,CSIdicom,pet_dicom)

print('done generate pet_dicom')

