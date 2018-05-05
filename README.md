# CC_Data_Reduction

Codes CanariCam (CC) polarimetry data reduction 

Requirements:

numpy, scipy, photutils, skimage, matplotlib

Functions

readimage.py: read in the FITS file generate by CC. It outputs four 2-dimensionary arry, corresponding to four wave-plate posirion angles.

cc_image_pol.py: contains five functions, implementing image alignment, image crop, Stokes parameters computation, polarization computation and the resampling of the results. 

