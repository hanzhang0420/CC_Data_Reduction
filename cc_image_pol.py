from scipy import ndimage as snd
import numpy as np
from photutils import centroid_com, centroid_1dg, centroid_2dg 
from photutils import make_source_mask
from astropy.stats import sigma_clipped_stats
import math
import matplotlib.pyplot as plt
from skimage import data 
from skimage.feature import register_translation
from skimage.feature.register_translation import _upsampled_dft
from scipy.ndimage import fourier_shift
import numpy as np
import matplotlib.pyplot as plt

def crop(image,sizex,sizey1,sizey2,plots=None): 
   # global make_source_mask
     #crop the images with separate o and e ray 
    o_ray=[] ; e_ray=[]; noise=[]
    for img in image:
        #print 'size of the image', img.shape
        o_ray.append(img[0:122,:])
        e_ray.append(img[110:239,:])
        mask = make_source_mask(img, snr=1.5, npixels=10, dilate_size=10)
        maskr= np.invert(mask)
        mean, median, std = sigma_clipped_stats(img, sigma=2.3, iters=8, mask=mask)
        noise.append(std)
        
    xc=[];yc=[];crop_image_o=[];crop_image_e=[]
    
   # show the result  #columns = 4 Calculate the centroid of a 2D array 
   # mask the region out of the object
    if plots is not None:
        fig=plt.figure()
        rows = 4; columns=1; counts=0; 
        
    print 'O ray images'
    for o in o_ray:
        mask = make_source_mask(o, snr=2.5, npixels=8, dilate_size=8); maskr= np.invert(mask)
        xc_o, yc_o = centroid_com(o, mask=maskr) 
        print 'center of the image',xc_o,yc_o#centroid to find the center of the image to do alignment 
        xc.append(xc_o); yc.append(yc_o)
        yc_o=int(round(yc_o));xc_o=int(round(xc_o))
        crop_image_o.append(o[(yc_o-2*sizey1):(yc_o+2*sizey2+1),(xc_o-2*sizex):(xc_o+2*sizex+1)]) #usually x axis is fine, y axis 
        #   needs two parameters just in case the object is too close to the edge of the mask 
        print 'size of the cropped image', o[(yc_o-2*sizey1):(yc_o+2*sizey2+1),(xc_o-2*sizex):(xc_o+2*sizex+1)].shape
        if plots is not None:
            counts += 1
            fig.add_subplot(4, 1, counts)
            fig.subplots_adjust(hspace=0.1)
            plt.imshow(o,origin='lower')
            plt.plot(xc_o,yc_o,'ro')
            
    if plots is not None:   
        plt.show()
        
    if plots is not None:
        fig=plt.figure()
        rows = 4; columns=1; counts=0; 
        
    print 'E ray images'
    for e in e_ray:
        mask = make_source_mask(e, snr=2.5, npixels=8, dilate_size=8); maskr= np.invert(mask)
        xc_e, yc_e = centroid_com(e,mask=maskr) 
        print 'center of the image',xc_e,yc_e
        xc.append(xc_e); yc.append(yc_e)
        yc_e=int(round(yc_e));xc_e=int(round(xc_e))
        crop_image_e.append(e[(yc_e-2*sizey1):(yc_e+2*sizey2+1),(xc_e-2*sizex):(xc_e+2*sizex+1)])
        print 'size of the cropped image', e[(yc_e-2*sizey1):(yc_e+2*sizey2+1),(xc_e-2*sizex):(xc_e+2*sizex+1)].shape
        if plots is not None:
            counts += 1
            fig.add_subplot(4, 1, counts)
            fig.subplots_adjust(hspace=0.1)
            plt.imshow(e,origin='lower')
            plt.plot(xc_e,yc_e,'ro')
  
    if plots is not None:   
        plt.show()      
  #  plt.plot([x for x in xc],[y for y in yc],'ro')
  #  plt.legend('centroid of the o and e ray images')
  #  plt.show()
    return crop_image_o,crop_image_e,noise


def align(ref,offset_image,plots=None): # if you want to show the cross-coorlation image, #image_s=image_align(ref,offset_image,2)
    shift, error, diffphase = register_translation(ref,offset_image, 100)
    print("Detected subpixel offset (y, x): {}".format(shift))
    #cross-correlation 
    image_product = np.fft.fft2(ref) * np.fft.fft2(offset_image).conj()
    cc_image = np.fft.fftshift(np.fft.ifft2(image_product))
    
    # shift the offset image according to reference 
    image_s = fourier_shift(np.fft.fftn(offset_image), shift)
    image_s = np.fft.ifftn(image_s)
    
    # show the new cross correlation after shift the image 
    image_product_s = np.fft.fft2(ref) * np.fft.fft2(image_s).conj()
    cc_image_s = np.fft.fftshift(np.fft.ifft2(image_product_s))
    
    if plots is not None:
        fig = plt.figure()
        ax1 = plt.subplot(3, 1, 1)
        ax2 = plt.subplot(3, 1, 2)
        ax3 = plt.subplot(3, 1, 3)
        ax1.set_title('Before Shift')
        ax1.imshow(ref,origin='lower')
        ax2.imshow(offset_image,origin='lower')
        ax3.imshow(cc_image.real,origin='lower')
        ax3.set_axis_off()
        plt.show()
        ##########################################
        fig = plt.figure()
        ax1 = plt.subplot(4, 1, 1)
        ax2 = plt.subplot(4, 1, 2)
        ax3 = plt.subplot(4, 1, 3)
        ax4 = plt.subplot(4, 1, 4)
    
        ax1.imshow(ref,origin='lower')
        ax2.imshow(image_s.real,origin='lower')
        ax3.imshow(cc_image_s.real,origin='lower')
        ax4.imshow(ref-image_s.real,origin='lower')
        #plt.set_title("Cross-correlation")
        ax1.set_title('After Shift')
        plt.show()
    return image_s.real
    
def stokes(image_o,image_e,noise):
    I_1=image_o[0]+image_e[0] ; Q_1=image_o[0]-image_e[0]
    I_2=image_o[1]+image_e[1] ; Q_2=image_o[1]-image_e[1]
    I_3=image_o[2]+image_e[2] ; U_1=image_o[2]-image_e[2]
    I_4=image_o[3]+image_e[3] ; U_2=image_o[3]-image_e[3]
    #decide the shape of input array 
    n=I_1.shape
    
    Q=(Q_1-Q_2)/2.0
    U=(U_1-U_2)/2.0
    I=(I_1+I_2+I_3+I_4)/4.0
    err_Q=(noise[0]+noise[1])/2.0+np.zeros((n[0],n[1])) #Q_1 noise =sqrt(2)*noise[0]  Q noise=Q_1 noise/sqrt(2)
    err_U=(noise[2]+noise[3])/2.0+np.zeros((n[0],n[1])) #same 
    err_I=np.mean(noise)/math.sqrt(2.0)+np.zeros((n[0],n[1]))
    return Q, U, I, err_Q, err_U, err_I

def p_theta(U, Q, err_U, err_Q):
    Ip=np.sqrt(U**2+Q**2)
    err_p=np.sqrt((U*err_U)**2+(Q*err_Q)**2)/Ip
    theta=0.5*np.arctan(U,Q)
    err_t=0.5*(err_p/Ip)
    return Ip,err_p,theta,err_t

def rebin(array, err_array, n):
    narray=array.shape ; c_loc=np.array([int(narray[0]/2.0), int(narray[1]/2.0)])
    m=np.array([int(narray[0]/n),int(narray[1]/n)])
    bin=(m-1)/2 ; index_start=(c_loc-1)-bin*n; index_end=(c_loc+1)+bin*n
    new_array=array[index_start[0]:index_end[0]+1,index_start[1]:index_end[1]+1]
    rebin_array=np.reshape(new_array,\
                   ((2*bin[0]+1),n,(2*bin[1]+1),n))

    rebin_array = np.mean(rebin_array,axis=3)
    rebin_array = np.mean(rebin_array,axis=1)
    rebin_err_array = err_array/n
    
    narray=new_array.shape
    arr_nan = np.empty((narray[0],narray[1],))
    arr_nan[:] = np.nan
    arra_re=rebin_array.repeat(3, axis=0).repeat(3, axis=1)
    
    for i in range(0,narray[0],3):
        for j in range(0,narray[1],3):
            arr_nan[i,j]=arra_re[i,j]
            
    return arr_nan, rebin_err_array
    
    

