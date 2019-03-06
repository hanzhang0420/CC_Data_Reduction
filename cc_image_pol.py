from scipy import ndimage as snd
import numpy as np
from photutils import centroid_com, centroid_1dg, centroid_2dg 
from photutils import make_source_mask
from astropy.stats import sigma_clipped_stats
import math
from skimage import data 
from skimage.feature import register_translation
from skimage.feature.register_translation import _upsampled_dft
from scipy.ndimage import fourier_shift

from astropy.io import fits
import matplotlib

from pylab import *
from matplotlib.ticker import MultipleLocator
from matplotlib import rc, font_manager

class CcData:
    
    if __name__ == "__main__":
        self.fits_file=fits_file
        
        # check file exists
        if not path.exists():
            self._index = []
        return

    
    def __init__(self,fits_file)
        self.fits_file=fits_file
    
        # check file exists
        if not path.exists():
            self._index = []
        return

        '''Readimage.py by Han Zhang (2018)
        
        Readimage.py is a function to read CanariCam polarimetry fits.gz files.
        Input: path of the fits file
        Output: images list (correspond to images at four different waveplate angles 0, 22.5, 67.5, 45) and instrument position angle
        
        '''
        image=fits.open(self.fits_file)
        #image.info()   # gives extensions
        nsave=image[0].header['NSAVSETS']
        nnods=image[0].header['NNODS']
        nnodsets=image[0].header['NNODSETS']#nod A or B different subtraction order  four different waveplate angles
        ZD1=image[0].header['ZD1'] ; ZD2=image[0].header['ZD2']
        RMA1=image[0].header['ROTATOR1'] ; RMA2=image[0].header['ROTATOR2']
        theta_instr=np.mean([ZD2,ZD2])-np.mean([RMA1,RMA2])
        obj=image[0].header['OBJECT']; f_name=image[0].header['FILTER1'] ; PA=image[0].header['INSTRPA']
        date=image[0].header['DATE']
        print('=============================================')
        print('This is the data of', obj, 'at', date, 'with filter',f_name, 'Position Angle', PA)
        print('=============================================')
        image_p1=[] ; image_p2=[] ; image_p3=[] ; image_p4=[]  #initial defite four list; list can contain various data types
    
        for i in range(1,nsave*nnods*nnodsets):
            if image[i].header['NOD'] == 'A':
                fir=0 ; sec=1
            elif image[i].header['NOD'] == 'B':
                fir=1 ; sec=0
            
            wplate=image[i].header['WPLATE']
            if wplate==float(0.0):
                image_p1.append(image[i].data[fir,:,:]-image[i].data[sec,:,:])
            elif wplate==float(45.0):
                image_p2.append(image[i].data[fir,:,:]-image[i].data[sec,:,:])
            elif wplate==float(22.5):
                image_p3.append(image[i].data[fir,:,:]-image[i].data[sec,:,:])
            elif wplate==float(67.5):
                image_p4.append(image[i].data[fir,:,:]-image[i].data[sec,:,:])

        print('for each waveplate, the number of frames', len(image_p1))#image_p1 is a list. len() returns the items in the list
    
        image_r=[]
        image_r.append(np.mean(image_p1,axis=0))
        image_r.append(np.mean(image_p2,axis=0))
        image_r.append(np.mean(image_p3,axis=0))
        image_r.append(np.mean(image_p4,axis=0))
    
        self.image = image_r
        self.theta_instr = theta_instr


    def crop(self,sizex,sizey1,sizey2,plots=None):
        '''
            global make_source_mask
            crop the images with separate o and e ray sizex is the cropped half xcoordinate sizey1,sizey2 the boundary in y direction
        '''

        o_ray=[] ; e_ray=[]; noise=[]
        for img in self.image:
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

        print('O ray images')
        
        for o in o_ray:
            mask = make_source_mask(o, snr=2.5, npixels=8, dilate_size=8); maskr= np.invert(mask)
            xc_o, yc_o = centroid_com(o, mask=maskr)
            print('center of the image',xc_o,yc_o)#centroid to find the center of the image to do alignment
            xc.append(xc_o); yc.append(yc_o)
            yc_o=int(round(yc_o));xc_o=int(round(xc_o))
            crop_image_o.append(o[(yc_o-2*sizey1):(yc_o+2*sizey2+1),(xc_o-2*sizex):(xc_o+2*sizex+1)]) #usually x axis is fine, y axis
                #   needs two parameters just in case the object is too close to the edge of the mask
            print('size of the cropped image', o[(yc_o-2*sizey1):(yc_o+2*sizey2+1),(xc_o-2*sizex):(xc_o+2*sizex+1)].shape)
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
        
        print('E ray images')

        for e in e_ray:
            mask = make_source_mask(e, snr=2.5, npixels=8, dilate_size=8); maskr= np.invert(mask)
            xc_e, yc_e = centroid_com(e,mask=maskr)
            print('center of the image',xc_e,yc_e)
            xc.append(xc_e); yc.append(yc_e)
            yc_e=int(round(yc_e));xc_e=int(round(xc_e))
            crop_image_e.append(e[(yc_e-2*sizey1):(yc_e+2*sizey2+1),(xc_e-2*sizex):(xc_e+2*sizex+1)])
            print('size of the cropped image', e[(yc_e-2*sizey1):(yc_e+2*sizey2+1),(xc_e-2*sizex):(xc_e+2*sizex+1)].shape)
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

    def align(self, ref,offset_image,plots=None): # if you want to show the cross-coorlation image, #image_s=image_align(ref,offset_image,2)
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
 
        #check
        if plots is not None:
            fig = plt.figure()
            ax1 = plt.subplot(3, 1, 1)
            ax2 = plt.subplot(3, 1, 2)
            ax3 = plt.subplot(3, 1, 3)
            ax1.set_title('Before Shift')
            ax1.imshow(ref,origin='lower')
            ax1.set_axis_off()
            ax2.imshow(offset_image,origin='lower')
            ax2.set_axis_off()
            ax3.imshow(cc_image.real,origin='lower')
        
            plt.show()

            fig = plt.figure()
            ax1 = plt.subplot(3, 1, 1)
            ax2 = plt.subplot(3, 1, 2)
            ax3 = plt.subplot(3, 1, 3)
            #ax4 = plt.subplot(4, 1, 4)
    
            ax1.imshow(ref,origin='lower')
            ax1.set_axis_off()
            ax2.imshow(image_s.real,origin='lower')
            ax2.set_axis_off()
            ax3.imshow(cc_image_s.real,origin='lower')
            # ax4.imshow(ref-image_s.real,origin='lower')
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

    def rebin_center(array, err_array, n):
        #make sure one vector is at the center of the image#############################################
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
        arr_nan = np.empty((narray[0],narray[1]))
        arr_nan[:] = np.nan
        arra_re=rebin_array.repeat(3, axis=0).repeat(3, axis=1)
    
        for i in range(0,narray[0],3):
            for j in range(0,narray[1],3):
                arr_nan[i,j]=arra_re[i,j]
            
        return arr_nan, rebin_err_array
    
    def instrp_correct(U,Q,p_instr,theta_instr,I):
        U_c=U-p_instr*cos(2.0*theta_instr*np.pi/180.0)*I
        Q_c=Q-p_instr*sin(2.0*theta_instr*np.pi/180.0)*I
    
        return U_c,Q_c





    def rebin(a, fshape):
        if len(a.shape)>1:
            sh = fshape[0],a.shape[0]//fshape[0],fshape[1],a.shape[1]//fshape[1]
            return a.reshape(sh).mean(-1).mean(1)
        else:
            sh = fshape,a.shape[0]//fshape
        return a.reshape(sh).mean(1)

    def cc_spec(self,width,p_instr,theta_instr):
        # global make_source_mask
        #crop the images with separate o and e ray sizex is the cropped half xcoordinate sizey1,sizey2 the boundary in y direction
        o_ray=[] ; e_ray=[]; noise=[]; ny=240
        for img in image:
            img=rebin(img[:,20:280],[240,130])
            #print 'size of the image', img.shape
            cen_o=int(round(101.0/(240/ny)))
            cen_e=int(round(135/(240/ny)))
            o_ray.append(img[(cen_o-width):(cen_o+width),:])
            e_ray.append(img[(cen_e-width):(cen_e+width),:])
            mask = make_source_mask(img, snr=1.5, npixels=10, dilate_size=10)
            maskr= np.invert(mask)
            mean, median, std = sigma_clipped_stats(img, sigma=2.3, iters=8, mask=mask)
            noise.append(std)
    
        xc=[];yc=[];crop_image_o=[];crop_image_e=[]
        I_1=np.mean(o_ray[0]+e_ray[0],axis=0); Q_1=np.mean(o_ray[0]-e_ray[0],axis=0)
        I_2=np.mean(o_ray[1]+e_ray[1],axis=0); Q_2=np.mean(o_ray[1]-e_ray[1],axis=0)
        I_3=np.mean(o_ray[2]+e_ray[2],axis=0); U_1=np.mean(o_ray[2]-e_ray[2],axis=0)
        I_4=np.mean(o_ray[3]+e_ray[3],axis=0); U_2=np.mean(o_ray[3]-e_ray[3],axis=0)
    
        Q=(Q_1-Q_2)/2.0
        U=(U_1-U_2)/2.0
        I=(I_1+I_2+I_3+I_4)/4.0
        n=I_1.shape
        ###Instrument correction
        U_c=U-p_instr*np.cos(2.0*theta_instr*np.pi/180.0)*I
        Q_c=Q-p_instr*np.sin(2.0*theta_instr*np.pi/180.0)*I
    
        factor=((2.0*width+1.0)**0.5)
        err_Q=(noise[0]+noise[1])/(2.0*factor)+np.zeros((n[0]))
        #Q_1 noise =sqrt(2)*noise[0]  Q noise=Q_1 noise/sqrt(2)
        err_U=(noise[2]+noise[3])/(2.0*factor)+np.zeros((n[0]))
        #same
        err_I=np.mean(noise)/((2.0)**0.5*factor)+np.zeros((n[0]))
        return Q, U, I, err_Q, err_U, err_I


class PlotData:
    def __self__(self,I,u,q,err_I,err_q,err_u):
        self.I = I
        self.u = u
        self.q = q
        self.err_I = err_I
        self.err_q = err_q
        self.err_u = err_u

    def plot_pol_decom(self, q_a,q_e,u_a,u_e,wavet,um):
        fig_size=plt.rcParams["figure.figsize"]
    
        fig_size[0]=10.0
        fig_size[1]=10.0
    
        plt.rcParams["figure.figsize"]=fig_size
        p,dp,theta,dt=p_theta(self.q,self.u,self.err_q,self.err_u)
        p_model=((u_a+u_e)**2.0+(q_a+q_e)**2.0)**(0.5)
        p_a = np.sqrt(q_a**2 + u_a**2)
        p_e = np.sqrt(q_e**2 + u_e**2)
    
        theta_a = 0.5*np.arctan(u_a/q_a)*180./np.pi
        theta_e = 0.5*np.arctan(u_e/q_e)*180./np.pi
    
        plt.figure()
    
        ax=subplot(221)
        plt.plot(wavet,q_a,'r--',label='absorption')
        plt.plot(wavet,q_e,'b:',label='emission')
        plt.plot(wavet,q_a+q_e,'k-',label='P')
        plt.errorbar(wavet,self.q,yerr=self.err_q,fmt='k.')
    
        plt.xlim([8.0,13])
        plt.ylim([-um,um])
        plt.setp(ax.get_xticklabels(), visible=False)
        text(11.0,0.8*um,'$q (Q/I) (\%)$',fontsize=13)
        plt.setp(ax.spines.values(), linewidth=2.0)
    
    
        ax=subplot(222)
        plt.plot(wavet,u_a,'r--',label='absorption')
        plt.plot(wavet,u_e,'b:',label='emission')
        plt.plot(wavet,u_a+u_e,'k-',label='P')
        plt.errorbar(wavet,self.u,yerr=self.err_u,fmt='k.')
    
        plt.xlim([8.0,13])
        plt.ylim([-um,um])
        text(11.0,0.8*um,'$u (U/I) (\%)$',fontsize=13)
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.yaxis.set_ticks_position('right')
        ax.yaxis.set_label_position('right')
        plt.setp(ax.spines.values(), linewidth=2.0)
    
        ax=subplot(223)
        plt.plot(wavet,p_a,'r--',label='absorption')
        plt.plot(wavet,p_e,'b:',label='emission')
        plt.plot(wavet,p_model,'k-',label='P')
        plt.errorbar(wavet,p,yerr=dp,fmt='k.')
    
        plt.xlabel("Wavelength ($\mu$m)")
        plt.xlim([8.0,13.])
        plt.ylim([0,um])
        text(12.,0.8*um,'$p (\%)$',fontsize=13)
        plt.setp(ax.spines.values(), linewidth=2.0)
    
        ax=subplot(224)
        plt.plot(wavet,theta_a,'r--',label='Absorption')
        plt.plot(wavet,theta_e,'b:',label='Emission')
    
        plt.errorbar(wavet,theta,yerr=dt,fmt='k.')
        plt.plot(wavet,theta,'k-',label='Total')
        plt.xlabel("Wavelength ($\mu$m)")
        plt.xlim([8.0,13])
        plt.ylim([-100,100])
        text(12.,80,'$\Theta (^{\circ})$',fontsize=13)
        ax.yaxis.set_label_position('right')
        ax.yaxis.set_ticks_position('right')
        legend = plt.legend(loc='upper left',prop={'size': 9})
    
        subplots_adjust(wspace=0.05,hspace=0.05)
        plt.setp(ax.spines.values(), linewidth=2.0)
        plt.savefig('output/Plot_decom.png',bbox_inches='tight')
        plt.show()
        return p_a,p_e,theta_a,theta_e


