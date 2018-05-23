# rebin an array (2d or 1d), works as the IDL function rebin
import numpy
def rebin(a, fshape):
    if len(a.shape)>1:
        sh = fshape[0],a.shape[0]//fshape[0],fshape[1],a.shape[1]//fshape[1]
        return a.reshape(sh).mean(-1).mean(1)
    else:
        sh = fshape,a.shape[0]//fshape
        return a.reshape(sh).mean(1)
