"""
@author: Mark
"""

#multiqtau in python.
import numpy as np
from tqdm import tqdm_notebook as tqdm
def insertimg(FD,n,plist,bind,nobind,nobuf,nolev,norm=1):
    """insertimg(FD,n,plist,bind,nobind,nobuf,nolev,norm=1)
       read and insert image, n, into accumulator buffers,
       plist is list of pixels, bind their bins and
       nobind is the number in each bind.
       then increments appropriate correlators for new image.
       Use function rdframe(n) to interface to the image reading routines.
       norm is 1 or SG used to normalize (from smoothed image).
    """
    global buf,cur,cts
    cur[0]=(1+cur[0])%nobuf
#     img=eliminate_module_beamstop(FD.rdframe(n))[1] #eliminate_module_beamsstop added by Peco
    img=FD.rdframe(n)
    buf[0,cur[0],:]= img.ravel()[plist]/norm
    process(0,cur[0],bind,nobind,nobuf)
    processing=1
    lev=1
    while(processing):
        if(cts[lev]):
            prev= (cur[lev-1]+nobuf-1)%nobuf
            cur[lev]= (cur[lev]+1)%nobuf
            buf[lev,cur[lev],:]=(buf[lev-1,prev,:]+buf[lev-1,cur[lev-1],:])/2.0
            cts[lev]=0
            process(lev,cur[lev],bind,nobind,nobuf)
            lev += 1
            processing= 1 if (lev<nolev) else 0
        else:
            cts[lev]=1
            processing=0

def process(lev,bufno,bind,nobind,nobuf):
    """ process(lev,bufno,bind,nobind,nobuf)
        The buffer bufno at level lev has just been added so update
        the correlation and intensity averages it affects. bind are
        bin indexes and nobind is histogram(bind)
    """
    global num,G,IAP,IAF
    num[lev]+=1
    imin=0 if lev==0 else nobuf//2

    for i in range(imin, min(num[lev],nobuf)):
        ptr=lev*nobuf//2+i
        delayno=(bufno-(i)+nobuf)%nobuf
        IP=buf[lev,delayno, :]
        IF=buf[lev,bufno,:]
        G[ptr,:]+= (np.bincount(bind,IF*IP)/nobind-  G[ptr,:])/(num[lev]-i)
        IAP[ptr,:] += (np.bincount(bind,IP)/nobind-IAP[ptr,:])/(num[lev]-i)
        IAF[ptr,:] += (np.bincount(bind,IF)/nobind-IAF[ptr,:])/(num[lev]-i)
    #print('proc: ', G[lev*nobuf//2+imin,:])

def autocor(FD,begframe,noframes,plist,bind,nobuf=8,nolev=8,skip=1,norm=1):
    """autocor(FD,begframe,noframes,plist,bind,nobuf=8,nolev=8,skip=1,norm=1)
        Function to drive the acquisition of data and to autocorrelate
        the data using multiple tau method. Uses variables in .des.i for info.
        begframe is frame number of first pattern to use.
        noframes is number of frames to process.
        plist is list of pixels to analyze (pixellist).
        bind is binning index for pixels.
        KEYWORD: norm=1 or SG  to normalize by Smooth image SG
                 nobuf=8 number of images at each level
                 nolev=8 number of levels (8 goes to 896 delays)
                 skip=1 analyze skip'th frame
    """
    global buf,cts,cur,G,IAP,IAF,num
    nobs=1+int(max(bind))
    buf=np.zeros((nolev,nobuf,plist.size),dtype=float)
    cts=np.zeros(nolev,dtype=int)
    cur=(nobuf-1)*np.ones(nolev,dtype=int)
    G=np.zeros([(nolev+1)*nobuf//2,nobs])
    IAP=np.zeros([(nolev+1)*nobuf//2,nobs])
    IAF=np.zeros([(nolev+1)*nobuf//2,nobs])
    num=np.zeros(nolev,dtype=int)
    nobind=np.bincount(bind)
    w=np.flatnonzero(nobind==0)
    if(w.size > 0):
        print("{:d} bins have no pixels in them.".format(w.size))
        nobind[w]=1 #prevent divide by 0s
    for n in tqdm(np.arange(0,noframes,skip)):
        insertimg(FD,begframe+n,plist,bind,nobind,nobuf,nolev,norm=norm)
    gmax=np.flatnonzero(IAP[:,0]==0)[0]
    if(gmax.size==0):gmax=IAP[:,0].size
    else: gmax=gmax.min()
    a=np.arange(nobuf//2)
    #generate times for G.
    tt=np.append(a,np.multiply.outer(\
            2**np.arange(1+(gmax-nobuf//2)//(nobuf//2)),a+nobuf//2).ravel())
    return(tt[:gmax],G[:gmax,:]/(IAP[:gmax,:]*IAF[:gmax,:]))
