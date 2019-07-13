"""
@author: Mark
"""
"""
New analysis code in python.
Master Module 
"""
import struct
import numpy as np
import time
import os


"""    Description of IMM files.

     This is code that Mark wrote to open the multifile format
    in compressed mode, translated to python.
    This seems to work for DALSA and FCCD in compressed mode
    should be included in the respective detector.i files
    Currently, this refers to the compression mode being '6'
    Each file is image descriptor files chunked together as follows:

    |--------------IMG N begin--------------|
    |        Header (1024 bytes)            |
    |---------------------------------------|
    |       Pixel positions (dlen*4 bytes   |
    |      (0 based indexing in file)       |
    |---------------------------------------|
    |    Pixel data(dlen*bytes bytes)       |
    |    (bytes is found in header          |
    |    at position 116)                   |
    |--------------IMG N end----------------|
    |--------------IMG N+1 begin------------|
    |----------------etc.....---------------|

    Here is the file layout as follows:
    0: mode (4 bytes)
    4: compression (4 bytes)
    8: date (32 bytes)
    40: prefix (16 bytes)
    56: number (4 bytes)
    60: suffix (16 bytes)
    76: monitor (4 bytes)
    80: shutter (4 bytes)
    84: row_beg (4 bytes)
    88: row_end (4 bytes)
    92: col_beg (4 bytes)
    96: col_end (4 bytes)
    100: row_bin (4 bytes)
    104: col_bin (4 bytes)
    108: rows (4 bytes)
    112: cols (4 bytes)
    116: bytes (4 bytes)
    120: kinetics (4 bytes)
    124: kinwinsize (4 bytes)
    128: elapsed (8 bytes)
    136: preset (8 bytes) in seconds
    144: topup (4 bytes)
    148: inject (4 bytes)
    152: dlen (4 bytes)
    156: roi_number (4 bytes)
    160: buffer_number (4 bytes)
    164: systick (4 bytes) in clock cycles
    608: imageserver (4 bytes)
    612: CPUspeed (4 bytes)
    616: immversion (4 bytes)
    620: corecotick (4 bytes) in microseconds
    624: cameratype (4 bytes)
    628: threshold (4 bytes)

"""

class Immfile:
    '''The class representing the uncompressed multifiles.
        The file is in 1 based numbering scheme (first record is 1)
        This assumes frame n is n frames into file for seeking.

        Member objects include:
        compression (4 bytes)
        108: rows (4 bytes)
        112: cols (4 bytes)
        116: byts (4 bytes, int)
        128: elapsed (8 bytes, double)
        136: preset (8 bytes, double)
        152: dlen (4 bytes)
        164: systick (8 bytes, double)
        620: corecotick (8 bytes, double)
        filepos (8 bytes)
        t0 (8 bytes, double)
        elapsed0 (8 bytes, double)
        recno (4 bytes, int)
        beg (4 bytes, int)
        end (4 bytes, int)
        filename (string)
    '''
    def __init__(self,filename,beg,end):
        '''Multifile initialization. Open the file.
            Here I use the read routine which returns byte objects
            (everything is an object in python). I use struct.unpack
            to convert the byte object to other data type (int object
            etc)
            NOTE: At each record n, the file cursor points to record n+1,
            so ready to read next image.
        '''
        self.FID = open(filename,"rb")
        self.beg = beg
        self.end = end
        self.number = end-beg+1
        self.filename = filename
        #br: bytes read
        br = self.FID.read(1024)
        self.recno = 1
        self.DK=0.0
        self.THRESHOLD=0.0
        # some initialization stuff
        (self.row_beg,self.row_end,self.col_beg,self.col_end,self.row_bin,self.col_bin,self.rows,self.cols,self.byts) = struct.unpack("@IIIIIIIII",br[84:120])
        if (self.byts==2):
            self.valtype = np.uint16
        elif (self.byts == 4):
            self.valtype = np.uint32
        #now convert pieces of these bytes to our data
        (self.elapsed,self.preset) = struct.unpack("@dd",br[128:144])
        #(self.dlen,) = struct.unpack("@I",br[152:156])
        (self.systick,) = struct.unpack("@d",br[164:172])
        # now read first image
        #print("Opened file. Bytes per data is {0}".format(self.byts))

    def _readImage(self):
        img=np.fromfile(self.FID,count=self.rows*self.cols,dtype=self.valtype).reshape(self.rows,self.cols)
        return(img)

    def seekimg(self,n):
        '''Position file to read the nth image.
        '''
        if((n<self.beg)|(n>self.end)):raise IOError('no record %d'%n)
        #pos=(self.rows*self.cols*self.byts+1024)*(n-self.beg)
        pos=(self.rows*self.cols*self.byts+1024)*(n-1)
        self.FID.seek(pos,0)
    def rdimg(self,n):
        '''raw read image n. Use rdframe(n) for processed image.
        '''
        self.seekimg(n)
        br = self.FID.read(1024)
        (self.elapsed,self.preset) = struct.unpack("@dd",br[128:144])
        (self.systick,) = struct.unpack("@d",br[164:172])
        return(self._readImage())
    def rdframe(self,n,tdk=True,tflag=True):
        ''' read processed images. Used by all analysis routines.
        '''
        img=self.rdimg(n).astype(dtype=np.float64) #newer numpy
        #img=self.rdimg(n).astype('float64') # older numpy
        if(tdk): img -= self.DK
        if(tflag): img *= (img>self.THRESHOLD)
        return(img)

#Multifiles
class Multifile:
    '''The class representing the compressed multifiles.
        The recno is in 1 based numbering scheme (first record is 1)
        This is efficient for reading in increasing order.
        Note: reading same image twice in a row is like reading an earlier
        numbered image and means the program starts from the beginning again.

        Member objects include:
        compression (4 bytes)
        108: rows (4 bytes)
        112: cols (4 bytes)
        116: byts (4 bytes, int)
        128: elapsed (8 bytes, double)
        136: preset (8 bytes, double)
        152: dlen (4 bytes)
        164: systick (8 bytes, double)
        620: corecotick (8 bytes, double)
        filepos (8 bytes)
        t0 (8 bytes, double)
        elapsed0 (8 bytes, double)
        recno (4 bytes, int)
        beg (4 bytes, int)
        end (4 bytes, int)
        filename (string)
    '''
    def __init__(self,filename,beg,end):
        '''Multifile initialization. Open the file.
            Here I use the read routine which returns byte objects
            (everything is an object in python). I use struct.unpack
            to convert the byte object to other data type (int object
            etc)
            NOTE: At each record n, the file cursor points to record n+1
        '''
        self.FID = open(filename,"rb")
#        self.FID.seek(0,os.SEEK_SET)
        self.beg = beg
        self.end = end
        self.number = end-beg+1
        self.filename = filename
        #br: bytes read
        br = self.FID.read(1024)
        self.imgread=0
        self.recno = 1
        # some initialization stuff
        (self.rows,self.cols,self.byts) = struct.unpack("@III",br[108:120])
        if (self.byts==2):
            self.valtype = np.uint16
        elif (self.byts == 4):
            self.valtype = np.uint32
        #now convert pieces of these bytes to our data
        (self.elapsed,self.preset) = struct.unpack("@dd",br[128:144])
        (self.dlen,) = struct.unpack("@I",br[152:156])
        #print 'i',self.dlen
        (self.systick,) = struct.unpack("@d",br[164:172])
        # now read first image
        #print "Opened file. Bytes per data is {0}".format(self.byts)

    def _readHeader(self):
        br = self.FID.read(1024)
        (self.elapsed,self.preset) = struct.unpack("@dd",br[128:144])
        (self.dlen,) = struct.unpack("@I",br[152:156])
        #print self.dlen
        (self.systick,) = struct.unpack("@d",br[164:172])

    def _readImageRaw(self):
        p = np.fromfile(self.FID,count=self.dlen,dtype=np.uint32)
        self.imgread=1
        return(p,np.fromfile(self.FID,count=self.dlen,dtype=self.valtype))

    def _readImage(self):
        (p,v)=self._readImageRaw()
        img = np.zeros((self.rows*self.cols))
        np.put(img,p, v)
        img.shape = (self.rows,self.cols)
        return(img)

    def seekimg(self,n=None):
        '''Position file pointer to read the nth image.
        '''
        # the logic involving finding the cursor position
        if (n==None):
            n = self.recno
        if (n < self.beg or n > self.end):
            print("Error, record out of range {:4d}".format(n))
            return -1

        if ((n == self.recno)  and (self.imgread==0)):
            pass # do nothing
        else:
            if (n <= self.recno): #ensure cursor less than search pos
                self.FID.seek(0,os.SEEK_SET) #go to beginning
                self.recno = 0
                self.imgread=1
            #have to iterate on seeking since dlen varies
            #remember for rec recno, cursor is always at recno+1
            if(self.imgread==0): #move to next header if need to
                self.FID.seek(self.dlen*(4+self.byts),os.SEEK_CUR)
            for i in range(self.recno+1,n):
                #the less seeks performed the faster
                br = self.FID.read(1024)
                (self.dlen,) = struct.unpack("@I",br[152:156])
                self.FID.seek(self.dlen*(4+self.byts),os.SEEK_CUR)

            # we are now at recno in file, read the header and data
            #self._clearImage()
            self._readHeader()
            self.imgread=0
            self.recno = n
    def rdframe(self,n):
        self.seekimg(n)
        return(self._readImage())

    def rdrawframe(self,n):
         self.seekimg(n)
         return(self._readImageRaw())

def rdmask(file,datadesc,pos=192,type=np.int64):
    ''' mask=rdmask(filename,pos=192,type=np.int64)
        read a mask from a binary file. This reads a set
        of data of type=type, starting at pos=pos.
        Note defaults for a yorick created mask.
        '''
    m=rdbin(file,datadesc,pos,type) 
    m=np.int16(m==0) #invert mask 0 and 1's
    return(m)

def rdbin(file,datadesc,pos=192,type=np.int64):
    ''' mask=rdbin(filename,pos=192,type=np.int64)
        read an image from a binary file. This reads an image
        of type=type, starting at pos=pos (size determined 
	from datadesc structure.
        Note default pos for a yorick binary file.
        '''
    f=open(file,'rb')
    f.seek(pos)
    m=np.fromfile(f,count=datadesc['rows']*datadesc['cols'],dtype=type)
    m.shape=(datadesc['rows'],datadesc['cols'])
    return(m)

#ldimgs
def ldimgs(FD,begframe=None,noframes=None):

    if(begframe==None):begframe=FD.beg
    if(noframes==None):noframes=FD.end-FD.beg+1
    ta = time.time()
    imgs=np.zeros([noframes, FD.row_end-FD.row_beg,FD.col_end-FD.col_beg])
    p=0
    for i in range(begframe,begframe+noframes):
        imgs[p,:,:]=FD.rdframe(i)
        p+=1
    tb = time.time()
    print("total time: {0}".format(tb-ta))
    return(imgs)
