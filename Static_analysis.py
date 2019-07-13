# -*- coding: utf-8 -*-
"""
Last modified date 02/27/2019

@author: Peco
"""
############################################################################################################
def get_intensity_for_pixel_box(yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe):
    Is=np.zeros((endframe-startframe,2))
    pos = 0
    for i in tqdm(range(startframe,endframe)):
        img = DD.FD.rdframe(i)
        Is[pos,0] = pos*datatakingtime #labelling time since it will be useful for fitting with log later .. use either i or pos here
        Is[pos,1] = (np.sum(img[yslicestart:ysliceend,xslicestart:xsliceend]))/((ysliceend-yslicestart)*(xsliceend-xslicestart))
        pos += 1
    return(Is)

def plot_intensity_for_pixel_box(Is):
    plt.figure(figsize=(7,7))
    plt.scatter(Is[:,0], Is[:,1], color = 'black', marker = '.', label='Data')
    plt.ylabel('Intensity')
    plt.xlabel('Time (per frame)')
    plt.legend(loc='best')
    plt.show()


def plot_intensity_adv_for_pixel_box(Is,yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe):
    """
    same as plot_intensity_for_pixel_box but can label things and save the output plot 
    plots the intensity that is given, plots the fitted function, save the image to the folder
    
    make sure to change the p0 guess when appropriate
    """
    plt.figure(figsize=(7,7))
    plt.scatter(Is[:,0], Is[:,1], color = 'black', marker = '.', label='Data')
    plt.ylabel('Intensity')
    plt.xlabel('Time (per frame)')
    plt.legend(loc='best')
    plt.suptitle(' \n yslicestart = {0}'.format(yslicestart) + 
                 ' \n ysliceend = {0}'.format(ysliceend) + 
                 ' \n xslicestart = {0}'.format(xslicestart) + 
                 ' \n xsliceend = {0}'.format(xsliceend) + 
                 ' \n startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe) , fontsize=8)
    
    plt.show()
    plt.savefig(' yslicestart = {0}'.format(yslicestart) + 
                 ' ysliceend = {0}'.format(ysliceend) + 
                 ' xslicestart = {0}'.format(xslicestart) + 
                 ' xsliceend = {0}'.format(xsliceend) + 
                 ' startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe))

def exponential_function(x, a, b, c):
    """
    used by functions for fitting R, b is R here. 
    """
    return a * np.exp(2*b*x) + c

def plot_intensity_andfit_for_pixel_box(Is,yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe,fitsign):
    """
    plots the intensity that is given, plots the fitted function, save the image to the folder
    
    make sure to change the p0 guess when appropriate
    """
    params, params_covariance = optimize.curve_fit(exponential_function, Is[:,0], Is[:,1], p0=[1.3, fitsign*0.04,0.1]) 
    plt.figure(figsize=(7,7))
    plt.scatter(Is[:,0], Is[:,1], color = 'black', marker = '.', label='Data')
    plt.plot(Is[:,0], exponential_function(Is[:,0], params[0], params[1],params[2]),label='Fitted function')
    plt.ylabel('Intensity')
    plt.xlabel('Time (per frame)')
    plt.legend(loc='best')
    plt.suptitle('parameters for a*exp(2*b*x)+c are {0}'.format(params) + ' \n yslicestart = {0}'.format(yslicestart) + 
                 ' \n ysliceend = {0}'.format(ysliceend) + 
                 ' \n xslicestart = {0}'.format(xslicestart) + 
                 ' \n xsliceend = {0}'.format(xsliceend) + 
                 ' \n startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe) , fontsize=8)
    
    plt.show()
    plt.savefig(' yslicestart = {0}'.format(yslicestart) + 
                 ' ysliceend = {0}'.format(ysliceend) + 
                 ' xslicestart = {0}'.format(xslicestart) + 
                 ' xsliceend = {0}'.format(xsliceend) + 
                 ' startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe))

############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxx    START    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################
def read_and_plot_intensity_pixel(yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe):
    Is = get_intensity_for_pixel_box(yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe)
    plot_intensity_for_pixel_box(Is)
    
def read_and_plot_fit_intensity_pixel(yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe):
    Is = get_intensity_for_pixel_box(yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe)
    plot_intensity_andfit_for_pixel_box(Is,yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe,fitsign)
    
def read_and_plot_fit_intensity_along_q_parallel(yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe,fitsign):
    loop = np.arange(xslicestart,xsliceend,delx)
    Rs = np.zeros((len(loop),2))
    pos = 0
    for i in tqdm(loop):
        Is = get_intensity_for_pixel_box(yslicestart,ysliceend,i,i+delx,startframe,endframe)
        plot_intensity_andfit_for_pixel_box(Is,yslicestart,ysliceend,i,i+delx,startframe,endframe,fitsign)

def read_and_plot_intensity_along_q_parallel(yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe):
    loop = np.arange(xslicestart,xsliceend,delx)
    Rs = np.zeros((len(loop),2))
    pos = 0
    for i in tqdm(loop):
        Is = get_intensity_for_pixel_box(yslicestart,ysliceend,i,i+delx,startframe,endframe)
        plot_intensity_adv_for_pixel_box(Is,yslicestart,ysliceend,i,i+delx,startframe,endframe)
    
############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxx     END     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################
def read_and_get_plot_R_along_q_parallel(yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe, fitsign):
    Rs = read_and_get_R_along_q_parallel(yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe, fitsign)
    plt.figure(figsize=(7,7))
    plt.scatter(Rs[:,0], Rs[:,1]*1E3, color = 'black', marker = '.', label='R')
    plt.ylabel(r'R($q_{y}$) ($10^{-3} s^{-1}$)')
    plt.xlabel(r'q ($nm^{-1}$)')
    plt.legend(loc='best')
    plt.suptitle('yslicestart = {0}'.format(yslicestart) + 
                 ' \n ysliceend = {0}'.format(ysliceend) + 
                 ' \n xslicestart = {0}'.format(xslicestart) + 
                 ' \n xsliceend = {0}'.format(xsliceend) + 
                 ' \n delx = {0}'.format(delx) + 
                 ' \n startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe) , fontsize=8)
    
    plt.show()
    plt.savefig(' yslicestart = {0}'.format(yslicestart) + 
                 ' ysliceend = {0}'.format(ysliceend) + 
                 ' xslicestart = {0}'.format(xslicestart) + 
                 ' xsliceend = {0}'.format(xsliceend) +
                '  delx = {0}'.format(delx) + 
                 ' startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe)) 

############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxx    START    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################
def read_and_get_R_along_q_parallel(yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe, fitsign):
    """
    make sure you run "read_and_plot_fit_intensity_pixel" first to make sure the fit is decent then run this
    
    also make sure to change p0 guess whenver appropriate so that the fitting converges
    """
    loop = np.arange(xslicestart,xsliceend,delx)
    Rs = np.zeros((len(loop),2))
    pos = 0
    params_init = [1.3, 0.04,0.1]
    for i in tqdm(loop):
        Is = get_intensity_for_pixel_box(yslicestart,ysliceend,i,i+delx,startframe,endframe)
        params, params_covariance = optimize.curve_fit(exponential_function, Is[:,0], Is[:,1], p0=[1.3, fitsign*0.04,0.1])
        Rs[pos,1] = params[1]
        Rs[pos,0] = convert_pixel_to_q(i+(delx)/2,(ysliceend-yslicestart)/2,ccdx,ccdz,x0,y0,ccdx0,ccdz0,d,E,alpha_i)[1]#(i+(delx)/2) #extract q... select the center pixel of the box to calculate q
        pos +=  1
    return(Rs)    
    
def stitch(a,b):
    """
    used to stich the R values generated by read_and_get_R_along_q_parallel
    """
    c = np.zeros((len(a)+len(b),2))
    c[0:len(a),:] = a
    c[len(a):,:]=b
    return(c)

def plotRs(Rs,yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe):
    plt.figure(figsize=(7,7))
    plt.scatter(Rs[:,0], Rs[:,1]*1E3, color = 'black', marker = '.', label='R')
    plt.ylabel(r'R($q_{y}$) ($10^{-3} s^{-1}$)')
    plt.xlabel(r'q ($nm^{-1}$)')
    plt.legend(loc='best')
    plt.suptitle('yslicestart = {0}'.format(yslicestart) + 
                 ' \n ysliceend = {0}'.format(ysliceend) + 
                 ' \n xslicestart = {0}'.format(xslicestart) + 
                 ' \n xsliceend = {0}'.format(xsliceend) + 
                 ' \n delx = {0}'.format(delx) + 
                 ' \n startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe) , fontsize=8)
    
    plt.show()
    plt.savefig(' yslicestart = {0}'.format(yslicestart) + 
                 ' ysliceend = {0}'.format(ysliceend) + 
                 ' xslicestart = {0}'.format(xslicestart) + 
                 ' xsliceend = {0}'.format(xsliceend) +
                '  delx = {0}'.format(delx) + 
                 ' startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe)) 
    
def plotRs_manuallabel(Rs,text,ylim):
    plt.figure(figsize=(7,7))
    plt.scatter(Rs[:,0], Rs[:,1]*1E3, color = 'black', marker = '.', label='R')
    plt.ylabel(r'R($q_{y}$) ($10^{-3} s^{-1}$)')
    plt.xlabel(r'q ($nm^{-1}$)')
    #plt.ylim(bottom=-10.5)
    plt.legend(loc='best')
    plt.suptitle(text, fontsize=8)
    
    plt.show()
    print('the figure is not saved automatically. Please save manually')
            
############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxx     END     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################

def convert_pixel_to_q(x,y,ccdx,ccdz,x0,y0,ccdx0,ccdz0,d,E,alpha_i):
    """Notice that x comes first here. sorry for inconsistency. normally Y is taken first because of 
    the way data is sliced. but for computing I would like to keep things in x and y order"""
    Psi = ((ccdx0 - ccdx)*(1E-3) + (x - x0)*(55E-6)) / (d*1E-3)
    alpha_f = ((ccdz - ccdz0)*(1E-3) + (y0 - y)*(55E-6)) / (d*1E-3)
    h = 6.62607 * 1E-34
    c = 3E8
    lambda_ = h*c/E
    q_x = ((2*np.pi)/lambda_)*(np.cos(alpha_f)*np.cos(Psi)-np.cos(alpha_i))
    q_y = ((2*np.pi)/lambda_)*(np.cos(alpha_f)*np.sin(Psi))
    q_z = ((2*np.pi)/lambda_)*(np.sin(alpha_i)+np.sin(alpha_f))
    
    #convert m^-1 to nm^-1
    q_x,q_y,q_z = q_x*1E-9,q_y*1E-9,q_z*1E-9
    #calculate q magnitude
    q = (q_x**2+q_y**2+q_z**2)**(0.5)
    
    return(q_x,q_y,q_z,q) 