"""
Last modified date 02/27/2019

@author: Peco
"""
def long_wave_form(q,Sy,B):
    return -B*(q**4)+Sy*(q**2)


def Rwithorchard(q,Sy,A,h0):
    Q = q*h0
    return Sy*(q**2)-A*((Q*(np.sinh(2*Q)-2*Q))/(1+2*(Q**2)+np.cosh(2*Q)))

def plotandfit_R(R,yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe):
    #Sy = 0
    params, params_covariance = optimize.curve_fit(long_wave_form, R[:,0], R[:,1], p0=[6.4,2]) #lambda q, B: long_wave_form(q, B, Sy), R[:,0], R[:,1], p0=[2]) 
    paramsOrchard, params_covarianceOrchard = optimize.curve_fit(Rwithorchard, R[:,0], R[:,1], p0=[0.5, 0.203,3.4]) 
    plt.figure(figsize=(7,7))
    plt.scatter(R[:,0], R[:,1]*1E3, color = 'black', marker = '.', label='data') #multiplied by 0.1*1e3/1e3 which mean it is per 1000 seconds
    plt.plot(R[:,0], Rwithorchard(R[:,0], paramsOrchard[0], paramsOrchard[1],paramsOrchard[2] )*1E3,label='Orchard Term fit')
    plt.plot(R[:,0], long_wave_form(R[:,0], params[0],params[1])*1E3,label=r'$S_{y}q^{2} - Bq^{4}$ fit') 
    plt.ylabel(r'R($q_{y}$) ($10^{-3} s^{-1}$)')
    plt.xlabel(r'q ($nm^{-1}$)')
    plt.legend(loc='best')
    plt.suptitle('yslicestart = {0}'.format(yslicestart) + 
                 ' ysliceend = {0}'.format(ysliceend) + 
                 ' \n xslicestart = {0}'.format(xslicestart) + 
                 ' xsliceend = {0}'.format(xsliceend) + 
                 ' delx = {0}'.format(delx) + 
                 ' \n startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe) + 
                 ' \n Sy, B = {0}'.format(params)+ 
                 ' Sy, A, h0 = {0}'.format(paramsOrchard), fontsize=8)

    plt.show()
    plt.savefig(DD.DD['filename'][:-4] + '_R_extracted')

def plotandfit_R_constraint(R,yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe):
    """
    here we constarint one parameter with brute force while fitting Rs
    for example, we can also set h0 as a certain value 
    """
    Sy = 0
    params, params_covariance = optimize.curve_fit(lambda q, B: long_wave_form(q, B, Sy), R[:,0], R[:,1], p0=[2]) 
    paramsOrchard, params_covarianceOrchard = optimize.curve_fit(lambda q, A, h0: Rwithorchard(q,Sy,A,h0), R[:,0], R[:,1], p0=[0.203,3.4]) 
    plt.figure(figsize=(7,7))
    plt.scatter(R[:,0], R[:,1]*1E3, color = 'black', marker = '.', label='data') #multiplied by 0.1*1e3/1e3 which mean it is per 1000 seconds
    plt.plot(R[:,0], Rwithorchard(R[:,0], Sy, paramsOrchard[0],paramsOrchard[1] )*1E3,label='Orchard Term fit')
    plt.plot(R[:,0], long_wave_form(R[:,0], params[0],Sy)*1E3,label=r'$S_{y}q^{2} - Bq^{4}$ fit') 
    plt.ylabel(r'R($q_{y}$) ($10^{-3} s^{-1}$)')
    plt.xlabel(r'q ($nm^{-1}$)')
    plt.legend(loc='best')
    plt.suptitle('yslicestart = {0}'.format(yslicestart) + 
                 ' ysliceend = {0}'.format(ysliceend) + 
                 ' \n xslicestart = {0}'.format(xslicestart) + 
                 ' xsliceend = {0}'.format(xsliceend) + 
                 ' delx = {0}'.format(delx) + 
                 ' \n startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe) + 
                 ' \n Sy = 0, B = {0}'.format(params)+ 
                 ' Sy = 0, A, h0 = {0}'.format(paramsOrchard), fontsize=8)

    plt.show()
    plt.savefig(DD.DD['filename'][-4] + '_R_fitted')