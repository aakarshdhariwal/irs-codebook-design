import numpy as np
import matplotlib.pyplot as plt
import pdb

############ To Do###########
## make distance vectro as a constant value to evaluate the user ata fixed point and 
## the use the theta and A_x and A-y thing from codebook_design code and plot the beam shape with radius also fied point.
#############################
#### Calculate and implement the matched filter influenced IRS gain management.
#### first we initialize all the required parameters such as IRS size, wavelength, position vectors such as p_r,p_n and p_i 
#### namely target user, IRS elements and point source BS. 
### In order to define the new phase shifts in the codebook design we implement artificial codeword design that 
#### emulates fixed user and finds the target in a region.


###Function to calculate scalr gain g_m for various target positions over the spherical grid.

def discrete_phase_shift(b,w_n):

    mul=2*np.pi/(2**b)
    phase=np.mod(w_n,2*np.pi)

    #phase=mul*np.mod(np.round_(w_n/mul),2*np.pi)
    #phase=np.pi*np.remainder(w_n,2)

    #w_n=phase
    w_n=mul*np.round_(phase/mul)
   
    if w_n==0:
        w_n=0
    elif w_n==np.pi:
        w_n=8*np.pi/9
    elif w_n==(np.pi/2):
        w_n=4*np.pi/9
    elif w_n==(1.5*np.pi):
        w_n=12*np.pi/9
    else:     
        w_n
      

    #print(w_n)
    return w_n


def calculate_gain(Q,wavelength,p_n,p_i,p_t,b,p_r):
    
    res=1
    #wave number
    k=1j*2*np.pi/wavelength
    d_x=wavelength/2
    d_y=wavelength/2
    g_bar=4*np.pi*d_x*d_y/(wavelength**2) # unit cell factor
    # distance between two concentric spheres radius for the receiver of point of interest
    #d=np.random.uniform(r,r+10,1)
    #p_r=np.array([d*np.sin(0),0,d*np.cos(0)])
    g_m=np.zeros((1,res),dtype='complex')
    for j in range(res):
        noise=np.random.uniform(-np.pi/12,np.pi/12,Q*Q)
        #noise=np.pi/3
        gaussian_noise=np.random.normal(0,0.1008,Q*Q)
        #g_m=np.zeros((1,100),dtype='complex')
        for i in range(Q*Q):
            w_n=(2*np.pi/wavelength)*(np.linalg.norm(p_i-p_n[i])+np.linalg.norm(p_t-p_n[i]))

            w_n=discrete_phase_shift(b,w_n)
            w_n_noisy=w_n+noise[i]
            # For Gaussian noise in the beam
            w_n_gaussian_noise=w_n+gaussian_noise[i]

            #g_m+=np.exp(k*(np.linalg.norm(p_i-p_n[i])+np.linalg.norm(p_r-p_n[i])))*np.exp(-1j*w_n_noisy)
            g_m[0,j]+=np.exp(k*(np.linalg.norm(p_i-p_n[i])+np.linalg.norm(p_r-p_n[i])))*np.exp(-1j*w_n_gaussian_noise)
        #g_m+=np.exp(k*(np.linalg.norm(p_i-p_n[i])+np.linalg.norm(p_r-p_n[i])))*np.exp(-k*(np.linalg.norm(p_i-p_n[i])+np.linalg.norm(p_t-p_n[i])))
    g_m=np.mean(g_m)
    gain_m=20*np.log10(np.abs(g_bar*g_m))
    #g[0,j]=g_m
    
    return gain_m


if __name__ == '__main__':
    #initialize parameters
    frequency=int(3e10)
    c=int(3e8)
    wavelength=c/frequency
    #wavelength=1
    # radius of the concentric crcles we need to evaluate for designing the phase shifts
    radius=np.array([30])*wavelength
    # direction and vector to define the evaluation distance of the users
    #distance_vec=np.linspace(0,100,100,dtype='int')*wavelength
    distance_vec=np.array([30])*wavelength
    # number of discrete points that define the spherical grid 
    
    Q=20
    M_x=30
    M_y=30
    theta,phi=np.mgrid[0.0:2*np.pi:30j,0.0:2*np.pi:30j] #20==M
    #theta,phi=np.mgrid[-np.pi/2:np.pi/2:20j,-np.pi/2:np.pi/2:20j] 
#####################################################################
    #theta[0:9,0:9]=0
    #theta[10:19,10:19]=np.pi
######################################################################

    ######################################
    res=int(1e2)
    theta_r=np.linspace(0,(np.pi/2),res)
    phi_r=np.array([0,np.pi])
    #Psi_i=np.array([theta_i,phi_i],dtype='object')
    Psi_r=np.array([theta_r,phi_r],dtype='object')
    #A_x_psi_i=np.array([np.sin(Psi_i[0])*np.cos(Psi_i[1])])

   # azimuth reflection angle vector from IRS to  UE
    A_x_psi_r=np.array([np.sin(Psi_r[0])*np.cos(j) for j in Psi_r[1] ])
    A_x_psi_r=np.concatenate((A_x_psi_r[1,::-1],A_x_psi_r[0,:]))
    #print(A_x_psi_r.shape)

    A_y_psi_r=np.array([np.sin(Psi_r[0])*np.sin(j) for j in Psi_r[1]])
    A_y_psi_r=np.concatenate((A_y_psi_r[1,::-1],A_y_psi_r[0,:]))
    #A_y_psi_r=A_y_psi_r.reshape((phi_r.shape[0]*res,))
    theta_rr=np.linspace(-np.pi/2,(np.pi/2),2*res)
    ###########################################################################
    k=1j*2*np.pi/wavelength
    # IRS element spacing
    d_x=wavelength/2
    d_y=wavelength/2
    p_n=np.zeros((Q*Q,3))
    #initialize BS as point source vector
    p_i=np.array([0,0,1000*wavelength])
    #origin of the IRS 
    p_n_origin=np.array([d_x/2,d_y/2,0])
    for i in range(0,Q): 
        for j in range(0,Q):

            p_n[i*Q+j,0]=p_n_origin[0]+(i-10)*d_x
            p_n[i*Q+j,1]=p_n_origin[1]+(j-10)*d_y
     
    bits=np.array([1,2,3,10])
    #gain=np.zeros((bits.shape[0],radius.shape[0],distance_vec.shape[0]))

    gain=np.zeros((bits.shape[0],2*res))
    #gain=np.zeros((1,res))
    # the position matrix is a variable that stores for each sphere M_x*M_y total grid vectors in[x,y,z] plane
    #p_r_artificial=np.zeros((radius.shape[0],M_x*M_y,3))
    p_r_artificial=np.zeros((radius.shape[0],3,theta.shape[0],theta.shape[1]))
   
    # construction of p_r positional vector whoese elements lie on the spherical grid
    
    for b in range(bits.shape[0]):
        
                x=radius[0]*np.sin(theta)*np.cos(phi) #x
                # might need to change y and z
                #x=np.zeros_like(x)
                y=radius[0]*np.sin(theta)*np.sin(phi) #y
                #y=np.zeros_like(y)
                #theta=np.zeros_like(theta)
                z=radius[0]*np.cos(theta) #z
                #z=radius[r]*np.ones_like(z)
                p_r_artificial[0]=np.array([x,y,z],dtype='object')
                p_target=p_r_artificial[0].reshape(3,M_x,M_y)
                
                for i_theta in range(A_x_psi_r.shape[0]):

                    ########################################
                    ### While calculating gain for a user we inherently fix a target point which is defined 
                    ## using p_r_artificial and then defining just a target vector for the codebook design scenario.
                    ## so basically when the target is moving and user/ receiver is fixed we compute gain and plot it against the distances 
                    ## this analysis works out for my every plot which I failed to understnd completely.
                    ## Now the goal is to change the direction of the users keeping fixed distance and fixed target point.
                    p_t=[p_target[0,0,0],p_target[1,0,0],p_target[2,0,0]]
                    #p_r=np.array([distance_vec[0]*np.sin(A_x_psi_r[i_theta]),0,distance_vec[0]*np.cos(A_x_psi_r[i_theta])])
                    p_r=np.array([distance_vec[0]*A_x_psi_r[i_theta],distance_vec[0]*A_y_psi_r[i_theta],distance_vec[0]*np.cos(theta_rr[i_theta])])
                    gain[b,i_theta]=calculate_gain(Q,wavelength,p_n,p_i,p_t,bits[b],p_r)
                    #########################################
            
        
            
                # for n_x in range(M_x):
                #     for n_y in range(M_y):
                #             #p_target=[x[n_x,n_y],y[n_x,n_y],z[n_x,n_y]]
                #                 p_t=[p_target[0,n_x,n_y],p_target[1,n_x,n_y],p_target[2,n_x,n_y]]
                #                 gain[0,n_x,n_y]=calculate_gain(Q,wavelength,p_n,p_i,p_t,distance_vec[0],bits[b])
                               
               
            #pdb.set_trace()

    #print(gain.shape)
    #gain_30_gaussian_noise=np.save('gain_30_gaussian_noise.npy',gain)
    #gain_30_noisy=np.save('gain_30_noisy.npy',gain)
    #gain_20_degrees_noisy_averaged=np.save('gain_20_degrees_noisy_averaged.npy',gain)
    #gain_r=gain.reshape(1,M_x*M_y)
    #gain_1bit_100res=np.save('gain_1bit_100res.npy',gain)
    #plot gain vs the distace vector
    x_1=28.288
    #y_1=gain[0,0,28]
    x_2=1.0
    #y_2=gain[1,0,1]
    #theta_r=theta.reshape(1,400)
    theta_r=np.linspace(-90,90,2*res)
    fig,ax=plt.subplots()
    
    #ax.plot(theta_r,gain.T)
    ax.plot(theta_r,gain[0].T,label=f'Gain for radius={int(radius[0]*100)} m and bits={bits[0]}')
    arrowprops={'arrowstyle': '-', 'ls':'--'}
    
   
    ax.plot(theta_r,gain[1].T,ls='--',label=f'Gain for radius={int(radius[0]*100)} m and bits={bits[1]}')
    
    
    ax.plot(theta_r,gain[2].T,ls='-.',label=f'Gain for radius={int(radius[0]*100)} m and bits={bits[2]}')
    
    ax.plot(theta_r,gain[3].T,ls=':',label=f'Gain for radius={int(radius[0]*100)} m and bits={bits[3]}')
    

    #ax.set_title(" IRS Response gain vs angle of reflection in azimuth")
    ax.set_xlabel(r'$\theta_r $',fontsize=12)
    #ax.set_xlabel(r'$\dfrac{distance_{vec}}{\lambda}$',fontsize=14)
    #ax.set_xlabel(r' Distance from the IRS')
    ax.set_ylabel(r'$10\log{||g_{ris}(\mathbf{p_i},\mathbf{p_r})||^2}$',fontsize=12)
    ax.set_xticks(np.arange(-90,100,10))
    #ax.set_xticks(np.arange(0,40,1))
    #ax.set_yticks(np.arange(0,100,10))
    #ax.set_ylim(-20,100
    ax.set_xlim(-90,90)
    ax.grid()
    ax.legend(loc="best")
   
    plt.show()