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
    elif w_n==np.pi/2:
        w_n=4*np.pi/9
    elif w_n==(1.5*np.pi):
        w_n=12*np.pi/9
    else:
        w_n=w_n
    
      

    #print(w_n)
    return w_n


def calculate_gain(Q,wavelength,p_n,p_i,p_t,d,b):
    
    #wave number
    k=1j*2*np.pi/wavelength
    d_x=wavelength/2
    d_y=wavelength/2
    g_bar=4*np.pi*d_x*d_y/(wavelength**2) # unit cell factor
    # distance between two concentric spheres radius for the receiver of point of interest
    #d=np.random.uniform(r,r+10,1)
    # maybe express p_r also using theta and phi
    p_r=np.array([d*np.sin(0),0,d*np.cos(0)])
    # Adding noise to the phase shifts using unifrom distribution and gaussiaon noise
    noise=np.random.uniform(-np.pi/12,np.pi/12,Q*Q)
    gaussian_noise=np.random.normal(0,0.1008,Q*Q)
    g_m=0
    for i in range(Q*Q):
        w_n=(2*np.pi/wavelength)*(np.linalg.norm(p_i-p_n[i])+np.linalg.norm(p_t-p_n[i]))

        #w_n=discrete_phase_shift(b,w_n)
        w_n_noisy=w_n+noise[i]
        ## For Gaussian noise in the beam
        w_n_gaussian_noise=w_n+gaussian_noise[i]

        g_m+=np.exp(k*(np.linalg.norm(p_i-p_n[i])+np.linalg.norm(p_r-p_n[i])))*np.exp(-1j*w_n)
        
        #g_m+=np.exp(k*(np.linalg.norm(p_i-p_n[i])+np.linalg.norm(p_r-p_n[i])))*np.exp(-k*(np.linalg.norm(p_i-p_n[i])+np.linalg.norm(p_t-p_n[i])))
    g_m=20*np.log10(np.abs(g_bar*g_m))
    #g[0,j]=g_m
    
    return g_m


if __name__ == '__main__':
    #initialize parameters
    frequency=int(3e10)
    c=int(3e8)
    wavelength=c/frequency
    #wavelength=1
    # radius of the concentric crcles we need to evaluate for designing the phase shifts
    radius=np.array([10,20,30,60,100])*wavelength
    # direction and vector to define the evaluation distance of the users
    distance_vec=np.linspace(0,100,1000,dtype='int')*wavelength
   
    # number of discrete points that define the spherical grid 
    
    Q=20
    M_x=20
    M_y=30
    theta,phi=np.mgrid[0.0:2*np.pi:20j,0.0:2*np.pi:30j] #20==M
    #theta,phi=np.mgrid[-np.pi/2:np.pi/2:20j,-np.pi/2:np.pi/2:20j] 
    
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
     
    bits=np.array([0])
    gain=np.zeros((bits.shape[0],radius.shape[0],distance_vec.shape[0]))

    
    
    # the position matrix is a variable that stores for each sphere M_x*M_y total grid vectors in[x,y,z] plane
    #p_r_artificial=np.zeros((radius.shape[0],M_x*M_y,3))
    p_r_artificial=np.zeros((radius.shape[0],3,theta.shape[0],theta.shape[1]))
   
    # construction of p_r positional vector whoese elements lie on the spherical grid
    
    for b in range(bits.shape[0]):
        
            for r in range(radius.shape[0]):
         
                x=radius[r]*np.sin(theta)*np.cos(phi) #x
                # might need to change y and z
                #x=np.zeros_like(x)
                y=radius[r]*np.sin(theta)*np.sin(phi) #y
                #y=np.zeros_like(y)
                #theta=np.zeros_like(theta)
                z=radius[r]*np.cos(theta) #z
                #z=radius[r]*np.ones_like(z)
                p_r_artificial[r]=np.array([x,y,z],dtype='object')
                p_target=p_r_artificial[r].reshape(3,M_x,M_y)
                
                #for i_theta in range(A_x_psi_r.shape[0]):
            
                for d in range(distance_vec.shape[0]):
                
                    #p_r=np.array([d*np.sin(A_x_psi_r[i_theta]),0,d*np.cos(A_x_psi_r[i_theta])])
                    n_x=0
                    n_y=1
                    #for n_x in range(M_x):
                     #   for n_y in range(M_y):
                            #p_target=[x[n_x,n_y],y[n_x,n_y],z[n_x,n_y]]
                    p_t=[p_target[0,n_x,n_y],p_target[1,n_x,n_y],p_target[2,n_x,n_y]]
                    gain[b,r,d]=calculate_gain(Q,wavelength,p_n,p_i,p_t,distance_vec[d],bits[b])
                                
            #pdb.set_trace()

    #print(gain.shape)
    #plot gain vs the distace vector
    x_1=28.288
    y_1=63
    x_2=1.0
    #y_2=gain[1,0,1]
    #gain_linear_extension=np.save('gain_linear_extension.npy',gain)
    #gain_10_bits=np.save('gain_10_bits',gain)
    fig,ax=plt.subplots()
    
    ax.plot(distance_vec/wavelength,gain[0,0].T,label=f'Gain for radius={int(radius[0]*100)} m ')
    arrowprops={'arrowstyle': '-', 'ls':'--'}
    
    plt.annotate(r'$d_B$', xy=(x_1,y_1), xytext=(x_1, 0), 
             textcoords=plt.gca().get_xaxis_transform(),
             arrowprops=arrowprops,
             va='top', ha='center')
    
    ax.plot(distance_vec/wavelength,gain[0,1].T,ls='--',label=f'Gain for radius={int(radius[1]*100)} m ')
    
    plt.annotate(r'$d_B$', xy=(x_1,y_1), xytext=(x_1, 0), 
             textcoords=plt.gca().get_xaxis_transform(),
             arrowprops=arrowprops,
             va='top', ha='center')
    ax.plot(distance_vec/wavelength,gain[0,2].T,ls='-.',label=f'Gain for radius={int(radius[2]*100)} m ')
    
    plt.annotate(r'$d_B$', xy=(x_1,y_1), xytext=(x_1, 0), 
             textcoords=plt.gca().get_xaxis_transform(),
             arrowprops=arrowprops,
            va='top', ha='center')

    ax.plot(distance_vec/wavelength,gain[0,3].T,label=f'Gain for radius={int(radius[3]*100)} m ')
    
    plt.annotate(r'$d_B$', xy=(x_1,y_1), xytext=(x_1, 0), 
              textcoords=plt.gca().get_xaxis_transform(),
              arrowprops=arrowprops,
             va='top', ha='center')
    
    ax.plot(distance_vec/wavelength,gain[0,4].T,label=f'Gain for radius={int(radius[4]*100)} m ')
    
    plt.annotate(r'$d_B$', xy=(x_1,y_1), xytext=(x_1, 0), 
             textcoords=plt.gca().get_xaxis_transform(),
             arrowprops=arrowprops,
            va='top', ha='center')
    
   

    #ax.set_title(" IRS Response gain vs DIstance from the center of the IRS")
    #ax.set_xlabel(r'$\theta_r$')
    ax.set_xlabel(r'$\dfrac{d}{\lambda}$',fontsize=12)
    #ax.set_xlabel(r' Distance from the IRS')
    ax.set_ylabel(r'$10\log{||g_{ris}(\mathbf{p_i},\mathbf{p_r})||^2}$',fontsize=12)
    #ax.set_xticks(np.arange(-90,100,10))
    #ax.set_xticks(np.arange(0,40,1))
    ax.set_yticks(np.arange(0,80,10))
    ax.set_ylim(0,80)
    ax.set_xlim(0,100)
    ax.grid()
    ax.legend(loc="best")
   
    plt.show()