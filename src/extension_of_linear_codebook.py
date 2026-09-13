import numpy as np
import matplotlib.pyplot as plt
import pdb
#### Calculate and implement the matched filter influenced IRS gain management.
#### first we initialize all the required parameters such as IRS size, wavelength, position vectors such as p_r,p_n and p_i 
#### namely target user, IRS elements and point source BS. 
### In order to define the new phase shifts in the codebook design we implement artificial codeword design that 
#### emulates fixed user and finds the target in a region.


###Function to calculate scalr gain g_m for various target positions over the spherical grid.

def calculate_gain(Q,r,wavelength,p_n,p_i,phase_shifts,A_x_psi_r):
    #wave number
    k=1j*2*np.pi/wavelength
    d_x=wavelength/2
    d_y=wavelength/2
    g_bar=4*np.pi*d_x*d_y/(wavelength**2) # unit cell factor
    # distance between two concentric spheres radius for the receiver of point of interest
    distance=np.linspace(0,1000,1001)*wavelength
    g=np.zeros(distance.shape)
    for j in range(distance.shape[0]):
        p_r=np.array([distance[j]*np.sin(A_x_psi_r),0,distance[j]*np.cos(A_x_psi_r)])
        g_m=0
        for i in range(Q*Q):
            g_m+=np.exp(k*(np.linalg.norm(p_i-p_n[i])+np.linalg.norm(p_r-p_n[i])))*np.exp(1j*phase_shifts[i])
        g_m=20*np.log10(np.abs(g_bar*g_m))
        g[j]=g_m
    return g


def gain_and_phase_shift_initializer(Q):

    frequency=int(3e10)
    c=int(3e8)
    wavelength=c/frequency
    # radius of the concentric crcles we need to evaluate for designing the phase shifts
    radius=np.array([100])*wavelength
    # number of discrete points that define the spherical grid 
    M_x=20
    M_y=20
    theta,phi=np.mgrid[0.0:2*np.pi:20j,0.0:2*np.pi:20j] #20==M
    #phase_shifts=np.zeros((radius.shape,Q*Q,Q*Q))
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
     #resoluton for points theta_r
    res = int(1e1)
    # initialization of elevation and azimuth angles respectively for incidence and reflection.
    theta_i=0
    phi_i=0
    phi_r=np.array([0,np.pi])
    theta_r=np.linspace(0,(np.pi/2),res)
    theta_r[:]=0

    Psi_i=np.array([theta_i,phi_i],dtype='object')
    Psi_r=np.array([theta_r,phi_r],dtype='object')
    A_x_psi_i=np.array([np.sin(Psi_i[0])*np.cos(Psi_i[1])])

   # azimuth reflection angle vector from IRS to  UE
    A_x_psi_r=np.array([np.sin(Psi_r[0])*np.cos(j) for j in Psi_r[1] ])
    A_x_psi_r=np.concatenate((A_x_psi_r[1,::-1],A_x_psi_r[0,:]))

    #incidence azimuth angle vector from bS to IRS
    A_y_psi_i=np.array([np.sin(Psi_i[0])*np.sin(Psi_i[1])])
   # azimuth reflection angle vector from IRS to  UE

    A_y_psi_r=np.array([np.sin(Psi_r[0])*np.sin(j) for j in Psi_r[1]])
    A_y_psi_r=A_y_psi_r.reshape((phi_r.shape[0]*res,))

    #gain=np.zeros((radius.shape[0],A_x_psi_r.shape[0]))
    gain=np.zeros((radius.shape[0],1001))
    p_r_artificial=np.zeros((radius.shape[0],M_x*M_y,3))
    phase_shifts=np.zeros((radius.shape[0],Q*Q))
    # construction p_r positional vector whoese elements lie on the spherical grid
    for r in range(radius.shape[0]):
         
         x=radius[r]*np.sin(theta)*np.cos(phi)
         x=np.zeros_like(x)
         # might need to change y and z
         y=radius[r]*np.sin(theta)*np.sin(phi)
         z=radius[r]*np.cos(theta)
         p_r_artificial[r]=np.array([x,y,z]).reshape(M_x*M_y,3)
         #pdb.set_trace()
         #create phase shift across points on the spherical grid
         for n in range(M_x*M_y):
            phase_shifts[r,n]=-2*np.pi*(np.linalg.norm(p_i-p_n)+np.linalg.norm(p_r_artificial[r,n]-p_n))/wavelength
         
         for  i_theta in range(A_x_psi_r.shape[0]):
            gain[r]=calculate_gain(Q,radius[r],wavelength,p_n,p_i,phase_shifts[r],A_x_psi_r[i_theta])
    #pdb.set_trace()
    return gain



if __name__ == '__main__':
    #initialize parameters
    frequency=int(3e10)
    c=int(3e8)
    wavelength=c/frequency
    #gain calu
    gain=gain_and_phase_shift_initializer(Q=20) 
    
    #plot gain vs the radius
    fig,ax=plt.subplots()
    #distance to evaluate the gain from the IRS to receivers
    distance=np.linspace(0,1000,1001)*wavelength
    #distance=np.mgrid[0:200:200j]
    ax.plot(distance/wavelength,gain[0].T,label='Gain for radius_10')
    ax.plot(distance/wavelength,gain[1].T,label='Gain for radius_20')
    ax.plot(distance/wavelength,gain[2].T,label='Gain for radius_30')
    ax.plot(distance/wavelength,gain[3].T,label='Gain for radius_100')
    

   
    
    ax.set_title(" IRS Response gain vs angle of reflection in azimuth")
    #ax.set_xlabel(r'$\theta_r$')
    #ax.set_xlabel(r'$\dfrac{d}{\lambda}$')
    ax.set_xlabel(r' Distance from the IRS')
    ax.set_ylabel(r'$||g_{ris}(\mathbf{p_i},\mathbf{p_r})||^2$')
    #ax.set_xticks(np.arange(-90,100,10))
    #ax.set_xticks(np.arange(0,40,1))
    ax.set_yticks(np.arange(0,100,10))
    ax.set_ylim(0,100)
    ax.grid()
    ax.legend(loc="best")
   
    plt.show()


