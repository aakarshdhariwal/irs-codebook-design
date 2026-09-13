import numpy as np
import matplotlib.pyplot as plt
import pdb
#plt.rcParams['text.usetex'] = True

# function to calculate the gain of the IRS response function
# To do: implement DFT, Linear and quadratic codebook designs and save each implementation results in a file to load that
# to compare results afterwards.
def calculate_g_m(A_x_BS_to_UE,A_y_BS_to_UE,wavelength,Q,m_x,m_y):
    #wave number
    k1=1j*2*np.pi/wavelength
    k2=1j*2*np.pi/(Q)
    n_x=np.linspace(0,Q-1,Q)
    n_y=np.linspace(0,Q-1,Q)
    # antenna spacing factors 
    d_x=wavelength/2
    d_y=wavelength/2
    g_bar=4*np.pi*d_x*d_y/(wavelength**2) # unit cell factor
    #Calculating the IRS response gain as a scalar value given the formula from Jamali's Power efficiency, overhead
    # and complexity tradeoff design for quadratic phase shift design.
    temp=np.sum(np.exp(-k2*n_x*m_x+k1*d_x*n_x*A_x_BS_to_UE)*np.sum(np.exp(-k2*n_y*m_y+k1*d_y*n_y*A_y_BS_to_UE)))
        
    g_m=20*np.log10(np.abs(g_bar*temp))
    
    return g_m



def calculate_IRS_gain(A_x_psi_i,A_x_psi_r,A_y_psi_i,A_y_psi_r,phase_shift,wavelength,Q):
    #wave number
    k1=1j*2*np.pi/wavelength
    # antenna spacing factors 
    d_x=wavelength/2
    d_y=wavelength/2
    g_bar=4*np.pi*d_x*d_y/(wavelength**2) # unit cell factor
    # mutlipath taps from BS to IRS
    L_i=1
    a_i=np.array([np.exp(k1*(d_x*A_x_psi_i*n_x+d_y*A_y_psi_i*n_y)) \
     for n_x in range(Q) for n_y in range(Q)])
    # multipath from IRS to UE
    L_r_k=1
    d_r_h=np.array([np.exp(k1*(d_x*A_x_psi_r*n_x+d_y*A_y_psi_r*n_y)) \
    for n_x in range(Q) for n_y in range(Q)])
    d_r_h=d_r_h.reshape((L_r_k,Q*Q))
    #g_m=np.dot(np.dot(d_r_h,phase_shift),a_i)
    g_m=np.sum((d_r_h*phase_shift)*a_i)
    g_m=20*np.log10(np.abs(g_bar*g_m))

    return g_m

def calculate_near_field_gain(A_x_psi_i,A_x_psi_r,A_y_psi_i,A_y_psi_r,phase_shift,wavelength,Q):
    # Initialize near field parameters such as position vectors p_n and p_r.
    k=1j*2*np.pi/wavelength
    # IRS element spacing
    d_x=wavelength/2
    d_y=wavelength/2
    g_bar=4*np.pi*d_x*d_y/(wavelength**2) # unit cell factor
    #distance to user d in meters
    d=1000*wavelength
    p_r=np.array([d*np.sin(A_x_psi_r),0,d*np.cos(A_x_psi_r)])
    p_n=np.zeros((Q*Q,3))
    #origin of the IRS 
    p_n_origin=np.array([d_x/2,d_y/2,0])
    for i in range(0,Q):
        for j in range(0,Q):
            # TODO think about this offset
            p_n[i*Q+j,0]=p_n_origin[0]+(i-10)*d_x
            p_n[i*Q+j,1]=p_n_origin[1]+(j-10)*d_y
    #initialize BS as point source vector
    #phase_shift[:] = 1
    p_i=np.array([0,0,1000*wavelength])
    g_m=0
    for n in range(Q*Q):
        g_m+=np.exp(k*(np.linalg.norm(p_i-p_n[n],ord=2)+np.linalg.norm(p_r-p_n[n],ord=2)))*phase_shift[n,n]
        # pdb.set_trace()
        #print(k*(np.linalg.norm(p_i-p_n[n],ord=2)+np.linalg.norm(p_r-p_n[n],ord=1)))
    g_m=20*np.log10(np.abs(g_bar*g_m))
    return g_m


def calculate_gain(codebook):

    # Calculate the IRS response g_m for the given IRS with QxQ IRS elements.
    # Define inputs angle of incident of arrival and angle of departure psi_i and psi_r respectively.
    # Here, we make M the size of DFT codebook as M=Q^2 and calculate the gain for varying theta_r i.e angle of reflection in elevation direction.dir
    # We define Inputs as:
    # -resolution: number of points to fine tune our gain against elevation reflection angle.
    # -theta_i: elevation angle of incidence
    # -psi_i: azimuth angle of incidence
    # -theta_r: elevation angle of reflection
    # -psi_r: azimuth angle of reflection
    # output: g_m IRS response gain for varying theta_r over the m_x --> 0.....Q-1 and m_y=1.

    #basic initializations
    Q=20
    frequency=int(3e10)
    c=int(3e8)
    wavelength=c/frequency
    #wavelength=1
    #resoluton for points theta_r
    res = int(1e2)
    # initialization of elevation and azimuth angles respectively for incidence and reflection.
    theta_i=0
    phi_i=0
    phi_r=np.array([0,np.pi])
    theta_r=np.linspace(0,(np.pi/2),res)
   # theta_r[:] = 0 # TODO remove
    #DFT codebook  indices
    m_x=1
    m_y=0
    # initialize IRS response function
    M = 5
    g_m=np.zeros((M,phi_r.shape[0]*res))
    # defining angle of arrival and angle of departure

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
    #print(A_y_psi_r.shape)
    #print(A_x_psi_i.shape, A_x_psi_r.shape, A_y_psi_i.shape, A_y_psi_r.shape)

    for m in range(M):
        if codebook=='DFT':
            k2=1j*2*np.pi/Q # wave number
            phase_shift=np.array([np.exp(-k2*(n_x*m+n_y*m_y)) for n_x in range(Q) for n_y in range(Q)])
            phase_shift=np.diag(phase_shift)

        elif codebook=='Linear':
            # phase shift gradient to model variable unit cell spacing between elements for a linear codebook 
            beta_bar=2
            M_x=5
            M_y=5
            d_x=wavelength/2
            d_y=wavelength/2
            k3=1j*2*np.pi*beta_bar/wavelength # wave number
            #phase shift calculation 
            phase_shift=np.array([np.exp(-k3*(d_x*n_x*m/(M_x)+d_y*n_y*m_y/(M_y))) \
            for n_x in range(Q) for n_y in range(Q)])
            phase_shift=np.diag(phase_shift)

        elif codebook=='Quadratic':
            beta_bar_x=2
            beta_bar_y=2
            M_x=5
            M_y=5 #Q/4
            # phase shift gradient offset udate terms from equation 11
            delta_beta_x=beta_bar_x/M_x
            delta_beta_y=beta_bar_y/M_y
            beta_bar_x_mx=m*delta_beta_x-1
            beta_bar_y_my=m_y*delta_beta_y
            d_x=wavelength/2
            d_y=wavelength/2

            k3=1j*2*np.pi/wavelength # wave number
            # intermediate variables to initialize delta beta 
            quadratic_x=delta_beta_x/(2*Q)
            quadratic_y=delta_beta_y/(2*Q)
            # phase shift design for quadratic codebook
            # import pdb
            # pdb.set_trace()
            phase_shift=np.array([np.exp(-k3*(d_x*(quadratic_x*(n_x**2)+beta_bar_x_mx*n_x)+d_y*(quadratic_y*(n_y**2)+beta_bar_y_my*n_y))) \
            for n_x in range(Q) for n_y in range(Q)])
            phase_shift=np.diag(phase_shift)

        for i_theta in range(A_x_psi_r.shape[0]):
                
        # beam steering angle vectors form BS to Ue
                A_x_BS_to_UE=A_x_psi_i+A_x_psi_r[i_theta]
                A_y_BS_to_UE=A_y_psi_i+A_y_psi_r[i_theta]
                #g_m[m,i_theta]=calculate_IRS_gain(A_x_psi_i,A_x_psi_r[i_theta],A_y_psi_i,A_y_psi_r[i_theta],phase_shift,wavelength,Q)
                g_m[m,i_theta]=calculate_near_field_gain(A_x_psi_i,A_x_psi_r[i_theta],A_y_psi_i,A_y_psi_r[i_theta],phase_shift,wavelength,Q)
    #print(g_m.shape)
    #print(g_m[0,0:100])
    
    return g_m

if __name__ == '__main__':
    #initialize parameters
    res=int(1e2)
    theta=np.linspace(-90,90,2*res)
    #codebook=['DFT','Linear','Quadratic']
    
    g_m=calculate_gain(codebook='Linear')
    #linear_array=np.save('Linear_array.npy',g_m)
    #g_m_Linear=npload('Linear_array.npy')
    fig,ax=plt.subplots()
    
    ax.plot(theta,g_m.T)
    ax.set_title(" IRS Response gain vs angle of reflection in azimuth")
    ax.set_xlabel('theta_r')
    ax.set_ylabel('IRS response gain')
    ax.set_xticks(np.arange(-90,100,10))
    ax.set_yticks(np.arange(0,80,10))
    ax.set_ylim(0,80)
    ax.grid()
    plt.show()
    