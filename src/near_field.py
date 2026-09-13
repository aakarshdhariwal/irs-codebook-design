### design the codebook for near field scenario where the distance from the antennas to the desired user is in 
# the range of 1lambda or 1 times the wavelength. Using the phase shift design we need to illuminate the 
# observation area with MU at the center. 

import numpy as np
import matplotlib.pyplot as plt


class near_field:

    ### Firstly initializing the constructor with parametres wavelength, frequency, aperture size of the RIS,
    # RIS gain factor g, total no. of RIS elements Q with Qz and Qy elemsts in z and y dimensions,
    #  angle of arrival and departures.
    def __init__(self):
        self.frequency=1e10
        self.c=3e8
        self.res=int(1e2)
        self.Qx=20
        self.Qy=20
        self.theta_i=0
        self.phi_i=0
        self.theta_r=np.linspace(0,np.pi/2,self.res)
        self.phi_r=np.array([0,np.pi])
        self.R_x=16 # in meters
        self.R_y=4 # in meters
        self.wavelength= self.c/self.frequency
        # L_z and L_y
        self.L_z=0.5 #in meters
        self.L_y=0.5
        self.aperture_size=20*self.wavelength
        
        #position vector of center of observation area
        self.p_b=np.array([30,30,1])
        #position vector for point source BS
        self.p_i=np.array([40, 5, 10])
        self.pris=np.array([0, 0, 0])
        self.d=self.wavelength/2
        self.g=4*np.pi*self.d**2/self.wavelength**2

    def phase_shift(self):
        ### calculate the angle of arrival and angle of departure for impinging waves
        Psi_i=np.array([self.theta_i,self.phi_i],dtype='object')
        Psi_r=np.array([self.theta_r,self.phi_r],dtype='object')
        A_x_psi_i=np.array([np.sin(Psi_i[0])*np.cos(Psi_i[1])])

        # azimuth reflection angle vector from IRS to  UE
        A_x_psi_r=self.aperture_size*np.array([np.sin(Psi_r[0])*np.cos(j) for j in Psi_r[1] ])
        A_x_psi_r=np.concatenate((A_x_psi_r[1,::-1],A_x_psi_r[0,:]))

        #incidence azimuth angle vector from bS to IRS
        A_y_psi_i=np.array([np.sin(Psi_i[0])*np.sin(Psi_i[1])])
        # azimuth reflection angle vector from IRS to  UE

        A_y_psi_r=self.aperture_size*np.array([np.sin(Psi_r[0])*np.sin(j) for j in Psi_r[1]])
        A_y_psi_r=A_y_psi_r.reshape((self.phi_r.shape[0]*self.res,))

        #Initialize the nth RIS position vector
        self.p_n=np.zeros(shape=(self.Qx*self.Qy,3))

        # Formulate position vector of observation area
        x=np.linspace(-self.R_x/2,self.R_x/2,self.Qx)
        y=np.linspace(-self.R_y/2,self.R_y/2,self.Qy)
        # z=np.zeros((self.Qx,1))
        #self.p_r=self.p_b+np.concatenate((x,y,z),axis=0)
        self.p_r_n=np.zeros_like(self.p_n)
        alpha=0.8

        g_ris=np.zeros((A_x_psi_r.shape[0],),dtype='complex_')
        for i_theta in range(A_x_psi_r.shape[0]):

            self.p_n[:,0]=A_x_psi_r[i_theta]
            #Mapping of p_n to the observation area
            delta_x=alpha*self.R_x/self.Qx
            delta_y=alpha*self.R_y/self.Qy
            M_p_n = np.zeros((self.Qx*self.Qy, 3))
            for i in range(self.Qx):
                for j in range(self.Qy):
                    self.p_r_n[i*self.Qy+j, 0]=x[i]
                    self.p_r_n[i*self.Qy+j, 1]=y[j]
                    x_w_x=i*self.R_x/self.Qx
                    y_w_y=j*self.R_y/self.Qy
                    M_p_n[i*self.Qy+j, 0] = (delta_x/self.L_z)*A_x_psi_r[i_theta]+x_w_x
                    M_p_n[i*self.Qy+j, 1] = (delta_y/self.L_y)*A_y_psi_r[i_theta]+y_w_y
            #print(M_p_n)
            #print(self.p_r_n)
            # x_w_x=np.linspace(0,self.Qx**2,self.Qx**2).reshape(self.Qx**2,1)
            # x_dim=(delta_x/self.L_z)*A_x_psi_r[i_theta]+x_w_x
            # y_w_y=np.linspace(0,self.Qy**2,self.Qy**2).reshape(self.Qy**2,1)
            # y_dim=(delta_y/self.L_y)*A_y_psi_r[i_theta]+y_w_y
            # z_dim=np.zeros((self.Qx**2,1))
            # self.p_tilda=np.concatenate((x_dim,y_dim,z_dim),axis=0)
            #mapping from p_n to observation area p_r_n
            self.map_p_n=self.p_r_n+M_p_n
            # intermediate calculations of phase shift terms
            x1=np.linalg.norm(self.map_p_n-self.p_n,ord=1,axis=1)
            x2=np.linalg.norm(self.map_p_n-self.pris,ord=1,axis=1)
            x3=np.linalg.norm(self.p_i-self.p_n,ord=1,axis=1)
            self.w_n=(-2*np.pi/self.wavelength)*(x1-x2+x3)
            #print(self.w_n.shape)
            for n in range(self.w_n.shape[0]):
                g_ris[i_theta]+=self.g*np.exp((1j*2*np.pi/self.wavelength)*(np.linalg.norm(self.p_i-self.p_n[n])+\
                np.linalg.norm(self.p_r_n[n]-self.p_n[n])))*np.exp(1j*self.w_n[n])
            
            
        return g_ris



if __name__ == '__main__':
    res=int(1e2)
    theta=np.linspace(-90,90,2*res)
    gris=near_field()
    ret = gris.phase_shift()
    g_ris=20*np.log10(np.abs(ret))
    fig,ax=plt.subplots()
    
    ax.plot(theta,g_ris.T)
    ax.set_title(" IRS Response gain vs angle of reflection in azimuth")
    ax.set_xlabel('theta_r')
    ax.set_ylabel('IRS response gain ')
    ax.set_xticks(np.arange(-90,100,10))
    ax.set_yticks(np.arange(-10,60,10))
    #ax.set_ylim(-40,20,10)
    ax.grid()
    plt.show()













