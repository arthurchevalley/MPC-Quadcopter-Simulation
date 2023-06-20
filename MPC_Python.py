#%%
import numpy as np
import math

import matplotlib.pyplot as plt

from scipy.linalg import solve_discrete_are
from scipy.integrate import solve_ivp

import polytope as pc
from cvxpy import Variable, Problem, Minimize, Parameter
import cvxpy as cp

import control

def plotState(full_command, full_state, simTime, dt=0.2, coords = None):
    
    x, y, z =[],[], []
    ref_x, ref_y, ref_z =[],[], []
    phi, theta, psi = [],[],[]
    vx,vy,vz = [],[],[]
    p,q,r = [],[],[]
    cmd1, cmd2, cmd3, cmd4 = [],[],[],[]
    id = 0
    for i in range(len(full_state)):

        cmd1.append(np.round(full_command[i][0],5))
        cmd2.append(np.round(full_command[i][1],5))
        cmd3.append(np.round(full_command[i][2],5))
        cmd4.append(np.round(full_command[i][3],5))

        x.append(np.round(full_state[i][9],5))
        y.append(np.round(full_state[i][10],5))
        z.append(np.round(full_state[i][11],5))
        if coords is not None:
            if i < len(coords) - 1:
                id += 1
            ref_y.append(coords[id][1])
            ref_x.append(coords[id][0])
            ref_z.append(coords[id][2])

        vx.append(np.round(full_state[i][6],5))
        vy.append(np.round(full_state[i][7],5))
        vz.append(np.round(full_state[i][8],5))

        phi.append(np.round(np.rad2deg(full_state[i][3]),5))
        theta.append(np.round(np.rad2deg(full_state[i][4]),5))
        psi.append(np.round(np.rad2deg(full_state[i][5]),5))

        p.append(np.round(np.rad2deg(full_state[i][0]),5))
        q.append(np.round(np.rad2deg(full_state[i][1]),5))
        r.append(np.round(np.rad2deg(full_state[i][2]),5))

    # position
    fig, ax = plt.subplots(3,1)
    fig.suptitle("Position in meters")
    ax[0].plot(x)
    ax[1].plot(y)
    ax[2].plot(z)
    if coords is not None:
        ax[1].plot(ref_y)
        ax[0].plot(ref_x)   
        ax[2].plot(ref_z)

    ax[2].set_xlabel('Simulation time')
    
    ax[0].set_ylabel('X')
    ax[1].set_ylabel('Y')
    ax[2].set_ylabel('Z')
    for i in range(2):
        ax[i].xaxis.set_tick_params(labelbottom=False)
        ax[i].set_xticks([])
    ax[2]
    # Velocities
    fig2, ax2 = plt.subplots(3,1)
    fig2.suptitle("Speed value in m/s")
    ax2[0].plot(vx)
    ax2[1].plot(vy)
    ax2[2].plot(vz)

    ax2[2].set_xlabel('Simulation time')
    
    ax2[0].set_ylabel('VX')
    ax2[1].set_ylabel('VY')
    ax2[2].set_ylabel('VZ')
    for i in range(2):
        ax2[i].xaxis.set_tick_params(labelbottom=False)
        ax2[i].set_xticks([])

    # Angles
    fig3, ax3 = plt.subplots(3,1)
    
    fig3.suptitle("Angles in degrees")
    ax3[0].plot(phi)
    ax3[1].plot(theta)
    ax3[2].plot(psi)

    ax3[2].set_xlabel('Simulation time')
    
    ax3[0].set_ylabel('Roll')
    ax3[1].set_ylabel('Pitch')
    ax3[2].set_ylabel('Yaw')
    for i in range(2):
        ax3[i].xaxis.set_tick_params(labelbottom=False)
        ax3[i].set_xticks([])


    # Rates
    fig4, ax4 = plt.subplots(3,1)
    fig4.suptitle("Rates in degrees/s")
    ax4[0].plot(p)
    ax4[1].plot(q)
    ax4[2].plot(r)

    ax4[2].set_xlabel('Simulation time')
    
    ax4[0].set_ylabel('Roll')
    ax4[1].set_ylabel('Pitch')
    ax4[2].set_ylabel('Yaw')
    for i in range(2):
        ax4[i].xaxis.set_tick_params(labelbottom=False)
        ax4[i].set_xticks([])

    # commadn
    fig5, ax5 = plt.subplots(4,1)
    fig5.suptitle("Commands")
    ax5[0].plot(cmd1)
    ax5[1].plot(cmd2)
    ax5[2].plot(cmd3)
    ax5[3].plot(cmd4)
    
    for i in range(3):
        ax5[i].xaxis.set_tick_params(labelbottom=False)
        ax5[i].set_xticks([])
        
    
    ax5[3].set_xlabel('Simulation time')
    
    ax6 = plt.figure().add_subplot(projection='3d')
    ax6.plot(x,y,z)
    plt.title("3D views ")
    
    fig7, ax7 = plt.subplots(1,3)
    ax7[0].set_xlabel('X Views')
    ax7[0].set_ylabel('Y Views')
    ax7[0].plot(x,y)
    if coords is not None:
        ax7[0].plot(ref_x, ref_y)
    
    ax7[1].set_xlabel('X Views')
    ax7[1].set_ylabel('Z Views')
    ax7[1].plot(x,z)
    if coords is not None:
        ax7[1].plot(ref_x, ref_z)
    
    ax7[2].set_xlabel('Y Views')
    ax7[2].set_ylabel('Z Views')
    ax7[2].plot(y, z)
    if coords is not None:
        ax7[2].plot(ref_y,ref_z)
    
    plt.show()

class MPC_Param:
    def __init__(self,
                 A=None,
                 B=None,
                 C=None,
                 Klqr=None,
                 Qf=None,
                 Qlqr=None,
                 Rlqr=None,
                 xs=None,
                 us=None,
                 R=None, 
                 M=None, 
                 m=None, 
                 F=None, 
                 f=None,
                 N = None):
        self.A = A
        self.B=B
        self.C=C
        self.Klqr=Klqr
        self.Qf=Qf
        self.Qlqr=Qlqr
        self.Rlqr=Rlqr
        self.xs=xs
        self.us=us
        self.u = None
        self.R= R
        self.M= M
        self.m= m
        self.F= F
        self.f=f
        self.Hu = None
        self.ku= None
        self.Hx = None
        self.kx= None
        self.N = N
        self.Lest = None
        self.x_selec = None

    def estimate(self, z_hat=None, poles=[0.4, 0.5, 0.9]):
        nx = self.A.shape[0]
        self.A_bar = []
        self.B_bar = []
        self.C_bar = []
        AB = np.concatenate([self.A, np.expand_dims(self.B, axis=1)], axis=1)
        C =  np.concatenate([np.zeros((1, nx)), np.ones((1,1))], axis=1)
        self.A_bar = np.concatenate([AB,C], axis=0)
        self.B_bar = np.append(self.B, np.zeros(1))
        self.C_bar = np.expand_dims(np.append(self.C, np.zeros(1)), axis=1)
        
        self.Lest = control.place(self.A_bar.T, self.C_bar, poles)[0]
      
        self.z_hat = z_hat
        
    def UpdateEstimate(self, state, input):
        aa = self.A_bar @ self.z_hat
        bb = self.B_bar * input
        cc = self.Lest * (self.C_bar.T @ self.z_hat - state[1])
        self.z_hat = aa + bb + cc
        d = 3
        
class Quad:
    def __init__(self,
                 x = np.zeros(12),
                 u = np.zeros(4),
                 mass=8,
                 thrustLimit = 1.5*np.array([[0, 0, 0, 0],
                                             [1, 1, 1, 1]]),
                 I = np.array([[10, 0, 0],
                               [0, 10, 0],
                               [0, 0, 15]]),
                 Kf = 28,
                 Km = 11,
                 nL = 0.2,
                 rad = 0.04,
                 bladeRad = 0.08,
                 samplingTime=0.1):
        
        # omega 0:3, 
        # theta 3:6 roll, pitch, yaw
        # speed 6:9, 
        # pos 9:12 x, y, z
        self.g = 9.81
        self.x = x
        self.u = u
        self.m = mass
        self.thrustLimit = thrustLimit
        self.I = I
        self.Kf = Kf
        self.Km = Km
        self.K = np.array([[Kf, Kf, Kf, Kf],
                            [0, Kf*nL, 0, -Kf*nL],
                            [-Kf*nL, 0, Kf*nL, 0],
                            [Km, -Km, Km, -Km]])
        
        self.T = self.K
        self.nL = nL
        self.rad = rad
        self.bladeRad = bladeRad
        self.Ts = samplingTime
        self.thrsutDir = np.array([[0,0,0,0],
                                   [0,0,0,0],
                                   [1,1,1,1]])
        
        self.L = np.array([[0.2,0,-0.2,0],
                           [0,0.2,0,-0.2],
                           [1,1,1,1]])
        self.mpc_x = None
        self.mpc_y = None
        self.mpc_z = None
        self.mpc_yaw = None
        
    def dynamics(self, t, x):
        self.u = np.clip(self.u, self.thrustLimit[0], self.thrustLimit[1])
        self.uTotal = self.K[0] @ self.u
        self.uMomemnts = self.K[1:] @ self.u
        self.R = np.array([[1, 0, 0],
                           [0, math.cos(x[3]), -math.sin(x[3])],
                           [0, math.sin(x[3]), math.cos(x[3])]])
        self.R = self.R @ np.array([[math.cos(x[4]), 0, math.sin(x[4])],
                                    [0, 1, 0],
                                    [-math.sin(x[4]), 0, math.cos(x[4])]])
        self.R = self.R @ np.array([[math.cos(x[5]), -math.sin(x[5]), 0],
                                    [math.sin(x[5]), math.cos(x[5]), 0],
                                    [0, 0, 1]])
        vDot = np.array([0,0,-self.m*9.81])+ self.uTotal*self.R@np.array([0,0,1])
        omegaDot = -np.cross(x[:3], self.I@x[:3]) + self.uMomemnts

        omegaDot /= np.array([self.I[0,0],self.I[1,1],self.I[2,2]])
        dx = np.append(omegaDot, [x[:3], vDot/self.m, x[6:9]])
        return dx
    
    def setStateSpace(self, xs=np.zeros(12), us=0.707*np.ones(4)):
        A = np.zeros((12,12))
        A[3,0] = 1
        A[4,1] = 1
        A[5,2] = 1
        A[7,3] = -self.g
        A[6,4] = self.g
        A[9,6] = 1
        A[10,7] = 1
        A[11,8] = 1

        B = np.zeros((12, 4))
        B[0] = [0, 0.1*self.Kf*self.nL, 0, -0.1*self.Kf*self.nL]
        B[1,:] = [-0.1*self.Kf*self.nL, 0, 0.1*self.Kf*self.nL, 0]
        B[2,:] = [self.Km/15, -self.Km/15, self.Km/15, -self.Km/15]
        B[8,:] = [self.Kf/self.m, self.Kf/self.m, self.Kf/self.m, self.Kf/self.m]

        Bp = B @ np.linalg.inv(self.T)
        
        C = np.eye(12)
        D = np.zeros((12,4))

        sys = control.StateSpace(A, Bp, C, D)
        
        sys_discrete = control.c2d(sys, self.Ts, method='zoh')
        
        self.A = np.array(sys_discrete.A)
        self.Bp = np.array(sys_discrete.B)
        self.C = np.array(sys_discrete.C)
        self.D = np.array(sys_discrete.D)
        
        
        
        self.xs = xs
        self.us = us

    def UpdateSubState(self, mpc_paramX=None, mpc_paramY=None, mpc_paramZ=None,mpc_paramYaw=None):
        if mpc_paramX is not None:
            mpc_paramX.x = self.x[mpc_paramX.x_selec].copy()
        if mpc_paramY is not None:
            mpc_paramY.x = self.x[mpc_paramY.x_selec].copy()
        if mpc_paramZ is not None:
            mpc_paramZ.x = self.x[mpc_paramZ.x_selec].copy()
        if mpc_paramYaw is not None:
            mpc_paramYaw.x = self.x[mpc_paramYaw.x_selec].copy()
        
    def decomposeSS(self, mpc_paramX, mpc_paramY, mpc_paramZ,mpc_paramYaw, non_linear=False):
        
        self.UpdateSubState(mpc_paramX, mpc_paramY, mpc_paramZ,mpc_paramYaw)
        
        
        mpc_paramX.A = self.A[[ 1, 4, 6, 9],:][:,[1,4,6,9]]
        mpc_paramY.A = self.A[[ 0, 3, 7, 10],:][:,[ 0, 3, 7, 10]]
        mpc_paramZ.A = self.A[[ 8, 11],:][:,[ 8, 11]]
        mpc_paramYaw.A = self.A[[ 2, 5],:][:,[ 2, 5]]

        mpc_paramX.B = self.Bp[[ 1, 4, 6, 9], 2]
        mpc_paramY.B = self.Bp[[ 0, 3, 7, 10], 1]
        mpc_paramZ.B = self.Bp[[ 8, 11], 0]
        mpc_paramYaw.B = self.Bp[[ 2, 5], 3]

        mpc_paramX.C = np.array([0, 0, 0, 1])
        mpc_paramY.C = np.array([0, 0, 0, 1])
        mpc_paramZ.C = np.array([0, 1])
        mpc_paramYaw.C = np.array([0, 1])

        if non_linear:
            mpc_paramX.Q = np.array([[ 25, 0, 0, 0],
                                     [ 0, 0, 0, 0],
                                     [ 0, 0, 5, 0],
                                     [ 0, 0, 0, 10]])
            mpc_paramY.Q = np.array([[ 10, 0, 0, 0],
                                     [ 0, 0, 0, 0],
                                     [ 0, 0, 10, 0],
                                     [ 0, 0, 0, 20]])
            mpc_paramZ.Q = np.array([[ 50, 0],
                                     [ 0, 400]])
            mpc_paramYaw.Q = np.array([[ 10, 0],
                                     [ 0, 40]])
            
            mpc_paramX.R = 20*np.ones(1)
            mpc_paramY.R = 20*np.ones(1)
            mpc_paramZ.R = 20*np.ones(1)
            mpc_paramYaw.R = 20*np.ones(1)
        else:
            mpc_paramX.Q = 10*np.eye(4)
            mpc_paramY.Q = 10*np.eye(4)
            mpc_paramZ.Q = 10*np.eye(2)
            mpc_paramYaw.Q = 10*np.eye(2)
            mpc_paramX.R = 8*np.ones(1)
            mpc_paramY.R = 8*np.ones(1)
            mpc_paramZ.R = 8*np.ones(1)
            mpc_paramYaw.R = 8*np.ones(1)

        mpc_paramX.Qlqr = np.eye(4)
        mpc_paramX.Rlqr = 8*np.ones(1)
        mpc_paramY.Qlqr = np.eye(4)
        mpc_paramY.Rlqr = 8*np.ones(1)
        mpc_paramZ.Qlqr = np.eye(2)
        mpc_paramZ.Rlqr = 8*np.ones(1)
        mpc_paramYaw.Qlqr = np.eye(2)
        mpc_paramYaw.Rlqr = 8*np.ones(1)

    def fuseCtrl(self, mpc_paramX, mpc_paramY, mpc_paramZ,mpc_paramYaw, ref):
        
        ref_x = ref[mpc_paramX.x_selec]
        opti_ctrX = self.ControlMPC(mpc_paramX, ref_x, x_0 = self.x[mpc_paramX.x_selec])
        
        ref_y = ref[mpc_paramY.x_selec]
        opti_ctrY = self.ControlMPC(mpc_paramY, ref_y, x_0 = self.x[mpc_paramY.x_selec])
        
        ref_z = ref[mpc_paramZ.x_selec]
        if mpc_paramZ.Lest is not None:
            opti_ctrZ = self.ControlMPC(mpc_paramZ, ref_z, d_est=mpc_paramZ.z_hat[-1], x_0=mpc_paramZ.z_hat[:2])
            mpc_paramZ.u = opti_ctrZ
        else:
            opti_ctrZ = self.ControlMPC(mpc_paramZ, ref_z, x_0 = self.x[mpc_paramZ.x_selec])
            
        ref_yaw = ref[mpc_paramYaw.x_selec]
        opti_ctrYaw = self.ControlMPC(mpc_paramYaw, ref_yaw, x_0 = self.x[mpc_paramYaw.x_selec])

        Tinv = np.linalg.inv(self.T)
        u_v = np.array([opti_ctrZ, opti_ctrY, opti_ctrX, opti_ctrYaw]).flatten()
        u = Tinv @ u_v
        return u, u_v
    
    def GetLQRGain(self, A, B, Q, R):

        # You can play with this number
        N = 20
    
        # Create a list of N + 1 elements
        P = [None] * (N + 1)
        
        Qf = Q
    
        # LQR via Dynamic Programming
        P[N] = Qf
    
        # Calculate the optimal feedback gain K
        B = np.expand_dims(B, axis=0).T

        P = solve_discrete_are(A, B, Q, R)
        K = (-np.linalg.inv(B.T @ P @ B + R) @ (B.T @ P @ A))[0]

        #K = -np.linalg.pinv(R + (B).T @ P[0] @ B) @ B.T @ P[0] @ A
        return K, Qf

    def GetSubSystemLQR(self):
        self.KlqrX, self.Qfx = self.GetLQRGain(self.Ax, self.Bx, self.Qlqr_X, self.Rlqr_X)
        self.KlqrY, self.Qfy  = self.GetLQRGain(self.Ay, self.By, self.Qlqr_Y, self.Rlqr_Y)
        self.KlqrZ, self.Qfz  = self.GetLQRGain(self.Az, self.Bz, self.Qlqr_Z, self.Rlqr_Z)
        self.KlqrYaw, self.Qfyaw  = self.GetLQRGain(self.Ayaw, self.Byaw, self.Qlqr_Yaw, self.Rlqr_Yaw)

    def GetOptiInputTracking(self, mpc_param, xs_val, us_val, d_est=None, x_0=None):

        if len(mpc_param.B.shape) > 1:
            n, m = mpc_param.B.shape
        else:
            n = mpc_param.B.shape[0]
            m=1
        # Define variables and parameters
        N = mpc_param.N  # Number of time steps
        x = Variable((n, N))  # State variables
        u = Variable((m, N))  # Input variables
        xs = Parameter(n)  # Reference state
        us = Parameter(m)  # Reference input
        x_init = Parameter(n)
        xs.value = xs_val
        us.value = us_val
        
        # Define constraints and objective

        constr = [x[:, 0] == x_init]
        if d_est is not None:
            constr += [x[:, 1] == mpc_param.A @ (x[:, 0]) + mpc_param.B *( u[:, 0]) + mpc_param.B * d_est]
        else:
            constr += [x[:, 1] == mpc_param.A @ (x[:, 0]) + mpc_param.B *( u[:, 0])]
        obj = cp.quad_form(u[:, 0], np.expand_dims(mpc_param.R, axis=1)) 
        if mpc_param.M.shape[0]:
                constr += [cp.multiply(mpc_param.M, u[:, 0]) <= m]
                
        for i in range(1,N-1):
            
            if d_est is not None:
                constr += [x[:, i+1] == mpc_param.A @ (x[:, i]) + mpc_param.B *( u[:, i]) + mpc_param.B * d_est]
            else:
                constr += [x[:, i + 1] == mpc_param.A @ (x[:, i]) + mpc_param.B * (u[:, i])]
            if mpc_param.f is not None: 
                AA = mpc_param.F @ x[:, i]
                constr += [AA <= (mpc_param.f[0])]
            
            if mpc_param.M is not None:
                constr += [cp.multiply(mpc_param.M, u[:, i]) <= m] # * problem
                
            obj += cp.quad_form((x[:, i] - xs), mpc_param.Q) + cp.quad_form((u[:, i] - us), np.expand_dims(mpc_param.R, axis=1)) 
        # Define the optimization problem
        obj += cp.quad_form((x[:, N-1] - xs), mpc_param.Q) + cp.quad_form((u[:, N-1] - us), np.expand_dims(mpc_param.R, axis=1)) 

        problem = Problem(Minimize(obj), constr)

        # Solve the optimization problem
        if x_0 is None:
            x_init.value = mpc_param.x
        else:
            x_init.value = x_0

        problem.solve(solver=cp.OSQP, warm_start=True, max_iter=100000)
        new_u = u[:,0].value
        return new_u

    def GetOptiInput(self, mpc_param):
        
        # Compute LQR controller for unconstrained system
        # Compute maximal invariant set
        if mpc_param.F.shape[0] and mpc_param.M.shape[0]:
            FM = np.vstack((mpc_param.F, mpc_param.M * mpc_param.Klqr))
            fm = np.vstack((mpc_param.f, mpc_param.m))
        elif mpc_param.F.shape[0]: 
            FM = mpc_param.F
            fm = mpc_param.f
        elif mpc_param.M.shape[0]:
            FM = mpc_param.M * mpc_param.Klqr
            fm = mpc_param.m
        else:
            raise Exception('Bad state or input contraints') 
        Xf = pc.Polytope(FM, fm)
        BB = np.expand_dims(mpc_param.B, axis=0).T
        Acl = mpc_param.A + (BB) * mpc_param.Klqr

        while True:

            prevXf = Xf
            T, t = Xf.A, Xf.b
            preXf = pc.Polytope(T @ Acl, t)
            Xf = Xf.intersect(preXf)
            if prevXf == Xf:
                break

        # Extract the vertices of Xf
        Ff, ff = Xf.A, Xf.b
        if len(mpc_param.B.shape) > 1:
            n, m = mpc_param.B.shape
        else:
            n = mpc_param.B.shape[0]
            m=1
        # Define variables and parameters
        N = mpc_param.N  # Number of time steps
        x = Variable((n, N))  # State variables
        u = Variable((m, N))  # Input variables
        
        x_init = Parameter(n)
        
        # Define constraints and objective

        constr = [x[:, 0] == x_init]
        constr += [x[:, 1] == mpc_param.A @ (x[:, 0]) + mpc_param.B *( u[:, 0])]
        obj = cp.quad_form(u[:, 0], np.expand_dims(mpc_param.R, axis=1)) 
        if mpc_param.M.shape[0]:
                constr += [cp.multiply(mpc_param.M, u[:, 0]) <= m]
                
        for i in range(1,N-1):
            constr += [x[:, i + 1] == mpc_param.A @ (x[:, i]) + mpc_param.B * (u[:, i])]
            
            if mpc_param.F.shape[0]: 
                AA = mpc_param.F @ x[:, i]
                constr += [AA <= (mpc_param.f[0])]
            
            if mpc_param.M.shape[0]:
                constr += [cp.multiply(mpc_param.M, u[:, i]) <= m] # * problem
                
            obj += cp.quad_form((x[:, i]), mpc_param.Q) + cp.quad_form((u[:, i]), np.expand_dims(mpc_param.R, axis=1)) 
        # Define the optimization problem
        
        constr += [Ff @ (x[:, N-1]) <= ff]
        obj += cp.quad_form((x[:, (N-1)]), mpc_param.Qf)

        problem = Problem(Minimize(obj), constr)

        # Solve the optimization problem
        x_init.value = mpc_param.x

        problem.solve(solver=cp.OSQP, warm_start=True, max_iter=100000)
        new_u = u[:,0].value
        return new_u

    
    def setup_steady_state_target(self,mpc_ctrl, ref_val=None, d_est=None):
        # Steady-state targets

        if len(mpc_ctrl.B.shape) > 1:
            n, m = mpc_ctrl.B.shape
        else:
            n = mpc_ctrl.B.shape[0]
            m=1
        xs = Variable(n)
        us = Variable(m)
        
        # Reference position (Ignore this before Todo 3.2)
        ref = Parameter(n)
        ref.value = ref_val
        if d_est is not None:
            d_est_p = Parameter(m)
            d_est_p.value = np.array([d_est])
            
        # Constraints
        
        
        Rs = np.eye(us.shape[0])
        if d_est is not None:
            BB = mpc_ctrl.B * d_est_p
            constr = [xs == mpc_ctrl.A @ xs + mpc_ctrl.B * us + BB]
            if ref_val is not None:
                constr += [mpc_ctrl.C @ ref == (mpc_ctrl.C @ xs)]
        
        else:
            constr = [xs == mpc_ctrl.A @ xs + mpc_ctrl.B * us]
            if ref_val is not None:
                constr += [mpc_ctrl.C @ ref == (mpc_ctrl.C @ xs)]
            
        if mpc_ctrl.f is not None: 
            constr += [mpc_ctrl.F @ xs <= mpc_ctrl.f.flatten()]
        if mpc_ctrl.m is not None:
            constr += [cp.multiply(mpc_ctrl.M, us) <= mpc_ctrl.m]
        obj = cp.quad_form(us, Rs)

        # Define the optimization problem

        # Compute the steady-state target
        problem = Problem(Minimize(obj), constr)

        # Solve the optimization problem
        problem.solve()

        return xs.value.flatten(), us.value.flatten()

    def ControlMPC(self, mpc_param, ref_val =None, d_est = None, x_0 = None):
        if ref_val is None:
            mpc_param.Klqr, mpc_param.Qf = self.GetLQRGain(mpc_param.A, mpc_param.B, mpc_param.Qlqr, mpc_param.Rlqr)
            opti_ctr = self.GetOptiInput(mpc_param)
        else:
            xs_val, us_val = self.setup_steady_state_target(mpc_param, ref_val=ref_val, d_est=d_est)
            opti_ctr = self.GetOptiInputTracking(mpc_param, xs_val, us_val, d_est=d_est, x_0=x_0)
        return opti_ctr
        
    def UpdateSS(self, mpc_param, u):
        nx = mpc_param.A @ (mpc_param.x) + mpc_param.B * (u)
        return nx
       
    def MPC_ref(self, coords, id, waypt_thresh = .05):
        if np.linalg.norm(self.x[11] - coords[id,2]) <= waypt_thresh/2 and id < (coords.shape[0]-1):
            if np.linalg.norm(self.x[9:11] - coords[id,:2]) <= waypt_thresh and id < (coords.shape[0]-1):
                id += 1

        xr = np.zeros(12)
        xr[9:] = coords[id,:]
        return xr, id

def MPC_Regulator():
    drone = Quad()

    #drone.dynamics()
    simulationTime = 15
    drone.Ts = .2

    
    xs = np.zeros(12)
    u_s = 0.707*np.ones(4)
    x0 = np.zeros(12)
    x0[9:] = [2,2,2]
    x0[5] = math.pi/4
    drone.x = -x0.copy()
    drone.mpc_x = MPC_Param(M = np.array([[1],[-1]]),
                            m = np.array([[0.3], [0.3]]),
                            F = np.array([[0, 1, 0, 0],[0, -1, 0, 0]]), 
                            f = np.array([[0.035],[0.035]]),
                            N = 20)
    drone.mpc_x.x_selec = [1,4,6,9]
    drone.mpc_x.xs = xs[drone.mpc_x.x_selec]
    drone.mpc_x.us = np.array([u_s[2]])
    
    drone.mpc_y = MPC_Param(M = np.array([[1],[-1]]),
                            m = np.array([[0.3], [0.3]]),
                            F = np.array([[0, 1, 0, 0],[0, -1, 0, 0]]), 
                            f = np.array([[0.035],[0.035]]),
                            N = 20)
    drone.mpc_y.x_selec = [0, 3, 7, 10]
    drone.mpc_y.xs = xs[drone.mpc_y.x_selec]
    drone.mpc_y.us = np.array([u_s[1]])

    drone.mpc_z = MPC_Param(M = np.array([[1],[-1]]),
                            m = np.array([[0.3], [0.2]]),
                            N = 20)
    drone.mpc_z.x_selec = [8, 11]
    drone.mpc_z.xs = xs[drone.mpc_z.x_selec]
    drone.mpc_z.us = np.array([u_s[0]])

    drone.mpc_yaw = MPC_Param(M = np.array([[1],[-1]]),
                              m = np.array([[0.2], [0.2]]),
                            N = 20)
    drone.mpc_yaw.x_selec = [2, 5]
    drone.mpc_yaw.xs = xs[drone.mpc_yaw.x_selec]
    drone.mpc_yaw.us = np.array([u_s[3]])

    
    drone.setStateSpace()
    drone.decomposeSS(drone.mpc_x, drone.mpc_y, drone.mpc_z, drone.mpc_yaw)

    state_hist = [-x0]
    input_hist = [u_s]
    for t in range(0, int(simulationTime/drone.Ts)):
        # X subsystem       
        opti_ctrX = drone.ControlMPC(drone.mpc_x)
        nx = drone.UpdateSS(drone.mpc_x, opti_ctrX)
        drone.mpc_x.x = nx

        # Y subsystem
        opti_ctrY = drone.ControlMPC(drone.mpc_y)
        ny = drone.UpdateSS(drone.mpc_y, opti_ctrY)
        drone.mpc_y.x = ny
        
        # Z subsystem
        opti_ctrZ = drone.ControlMPC(drone.mpc_z)
        nz = drone.UpdateSS(drone.mpc_z, opti_ctrZ)
        drone.mpc_z.x = nz
        
        # Yaw subsystem
        opti_ctrYaw = drone.ControlMPC(drone.mpc_yaw)
        nyaw = drone.UpdateSS(drone.mpc_yaw, opti_ctrYaw)
        drone.mpc_yaw.x = nyaw

        drone.x[drone.mpc_x.x_selec] = nx
        drone.x[drone.mpc_y.x_selec] = ny
        drone.x[drone.mpc_z.x_selec] = nz
        drone.x[drone.mpc_yaw.x_selec] = nyaw 

        u_sys = np.array([opti_ctrX,opti_ctrY,opti_ctrZ,opti_ctrYaw])

        
        input_hist.append(u_sys)
        state_hist.append(drone.x.copy())
        #print(f'{t} second have been simulated')
        if not((t+1)*drone.Ts % 1):
            print(f'{(t+1)*drone.Ts} second have been simulated')
    plotState(input_hist, state_hist, simulationTime, drone.Ts)    
    print("DOne")

def MPC_Tracking(ref):
    drone = Quad()

    #drone.dynamics()
    simulationTime = 80
    drone.Ts = .2

    
    xs = np.zeros(12)
    u_s = 0.707*np.ones(4)
    x0 = np.zeros(12)
    drone.x = x0.copy()
    drone.mpc_x = MPC_Param(M = np.array([[1],[-1]]),
                            m = np.array([[0.3], [0.3]]),
                            F = np.array([[0, 1, 0, 0],[0, -1, 0, 0]]), 
                            f = np.array([[0.035],[0.035]]),
                            N = 20)
    drone.mpc_x.x_selec = [1,4,6,9]
    drone.mpc_x.xs = xs[drone.mpc_x.x_selec]
    drone.mpc_x.us = np.array([u_s[2]])
    
    drone.mpc_y = MPC_Param(M = np.array([[1],[-1]]),
                            m = np.array([[0.3], [0.3]]),
                            F = np.array([[0, 1, 0, 0],[0, -1, 0, 0]]), 
                            f = np.array([[0.035],[0.035]]),
                            N = 20)
    drone.mpc_y.x_selec = [0, 3, 7, 10]
    drone.mpc_y.xs = xs[drone.mpc_y.x_selec]
    drone.mpc_y.us = np.array([u_s[1]])

    drone.mpc_z = MPC_Param(M = np.array([[1],[-1]]),
                            m = np.array([[0.3], [0.2]]),
                            F = np.array([[1,0],[-1,0]]), 
                            f = np.array([[15],[15]]),
                            N = 20)
    drone.mpc_z.x_selec = [8, 11]
    drone.mpc_z.xs = xs[drone.mpc_z.x_selec]
    drone.mpc_z.us = np.array([u_s[0]])

    drone.mpc_yaw = MPC_Param(M = np.array([[1],[-1]]),
                              m = np.array([[0.2], [0.2]]),
                              F = np.array([[1,0],[-1,0]]), 
                              f = np.array([[15],[15]]),
                              N = 20)
    drone.mpc_yaw.x_selec = [2, 5]
    drone.mpc_yaw.xs = xs[drone.mpc_yaw.x_selec]
    drone.mpc_yaw.us = np.array([u_s[3]])

    
    drone.setStateSpace()
    drone.decomposeSS(drone.mpc_x, drone.mpc_y, drone.mpc_z, drone.mpc_yaw)

    coords = np.array([[0, 0, 1],[0, 2, 1],[1, 1, 1],[2, 2, 1],[2, 0, 1],
                           [3, 0, 1.5],[3, 2, 1.5],[4, 2, 1.5],[4, 1, 1.5],[3, 1, 1.5],
                           [3, 0, 2],[7, 0, 2],[5, 0, 2],[5, 2, 2],[7, 2, 2]])
    
    state_hist = [x0]
    input_hist = [u_s]
    wypt_id = 0
    for t in np.arange(0, simulationTime, drone.Ts):
        # X subsystem       
        ref, wypt_id = drone.MPC_ref(coords, wypt_id)
        ref_x = ref[drone.mpc_x.x_selec]
        opti_ctrX = drone.ControlMPC(drone.mpc_x, ref_x)
        nx = drone.UpdateSS(drone.mpc_x, opti_ctrX)
        drone.mpc_x.x = nx

        # Y subsystem
        ref_y = ref[drone.mpc_y.x_selec]
        opti_ctrY = drone.ControlMPC(drone.mpc_y, ref_y)
        ny = drone.UpdateSS(drone.mpc_y, opti_ctrY)
        drone.mpc_y.x = ny
        
        # Z subsystem
        ref_z = ref[drone.mpc_z.x_selec]
        opti_ctrZ = drone.ControlMPC(drone.mpc_z, ref_z)
        nz = drone.UpdateSS(drone.mpc_z, opti_ctrZ)
        drone.mpc_z.x = nz
        
        # Yaw subsystem
        ref_yaw = ref[drone.mpc_yaw.x_selec]
        opti_ctrYaw = drone.ControlMPC(drone.mpc_yaw, ref_yaw)
        nyaw = drone.UpdateSS(drone.mpc_yaw, opti_ctrYaw)
        drone.mpc_yaw.x = nyaw

        drone.x[drone.mpc_x.x_selec] = nx
        drone.x[drone.mpc_y.x_selec] = ny
        drone.x[drone.mpc_z.x_selec] = nz
        drone.x[drone.mpc_yaw.x_selec] = nyaw 

        u_sys = np.array([opti_ctrX,opti_ctrY,opti_ctrZ,opti_ctrYaw])

        input_hist.append(u_sys)
        state_hist.append(drone.x.copy())

        if not(t % 1):
            print(f'{t} second have been simulated')
    plotState(input_hist, state_hist, simulationTime, drone.Ts, coords)    
    print("DOne")

def MPC_Tracking_Estimate(ref):
    drone = Quad()

    #drone.dynamics()
    simulationTime = 80
    drone.Ts = .2

    
    xs = np.zeros(12)
    u_s = 0.707*np.ones(4)
    x0 = np.zeros(12)
    drone.x = x0.copy()
    drone.mpc_x = MPC_Param(M = np.array([[1],[-1]]),
                            m = np.array([[0.3], [0.3]]),
                            F = np.array([[0, 1, 0, 0],[0, -1, 0, 0]]), 
                            f = np.array([[0.035],[0.035]]),
                            N = 20)
    drone.mpc_x.x_selec = [1,4,6,9]
    drone.mpc_x.xs = xs[drone.mpc_x.x_selec]
    drone.mpc_x.us = np.array([u_s[2]])
    
    drone.mpc_y = MPC_Param(M = np.array([[1],[-1]]),
                            m = np.array([[0.3], [0.3]]),
                            F = np.array([[0, 1, 0, 0],[0, -1, 0, 0]]), 
                            f = np.array([[0.035],[0.035]]),
                            N = 20)
    drone.mpc_y.x_selec = [0, 3, 7, 10]
    drone.mpc_y.xs = xs[drone.mpc_y.x_selec]
    drone.mpc_y.us = np.array([u_s[1]])

    drone.mpc_z = MPC_Param(M = np.array([[1],[-1]]),
                            m = np.array([[0.3], [0.2]]),
                            F = np.array([[1,0],[-1,0]]), 
                            f = np.array([[15],[15]]),
                            N = 20)
    drone.mpc_z.x_selec = [8, 11]
    drone.mpc_z.xs = xs[drone.mpc_z.x_selec]
    drone.mpc_z.us = np.array([u_s[0]])

    drone.mpc_yaw = MPC_Param(M = np.array([[1],[-1]]),
                              m = np.array([[0.2], [0.2]]),
                              F = np.array([[1,0],[-1,0]]), 
                              f = np.array([[15],[15]]),
                            N = 20)
    drone.mpc_yaw.x_selec = [2, 5]
    drone.mpc_yaw.xs = xs[drone.mpc_yaw.x_selec]
    drone.mpc_yaw.us = np.array([u_s[3]])

    
    drone.setStateSpace()
    drone.decomposeSS(drone.mpc_x, drone.mpc_y, drone.mpc_z, drone.mpc_yaw)

    coords = np.array([[0, 0, 1],[0, 2, 1],[1, 1, 1],[2, 2, 1],[2, 0, 1],
                           [3, 0, 1.5],[3, 2, 1.5],[4, 2, 1.5],[4, 1, 1.5],[3, 1, 1.5],
                           [3, 0, 2],[7, 0, 2],[5, 0, 2],[5, 2, 2],[7, 2, 2]])/2
    
    d_est = np.array([-0.1])
    
    drone.mpc_z.estimate(z_hat=np.array([0,0,d_est[0]]))
    state_hist = [x0]
    input_hist = [u_s]
    
    wypt_id = 0
    for t in np.arange(0, simulationTime, drone.Ts):
        # X subsystem       
        ref, wypt_id = drone.MPC_ref(coords, wypt_id)
        ref_x = ref[drone.mpc_x.x_selec]
        opti_ctrX = drone.ControlMPC(drone.mpc_x, ref_x)
        nx = drone.UpdateSS(drone.mpc_x, opti_ctrX)
        drone.mpc_x.x = nx

        # Y subsystem
        ref_y = ref[drone.mpc_y.x_selec]
        opti_ctrY = drone.ControlMPC(drone.mpc_y, ref_y)
        ny = drone.UpdateSS(drone.mpc_y, opti_ctrY)
        drone.mpc_y.x = ny
        
        # Z subsystem
        ref_z = ref[drone.mpc_z.x_selec]
        opti_ctrZ = drone.ControlMPC(drone.mpc_z, ref_z, d_est=drone.mpc_z.z_hat[2], x_0 = drone.mpc_z.z_hat[:2])
        drone.mpc_z.UpdateEstimate(drone.mpc_z.x, opti_ctrZ)

        drone.mpc_z.x = drone.mpc_z.z_hat[:2].copy()
        
        
        
        # Yaw subsystem
        ref_yaw = ref[drone.mpc_yaw.x_selec]
        opti_ctrYaw = drone.ControlMPC(drone.mpc_yaw, ref_yaw)
        nyaw = drone.UpdateSS(drone.mpc_yaw, opti_ctrYaw)
        drone.mpc_yaw.x = nyaw

        drone.x[drone.mpc_x.x_selec] = nx
        drone.x[drone.mpc_y.x_selec] = ny
        drone.x[drone.mpc_z.x_selec] = drone.mpc_z.z_hat[:2].copy()
        drone.x[drone.mpc_yaw.x_selec] = nyaw 

        u_sys = np.array([opti_ctrX,opti_ctrY,opti_ctrZ,opti_ctrYaw])

        input_hist.append(u_sys)
        state_hist.append(drone.x.copy())

        if not(t % 1):
            print(f'{t} second have been simulated')
    plotState(input_hist, state_hist, simulationTime, drone.Ts, coords)    
    print("DOne")

def NonLinear_MPC_Tracking():
    drone = Quad()

    simulationTime = 20
    drone.Ts = .2

    
    xs = np.zeros(12)
    u_s = 0.707*np.ones(4)
    x0 = np.zeros(12)
    drone.x = x0.copy()
    drone.mpc_x = MPC_Param(M = np.array([[1],[-1]]),
                            m = np.array([[0.3], [0.3]]),
                            F = np.array([[0, 1, 0, 0],[0, -1, 0, 0]]), 
                            f = np.array([[0.035],[0.035]]),
                            N = 30)
    drone.mpc_x.x_selec = [1,4,6,9]
    drone.mpc_x.xs = xs[drone.mpc_x.x_selec]
    drone.mpc_x.us = np.array([u_s[2]])
    
    drone.mpc_y = MPC_Param(M = np.array([[1],[-1]]),
                            m = np.array([[0.3], [0.3]]),
                            F = np.array([[0, 1, 0, 0],[0, -1, 0, 0]]), 
                            f = np.array([[0.035],[0.035]]),
                            N = 30)
    drone.mpc_y.x_selec = [0, 3, 7, 10]
    drone.mpc_y.xs = xs[drone.mpc_y.x_selec]
    drone.mpc_y.us = np.array([u_s[1]])

    drone.mpc_z = MPC_Param(M = np.array([[1],[-1]]),
                            m = np.array([[0.3], [0.2]]),
                            N = 30)
    drone.mpc_z.x_selec = [8, 11]
    drone.mpc_z.xs = xs[drone.mpc_z.x_selec]
    drone.mpc_z.us = np.array([u_s[0]])

    drone.mpc_yaw = MPC_Param(M = np.array([[1],[-1]]),
                              m = np.array([[0.2], [0.2]]),
                              N = 30)
    drone.mpc_yaw.x_selec = [2, 5]
    drone.mpc_yaw.xs = xs[drone.mpc_yaw.x_selec]
    drone.mpc_yaw.us = np.array([u_s[3]])

    
    drone.setStateSpace()
    drone.decomposeSS(drone.mpc_x, drone.mpc_y, drone.mpc_z, drone.mpc_yaw, non_linear=True)

    state_hist = [x0]
    input_hist = [u_s]
    coords = np.array([[0,0,0],[0, 0, 1],[0, 2, 1],[1, 1, 1],[2, 2, 1],[2, 0, 1],
                           [3, 0, 1.5],[3, 2, 1.5],[4, 2, 1.5],[4, 1, 1.5],[3, 1, 1.5],
                           [3, 0, 2],[7, 0, 2],[5, 0, 2],[5, 2, 2],[7, 2, 2]])/2
    
    ref_hist = []
    wypt_id = 0
    for t in np.arange(0, simulationTime, drone.Ts):
        tspan = (t, t+drone.Ts)
        ref, wypt_id = drone.MPC_ref(coords, wypt_id)

        u_input, u_v = drone.fuseCtrl(drone.mpc_x, drone.mpc_y, drone.mpc_z, drone.mpc_yaw, ref)
        drone.u =  u_input + u_s
        x_ode = solve_ivp(drone.dynamics, t_span=tspan, y0=drone.x)
        drone.x = x_ode.y[:,-1]
        #drone.UpdateSubState(drone.mpc_x, drone.mpc_y, drone.mpc_z,drone.mpc_yaw)
        
        ref_hist.append(coords[wypt_id, :].copy())
        input_hist.append(drone.u.copy())
        state_hist.append(drone.x.copy())

        if not t % 1:
            print(f'{t} second have been simulated')
    plotState(input_hist, state_hist, simulationTime,drone.Ts ,ref_hist)    
    print("DOne")



if __name__ == '__main__':

    mpc_type = ['regulator', 'tracking','MPC_Tracking_Estimate','non-linear tracking']
    mpc_choice = mpc_type[3]
    
    match mpc_choice:
        case 'regulator':
            MPC_Regulator()
        case 'tracking':
            ref = np.array([0, 0, 0, 0, 0, math.pi/4,0,0,0,-2,-2,-2])
            MPC_Tracking(ref)
        case 'MPC_Tracking_Estimate':
            ref = np.array([0, 0, 0, 0, 0, math.pi/4,0,0,0,-2,-2,-2])
            MPC_Tracking_Estimate(ref)
        case 'non-linear tracking':
            NonLinear_MPC_Tracking()
