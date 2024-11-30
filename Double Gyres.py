"""
@ A program for modeling double gyres
@ author ZimoJupiter
@ w.zimo@outlook.com
@ date 30 Nov 2024
@ license MIT License
"""
import numpy as np
import pandas as pd
from numpy import pi, exp, sqrt, sin, cos
import scipy
import matplotlib.pyplot as plt
import imageio
import glob
plt.rcParams['font.weight'] = 'normal'
plt.rcParams["figure.figsize"] = (3.2, 3.2*3/4)
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'

def meshing():
    Lx = 2
    Ly = 1
    dx = 0.01
    dy = 0.01
    x = np.arange(0, Lx+dx, dx)
    y = np.arange(0, Ly+dx, dy)
    mesh_X, mesh_Y = np.meshgrid(x, y)
    # plt.figure()
    # plt.scatter(X, Y, s=0.01)
    # plt.show()
    return mesh_X, mesh_Y

def DoubleGyres(x, y, t):
    A = 0.1
    epsilon = 0.25
    omega = 2*pi/10

    a = epsilon*sin(omega*t)
    b = 1 - 2*epsilon*sin(omega*t)
    f = a*x**2 + b*x
    dfdx = 2*a*x + b

    u = -pi*A*sin(pi*f)*cos(pi*y)
    v = pi*A*cos(pi*f)*sin(pi*y)*dfdx
    return np.array([u, v])

def compute_Q_criterion_2D(u, v, x, y, t):
    x_1d = x[0, :]
    y_1d = y[:, 0]

    dudx = np.gradient(u, x[0, :], axis=1)
    dudy = np.gradient(u, y[:, 0], axis=0)
    dvdx = np.gradient(v, x[0, :], axis=1)
    dvdy = np.gradient(v, y[:, 0], axis=0)
    
    S2 = (dudx**2 + dvdy**2 + (dudy + dvdx)**2) / 2
    Omega2 = ((dudy - dvdx)**2) / 2
    Q = 0.5 * (Omega2 - S2)

    plt.figure()
    plt.contourf(x, y, Q, levels=np.arange(-2,2+0.1,0.1), cmap='jet')
    plt.colorbar(label='Q Criterion')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig(f'Figures/Q_criterion_{t}.png')

def Trajectories():
    dt = 0.02
    t = np.arange(0, 10+2*dt, dt)
    X = np.zeros((t.shape[0], mesh_X.shape[0], mesh_X.shape[1]))
    Y = np.zeros((t.shape[0], mesh_Y.shape[0], mesh_Y.shape[1]))
    X[0] = mesh_X
    Y[0] = mesh_Y
    Tr_A2 = np.zeros((t.shape[0]))
    a2 = np.zeros((t.shape[0]))
    d0_index = np.arange(5, 201, 5)
    d0_index = np.array([1, 5, 10, 20])
    d0 = 0.01*d0_index
    Tr_R2 = np.zeros((t.shape[0], d0_index.shape[0]))
    r2 = np.zeros((t.shape[0], d0_index.shape[0]))

    U = np.zeros((t.shape[0], mesh_X.shape[0], mesh_X.shape[1]))
    V = np.zeros((t.shape[0], mesh_Y.shape[0], mesh_Y.shape[1]))

    k1 = np.zeros((2, mesh_X.shape[0], mesh_Y.shape[1]))
    k2 = np.zeros((2, mesh_X.shape[0], mesh_Y.shape[1]))
    k3 = np.zeros((2, mesh_X.shape[0], mesh_Y.shape[1]))
    k4 = np.zeros((2, mesh_X.shape[0], mesh_Y.shape[1]))

    for t_i in range(t.shape[0]-1):
        U[t_i], V[t_i] = DoubleGyres(mesh_X, mesh_Y, t[t_i])

        k1 = dt*DoubleGyres(X[t_i], Y[t_i], t[t_i])
        k2 = dt*DoubleGyres(X[t_i]+k1[0]/2, Y[t_i]+k1[1]/2, t[t_i]+dt/2)
        k3 = dt*DoubleGyres(X[t_i]+k2[0]/2, Y[t_i]+k2[1]/2, t[t_i]+dt/2)
        k4 = dt*DoubleGyres(X[t_i]+k3[0], Y[t_i]+k3[1], t[t_i]+dt)

        X[t_i+1] = X[t_i] + (k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
        Y[t_i+1] = Y[t_i] + (k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
        
        plt.figure()
        plt.scatter(X[t_i,::5,::5], Y[t_i,::5,::5], facecolor='r', edgecolor='None', s=10, marker='.')
        plt.quiver(mesh_X[::10,::10], mesh_Y[::10,::10], U[t_i,::10,::10], V[t_i,::10,::10], \
                color='b', headaxislength=5, minshaft=5, scale=3)
        plt.text(0.75, 0.85, f't = {t[t_i]:.2f} s', transform=plt.gca().transAxes, \
            bbox=dict(facecolor='white', edgecolor='white', alpha=0.8))
        plt.xlabel('X', labelpad=0)
        plt.ylabel('Y', labelpad=0)
        plt.tight_layout()
        plt.savefig(f'Figures/DG/DoubleGyres_velocity_{t[t_i]:.2f}.png')

        Tr_A2[t_i] = np.sum((X[t_i]-X[0])**2+(Y[t_i]-Y[0])**2)/(2*mesh_X.shape[0]*mesh_X.shape[1])
        a2[t_i] = (np.sum((X[t_i]-X[0])*U[t_i]+(Y[t_i]-Y[0])*V[t_i])/(mesh_X.shape[0]*mesh_X.shape[1]))
        
        for d0_i in range(d0_index.shape[0]):
            X_diff_x = X[t_i, :-d0_index[d0_i], :] - X[t_i, d0_index[d0_i]:, :]
            Y_diff_x = Y[t_i, :-d0_index[d0_i], :] - Y[t_i, d0_index[d0_i]:, :]
            X_diff_y = X[t_i, :, :-d0_index[d0_i]] - X[t_i, :, d0_index[d0_i]:]
            Y_diff_y = Y[t_i, :, :-d0_index[d0_i]] - Y[t_i, :, d0_index[d0_i]:]

            U_diff_x = U[t_i, :-d0_index[d0_i], :] - U[t_i, d0_index[d0_i]:, :]
            V_diff_x = V[t_i, :-d0_index[d0_i], :] - V[t_i, d0_index[d0_i]:, :]
            U_diff_y = U[t_i, :, :-d0_index[d0_i]] - U[t_i, :, d0_index[d0_i]:]
            V_diff_y = V[t_i, :, :-d0_index[d0_i]] - V[t_i, :, d0_index[d0_i]:]

            X_sum = np.sum(X_diff_x**2) + np.sum(X_diff_y**2)
            Y_sum = np.sum(Y_diff_x**2) + np.sum(Y_diff_y**2)
            XU_sum = np.sum(X_diff_x*U_diff_x) + np.sum(X_diff_y*U_diff_y)
            YV_sum = np.sum(Y_diff_x*V_diff_x) + np.sum(Y_diff_y*V_diff_y)

            Tr_R2[t_i, d0_i] = (X_sum + Y_sum) \
                                /(2*(X_diff_x.shape[0])*(X_diff_y.shape[1]))
            r2[t_i, d0_i] = (XU_sum + YV_sum) \
                                /((X_diff_x.shape[0])*(X_diff_y.shape[1]))
        
        if t_i == 0 or t_i == 125 or t_i == 250 or t_i == 375 or t_i == 500:
            compute_Q_criterion_2D(U[t_i], V[t_i], mesh_X, mesh_Y, t[t_i])

    file_list = sorted(glob.glob('Figures/DG/DoubleGyres_velocity_*.png'))
    file_list = sorted(file_list, key=lambda x: int(x.split('_')[-1].split('.')[0]))
    gif_images = []
    for file in file_list:
        image = imageio.imread(file)
        gif_images.append(image)
    imageio.mimsave('Figures/DG/DoubleGyres.gif', gif_images, loop=0, duration=5)

    plt.figure()
    plt.plot(t[:-1], Tr_A2[:-1], 'b', linewidth=0.5)
    plt.xlabel(r'Time [$s$]')
    plt.ylabel(r'Tr($A^2$) [$m^2$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.tight_layout()
    # plt.show()
    plt.savefig('Figures/Absolute dispersion coefficient.png')

    plt.figure()
    plt.plot(t[:-1], a2[:-1], 'b', linewidth=0.5)
    plt.xlabel(r'Time [$s$]')
    plt.ylabel(r'Absolute dispersion rate [$m^2/s$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.tight_layout()
    # plt.show()
    plt.savefig('Figures/Absolute dispersion rate.png')

    plt.figure()
    plt.plot(t[:-1], r2[:-1, 0], 'b', linewidth=0.5, label=f'$d0 = {d0[0]} m$')
    plt.plot(t[:-1], r2[:-1, 1], 'r', linewidth=0.5, label=f'$d0 = {d0[1]} m$')
    plt.plot(t[:-1], r2[:-1, 2], 'g', linewidth=0.5, label=f'$d0 = {d0[2]} m$')
    plt.plot(t[:-1], r2[:-1, 3], 'orange', linewidth=0.5, label=f'$d0 = {d0[3]} m$')
    plt.xlabel(r'Time [$s$]')
    plt.ylabel(r'Relative dispersion rate [$m^2/s$]')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlim([5e-2, 1])
    # plt.ylim([1e-5, 1])
    plt.grid(True)
    plt.legend(loc='best',facecolor='white', edgecolor='white', ncol=1, handletextpad=0.2, borderpad=0, columnspacing=0.5)
    plt.tight_layout()
    # plt.show()
    plt.savefig('Figures/Relative dispersion rate.png')

    plt.figure()
    plt.plot(t[:-1], Tr_R2[:-1, 0], 'b', linewidth=0.5, label=f'$d0 = {d0[0]} m$')
    plt.plot(t[:-1], Tr_R2[:-1, 1], 'r', linewidth=0.5, label=f'$d0 = {d0[1]} m$')
    plt.plot(t[:-1], Tr_R2[:-1, 2], 'g', linewidth=0.5, label=f'$d0 = {d0[2]} m$')
    plt.plot(t[:-1], Tr_R2[:-1, 3], 'orange', linewidth=0.5, label=f'$d0 = {d0[3]} m$')
    plt.xlabel(r'Time [$s$]')
    plt.ylabel(r'Tr($R^2$) [$m^2$]')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlim([5e-2, 1])
    # plt.ylim([1e-5, 1])
    plt.grid(True)
    plt.legend(loc='best',facecolor='white', edgecolor='white', ncol=1, handletextpad=0.2, borderpad=0, columnspacing=0.5)
    plt.tight_layout()
    # plt.show()
    plt.savefig('Figures/Relative dispersion.png')

    plt.figure()
    plt.scatter(Tr_R2[:-1, 0], r2[:-1, 0], s=0.5, c='b', label=f'$d0 = {d0[0]} m$')
    plt.scatter(Tr_R2[:-1, 1], r2[:-1, 1], s=0.5, c='r', label=f'$d0 = {d0[1]} m$')
    plt.scatter(Tr_R2[:-1, 2], r2[:-1, 2], s=0.5, c='g', label=f'$d0 = {d0[2]} m$')
    plt.scatter(Tr_R2[:-1, 3], r2[:-1, 3], s=0.5, c='orange', label=f'$d0 = {d0[3]} m$')
    plt.xlabel(r'Tr($R^2$) [$m^2$]')
    plt.ylabel(r'Relative dispersion rate [$m^2/s$]')
    plt.xscale('log')
    plt.yscale('log')
    # plt.xlim([5e-2, 1])
    plt.ylim([1e-5, 1])
    plt.grid(True)
    plt.legend(loc='best',facecolor='white', edgecolor='white', ncol=1, handletextpad=0.2, borderpad=0, columnspacing=0.5)
    plt.tight_layout()
    # plt.show()
    plt.savefig('Figures/Tr vs rate.png')
    
if __name__ == '__main__':
    mesh_X, mesh_Y = meshing()
    # DoubleGyres(X, Y, 0)
    Trajectories()

    breakpoint
