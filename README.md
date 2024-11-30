# Double-gyres
Single and multiple particles statistics for Double Gyres diffusion coefficients

## Formulation
Double Gyres model is a simple model for ocean circulation. The stream function is given by:

$$\phi(x,y,t) = A\sin(\pi f(x,t))\sin(\pi y),$$

with:

$$f(x,t) = \epsilon\sin(\omega t) x^2 + [1-2\epsilon\sin(\omega t)] x,$$

where $$A = 0.1$$, $$\epsilon = 0.25$$, $$\omega = 2\pi/10$$.

The velocity field is given by:
$$u = -\frac{\partial \phi}{\partial y} = -\pi A \sin(\pi f(x))\cos(\pi y),$$

$$v = \frac{\partial \phi}{\partial x} = \pi A \cos(\pi f(x))\sin(\pi y)\frac{\partial f}{\partial x},$$

with:
$$\frac{\partial f}{\partial x} = 2 \epsilon\sin(\omega t) x + [1-2\epsilon\sin(\omega t)].$$

The Q-criterion is a fluid dynamics criterion based on the vorticity tensor, used to identify vortex cores in flow fields. Its basic idea is to determine the existence of vortex cores by comparing the relative strength of rotational (vortical) and strain (shear) components of the flow. The definition of the Q-criterion is as follows:

$$Q = \frac{1}{2} \left( \parallel\pmb{\Omega}\parallel^2 - \parallel\pmb{S}\parallel^2 \right),$$

where: $$\left.\pmb{\Omega}\right.$$ is the antisymmetric part of the velocity gradient tensor (the vorticity tensor), defined as:

$$\left. \pmb{\Omega} = \frac{1}{2} \left( \nabla \pmb{u} - \nabla \pmb{u}^\top \right). \right.$$

$$\left.\pmb{S}\right.$$ is the symmetric part of the velocity gradient tensor (the strain-rate tensor), defined as:

$$\left.\pmb{S} = \frac{1}{2} \left( \nabla \pmb{u} + \nabla \pmb{u}^\top \right). \right.$$

$$\left.\parallel\pmb{\Omega}\parallel^2\right.$$ and $$\left.\parallel\pmb{S}\parallel^2\right.$$ are the second invariants of the tensors, defined as:

$$\left.\parallel\pmb{\Omega}\parallel^2 = \text{Tr}(\pmb{\Omega} \cdot \pmb{\Omega}^\top), \quad \parallel\pmb{S}\parallel^2 = \text{Tr}(\pmb{S} \cdot \pmb{S}^\top). \right.$$

The Runge-Kutta method is used to solve the system of differential equations for the partical trajectories. Define four slopes: $k_1$ is the slope at the beginning of the interval; $$k_2$$ is the slope at the midpoint of the interval; $$k_3$$ is again the slope at the midpoint; $k_4$ is the slope at the end of the interval. The specific formulations are as follows:

$$\left.
    \begin{aligned}
    \pmb{k_1} &= \Delta t \cdot \left( u\left[x(t_{i}), y(t_{i}), t(i)\right], v\left[x(t_{i}), y(t_{i}), t(i)\right] \right); \\
    \pmb{k_2} &= \Delta t \cdot \left( u\left[x(t_{i}) + \frac{\pmb{k_1}(1)}{2}, y(t_{i}) + \frac{\pmb{k_1}(2)}{2}, t_{i} + \frac{\Delta t}{2}\right], v\left[x(t_{i}) + \frac{\pmb{k_1}(1)}{2}, y(t_{i}) + \frac{\pmb{k_1}(2)}{2}, t_{i} + \frac{\Delta t}{2}\right] \right); \\
    \pmb{k_3} &= \Delta t \cdot \left( u\left[x(t_{i}) + \frac{\pmb{k_2}(1)}{2}, y(t_{i}) + \frac{\pmb{k_2}(2)}{2}, t_{i} + \frac{\Delta t}{2}\right], v\left[x(t_{i}) + \frac{\pmb{k_2}(1)}{2}, y(t_{i}) + \frac{\pmb{k_2}(2)}{2}, t_{i} + \frac{\Delta t}{2}\right] \right); \\
    \pmb{k_4} &= \Delta t \cdot \left( u\left[x(t) + \pmb{k_3}(1), y(t) + \pmb{k_3}(2), t_{i} + \Delta t\right], v\left[x(t) + \pmb{k_3}(1), y(t) + \pmb{k_3}(2), t_{i} + \Delta t\right] \right). \\
    \end{aligned}
\right.$$

Here, $$\Delta t$$ is the time step, and $$t_i$$ represents the current time step index. The final values of $$X$$ and $$Y$$ at the next time step $$t_{i+1}$$ are computed as a weighted average of these intermediate values:

$$\left[x(t_{i+1}), y(t_{i+1})\right] = \left[x(t_{i}), y(t_{i})\right] + \frac{1}{6} \left(\pmb{k_1} + 2\pmb{k_2} + 2\pmb{k_3} + \pmb{k_4}\right).$$

The covariance statistics of particle displacement is defined as:

$$A^2_{ij}(t,t_0) = \frac{1}{M}\sum_{m=1}^{M}[x_i^m(t)-x_i^m(t_0)][x_j^m(t)-x_j^m(t_0)],$$

where $$M$$ is the number of particles, $$x_i^m(t)$$ is the $$i$$-th component of displacement of the $$m$$-th particle at time $$t$$.
The covariance matrix for all directions is:

$$ 
A^2(t,t_0) = \begin{bmatrix}
        A^2_{xx}(t,t_0) & A^2_{xy}(t,t_0) \\
        A^2_{yx}(t,t_0) & A^2_{yy}(t,t_0)
    \end{bmatrix}.
$$
    
The trace of covariance matrix characterizes the level of diffusion in both directions:

$$\text{Tr}[A^2(t,t_0)] = A^2_{xx}(t,t_0) + A^2_{yy}(t,t_0).$$

Then, the absolute diffusion rate is defined as:

$$
\begin{aligned}
    a^2(t) &= \frac{1}{2}\frac{\text{d}}{\text{d}t}\{\text{Tr}[A^2(t,t_0)]\} \\
    &= \frac{1}{2M}\sum_{m=1}^{M}\frac{\text{d}}{\text{d}t}\{[x^m(t)-x^m(t_0)]^2 + [y^m(t)-y^m(t_0)]^2\} \\
    &= \frac{1}{M}\sum_{m=1}^{M}\{u^m(t)[x^m(t)-x^m(t_0)] + v^m(t)[y^m(t)-y^m(t_0)]\}.
    \end{aligned}
$$

Similarly, The covariance statistics of particle relative displacement is defined as:

$$R^2_{ij}(t,d_0) = \frac{1}{M-1}\sum_{m=1}^{M-1}[x_i^m(t)-x_i^{m+1}(t)][x_j^m(t)-x_j^{m+1}(t)],$$

where $$d_0$$ is the initial distance between two particles. The Relative diffusion tensor is:

$$
    R^2(t,d_0) = \begin{bmatrix}
        R^2_{xx}(t,d_0) & R^2_{xy}(t,d_0) \\
        R^2_{yx}(t,d_0) & R^2_{yy}(t,d_0)
    \end{bmatrix}.
$$

The trace of relative diffusion tensor characterizes the level of mean squared relative distance in all directions:

$$\text{Tr}[R^2(t,d_0)] = R^2_{xx}(t,d_0) + R^2_{yy}(t,d_0).$$

Then, the relatice diffusion rate is defined as:

$$
    \begin{aligned}
    r^2(t) &= \frac{1}{2}\frac{\text{d}}{\text{d}t}\{\text{Tr}[R^2(t,d_0)]\} \\
    &= \frac{1}{2(M-1)}\sum_{m=1}^{M-1}\frac{\text{d}}{\text{d}t}\{[x^m(t)-x^{m+1}(t)]^2 + [y^m(t)-y^{m+1}(t)]^2\} \\
    &= \frac{1}{M-1}\sum_{m=1}^{M-1}\{[u^m(t)-u^{m+1}(t)][x^m(t)-x^{m+1}(t)] + [v^m(t)-v^{m+1}(t)][y^m(t)-y^{m+1}(t)]\}.
    \end{aligned}
$$





## Results

### Animation
![image](https://github.com/ZimoJupiter/Double-gyres/blob/main/Figures/DG/DoubleGyres.gif)

### Q criterion
#### t = 0.0s
![image](https://github.com/ZimoJupiter/Double-gyres/blob/main/Figures/Q_criterion_0.0.png)
#### t = 2.5
![image](https://github.com/ZimoJupiter/Double-gyres/blob/main/Figures/Q_criterion_2.5.png)
#### t = 5.0s
![image](https://github.com/ZimoJupiter/Double-gyres/blob/main/Figures/Q_criterion_5.0.png)
#### t = 7.5s
![image](https://github.com/ZimoJupiter/Double-gyres/blob/main/Figures/Q_criterion_7.5.png)
#### t = 10.0s
![image](https://github.com/ZimoJupiter/Double-gyres/blob/main/Figures/Q_criterion_10.0.png)

### Absolut diffusion coefficient and rate
![image](https://github.com/ZimoJupiter/Double-gyres/blob/main/Figures/Absolute%20dispersion%20coefficient.png)
![image](https://github.com/ZimoJupiter/Double-gyres/blob/main/Figures/Absolute%20dispersion%20rate.png)

### Relative diffusion coefficient and rate
![image](https://github.com/ZimoJupiter/Double-gyres/blob/main/Figures/Relative%20dispersion.png)
![image](https://github.com/ZimoJupiter/Double-gyres/blob/main/Figures/Absolute%20dispersion%20rate.png)
![image](https://github.com/ZimoJupiter/Double-gyres/blob/main/Figures/Tr%20vs%20rate.png)
