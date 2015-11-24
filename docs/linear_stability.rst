Linear Stability Analysis
=========================

One-Layer Shallow Water
-----------------------

In this subsection we consider the one-layer reduced gravity RSW model
with topography below. We define the following:

-  :math:`H`: mean depth of the layer

-  :math:`z=\eta`: height of the free surface

-  :math:`z=-H + \eta_B`: height of the topography.

-  :math:`h = H + \eta - \eta_B`: total depth of layer

-  :math:`(u,v)`: horizontal velocity

-  :math:`g'`: reduced gravity

-  :math:`\rho_0`: reference density

The governing nonlinear equations are,

.. math::

   \begin{aligned}
   \frac{\partial u}{\partial t} + {\vec u} \cdot \vec \nabla u - f v & 
   = - g \frac{\partial}{\partial x} \left( h + \eta_B \right) , \\
    \frac{\partial v}{\partial t} + {\vec u} \cdot \vec \nabla v + f u & 
   = - g \frac{\partial}{\partial y} \left( h + \eta_B \right) , \\
   \frac{\partial h}{\partial t} + \vec\nabla \cdot \left( h \vec u_1 \right) & = 0.\end{aligned}

Basic State
-----------

To study shear flows in a meridional channel we consider solutions of
the form,

.. math::

   \begin{aligned}
   u & = U_B(y), \\
   v & = 0,\\
   h & = H_B(y).\end{aligned}

For this to be an exact solution we require that the flow is in
geostrophic balance,

.. math:: f U_B = - g \frac{d}{dy}\left( H_B + \eta_B \right).

Perturbation
------------

We perturb the basic state with infinitesimal quantities,

.. math::

   \begin{aligned}
   u & = U_B(y) + u', \\
   v & = 0     + v',\\
   h & = H_B(y) + h'.\end{aligned}

We substitute our perturbation into the governing equations and drop the
primes (for brevity) and cancelling out the geostrophic terms

.. math::

   \begin{aligned}
   \frac{\partial u}{\partial t} + (u + U_B) \frac{\partial u}{\partial x} + v \frac{\partial}{\partial y} \left(u + U_B \right)  - f v & = - g \frac{\partial h}{\partial x}, \\
    \frac{\partial v}{\partial t}   + (u + U_B) \frac{\partial v}{\partial x} + v \frac{\partial v }{\partial y} + f u 
    & = - g \frac{\partial h}{\partial y}, \\
   \frac{\partial h}{\partial t}  + (u + U_B) \frac{\partial h}{\partial x}   + v \frac{\partial}{\partial y}(H_B + h)
   + (H_B + h)  & \left( \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} \right) =  0.\end{aligned}

Now we neglect the quadratic terms to obtain the linearized equations,

.. math::

   \begin{aligned}
   \frac{\partial u}{\partial t}
   & = - U_B \frac{\partial u}{\partial x} + \left( f - \frac{d U_B}{d y}  \right) v  - g \frac{\partial h_1}{\partial x}, \\
    \frac{\partial v}{\partial t}    & =  - f u  - U_B \frac{\partial v}{\partial x}  - g \frac{\partial h_1}{\partial y}, \\
   \frac{\partial h}{\partial t}   & = - H_B \frac{\partial u}{\partial x}    - v \frac{d H_B}{d y}
    - H_B \frac{\partial v}{\partial y} - U_B \frac{\partial h}{\partial x} .\end{aligned}

Finally, we assume a normal mode decomposition in the zonal direction
and time,

.. math::

   \begin{aligned}
   [u, v, h] = \mbox{Re}\left\{ e^{ik(x - c t)} [\hat u, ik \hat v, \hat h] \right\},\end{aligned}

which we can substitute into the above equations to yield in the
inviscid limit

.. math::

   \begin{aligned}
   c \hat u  &  = U_B \hat u - (f - \frac{dU_B}{dy}) \hat v + g \hat h, \\
   c \hat v  &  =- \frac{f}{k^2} \hat u + U_B \hat v -\frac{ g}{k^2} \frac{d h}{d y}, \\
   c \hat h  &=   H_B \hat u  +  \frac{d}{d y}\left(H_B \hat v\right) +  U_B \hat h.\end{aligned}


