Layered quasigeostrophic model
==============================

The :math:`{\mathrm{N}}`-layer quasigeostrophic (QG) potential vorticity
is

.. math::

   \begin{aligned}
   {q_1} &= {\nabla^2}\psi_1 + \frac{f_0^2}{H_1} \left(\frac{\psi_{2}-\psi_1}{g'_{1}}\right)\,,  \qquad & i =1{\, ,}\nonumber \\
   {q_n} &= {\nabla^2}\psi_n + \frac{f_0^2}{H_n} \left(\frac{\psi_{n-1}-\psi_n}{g'_{n-1}}  - \frac{\psi_{i}-\psi_{n+1}}{g'_{i}}\right)\,,  \qquad &i = 2,{\mathrm{N}}-1 {\, ,}\nonumber \\
   {q_{\mathrm{N}}} &= {\nabla^2}\psi_{\mathrm{N}}+ \frac{f_0^2}{H_{\mathrm{N}}} \left(\frac{\psi_{\textsf{N}-1}-\psi_{\mathrm{N}}}{g'_{{\mathrm{N}}-1}}\right) + \frac{f_0}{H_{\mathrm{N}}}h_b (x,y)\,,  \qquad & i ={\mathrm{N}}\,,\end{aligned}

where :math:`q_n` is the n’th layer QG potential vorticity, and
:math:`\psi_n` is the streamfunction, :math:`f_0` is the inertial
frequency, n’th :math:`H_n` is the layer depth, and :math:`h_b` is the
bottom topography. (Note that in QG :math:`h_b/H_{\mathrm{N}}<< 1`.)
Also the n’th buoyancy jump (reduced gravity) is

.. math:: g'_n \equiv g \frac{\rho_{n}-\rho_{n+1}}{\rho_n}{\, ,}

where :math:`g` is the acceleration due to gravity and :math:`\rho_n` is
the layer density.

The dynamics of the system is given by the evolution of PV. We introduce
a background flow that can vary in the horizontal. The streamfunction
associated with this flow can be denoted with :math:`\Psi_n(x,y)` for
each layer and geostrophy yields its corresponding velocity
:math:`\vec{V_n} = (U_n(x,y),V_n(x,y))` where :math:`\Psi_{ny} = - U_n`
and :math:`\Psi_{nx} = V_n`. We can perturb the stream function in each
layer into a background flow and deviations from that flow as,

.. math::

   \begin{aligned}
   \psi_n^{{\text{tot}}} = \Psi_n + \psi_n.\end{aligned}

With this basic decomposition we can than write out the corresponding
decompositions in velocity

.. math::

   \begin{aligned}
   \label{eq:Uequiv}
   u_n^{{{\text{tot}}}} = U_n - \psi_{n y}{\, ,}\nonumber \\
   v_n^{{\text{tot}}} = V_n + \psi_{n x} {\, ,}\end{aligned}

and

.. math:: q_n^{{\text{tot}}} = Q_n + \delta_{n{\mathrm{N}}}\frac{f_0}{H_{\mathrm{N}}}h_b + q_n {\, ,}

where :math:`Q_n + \delta_{n{\mathrm{N}}}\frac{f_0}{H_{\mathrm{N}}}h_b`
is n’th layer background PV, we obtain the evolution equations

.. math::

   \begin{aligned}
   \label{eq:qg_dynamics}
   {q_n}_t + \mathsf{J}(\psi_n,q_n + \delta_{n {\mathrm{N}}} \frac{f_0}{H_{\mathrm{N}}}h_b )& + U_n ({q_n}_x + \delta_{n {\mathrm{N}}} \frac{f_0}{H_{\mathrm{N}}}h_{bx}) + V_n ({q_n}_y + \delta_{n {\mathrm{N}}} \frac{f_0}{H_{\mathrm{N}}}h_{by})+ \nonumber
   \\ & {Q_n}_y {\psi_n}_x - {Q_n}_x {\psi_n}_y = {\text{ssd}}- r_{ek} \delta_{n{\mathrm{N}}} {\nabla^2}\psi_n {\, ,}\qquad n = 1,{\mathrm{N}}{\, ,}\end{aligned}

where :math:`{\text{ssd}}` is stands for small scale dissipation, which
is achieved by an spectral exponential filter or hyperviscosity, and
:math:`r_{ek}` is the linear bottom drag coefficient. The Dirac delta,
:math:`\delta_{nN}`, indicates that the drag is only applied in the
bottom layer.

Linear Stability Analysis
-------------------------

In order to study the stability of a jet in the context of our
:math:`n`-layer QG model we focus our attention on basic states that
consist of zonal flows. i.e. :math:`\Psi_n(y)` only. If we assume that
the quadratic quantities we can then linearize to obtain in the
conservative limit over a flat bottom,

.. math::

   \begin{aligned}
   \label{eq:qglin_dynamics}
   {q_n}_t  + U_n {q_n}_x + {Q_n}_y {\psi_n}_x  = 0,\end{aligned}

for :math:`n = 1, \cdots, N`.

We assume that the perturbations are normal modes in the zonal direction
and time,

.. math:: \psi_n  = {\mathrm{Re}}[ \hat \psi_n e^{i(kx - \omega t)} ].

This implies that the PV will be modified appropriately and we denote it
with :math:`\hat q_n`.

We substitute this into the linear equations and then divide by the
exponential to obtain,

.. math::

   \begin{aligned}
   c  {\hat q_n}  =  U_n {\hat q_n} +  {Q_n}_y {\hat \psi_n} ,\end{aligned}

where the basic state only depends on :math:`y`, and layer of course,
and we have introduced the phase speed :math:`c=\omega/k`. Note that the
actual PVs are

.. math::

   \begin{aligned}
   {\hat q_1} &= (\partial_{yy} - k^2) \hat \psi_1 + \frac{f_0^2}{H_1} \left(\frac{\hat \psi_{2}-\hat \psi_1}{g'_{1}}\right)\,,  \qquad & i =1{\, ,}\nonumber \\
   {\hat q_n} &= (\partial_{yy} - k^2)\psi_n + \frac{f_0^2}{H_n} \left(\frac{\hat \psi_{n-1}-
   \hat \psi_n}{g'_{n-1}}  - \frac{\hat \psi_{i}-\hat \psi_{n+1}}{g'_{i}}\right)\,,  \qquad &i = 2,{\mathrm{N}}-1 {\, ,}\nonumber \\
   {\hat q_{\mathrm{N}}} &= (\partial_{yy} - k^2)\hat \psi_{\mathrm{N}}+ \frac{f_0^2}{H_{\mathrm{N}}} \left(\frac{\hat \psi_{\textsf{N}-1} - \hat \psi_{\mathrm{N}}}{g'_{{\mathrm{N}}-1}}\right)\,,  \qquad & i ={\mathrm{N}}\,,\end{aligned}

Special case: one-layer model
=============================

In the one-layer case we have

.. math::

   \begin{aligned}
   c  {\hat q_1}  =  U_1 {\hat q_1} +  {Q_1}_y {\hat \psi_1} ,\end{aligned}

.. math::

   \begin{aligned}
   {\hat q_1} = \left[ \partial_{yy} - k^2 -  \frac{f_0^2}{g'_1 H_1} \right] \hat \psi_1.\end{aligned}

Special case: two-layer model
=============================

In the two-layer case we have

.. math::

   \begin{aligned}
   c  {\hat q_n}  =  U_n {\hat q_n} +  {Q_n}_y {\hat \psi_n} ,\end{aligned}

.. math::

   \begin{aligned}
   {\hat q_1} &= \left[ \partial_{yy} - k^2 -  \frac{f_0^2}{g'_1 H_1}\right] \hat \psi_1 + \frac{f_0^2}{g'_1 H_1} \hat \psi_{2}, \\
   {\hat q_2} &= \frac{f_0^2}{g'_1 H_2}\hat \psi_1+ \left[ \partial_{yy} - k^2 -  \frac{f_0^2}{g'_1 H_2} \right] \hat \psi_2 
   .\end{aligned}


