.. figure:: https://github.com/PyRsw/PyRsw/blob/master/docs/_static/vort_jet.png
   :alt:

PyRsw: Python Rotating Shallow Water Model
===================================

PyRsw is a python solver of the Rotating Shallow Water Model.  
This is an approximate model that can be derived from the full
Bossinessq equations in the limit of constant density and strictly
hydrostatic motions. 

It is more general than the one-layer 
quasigeostrophic model in that it allows for variable Rossby numbers
and can include gravity waves. It is a useful model to study
the oceans, the atmosphere and other planetary systems as well.

We hope that this will be of interest both to students and 
researchers in the field. The aim is to create examples that 
can illustrate the basic physical processes and documentation 
that explains it.  Then the user can modify it to study other processes.

PyRsw is a work in progress and we greatly value peoples comments and
contributions.  

Presently, PyRsw is strictly pseudo-spectral but there are plans of including
a Finite Volume option.

PyRws is threaded if you have pyfftw installed.  If this is not present than
the fft’s from numpy will be used instead.   Extending this to include
mpi would make things run even faster.


Links
-----

-  HTML documentation: http://pyrsw.readthedocs.org/en/latest/
-  Issue tracker: https://github.com/PyRsw/PyRsw/issues
-  Source code: https://github.com/PyRsw/PyRsw/tree/master/src
-  pyfftw: http://github.com/hgomersall/pyFFTW

