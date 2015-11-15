all:
	f2py-2.7 -c -m WENOx_periodic FV_x_WENO_periodic.f
	f2py-2.7 -c --f77flags=-fPIC -m WENOy_periodic FV_y_WENO_periodic.f

clean:
	rm *.so *.pyc
