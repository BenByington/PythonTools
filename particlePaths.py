from scipy.interpolate import RectBivariateSpline
import numpy as np

#assumes a periodic field
def integrate(xloc, zloc, xfield, zfield, duration, xmx=1.0, zmx=1.0):


	xloc = np.array(xloc)
	zloc = np.array(zloc)

	time = 0
	dt = duration / 1000.
	s = xfield.shape
	nz = s[0]
	nx = s[2]
	tol = np.min([xmx/nx,zmx/nz])/10.

	xpath = np.array(xloc,copy=True,ndmin=2)
	zpath = np.array(zloc,copy=True,ndmin=2)

	tempxf = __wrapField__(xfield)
	tempzf = __wrapField__(zfield)
	xc,zc = __getCoords__(xfield,xmx,zmx)
#	xc = xmx * np.array(range(nx))/nx
#	zc = zmx * np.array(range(nz))/nz

#	xfunc = RectBivariateSpline(zc,xc,xfield[:,0,:])
#	zfunc = RectBivariateSpline(zc,xc,zfield[:,0,:])
	print zc.shape
	print xc.shape
	print tempxf.shape
	xfunc = RectBivariateSpline(zc,xc,tempxf[:,0,:])
	zfunc = RectBivariateSpline(zc,xc,tempzf[:,0,:])

#	while (time <  duration):
	for i in range(10000):
		print i
		sx,sz = __rungeStep__(xloc,zloc,xfunc,zfunc,dt,xmx,zmx)
		dx,dz = __rungeStep__(xloc,zloc,xfunc,zfunc,dt/2.0,xmx,zmx)
		dx,dz = __rungeStep__(dx,dz,xfunc,zfunc,dt/2.0,xmx,zmx)

		error = np.max(np.sqrt((sx - dx)**2 + (sz - dz)**2))

		if(error == 0.0):
			mult = 10
		else:
			mult = (tol / error)**(2.0/3.0)
		if mult > 100:
			dt = dt * 100
		else:
			dt = dt * mult

		if error < tol:
			xloc = dx
			zloc = dz
			time = time + dt
		 	xpath = np.append(xpath,np.expand_dims(xloc,axis=0),axis=0)
			zpath = np.append(zpath,np.expand_dims(zloc,axis=0),axis=0)

	return [xpath,zpath]

def __rungeStep__(xloc, zloc, xfunc, zfunc, dt, xmx, zmx):
    
	k1 = []
	k1.append(xfunc.ev(zloc,xloc))
	k1.append(zfunc.ev(zloc,xloc))
    
    #move to the middle using the derivative at the beginging, and take the derivative there
	tempx = xloc + dt *  k1[0] / 2.0
	tempz = zloc + dt *  k1[1] / 2.0;

	tempx[tempx > xmx] -= xmx
	tempx[tempx < 0] += xmx
	tempz[tempz > zmx] -= zmx
	tempz[tempz < 0] += zmx

	k2 = []
	k2.append(xfunc.ev(tempz,tempx))
	k2.append(zfunc.ev(tempz,tempx))
    
    #move to the middle again using the new derivative found, and take the derivative there
	tempx = xloc + dt * k2[0] / 2.0;
	tempz = zloc + dt * k2[1] / 2.0;

	tempx[tempx > xmx] -= xmx
	tempx[tempx < 0] += xmx
	tempz[tempz > zmx] -= zmx
	tempz[tempz < 0] += zmx

	k3 = []
	k3.append(xfunc.ev(tempz,tempx))
	k3.append(zfunc.ev(tempz,tempx))

    #move to the end using the newest derivative, and take the derivative there
	tempx = xloc + dt * k3[0];
	tempz = zloc + dt * k3[1];

	tempx[tempx > xmx] -= xmx
	tempx[tempx < 0] += xmx
	tempz[tempz > zmx] -= zmx
	tempz[tempz < 0] += zmx

	k4 = []
	k4.append(xfunc.ev(tempz,tempx))
	k4.append(zfunc.ev(tempz,tempx))

    #use a weighted average of the derivatives to make the integration step
	retx = xloc + dt * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6.0
	retz = zloc + dt * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6.0
    
	retx[retx > xmx] -= xmx
	retx[retx < 0] += xmx
	retz[retz > zmx] -= zmx
	retz[retz < 0] += zmx

	return [retx,retz]

def __wrapField__(field):
	s = field.shape
	nz = s[0]
	nx = s[2]

	xl = field[:,:,0:5]
	xr = field[:,:,nx-5:]

	tempf = np.append(xr,field,2)
	tempf = np.append(tempf,xl,2)

	zd = tempf[0:5,:,:]
	zu = tempf[nz-5:,:,:]

	tempf = np.append(zu,tempf,0)
	tempf = np.append(tempf,zd,0)

	return tempf

def __getCoords__(field,xmx,zmx):
	
	s = field.shape
	nz = s[0]
	nx = s[2]

	xc = xmx * np.array(range(nx+10),dtype=np.float64)/nx
	xc = xc - xc[5]
	zc = zmx * np.array(range(nz+10),dtype=np.float64)/nz
	zc = zc - zc[5]

	return xc,zc
