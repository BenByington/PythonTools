import numpy as np
import numpy.fft as fft
try:
	from scipy import ndimage
except ImportError:
	print "SciPy is not installed on this machine.  Any calls to spec will fail gracelessly"

def populateSphere(rad,num):

	out = np.zeros((num,3))

	inc = np.pi * (3. - np.sqrt(5.))
	off = 2. / num

	k = np.array(range(num))
	out[:,1] = k * off - 1. + (off / 2.)

	r = np.sqrt(1 - out[:,1]**2)
	phi = k * inc
	out[:,2] = np.cos(phi)*r
	out[:,0] = np.sin(phi)*r

	def filt(a):
		if a[0] >= 0.0 and a[1] >= 0.0 and a[2] >= 0.0:
			return True
		else:
			return False

	return rad*np.array(filter(filt,out))

def sphereAvg(a,limiter,factors):

	avg = np.zeros(a.shape[limiter])
	avg[0] = a[0,0,0]

	for j in range(1,a.shape[limiter]):
		points = populateSphere(j,2*(j+3)**2)
		points[:,0] *= factors[0]
		points[:,1] *= factors[1]
		points[:,2] *= factors[2]
		vals = ndimage.map_coordinates(a,np.transpose(points),order=1)
		avg[j] = np.average(vals)*j*j

	return avg

def populateCircle(rad,num):

	out = np.zeros((num,2))

	k = np.array(range(num))
	da = 0.5*np.pi / (num-1)

	out[:,0] = rad * np.sin(da * k)
	out[:,1] = rad * np.cos(da * k)

#	for a in out:
#		if a[0] < 0.0 or a[1] < 0.0:
#			print "Broken populateCircle!"
#			print a[0], a[1]
	
	return out

def circleAvg(a,limiter,factors):

	avg = np.zeros(a.shape[limiter])
	avg[0] = a[0,0]

	for j in range(1,a.shape[limiter]):
		points = populateCircle(j,j+3)
		points[:,0] *= factors[0]
		points[:,1] *= factors[1]
		vals = ndimage.map_coordinates(a,np.transpose(points),order=1)
		avg[j] = np.average(vals)*j

	return avg

def spec_comp(u,v,w,r,dims=(1.0,1.0,1.0)):
	u = u * np.sqrt(r)
	v = v * np.sqrt(r)
	w = w * np.sqrt(r)
	return spec(u,v,w,dims=dims)

def spec(u,v,w,dims=(1.0,1.0,1.0)):
	if u.shape != v.shape or u.shape != w.shape:
		print "Incompatible input shapes!"
		print u.shape,v.shape,w.shape
		return

	if len(u.shape) == 3:
		return spec3d(u,v,w,dims=dims)
	elif len(u.shape) == 2:
		return spec2d(u,v,w,dims=dims)
	else:
		print "No spec routine for dimension: ",len(u.shape)

def spec3d(u,v,w,dims=(1.0,1.0,1.0)):
	nx = u.shape[2]
	ny = u.shape[1]
	nz = u.shape[0]

	lens = np.array(dims)
	res = np.array(u.shape)
	minres = res.min()
	limiter = (res/lens).argmin()
	factors = lens / lens[limiter]

	fu = np.abs(fft.ifftn(u))**2 + np.abs(fft.ifftn(v))**2 + np.abs(fft.ifftn(w))**2
	fu = sphereAvg(fu[0:nz/2,0:ny/2,0:nx/2],limiter,factors)

	return fu

def spec2d(u,v,w,dims=(1.0,1.0)):
	ny = u.shape[0]
	nx = u.shape[1]

	lens = np.array(dims)
	res = np.array(u.shape)
	minres = res.min()
	limiter = (res/lens).argmin()
	factors = lens / lens[limiter]

	fu = np.abs(fft.ifftn(u))**2 + np.abs(fft.ifftn(v))**2 + np.abs(fft.ifftn(w))**2
	fu = circleAvg(fu[0:ny/2,0:nx/2],limiter,factors)

	return fu
