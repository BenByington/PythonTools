import numpy as np
import gc
from bashcmd import *

def get_rec(name,shape,type,iswap=False):
	dat = np.fromfile(name,dtype=type,count=-1)

	if shape != None:
		len = dat.size
		slen = 1
		for d in shape:
			slen *= d
		if len != slen:
			raise NameError("Incompatible Dimensions for Data File " \
                + name + ": " + str(shape) + " " + str(len))      
		#Can run into some memory issues when a few large things are loaded
		#into memory.  This is a stopgap to try and mitigate that
		if(len > 2097152):
			gc.collect()
   
 		dat.shape = shape

	if iswap == True:
		return dat.byteswap()
	else:
		return dat

def save_rec(name,data):
	data.tofile(name)

def get_b(nx=96,ny=96,nz=94,type=np.float32,iswap=False):
	bx = get_rec('Bx',(nz,ny,nx),type,iswap)
	by = get_rec('By',(nz,ny,nx),type,iswap)
	bz = get_rec('Bz',(nz,ny,nx),type,iswap)
	b2 = bx**2 + by**2 + bz**2

	return (bx,by,bz,b2)

def Iget_b(nx=96,ny=96,nz=96,type=np.float64,iswap=False):
	return get_b(nx=nx,ny=ny,nz=nz,type=type,iswap=iswap)

def Iget_v(nx=96,ny=96,nz=96,t=False,type=np.float64,iswap=False):
	u = get_rec('u',(nz,ny,nx),type,iswap)
	v = get_rec('v',(nz,ny,nx),type,iswap)
	w = get_rec('w',(nz,ny,nx),type,iswap)

	if t == True:
		t = get_rec('T',(nz,ny,nx),type,iswap)
		return (u,v,w,t)
	else:
		return (u,v,w)

def get_v(nx=96,ny=96,nz=94,type=np.float32,iswap=False):
	u = get_rec('u',(nz,ny,nx),type,iswap)
	v = get_rec('v',(nz,ny,nx),type,iswap)
	w = get_rec('w',(nz,ny,nx),type,iswap)
	t = get_rec('temperature',(nz,ny,nx),type,iswap)
	r = get_rec('density',(nz,ny,nx),type,iswap)

	return (u,v,w,t,r)

def make_scalar_accum(nstart=None,nstop=None,nstep=None,nfreq=50,nrec=None, \
					  n_scalars=40,list=None, type=np.float32,zcount=6, \
					  forceNew = False):

	if list == None:
		if(nstart == None or nstop == None or nstep == None):
			list = get_data_files(nz=zcount)
			if nrec == None:
				nrec = (int(list[1]) - int(list[0]))/nfreq
		else:
			list=range(nstart,nstop+nstep,nstep)
			for i in range(len(list)):
				list[i] = str(list[i]).zfill(zcount)
			if nrec == None:
				nrec = nstep / nfreq
	name = str(list[0]) + "_" + str(list[-1]) + ".accm"
	shape = (nrec*len(list),n_scalars)

	if(forceNew == False):
		try:
			ret = get_rec(name,shape,type)
			return ret
		except NameError:
			pass
		except IOError:
			pass
		
	ret = np.zeros(shape,dtype=np.float32)
	index = int(0)
	for f in list:
		print "Doing file: " + f + " of " + list[-1]
		temp = get_rec(f,(nrec,n_scalars),type)
		ret[index*nrec:(index+1)*nrec,:] = temp
		index += 1
	save_rec(name,ret)
	return ret

def make_flux_accum(nstart=None,nstop=None,nstep=None,nfreq=50,nrec=None, \
					  nz=94, n_fluxes=40,list=None,type=np.float32,zcount=6, \
					forceNew = False):

	if list == None:
		if(nstart == None or nstop == None or nstep == None):
			list = get_data_files(nz=zcount)
			if nrec == None:
				nrec = int(list[0])/nfreq
		else:
			list=range(nstart,nstop+nstep,nstep)
			for i in range(len(list)):
				list[i] = str(list[i]).zfill(zcount)
			if nrec == None:
				nrec = nstep / nfreq
	name = str(list[0]) + "_" + str(list[-1]) + ".accm"
	shape = (nrec * len(list),n_fluxes, nz)

	if(forceNew == False):
		try:
			ret = get_rec(name,shape,type)
			return ret
		except NameError:
			pass
		except IOError:
			pass

	ret = np.zeros(shape,dtype=type)
	index = int(0)
	for f in list:
		print "Doing file: " + f + " of " + list[-1]
		temp = get_rec(f,(nrec,6,n_fluxes,nz),type)
		ret[index*nrec:(index+1)*nrec,:] = temp[:,0,:,:]
		index += 1

	save_rec(name,ret)
	return ret

