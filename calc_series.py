import numpy as np
try:
	import matplotlib.pyplot as plt
	import colormaps as cm
except:
	print "matplotlib not available on this machine.  Prepare to fail gracelessly if used!"

import os
from spec import spec
from filter import filter_conv
from filter import filter_single
from physics import *
#from tubeDynamics import tDynamics
#from tubeDynamics2 import tDynamics2


class sim_operation(object):
"""This file is a cluttered mess, but the principle isn't too tough.  Any time
   you wish to do a calculation or make a movie over a set of 3D data snapshots
   from a simulation, you need to use/define one of these objects here.  This
   is the base class, which has the various methods you'll need to impliment
   in the child class in order for the imhd object to interact with it as 
   desired.  
"""
	sim = None

	def __init__(self):
		self.sim = None
		self.suffix = None

	def __call__(self):
"""EVERY child MUST have a call routine implimented.  The first argument will
   always be self, after that it is whatever you need (e.g. bx,by and bz).  This
   routine will either make calls to matplotlib to create a frame when doing an
   animation, or return a scalar or small vector of data when doing a
   calculation.
"""
		raise NotImplementedError()

	def get_params(self):
"""This is the most fragile part of this operation.  Here you must define a 
   routine that returns a string describing the arguments needed to call the
   __call__ function.  See below for examples on how to do this.
"""
		raise NotImplementedError()

	def get_ret_dim(self):
"""IF defining an object to do caluclations, here you must define how many 
   elements are in your return vector.  e.g. if you calculate the mean emf in
   the box, then you would return a 3.
""" 
		raise NotImplementedError()

class getIteration(sim_operation): 
	def get_params(self):
		return ['iteration']
	def __call__(self,iteration):
		return iteration
	def get_ret_dim(self):
		return 1
	

class get_x_var(sim_operation): 
	def get_params(self):
		return ['iteration','u']
	def __call__(self,iteration,u):
		return [self.sim.get_time(int(iteration)),np.max(self.sim.ddx(u))]
	def get_ret_dim(self):
		return 2

class instability_hunt(sim_operation): 
	def get_params(self):
		return ['u','v','w','b2']
	def get_ret_dim(self):
		return 11
	def __call__(self,u,v,w,b2):

		def center_mass(field):
			(nz,ny) = field.shape
		
			cz = 0
			cy = 0
			for i in range(nz):
				cz += (field[i,:]*i).sum()
			for i in range(ny):
				cy += (field[:,i]*i).sum()
			total = field.sum()
			return [cz/total, cy/total]

		def center_extent(field,cz,cy):
			(nz,ny) = field.shape

			m2 = 0
			for i in range(ny):
				for j in range(nz):
					dist = np.sqrt((i-cy)**2.+ (j-cz)**2.)
					m2 += dist * field[j,i]
			return m2/field.sum()
	
		def average_circle(field,cz,cy,rad):
			(nz,ny) = field.shape
			avg = 0
			count = 0
			for i in range(ny):
				for j in range(nz):
					dist = np.sqrt((i-cy)**2. + (j-cz)**2.)
					if dist < rad:
						avg += field[j,i]
						count += 1
			return avg / count

		avgv = np.average(v,2)
		avgw = -np.average(w,2)
		avgb2 = np.average(b2,2)
		xd = self.sim.ddx(u)
		xd2d = np.max(xd,2)
		mvar = np.max(xd2d)

		yshear = np.max(np.abs(self.sim.ddy(u)),2)
		zshear = np.max(np.abs(self.sim.ddz(u)),2)
		tshear = np.max(self.sim.ddz(u)**2 + self.sim.ddy(u)**2,2)
		div = np.average(self.sim.ddz(w) + self.sim.ddy(v),2)

		[varz,vary] = center_mass(xd2d[40:80,:])
		varSig = center_extent(xd2d[40:80,:],varz,vary)
		varz += 40

		mv = average_circle(avgv,varz,vary,varSig)
		mw = average_circle(avgw,varz,vary,varSig)
		mdiv = average_circle(div,varz,vary,varSig)
		mshy = average_circle(yshear,varz,vary,varSig)
		mshz = average_circle(zshear,varz,vary,varSig)
		msht = average_circle(tshear,varz,vary,varSig)
		mb2 = average_circle(avgb2,varz,vary,varSig)

		return [mvar, varz, vary, varSig, mv, mw, mdiv, mshy, mshz, msht, mb2]

class xavg_mag_gif(sim_operation):
	x = None
	y = None

	def __init__(self):
		x = None
		y = None

	def get_params(self):
		return ["v", "w", "bx", "by", "bz","b2"]

	def __call__(self,v,w,bx,by,bz,b2):
		q = 4

		map = cm.red_blue()
		if self.x == None:
			ny = v.shape[1]
			nz = v.shape[0]
			self.x,self.y = np.meshgrid(range(ny),range(nz))
		x,y = self.x,self.y

		avgv = np.average(v,2)
		avgw = -np.average(w,2)
		avgbx = np.average(bx,2)/np.sqrt(np.max(b2))
		avgby = np.average(by,2)/np.sqrt(np.max(b2))
		avgbz = -np.average(bz,2)
		avgb2 = np.average(b2,2)

		avgbx[0,0] = -np.max(avgbx)
		avgbx[0,1] = -np.min(avgbx)
		avgby[0,0] = -np.max(avgby)
		avgby[0,1] = -np.min(avgby)
		avgbz[0,0] = -np.max(avgbz)
		avgbz[0,1] = -np.min(avgbz)

		plt.subplot(221)
		plt.imshow(avgb2,cmap=map)
		plt.colorbar()
		plt.quiver(x[::q,::q],y[::q,::q],avgv[::q,::q],avgw[::q,::q])
		plt.title('B2')
		plt.axis("tight")

		plt.subplot(222)
		plt.imshow(avgbx,cmap=map)
		plt.colorbar()
		plt.quiver(x[::q,::q],y[::q,::q],avgv[::q,::q],avgw[::q,::q])
		plt.title('Bx')
		plt.axis("tight")

		plt.subplot(223)
		plt.imshow(avgby,cmap=map)
		plt.colorbar()
		plt.quiver(x[::q,::q],y[::q,::q],avgv[::q,::q],avgw[::q,::q])
		plt.title('By')
		plt.axis("tight")

		plt.subplot(224)
		plt.imshow(avgbz,cmap=map)
		plt.colorbar()
		plt.quiver(x[::q,::q],y[::q,::q],avgv[::q,::q],avgw[::q,::q])
		plt.title('Bz')
		plt.axis("tight")

class yavg_mag_gif(sim_operation):
	x = None
	y = None

	def __init__(self):
		x = None
		y = None

	def get_params(self):
		return ["u", "w", "bx", "by", "bz","b2"]

	def __call__(self,u,w,bx,by,bz,b2):
		q = 16

		map = cm.red_blue()
		if self.x == None:
			nx = u.shape[2]
			nz = u.shape[0]
			self.x,self.y = np.meshgrid(range(nx),range(nz))
		x,y = self.x,self.y

		avgu = np.average(u,1)
		avgw = np.average(w,1)
		avgbx = np.average(bx,1)
		avgby = np.average(by,1)
		avgbz = np.average(bz,1)
		avgb2 = np.average(b2,1)

		plt.subplot(221)
		plt.imshow(avgb2,cmap=map,origin='lower')
		plt.colorbar()
		plt.quiver(x[::q,::q],y[::q,::q],avgu[::q,::q],avgw[::q,::q])
		plt.title('B2')
		plt.axis("tight")

		plt.subplot(222)
		plt.imshow(avgbx,cmap=map,origin='lower')
		plt.colorbar()
		plt.quiver(x[::q,::q],y[::q,::q],avgu[::q,::q],avgw[::q,::q])
		plt.title('Bx')
		plt.axis("tight")

		plt.subplot(223)
		plt.imshow(avgby,cmap=map,origin='lower')
		plt.colorbar()
		plt.quiver(x[::q,::q],y[::q,::q],avgu[::q,::q],avgw[::q,::q])
		plt.title('By')
		plt.axis("tight")

		plt.subplot(224)
		plt.imshow(avgbz,cmap=map,origin='lower')
		plt.colorbar()
		plt.quiver(x[::q,::q],y[::q,::q],avgu[::q,::q],avgw[::q,::q])
		plt.title('Bz')
		plt.axis("tight")

class flux_gif(sim_operation):
	x = None
	y = None

	def __init__(self,q=8):
		x = None
		y = None
		self.q = q
		self.suffix = "flux"

	def get_params(self):
		return ["u", "w", "by"]

	def __call__(self,u,w,by):
		q = self.q

		map = cm.red_blue()
		if self.x == None:
			nx = u.shape[2]
			nz = u.shape[0]
			self.x,self.y = np.meshgrid(np.array(range(nx))*self.sim.xmx/nx,np.array(range(nz))*self.sim.zmx/nz)
		x,y = self.x,self.y

		avgu = np.average(u,1)
		avgw = np.average(w,1)
		avgby = np.average(by,1)

		plt.imshow(avgby,cmap=map,origin='lower',extent=[0,self.sim.xmx,0,self.sim.zmx])
		plt.quiver(x[::q,::q],y[::q,::q],avgu[::q,::q],avgw[::q,::q])
		plt.axis("tight")

class inc_flux_gif(sim_operation):
	x = None
	y = None

	def __init__(self,q=8):
		x = None
		y = None
		self.q = q

	def get_params(self):
		return ["u", "w", "bx", "by", "bz","b2"]

	def __call__(self,u,w,bx,by,bz,b2):
		q = self.q

		map = cm.red_blue()
		if self.x == None:
			nx = u.shape[2]
			nz = u.shape[0]
			self.x,self.y = np.meshgrid(range(nx),range(nz))
		x,y = self.x,self.y

		avgu = np.average(u,1)
		avgw = np.average(w,1)
		avgbx = np.average(bx,1)
		avgby = np.average(by,1)
		avgbz = np.average(bz,1)
		avgb2 = np.average(b2,1)

		plt.subplot(121)
		plt.imshow(avgby,cmap=map,origin='lower')
		plt.colorbar()
		plt.quiver(x[::q,::q],y[::q,::q],avgu[::q,::q],avgw[::q,::q])
		plt.title('By-Vel')
		plt.axis("tight")

		plt.subplot(122)
		plt.imshow(avgby,cmap=map,origin='lower')
		plt.colorbar()
		plt.quiver(x[::q,::q],y[::q,::q],avgbx[::q,::q],avgbz[::q,::q])
		plt.title('By-Twist')
		plt.axis("tight")

class tracer_flux_gif(sim_operation):
	x = None
	y = None

	def __init__(self):
		x = None
		y = None
		self.suffix = None

	def get_params(self):
		return ["u", "w", "bx", "by", "bz","b2","t"]

	def __call__(self,u,w,bx,by,bz,b2,t):
		q = 8

		map = cm.red_blue()
		if self.x == None:
			nx = u.shape[2]
			nz = u.shape[0]
			self.x,self.y = np.meshgrid(range(nx),range(nz))
		x,y = self.x,self.y

		avgu = np.average(u,1)
		avgw = np.average(w,1)
		avgbx = np.average(bx,1)
		avgby = np.average(by,1)
		avgbz = np.average(bz,1)
		avgb2 = np.average(b2,1)
		avgt = np.average(t,1)

		plt.subplot(121)
		plt.imshow(avgt,cmap=map,origin='lower')
		plt.colorbar()
		plt.quiver(x[::q,::q],y[::q,::q],avgu[::q,::q],avgw[::q,::q])
		plt.title('Tracer-Vel')
		plt.axis("tight")

		plt.subplot(122)
		plt.imshow(avgby,cmap=map,origin='lower')
		plt.colorbar()
		plt.quiver(x[::q,::q],y[::q,::q],avgbx[::q,::q],avgbz[::q,::q])
		plt.title('By-Twist')
		plt.axis("tight")

class vort_forces_gif(sim_operation):
	x = None
	y = None

	def __init__(self,v,vh,bulk_freq):
		self.x = None
		self.y = None
		self.v = v
		self.vh = vh
		self.bf = bulk_freq

	def get_params(self):
		return ["iteration","u", "w", "bx", "by", "bz","b2","t"]

	def __call__(self,iteration,u,w,bx,by,bz,b2,t):
		q = 8

		index = int(iteration)/self.bf - 1
		dat = self.v[index]
		datH = self.vh[index]

		norm = max([dat[8],datH[8]])

		map = cm.red_blue()
		if self.x == None:
			nx = u.shape[2]
			nz = u.shape[0]
			self.x,self.y = np.meshgrid(range(nx),range(nz))
		x,y = self.x*self.sim.xmx / self.sim.nx,self.y * self.sim.zmx / self.sim.nz

		avgu = np.average(u,1)
		avgw = np.average(w,1)
		avgbx = np.average(bx,1)
		avgby = np.average(by,1)
		avgbz = np.average(bz,1)
		avgb2 = np.average(b2,1)
		avgt = np.average(t,1)

#		plt.subplot(121)
#		plt.imshow(avgt,cmap=map,origin='lower',extent=(0,self.sim.xmx,0,self.sim.zmx))
#		plt.colorbar()
#		plt.quiver(x[::q,::q],y[::q,::q],avgu[::q,::q],avgw[::q,::q])
#		plt.title('Tracer-Vel')
#		plt.axis("tight")

#		plt.subplot(122)
		plt.imshow(avgby,cmap=map,origin='lower',extent=(0,self.sim.xmx,0,self.sim.zmx))
		plt.arrow(dat[1]+self.sim.xmx/50.,dat[2]+self.sim.zmx/50,datH[6]/norm/3,datH[7]/norm/3,width=.02,head_width=.06,head_length=.02,length_includes_head=True,color='green')
		plt.arrow(dat[1],dat[2],dat[6]/norm/3,dat[7]/norm/3,width=.02,head_width=.06,head_length=.02,length_includes_head=True,color="black")
		plt.colorbar()
		plt.title('Vortex Forces')
		plt.axis("tight")

class shear_gif(sim_operation):
	x = None
	y = None
	iterations = None
	xvar = None

	def __init__(self):
		x = None
		y = None
		iterations = None
		xvar = None

	def get_params(self):
		return ["u","v","w","iteration"]

	def __call__(self,u,v,w,iteration):
		q = 4

		plt.cool()
		if self.x == None:
			ny = v.shape[1]
			nz = v.shape[0]
			self.x,self.y = np.meshgrid(range(ny),range(nz))
		x,y = self.x,self.y

		if self.iterations == None:
			self.iterations = self.sim.bulk_calc(getIteration())
		all_itr = self.iterations

		if self.xvar == None:
			class temp(sim_operation):
				def get_params(self):
					return ["u"]
				def __call__(self,u):
					return np.max(self.sim.ddx(u))

			self.xvar = self.sim.bulk_calc(temp())
		xvar_series = self.xvar

		min = np.min(xvar_series)
		max = np.max(xvar_series)
		if min <= 0:
			min = 0.000001
		if max <= min:
			max = 0.00001
	
		avgv = np.average(v,2)
		avgw = -np.average(w,2)
		xd = self.sim.ddx(u)
		xd2d = np.max(xd,2)
		yshear = np.max(np.abs(self.sim.ddy(u)),2)
		zshear = np.max(np.abs(self.sim.ddz(u)),2)

		plt.subplot(221)
		plt.imshow(yshear)
		plt.colorbar()
		plt.title('du/dy')
		plt.axis("tight")

		plt.subplot(222)
		plt.imshow(xd2d)
		plt.colorbar()
		plt.quiver(x[::q,::q],y[::q,::q],avgv[::q,::q],avgw[::q,::q])
		plt.title('Max x Variation (y-z)')
		plt.axis("tight")

		plt.subplot(223)
		plt.imshow(zshear)
		plt.colorbar()
		plt.title('du/dz')
		plt.axis("tight")

		plt.subplot(224)
		plt.plot(all_itr,xvar_series, '--')
		plt.plot([iteration,iteration],[min,max])
		plt.semilogy()
		plt.title('Max x Variation (t)')
		plt.axis("tight")


class instab_gif(sim_operation):
	x = None
	y = None
	iterations = None
	xvar = None

	def __init__(self):
		x = None
		y = None
		iterations = None
		xvar = None

	def get_params(self):
		return ["u","v","w","iteration"]

	def __call__(self,u,v,w,iteration):
		q = 4

		plt.cool()
		if self.x == None:
			ny = v.shape[1]
			nz = v.shape[0]
			self.x,self.y = np.meshgrid(range(ny),range(nz))
		x,y = self.x,self.y

		if self.iterations == None:
			self.iterations = self.sim.bulk_calc(getIteration())
		all_itr = self.iterations

		if self.xvar == None:
			class temp(sim_operation):
				def get_params(self):
					return ["u"]
				def __call__(self,u):
					return np.max(self.sim.ddx(u))

			self.xvar = self.sim.bulk_calc(temp())
		xvar_series = self.xvar

		min = np.min(xvar_series)
		max = np.max(xvar_series)
		if min <= 0:
			min = 0.000001
		if max <= min:
			max = 0.00001
	
		avgu = np.average(u,2)
		avgv = np.average(v,2)
		avgw = -np.average(w,2)
		xd = self.sim.ddx(u)
		xd2d = np.max(xd,2)
		xd1d = np.max(xd2d,1)

		plt.subplot(221)
		plt.imshow(avgu)
		plt.quiver(x[::q,::q],y[::q,::q],avgv[::q,::q],avgw[::q,::q])
		plt.title('Avg u')
		plt.axis("tight")

		plt.subplot(222)
		plt.imshow(xd2d)
		plt.title('Max x Variation (y-z)')
		plt.axis("tight")

		plt.subplot(223)
		plt.plot(xd1d)
		plt.title('Max x Variation (z)')
		plt.axis("tight")

		plt.subplot(224)
		plt.plot(all_itr,xvar_series, '--')
		plt.plot([iteration,iteration],[min,max])
		plt.semilogy()
		plt.title('Max x Variation (t)')
		plt.axis("tight")

class shear_steep_gif(sim_operation):
	iterations = None
	xvar = None

	def __init__(self):
		iterations = None
		xvar = None

	def get_params(self):
		return ["u","iteration"]

	def __call__(self,u,iteration):

		plt.cool()

		if self.iterations == None:
			self.iterations = self.sim.bulk_calc(getIteration())
		all_itr = self.iterations

		if self.xvar == None:
			class temp(sim_operation):
				def get_params(self):
					return ["u"]
				def __call__(self,u):
					return np.max(self.sim.ddx(u))

			self.xvar = self.sim.bulk_calc(temp())
		xvar_series = self.xvar

		min = np.min(xvar_series)
		max = np.max(xvar_series)
		if min <= 0:
			min = 0.000001
		if max <= min:
			max = 0.00001
	
		uavg = np.average(u,2)
		xd = self.sim.ddx(u)
		xd2d = np.max(xd,2)

		plt.subplot(241)
		plt.imshow(xd2d)
		plt.title('Max x Variation (y-z)')
		plt.axis("tight")

		plt.subplot(242)
		plt.plot(uavg[58,:])
		plt.title('z = 58')
		plt.axis("tight")

		plt.subplot(243)
		plt.plot(uavg[60,:])
		plt.title('z = 60')
		plt.axis("tight")

		plt.subplot(244)
		plt.plot(uavg[62,:])
		plt.title('z = 62')
		plt.axis("tight")

		plt.subplot(245)
		plt.plot(all_itr,xvar_series, '--')
		plt.plot([iteration,iteration],[min,max])
		plt.semilogy()
		plt.title('Max x Variation (t)')
		plt.axis("tight")

		plt.subplot(246)
		plt.plot(uavg[64,:])
		plt.title('z = 64')
		plt.axis("tight")

		plt.subplot(247)
		plt.plot(uavg[64,:])
		plt.title('z = 64')
		plt.axis("tight")

		plt.subplot(248)
		plt.plot(uavg[66,:])
		plt.title('z = 66')
		plt.axis("tight")

class analyzeJ(sim_operation):
	def __init__(self,sim,divide=2./3):
		self.nd = np.floor(sim.nz*divide)
 
	def get_params(self):
		return ['iteration','bx','by','bz','b2']
	def get_ret_dim(self):
		return 7
	def __call__(self,iteration,bx,by,bz,b2):
		x,y,z = self.sim.curl(bx,by,bz)
		jmag = np.sqrt(x**2 + y**2 + z**2)
		bmag = np.sqrt(b2)

		nz = self.sim.nz

		lenAll = bmag.sum() / jmag.sum()
		lenTop = bmag[0:self.nd,:,:].sum() / jmag[0:self.nd,:,:].sum()
		lenBot = bmag[self.nd:nz-1,:,:].sum() / jmag[self.nd:nz-1,:,:].sum()

		avgJ = np.average(jmag)
		layerJ = np.average(np.average(jmag,axis=2),axis=1)
		weight = np.array(range(nz))
		
		cm = (layerJ * weight).sum()/layerJ.sum()
		stddev = (layerJ * weight * weight).sum()/layerJ.sum()
		stddev = np.sqrt(stddev - cm*cm)

		return [self.sim.get_time(iteration), lenAll, lenTop, lenBot, avgJ, cm, stddev]
		
class vert_tense(sim_operation):
	def __init__(self,coef = 1.0):
		self.coef = coef
	def get_params(self):
		return ['iteration','bx','by','bz']
	def get_ret_dim(self):
		return 1
	def __call__(self,iteration,bx,by,bz):

		lz = self.sim.ddx(bx*bz)
		lz += self.sim.ddy(by*bz)
		lz += self.sim.ddz(bz*bz)

		nz = self.sim.nz

		return [np.average(np.average(lz,1),1)]
		
"""Next time you want this, edit spec_gif to handle imhd"""		
class Spec_2D_gif(sim_operation):

	def get_params(self):
		return ["u","v", "w", "bx", "by", "bz"]

	def __call__(self,u,v,w,bx,by,bz):
		q = 4

		u = u[:,0,:]
		v = v[:,0,:]
		w = w[:,0,:]
		bx = bx[:,0,:]
		by = by[:,0,:]
		bz = bz[:,0,:]

		d = (self.sim.xmz, self.sim.xmx)
		vspec = spec(u,v,w,dims=d)
		bspec = spec(bx,by,bz,dims=d)

		plt.subplot(221)
		plt.plot(vspec)
		plt.semilogy()
		plt.title('v-spec')
		plt.axis("tight")
		plt.ylim([1e-12,1e-2])

		plt.subplot(222)
		plt.plot(vspec)
		plt.loglog()
		plt.title('v-spec')
		plt.axis("tight")
		plt.ylim([1e-12,1e-2])

		plt.subplot(223)
		plt.plot(bspec)
		plt.semilogy()
		plt.title('b-spec')
		plt.axis("tight")
		plt.ylim([1e-12,1e-2])

		plt.subplot(224)
		plt.plot(bspec)
		plt.loglog()
		plt.title('b-spec')
		plt.axis("tight")
		plt.ylim([1e-12,1e-2])

"""Next time you want this, edit spec_gif to handle imhd"""		
class Spec_3D_gif(sim_operation):

	def get_params(self):
		return ["u","v", "w", "bx", "by", "bz"]

	def __call__(self,u,v,w,bx,by,bz):
		q = 4

		d = (self.sim.zmx, self.sim.ymx, self.sim.xmx)
		vspec = spec(u,v,w,dims=d)
		bspec = spec(bx,by,bz,dims=d)

		plt.subplot(221)
		plt.plot(vspec)
		plt.semilogy()
		plt.title('v-spec')
		plt.axis("tight")
		plt.ylim([1e-12,1e-2])

		plt.subplot(222)
		plt.plot(vspec)
		plt.loglog()
		plt.title('v-spec')
		plt.axis("tight")
		plt.ylim([1e-12,1e-2])

		plt.subplot(223)
		plt.plot(bspec)
		plt.semilogy()
		plt.title('b-spec')
		plt.axis("tight")
		plt.ylim([1e-12,1e-2])

		plt.subplot(224)
		plt.plot(bspec)
		plt.loglog()
		plt.title('b-spec')
		plt.axis("tight")
		plt.ylim([1e-12,1e-2])

class hyperbolic_leak(sim_operation):
	def __init__(self,focused=True):
		self.suffix = "LeakStructs"
		self.focused = focused
	def get_params(self):
		return ["iteration","by"]

	def __call__(self,iteration,by):

		try:
			unstable = self.sim.get_file(iteration, "unstable", (self.sim.nz, self.sim.nx))
		except:
			print "filling unstable " + iteration + " with zeros"
			unstable = np.zeros((self.sim.nz,self.sim.nx),self.sim.dtype)
		try:
			stable = self.sim.get_file(iteration, "stable", (self.sim.nz, self.sim.nx))
		except:
			stable = np.zeros((self.sim.nz,self.sim.nx),self.sim.dtype)
			print "filling stable " + iteration + " with zeros"
		try:
			struct = self.sim.get_file(iteration, "vortexLeak", (self.sim.nz, self.sim.nx))
		except:
			struct = np.zeros((self.sim.nz,self.sim.nx),self.sim.dtype)
			print "filling struct " + iteration + " with zeros"

		hyper = unstable + stable

		if(self.focused == True):
			centerz  = np.unravel_index(np.argmax(by),by.shape)[0]
			width = by.shape[2];
			if(centerz < width+1):
				hyper = hyper[0:2*width,:]
				struct = struct[0:2*width,:]
			else:
				hyper = hyper[centerz-width:centerz+width,:]
				struct = struct[centerz-width:centerz+width,:]

		plt.subplot(121)
		plt.imshow(hyper,origin=0)
		plt.colorbar()
		plt.title('Hyperbolic Structures')
		plt.axis("tight")

		plt.subplot(122)
		plt.imshow(struct,origin=0,interpolation='nearest')
		plt.colorbar()
		plt.title('Leak Regions')
		plt.axis("tight")

class hyperbolic_mask(sim_operation):
	def get_params(self):
		return ["iteration"]

	def get_ret_dim(self):
		return 1

	def __call__(self,iteration):

		try:
			unstable = self.sim.get_file(iteration, "unstable")
		except:
			print "filling unstable " + iteration + " with zeros"
			unstable = np.zeros((self.sim.nz,1,self.sim.nx),self.sim.dtype)
		try:
			stable = self.sim.get_file(iteration, "stable")
		except:
			stable = np.zeros((self.sim.nz,1,self.sim.nx),self.sim.dtype)
			print "filling stable " + iteration + " with zeros"

		hyper = np.log(unstable + stable)

		plt.figure(1)
		plt.clf()
		plt.imshow(hyper[:,0,:],origin=0)
		plt.colorbar()
		plt.title('Hyperbolic Structures')
		plt.axis("tight")

		print "Select Bounding Box"
		trash = raw_input()

		xmin,xmax = plt.xlim()  
		zmin,zmax = plt.ylim()  

		print xmin,xmax
		print zmin,zmax

		mask = np.zeros_like(hyper)
		mask[zmin:zmax,:,xmin:xmax] = 1
		self.sim.write_file(iteration,"mask",mask)

		return 0
	
class hyperbolic(sim_operation):

	def get_params(self):
		return ["iteration"]

	def __call__(self,iteration):

		try:
			unstable = self.sim.get_file(iteration, "unstable", (self.sim.nz, self.sim.nx))
		except:
			print "filling unstable " + iteration + " with zeros"
			unstable = np.zeros((self.sim.nz,self.sim.nx),self.sim.dtype)
		try:
			stable = self.sim.get_file(iteration, "stable", (self.sim.nz, self.sim.nx))
		except:
			stable = np.zeros((self.sim.nz,self.sim.nx),self.sim.dtype)
			print "filling stable " + iteration + " with zeros"

		plt.subplot(121)
		plt.imshow(unstable,origin=0)
		plt.colorbar()
		plt.title('unstable')
		plt.axis("tight")

		plt.subplot(122)
		plt.imshow(stable,origin=0)
		plt.colorbar()
		plt.title('stable')
		plt.axis("tight")

class by_flux(sim_operation): 
	def __init__(self,dc):
		self.dc = dc
	def get_params(self):
		return ['w','by']
	def get_ret_dim(self):
		return 3
	def __call__(self,w,by):

		mean = np.average(np.average(by,1),1)
		adv = np.average(np.average(by*w,1),1)
		diff = -np.average(np.average(self.sim.ddz(by),1),1)*self.dc

		return [mean, adv, diff]

class b2_flux(sim_operation): 
	def get_params(self):
		return ['iteration','w','b2']
	def get_ret_dim(self):
		return 1
	def __call__(self,iteration,w,b2):

		adv = np.average(np.average(b2*w,1),1)

		return [adv]

class Spec_gif(sim_operation):
	def __init__(self, heights, mag=False, loglog=False):
		self.mag = mag
		self.heights = heights
		self.loglog = loglog

	def get_params(self):
		return ["iteration"]

	def __call__(self, iteration):
		self.sim.plot_spec(int(iteration),mag=self.mag,loglog=self.loglog,heights=self.heights)


class center_mass(sim_operation): 
	def __init__(self,field='b2'):
		self.field = field
	def get_params(self):
		return ['iteration',self.field]
	def get_ret_dim(self):
		return 3
	def __call__(self,iteration,field):
		
		field = np.average(field,2)

		(nz,ny) = field.shape
		
		cz = 0
		cy = 0
		for i in range(nz):
			cz += (field[i,:]*i).sum()
		for i in range(ny):
			cy += (field[:,i]*i).sum()
		total = field.sum()
		return [self.sim.get_time(iteration),cz/total, cy/total]

class filters_gif(sim_operation):

	def __init__(self,heights):
		self.heights = heights

	def get_params(self):
		return ["u","v","w","r","b2"]

	def __call__(self,u,v,w,r,b2):

		eng = r*(u*u + v*v + w*w)

		plt.subplot(211)
		for i in self.heights:
			layer = eng[i,:,:]
			vals = filter_conv(layer)
			plt.plot(vals,label="z = " + str(i))
		plt.title('Kinetic')
		plt.semilogy()
		plt.legend()

		plt.subplot(212)
		for i in self.heights:
			layer = b2[i,:,:]
			vals = filter_conv(layer)
			plt.plot(vals,label="z = " + str(i))
		plt.title('Magnetic')
		plt.xlabel('k')
		plt.semilogy()
		plt.legend()


class filters2_gif(sim_operation):

	def __init__(self,heights):
		self.heights = heights

	def get_params(self):
		return ["u","v","w","r","bx","by","bz"]

	def __call__(self,u,v,w,r,bx,by,bz):

		freq = self.sim.nx/2

		plt.subplot(211)
		for i in self.heights:
			vals = filter_conv(u[i,:,:],v[i,:,:],w[i,:,:],r[i,:,:])
			plt.plot(vals,label="z = " + str(i))
		plt.title('Kinetic')
		plt.semilogy()
		plt.legend()

		plt.subplot(212)
		for i in self.heights:
			vals = filter_conv(bx[i,:,:],by[i,:,:],bz[i,:,:])
			plt.plot(vals,label="z = " + str(i))
		plt.title('Magnetic')
		plt.xlabel('k')
		plt.semilogy()
		plt.legend()

class tube_field_gif(sim_operation):
	def __init__(self, focused=True):
		self.suffix="fieldlines"
		self.focused = focused

	def get_params(self):
		return ["bz","b2"]

	def __call__(self,bz,b2):

		phi = np.zeros_like(bz)
		dx = self.sim.xmx / self.sim.nx

		for i in range(1,self.sim.nx-1):
			phi[:,:,i] = phi[:,:,i-2] - 2 * dx * bz[:,:,i-1]

		if(self.focused == True):
			centerz  = np.unravel_index(np.argmax(b2),b2.shape)[0]
			width = phi.shape[2];
			if(centerz < width/2+1):
				phi = phi[0:width,:,:]
			else:
				phi = phi[centerz-width/2:centerz+width/2,:,:]

		phi = phi - phi.min()
		m = np.log(np.max(phi))
		levels = np.arange(m,m-3,-3./30)
		plt.contour(np.log(phi[:,0,:]),levels)

		#plt.imshow(phi[:,0,:],origin=0)
		#plt.colorbar()
		
		plt.axis('tight')

class tube_lorentz_gif(sim_operation):
	def __init__(self, t,focused=True):
		self.suffix="fieldlines"
		self.focused = focused
		self.t = tDyn

	def get_params(self):
		return ["bx","by","bz"]

	def __call__(self,bx,by,bz):

		dx = self.sim.xmx / self.sim.nx

		for i in range(1,self.sim.nx-1):
			phi[:,:,i] = phi[:,:,i-2] - 2 * dx * bz[:,:,i-1]

		if(self.focused == True):
			centerz  = np.unravel_index(np.argmax(b2),b2.shape)[0]
			width = phi.shape[2];
			if(centerz < width/2+1):
				phi = phi[0:width,:,:]
			else:
				phi = phi[centerz-width/2:centerz+width/2,:,:]

		phi = phi - phi.min()
		m = np.log(np.max(phi))
		levels = np.arange(m,m-3,-3./30)
		plt.contour(np.log(phi[:,0,:]),levels)

		#plt.imshow(phi[:,0,:],origin=0)
		#plt.colorbar()
		
		plt.axis('tight')

class var_gif(sim_operation):

	def __init__(self,variable_name,x=10):
		self.name = variable_name
		self.x = x

	def get_params(self):
		return [self.name]

	def __call__(self,var):

		plt.imshow(var[:,:,self.x])
		plt.colorbar()
		plt.axis('tight')

class vortDynamics(sim_operation): 
	def __init__(self,tDyn):
		self.t = tDyn
	def get_params(self):
		return ['iteration']
	def get_ret_dim(self):
		return 1
	def __call__(self,iteration):
		
		time = self.sim.get_time(iteration)
		v = self.t.calcTrajectory(iteration, xr=[0,self.sim.nx/2])

		return [time, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]]

class vortDynamics2(sim_operation): 
	def __init__(self,tDyn,type):
		self.t = tDyn
		self.type = type
	def get_params(self):
		return ['iteration']
	def get_ret_dim(self):
		return 1
	def __call__(self,iteration):
		
		vels = self.t.calcTrajectory(iteration, xr=[0,self.sim.nx/2])

		return vels[self.type]

class pressureGif(sim_operation): 
	def __init__(self,tDyn):
		self.t = tDyn
		self.suffix="pressure"
	def get_params(self):
		return ['iteration']
	def get_ret_dim(self):
		return 1
	def __call__(self,iteration):
		
		p = self.t.calcPress(iteration,pOnly=True)

		plt.imshow(p[:,0,:],origin=0)
		plt.axis('tight')
		plt.colorbar()

class twistStrength(sim_operation): 
	def __init__(self,focused=True):
		self.suffix="twist_strength"
		self.focused = focused
	def get_params(self):
		return ['bx','by','bz']
	def __call__(self,bx,by,bz):
		
		bt = np.sqrt(bx*bx + bz*bz)

		if(self.focused == True):
			centerz  = np.unravel_index(np.argmax(by),by.shape)[0]
			width = bt.shape[2];
			if(centerz < width/2+1):
				bt = bt[0:width,:,:]
			else:
				bt = bt[centerz-width/2:centerz+width/2,:,:]

		plt.imshow((bt/by.max())[:,0,:],origin=0)
		plt.colorbar()
		plt.axis('tight')

class calc_leaks(sim_operation): 
	def get_params(self):
		return ['iteration','t','by']
	def __call__(self,iteration,t,by):
		try:
			struct = self.sim.get_file(iteration, "vortexLeak")
		except:
			print "filling structure masks " + iteration + " with zeros"
			struct = np.zeros_like(t)

		dv = self.sim.xmx*self.sim.zmx/self.sim.nx/self.sim.nz

		insideBoth = np.sum(t[struct > 2.5])*dv
		insideNew = np.sum(t[np.logical_and(struct > 1.5, struct < 2.5)])*dv
		insideOld = np.sum(t[np.logical_and(struct > .5,struct < 1.5)])*dv
		insideMag = np.sum(t[by > .1*by.max()])*dv
		total = np.sum(t)*dv

		time = self.sim.get_time(iteration)
		print insideBoth + insideNew + insideOld, np.sum(t[struct > .5])*dv
		ret = [time,(insideBoth + insideNew)/total, insideNew/total, insideOld/total, insideMag/total]
		return ret
		
	def get_ret_dim(self):
		return 5

class bScales(sim_operation):
	def __init__(self, wave):
		self.wave = wave

	def get_params(self):
		return ['iteration','u','v','w','bx','by','bz']

	def get_ret_dim(self):
		return 9

	def __call__(self,iteration, u, v, w, bx, by, bz):

		uf = filter_single(self.wave,u,calcError=0)
		vf = filter_single(self.wave,v,calcError=0)
		wf = filter_single(self.wave,w,calcError=0)
		bxf = filter_single(self.wave,bx,calcError=0)
		byf = filter_single(self.wave,by,calcError=0)
		bzf = filter_single(self.wave,bz,calcError=0)

		uf[1] = u - uf[0]
		vf[1] = v - vf[0]
		wf[1] = w - wf[0]
		bxf[1] = bx - bxf[0]
		byf[1] = by - byf[0]
		bzf[1] = bz - bzf[0]

		lldbx,lldby,lldbz = self.__dB__(bxf[0],byf[0],bzf[0],uf[0],vf[0],wf[0])
		lsdbx,lsdby,lsdbz = self.__dB__(bxf[0],byf[0],bzf[0],uf[1],vf[1],wf[1])
		sldbx,sldby,sldbz = self.__dB__(bxf[1],byf[1],bzf[1],uf[0],vf[0],wf[0])
		ssdbx,ssdby,ssdbz = self.__dB__(bxf[1],byf[1],bzf[1],uf[1],vf[1],wf[1])
		fdbx,fdby,fdbz = self.__dB__(bx,by,bz,u,v,w)

		lll = (lldbx*bxf[0] + lldby*byf[0] + lldbz*bzf[0]).mean()
		lsl = (lsdbx*bxf[0] + lsdby*byf[0] + lsdbz*bzf[0]).mean()
		sll = (sldbx*bxf[0] + sldby*byf[0] + sldbz*bzf[0]).mean()
		ssl = (ssdbx*bxf[0] + ssdby*byf[0] + ssdbz*bzf[0]).mean()
		lls = (lldbx*bxf[1] + lldby*byf[1] + lldbz*bzf[1]).mean()
		lss = (lsdbx*bxf[1] + lsdby*byf[1] + lsdbz*bzf[1]).mean()
		sls = (sldbx*bxf[1] + sldby*byf[1] + sldbz*bzf[1]).mean()
		sss = (ssdbx*bxf[1] + ssdby*byf[1] + ssdbz*bzf[1]).mean()

		full = (fdbx*bx + fdby*by + fdbz*bz).mean()
	#	print (full - lll -lsl - sll - ssl - lls - lss - sls - sss) / full 
#		print np.max([uf[2],vf[2],wf[2],bxf[2],byf[2],bzf[2]])
		ret = [self.sim.get_time(iteration),lll,lsl,sll,ssl,lls,lss,sls,sss]
#		print ret
		return ret
	
	def __dB__(self,bx,by,bz,u,v,w):
		ix,iy,iz = self.sim.curl(*cross(u,v,w,bx,by,bz))

		return ix,iy,iz

