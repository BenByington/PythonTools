import inspect
from load_data import *
from bashcmd import *
import os
import sys, traceback
from spec import spec
from spec import spec_comp

import matplotlib.pyplot as plt

class HPS:
	_ratt = ["directory","nx","ny","nz", "xmx", "ymx", "zmx","scalar_f",\
		     "flux_f"]

	directory = None

	def __init__(self, loc='.', nx=None, ny=None, nz=None, \
		  xmx = None, ymx=None, zmx = None,\
				  scalar_f = None, flux_f = None):
		failed = ValueError("Could not create HPS object")
		try:
			curr = cd(loc,save=True)
		except OSError:
			raise failed
		self.directory = pwd()
		if nx != None:
			self.nx = int(nx)
		if ny != None:
			self.ny = int(ny)
		if nz != None:
			self.nz = int(nz)
		if xmx != None:
			self.xmx = float(xmx)
		if ymx != None:
			self.ymx = float(ymx)
		if zmx != None:
			self.zmx = float(zmx)
		if scalar_f != None:
			self.scalar_f = int(scalar_f)
		if flux_f != None:
			self.flux_f = int(flux_f)
		self.parse_start("start.in")
		cd(curr)
		
		for a in self._ratt:
			if hasattr(self,a) == False:
				cd(curr)
				print "Attribute " + a + " does not exist"
				raise failed

	def parse_start(self, file):
		try:
			file = open(file)
			while 1:
				line = file.readline()
				if not line:
					break
				words = line.replace('=',' ').split()
				if len(words) == 0:
					continue
				if words[0] == 'nx' and not hasattr(self,'nx'):
					self.nx = int(words[1])
				if words[0] == 'ny' and not hasattr(self,'ny'):
					self.ny = int(words[1])
				if words[0] == 'nz' and not hasattr(self,'nz'):
					self.nz = int(words[1])
				if words[0] == 'x_max' and not hasattr(self,'xmx'):
					self.xmx = float(words[1])
				if words[0] == 'y_max' and not hasattr(self,'ymx'):
					self.ymx = float(words[1])
				if words[0] == 'height' and not hasattr(self,'zmx'):
					self.zmx = float(words[1])
				if words[0] == 'scalar%frequency' and not hasattr(self,'scalar_f'):
					self.scalar_f = int(words[1])
				if words[0] == 'average%frequency' and not hasattr(self,'flux_f'):
					self.flux_f = int(words[1])
		except IOError:
			print "Warning, could not parse: File " + file + " does not exist in " + pwd()

	def load_scalars(self, nstart=None,nstop=None,nstep=None, \
					 nrec=None,list=None):
		curr = cd(self.directory + "/Scalar",save=True)
		self.scalar_data = make_scalar_accum(nstart=nstart,nstop=nstop,\
											 nstep=nstep,\
										 nfreq=self.scalar_f,nrec=nrec, \
											 n_scalars=40,list=list)
		cd(curr)
	
	def load_flux(self, nstart=None,nstop=None,nstep=None, \
				  nrec=None,list=None):
		curr = cd(self.directory + "/Flux",save=True)
		self.flux_data = make_flux_accum(nstart=nstart,nstop=nstop,\
											 nstep=nstep,\
										 nfreq=self.flux_f,nrec=nrec, \
											 nz = self.nz, \
											 n_fluxes=40,list=list)
		cd(curr)

	def clear_flux(self):
		if hasattr(self,'flux_data'):
			del self.flux_data
	
	def clear_scalars(self):
		if hasattr(self,'scalar_data'):
			del self.scalar_data
	
	def scalars(self,var=None):
		if hasattr(self,'scalar_data') == False:
			self.load_scalars()
		if var == None:
			return self.scalar_data
		else:
			return self.scalar_data[:,var]

	def flux(self):
		if hasattr(self,'flux_data') == False:
			self.load_flux()
		return self.flux_data

	def get_b(self, itr):
		try:
			if "start" in str(itr).lower():
				curr = cd(self.directory + "/Start/",save=True)
			elif "end" in str(itr).lower() or "last" in str(itr).lower():
				dir = np.array(map(int,get_data_files(path=self.directory + '/Bulk/'))).max()
				print "last snapshot is iteration ",dir
				curr = cd(self.directory + "/Bulk/" + str(dir).zfill(6),save=True)
			else:
				curr = cd(self.directory + "/Bulk/" + str(itr).zfill(6),save=True)
			ret = get_b(nx=self.nx,ny=self.ny,nz=self.nz)
			cd(curr)
			return ret
		except OSError:
			print "Bulk record " + str(itr) + " does not exist"
			
	def get_v(self, itr):
		try:
			if "start" in str(itr).lower():
				curr = cd(self.directory + "/Start/",save=True)
			elif "end" in str(itr).lower() or "last" in str(itr).lower():
				dir = np.array(map(int,get_data_files(path=self.directory + '/Bulk/'))).max()
				print "last snapshot is iteration ",dir
				curr = cd(self.directory + "/Bulk/" + str(dir).zfill(6),save=True)
			else:
				curr = cd(self.directory + "/Bulk/" + str(itr).zfill(6),save=True)
			ret = get_v(nx=self.nx,ny=self.ny,nz=self.nz)
			cd(curr)
			return ret
		except OSError:
			print "Bulk record " + str(itr) + " does not exist"
			
	def get_time(self, itr):
		try:
			itr = int(itr)
		except ValueError:
			print itr + "is not a valid integer"
			raise
		try:
			if "start" in str(itr).lower():
				curr = cd(self.directory + "/Start/",save=True)
			elif "end" in str(itr).lower() or "last" in str(itr).lower():
				dir = np.array(map(int,get_data_files(path=self.directory + '/Bulk/'))).max()
				print "last snapshot is iteration ",dir
				curr = cd(self.directory + "/Bulk/" + str(dir).zfill(6),save=True)
			else:
				curr = cd(self.directory + "/Bulk/" + str(itr).zfill(6),save=True)
			file = open('bulk_info')
			for i in range(7):
				info_str = file.readline()
			cd(curr)
			return float(info_str[7:-2])
		except OSError:
			print "Info record " + str(itr) + " does not exist"
			cd(curr)
		except ValueError:
			print "Could not find time for record " + str(itr)
			return 0
			

	def bulk_itr(self,vars=[],starti=None,endi=None,bstep=None):
		
		validV = {'bx' : 'Bx', 'by' : 'By', 'bz' : 'Bz', 'b2' : 'B2', \
		  'u' : 'u', 'v' : 'v', 'w' : 'w', \
		   't' : 'temperature', 'r' : 'density'}
		shape = (self.nz, self.ny, self.nx)
		if 'b2' in vars:
			allB = True
		else:
			allB = False

		if(starti == None or endi == None or bstep == None):
			curr = cd(self.directory + "/Bulk",save=True)
			list = get_data_files()
			cd(curr)
		else:
			list = range(starti+bstep,endi+bstep,bstep)
			for i in range(len(list)):
				list[i] = str(list[i]).zfill(6)

		for f in list:
			print "Grabbing Bulk " + f + " out of " + list[-1]
			try:
				curr = cd(self.directory + "/Bulk/" + f,save=True)
				if allB == True:
					bx,by,bz,b2 = get_b(nx=self.nx,ny=self.ny,nz=self.nz)
				
				ret = []
				for v in vars:
					if allB == True and v is 'bx':
						ret.append(bx)
					elif allB == True and v is 'by':
						ret.append(by)
					elif allB == True and v is 'bz':
						ret.append(bz)
					elif allB == True and v is 'b2':
						ret.append(b2)
					elif v is 'iteration':
						ret.append(f.zfill(len(list[-1])))
					else:
						ret.append(get_rec(validV[v],shape,np.float32))

				cd(curr)
				yield ret
			except OSError:
				print str(OSError)
				print "Trying to iterate over invalid Bulk files. "\
				+ f + " does not exist!  ...Skipping..."
			except IOError:
				print "File " + validV[v] + " in record " + f + " does not exist!  Skipping this snapshot..."
			except KeyError:
				print "Invalid argument list " + str(vars) + " .  Valid arugments are " + str(validV) + " and iteration"
				raise
				

	def make_gif(self, draw, starti=None, endi=None, bstep=None, dir=None, 
				 size=(8,6), suffix=None):
		if dir == None:
			dir = self.directory

		words = dir.split("/")
		if suffix == None:
			suffix = words[-1]
		else:
			suffix = words[-1] + "_" + suffix

		args = draw.get_params()
		args.append('iteration')
	
		draw.sim = self

		try:
			curr = cd(dir,save=True)
		except OSError:
			curr = cd(self.directory,save=True)
			print "Warning: " + dir + " does not exist.  Working in directory "\
			+ self.directory

		interactive = plt.isinteractive()
		if interactive == True:
			plt.ioff()

		figure = plt.figure(figsize=size)
		figure.rastorize = True
		ax = plt.gca()
		ax.ticklabel_format(style='sci',axis='x')
		ax.ticklabel_format(style='sci',axis='y')
		ax.xaxis.major.formatter.set_powerlimits((-3,4))
		ax.yaxis.major.formatter.set_powerlimits((-3,4))

		try:
			for b in self.bulk_itr(vars=args,starti=starti, \
							   endi=endi, bstep=bstep):
				plt.clf()
				draw(*b[0:-1])
				figure.savefig(b[-1] + suffix + ".png")
		
			os.system('tar -cvf ' + suffix + '.tar *' + suffix + '.png')
			os.system("convert -verbose *" + suffix + ".png " + suffix + ".gif")
			os.system("rm *" + suffix + ".png")
					
	
		except ValueError:
			print "Invalid Bulk folders. Aborting animation"
			print '-'*60
			traceback.print_exc(file=sys.stdout)
			print '-'*60
		
		if interactive == True:
			plt.ion()

		cd(curr)

	def show_frame(self, draw, iteration, size=(8,6)):

		args = draw.get_params()
		args.append('iteration')
	
		draw.sim = self

		figure = plt.figure(figsize=size)
		figure.rastorize = True
		ax = plt.gca()
		ax.ticklabel_format(style='sci',axis='x')
		ax.ticklabel_format(style='sci',axis='y')
		ax.xaxis.major.formatter.set_powerlimits((-3,4))
		ax.yaxis.major.formatter.set_powerlimits((-3,4))

		try:
			for b in self.bulk_itr(vars=args,starti=int(iteration)-1, \
							   endi=int(iteration), bstep=1):
				plt.clf()
				draw(*b[0:-1])
		
	
		except ValueError:
			print "Invalid Bulk folder. Aborting animation"
			print '-'*60
			traceback.print_exc(file=sys.stdout)
			print '-'*60
		
	def bulk_calc(self, calc, starti=None, endi=None, bstep=None):

		args = calc.get_params()
		ret = []
		s = calc.get_ret_dim()
		for i in range(s):
			ret.append([])
	
		calc.sim = self

		try:
			for b in self.bulk_itr(vars=args,starti=starti, \
								   endi=endi, bstep=bstep):
				temp = calc(*b)
				if s > 1:
					for i in range(s):
						ret[i].append(temp[i])
				else:
					ret[0].append(temp)

			if s > 1:
				return np.array(ret,np.float32)
			else:
				return np.array(ret[0],np.float32)
	
		except ValueError:
			print "Invalid Bulk folders. Aborting calculation"
			print '-'*60
			traceback.print_exc(file=sys.stdout)
			print '-'*60
		

	def ddx(self, data,xmx=None):
		nx = data.shape[2]
		if xmx == None:
			xmx = self.xmx

		dx = xmx / float(nx)

		return (np.roll(data,-1,2) - np.roll(data,1,2))/2/dx

	def ddy(self, data,ymx=None):
		ny = data.shape[1]
		if ymx == None:
			ymx = self.ymx

		dy = ymx / float(ny)

		return (np.roll(data,-1,1) - np.roll(data,1,1))/2/dy

	def ddz(self, data,zmx=None):
		nz = data.shape[0]
		if zmx == None:
			zmx = self.zmx

		dz = zmx / float(nz-1)

		ret =  (np.roll(data,-1,0) - np.roll(data,1,0))/2/dz
		ret[0,:,:] = (-3*data[0,:,:] + 4*data[1,:,:] - data[2,:,:])/2/dz
		ret[-1,:,:] = (3*data[-1,:,:] - 4*data[-2,:,:] + data[-3,:,:])/2/dz

		return ret

	def curl(self, x, y, z, xmx=None, ymx=None, zmx=None):
		if xmx == None:
			xmx = self.xmx
		if ymx == None:
			ymx = self.ymx
		if zmx == None:
			zmx = self.zmx

		retx = self.ddy(z,ymx=ymx) - self.ddz(y,zmx=zmx)
		rety = self.ddz(x,zmx=zmx) - self.ddx(z,xmx=xmx)
		retz = self.ddx(y,xmx=xmx) - self.ddy(x,ymx=ymx)

		return [retx, rety, retz]

	def diffusion(self, x, y, z, xmx=None, ymx=None, zmx=None, coef=1):
		if xmx == None:
			xmx = self.xmx
		if ymx == None:
			ymx = self.ymx
		if zmx == None:
			zmx = self.zmx

		retx = self.ddx(self.ddx(x,xmx=xmx),xmx=xmx)
		retx = retx + self.ddy(self.ddy(x,ymx=ymx),ymx=ymx)
		retx = retx + self.ddz(self.ddz(x,zmx=zmx),zmx=zmx)

		rety = self.ddx(self.ddx(y,xmx=xmx),xmx=xmx)
		rety = rety + self.ddy(self.ddy(y,ymx=ymx),ymx=ymx)
		rety = rety + self.ddz(self.ddz(y,zmx=zmx),zmx=zmx)

		retz = self.ddx(self.ddx(z,xmx=xmx),xmx=xmx)
		retz = retz + self.ddy(self.ddy(z,ymx=ymx),ymx=ymx)
		retz = retz + self.ddz(self.ddz(z,zmx=zmx),zmx=zmx)

		return coef * np.array([retx, rety, retz])

	def spec(self,iteration, mag=False,heights=0):
		if type(heights) is not list: heights = [heights]

		if (mag):
			a,b,c,trash = self.get_b(iteration)
		else:
			a,b,c,d,trash = self.get_v(iteration)

		ret = []
		dims = (self.ymx,self.xmx)

		if(mag):
			for i in heights:
				ret.append(spec(a[i,:,:],b[i,:,:],c[i,:,:],dims=dims))
		else:
			for i in heights:
				ret.append(spec_comp(a[i,:,:],b[i,:,:],c[i,:,:],d[i,:,:],dims=dims))

		return ret

	def plot_spec(self,iteration,mag=False,loglog=False,heights=0):
		if type(heights) is not list: heights = [heights]
	
		if(mag):
			title = "Magnetic"
		else:
			title = "Momentum"

		title = title + " Power Spectrum: Iteration " + str(iteration)

		powers = self.spec(iteration,mag=mag,heights=heights)
		plt.clf()
		for i in range(len(heights)):
			plt.plot(powers[i],label="z = " + str(heights[i]))
		if(loglog):
			plt.loglog()
		else:
			plt.semilogy()
		plt.legend()
		plt.title(title)
		plt.xlabel("z")
		return
		
