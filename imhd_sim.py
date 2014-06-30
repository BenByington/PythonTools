import inspect
from load_data import *
from bashcmd import *
import os
import sys, traceback
from spec import spec

try:
	import matplotlib.pyplot as plt
except:
	print "no matplotlib available on this machine.  Prepare to fail gracelessly if used!"

class IMHD:
	_ratt = ["directory","nx","ny","nz", "xmx", "ymx", "zmx","scalar_f",
	      "magnetic","temperature"]

	directory = None

	def __init__(self, loc='.', nx=None, ny=None, nz=None, \
		  xmx = None, ymx=None, zmx = None,\
				  scalar_f = None, dtype=np.float64, \
				  magnetic = None, temperature = None):

		self.dtype = dtype
		failed = ValueError("Could not create IMHD object")
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
		if magnetic != None:
			self.magnetic = magnetic
		if temperature != None:
			self.temperature == temperature
		self.parse_start("start.in")
		cd(curr)
		
		for a in self._ratt:
			if hasattr(self,a) == False:
				cd(curr)
				print "Attribute " + a + " does not exist"
				raise failed

		self.dkx = None
		self.dky = None
		self.dkz = None

		self.kxmin = (self.nx/2)*2/3+1
		self.kymin = (self.ny/2)*2/3+1
		self.kzmin = (self.nz/2)*2/3+1

		self.kxmax = self.nx - 1
		self.kymax = self.ny - self.kymin
		self.kzmax = self.nz - self.kzmin


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
				if words[0] == 'xmx' and not hasattr(self,'xmx'):
					self.xmx = float(words[1])
				if words[0] == 'ymx' and not hasattr(self,'ymx'):
					self.ymx = float(words[1])
				if words[0] == 'zmx' and not hasattr(self,'zmx'):
					self.zmx = float(words[1])
				if words[0] == 'scalarRate' and not hasattr(self,'scalar_f'):
					self.scalar_f = int(words[1])
				if words[0] == 'magneticEQ' and not hasattr(self,'magnetic'):
					if words[1] == 'on':
						self.magnetic = True
					else:
						self.magnetic = False
				if words[0] == 'temperatureEQ' and not hasattr(self,'temperature'):
					if words[1] == 'on':
						self.temperature = True
					else:
						self.temperature = False
		except IOError:
			print "Warning, could not parse: File " + file + " does not exist in " + pwd()

	def load_scalars(self, nstart=None,nstop=None,nstep=None, \
					 nrec=None,list=None):
		if self.magnetic:
			nscalars = 24
		else:
			nscalars = 13
		curr = cd(self.directory + "/Scalars",save=True)
		self.scalar_data = make_scalar_accum(nstart=nstart,nstop=nstop,\
											 nstep=nstep,\
										 nfreq=self.scalar_f,nrec=nrec, \
											 n_scalars=nscalars,list=list,\
											 type = self.dtype, \
											 zcount = 8)
		cd(curr)
	
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

	def get_b(self, itr):
		if self.magnetic == False:
			print "ERROR: This simulation has no magnetic fields!!"
		try:
			if "start" in str(itr).lower():
				curr = cd(self.directory + "/Start/",save=True)
			elif "end" in str(itr).lower() or "last" in str(itr).lower():
				dir = np.array(map(int,get_data_files(path=self.directory + '/Spatial/'))).max()
				print "last snapshot is iteration ",dir
				curr = cd(self.directory + "/Spatial/" + str(dir).zfill(8),save=True)
			else:
				curr = cd(self.directory + "/Spatial/" + str(itr).zfill(8),save=True)
			bx,by,bz,b2 = get_b(nx=self.nx,ny=self.ny,nz=self.nz,type=self.dtype)
			cd(curr)
			return bx,by,bz,b2
		except OSError:
			print "Spatial record " + str(itr) + " does not exist"
			
	def get_time(self, itr):
		try:
			if "start" in str(itr).lower():
				curr = cd(self.directory + "/Start/",save=True)
			elif "end" in str(itr).lower() or "last" in str(itr).lower():
				dir = np.array(map(int,get_data_files(path=self.directory + '/Spatial/'))).max()
				print "last snapshot is iteration ",dir
				curr = cd(self.directory + "/Spatial/" + str(dir).zfill(6),save=True)
			else:
				curr = cd(self.directory + "/Spatial/" + str(itr).zfill(8),save=True)
			info_str = open('info').read()
			cd(curr)
			return float(info_str[6:])
		except OSError:
			print "Info record " + str(itr) + " does not exist"
			cd(curr)
			
	def get_v(self, itr):
		try:
			if "start" in str(itr).lower():
				curr = cd(self.directory + "/Start/",save=True)
			elif "end" in str(itr).lower() or "last" in str(itr).lower():
				dir = np.array(map(int,get_data_files(path=self.directory + '/Spatial/'))).max()
				print "last snapshot is iteration ",dir
				curr = cd(self.directory + "/Spatial/" + str(dir).zfill(8),save=True)
			else:
				curr = cd(self.directory + "/Spatial/" + str(itr).zfill(8),save=True)
			ret = Iget_v(t=self.temperature,nx=self.nx,ny=self.ny,nz=self.nz,type=self.dtype)
			cd(curr)
			return ret
		except OSError:
			print "Spatial record " + str(itr) + " does not exist"
			cd(curr)
			
	def get_file(self, itr, name, shape=None):
		try:
			if "start" in str(itr).lower():
				curr = cd(self.directory + "/Start/",save=True)
			elif "end" in str(itr).lower() or "last" in str(itr).lower():
				dir = np.array(map(int,get_data_files(path=self.directory + '/Spatial/'))).max()
				print "last snapshot is iteration ",dir
				curr = cd(self.directory + "/Spatial/" + str(dir).zfill(8),save=True)
			else:
				curr = cd(self.directory + "/Spatial/" + str(itr).zfill(8),save=True)
			if(shape == None):
				shape = (self.nz, self.ny, self.nx)
			ret = get_rec(name,shape,self.dtype)

			cd(curr)
			return ret
		except OSError:
			print "Spatial record " + str(itr) + " does not exist or does not have file " + name
			cd(curr)
			raise
		except:
			cd(curr)
			raise
			

	def write_file(self, itr, name, data):
		try:
			if "start" in str(itr).lower():
				curr = cd(self.directory + "/Start/",save=True)
			else:
				curr = cd(self.directory + "/Spatial/" + str(itr).zfill(8),save=True)
			save_rec(name,data)

			cd(curr)
		except OSError:
			print "Spatial record " + str(itr) + " does not exist or does not have file " + name
			cd(curr)
			raise
		except:
			cd(curr)
			raise
			

	def bulk_itr(self,vars=[],starti=None,endi=None,bstep=None):
		validV = {'u' : 'u', 'v' : 'v', 'w' : 'w'}
		if self.magnetic == True:
			validV.update({'bx' : 'Bx', 'by' : 'By', 'bz' : 'Bz', 'b2' : 'B2'})
		if self.temperature == True:
			validV.update({'t' : 'T'})

		shape = (self.nz, self.ny, self.nx)
		if 'b2' in vars and self.magnetic == True:
			allB = True
		else:
			allB = False

		if(starti == None or endi == None or bstep == None):
			curr = cd(self.directory + "/Spatial",save=True)
			list = get_data_files(nz=8)
			cd(curr)
		else:
			list = range(starti+bstep,endi+bstep,bstep)
			for i in range(len(list)):
				list[i] = str(list[i]).zfill(8)

		for f in list:
			print "Grabbing Spatial " + f + " out of " + list[-1]
			try:
				curr = cd(self.directory + "/Spatial/" + f,save=True)
				if allB == True:
					bx,by,bz,b2 = get_b(nx=self.nx,ny=self.ny,nz=self.nz,type=self.dtype)
				
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
						ret.append(get_rec(validV[v],shape,self.dtype))

				cd(curr)
				yield ret
			except OSError:
				print str(OSError)
				print "Trying to iterate over invalid Bulk files. "\
				+ f + " does not exist! ...Skipping..."
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
		if suffix == None and draw.suffix == None:
			suffix = words[-1]
		elif suffix == None and draw.suffix != None:
			suffix = words[-1] + "_" + draw.suffix	
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
				return np.array(ret,self.dtype)
			else:
				return np.array(ret[0],self.dtype)
	
		except ValueError:
			print "Invalid Bulk folders. Aborting calculation"
			print '-'*60
			traceback.print_exc(file=sys.stdout)
			print '-'*60

	def dealias(self,data):
		data[self.kzmin:self.kzmax,:,:] = 0
		data[:,self.kymin:self.kymax,:] = 0
		data[:,:,self.kxmin:self.kxmax] = 0

	def assertDkx(self,data):
		if (self.dkx == None):
			self.dkx = np.ones_like(data)
			for x in range(0,self.kxmin+2,1):
				self.dkx[:,:,x] = 1j * 2 * np.pi * x / self.xmx

	def assertDky(self,data):
		if (self.dky == None):
			if (self.ny == 1):
				self.dky = np.zeros_like(data)
			else:
				self.dky = np.ones_like(data)
				for y in range(-self.kymin-2,self.kymin+2,1):
					self.dky[:,y,:] = 1j * 2 * np.pi * y / self.ymx

	def assertDkz(self,data):
		if (self.dkz == None):
			self.dkz = np.ones_like(data)
			for z in range(-self.kzmin-2,self.kzmin+2,1):
				self.dkz[z,:,:] = 1j * 2 * np.pi * z / self.zmx

	def makeSpect(self, a):
		ret = np.fft.rfftn(a)
		self.dealias(ret)
		return ret

	def mulSpect(self, a, b):
		if (np.isrealobj(a)):
			return a*b
		a = np.fft.irfftn(a)
		b = np.fft.irfftn(b)
		ret = a*b
		return self.makeSpect(ret)

	def ddx(self, data):

		transform = False
		if (np.isrealobj(data)):
			data = self.makeSpect(data)
			transform = True
		
		self.assertDkx(data)        

		ret = self.dkx * data

		if (transform == True):
			ret = np.fft.irfftn(ret)
    
 		return ret

	def ddy(self, data):

		transform = False
		if (np.isrealobj(data)):
			data = self.makeSpect(data)
			transform = True
        
		self.assertDky(data)

		ret = self.dky * data

		if (transform == True):
			ret = np.fft.irfftn(ret)
    
 		return ret

	def ddz(self, data):

		transform = False
		if (np.isrealobj(data)):
			data = self.makeSpect(data)
			transform = True
        
		self.assertDkz(data)

		ret = self.dkz * data

		if (transform == True):
			ret = np.fft.irfftn(ret)
    
 		return ret

	def inv_ddx(self, data):

		transform = False
		if (np.isrealobj(data)):
			data = self.makeSpect(data)
			transform = True
       
		self.assertDkx(data) 

		ret = np.copy(data)
		ret[:,:,1:] /= self.dkx[:,:,1:]
		ret[:,:,0] *= self.xmx

		if (transform == True):
			ret = np.fft.irfftn(ret)
    
 		return ret

	def inv_ddy(self, data):

		transform = False
		if (np.isrealobj(data)):
			data = self.makeSpect(data)
			transform = True
        
		self.assertDky(data)

		ret = np.copy(data)
		ret[self.dky != 0] /= self.dky[self.dky != 0]

		if (transform == True):
			ret = np.fft.irfftn(ret)
    
 		return ret

	def inv_ddz(self, data):

		transform = False
		if (np.isrealobj(data)):
			data = self.makeSpect(data)
			transform = True
        
		self.assertDkz(data)

		ret = np.copy(data)
		ret[self.dkz != 0] /= self.dkz[self.dkz != 0]

		if (transform == True):
			ret = np.fft.irfftn(ret)
    
 		return ret

	def fast_ddx(self, data,xmx=None):
		nx = data.shape[2]
		if xmx == None:
			xmx = self.xmx

		dx = xmx / float(nx)

		return (np.roll(data,-1,2) - np.roll(data,1,2))/2/dx

	def fast_ddy(self, data,ymx=None):
		ny = data.shape[1]
		if ymx == None:
			ymx = self.ymx

		dy = ymx / float(ny)

		return (np.roll(data,-1,1) - np.roll(data,1,1))/2/dy

	def fast_ddz(self, data,zmx=None):
		nz = data.shape[0]
		if zmx == None:
			zmx = self.zmx

		dz = zmx / float(nz)

		return (np.roll(data,-1,0) - np.roll(data,1,0))/2/dz

	def curl(self, x, y, z):

		transform = False
		if (np.isrealobj(x)):
			x = self.makeSpect(x)
			y = self.makeSpect(y)
			z = self.makeSpect(z)
			transform = True

		retx = self.ddy(z) - self.ddz(y)
		rety = self.ddz(x) - self.ddx(z)
		retz = self.ddx(y) - self.ddy(x)

		if(transform == True):
			retx = np.fft.irfftn(retx)
			rety = np.fft.irfftn(rety)
			retz = np.fft.irfftn(retz)

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

		return ceof * [retx, rety, retz]
	
	def spec(self,iteration, magnetic=False, cutDim=None, slices=0):
		if magnetic == True:
			temp = self.get_b(iteration)
		else:
			temp = self.get_v(iteration)

		a=temp[0]
		b=temp[1]
		c=temp[2]

		ret  = []
		if cutDim == 0:
			a = a[slices,:,:]
			b = b[slices,:,:]
			c = c[slices,:,:]
			dims = (self.ymx,self.xmx)
			for i in range(len(slices)):
				ret.append(spec(a[i,:,:],b[i,:,:],c[i,:,:],dims=dims))
		elif cutDim == 1:
			a = a[:,slices,:]
			b = b[:,slices,:]
			c = c[:,slices,:]
			dims = (self.zmx,self.xmx)
			for i in range(len(slices)):
				ret.append(spec(a[:,i,:],b[:,i,:],c[:,i,:],dims=dims))
		elif cutDim == 2:
			a = a[:,:,slices]
			b = b[:,:,slices]
			c = c[:,:,slices]
			dims = (self.zmx,self.ymx)
			for i in range(len(slices)):
				ret.append(spec(a[:,:,i],b[:,:,i],c[:,:,i],dims=dims))
		elif cutDim == None:
			dims = (self.zmx,self.ymx,self.xmx)
			return spec(a,b,c,dims)
		else:
			print "nonsensical cutDim! ",cutDim

		return ret


	#currently only handles full spectral decomposition.  No slices here!
	def plot_spec(self,iteration,mag=False,loglog=False,cutDim=None, slices=0):
		if(mag):
			title = "Magnetic"
		else:
			title = "Momentum"

		title = title + " Power Spectrum: Iteration " + str(iteration)

		powers = self.spec(iteration,magnetic=mag,cutDim=cutDim,slices=slices)
		plt.clf()
		if(cutDim == None):
			plt.plot(powers)
		else:
			if cutDim == 0:
				axis = "z"
			elif cutDim == 1:
				axis = "y"
			elif cutDim == 2:
				axis = "x"
			for i in range(len(slices)):
				plt.plot(powers[i],label=axis+" = " + str(slices[i]))
			plt.legend()
		if(loglog):
			plt.loglog()
		else:
			plt.semilogy()
		plt.title(title)
		return
