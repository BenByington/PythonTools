import numpy as np
try:
	import matplotlib.pyplot as plt
except:
	"matplotlib is not on this machine.  tDynamics will fali if drawDp is used!"
from physics import cross

class tDynamics:

	def __init__(self,sim, alpha, pr, mpr):
		self.sim = sim
		self.alpha = alpha
		self.pr = pr
		self.mpr = mpr
		self.magBuoy = True
		self.lorentz = True
		self.viscosity = True
		self.mDiff = True

	def drawDp(self,iteration,pres=None,dp=None):
		if (pres == None or dp == None):
			pres,dp,ddp = self.calcPress(iteration)

		plt.figure()
		plt.imshow(dp[:,0,:],origin=0)
		plt.title("dP/dt")
		plt.colorbar()

		plt.figure()
		plt.imshow(pres[:,0,:],origin=0)
		plt.title("P")
		plt.colorbar()

		plt.figure()
		plt.imshow(ddp[:,0,:],origin=0)
		plt.title("ddP/ddt")
		plt.colorbar()


	def calcTrajectory(self, iteration, pres=None, dp=None, ddp=None, xr=None,zr=None):
		if (pres == None or dp == None):
			pres,dp,ddp = self.calcPress(iteration)

		if (xr == None):
			xr = [0,dp.shape[2]]
		if (zr == None):
			zr = [0,dp.shape[0]]

		temp = pres[zr[0]:zr[1],:,xr[0]:xr[1]]
		ind = np.unravel_index(np.argmin(temp),temp.shape)
		ind = (ind[0] + zr[0], ind[1], ind[2] + xr[0])

		vel = calcPVel(self.sim,pres,dp,ddp,ind)

		return [ind[2]*self.sim.xmx/self.sim.nx, ind[0]*self.sim.zmx/self.sim.nz, vel[0], vel[1], vel[2], vel[3], vel[4], vel[5]]

	def calcPress(self,iteration,pOnly=False):
		sim = self.sim
		alpha = self.alpha
		pr = self.pr
		mpr = self.mpr

		vel = sim.get_v(iteration)
		u=vel[0]
		v=vel[1]
		w=vel[2]
		bx,by,bz,b2 = sim.get_b(iteration)

		dp = None

		temp = advectionAccel(sim,u,v,w)
		pres = temp[3]
		ax = temp[0]
		ay = temp[1]
		az = temp[2]
	
		if(pOnly == True):
			return pres

		if (self.viscosity == True):
			temp = diffAccel(sim,u,v,w,pr)
			pres = pres + temp[3]
			ax = ax + temp[0]
			ay = ay + temp[1]
			az = az + temp[2]

		if(self.magBuoy == True):
			temp = mBuoyAccel(sim,b2,alpha)
			pres = pres + temp[3]
			ax = ax + temp[0]
			ay = ay + temp[1]
			az = az + temp[2]

		if(self.lorentz == True):
			temp = lorentzAccel(sim,bx,by,bz,alpha)
			pres = pres + temp[3]
			ax = ax + temp[0]
			ay = ay + temp[1]
			az = az + temp[2]

		if (self.magBuoy == True or self.lorentz == True):
			temp = induct(sim,u,v,w,bx,by,bz)
			dbx = temp[0]
			dby = temp[1]
			dbz = temp[2]

			if (self.mDiff == True):
				temp = magDiff(sim,bx,by,bz,pr/mpr)
				dbx = dbx + temp[0]
				dby = dby + temp[1]
				dbz = dbz + temp[2]

			temp = dInduct(sim,u,v,w,ax,ay,az,bx,by,bz,dbx,dby,dbz)
			ddbx = temp[0]
			ddby = temp[1]
			ddbz = temp[2]
		
			if (self.mDiff == True):
				temp = magDiff(sim,dbx,dby,dbz,pr/mpr)
				ddbx = ddbx + temp[0]
				ddby = ddby + temp[1]
				ddbz = ddbz + temp[2]

		temp = dAdvectionAccel(sim,u,v,w,ax,ay,az)
		dax = temp[0]
		day = temp[1]
		daz = temp[2]
		dp =  temp[3]

		if (self.viscosity == True):
			temp = diffAccel(sim,ax,ay,az,pr)
			dax = dax + temp[0]
			day = day + temp[1]
			daz = daz + temp[2]
			dp = dp + temp[3]

		if (self.lorentz == True):
			temp = dLorentzAccel(sim,bx,by,bz,dbx,dby,dbz,alpha)
			dax = dax + temp[0]
			day = day + temp[1]
			daz = daz + temp[2]
			dp = dp + temp[3]

		if (self.magBuoy == True):
			temp = dMBuoyAccel(sim,bx,by,bz,dbx,dby,dbz,alpha)
			dax = dax + temp[0]
			day = day + temp[1]
			daz = daz + temp[2]
			dp = dp + temp[3]

		if (dp == None):
			dp = dpAdvect(sim,u,w,ax,az)
			if (self.magBuoy == True):
				dp = dp + dpBuoy(sim,bx,by,bz,dbx,dby,dbz,alpha)
			if (self.lorentz == True):
				dp = dp + dpLorentz(sim,bx,bz,dbx,dbz)
			
		ddp = ddpAdvect(sim,u,w,ax,az,dax,daz)

		if (self.magBuoy == True):
			ddp = ddp + ddpBuoy(sim, bx, by, bz, dbx, dby, dbz, ddbx, ddby, ddbz,alpha)
		if (self.lorentz == True):
			ddp = ddp + ddpLorentz(sim,bx,bz,dbx,dbz,ddbx,ddbz,alpha)

		return [pres,dp,ddp]

	def enablePureVort(self):
		self.magBuoy = False
		self.lorentz = False
		self.mDiff = False 
		self.viscosity = True
	def enableAllVort(self):
		self.magBuoy = True
		self.lorentz = True
		self.mDiff = True
		self.viscosity = True
	def enableUntwisted(self):
		self.magBuoy = True
		self.lorentz = False
		self.mDiff = True
		self.viscosity = True
	def setViscosity(bool):
		self.viscosity = bool
	def setLorentz(bool):
		self.lorentz = bool
	def setMagBuoy(bool):
		self.magBuoy = bool
	def setMDiff(bool):
		self.mdiff = bool

def dpBuoy(sim, bx, by, bz, dbx, dby, dbz, alpha):
	nl1 = sim.makeSpect(bx * dbx + by*dby + bz*dbz)
	div = (2*alpha)*sim.ddz(nl1)
	return np.fft.irfftn(inverseLaplace(sim,div))

def ddpBuoy(sim, bx, by, bz, dbx, dby, dbz, ddbx, ddby, ddbz,alpha):
	nl1 = sim.makeSpect(bx*ddbx + by*ddby + bz*ddbz + dbx**2 + dby**2 + dbz**2)
	div = (2*alpha)*sim.ddz(nl1)
	return np.fft.irfftn(inverseLaplace(sim,div))

def dpAdvect(sim, u, w, ax, az):
	nl1 = sim.makeSpect(u*ax)
	nl2 = sim.makeSpect(u*az + w*ax)
	nl3 = sim.makeSpect(w*az)
	div = sim.ddx(sim.ddx(nl1)) + sim.ddx(sim.ddz(nl2)) + sim.ddz(sim.ddz(nl3))
	return np.fft.irfftn(-2 * inverseLaplace(sim,div))

def ddpAdvect(sim, u, w, ax, az, dax, daz):
	nl1 = sim.makeSpect(ax**2 + u*dax)
	nl2 = sim.makeSpect(u*daz + 2*ax*az + w*dax)
	nl3 = sim.makeSpect(w*daz + az**2)
	div = sim.ddx(sim.ddx(nl1)) + sim.ddx(sim.ddz(nl2)) + sim.ddz(sim.ddz(nl3))
	return np.fft.irfftn(-2 * inverseLaplace(sim,div))

def dpLorentz(sim, bx, bz, dbx, dbz,alpha):
	nl1 = sim.makeSpect(bx*dbx)
	nl2 = sim.makeSpect(bx*dbz + bz*dbx)
	nl3 = sim.makeSpect(bz*dbz)
	div = sim.ddx(sim.ddx(nl1)) + sim.ddx(sim.ddz(nl2)) + sim.ddz(sim.ddz(nl3))
	return np.fft.irfftn(2 * alpha * inverseLaplace(sim,div))

def ddpLorentz(sim, bx, bz, dbx, dbz, ddbx, ddbz,alpha):
	nl1 = sim.makeSpect(bx * ddbx + dbx**2)
	nl2 = sim.makeSpect(bx*ddbz + 2*dbx*dbz + bz*ddbx)
	nl3 = sim.makeSpect(bz*ddbz + dbz**2)
	div = sim.ddx(sim.ddx(nl1)) + sim.ddx(sim.ddz(nl1)) + sim.ddz(sim.ddz(nl3))
	return np.fft.irfftn(2 * alpha * inverseLaplace(sim,div))

def advectionAccel(sim,u,v,w):
	ax = -sim.ddx(sim.makeSpect(u*u))
	az = -sim.ddz(sim.makeSpect(w*w))

	temp = sim.makeSpect(u*w)
	ax = ax - sim.ddz(temp)
	az = az - sim.ddx(temp)

	ay = sim.ddx(sim.makeSpect(u*v))
	ay = ay + sim.ddz(sim.makeSpect(w*v))

	ap = pressureResponse(sim,ax,az)

	retx = np.fft.irfftn(ax + ap[0])
	rety = np.fft.irfftn(ay)
	retz = np.fft.irfftn(az + ap[1])
	retp = np.fft.irfftn(ap[2])
	return [retx, rety, retz, retp]

def dAdvectionAccel(sim,u,v,w,du,dv,dw):
	ax = -2*sim.ddx(sim.makeSpect(u*du))
	az = -2*sim.ddz(sim.makeSpect(w*dw))

	temp = sim.makeSpect(u*dw + w*du)
	ax = ax - sim.ddz(temp)
	az = az - sim.ddx(temp)

	ay = -sim.ddx(sim.makeSpect(du*v + u*dv))
	ay = ay - sim.ddz(sim.makeSpect(dw*v + w*dv))

	ap = pressureResponse(sim,ax,az)

	retx = np.fft.irfftn(ax + ap[0])
	rety = np.fft.irfftn(ay)
	retz = np.fft.irfftn(az + ap[1])
	retp = np.fft.irfftn(ap[2])
	return [retx, rety, retz, retp]

def mBuoyAccel(sim,b2,alpha):
	az = sim.makeSpect(alpha*b2)
	ax = np.zeros_like(az)

	bp = pressureResponse(sim,ax,az)

	retx = np.fft.irfftn(ax + bp[0])
	retz = np.fft.irfftn(az + bp[1])
	retp = np.fft.irfftn(bp[2])
	return [retx, np.zeros_like(retx), retz, retp]

def dMBuoyAccel(sim,bx,by,bz,dbx,dby,dbz,alpha):
	az = 2*sim.makeSpect(alpha*(bx*dbx + by*dby + bz*dbz))
	ax = np.zeros_like(az)

	bp = pressureResponse(sim,ax,az)

	retx = np.fft.irfftn(ax + bp[0])
	retz = np.fft.irfftn(az + bp[1])
	retp = np.fft.irfftn(bp[2])
	return [retx, np.zeros_like(retx), retz, retp]

def lorentzAccel(sim,bx,by,bz,alpha):
	ax = sim.ddx(sim.makeSpect(bx*bx*alpha))
	az = sim.ddz(sim.makeSpect(bz*bz*alpha))

	temp = sim.makeSpect(bx*bz*alpha)
	ax = ax + sim.ddz(temp)
	az = az + sim.ddx(temp)

	ay = sim.ddx(sim.makeSpect(bx*by*alpha))
	ay = ay + sim.ddz(sim.makeSpect(bz*by*alpha))
	
	lp = pressureResponse(sim,ax,az)

	retx = np.fft.irfftn(ax + lp[0])
	rety = np.fft.irfftn(ay)
	retz = np.fft.irfftn(az + lp[1])
	retp = np.fft.irfftn(lp[2])
	return [retx, rety, retz, retp]

def dLorentzAccel(sim,bx,by,bz,dbx,dby,dbz,alpha):
	
	ax = 2*sim.ddx(sim.makeSpect(bx*dbx*alpha))
	az = 2*sim.ddz(sim.makeSpect(bz*dbz*alpha))

	temp = sim.makeSpect((bx*dbz + dbx*bz)*alpha)
	ax = ax + sim.ddz(temp)
	az = az + sim.ddx(temp)

	ay = sim.ddx(sim.makeSpect(by*dbx + dby * bx))*alpha
	ay = ay + sim.ddz(sim.makeSpect(by*dbz + dby * bz))*alpha
	
	lp = pressureResponse(sim,ax,az)

	retx = np.fft.irfftn(ax + lp[0])
	rety = np.fft.irfftn(ay)
	retz = np.fft.irfftn(az + lp[1])
	retp = np.fft.irfftn(lp[2])
	return [retx, rety, retz, retp]

def diffAccel(sim,u,v,w,pr):
	u = sim.makeSpect(u)
	v = sim.makeSpect(v)
	w = sim.makeSpect(w)
	ax = (sim.ddx(sim.ddx(u)) + sim.ddz(sim.ddz(u)))*pr
	ay = (sim.ddx(sim.ddx(v)) + sim.ddz(sim.ddz(v)))*pr
	az = (sim.ddx(sim.ddx(w)) + sim.ddz(sim.ddz(w)))*pr

	dp = pressureResponse(sim,ax,az)

	retx = np.fft.irfftn(ax + dp[0])
	rety = np.fft.irfftn(ay)
	retz = np.fft.irfftn(az + dp[1])
	retp = np.fft.irfftn(dp[2])
	return [retx, rety, retz, retp]

def induct(sim,u,v,w,bx,by,bz):
	ix,iy,iz = cross(u,v,w,bx,by,bz)
	return sim.curl(ix,iy,iz)

def dInduct(sim,u,v,w,du,dv,dw,bx,by,bz,dbx,dby,dbz):
	t1x,t1y,t1z = cross(du,dv,dw,bx,by,bz)
	t2x,t2y,t2z = cross(u,v,w,dbx,dby,dbz)

	return sim.curl(t1x+t2x,t1y+t2y,t1z+t2z)

def magDiff(sim,bx,by,bz,diff):
	bx = sim.makeSpect(bx)
	by = sim.makeSpect(by)
	bz = sim.makeSpect(bz)

	dbx = (sim.ddx(sim.ddx(bx)) + sim.ddz(sim.ddz(bx)))*diff
	dby = (sim.ddx(sim.ddx(by)) + sim.ddz(sim.ddz(by)))*diff
	dbz = (sim.ddx(sim.ddx(bz)) + sim.ddz(sim.ddz(bz)))*diff

	dbx = np.fft.irfftn(dbx)
	dby = np.fft.irfftn(dby)
	dbz = np.fft.irfftn(dbz)
	return [dbx,dby,dbz]

def inverseLaplace(sim, force):

	nz = sim.nz
	nx = sim.nx

	transform = False
	if (np.isrealobj(force)):
		sim.makeSpect(force)
		transform = True

	sim.assertDkx(force)
	sim.assertDkz(force)

	pres = force / (sim.dkz**2 + sim.dkx**2)
	pres[0,0,0] = 0

	if (transform == True):
		pres = np.fft.irrfn(pres)

	return pres

def pressureResponse(sim, forcex, forcez):

	nz = sim.nz
	nx = sim.nx

	transform = False
	if(np.isrealobj(forcex)):
		forcex = sim.makeSpect(forcex)
		forcez = sim.makeSpect(forcez)
		transform = True

	sim.assertDkx(forcex)
	sim.assertDkz(forcex)

	div = forcex * sim.dkx + forcez * sim.dkz

	pres = div / (sim.dkz**2 + sim.dkx**2)
	pres[0,0,0] = 0
 	solx = -pres * sim.dkx
	solz = -pres * sim.dkz

	if (transform == True):
		solx = np.fft.irfftn(solx)
		solz = np.fft.irfftn(solz)
		pres = np.fft.irfftn(pres)

	return [solx,solz,pres]

def calcPVel(sim,p,dp,ddp,ind):
	# _c refers to complex variables
	# _st refers to full derivatives (d/dt + v*del)

	#easier just to have everyting spectral until explicitly want nl term
	if(np.isrealobj(p)):
		p = sim.makeSpect(p)
		dp = sim.makeSpect(dp)
		ddp = sim.makeSpect(ddp)

	#grab terms for H(p), and nab gradient for free along the way
	dxP_c = sim.ddx(p)
	dxzP_c = sim.ddz(dxP_c)
	dxxP_c = sim.ddx(dxP_c)

	dzP_c = sim.ddz(p)
	dzzP_c = sim.ddz(dzP_c)

	#no more differentiating these.  Make them real, and ditch full array
	dxxP = np.fft.irfftn(dxxP_c)[ind]
	dxzP = np.fft.irfftn(dxzP_c)[ind]
	dzzP = np.fft.irfftn(dzzP_c)[ind]

	#grab gradient of dp
	dtxP_c = sim.ddx(dp)
	dtzP_c = sim.ddz(dp)

	#set up array of test thetas
	theta = np.array(range(1000),dtype=np.float64)*2*np.pi/1000
	x = np.cos(theta)
	z = np.sin(theta)

	#get spatial versions of grad dp
	dtxP = np.fft.irfftn(dtxP_c)[ind]
	dtzP = np.fft.irfftn(dtzP_c)[ind]

	#calculate V for the full circle of directions, and grab the best one
	T = x*dtxP  + z*dtzP
	B = dxxP*x*x + 2*dxzP*x*z + dzzP*z*z
	best = np.argmax(-T/B)
	
	#no need for the full array of directions any more
	x = x[best]
	z = z[best]
	T = T[best]
	B = B[best]
	speed = -T/B

	P_cst = dp + x*speed*dxP_c + z*speed*dzP_c
	dtP_cst = ddp + x*speed*dtxP_c + z*speed*dtzP_c

	#gradient of dp_Star 
	dtxP_cst = sim.ddx(dtP_cst)
	dtzP_cst = sim.ddz(dtP_cst)

	dtxP_st = np.fft.irfftn(dtxP_cst)[ind]
	dtzP_st = np.fft.irfftn(dtzP_cst)[ind]

	#terms for the hessian of P star
	dxxP_cst = sim.ddx(P_cst)
	dxzP_cst = sim.ddz(dxxP_cst)
	dxxP_cst = sim.ddx(dxxP_cst)
	dzzP_cst = sim.ddz(sim.ddz(P_cst))

	dxxP_st = np.fft.irfftn(dxxP_cst)[ind]
	dxzP_st = np.fft.irfftn(dxzP_cst)[ind]
	dzzP_st = np.fft.irfftn(dzzP_cst)[ind]

	#calculate all of the D terms:

	D1 = x*dtzP - z*dtxP
	D2 = x*dtxP_st + z*dtzP_st
	D3 = z*x*(dzzP-dxxP) + (x*x - z*z)*dxzP
	D4 = x*x*dxxP_st + 2*x*z*dxzP_st + z*z*dzzP_st
	D5 = z*x*(dzzP_st-dxxP_st) + (x*x - z*z)*dxzP_st
	D6 = x*dtzP_st - z*dtxP_st
	D7 = z*z*dxxP - 2*x*z*dxzP + x*x*dzzP
	
	d_th = (2*D5*T + 2*D3*D2 - B*D6 - D4*D1)/(B*T - 2*D7*T)

	d_mag_v = -((d_th * D1 + D2)*B - T*(2*d_th*D3 + D4))/(B**2)

	a_x = x * d_mag_v - z * d_th * speed
	a_z = z * d_mag_v + x * d_th * speed

	return [x*speed, z*speed, speed, a_x, a_z, np.sqrt(a_x**2 + a_z**2)] 
	
