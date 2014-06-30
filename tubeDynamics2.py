import numpy as np
try:
	import matplotlib.pyplot as plt
except:
	"matplotlib is not on this machine.  tDynamics will fali if drawDp is used!"
from physics import cross

class tDynamics2:

	def __init__(self,sim, alpha, pr, mpr):
		self.sim = sim
		self.alpha = alpha
		self.pr = pr
		self.mpr = mpr

		self.saved = {}

	def calcTrajectory(self, iteration, xr=None,zr=None):

		try:
			ret = self.saved[iteration]
		except KeyError:

			p,dp,ddp = self.calcPress(iteration)

			if (xr == None):
				xr = [0,dp[0].shape[2]]
			if (zr == None):
				zr = [0,dp[0].shape[0]]

			temp = p[0][zr[0]:zr[1],:,xr[0]:xr[1]]
			ind = np.unravel_index(np.argmin(temp),temp.shape)
			ind = (ind[0] + zr[0], ind[1], ind[2] + xr[0])

			xp = ind[2]*self.sim.xmx/self.sim.nx
			zp = ind[0]*self.sim.zmx/self.sim.nz
			time = self.sim.get_time(iteration)

			temp = calcPVel(self.sim,p[0],dp[0],dp[0],ddp[0],ind)
			vel = [temp[3],temp[4]]

			ret = { 'full' : [time, xp , zp]  + temp }
			ret['advect'] = [time, xp, zp] + calcPVel(self.sim,p[0],dp[0],dp[1],ddp[1],ind)
			ret['diff'] = [time, xp, zp] + calcPVel(self.sim,p[0],dp[0],dp[2],ddp[2],ind)
			ret['buoy'] = [time, xp, zp] + calcPVel(self.sim,p[0],dp[0],dp[3],ddp[3],ind)
			ret['lor'] = [time, xp, zp] + calcPVel(self.sim,p[0],dp[0],dp[4],ddp[4],ind)
			self.saved[iteration] = ret

		return ret

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

		#Calculate the different components of the acceleration.
		adv_ax,adv_ay,adv_az,adv_p = advectionAccel(sim,u,v,w)
		diff_ax,diff_ay,diff_az = diffAccel(sim,u,v,w,pr)
		buoy_ax,buoy_ay,buoy_az,buoy_p = mBuoyAccel(sim,b2,alpha)
		lor_ax,lor_ay,lor_az,lor_p = lorentzAccel(sim,bx,by,bz,alpha)

		ax = adv_ax + diff_ax + buoy_ax + lor_ax
		ay = adv_ay + diff_ay + buoy_ay + lor_ay
		az = adv_az + diff_az + buoy_az + lor_az
		p = adv_p + buoy_p + lor_p
		if(pOnly == True):
			return p
		#grab the derivative of the B field

		dbx,dby,dbz = induct(sim,u,v,w,bx,by,bz)

		temp = magDiff(sim,bx,by,bz,pr/mpr)
		dbx += temp[0]
		dby += temp[1]
		dbz += temp[2]

		#start calculating dp.  Each version of DP has some portion of the derivative of the advenction term, so each of them contributes to adv_dau.  Each term contributes to a different pressure, so be careful.
		temp  = dAdvectionAccel(sim,u,v,w,adv_ax, adv_ay, adv_az)
		adv_dax = temp[0]
		adv_day = temp[1]
		adv_daz = temp[2]
		adv_dp = temp[3]

		temp = dAdvectionAccel(sim,u,v,w,diff_ax,diff_ay,diff_az)
		adv_dax += temp[0]
		adv_day += temp[1]
		adv_daz += temp[2]
		diff_dp = temp[3]

		temp = dAdvectionAccel(sim,u,v,w,buoy_ax,buoy_ay,buoy_az)
		adv_dax += temp[0]
		adv_day += temp[1]
		adv_daz += temp[2]
		buoy_dp = temp[3]

		temp  = dAdvectionAccel(sim,u,v,w,lor_ax, lor_ay, lor_az)
		adv_dax += temp[0]
		adv_day += temp[1]
		adv_daz += temp[2]
		lor_dp = temp[3]

		#grab the remaining bits for lorents and buoyancy
		temp = dMBuoyAccel(sim,bx,by,bz,dbx,dby,dbz,alpha)
		buoy_dax = temp[0]
		buoy_day = temp[1]
		buoy_daz = temp[2]
		buoy_dp += temp[3]

		temp = dLorentzAccel(sim,bx,by,bz,dbx,dby,dbz,alpha)
		lor_dax = temp[0]
		lor_day = temp[1]
		lor_daz = temp[2]
		lor_dp += temp[3]

		diff_dax,diff_day,diff_daz = diffAccel(sim,ax,ay,az,pr)

		dp = adv_dp + diff_dp + buoy_dp + lor_dp

		#We need the second derivative of the magnetic field to finish things off 
		ddbx,ddby,ddbz = dInduct(sim,u,v,w,ax,ay,az,bx,by,bz,dbx,dby,dbz)
		
		temp = magDiff(sim,dbx,dby,dbz,pr/mpr)
		ddbx += temp[0]
		ddby += temp[1]
		ddbz += temp[2]

		#every portion of ddp has a contribution from the advection term
		adv_ddp = ddpAdvect(sim,u,w,ax,az,adv_ax,adv_az,adv_dax,adv_daz)		
		diff_ddp = ddpAdvect(sim,u,w,ax,az,diff_ax,diff_az,diff_dax,diff_daz)		
		buoy_ddp = ddpAdvect(sim,u,w,ax,az,buoy_ax,buoy_az,buoy_dax,buoy_daz)		
		lor_ddp = ddpAdvect(sim,u,w,ax,az,lor_ax,lor_az,lor_dax,lor_daz)

		#add in the remaining contributions from lorentz and buoyancy
		buoy_ddp += ddpBuoy(sim,bx,by,bz,dbx,dby,dbz,ddbx,ddby,ddbz,alpha)
		lor_ddp += ddpLorentz(sim,bx,bz,dbx,dbz,ddbx,ddbz,alpha)		

		ddp = adv_ddp + diff_ddp + buoy_ddp + lor_ddp	

		#prepare return values
		retp = [p, adv_p, buoy_p, lor_p]
		retdp = [dp, adv_dp, diff_dp, buoy_dp, lor_dp]
		retddp = [ddp, adv_ddp, diff_ddp, buoy_ddp, lor_ddp]

		return [retp, retdp, retddp]

def ddpBuoy(sim, bx, by, bz, dbx, dby, dbz, ddbx, ddby, ddbz,alpha):
	nl1 = sim.makeSpect(bx*ddbx + by*ddby + bz*ddbz + dbx**2 + dby**2 + dbz**2)
	div = (2*alpha)*sim.ddz(nl1)
	return np.fft.irfftn(inverseLaplace(sim,div))

def ddpAdvect(sim, u, w, fax, faz, pax, paz, dpax, dpaz):
	nl1 = sim.makeSpect(fax*pax + u*dpax)
	nl2 = sim.makeSpect(u*dpaz + fax*paz + faz*pax + w*dpax)
	nl3 = sim.makeSpect(w*dpaz + faz*paz)
	div = sim.ddx(sim.ddx(nl1)) + sim.ddx(sim.ddz(nl2)) + sim.ddz(sim.ddz(nl3))
	return np.fft.irfftn(-2 * inverseLaplace(sim,div))

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

	retx = np.fft.irfftn(ax)
	rety = np.fft.irfftn(ay)
	retz = np.fft.irfftn(az)
	return [retx, rety, retz]

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

def calcPVel(sim,p,fdp,pdp, pddp,ind,vel=None):
	# _c refers to complex variables
	# _st refers to full derivatives (d/dt + v*del)

	#dfp is the full time derivative of p
	#pddp is the full time derivative of pdp
	#pdp is NOT the full time derivative of p.  It is the part of dp 
	#resulting from some process of interest (i.e. mag buoy)

	#easier just to have everyting spectral until explicitly want nl term
	if(np.isrealobj(p)):
		p = sim.makeSpect(p)
		fdp = sim.makeSpect(fdp)
		pdp = sim.makeSpect(pdp)
		pddp = sim.makeSpect(pddp)

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

	#grab gradient of pdp
	pdtxP_c = sim.ddx(pdp)
	pdtzP_c = sim.ddz(pdp)

	#set up array of test thetas
	theta = np.array(range(1000),dtype=np.float64)*2*np.pi/1000
	x = np.cos(theta)
	z = np.sin(theta)

	#get spatial versions of grad dp
	pdtxP = np.fft.irfftn(pdtxP_c)[ind]
	pdtzP = np.fft.irfftn(pdtzP_c)[ind]

	#calculate V for the full circle of directions, and grab the best one
	T = x*pdtxP  + z*pdtzP
	B = dxxP*x*x + 2*dxzP*x*z + dzzP*z*z
	best = np.argmax(-T/B)
	
	#no need for the full array of directions any more
	x = x[best]
	z = z[best]
	T = T[best]
	B = B[best]
	speed = -T/B

	if(vel == None):
		vel = [x*speed,z*speed]

	P_cst = fdp + vel[0]*dxP_c + vel[1]*dzP_c
	dtP_cst = pddp + vel[0]*pdtxP_c + vel[1]*pdtzP_c

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

	D1 = x*pdtzP - z*pdtxP
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
	
