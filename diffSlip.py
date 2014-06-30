from tubeDynamics2 import diffAccel
import numpy as np

def calcSlip(s,bx,bz,zeta):
	dbx = (s.ddx(s.ddx(bx)) + s.ddz(s.ddz(bx)))*zeta
	dbz = (s.ddx(s.ddx(bz)) + s.ddz(s.ddz(bz)))*zeta
	

	bnorm = np.sqrt(bx*bx + bz*bz)
	bnx = np.zeros_like(bx)
	bnz = np.zeros_like(bz)
	valid = bnorm > 0
	bnx[valid] = bx[valid] / bnorm[valid]
	bnz[valid] = bz[valid] / bnorm[valid]

	temp = (dbx*bnx + dbz*bnz)/bnorm
	dbnx = dbx/bnorm - temp*bnx
	dbnz = dbz/bnorm - temp*bnz

	ax = s.ddx(bnx*bnx) + s.ddz(bnx*bnz)
	az = s.ddz(bnz*bnz) + s.ddx(bnx*bnz)

	dax = s.ddx(2*bnx*dbnx) + s.ddz(bnx*dbnz + bnz*dbnx)
	daz = s.ddz(2*bnz*dbnz) + s.ddx(bnx*dbnz + bnz*dbnx)

	anorm = np.sqrt(ax*ax + az*az)
	rad = 1/anorm

	drad = (ax*dax + ax*daz) / (anorm ** 3)

	gradx = s.ddx(drad)
	gradz = s.ddz(drad)

	dvdensity = gradx*bnx + gradz*bnz

	return {"curv" : [ax, az, anorm], "rad":rad, "drad":drad, "dvdensity":dvdensity}
