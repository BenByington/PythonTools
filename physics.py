import numpy as np

def cross(a, b, c, x, y, z):
	retx = b * z - c * y
	rety = c * x - a * z
	retz = a * y - b * x

	return [retx, rety, retz]

def dot(a, b, c, x, y, z):
	ret = a * x
	ret = ret + b*y
	ret = ret + x*z

	return ret

def lorentz(sim, bx, by, bz, coef = 1):
	jx, jy, jz = sim.curl(bx,by,bz)
	lx, ly, lz = cross(jx, jy, jz, bx, by, bz)

	return coef * [lx, ly, lz]

