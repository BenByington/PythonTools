import numpy as np

from load_data import save_rec

def generate(nx=128,ny=8,nz=1024,q=1.0,rad=0.25/2,cx=0.5,cz=0.05):

	bx = np.zeros((nz,ny,nx),np.float64)
	by = np.zeros((nz,ny,nx),np.float64)
	bz = np.zeros((nz,ny,nx),np.float64)
	t = np.zeros((nz,ny,nx),np.float64)
	u = np.zeros((nz,ny,nx),np.float64)

	save_rec("u", u)
	save_rec("v", u)
	save_rec("w", u)

	midx = nx*cx
	midz = nz*cz

	for i in range(nz):
		for j in range(nx):
			r = np.sqrt((j - midx)**2 + (i-midz)**2)
			if r < nx*rad:
				by[i,:,j] = 1 - (r/nx/rad)**2
				t[i,:,j] = 1.0

	save_rec("By",by)
	save_rec("T",t)

	for i in range(nz):
		for j in range(nx):
			bx[i,:,j] = -2*q*(i-midz)*by[i,:,j] / rad / nx
			bz[i,:,j] = 2*q*(j-midx)*by[i,:,j] / rad / nx

	save_rec("Bx", bx)
	save_rec("Bz", bz)
