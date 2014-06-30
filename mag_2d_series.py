import numpy as np
import matplotlib.pyplot as plt
import os

def x_avg_plots(data):
	
	v = data[2]
	w = data[3]
	bx = data[6]
	by = data[7]
	bz = data[8]
	b2 = data[9]
	q = 4

	plt.cool()
	if hasattr(x_avg_plots,'x') == False:
		ny = v.shape[1]
		nz = v.shape[0]
		x_avg_plots.x,x_avg_plots.y = np.meshgrid(range(ny),range(nz))
	x,y = x_avg_plots.x,x_avg_plots.y
	
	avgv = np.average(v,2)
	avgw = -np.average(w,2)
	avgbx = np.average(bx,2)
	avgby = np.average(by,2)
	avgbz = -np.average(bz,2)
	avgb2 = np.average(b2,2)

	plt.subplot(221)
	plt.imshow(avgb2)
	plt.quiver(x[::q,::q],y[::q,::q],avgv[::q,::q],avgw[::q,::q])
	plt.title('B2')
	plt.axis("tight")

	plt.subplot(222)
	plt.imshow(avgbx)
	plt.quiver(x[::q,::q],y[::q,::q],avgv[::q,::q],avgw[::q,::q])
	plt.title('Bx')
	plt.axis("tight")

	plt.subplot(223)
	plt.imshow(avgby)
	plt.quiver(x[::q,::q],y[::q,::q],avgv[::q,::q],avgw[::q,::q])
	plt.title('By')
	plt.axis("tight")

	plt.subplot(224)
	plt.imshow(avgbz)
	plt.quiver(x[::q,::q],y[::q,::q],avgv[::q,::q],avgw[::q,::q])
	plt.title('Bz')
	plt.axis("tight")

