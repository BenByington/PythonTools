from vorts_template import *
import matplotlib.pyplot as plt


def plots(pm, pr, q):
	
	plt.ion()
	vorts = vorts_template(pm,pr,q).getData()
	vortsH = vorts_template(pm,pr,q,hydro=True).getData()

	plt.figure()
	plt.plot((vorts[:,6]*vortsH[:,6] + vorts[:,7]*vortsH[:,7])/(vorts[:,8]*vortsH[:,8]))
	plt.title("Cos_th")

	plt.figure()
	plt.plot(vorts[:,8] / vortsH[:,8])
	plt.title("Accel Ratio")

	plt.figure()
	plt.plot(vorts[:,5] / vortsH[:,5])
	plt.title("Vel Ratio")
