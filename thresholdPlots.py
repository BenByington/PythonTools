import numpy as np
import matplotlib.pyplot as plt

def thresholdPlots(thresh, b2, by, t, printOnly = False):
	
	thresh.sort()
	print thresh
	thresh = [0] + thresh + [1]
	thresh = np.array(thresh)*np.max(b2)

	ht,e = np.histogram(b2,thresh,weights=t)
	hb1,e = np.histogram(b2,thresh,weights=b2)
	hb2,e = np.histogram(b2,thresh,weights=np.sqrt(b2))
	hb3,e = np.histogram(b2,thresh,weights=by)
	ha,e = np.histogram(b2,thresh)
	s = np.sum(t)
	print "Percent Mag Eng: ",hb1 / np.sum(b2)
	print "Percent Mag Amp: ", hb2 / np.sum(np.sqrt(b2))
	print "Percent By: ", hb3 / np.sum(by)
	print "Percent tracer: ", ht / s
	ret = ht/s

	if printOnly == True:
		return ret

	state = plt.isinteractive()
	if(state == False):
		plt.ion()

	for i in range(len(thresh)-1):
		plt.figure(i)
		plt.clf()
		plt.imshow(b2[:,0,:],origin=0,vmin=thresh[i],vmax=thresh[i+1])
		plt.colorbar()
		plt.axis('tight')
		plt.title("Range: " + str(thresh[i]) + " :: " + str(thresh[i+1]) + "\n " + str(100*ret[i]) + "%")

#	plt.figure(i+1)
#	plt.clf()
#	plt.imshow(np.log10(b2[:,0,:]),origin=0,vmin=np.log10(thresh[1]/100))
#	plt.axis('tight')
#	plt.colorbar()

#	plt.figure(i+2)
#	plt.clf()
#	plt.imshow(b2[:,0,:],origin=0)
#	plt.axis('tight')
#	plt.colorbar()

	if(state == False):
		plt.ioff()

	return ret

