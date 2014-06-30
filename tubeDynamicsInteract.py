import matplotlib.pyplot as plt

from vorts_template2 import *
from vorts_template_leaky2 import *
from vorts_template_detail2 import *

def loadVorts(pm,pr,q):
	full = vorts_template2(pm,pr,q,key="full").getData()
	adv = vorts_template2(pm,pr,q,key="advect").getData()
	diff = vorts_template2(pm,pr,q,key="diff").getData()
	buoy = vorts_template2(pm,pr,q,key="buoy").getData()
	lor = vorts_template2(pm,pr,q,key="lor").getData()

	return {'full':full,'adv':adv,'diff':diff,'buoy':buoy,'lor':lor}

def loadVortsLeaky(pm,pr,q,back):
	full = vorts_template_leaky2(pm,pr,q,back,key="full").getData()
	adv = vorts_template_leaky2(pm,pr,q,back,key="advect").getData()
	diff = vorts_template_leaky2(pm,pr,q,back,key="diff").getData()
	buoy = vorts_template_leaky2(pm,pr,q,back,key="buoy").getData()
	lor = vorts_template_leaky2(pm,pr,q,back,key="lor").getData()

	return {'full':full,'adv':adv,'diff':diff,'buoy':buoy,'lor':lor}

def loadVortsDetail(pm,pr,q,id):
	full = vorts_template_detail2(pm,pr,q,ident=id,key="full").getData()
	adv = vorts_template_detail2(pm,pr,q,ident=id,key="advect").getData()
	diff = vorts_template_detail2(pm,pr,q,ident=id,key="diff").getData()
	buoy = vorts_template_detail2(pm,pr,q,ident=id,key="buoy").getData()
	lor = vorts_template_detail2(pm,pr,q,ident=id,key="lor").getData()

	return {'full':full,'adv':adv,'diff':diff,'buoy':buoy,'lor':lor}

def assignVars(variables,vorts):
	variables.update(vorts)

def plotVars(v,i,title=None,f=None):
	full = v['full']
	adv = v['adv']
	buoy = v['buoy']
	lor = v['lor']
	diff = v['diff']

	titles = {0 : "Time", 1 : "x Pos", 2:"z Pos", 3:"x Vel", 4:"z Vel", 5:"Total Vel", 6:"x Accel", 7:"z Accel", 8:"Total Accel" }

	if(f != None):
		plt.clf()
	else:
		f = plt.figure()
	plt.plot(full[:,0],full[:,i],label="full")
	plt.plot(adv[:,0],adv[:,i],label="adv")
	plt.plot(buoy[:,0],buoy[:,i],label="buoy")
	plt.plot(lor[:,0],lor[:,i],label="lor")
	plt.plot(diff[:,0],diff[:,i],label="diff")
	plt.legend()
	if(title == None):
		plt.title(titles[i])
	else:
		plt.title(title)
	
	return f
	
def setXR(f,xr):
	try:
		for fig in f:
			plt.figure(fig.number)
			plt.xlim(xr)
	except TypeError:
		plt.figure(f.number)
		plt.xlim(xr)

def setYR(f,yr):
	try:
		for fig in f:
			plt.figure(fig.number)
			plt.ylim(yr)
	except TypeError:
		plt.figure(f.number)
		plt.ylim(yr)
