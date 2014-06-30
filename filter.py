import numpy as np

def filter_conv(fieldx, fieldy=None, fieldz=None, rho=None):

	freq = fieldx.shape[0]/2

	ret = []

	if fieldy == None or fieldz == None:
		spec = np.fft.rfftn(fieldx)
		scale = np.max(np.abs(fieldx))
		for j in range(freq):
			temp = filter(spec,j)
			error = np.fft.irfftn(temp[1])
			ret.append(np.max(np.abs(error))/scale)
	else:
		if rho != None:
			temp = np.sqrt(rho)
			specx = np.fft.rfftn(fieldx*temp)
			specy = np.fft.rfftn(fieldy*temp)
			specz = np.fft.rfftn(fieldz*temp)
			engScale = np.max(rho*(fieldx**2 + fieldy**2 + fieldz**2))
		else:
			specx = np.fft.rfftn(fieldx)
			specy = np.fft.rfftn(fieldy)
			specz = np.fft.rfftn(fieldz)
			engScale = np.max(fieldx**2 + fieldy**2 + fieldz**2)
			
			
		for j in range(freq):
			a = filter(specx,j)
			b = filter(specy,j)
			c = filter(specz,j)
			error = np.fft.irfftn(a[1])**2
			error += np.fft.irfftn(b[1])**2
			error += np.fft.irfftn(c[1])**2
			ret.append(np.max(error)/engScale)

	return ret

def filter_single(wn, fieldx, fieldy=None, fieldz=None, rho=None, calcError = 0):

	if fieldy == None or fieldz == None:
		spec = np.fft.rfftn(fieldx)
		temp = filter(spec,wn)
		retl = np.fft.irfftn(temp[0])
		rets = np.fft.irfftn(temp[1])
		if(calcError != 0):
			error = np.max(np.abs(fieldx - retl - rets)) / np.average(np.abs(fieldx))
			error2 = (np.mean(np.abs(fieldx)) - np.mean(np.abs(retl+rets)))/np.mean(np.abs(fieldx))
			print error
			print error2
			return [retl,rets,error]
		else:
			return [retl,rets]
	else:
		if rho != None:
			tmp = np.sqrt(rho)
			tempx = fieldx*tmp
			tempy = fieldy*tmp
			tempz = fieldz*tmp
		else:
			tempx = fieldx
			tempy = fieldy
			tempz = fieldz
			
		specx = np.fft.rfftn(tempx)
		specy = np.fft.rfftn(tempy)
		specz = np.fft.rfftn(tempz)

		a = filter(specx,wn)
		b = filter(specy,wn)
		c = filter(specz,wn)
		a=[np.fft.irfftn(a[0]),np.fft.irfftn(a[1])]
		b=[np.fft.irfftn(b[0]),np.fft.irfftn(b[1])]
		c=[np.fft.irfftn(c[0]),np.fft.irfftn(c[1])]
		return [a[0]**2 + b[0]**2  + c[0]**2, a[1]**2 + b[1]**2 + c[1]**2]

def filter(field, wn):
	if wn <= 0:
		wn = 1

	lfield = np.copy(field)
	sfield = np.zeros_like(field)

	size = field.shape
	if len(size) == 2:
		try:
			w = [wn[0],wn[1]]
		except TypeError:
			w = [wn,wn]
		sfield[w[0]:size[0]-w[0]+1,:] = lfield[w[0]:size[0]-w[0]+1,:]
		lfield[w[0]:size[0]-w[0]+1,:] = 0
		sfield[:,w[1]:] = lfield[:,w[1]:]
		lfield[:,w[1]:] = 0
	elif len(size) == 3:
		try:
			w = [wn[0],wn[1],wn[2]]
		except TypeError:
			w = [wn,wn,wn]
		sfield[w[0]:size[0]-w[0]+1,:,:] = lfield[w[0]:size[0]-w[0]+1,:,:]
		lfield[w[0]:size[0]-w[0]+1,:,:] = 0
		sfield[:,w[1]:size[1]-w[1]+1,:] = lfield[:,w[1]:size[1]-w[1]+1,:]
		lfield[:,w[1]:size[1]-w[1]+1,:] = 0
		sfield[:,:,w[2]:] = lfield[:,:,w[2]:]
		lfield[:,:,w[2]:] = 0
	else:
		print "unsupported data dimension!"

	return [lfield, sfield]
