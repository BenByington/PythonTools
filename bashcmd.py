import os

def cd(loc,save=False):
	try:
		ret = os.getcwd()
		os.chdir(loc)
		if save == True:
			return ret
	except OSError:
		print "Directory '" + loc + "' does not exist."
		raise

def pwd():
	return os.getcwd()


def mkdir(path):
	return os.mkdir(path)

def ls(path='.',all=False):
	temp = os.listdir(path)
	ret = []
	if all == False:
		for f in temp:
			if f.startswith('.'):
				pass
			else:
				ret.append(f)
	else:
		ret = temp

	ret.sort()
	return ret
	

def get_data_files(nz=6,path='.'):
	list = ls(path=path)
	for i in range(len(list)):
		try:
			list[i] = int(list[i])
		except ValueError:
			pass
	list = filter(lambda x: type(x) is int, list)
	list.sort()
	for i in range(len(list)):
		list[i] = str(list[i]).zfill(nz)
	return list

