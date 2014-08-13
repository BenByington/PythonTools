import os

def cd(loc,save=False):
	"""Change directories the python shell works in

       Input :
         loc  -- Directory (relative or absolute) to switch to

       Output:
         save -- Optional parameter used to save the current directory
                 to facilite programatically switching back
       
       Side Effects: 
	The working directory is changed
	"""
	try:
		ret = os.getcwd()
		os.chdir(loc)
		if save == True:
			return ret
	except OSError:
		print "Directory '" + loc + "' does not exist."
		raise

def pwd():
	"""Prints the current working directory for the python shell

       Input : NA

       Output: NA
       
       Side Effects: 
         Prints to the console the current directory
	"""
	return os.getcwd()


def mkdir(path):
	"""Create a new directory

       Input :
         path  -- New directory to create

       Output: NA

       Side Effects:
         A new directory will be created in the file system (permissions 
         allowing)
	"""
	return os.mkdir(path)

def ls(path='.',all=False):
	"""Lists files and folders in the current directory

       Input  :
         path : The path of the directory you wish to list the contents of.  
                Defaults to current directory
         all  : a boolean value that determines if hidden files (beggining with
                a '.') are listed.  Defaults to false

       Output: Returns an alphabetic list of strings containing all the 
               files/folders found.
       
       Side Effects: NA
         
	"""
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
	"""Lists all the data files in a given folder created by the IMHD or HPS 
       codes.  Useful for getting a list of files when some files are missing 
       or for some reason there are different numbers of iterations between 
       file outputs.  These files are expected to have names containing only 
       numeric characters, and with a minimum width of nz.

       Input  : 
         nz   -- The zero-fill for the numbers used to label files.  i.e. if 
                 the first file is called 01 then nz should be 2.  Default 
                 value is 6 which is for the HPS code.  IMHD is currently 
                 configured to use 8.
         path -- The path to the folder we wish to examine.  Defaults to 
                 current working directory. 

       Output: Returns a sorted list of strings containing all the data files 
               found in the specified folder
       
       Side Effects: NA
         
	"""
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

