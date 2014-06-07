# Shared functions
import gzip

################################################################################

def open2(file, mode="r"):
	"""
	Open a file, or a gzipped file if it ends in .gz
	"""
	if file[-3:]==".gz":
		return gzip.open(file, mode)
	else:
		return open(file, mode)

################################################################################

