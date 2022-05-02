#Fisher-py does not include a method for closing open file handles in the 'RawFile' class. These functions remedy that.

from fisher_py.raw_file_reader.raw_file_access import RawFileAccess
from fisher_py import RawFile

def close_raw_file(self):
	"""close .raw files"""
	self._raw_file_access.dispose()
	
def assign_close_method():
	"""assign method 'close_raw_file' to RawFile class"""
	RawFile.close_raw_file = close_raw_file


