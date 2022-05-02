#Raw Class has been adapted from the RawFile class in fisher_py. Primarily so that .raw files can be closed after opening them

from fisher_py.raw_file_reader.raw_file_access import RawFileAccess
from fisher_py import RawFile


def close_raw_file(self):
	self._raw_file_access.dispose()
	
def assign_close_method():
	RawFile.close_raw_file = close_raw_file


