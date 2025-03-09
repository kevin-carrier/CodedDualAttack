import ctypes
import numpy as np

class stat_t(ctypes.Structure):
    _fields_ = [('mean', ctypes.c_double),    # use correct types
                ('std_dev',ctypes.c_double),
                ('confidence', ctypes.c_double)]

lib_polar = ctypes.CDLL('./polar.o', mode=ctypes.RTLD_GLOBAL)

#PolarHandle = ctypes.POINTER(ctypes.c_char)
PolarHandle = ctypes.c_void_p

lib_polar._random_polar_code.argtypes = [ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_ulong,ctypes.c_ulong,ctypes.c_int,ctypes.c_int]
lib_polar._random_polar_code.restype = PolarHandle


lib_polar._copy_polar.argtypes = [PolarHandle]
lib_polar._copy_polar.restype = PolarHandle



lib_polar._init_random_seed.argtypes = [ctypes.c_uint]
lib_polar._init_random_seed.restype = None

lib_polar._init_random.argtypes = []
lib_polar._init_random.restype = None


lib_polar._free_polar_code.argtypes = [PolarHandle]
lib_polar._free_polar_code.restype = None

lib_polar._mean_error.argtypes = [PolarHandle,ctypes.c_ulong,ctypes.c_ulong]
#lib_polar._mean_error.restype = ctypes.c_double
lib_polar._mean_error.restype = stat_t

lib_polar._print_polar.argtypes = [PolarHandle]
lib_polar._print_polar.restype = None

lib_polar._gen_random_codeword.argtypes = [PolarHandle, ctypes.POINTER(ctypes.c_long)]
lib_polar._gen_random_codeword.restype = None

lib_polar._decode.argtypes = [PolarHandle, ctypes.POINTER(ctypes.c_long),ctypes.POINTER(ctypes.c_long),ctypes.c_ulong]
lib_polar._decode.restype = None

lib_polar._generator_matrix.argtypes = [PolarHandle, ctypes.POINTER(ctypes.c_long)]
lib_polar._generator_matrix.restype = None


lib_polar._init_random()

def polar_copy(C_polar):
	return lib_polar._copy_polar(C_polar)
def polar_print(C_polar):
	lib_polar._print_polar(C_polar)
def polar_random(q,n,k, nb_test_code = 5, nb_samples_polarization = 100, list_size = 2, nb_test_mean = 100):
	return lib_polar._random_polar_code(q,n,k,nb_test_code,nb_samples_polarization,list_size,nb_test_mean)
def polar_free(C_polar):
	lib_polar._free_polar_code(C_polar)
def polar_mean_error(C_polar, list_size = 2, nb_test = 1):
	res = lib_polar._mean_error(C_polar,list_size, nb_test)
	return res.mean
def polar_stats(C_polar, list_size = 2, nb_test = 1):
	res = lib_polar._mean_error(C_polar,list_size, nb_test)
	return [res.mean,res.std_dev,res.confidence]
def polar_decode(C_polar, noisy,list_size = 2):
	n = len(noisy)
	decoded = [0]*n
	C_noisy = (ctypes.c_long*n)(*list(noisy))
	C_decoded = (ctypes.c_long*n)(*decoded)
	lib_polar._decode(C_polar, C_noisy,C_decoded,list_size)
	return list(C_decoded)
def polar_random_codeword(C_polar,n):
	codeword = [0]*n
	C_codeword = (ctypes.c_long*n)(*codeword)
	lib_polar._gen_random_codeword(C_polar, C_codeword)
	return list(C_codeword)



def polar_str_repr(C_polar):
	from contextlib import contextmanager
	from contextlib import redirect_stdout
	import io
	import os, sys
	import tempfile
	libc = ctypes.CDLL(None)
	c_stdout = ctypes.c_void_p.in_dll(libc, 'stdout')
	@contextmanager
	def stdout_redirector(stream):
		# The original fd stdout points to. Usually 1 on POSIX systems.
		original_stdout_fd = sys.stdout.fileno()

		def _redirect_stdout(to_fd):
			"""Redirect stdout to the given file descriptor."""
			# Flush the C-level buffer stdout
			libc.fflush(c_stdout)
			# Flush and close sys.stdout - also closes the file descriptor (fd)
			sys.stdout.close()
			# Make original_stdout_fd point to the same file as to_fd
			os.dup2(to_fd, original_stdout_fd)
			# Create a new sys.stdout that points to the redirected fd
			sys.stdout = io.TextIOWrapper(os.fdopen(original_stdout_fd, 'wb'))

		# Save a copy of the original stdout fd in saved_stdout_fd
		saved_stdout_fd = os.dup(original_stdout_fd)
		try:
			# Create a temporary file and redirect stdout to it
			tfile = tempfile.TemporaryFile(mode='w+b')
			_redirect_stdout(tfile.fileno())
			# Yield to caller, then redirect stdout back to the saved fd
			yield
			_redirect_stdout(saved_stdout_fd)
			# Copy contents of temporary file to the given stream
			tfile.flush()
			tfile.seek(0, io.SEEK_SET)
			stream.write(tfile.read())
		finally:
			tfile.close()
			os.close(saved_stdout_fd)
	f = io.BytesIO()
	with stdout_redirector(f):
		polar_print(C_polar)
		#os.system('echo and this is from echo')
	st = "{0}".format(f.getvalue().decode('utf-8'))
	return st