import numpy as np

def _is_blank(name, filename_categories_vocab):

    blank_search = filename_categories_vocab[2]
    if blank_search.upper() in name.upper():
        is_blank = True
    else:
        is_blank = False

    return is_blank

def print_progress_bar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):

    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def ppm_diff(observed, theoretical):
	"""Take two numbers as arguments, return PPM difference."""
	return (theoretical - observed) / theoretical * 1000000
