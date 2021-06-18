import sys
import subprocess

py_exec = sys.executable
py_prefix = sys.exec_prefix
# ensure pip is installed & update
subprocess.call([str(py_exec), "-m", "ensurepip", "--user"])
subprocess.call([str(py_exec), "-m", "pip", "install", "--target={}".format(py_prefix), "--upgrade", "pip"])
# install dependencies using pip
# dependencies such as 'numpy' could be added to the end of this command's list
subprocess.call([str(py_exec),"-m", "pip", "install", "--target={}".format(py_prefix), "scipy"])