from distutils.core import setup
from distutils.extension import Extension
#from Cython.Distutils import build_ext
from Cython.Build import cythonize

#import os
#path_src = os.path.dirname(os.path.abspath(__file__)) + '/src'

#os.system('export CFLAGS="-I /opt/aci/sw/python/3.6.3_anaconda-5.0.1/lib/python3.6/site-packages/numpy/core/include $CFLAGS"')
extensions = [
    Extension("TGA",
               sources = ["TGA.pyx","Phase.cpp","Species.cpp", "Reactions.cpp"],
               language = "c++",
               include_dirs=[]),
]      

setup(
      name = "TGA_ext",
      ext_modules = cythonize(extensions,annotate=True),   
)

extensions = [
    Extension("CRT",
               sources = ["CRT.pyx","Phase.cpp","Species.cpp", "Reactions.cpp"],
               language = "c++",
               include_dirs=[]),
]      

setup(
      name = "CRT_ext",
      ext_modules = cythonize(extensions,annotate=True),   
)

