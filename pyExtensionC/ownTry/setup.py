from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(
    ext_modules = cythonize(Extension("cos_module", ["lib/dblcos.c", "lib/cos_module.pyx"]))
)

