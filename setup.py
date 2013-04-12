from distutils.core import setup, Extension
import numpy

nativemodule = Extension('diffsim',
                         sources=['bindings.c', 'calcsim.c'],
                         extra_compile_args=['-std=c99'])
setup(name='diffsim',
      version='1.0',
      description='Diffusion simulation',
      include_dirs=[numpy.get_include()],
      ext_modules=[nativemodule])