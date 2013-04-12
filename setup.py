from distutils.core import setup, Extension

nativemodule = Extension('diffsim',
                         sources=['bindings.c', 'calcsim.c'])
setup(name='diffsim',
      version='1.0',
      description='Diffusion simulation',
      ext_modules=[nativemodule])