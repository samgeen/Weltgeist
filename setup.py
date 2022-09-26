import setuptools
from skbuild import setup 
#import setuptools 

#from numpy.distutils.core import setup
#from numpy.distutils.extension import Extension

#rayext = Extension('raytracing', ['Radiation/raytracing.f90'], extra_compile_args=['-O2'])
#stellar = Extension('singlestar',sources=['StellarSources/Fortran/singlestar_f2py.f90'],libraries=['StellarSources/Fortran/table_1d_module.f90',
#                                                                                                'StellarSources/Fortran/singlestar_module.f90',],
#                                                                                        language="f90",
#                                                                                        f2py_options=["-m"])

#def configuration(parent_package='weltgeist',top_path=None):
#    from numpy.distutils.misc_util import Configuration
#    config = Configuration('mypackage',parent_package,top_path)
#    return config

setup(
    name='weltgeist',
    version='0.9',
    description='1D hydrodynamic code for spherical astrophysics',
    license=["MIT"],
    packages=['weltgeist'],
    ext_modules = [rayext],
)