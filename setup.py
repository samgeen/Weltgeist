from skbuild import setup 

setup(
    name='weltgeist',
    version='0.9',
    description='A 1D hydrodynamic code for spherical astrophysics',
    license="MIT",
    packages=['weltgeist'],
    cmake_args=['-DSKBUILD=ON']
)
