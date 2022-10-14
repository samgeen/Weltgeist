from skbuild import setup 

setup(
    name='weltgeist',
    version='0.911',
    description='A 1D hydrodynamic code for spherical astrophysics',
    license="MIT",
    packages=['weltgeist'],
    package_data={'weltgeist': ['Tgas.csv']},
    cmake_args=['-DSKBUILD=ON -fPIC']
)
