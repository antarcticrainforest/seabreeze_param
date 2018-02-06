from numpy.distutils.core import setup, Extension


meta=dict(\
        description= "A sub grid seabreeze parametrization for weather and climate models",
        url = 'https://github.com/antarcticrainforest/seabreezediag',
        author = 'Martin Bergemann',
        author_email = 'martin.bergemann@met.fu-berlin.de',
        license = 'GPL',
        version = '0.1')

wrapper = Extension('seabreeze', sources=['seabreezediag/seabreeze_diag_python.f90','seabreezediag/sobel.f90'],
        extra_f90_compile_args=["-g -O3 -std=f2003  -fsign-zero -fbounds-check -Wpedantic -shared -fPIC -fopenmp"],
        library_dirs=['/usr/lib'],
        include_dirs=['/usr/include'],
        libraries=['gomp'])

setup(name='seabreezediag',ext_modules=[wrapper] , packages=['seabreezediag'],**meta)



