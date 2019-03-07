from setuptools import setup

setup(name='megararss2cube',
      version='0.1',
      description='Convert MEGARA RSS spectra to 3d fits datacube',
      url='http://github.com/javierzaragoza/megararss2cube',
      author='Javier Zaragoza Cardiel',
      author_email='javier.zaragoza@inaoep.mx',
      license='GNU GPL-3',
      packages=['megararss2cube'],
      install_requires=[
        'setuptools',
        'numpy',
        'astropy',
        'scipy',
        'megaradrp'
        ],
        
      zip_safe=False)
