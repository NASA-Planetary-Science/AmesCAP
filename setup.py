from setuptools import setup, find_packages

setup(name='amesgcm',
      version='0.1',
      description='Analysis pipeline for the NASA Ames MGCM',
      url='http://github.com/alex-kling/amesgcm',
      author='Mars Climate Modeling Center',
      author_email='alexandre.m.kling@nasa.gov',
      license='TBD',
      scripts=['bin/MarsPull.py','bin/MarsDocumentation.sh','bin/MarsPlot.py','bin/MarsVars.py','bin/MarsFiles.py'],
      install_requires=['requests','netCDF4','numpy','matplotlib'],
      packages=['amesgcm'],
      data_files = [('mars_data', ['mars_data/Legacy.fixed.nc'])],
      include_package_data=True,
      zip_safe=False)
