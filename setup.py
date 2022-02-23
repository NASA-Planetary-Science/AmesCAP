from setuptools import setup, find_packages

setup(name='amesgcmcb',
      version='0.2',
      description='Analysis pipeline for the NASA Ames MGCM - working repository',
      url='http://github.com/falconstryker/amesgcm',
      author='Mars Climate Modeling Center',
      author_email='cmlbatterson@gmail.com',
      license='TBD',
      scripts=['bin/MarsPull.py','bin/MarsInterp.py','bin/MarsPlot.py','bin/MarsVars.py','bin/MarsFiles.py','bin/MarsViewer.py'],
      install_requires=['requests','netCDF4','numpy','matplotlib','scipy'],
      packages=['amesgcm'],
      data_files = [('mars_data', ['mars_data/Legacy.fixed.nc']),('mars_templates', ['mars_templates/legacy.in','mars_templates/amesgcm_profile'])],
      include_package_data=True,
      zip_safe=False)
