from setuptools import setup, find_packages

setup(name='amesgcm',
      version='0.2',
      description='Analysis pipeline for the NASA Ames MGCM',
      url='http://github.com/alex-kling/amesgcm',
      author='Mars Climate Modeling Center',
      author_email='alexandre.m.kling@nasa.gov',
      license='TBD',
      scripts=['bin/MarsPull.py','bin/MarsInterp.py','bin/MarsPlot.py','bin/MarsVars.py','bin/MarsFiles.py','bin/MarsViewer.py'],
      install_requires=['requests','netCDF4','numpy>=1.18.0','matplotlib','scipy'],
      packages=['amesgcm'],
      data_files = [('mars_data', ['mars_data/Legacy.fixed.nc']),('mars_templates', ['mars_templates/legacy.in','mars_templates/amesgcm_profile'])],
      include_package_data=True,
      zip_safe=False)
