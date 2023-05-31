from setuptools import setup, find_packages

setup(name='amescap',
      version='0.3',
      description='Analysis pipeline for the NASA Ames MGCM',
      url='https://github.com/NASA-Planetary-Science/AmesCAP',
      author='Mars Climate Modeling Center',
      author_email='alexandre.m.kling@nasa.gov',
      license='TBD',
      scripts=['bin/MarsPull.py','bin/MarsInterp.py','bin/MarsPlot.py','bin/MarsVars.py','bin/MarsFiles.py','bin/MarsViewer.py'],
      install_requires=['requests','netCDF4','numpy','matplotlib','scipy'],
      packages=['amescap'],
      data_files = [('mars_data', ['mars_data/Legacy.fixed.nc']),('mars_templates', ['mars_templates/legacy.in','mars_templates/amescap_profile'])],
      include_package_data=True,
      zip_safe=False)
