from setuptools import setup

setup(name='amesgcm',
      version='0.1',
      description='Analysis pipeline for the NASA Ames MGCM',
      url='http://github.com/alex-kling/amesgcm',
      author='Mars Climate Modeling Center',
      author_email='alexandre.m.kling@nasa.gov',
      license='TBD',
      scripts=['bin/MarsPull.py'],
      install_requires=['requests',],
      packages=['amesgcm'],
      zip_safe=False)
