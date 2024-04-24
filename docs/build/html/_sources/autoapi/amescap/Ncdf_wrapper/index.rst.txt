:orphan:

:py:mod:`amescap.Ncdf_wrapper`
==============================

.. py:module:: amescap.Ncdf_wrapper


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   amescap.Ncdf_wrapper.Ncdf
   amescap.Ncdf_wrapper.Fort




.. py:class:: Ncdf(filename=None, description_txt='', action='w', ncformat='NETCDF4_CLASSIC')


   Bases: :py:obj:`object`

   Alex K.
   NetCdf wrapper for quick archiving of data into netcdf format

   USAGE:

   from netcdf_wrapper import Ncdf

   Fgeo= 0.03 #W/m2, a constant
   TG=np.ones((24,8)) #ground temperature

   #---create file---
   filename="/lou/s2n/mkahre/MCMC/analysis/working/myfile.nc"
   description="results from new simulation, Alex 01-01-19"
   Log=Ncdf(filename,description)

   #---Save the constant to the file---
   Log.add_constant('Fgeo',Fgeo,"geothermal flux","W/m2")

   #---Save the TG array to the file---
   Log.add_dimension('Nx',8)
   Log.add_dimension('time',24)

   Log.log_variable('TG',TG,('time','Nx'),'soil temperature','K')

   Log.close()




.. py:class:: Fort(filename=None, description_txt='')


   Bases: :py:obj:`object`

   Alex K.
   A class that generate an object from fort.11 ... with similar attributes as a Netcdf file, e.g.:
   >>  f.variables.keys()
   >>  f.variables['var'].long_name
   >>  f.variables['var'].units
   >>  f.variables['var'].dimensions

   Create a Fort object using the following:
   f=Fort('/Users/akling/test/fort.11/fort.11_0684')

   PUBLIC METHODS:
   >> f.write_to_fixed(), f.write_to_average()  f.write_to_daily()  and f.write_to_diurn() can be used to generate FV3-like netcdf files

   .. py:class:: Fort_var(input_vals, name_txt, long_name_txt, units_txt, dimensions_tuple)


      Bases: :py:obj:`numpy.ndarray`

      Sub-class that emulate a netcdf-like variable by adding the name, long_name, units, dimensions attribute to a numpy array. [A. Kling]
      *** NOTE***
      A useful resource on subclassing in available at:
      https://numpy.org/devdocs/reference/arrays.classes.html

      Note that because we use an existing numpy.ndarray to define the object, we do not use a call to __array_finalize__(self, obj)


   .. py:method:: write_to_fixed()

      Create 'fixed' file, i.e.  all static variables


   .. py:method:: write_to_daily()

      Create daily file, e.g. contineuous time serie


   .. py:method:: write_to_average(day_average=5)

      Create average file, e.g. N day averages (typically 5)


   .. py:method:: write_to_diurn(day_average=5)

      Create diurn file, e.g.variable are organized by time of day. Additionally, the data is also binned  (typically 5)



