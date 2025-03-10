:py:mod:`amescap.Ncdf_wrapper`
==============================

.. py:module:: amescap.Ncdf_wrapper

.. autoapi-nested-parse::

   Ncdf_wrapper archives data into netCDF format. It serves as a wrapper
   for creating netCDF files.

   Third-party Requirements:
       * ``numpy``
       * ``amescap.FV3_utils``
       * ``scipy.io``
       * ``netCDF4``
       * ``os``
       * ``datetime



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   amescap.Ncdf_wrapper.Fort
   amescap.Ncdf_wrapper.Ncdf




.. py:class:: Fort(filename=None, description_txt='')


   Bases: :py:obj:`object`

   A class that generates an object from a fort.11 file. The new file
   will have netCDF file attributes. Alex Kling.

   EX::

       f.variables.keys()
       f.variables['var'].long_name
       f.variables['var'].units
       f.variables['var'].dimensions

   Create a Fort object using the following::

       f=Fort('/Users/akling/test/fort.11/fort.11_0684')

   Public methods can be used to generate FV3-like netCDF files::

       f.write_to_fixed()
       f.write_to_average()
       f.write_to_daily()
       f.write_to_diurn()

   :param object: _description_
   :type object: _type_
   :return: _description_
   :rtype: _type_

   .. py:class:: Fort_var(input_vals, name_txt, long_name_txt, units_txt, dimensions_tuple)


      Bases: :py:obj:`numpy.ndarray`

      Sub-class that emulates a netCDF-like variable by adding the
      ``name``, ``long_name``, ``units``, and ``dimensions``
      attributes to a numpy array. Inner class for
      ``fortran_variables`` (Fort_var) that comprise the Fort file.
      Alex Kling

      A useful resource on subclassing is available at:
      https://numpy.org/devdocs/reference/arrays.classes.html

      .. NOTE:: Because we use an existing ``numpy.ndarray`` to define
          the object, we do not call ``__array_finalize__(self, obj)``

      :param np.ndarray: _description_
      :type np.ndarray: _type_
      :return: _description_
      :rtype: _type_

      .. py:method:: all(axis=None, out=None, keepdims=False)


      .. py:method:: any(axis=None, out=None, keepdims=False)


      .. py:method:: argmax(axis=None, out=None)


      .. py:method:: argmin(axis=None, out=None)


      .. py:method:: argpartition(kth, axis=-1, kind='introselect', order=None)


      .. py:method:: argsort(axis=-1, kind='quicksort', order=None)


      .. py:method:: astype(dtype, order='K', casting='unsafe', subok=True, copy=True)


      .. py:method:: byteswap(inplace=False)


      .. py:method:: choose(choices, out=None, mode='raise')


      .. py:method:: clip(min=None, max=None, out=None)


      .. py:method:: compress(condition, axis=None, out=None)


      .. py:method:: conj()


      .. py:method:: conjugate()


      .. py:method:: copy(order='C')


      .. py:method:: cumprod(axis=None, dtype=None, out=None)


      .. py:method:: cumsum(axis=None, dtype=None, out=None)


      .. py:method:: diagonal(offset=0, axis1=0, axis2=1)


      .. py:method:: dot(b, out=None)


      .. py:method:: dump(file)


      .. py:method:: dumps()


      .. py:method:: fill(value)


      .. py:method:: flatten(order='C')


      .. py:method:: getfield(dtype, offset=0)


      .. py:method:: item(*args)


      .. py:method:: itemset(*args)


      .. py:method:: max(axis=None, out=None)


      .. py:method:: mean(axis=None, dtype=None, out=None, keepdims=False)


      .. py:method:: min(axis=None, out=None, keepdims=False)


      .. py:method:: newbyteorder(new_order='S')


      .. py:method:: nonzero()


      .. py:method:: partition(kth, axis=-1, kind='introselect', order=None)


      .. py:method:: prod(axis=None, dtype=None, out=None, keepdims=False)


      .. py:method:: ptp(axis=None, out=None)


      .. py:method:: put(indices, values, mode='raise')


      .. py:method:: ravel(order='C')


      .. py:method:: repeat(repeats, axis=None)


      .. py:method:: reshape(shape, order='C')


      .. py:method:: resize(new_shape, refcheck=True)


      .. py:method:: round(decimals=0, out=None)


      .. py:method:: searchsorted(v, side='left', sorter=None)


      .. py:method:: setfield(val, dtype, offset=0)


      .. py:method:: setflags(write=None, align=None, uic=None)


      .. py:method:: sort(axis=-1, kind='quicksort', order=None)


      .. py:method:: squeeze(axis=None)


      .. py:method:: std(axis=None, dtype=None, out=None, ddof=0, keepdims=False)


      .. py:method:: sum(axis=None, dtype=None, out=None, keepdims=False)


      .. py:method:: swapaxes(axis1, axis2)


      .. py:method:: take(indices, axis=None, out=None, mode='raise')


      .. py:method:: tobytes(order='C')


      .. py:method:: tofile(fid, sep='', format='%s')


      .. py:method:: tolist()


      .. py:method:: tostring(order='C')


      .. py:method:: trace(offset=0, axis1=0, axis2=1, dtype=None, out=None)


      .. py:method:: transpose(*axes)


      .. py:method:: var(axis=None, dtype=None, out=None, ddof=0, keepdims=False)


      .. py:method:: view(dtype=None, type=None)



   .. py:method:: close()


   .. py:method:: write_to_average(day_average=5)

      Create average file (e.g., N-day averages [N=5 usually])


   .. py:method:: write_to_daily()

      Create daily file (continuous time series)


   .. py:method:: write_to_diurn(day_average=5)

      Create diurn file (variables organized by time of day & binned
      (typically 5-day bins)


   .. py:method:: write_to_fixed()

      Create ``fixed`` file (all static variables)



.. py:class:: Ncdf(filename=None, description_txt='', action='w', ncformat='NETCDF4_CLASSIC')


   Bases: :py:obj:`object`

   netCDF wrapper for archiving data in netCDF format. Alex Kling.

   Usage::

       from netcdf_wrapper import Ncdf

       Fgeo = 0.03 # W/m2, a constant
       sfcT = np.ones((24,8)) # surface temperature

       # Create file
       filename = "/path/to/myfile.nc"
       description = "results from new simulation, Alex 01-01-19"
       Log = Ncdf(filename, description)

       # Save the constant (``Fgeo``) to the file
       Log.add_constant('Fgeo', Fgeo, "geothermal flux", "W/m2")

       # Save the sfcT array to the file
       Log.add_dimension('Nx', 8)
       Log.add_dimension('time', 24)

       Log.log_variable('sfcT', sfcT, ('time', 'Nx'),
                        'soil temperature', 'K')

       Log.close()

   :param object: _description_
   :type object: _type_
   :return: netCDF file

   .. py:method:: add_constant(variable_name, value, longname_txt='', units_txt='')


   .. py:method:: add_dim_with_content(dimension_name, DATAin, longname_txt='', units_txt='', cart_txt='')

      Function to define a dimension and add a variable at the
      same time. Equivalent to ``add_dimension()`` followed by
      ``log_axis1D()``::

          lon_array = np.linspace(0, 360)

      EX::

          Log.add_dim_with_content("lon", lon_array, "longitudes",
                                   "degree", "X")


   .. py:method:: add_dimension(dimension_name, length)


   .. py:method:: close()


   .. py:method:: copy_Ncaxis_with_content(Ncdim_var)

      Copy a netCDF DIMENSION variable (e.g.,
      ``Ncdim = f.variables["lon"]``). If the dimension does not exist
      yet, it will be created


   .. py:method:: copy_Ncvar(Ncvar, swap_array=None)

      Copy a netCDF variable from another file (e.g.,
      ``Ncvar = f.variables["ucomp"]``). All dimensions must already
      exist. If ``swap_array`` is provided, the original values are
      swapped with this array.


   .. py:method:: copy_all_dims_from_Ncfile(Ncfile_in, exclude_dim=[], time_unlimited=True)

      Copy all variables, dimensions, and attributes from another
      netCDF file


   .. py:method:: copy_all_vars_from_Ncfile(Ncfile_in, exclude_var=[])


   .. py:method:: log_axis1D(variable_name, DATAin, dim_name, longname_txt='', units_txt='', cart_txt='')

      EX::

          Log.log_axis1D("areo", areo, "time", "degree", "T")


   .. py:method:: log_variable(variable_name, DATAin, dim_array, longname_txt='', units_txt='')

      EX::

          Log.log_variable("sfcT", sfcT, ("time", "Nx"),
                           "soil temperature", "K")


   .. py:method:: merge_files_from_list(Ncfilename_list, exclude_var=[])


   .. py:method:: print_dimensions()


   .. py:method:: print_variables()



