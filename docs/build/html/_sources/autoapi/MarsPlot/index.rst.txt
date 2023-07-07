:orphan:

:py:mod:`MarsPlot`
==================

.. py:module:: MarsPlot


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   MarsPlot.CustomTicker



Functions
~~~~~~~~~

.. autoapisummary::

   MarsPlot.shift_data
   MarsPlot.MY_func
   MarsPlot.get_lon_index
   MarsPlot.get_lat_index
   MarsPlot.get_tod_index
   MarsPlot.get_level_index
   MarsPlot.get_time_index
   MarsPlot.filter_input
   MarsPlot.rT
   MarsPlot.read_axis_options
   MarsPlot.split_varfull
   MarsPlot.remove_whitespace
   MarsPlot.clean_comma_whitespace
   MarsPlot.get_list_varfull
   MarsPlot.get_overwrite_dim_2D
   MarsPlot.get_overwrite_dim_1D
   MarsPlot.fig_layout
   MarsPlot.namelist_parser
   MarsPlot.get_figure_header
   MarsPlot.format_lon_lat
   MarsPlot.get_Ncdf_num
   MarsPlot.select_range
   MarsPlot.create_name
   MarsPlot.path_to_template
   MarsPlot.progress
   MarsPlot.prep_file



.. py:function:: shift_data(lon, data)

   This function shift the longitude and data from a 0->360 to a -180/+180 grid.
   Args:
       lon: 1D array of longitude 0->360
       data: 2D array with last dimension being the longitude
   Returns:
       lon: 1D array of longitude -180/+180
       data: shifted data
   Note: Use np.ma.hstack instead of np.hstack to keep the masked array properties


.. py:function:: MY_func(Ls_cont)

   This function return the Mars Year
   Args:
       Ls_cont: solar longitude, contineuous
   Returns:
       MY : int the Mars year


.. py:function:: get_lon_index(lon_query_180, lons)

   Given a range of requested longitudes, return the indexes to extract data from the netcdf file
   Args:
       lon_query_180: requested longitudes in -180/+180 units: value, [min, max] or None
       lons:          1D array of longitude in the native coordinates (0->360)
   Returns:
       loni: 1D array of file indexes
       txt_lon: text descriptor for the extracted longitudes
   *** Note that the keyword 'all' is passed as -99999 by the rT() functions


.. py:function:: get_lat_index(lat_query, lats)

   Given a range of requested latitudes, return the indexes to extract data from the netcdf file
   Args:
       lat_query: requested latitudes -90/+90
       lats:      1D array of latitudes in the native coordinates
   Returns:
       lati: 1D array of file indexes
       txt_lat: text descriptor for the extracted latitudes
   *** Note that the keyword 'all' is passed as -99999 by the rT() functions


.. py:function:: get_tod_index(tod_query, tods)

   Given a range of requested time of day, return the indexes to extract data from the netcdf file
   Args:
       tod_query: requested time of day, 0-24
       tods:      1D array of time of days in the native coordinates
   Returns:
       todi: 1D array of file indexes
       txt_tod: text descriptor for the extracted time of day
   *** Note that the keyword 'all' is passed as -99999 by the rT() functions


.. py:function:: get_level_index(level_query, levs)

   Given a range of requested pressures (resp. depth for 'zgrid'), return the indexes to extract data from the netcdf file
   Args:
       level_query: requested  pressure in [Pa] (resp. depth in [m])
       levs:         1D array of levels in the native coordinates [Pa] (resp. [m])
   Returns:
       levi: 1D array of file indexes
       txt_lev: text descriptor for the extracted pressure (resp. depth)
   *** Note that the keyword 'all' is passed as -99999 by the rT() functions


.. py:function:: get_time_index(Ls_query_360, Ls)

   Given a range of requested solar longitude [0-360], return the indexes to extract data from the netcdf file.
   First try the Mars year of the last timestep, then try the year before then picks whichever Ls period is closest to the requested date.

   Args:
       Ls_query_360: requested  solar longitudes
       Ls_c:         1D array of continueous solar longitudes
   Returns:
       ti: 1D array of file indexes
       txt_time: text descriptor for the extracted solar longitudes
   *** Note that the keyword 'all' is passed as -99999 by the rT() functions


.. py:function:: filter_input(txt, typeIn='char')

   Read Template for the type of data expected
   Args:
       txt: a string, typical the right-hand sign of an equal sign '3', '3,4', or 'all'
       typeIn: type of data expected: 'char', 'float', 'int', 'bool'
   Returns:
       out: float or 1D array [val1,val2] in the expected format



.. py:function:: rT(typeIn='char')

   Read Template for the type of data expected
   Args:
       typeIn: type of data expected: 'char', 'float', 'int', 'bool'
   Returns:
       out: float or 1D array [val1,val2] in the expected format



.. py:function:: read_axis_options(axis_options_txt)

   Return axis customization options
   Args:
       axis_options_txt: One liner string: 'Axis Options  : lon = [5,8] | lat = [None,None] | cmap = jet | scale= lin | proj = cart'
   Returns:
       Xaxis: X-axis bounds as a numpy array or None if undedefined
       Yaxis: Y-axis bounds as a numpy array or None if undedefined
       custom_line1: string, i.e colormap ('jet', 'nipy_spectral') or line options, e.g '--r' for dashed red
       custom_line2: lin or log
       custom_line3: None of string for projections, e.g 'ortho -125,45'



.. py:function:: split_varfull(varfull)

   Split  the varfull object into its different components.
   Args:
       varfull: a varfull object, for example 'atmos_average@2.zsurf',
                                              '02400.atmos_average@2.zsurf'
   Returns:
       sol_array: a sol number e.g 2400 or None if none is provided
       filetype:  file type, i.e 'atmos_average'
       var:       variable of interest, i.e 'zsurf'
       simuID:    int, simulation ID = 2-1= 1 as Python indexes start at zero



.. py:function:: remove_whitespace(raw_input)

   Remove the white space inside an expression. This is different from the '.strip()' method that only remove white spaces at the edges of the string
   Args:
       raw_input: a string, e.g '[atmos_average.temp] +  2'
   Returns:
       processed_input  the string without white spaces, e.g [atmos_average.temp]+2'



.. py:function:: clean_comma_whitespace(raw_input)

   Remove the commas and white spaces inside an expression.
   Args:
       raw_input: a string, e.g 'lat=3. ,'
   Returns:
       processed_input  the string without white spaces or commas e.g 'lat=3.lon=2lev=10.'



.. py:function:: get_list_varfull(raw_input)

   Given an expression object with '[]' return the different variable needed
   Args:
       raw_input: a complex varfull object, for example '2*[atmos_average.temp]+[atmos_average2.ucomp]*1000'
   Returns:
       var_list  a list of variable to load, e.g ['atmos_average.temp', 'atmos_average2.ucomp']



.. py:function:: get_overwrite_dim_2D(varfull_bracket, plot_type, fdim1, fdim2, ftod)

   Given a single varfull object with '{}' return the new dimensions to overwrite the default dimensions
   Args:
       varfull_bracket: a  varfull object with any of the following atmos_average.temp{lev=10;ls=350;lon=155;lat=25} (brackets and semi-colons separated)
       plot_type: the type of plot

   Returns:
       varfull: the varfull without brackets: e.g 'atmos_average.temp'
       fdim_out1,fdim_out1,ftod_out: the dimensions to update
   NOTE:
   2D_lon_lat:   fdim1=ls
                 fdim2=lev

   2D_lat_lev: fdim1=ls
                 fdim2=lon

   2D_time_lat: fdim1=lon
                 fdim2=lev

   2D_lon_lev: fdim1=ls
                 fdim2=lat

   2D_time_lev:fdim1=lat
                 fdim2=lon

   2D_lon_time:  fdim1=lat
                 fdim2=lev


.. py:function:: get_overwrite_dim_1D(varfull_bracket, t_in, lat_in, lon_in, lev_in, ftod_in)

   Given a single varfull object with '{}' return the new dimensions to overwrite the default dimensions
   Args:
       varfull_bracket: a  varfull object with any of the following atmos_average.temp{lev=10;ls=350;lon=155;lat=25;tod=15}
       t_in,lat_in,lon_in,lev_in,ftod_in: the variables as defined by self.t ,self.lat,self.lon,self.lev,self.ftod

   Returns:
       varfull the varfull without brackets: e.g 'atmos_average.temp'
       t_out,lat_out,lon_out,lev_out,ftod_out: the dimensions to update
   NOTE:



.. py:function:: fig_layout(subID, nPan, vertical_page=False)

   Return figure layout
   Args:
       subID:    integer, current subplot number
       nPan : integer, number of panels desired on the figure, up to 64 (8x8 panel)
       vertical_page: if True, reverse the tuple for vertical layout
   Returns:
       out: tuple with approriate layout: plt.subplot(nrows=out[0],ncols=out[1],plot_number=out[2])


.. py:function:: namelist_parser(Custom_file)

   Parse a template
   Args:
       Custom_file: full path to Custom.in file
   Actions:
       Update  global variableFigLayout, objectList


.. py:function:: get_figure_header(line_txt)

   This function return the type of a figure and tells us if wanted
   Args:
       line_txt: string, figure header from Custom.in, i.e '<<<<<<<<<<<<<<| Plot 2D lon X lat = True |>>>>>>>>>>>>>'
   Returns:
       figtype : string, figure type, i.e:  Plot 2D lon X lat
       boolPlot: boolean, is the plot wanted?


.. py:function:: format_lon_lat(lon_lat, type)

   Format latitude and longitude as labels, e.g. 30S , 30N, 45W, 45E
   Args:
       lon_lat (float): latitude or longitude +/-180
       type (string) : 'lat' or 'lon'
   Returns:
       lon_lat_label : (string), formatted label


.. py:function:: get_Ncdf_num()

   Get the sol numbers of all the netcdf files in directory
   This test is based on the existence of a least one  XXXXX.fixed.nc in the current directory.
   Args:
       None
   Returns:
       Ncdf_num: a sorted array of sols


.. py:function:: select_range(Ncdf_num, bound)

   Args:
       Ncdf_num:  a sorted array of sols
       bound: a integer representing a date (e.g. 0350) or an array containing the sol bounds (e.g [min max])
   Returns:
       Ncdf_num: a sorted array of sols within the prescribed bounds


.. py:function:: create_name(root_name)

   Create a file name based on its existence in the current directory.
   Args:
       root_name: desired name for the file: "/path/custom.in" or "/path/figure.png"
   Returns:
       new_name: new name if the file already exists: "/path/custom_01.in" or "/path/figure_01.png"


.. py:function:: path_to_template(custom_name)

   Create a file name based on its existence in the current directory.
   Args:
       custom_name: custom file name, accepted formats are my_custom or my_custom.in
   Returns:
       full_path: full_path to /u/user/FV3/templates/my_custom.in

        If file not found, try shared directory


.. py:function:: progress(k, Nmax, txt='', success=True)

   Display a progress bar to monitor heavy calculations.
   Args:
       k: current iteration of the outer loop
       Nmax: max iteration of the outer loop
   Returns:
       Running... [#---------] 10.64 %


.. py:function:: prep_file(var_name, file_type, simuID, sol_array)

   Given the different information, open the file as a Dataset or MFDataset object.
   Note that the input arguments are typically extracted  from a varfull object, e.g.  '03340.atmos_average.ucomp',
   not a file from those the existence on the disk is known beforehand
   Args:
       var_name: variable to extract, e.g. 'ucomp'
       file_type: 'fixed', atmos_average_pstd
       simuID:    e.g 2 for 2nd simulation
       sol_array: e.g [3340,4008]

   Returns:
       f: Dataset or MFDataset object
       var_info: longname and units
       dim_info: dimensions e.g ('time', 'lat','lon')
       dims:    shape of the array e.g [133,48,96]


.. py:class:: CustomTicker(base=10.0, labelOnlyBase=False, minor_thresholds=None, linthresh=None)


   Bases: :py:obj:`matplotlib.ticker.LogFormatterSciNotation`

   Format values following scientific notation in a logarithmic axis.

   .. py:method:: base(base)

      Change the *base* for labeling.

      .. warning::
         Should always match the base used for :class:`LogLocator`


   .. py:method:: label_minor(labelOnlyBase)

      Switch minor tick labeling on or off.

      Parameters
      ----------
      labelOnlyBase : bool
          If True, label ticks only at integer powers of base.


   .. py:method:: set_locs(locs=None)

      Use axis view limits to control which ticks are labeled.

      The *locs* parameter is ignored in the present algorithm.


   .. py:method:: format_data(value)

      Return the full string representation of the value with the
      position unspecified.


   .. py:method:: format_data_short(value)

      Return a short string version of the tick value.

      Defaults to the position-independent long value.


   .. py:method:: format_ticks(values)

      Return the tick labels for all the ticks at once.


   .. py:method:: fix_minus(s)
      :staticmethod:

      Some classes may want to replace a hyphen for minus with the proper
      unicode symbol (U+2212) for typographical correctness.  This is a
      helper method to perform such a replacement when it is enabled via
      :rc:`axes.unicode_minus`.



