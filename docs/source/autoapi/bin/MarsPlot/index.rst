:py:mod:`bin.MarsPlot`
======================

.. py:module:: bin.MarsPlot

.. autoapi-nested-parse::

   The MarsPlot executable is for generating plots from Custom.in template
   files. It sources variables from netCDF files in a specified directory.

   The executable requires:

       * ``[-template --generate_template]`` Generates a Custom.in template
       * ``[-i --inspect]``         Triggers ncdump-like text to console
       * ``[Custom.in]``            To create plots in Custom.in template

   Third-party Requirements:

       * ``numpy``
       * ``netCDF4``
       * ``sys``
       * ``argparse``
       * ``os``
       * ``warnings``
       * ``subprocess``
       * ``matplotlib``



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   bin.MarsPlot.CustomTicker
   bin.MarsPlot.Fig_1D
   bin.MarsPlot.Fig_2D
   bin.MarsPlot.Fig_2D_lat_lev
   bin.MarsPlot.Fig_2D_lon_lat
   bin.MarsPlot.Fig_2D_lon_lev
   bin.MarsPlot.Fig_2D_lon_time
   bin.MarsPlot.Fig_2D_time_lat
   bin.MarsPlot.Fig_2D_time_lev



Functions
~~~~~~~~~

.. autoapisummary::

   bin.MarsPlot.MY_func
   bin.MarsPlot.clean_comma_whitespace
   bin.MarsPlot.create_exec
   bin.MarsPlot.create_name
   bin.MarsPlot.fig_layout
   bin.MarsPlot.filter_input
   bin.MarsPlot.format_lon_lat
   bin.MarsPlot.get_Ncdf_num
   bin.MarsPlot.get_figure_header
   bin.MarsPlot.get_lat_index
   bin.MarsPlot.get_level_index
   bin.MarsPlot.get_list_varfull
   bin.MarsPlot.get_lon_index
   bin.MarsPlot.get_overwrite_dim_1D
   bin.MarsPlot.get_overwrite_dim_2D
   bin.MarsPlot.get_time_index
   bin.MarsPlot.get_tod_index
   bin.MarsPlot.give_permission
   bin.MarsPlot.main
   bin.MarsPlot.make_template
   bin.MarsPlot.mean_func
   bin.MarsPlot.namelist_parser
   bin.MarsPlot.prep_file
   bin.MarsPlot.progress
   bin.MarsPlot.rT
   bin.MarsPlot.read_axis_options
   bin.MarsPlot.remove_whitespace
   bin.MarsPlot.select_range
   bin.MarsPlot.shift_data
   bin.MarsPlot.split_varfull



Attributes
~~~~~~~~~~

.. autoapisummary::

   bin.MarsPlot.add_sol_time_axis
   bin.MarsPlot.args
   bin.MarsPlot.current_version
   bin.MarsPlot.degr
   bin.MarsPlot.include_NaNs
   bin.MarsPlot.lon_coord_type
   bin.MarsPlot.namespace
   bin.MarsPlot.parser


.. py:class:: CustomTicker(base=10.0, labelOnlyBase=False, minor_thresholds=None, linthresh=None)


   Bases: :py:obj:`matplotlib.ticker.LogFormatterSciNotation`

   Format values following scientific notation in a logarithmic axis.

   .. py:attribute:: axis

      

   .. py:attribute:: locs
      :value: []

      

   .. py:method:: __call__(x, pos=None)

      Return the format for tick value *x* at position pos.
      ``pos=None`` indicates an unspecified location.


   .. py:method:: create_dummy_axis(**kwargs)


   .. py:method:: fix_minus(s)
      :staticmethod:

      Some classes may want to replace a hyphen for minus with the proper
      Unicode symbol (U+2212) for typographical correctness.  This is a
      helper method to perform such a replacement when it is enabled via
      :rc:`axes.unicode_minus`.


   .. py:method:: format_data(value)

      Return the full string representation of the value with the
      position unspecified.


   .. py:method:: format_data_short(value)

      Return a short string version of the tick value.

      Defaults to the position-independent long value.


   .. py:method:: format_ticks(values)

      Return the tick labels for all the ticks at once.


   .. py:method:: get_offset()


   .. py:method:: set_axis(axis)


   .. py:method:: set_base(base)

      Change the *base* for labeling.

      .. warning::
         Should always match the base used for :class:`LogLocator`


   .. py:method:: set_label_minor(labelOnlyBase)

      Switch minor tick labeling on or off.

      Parameters
      ----------
      labelOnlyBase : bool
          If True, label ticks only at integer powers of base.


   .. py:method:: set_locs(locs=None)

      Use axis view limits to control which ticks are labeled.

      The *locs* parameter is ignored in the present algorithm.



.. py:class:: Fig_1D(varfull='atmos_average.ts', doPlot=True)


   Bases: :py:obj:`object`

   .. py:method:: data_loader_1D(varfull, plot_type)


   .. py:method:: do_plot()


   .. py:method:: exception_handler(e, ax)


   .. py:method:: fig_init()


   .. py:method:: fig_save()


   .. py:method:: get_plot_type()

      Note that the ``self.t == "AXIS" test`` and the
      ``self.t = -88888`` assignment are only used when MarsPlot is
      not passed a template.

      :return: type of 1D plot to create (1D_time, 1D_lat, etc.)



   .. py:method:: make_template()


   .. py:method:: read_NCDF_1D(var_name, file_type, simuID, sol_array, plot_type, t_req, lat_req, lon_req, lev_req, ftod_req)

      Parse a Main Variable expression object that includes a square
      bracket [] (for variable calculations) for the variable to
      plot.

      :param var_name: variable name (e.g., ``temp``)
      :type var_name: str

      :param file_type: MGCM output file type. Must be ``fixed`` or
          ``average``
      :type file_type: str

      :param simuID: number identifier for netCDF file directory
      :type simuID: str

      :param sol_array: sol if different from default
          (e.g., ``02400``)
      :type sol_array:  str

      :param plot_type: ``1D_lon``, ``1D_lat``, ``1D_lev``, or
          ``1D_time``
      :type plot_type: str

      :param t_req: Ls requested
      :type t_req: str

      :param lat_req: lat requested
      :type lat_req: str

      :param lon_req: lon requested
      :type lon_req: str

      :param lev_req: level [Pa/m] requested
      :type lev_req: str

      :param ftod_req: time of day requested
      :type ftod_req: str

      :return: (dim_array) the axis (e.g., an array of longitudes),
               (var_array) the variable extracted



   .. py:method:: read_template()



.. py:class:: Fig_2D(varfull='fileYYY.XXX', doPlot=False, varfull2=None)


   Bases: :py:obj:`object`

   .. py:method:: data_loader_2D(varfull, plot_type)


   .. py:method:: exception_handler(e, ax)


   .. py:method:: fig_init()


   .. py:method:: fig_save()


   .. py:method:: filled_contour(xdata, ydata, var)


   .. py:method:: make_colorbar(levs)


   .. py:method:: make_template(plot_txt, fdim1_txt, fdim2_txt, Xaxis_txt, Yaxis_txt)


   .. py:method:: make_title(var_info, xlabel, ylabel)


   .. py:method:: plot_dimensions()


   .. py:method:: read_NCDF_2D(var_name, file_type, simuID, sol_array, plot_type, fdim1, fdim2, ftod)


   .. py:method:: read_template()


   .. py:method:: return_norm_levs()


   .. py:method:: solid_contour(xdata, ydata, var, contours)



.. py:class:: Fig_2D_lat_lev(varfull='fileYYY.XXX', doPlot=False, varfull2=None)


   Bases: :py:obj:`Fig_2D`

   .. py:method:: data_loader_2D(varfull, plot_type)


   .. py:method:: do_plot()


   .. py:method:: exception_handler(e, ax)


   .. py:method:: fig_init()


   .. py:method:: fig_save()


   .. py:method:: filled_contour(xdata, ydata, var)


   .. py:method:: make_colorbar(levs)


   .. py:method:: make_template()


   .. py:method:: make_title(var_info, xlabel, ylabel)


   .. py:method:: plot_dimensions()


   .. py:method:: read_NCDF_2D(var_name, file_type, simuID, sol_array, plot_type, fdim1, fdim2, ftod)


   .. py:method:: read_template()


   .. py:method:: return_norm_levs()


   .. py:method:: solid_contour(xdata, ydata, var, contours)



.. py:class:: Fig_2D_lon_lat(varfull='fileYYY.XXX', doPlot=False, varfull2=None)


   Bases: :py:obj:`Fig_2D`

   .. py:method:: data_loader_2D(varfull, plot_type)


   .. py:method:: do_plot()


   .. py:method:: exception_handler(e, ax)


   .. py:method:: fig_init()


   .. py:method:: fig_save()


   .. py:method:: filled_contour(xdata, ydata, var)


   .. py:method:: get_topo_2D(varfull, plot_type)

      This function returns the longitude, latitude, and topography
      to overlay as contours in a ``2D_lon_lat`` plot. Because the
      main variable requested may be complex
      (e.g., ``[00668.atmos_average_psdt2.temp]/1000.``), we will
      ensure to load the matching topography (here ``00668.fixed.nc``
      from the 2nd simulation). This function essentially does a
      simple task in a complicated way. Note that a great deal of
      the code is borrowed from the ``data_loader_2D()`` function.

      :param varfull: variable input to main_variable in Custom.in
          (e.g., ``03340.atmos_average.ucomp``)
      :type varfull: str

      :param plot_type: plot type (e.g.,
          ``Plot 2D lon X time``)
      :type plot_type: str

      :return: topography or ``None`` if no matching ``fixed`` file



   .. py:method:: make_colorbar(levs)


   .. py:method:: make_template()


   .. py:method:: make_title(var_info, xlabel, ylabel)


   .. py:method:: plot_dimensions()


   .. py:method:: read_NCDF_2D(var_name, file_type, simuID, sol_array, plot_type, fdim1, fdim2, ftod)


   .. py:method:: read_template()


   .. py:method:: return_norm_levs()


   .. py:method:: solid_contour(xdata, ydata, var, contours)



.. py:class:: Fig_2D_lon_lev(varfull='fileYYY.XXX', doPlot=False, varfull2=None)


   Bases: :py:obj:`Fig_2D`

   .. py:method:: data_loader_2D(varfull, plot_type)


   .. py:method:: do_plot()

      Create figure



   .. py:method:: exception_handler(e, ax)


   .. py:method:: fig_init()


   .. py:method:: fig_save()


   .. py:method:: filled_contour(xdata, ydata, var)


   .. py:method:: make_colorbar(levs)


   .. py:method:: make_template()

      Calls method from parent class



   .. py:method:: make_title(var_info, xlabel, ylabel)


   .. py:method:: plot_dimensions()


   .. py:method:: read_NCDF_2D(var_name, file_type, simuID, sol_array, plot_type, fdim1, fdim2, ftod)


   .. py:method:: read_template()


   .. py:method:: return_norm_levs()


   .. py:method:: solid_contour(xdata, ydata, var, contours)



.. py:class:: Fig_2D_lon_time(varfull='fileYYY.XXX', doPlot=False, varfull2=None)


   Bases: :py:obj:`Fig_2D`

   .. py:method:: data_loader_2D(varfull, plot_type)


   .. py:method:: do_plot()


   .. py:method:: exception_handler(e, ax)


   .. py:method:: fig_init()


   .. py:method:: fig_save()


   .. py:method:: filled_contour(xdata, ydata, var)


   .. py:method:: make_colorbar(levs)


   .. py:method:: make_template()


   .. py:method:: make_title(var_info, xlabel, ylabel)


   .. py:method:: plot_dimensions()


   .. py:method:: read_NCDF_2D(var_name, file_type, simuID, sol_array, plot_type, fdim1, fdim2, ftod)


   .. py:method:: read_template()


   .. py:method:: return_norm_levs()


   .. py:method:: solid_contour(xdata, ydata, var, contours)



.. py:class:: Fig_2D_time_lat(varfull='fileYYY.XXX', doPlot=False, varfull2=None)


   Bases: :py:obj:`Fig_2D`

   .. py:method:: data_loader_2D(varfull, plot_type)


   .. py:method:: do_plot()


   .. py:method:: exception_handler(e, ax)


   .. py:method:: fig_init()


   .. py:method:: fig_save()


   .. py:method:: filled_contour(xdata, ydata, var)


   .. py:method:: make_colorbar(levs)


   .. py:method:: make_template()


   .. py:method:: make_title(var_info, xlabel, ylabel)


   .. py:method:: plot_dimensions()


   .. py:method:: read_NCDF_2D(var_name, file_type, simuID, sol_array, plot_type, fdim1, fdim2, ftod)


   .. py:method:: read_template()


   .. py:method:: return_norm_levs()


   .. py:method:: solid_contour(xdata, ydata, var, contours)



.. py:class:: Fig_2D_time_lev(varfull='fileYYY.XXX', doPlot=False, varfull2=None)


   Bases: :py:obj:`Fig_2D`

   .. py:method:: data_loader_2D(varfull, plot_type)


   .. py:method:: do_plot()


   .. py:method:: exception_handler(e, ax)


   .. py:method:: fig_init()


   .. py:method:: fig_save()


   .. py:method:: filled_contour(xdata, ydata, var)


   .. py:method:: make_colorbar(levs)


   .. py:method:: make_template()


   .. py:method:: make_title(var_info, xlabel, ylabel)


   .. py:method:: plot_dimensions()


   .. py:method:: read_NCDF_2D(var_name, file_type, simuID, sol_array, plot_type, fdim1, fdim2, ftod)


   .. py:method:: read_template()


   .. py:method:: return_norm_levs()


   .. py:method:: solid_contour(xdata, ydata, var, contours)



.. py:function:: MY_func(Ls_cont)

   Returns the Mars Year

   :param Ls_cont: solar longitude (``areo``; continuous)
   :type Ls_cont: array [areo]

   :return: the Mars year
   :rtype: int



.. py:function:: clean_comma_whitespace(raw_input)

   Remove commas and whitespaces inside an expression.

   :param raw_input: dimensions specified by user input to Variable
       (e.g., ``lat=3. , lon=2 , lev = 10.``)
   :type raw_input: str

   :return: raw_input without whitespaces (e.g.,
       ``lat=3.,lon=2,lev=10.``)
   :rtype: str



.. py:function:: create_exec(raw_input, varfull_list)


.. py:function:: create_name(root_name)

   Modify file name if a file with that name already exists.

   :param root_name: path + default name for the file type (e.g.,
       ``/path/custom.in`` or ``/path/figure.png``)
   :type root_name: str

   :return: the modified name if the file already exists
       (e.g., ``/path/custom_01.in`` or ``/path/figure_01.png``)
   :rtype: str



.. py:function:: fig_layout(subID, nPan, vertical_page=False)

   Return figure layout.

   :param subID: current subplot number
   :type subID: int

   :param nPan: number of panels desired on page (max = 64, 8x8)
   :type nPan: int

   :param vertical_page: reverse the tuple for portrait format if
       ``True``
   :type vertical_page: bool

   :return: plot layout (e.g., ``plt.subplot(nrows = out[0], ncols =
       out[1], plot_number = out[2])``)
   :rtype: tuple



.. py:function:: filter_input(txt, typeIn='char')

   Read template for the type of data expected

   :param txt: text input into ``Custom.in`` to the right of an equal
       sign
   :type txt: str

   :param typeIn: type of data expected: ``char``, ``float``, ``int``,
       ``bool``, defaults to ``char``
   :type typeIn: str, optional

   :return: text input reformatted to ``[val1, val2]``
   :rtype: float or array



.. py:function:: format_lon_lat(lon_lat, type)

   Format latitude and longitude as labels (e.g., 30°S, 30°N, 45°W,
   45°E)

   :param lon_lat: latitude or longitude (+180/-180)
   :type lon_lat: float

   :param type: ``lat`` or ``lon``
   :type type: str

   :return: formatted label
   :rtype: str



.. py:function:: get_Ncdf_num()

   Return the prefix numbers for the netCDF files in the directory.
   Requires at least one ``fixed`` file in the directory.

   :return: a sorted array of sols
   :rtype: array



.. py:function:: get_figure_header(line_txt)

   Returns the plot type by confirming that template = ``True``.

   :param line_txt: template header from Custom.in (e.g.,
       ``<<<<<<<<<| Plot 2D lon X lat = True |>>>>>>>>``)
   :type line_txt: str

   :return: (figtype) figure type (e.g., ``Plot 2D lon X lat``)
   :rtype: str

   :return: (boolPlot) whether to plot (``True``) or skip (``False``)
       figure
   :rtype: bool



.. py:function:: get_lat_index(lat_query, lats)

   Returns the indices that will extract data from the netCDF file
   according to a range of *latitudes*.

   :param lat_query: requested latitudes (-90/+90)
   :type lat_query: list

   :param lats: latitude
   :type lats: array [lat]

   :return: 1d array of file indices
   :rtype: text descriptor for the extracted longitudes

   .. note::T
       The keyword ``all`` passed as ``-99999`` by the ``rt()``
       function



.. py:function:: get_level_index(level_query, levs)

   Returns the indices that will extract data from the netCDF file
   according to a range of *pressures* (resp. depth for ``zgrid``).

   :param level_query: requested pressure [Pa] (depth [m])
   :type level_query: float

   :param levs: levels (in the native coordinates)
   :type levs: array [lev]

   :return: file indices
   :rtype: array

   :return: descriptor for the extracted pressure (depth)
   :rtype: str

   .. note::
       The keyword ``all`` is passed as ``-99999`` by the ``rT()``
       functions



.. py:function:: get_list_varfull(raw_input)

   Return requested variable from a complex ``varfull`` object with ``[]``.

   :param raw_input: complex user input to Variable (e.g.,
       ``2*[atmos_average.temp]+[atmos_average2.ucomp]*1000``)
   :type raw_input: str

   :return: list required variables (e.g., [``atmos_average.temp``,
       ``atmos_average2.ucomp``])
   :rtype: str



.. py:function:: get_lon_index(lon_query_180, lons)

   Returns the indices that will extract data from the netCDF file
   according to a range of *longitudes*.

   :param lon_query_180: longitudes in -180/180: value,
       ``[min, max]``, or `None`
   :type lon_query_180: list

   :param lons: longitude in 0-360
   :type lons: array [lon]

   :return: 1D array of file indices
   :rtype: array

   :return: text descriptor for the extracted longitudes
   :rtype: str

   .. note::
       The keyword ``all`` passed as ``-99999`` by the rT() functions



.. py:function:: get_overwrite_dim_1D(varfull_bracket, t_in, lat_in, lon_in, lev_in, ftod_in)

   Return new dimensions that will overwrite default dimensions for a
   varfull object with ``{}`` for a 1D plot.

   :param varfull_bracket: a ``varfull`` object with ``{}`` (e.g.,
       ``atmos_average.temp{lev=10;ls=350;lon=155;lat=25}``)
   :type varfull_bracket: str

   :param t_in: self.t variable
   :type t_in: array [time]

   :param lat_in: self.lat variable
   :type lat_in: array [lat]

   :param lon_in: self.lon variable
   :type lon_in: array [lon]

   :param lev_in: self.lev variable
   :type lev_in: array [lev]

   :param ftod_in: self.ftod variable
   :type ftod_in: array [tod]

   :return: ``varfull`` object without brackets (e.g.,
       ``atmos_average.temp``);
       :return: (t_out) dimension to update;
       :return: (lat_out) dimension to update;
       :return: (lon_out) dimension to update;
       :return: (lev_out) dimension to update;
       :return: (ftod_out) dimension to update;



.. py:function:: get_overwrite_dim_2D(varfull_bracket, plot_type, fdim1, fdim2, ftod)

   Return new dimensions that will overwrite default dimensions for a
   varfull object with ``{}`` on a 2D plot.

   ``2D_lon_lat:  fdim1 = ls,  fdim2 = lev``
   ``2D_lat_lev:  fdim1 = ls,  fdim2 = lon``
   ``2D_time_lat: fdim1 = lon, fdim2 = lev``
   ``2D_lon_lev:  fdim1 = ls,  fdim2 = lat``
   ``2D_time_lev: fdim1 = lat, fdim2 = lon``
   ``2D_lon_time: fdim1 = lat, fdim2 = lev``

   :param varfull_bracket: a ``varfull`` object with ``{}`` (e.g.,
       ``atmos_average.temp{lev=10;ls=350;lon=155;lat=25}``)
   :type varfull_bracket: str

   :param plot_type: the type of the plot template
   :type plot_type: str

   :param fdim1: X axis dimension for plot
   :type fdim1: str

   :param fdim2: Y axis dimension for plot
   :type fdim2: str

   :return: (varfull) required file and variable (e.g.,
       ``atmos_average.temp``);
       (fdim_out1) X axis dimension for plot;
       (fdim_out2) Y axis dimension for plot;
       (ftod_out) if X or Y axis dimension is time of day



.. py:function:: get_time_index(Ls_query_360, LsDay)

   Returns the indices that will extract data from the netCDF file
   according to a range of solar longitudes [0-360].

   First try the Mars Year of the last timestep, then try the year
   before that. Use whichever Ls period is closest to the requested
   date.

   :param Ls_query_360: requested solar longitudes
   :type Ls_query_360: list

   :param LsDay: continuous solar longitudes
   :type LsDay: array [areo]

   :return: file indices
   :rtype: array

   :return: descriptor for the extracted solar longitudes
   :rtype: str

   .. note::
       The keyword ``all`` is passed as ``-99999`` by the ``rT()``
       function



.. py:function:: get_tod_index(tod_query, tods)

   Returns the indices that will extract data from the netCDF file
   according to a range of *times of day*.

   :param tod_query: requested time of day (0-24)
   :type tod_query: list

   :param tods: times of day
   :type tods: array [tod]

   :return: file indices
   :rtype: array [tod]

   :return: descriptor for the extracted time of day
   :rtype: str

   .. note::
       The keyword ``all`` is passed as ``-99999`` by the ``rT()``
       function



.. py:function:: give_permission(filename)

   Sets group permissions for files created on NAS.

   :param filename: name of the file
   :type filename: str



.. py:function:: main()


.. py:function:: make_template()

   Generate the ``Custom.in`` template file.

   :return: Custom.in blank template



.. py:function:: mean_func(arr, axis)

   This function calculates a mean over the selected axis, ignoring or
   including NaN values as specified by ``show_NaN_in_slice`` in
   ``amescap_profile``.

   :param arr: the array to be averaged
   :type arr: array

   :param axis: the axis over which to average the array
   :type axis: int

   :return: the mean over the time axis



.. py:function:: namelist_parser(Custom_file)

   Parse a ``Custom.in`` template.

   :param Custom_file: full path to ``Custom.in`` file
   :type Custom_file: str

   :return: updated global variables, ``FigLayout``, ``objectList``



.. py:function:: prep_file(var_name, file_type, simuID, sol_array)

   Open the file as a Dataset or MFDataset object depending on its
       status on Lou. Note that the input arguments are typically
       extracted from a ``varfull`` object (e.g.,
       ``03340.atmos_average.ucomp``) and not from a file whose disk
       status is known beforehand.

   :param var_name: variable to extract (e.g., ``ucomp``)
   :type var_name: str

   :param file_type: MGCM output file type (e.g., ``average``)
   :type file_name: str

   :param simuID: simulation ID number (e.g., 2 for 2nd simulation)
   :type simuID: int

   :param sol_array: date in file name (e.g., [3340,4008])
   :type sol_array: list

   :return: Dataset or MFDataset object;
       (var_info) longname and units;
       (dim_info) dimensions e.g., (``time``, ``lat``,``lon``);
       (dims) shape of the array e.g., [133,48,96]



.. py:function:: progress(k, Nmax, txt='', success=True)

   Display a progress bar when performing heavy calculations.

   :param k: current iteration of the outer loop
   :type k: float

   :param Nmax: max iteration of the outer loop
   :type Nmax: float

   :return: progress bar (EX: ``Running... [#---------] 10.64 %``)



.. py:function:: rT(typeIn='char')

   Read template for the type of data expected. Returns value to
   ``filter_input()``.

   :param typeIn: type of data expected: ``char``, ``float``, ``int``,
       ``bool``, defaults to ``char``
   :type typeIn: str, optional

   :return: text input reformatted to ``[val1, val2]``
   :rtype: float or array



.. py:function:: read_axis_options(axis_options_txt)

   Return axis customization options.

   :param axis_options_txt: a copy of the last line ``Axis Options``
       in ``Custom.in`` templates
   :type axis_options_txt: str

   :return: X-axis bounds as a numpy array or ``None`` if undedefined
   :rtype: array or None

   :return: Y-axis bounds as a numpy array or ``None`` if undedefined
   :rtype: array or None

   :return: colormap (e.g., ``jet``, ``nipy_spectral``) or line
       options (e.g., ``--r`` for dashed red)
   :rtype: str

   :return: linear (``lin``) or logarithmic (``log``) color scale
   :rtype: str

   :return: projection (e.g., ``ortho -125,45``)
   :rtype: str



.. py:function:: remove_whitespace(raw_input)

   Remove whitespace inside an expression.

   This is different from the ``.strip()`` method, which only removes
   whitespaces at the edges of a string.

   :param raw_input: user input for variable, (e.g.,
       ``[atmos_average.temp] + 2)``
   :type raw_input: str

   :return: raw_input without whitespaces (e.g.,
       ``[atmos_average.temp]+2)``
   :rtype: str



.. py:function:: select_range(Ncdf_num, bound)

   Return the prefix numbers for the netCDF files in the directory
   within the user-defined range.

   :param Ncdf_num: a sorted array of sols
   :type Ncdf_num: array

   :param bound: a sol (e.g., 0350) or range of sols ``[min max]``
   :type bound: int or array

   :return: a sorted array of sols within the bounds
   :rtype: array



.. py:function:: shift_data(lon, data)

   Shifts the longitude data from 0-360 to -180/180 and vice versa.

   :param lon: 1D array of longitude
   :type lon: array [lon]

   :param data: 2D array with last dimension = longitude
   :type data: array [1,lon]

   :raises ValueError: Longitude coordinate type is not recognized.

   :return: longitude (-180/180)
   :rtype: array [lon]

   :return: shifted data
   :rtype: array [1,lon]

   .. note::
       Use ``np.ma.hstack`` instead of ``np.hstack`` to keep the
       masked array properties



.. py:function:: split_varfull(varfull)

   Split ``varfull`` object into its component parts

   :param varfull: a ``varfull`` object (e.g,
       ``atmos_average@2.zsurf``, ``02400.atmos_average@2.zsurf``)
   :type varfull: str

   :return: (sol_array) a sol number or ``None`` (if none provided)
   :rtype: int or None

   :return: (filetype) file type (e.g, ``atmos_average``)
   :rtype: str

   :return: (var) variable of interest (e.g, ``zsurf``)
   :rtype: str

   :return: (``simuID``) simulation ID (Python indexing starts at 0)
   :rtype: int



.. py:data:: add_sol_time_axis

   

.. py:data:: args

   

.. py:data:: current_version
   :value: 3.5

   

.. py:data:: degr
   :value: '°'

   

.. py:data:: include_NaNs

   

.. py:data:: lon_coord_type

   

.. py:data:: namespace

   

.. py:data:: parser

   

