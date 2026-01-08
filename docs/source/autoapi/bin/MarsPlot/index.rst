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
       * ``pypdf``



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
   bin.MarsPlot.debug_wrapper
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
   bin.MarsPlot.debug
   bin.MarsPlot.degr
   bin.MarsPlot.exit_code
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

   Fig_1D is a parent class for generating and handling 1D plots of
   Mars atmospheric data.

   Attributes:
       title : str
           Title of the plot.
       legend : str
           Legend label for the plot.
       varfull : str
           Full variable specification, including file and variable
           name.
       t : str or float
           Time axis or identifier for the varying dimension.
       lat : float or str
           Latitude value or identifier.
       lon : float or str
           Longitude value or identifier.
       lev : float or str
           Vertical level value or identifier.
       ftod : float or str
           Time of day requested.
       hour : float or str
           Hour of day, used for diurnal plots.
       doPlot : bool
           Whether to generate the plot.
       plot_type : str
           Type of 1D plot (e.g., "1D_time", "1D_lat").
       sol_array : str
           Sol array extracted from varfull.
       filetype : str
           File type extracted from varfull.
       var : str
           Variable name extracted from varfull.
       simuID : str
           Simulation ID extracted from varfull.
       nPan : int
           Number of panels in the plot.
       subID : int
           Subplot ID.
       addLine : bool
           Whether to add a line to an existing plot.
       layout : list or None
           Page layout for multipanel plots.
       fdim_txt : str
           Annotation for free dimensions.
       success : bool
           Indicates if the plot was successfully created.
       vert_unit : str
           Vertical unit, either "m" or "Pa".
       Dlim : list or None
           Dimension limits for the axis.
       Vlim : list or None
           Variable limits for the axis.
       axis_opt1 : str
           Line style or axis option.
       axis_opt2 : str
           Additional axis option (optional).

   Methods:
       make_template():
           Writes a template for the plot configuration to a file.
       read_template():
           Reads plot configuration from a template file.
       get_plot_type():
           Determines the type of 1D plot to create based on which
           dimension is set to "AXIS" or -88888.
       data_loader_1D(varfull, plot_type):
           Loads 1D data for plotting, handling variable expressions
           and dimension overwrites.
       read_NCDF_1D(var_name, file_type, simuID, sol_array, plot_type,
       t_req, lat_req, lon_req, lev_req, ftod_req):
           Reads and processes 1D data from a NetCDF file for the
           specified variable and dimensions.
       exception_handler(e, ax):
           Handles exceptions during plotting, displaying an error
           message on the plot.
       fig_init():
           Initializes the figure and subplot for plotting.
       fig_save():
           Saves the generated figure to disk.
       do_plot():
           Main method to generate the 1D plot, handling all plotting
           logic and exceptions.
           

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
      :type  var_name: str
      :param file_type: MGCM output file type. Must be ``fixed`` or
          ``average``
      :type  file_type: str
      :param simuID: number identifier for netCDF file directory
      :type  simuID: str
      :param sol_array: sol if different from default
          (e.g., ``02400``)
      :type  sol_array:  str
      :param plot_type: ``1D_lon``, ``1D_lat``, ``1D_lev``, or
          ``1D_time``
      :type  plot_type: str
      :param t_req: Ls requested
      :type  t_req: str
      :param lat_req: lat requested
      :type  lat_req: str
      :param lon_req: lon requested
      :type  lon_req: str
      :param lev_req: level [Pa/m] requested
      :type  lev_req: str
      :param ftod_req: time of day requested
      :type  ftod_req: str
      :return: (dim_array) the axis (e.g., an array of longitudes),
               (var_array) the variable extracted


   .. py:method:: read_template()



.. py:class:: Fig_2D(varfull='fileYYY.XXX', doPlot=False, varfull2=None)


   Bases: :py:obj:`object`

   Base class for 2D figures. This class is not intended to be
   instantiated directly. Instead, it is used as a base class for
   specific 2D figure classes (e.g., ``Fig_2D_lon_lat``, ``Fig_2D_time_lat``,
   ``Fig_2D_lat_lev``, etc.). It provides common attributes and methods
   for all 2D figures, such as the variable name, file type, simulation
   ID, and plotting options. The class also includes methods for
   creating a template for the figure, reading the template from a
   file, and loading data for 2D plots. The class is designed to be
   extended by subclasses that implement specific plotting
   functionality for different types of 2D figures.

   :param varfull: full variable name (e.g., ``fileYYY.XXX``)
   :type  varfull: str
   :param doPlot: whether to plot the figure (default: ``False``)
   :type  doPlot: bool
   :param varfull2: second variable name (default: ``None``)
   :type  varfull2: str
   :return: None
   :rtype:  None
   :raises ValueError: If the input varfull is not a valid type
       for variable name.
   :raises TypeError: If the input doPlot is not a valid type
       for plotting.
   :raises Exception: If the input varfull2 is not a valid type
       for variable name.

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

   A subclass of Fig_2D for generating 2D plots with latitude and
   vertical level (pressure or altitude) axes.

   This class customizes the plotting template and plotting logic for
   visualizing data as a function of latitude and vertical level.
   It supports filled contour plots for a primary variable, and
   optionally overlays solid contour lines for a secondary variable.

   Methods:
       make_template():
           Sets up the plot template with appropriate titles and axis
           labels for latitude vs. level plots.

       do_plot():
           Loads data, creates a filled contour plot of the primary
           variable, optionally overlays contours of a secondary
           variable, configures axis scaling and formatting (including
           logarithmic pressure axis if needed), sets axis limits and
           tick formatting, handles exceptions, and saves the resulting
           figure.

   Attributes (inherited and/or used):
       varfull : str
           Name of the primary variable to plot.
       varfull2 : str or None
           Name of the secondary variable to overlay as contours, if
           any.
       plot_type : str
           Type of plot or data selection.
       vert_unit : str
           Unit for the vertical axis ("Pa" for pressure, otherwise
           altitude in meters).
       Xlim : tuple or None
           Limits for the x-axis (latitude).
       Ylim : tuple or None
           Limits for the y-axis (level).
       contour2 : list or None
           Contour levels for the secondary variable.
       nPan : int
           Number of panels in the plot (affects tick label size).
       success : bool
           Indicates if the plot was successfully created.

   .. py:method:: data_loader_2D(varfull, plot_type)


   .. py:method:: do_plot()

      Generates a 2D latitude-level plot for the specified variable(s).
      This method initializes the figure, loads the required data, and creates a filled contour plot
      of the primary variable. If a secondary variable is specified, it overlays solid contours for
      that variable. The y-axis is set to logarithmic scale and inverted if the vertical unit is pressure.
      Axis limits, labels, and tick formatting are applied as specified by the instance attributes.
      The plot title is generated based on the variable information. Handles exceptions during plotting
      and saves the resulting figure.

      Raises:
          Exception: Any exception encountered during plotting is handled and logged.


   .. py:method:: exception_handler(e, ax)


   .. py:method:: fig_init()


   .. py:method:: fig_save()


   .. py:method:: filled_contour(xdata, ydata, var)


   .. py:method:: make_colorbar(levs)


   .. py:method:: make_template()

      Creates and configures a plot template for a 2D latitude versus level plot.
      This method calls the parent class's `make_template` method with predefined
      titles and axis labels suitable for a plot displaying latitude against atmospheric
      level data.
      The plot is labeled as "Plot 2D lat X lev" with the following axis labels:
          - X-axis: "Ls 0-360 "
          - Y-axis: "Lon +/-180"
          - Additional axes: "Lat", "Level[Pa/m]"
      Returns:
          None


   .. py:method:: make_title(var_info, xlabel, ylabel)


   .. py:method:: plot_dimensions()


   .. py:method:: read_NCDF_2D(var_name, file_type, simuID, sol_array, plot_type, fdim1, fdim2, ftod)


   .. py:method:: read_template()


   .. py:method:: return_norm_levs()


   .. py:method:: solid_contour(xdata, ydata, var, contours)



.. py:class:: Fig_2D_lon_lat(varfull='fileYYY.XXX', doPlot=False, varfull2=None)


   Bases: :py:obj:`Fig_2D`

   Fig_2D_lon_lat is a class for creating 2D longitude-latitude plots.

   Fig_2D_lon_lat is a subclass of Fig_2D designed for generating 2D
   plots of longitude versus latitude, primarily for visualizing Mars
   climate data. It provides methods for figure creation, data loading,
   plotting, and overlaying topography contours, with support for
   various map projections and customization options.

   Attributes:
       varfull (str): Full variable name (e.g., "fileYYY.XXX") to plot.
       doPlot (bool): Whether to plot the figure (default: False).
       varfull2 (str, optional): Second variable name for overlaying
           contours (default: None).
       plot_type (str): Type of plot (default: "2D_lon_lat").
       fdim1 (str, optional): First free dimension (default: None).
       fdim2 (str, optional): Second free dimension (default: None).
       ftod (str, optional): Time of day (default: None).
       axis_opt1 (str, optional): First axis option, e.g., colormap
           (default: None).
       axis_opt2 (str, optional): Second axis option (default: None).
       axis_opt3 (str, optional): Projection type (e.g., "cart",
           "robin", "moll", "Npole", "Spole", "ortho").
       Xlim (tuple, optional): Longitude axis limits.
       Ylim (tuple, optional): Latitude axis limits.
       range (bool, optional): Whether to use a specified range for
           color levels.
       contour2 (float or list, optional): Contour levels for the
           second variable.
       title (str, optional): Custom plot title.
       nPan (int): Number of panels (for multi-panel plots).
       fdim_txt (str): Text describing free dimensions.
       success (bool): Status flag indicating if plotting succeeded.

   Methods:
       make_template():
           Sets up the plot template with appropriate axis labels and
           titles.

       get_topo_2D(varfull, plot_type):
           Loads and returns topography data (zsurf) for overlaying as
           contours, matching the simulation and file type of the main variable.

       do_plot():
           Main plotting routine. Loads data, applies projection,
           overlays topography and optional second variable contours,
           customizes axes, and saves the figure. Handles both
           standard and special map projections (cartesian, Robinson,
           Mollweide, polar, orthographic).

   Usage:
       This class is intended to be used within the MarsPlot software
       for visualizing Mars climate model outputs as longitude-latitude
       maps, with optional overlays and advanced projection support.

   .. py:method:: data_loader_2D(varfull, plot_type)


   .. py:method:: do_plot()

      Generate a 2D longitude-latitude plot with various projection options and optional overlays.

      This method creates a 2D plot of a variable (and optionally a second variable as contours)
      on a longitude-latitude grid. It supports multiple map projections, including cartesian,
      Robinson, Mollweide, and azimuthal (north pole, south pole, orthographic) projections.
      Topography contours can be added if available. The method handles axis formatting,
      colorbars, titles, and annotation of meridians and parallels.

      The plotting behavior is controlled by instance attributes such as:
          - self.varfull: Main variable to plot.
          - self.varfull2: Optional second variable for contour overlay.
          - self.plot_type: Type of plot to generate.
          - self.axis_opt1: Colormap or colormap option.
          - self.axis_opt3: Projection type.
          - self.contour2: Contour levels for the second variable.
          - self.Xlim, self.Ylim: Axis limits for cartesian projection.
          - self.range: Whether to use a specific range for color levels.
          - self.title: Custom plot title.
          - self.fdim_txt: Additional dimension text for the title.
          - self.nPan: Panel index for multi-panel plots.

      The method handles exceptions and saves the figure upon completion.

      Raises:
          Exception: Any error encountered during plotting is handled and reported.


   .. py:method:: exception_handler(e, ax)


   .. py:method:: fig_init()


   .. py:method:: fig_save()


   .. py:method:: filled_contour(xdata, ydata, var)


   .. py:method:: get_topo_2D(varfull, plot_type)

      This function returns the longitude, latitude, and topography
      to overlay as contours in a ``2D_lon_lat`` plot. Because the
      main variable requested may be complex
      (e.g., ``[01336.atmos_average_psdt2.temp]/1000.``), we will
      ensure to load the matching topography (here ``01336.fixed.nc``
      from the 2nd simulation). This function essentially does a
      simple task in a complicated way. Note that a great deal of
      the code is borrowed from the ``data_loader_2D()`` function.

      :param varfull: variable input to main_variable in Custom.in
          (e.g., ``03340.atmos_average.ucomp``)
      :type  varfull: str
      :param plot_type: plot type (e.g.,
          ``Plot 2D lon X time``)
      :type  plot_type: str
      :return: topography or ``None`` if no matching ``fixed`` file


   .. py:method:: make_colorbar(levs)


   .. py:method:: make_template()

      Creates and configures a plot template for 2D longitude vs latitude data.
      This method calls the parent class's `make_template` method with predefined
      parameters to set up the plot title and axis labels specific to a 2D longitude-latitude plot.
      The template includes:
          - Title: "Plot 2D lon X lat"
          - X-axis label: "Ls 0-360"
          - Y-axis label: "Level Pa/m"
          - Additional axis labels: "Lon" (longitude), "Lat" (latitude)


   .. py:method:: make_title(var_info, xlabel, ylabel)


   .. py:method:: plot_dimensions()


   .. py:method:: read_NCDF_2D(var_name, file_type, simuID, sol_array, plot_type, fdim1, fdim2, ftod)


   .. py:method:: read_template()


   .. py:method:: return_norm_levs()


   .. py:method:: solid_contour(xdata, ydata, var, contours)



.. py:class:: Fig_2D_lon_lev(varfull='fileYYY.XXX', doPlot=False, varfull2=None)


   Bases: :py:obj:`Fig_2D`

   A subclass of Fig_2D for generating 2D plots with longitude and
   vertical level (pressure or altitude) axes.

   This class customizes the template and plotting routines to
   visualize data as a function of longitude and vertical level.
   It supports plotting filled contours for a primary variable and
   optional solid contours for a secondary variable.
   The vertical axis can be displayed in pressure (Pa, logarithmic
   scale) or altitude (m).

   Methods:
       make_template():
           Sets up the plot template with appropriate titles and axis
           labels for longitude vs. level plots.

       do_plot():
           Loads data, applies longitude shifting, creates filled and
           optional solid contour plots,
           configures axis scales and labels, and handles exceptions
           during plotting.

   .. py:method:: data_loader_2D(varfull, plot_type)


   .. py:method:: do_plot()

      Generates a 2D plot of a variable as a function of longitude
      and vertical level (pressure or altitude).

      This method initializes the figure, loads the required data,
      applies longitude shifting, and creates filled and/or solid
      contour plots.

      It handles plotting of a secondary variable if specified, sets
      axis scales and labels based on the vertical coordinate unit,
      applies axis limits if provided, customizes tick formatting and
      font sizes, and manages exceptions during plotting.
      The resulting figure is saved to file.

      Raises:
          Exception: If any error occurs during the plotting process,
          it is handled and logged by the exception handler.


   .. py:method:: exception_handler(e, ax)


   .. py:method:: fig_init()


   .. py:method:: fig_save()


   .. py:method:: filled_contour(xdata, ydata, var)


   .. py:method:: make_colorbar(levs)


   .. py:method:: make_template()

      Creates and configures a plot template for 2D lon x lev data.

      This method sets up the plot with predefined titles and axis
      labels:
      - Title: "Plot 2D lon X lev"
      - X-axis: "Ls 0-360"
      - Y-axis: "Latitude"
      - Additional labels: "Lon +/-180" and "Level[Pa/m]"

      Overrides the base class method to provide specific
      configuration for this plot type.


   .. py:method:: make_title(var_info, xlabel, ylabel)


   .. py:method:: plot_dimensions()


   .. py:method:: read_NCDF_2D(var_name, file_type, simuID, sol_array, plot_type, fdim1, fdim2, ftod)


   .. py:method:: read_template()


   .. py:method:: return_norm_levs()


   .. py:method:: solid_contour(xdata, ydata, var, contours)



.. py:class:: Fig_2D_lon_time(varfull='fileYYY.XXX', doPlot=False, varfull2=None)


   Bases: :py:obj:`Fig_2D`

   A specialized 2D plotting class for visualizing data as a function
   of longitude and time (Ls).

   Inherits from: Fig_2D

   Methods:
       make_template():
           Sets up the plot template with appropriate titles and axis
           labels for longitude vs. time plots.

       do_plot():
           Generates a 2D plot with longitude on the x-axis and solar
           longitude (Ls) on the y-axis.
           Loads and processes data, applies shifting if necessary,
           and creates filled and/or solid contours.
           Handles axis formatting, tick labeling (including optional
           sol time annotation), and plot saving.
           Catches and handles exceptions during plotting.

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

   A 2D plotting class for visualizing data as a function of time (Ls)
   and latitude. Inherits from: Fig_2D

   Methods:
       make_template():
           Sets up the plot template with appropriate titles and axis
           labels for a 2D time vs latitude plot.
       do_plot():
           Loads 2D data (time and latitude), creates a filled contour
           plot of the primary variable, and optionally overlays a
           solid contour of a secondary variable.
           Formats axes, customizes tick labels to show both Ls and
           sol time (if enabled), and applies axis limits if specified.
           Handles exceptions during plotting and saves the resulting
           figure.

   Attributes (inherited and used):
       varfull : str
           Name of the primary variable to plot.
       varfull2 : str or None
           Name of the secondary variable to overlay as contours
           (optional).
       plot_type : str
           Type of plot/data to load.
       Xlim : tuple or None
           Limits for the x-axis (sol time).
       Ylim : tuple or None
           Limits for the y-axis (latitude).
       contour2 : list or None
           Contour levels for the secondary variable.
       nPan : int
           Number of panels (used for label sizing).
       success : bool
           Indicates if the plot was successfully created.

   .. py:method:: data_loader_2D(varfull, plot_type)


   .. py:method:: do_plot()

      Generates a 2D time-latitude plot for the specified variable(s).
      This method initializes the figure, loads the required 2D data arrays (time and latitude),
      and creates a filled contour plot of the primary variable. If a secondary variable is specified,
      it overlays solid contours for that variable. The method also formats the axes, including
      custom tick labels for solar longitude (Ls) and optionally sol time, and applies axis limits
      if specified. Additional plot formatting such as tick intervals and font sizes are set.
      The plot is saved at the end of the method. Any exceptions encountered during plotting
      are handled and reported.

      Raises:
          Exception: If any error occurs during the plotting process, it is handled and reported.


   .. py:method:: exception_handler(e, ax)


   .. py:method:: fig_init()


   .. py:method:: fig_save()


   .. py:method:: filled_contour(xdata, ydata, var)


   .. py:method:: make_colorbar(levs)


   .. py:method:: make_template()

      Creates and configures a plot template for a 2D time versus latitude figure.
      This method calls the superclass's `make_template` method with predefined
      titles and axis labels suitable for a plot displaying data across longitude,
      level, solar longitude (Ls), and latitude.

      Returns:
          None


   .. py:method:: make_title(var_info, xlabel, ylabel)


   .. py:method:: plot_dimensions()


   .. py:method:: read_NCDF_2D(var_name, file_type, simuID, sol_array, plot_type, fdim1, fdim2, ftod)


   .. py:method:: read_template()


   .. py:method:: return_norm_levs()


   .. py:method:: solid_contour(xdata, ydata, var, contours)



.. py:class:: Fig_2D_time_lev(varfull='fileYYY.XXX', doPlot=False, varfull2=None)


   Bases: :py:obj:`Fig_2D`

   A specialized 2D plotting class for visualizing data as a function
   of time (Ls) and vertical level (pressure or altitude).

   Inherits from: Fig_2D

   Methods:
       make_template():
           Sets up the plot template with appropriate axis labels and
           titles for 2D time vs. level plots.

       do_plot():
           Loads data and generates a filled contour plot of the
           primary variable as a function of solar longitude (Ls) and
           vertical level. Optionally overlays a solid contour of a
           secondary variable.
           Handles axis formatting, tick labeling (including optional
           sol time axis), and y-axis scaling (logarithmic for
           pressure). Sets plot titles and saves the figure. Catches
           and handles exceptions during plotting.

   Attributes (inherited and/or used):
       varfull : str
           Name of the primary variable to plot.
       varfull2 : str or None
           Name of the secondary variable to overlay as contours
           (optional).
       plot_type : str
           Type of plot/data selection.
       Xlim : tuple or None
           Limits for the x-axis (solar day).
       Ylim : tuple or None
           Limits for the y-axis (vertical level).
       vert_unit : str
           Unit for the vertical axis ("Pa" for pressure or other for
           altitude).
       nPan : int
           Number of panels/subplots (affects label size).
       contour2 : list or None
           Contour levels for the secondary variable.
       success : bool
           Indicates if the plot was successfully generated.
       

   .. py:method:: data_loader_2D(varfull, plot_type)


   .. py:method:: do_plot()

      Generates a 2D time-level plot for Mars atmospheric data.

      This method initializes the figure, loads the required data, and creates a filled contour plot
      of the primary variable over solar longitude (Ls) and pressure or altitude. If a secondary variable
      is specified, it overlays solid contours for that variable. The method also formats axes, applies
      custom tick labels (optionally including sol time), and adjusts axis scales and labels based on
      the vertical unit (pressure or altitude). The plot is titled and saved to file.
      Handles exceptions by invoking a custom exception handler and always attempts to save the figure.

      Attributes used:
          varfull (str): Name of the primary variable to plot.
          plot_type (str): Type of plot/data to load.
          varfull2 (str, optional): Name of the secondary variable for contour overlay.
          contour2 (list, optional): Contour levels for the secondary variable.
          Xlim (tuple, optional): Limits for the x-axis (solar day).
          Ylim (tuple, optional): Limits for the y-axis (pressure or altitude).
          vert_unit (str): Vertical axis unit, either "Pa" for pressure or other for altitude.
          nPan (int): Number of panels (affects label size).
          success (bool): Set to True if plotting succeeds.

      Raises:
          Handles all exceptions internally and logs them via a custom handler.


   .. py:method:: exception_handler(e, ax)


   .. py:method:: fig_init()


   .. py:method:: fig_save()


   .. py:method:: filled_contour(xdata, ydata, var)


   .. py:method:: make_colorbar(levs)


   .. py:method:: make_template()

      Creates and configures a plot template for 2D time versus level visualization.
      This method calls the superclass's `make_template` method with predefined
      titles and axis labels suitable for plotting data with latitude, longitude,
      solar longitude (Ls), and atmospheric level (in Pa/m).

      Returns:
          None


   .. py:method:: make_title(var_info, xlabel, ylabel)


   .. py:method:: plot_dimensions()


   .. py:method:: read_NCDF_2D(var_name, file_type, simuID, sol_array, plot_type, fdim1, fdim2, ftod)


   .. py:method:: read_template()


   .. py:method:: return_norm_levs()


   .. py:method:: solid_contour(xdata, ydata, var, contours)



.. py:function:: MY_func(Ls_cont)

   Returns the Mars Year.

   :param Ls_cont: solar longitude (``areo``; continuous)
   :type  Ls_cont: array [areo]
   :return: the Mars year
   :rtype:  int
   :raises ValueError: If the input Ls_cont is not a valid type for
       year calculation.


.. py:function:: clean_comma_whitespace(raw_input)

   Remove commas and whitespaces inside an expression.

   :param raw_input: dimensions specified by user input to Variable
       (e.g., ``lat=3. , lon=2 , lev = 10.``)
   :type  raw_input: str
   :return: raw_input without whitespaces (e.g.,
       ``lat=3.,lon=2,lev=10.``)
   :rtype:  str


.. py:function:: create_exec(raw_input, varfull_list)


.. py:function:: create_name(root_name)

   Modify file name if a file with that name already exists.

   :param root_name: path + default name for the file type (e.g.,
       ``/path/custom.in`` or ``/path/figure.png``)
   :type  root_name: str
   :return: the modified name if the file already exists
       (e.g., ``/path/custom_01.in`` or ``/path/figure_01.png``)
   :rtype:  str
   :raises ValueError: If the input root_name is not a valid type
       for file name.
   :raises TypeError: If the input root_name is not a valid type
       for file name.
   :raises Exception: If the file name creation fails for any
       reason.


.. py:function:: debug_wrapper(func)

   A decorator that wraps a function with error handling
   based on the --debug flag.
   If the --debug flag is set, it prints the full traceback
   of any exception that occurs. Otherwise, it prints a
   simplified error message.

   :param func: The function to wrap.
   :type   func: function
   :return: The wrapped function.
   :rtype:  function
   :raises Exception: If an error occurs during the function call.
   :raises TypeError: If the function is not callable.
   :raises ValueError: If the function is not found.
   :raises NameError: If the function is not defined.
   :raises AttributeError: If the function does not have the
       specified attribute.
   :raises ImportError: If the function cannot be imported.
   :raises RuntimeError: If the function cannot be run.
   :raises KeyError: If the function does not have the
       specified key.
   :raises IndexError: If the function does not have the
       specified index.
   :raises IOError: If the function cannot be opened.
   :raises OSError: If the function cannot be accessed.
   :raises EOFError: If the function cannot be read.
   :raises MemoryError: If the function cannot be allocated.
   :raises OverflowError: If the function cannot be overflowed.
   :raises ZeroDivisionError: If the function cannot be divided by zero.
   :raises StopIteration: If the function cannot be stopped.
   :raises KeyboardInterrupt: If the function cannot be interrupted.
   :raises SystemExit: If the function cannot be exited.
   :raises AssertionError: If the function cannot be asserted.


.. py:function:: fig_layout(subID, nPan, vertical_page=False)

   Return figure layout.

   :param subID: current subplot number
   :type  subID: int
   :param nPan: number of panels desired on page (max = 64, 8x8)
   :type  nPan: int
   :param vertical_page: reverse the tuple for portrait format if
       ``True``
   :type  vertical_page: bool
   :return: plot layout (e.g., ``plt.subplot(nrows = out[0], ncols =
       out[1], plot_number = out[2])``)
   :rtype:  tuple
   :raises ValueError: If the input subID is not a valid type for
       subplot number.
   :raises TypeError: If the input nPan is not a valid type for
       subplot number.
   :raises Exception: If the input vertical_page is not a valid type
       for subplot number.
   :raises Exception: If the figure layout calculation fails for any
       reason.


.. py:function:: filter_input(txt, typeIn='char')

   Read template for the type of data expected.

   Returns value to ``rT()``.

   :param txt: text input into ``Custom.in`` to the right of an equal
       sign
   :type  txt: str
   :param typeIn: type of data expected: ``char``, ``float``, ``int``,
       ``bool``, defaults to ``char``
   :type  typeIn: str, optional
   :return: text input reformatted to ``[val1, val2]``
   :rtype:  float or array
   :raises ValueError: If the input txt is not a valid type for
       filtering.
   :raises TypeError: If the input typeIn is not a valid type for
       filtering.
   :raises Exception: If the filtering operation fails for any reason.


.. py:function:: format_lon_lat(lon_lat, type)

   Format latitude and longitude as labels (e.g., 30°S, 30°N, 45°W,
   45°E)

   :param lon_lat: latitude or longitude (+180/-180)
   :type  lon_lat: float
   :param type: ``lat`` or ``lon``
   :type  type: str
   :return: formatted label
   :rtype:  str
   :raises ValueError: If the input lon_lat is not a valid type for
       latitude or longitude.
   :raises TypeError: If the input type is not a valid type for
       latitude or longitude.
   :raises Exception: If the formatting fails for any reason.


.. py:function:: get_Ncdf_num()

   Return the prefix numbers for the netCDF files in the directory.
   Requires at least one ``fixed`` file in the directory.

   :return: a sorted array of sols
   :rtype:  array
   :raises ValueError: If the input input_paths is not a valid type
       for file name.


.. py:function:: get_figure_header(line_txt)

   Returns the plot type by confirming that template = ``True``.

   :param line_txt: template header from Custom.in (e.g.,
       ``<<<<<<<<<| Plot 2D lon X lat = True |>>>>>>>>``)
   :type  line_txt: str
   :return: (figtype) figure type (e.g., ``Plot 2D lon X lat``)
   :rtype:  str
   :return: (boolPlot) whether to plot (``True``) or skip (``False``)
       figure
   :rtype:  bool
   :raises ValueError: If the input line_txt is not a valid type for
       figure header.
   :raises TypeError: If the input line_txt is not a valid type for
       figure header.
   :raises Exception: If the figure header parsing fails for any
       reason.


.. py:function:: get_lat_index(lat_query, lats)

   Returns the indices for a range of latitudes in a file.

   :param lat_query: requested latitudes (-90/+90)
   :type  lat_query: list
   :param lats: latitude
   :type  lats: array [lat]
   :return: 1d array of file indices
   :rtype:  text descriptor for the extracted longitudes
   :rtype:  str
   :raises ValueError: If the input lat_query is not a valid type for
       latitude calculation.

   .. note::T
       The keyword ``all`` passed as ``-99999`` by the ``rt()``
       function


.. py:function:: get_level_index(level_query, levs)

   Returns the indices for a range of pressures in a file.

   :param level_query: requested pressure [Pa] (depth [m])
   :type  level_query: float
   :param levs: levels (in the native coordinates)
   :type  levs: array [lev]
   :return: file indices
   :rtype:  array
   :return: descriptor for the extracted pressure (depth)
   :rtype:  str
   :raises ValueError: If the input level_query is not a valid type for
       level calculation.

   .. note::
       The keyword ``all`` is passed as ``-99999`` by the ``rT()``
       functions


.. py:function:: get_list_varfull(raw_input)

   Return requested variable from a complex ``varfull`` object with ``[]``.

   :param raw_input: complex user input to Variable (e.g.,
       ``2*[atmos_average.temp]+[atmos_average2.ucomp]*1000``)
   :type  raw_input: str
   :return: list required variables (e.g., [``atmos_average.temp``,
       ``atmos_average2.ucomp``])
   :rtype:  str
   :raises ValueError: If the input raw_input is not a valid type for
       variable extraction.


.. py:function:: get_lon_index(lon_query_180, lons)

   Returns the indices for a range of longitudes in a file.

   :param lon_query_180: longitudes in -180/180: value,
       ``[min, max]``, or `None`
   :type  lon_query_180: list
   :param lons: longitude in 0-360
   :type  lons: array [lon]
   :return: 1D array of file indices
   :rtype:  array
   :return: text descriptor for the extracted longitudes
   :rtype:  str
   :raises ValueError: If the input lon_query_180 is not a valid type
       for longitude calculation.

   .. note::
       The keyword ``all`` passed as ``-99999`` by the rT() functions


.. py:function:: get_overwrite_dim_1D(varfull_bracket, t_in, lat_in, lon_in, lev_in, ftod_in)

   1D plot: overwrite dimensions in ``varfull`` object with ``{}``.

   (e.g., ``atmos_average.temp{lev=10;ls=350;lon=155;lat=25}``)
   This function is used to overwrite the default dimensions in a
   ``varfull`` object with ``{}`` (e.g., ``atmos_average.temp{lev=10;
   ls=350;lon=155;lat=25}``) for a 1D plot. The function will return
   the new dimensions that will overwrite the default dimensions for
   the ``varfull`` object. The function will also return the required
   file and variable (e.g., ``atmos_average.temp``) and the X and Y
   axis dimensions for the plot.

   :param varfull_bracket: a ``varfull`` object with ``{}`` (e.g.,
       ``atmos_average.temp{lev=10;ls=350;lon=155;lat=25}``)
   :type  varfull_bracket: str
   :param t_in: self.t variable
   :type  t_in: array [time]
   :param lat_in: self.lat variable
   :type  lat_in: array [lat]
   :param lon_in: self.lon variable
   :type  lon_in: array [lon]
   :param lev_in: self.lev variable
   :type  lev_in: array [lev]
   :param ftod_in: self.ftod variable
   :type  ftod_in: array [tod]
   :return: ``varfull`` object without brackets (e.g.,
       ``atmos_average.temp``);
   :return: (t_out) dimension to update;
   :return: (lat_out) dimension to update;
   :return: (lon_out) dimension to update;
   :return: (lev_out) dimension to update;
   :return: (ftod_out) dimension to update;
   :rtype:  str
   :raises ValueError: If the input varfull_bracket is not a valid
       type for variable extraction.
   :raises TypeError: If the input t_in, lat_in, lon_in, lev_in,
       ftod_in are not valid types for variable extraction.
   :raises Exception: If the variable extraction fails for any reason.

   .. note:: This function is used for 1D plots only. The function
       will return the new dimensions that will overwrite the default
       dimensions for the ``varfull`` object. The function will also
       return the required file and variable (e.g.,
       ``atmos_average.temp``) and the X and Y axis dimensions for the
       plot.


.. py:function:: get_overwrite_dim_2D(varfull_bracket, plot_type, fdim1, fdim2, ftod)

   2D plot: overwrite dimensions in ``varfull`` object with ``{}``.

   (e.g., ``atmos_average.temp{lev=10;ls=350;lon=155;lat=25}``)

   This function is used to overwrite the default dimensions in a
   ``varfull`` object with ``{}`` (e.g., ``atmos_average.temp{lev=10;
   ls=350;lon=155;lat=25}``) for a 2D plot. The function will return
   the new dimensions that will overwrite the default dimensions for
   the ``varfull`` object. The function will also return the required
   file and variable (e.g., ``atmos_average.temp``) and the X and Y
   axis dimensions for the plot.

   ``2D_lon_lat:  fdim1 = ls,  fdim2 = lev``
   ``2D_lat_lev:  fdim1 = ls,  fdim2 = lon``
   ``2D_time_lat: fdim1 = lon, fdim2 = lev``
   ``2D_lon_lev:  fdim1 = ls,  fdim2 = lat``
   ``2D_time_lev: fdim1 = lat, fdim2 = lon``
   ``2D_lon_time: fdim1 = lat, fdim2 = lev``

   :param varfull_bracket: a ``varfull`` object with ``{}`` (e.g.,
       ``atmos_average.temp{lev=10;ls=350;lon=155;lat=25}``)
   :type  varfull_bracket: str
   :param plot_type: the type of the plot template
   :type  plot_type: str
   :param fdim1: X axis dimension for plot
   :type  fdim1: str
   :param fdim2: Y axis dimension for plot
   :type  fdim2: str
   :return: (varfull) required file and variable (e.g.,
       ``atmos_average.temp``);
       (fdim_out1) X axis dimension for plot;
       (fdim_out2) Y axis dimension for plot;
       (ftod_out) if X or Y axis dimension is time of day
   :rtype:  str
   :raises ValueError: If the input varfull_bracket is not a valid
       type for variable extraction.
   :raises TypeError: If the input plot_type is not a valid type for
       variable extraction.
   :raises Exception: If the variable extraction fails for any reason.


.. py:function:: get_time_index(Ls_query_360, LsDay)

   Returns the indices for a range of solar longitudes in a file.

   First try the Mars Year of the last timestep, then try the year
   before that. Use whichever Ls period is closest to the requested
   date.

   :param Ls_query_360: requested solar longitudes
   :type  Ls_query_360: list
   :param LsDay: continuous solar longitudes
   :type  LsDay: array [areo]
   :return: file indices
   :rtype:  array
   :return: descriptor for the extracted solar longitudes
   :rtype:  str
   :raises ValueError: If the input Ls_query_360 is not a valid type
       for solar longitude calculation.
   :raises TypeError: If the input LsDay is not a valid type for
       solar longitude calculation.
   :raises Exception: If the time index calculation fails for any
       reason.

   .. note::
       The keyword ``all`` is passed as ``-99999`` by the ``rT()``
       function


.. py:function:: get_tod_index(tod_query, tods)

   Returns the indices for a range of times of day in a file.

   :param tod_query: requested time of day (0-24)
   :type  tod_query: list
   :param tods: times of day
   :type  tods: array [tod]
   :return: file indices
   :rtype:  array [tod]
   :return: descriptor for the extracted time of day
   :rtype:  str
   :raises ValueError: If the input tod_query is not a valid type for
       time of day calculation.

   .. note::
       The keyword ``all`` is passed as ``-99999`` by the ``rT()``
       function


.. py:function:: give_permission(filename)

   Sets group permissions for files created on NAS.

   :param filename: name of the file
   :type  filename: str
   :raises ValueError: If the input filename is not a valid type
       for file name.


.. py:function:: main()

   Main entry point for the MarsPlot script.

   Handles argument parsing, global variable setup, figure object
   initialization, and execution of the main plotting workflow.
   Depending on the provided arguments, this function can:

       - Inspect the contents of a NetCDF file and print variable
       information or statistics.
       - Generate a template configuration file.
       - Parse a provided template file, select data based on optional
       date bounds, and generate
           diagnostic plots as individual files or as a merged
           multipage PDF.
       - Manage output directories and file naming conventions.
       - Display progress and handle debug output.

   Global variables are set for configuration and figure formatting.
   The function also manages error handling and user feedback for
   invalid arguments or file operations.


.. py:function:: make_template()

   Generate the ``Custom.in`` template file.

   :return: Custom.in blank template
   :rtype:  file
   :raises ValueError: If the input customFileIN is not a valid type
       for template generation.
   :raises TypeError: If the input customFileIN is not a valid type
       for template generation.
   :raises Exception: If the template generation fails for any
       reason.


.. py:function:: mean_func(arr, axis)

   Calculate the mean of an array along a specified axis.

   This function calculates a mean over the selected axis, ignoring or
   including NaN values as specified by ``show_NaN_in_slice`` in
   ``amescap_profile``.

   :param arr: the array to be averaged
   :type  arr: array
   :param axis: the axis over which to average the array
   :type  axis: int
   :return: the mean over the time axis
   :rtype:  array
   :raises ValueError: If the array is empty or the axis is out of bounds.
   :raises RuntimeWarning: If the mean calculation encounters NaN values.
   :raises TypeError: If the input array is not a valid type for mean
       calculation.
   :raises Exception: If the mean calculation fails for any reason.


.. py:function:: namelist_parser(Custom_file)

   Parse a ``Custom.in`` template.

   :param Custom_file: full path to ``Custom.in`` file
   :type  Custom_file: str
   :return: updated global variables, ``FigLayout``, ``objectList``
       ``panelList``, ``subplotList``, ``addLineList``, ``layoutList``
   :rtype:  list
   :raises ValueError: If the input Custom_file is not a valid type
       for file name.


.. py:function:: prep_file(var_name, file_type, simuID, sol_array)

   Open the file as a Dataset or MFDataset object depending on its
   status on Lou. Note that the input arguments are typically
   extracted from a ``varfull`` object (e.g.,
   ``03340.atmos_average.ucomp``) and not from a file whose disk status
   is known beforehand.

   :param var_name: variable to extract (e.g., ``ucomp``)
   :type  var_name: str
   :param file_type: MGCM output file type (e.g., ``average``)
   :type  file_name: str
   :param simuID: simulation ID number (e.g., 2 for 2nd simulation)
   :type  simuID: int
   :param sol_array: date in file name (e.g., [3340,4008])
   :type  sol_array: list
   :return: Dataset or MFDataset object;
   :return: (var_info) longname and units;
   :return: (dim_info) dimensions e.g., (``time``, ``lat``,``lon``);
   :return: (dims) shape of the array e.g., [133,48,96]
   :rtype:  Dataset or MFDataset object, str, tuple, list
   :raises ValueError: If the input var_name is not a valid type
       for variable name.
   :raises TypeError: If the input file_type is not a valid type
       for file type.
   :raises Exception: If the file preparation fails for any
       reason.
   :raises IOError: If the file is not found or cannot be opened.


.. py:function:: progress(k, Nmax, txt='', success=True)

   Display a progress bar when performing heavy calculations.

   :param k: current iteration of the outer loop
   :type  k: float
   :param Nmax: max iteration of the outer loop
   :type  Nmax: float
   :return: progress bar (EX: ``Running... [#---------] 10.64 %``)
   :rtype:  str
   :raises ValueError: If the input k is not a valid type for
       progress bar.
   :raises TypeError: If the input Nmax is not a valid type for
       progress bar.
   :raises Exception: If the progress bar creation fails for any
       reason.


.. py:function:: rT(typeIn='char')

   Read template for the type of data expected.

   Returns value to
   ``filter_input()``.

   :param typeIn: type of data expected: ``char``, ``float``, ``int``,
       ``bool``, defaults to ``char``
   :type  typeIn: str, optional
   :return: text input reformatted to ``[val1, val2]``
   :rtype:  float or array
   :raises ValueError: If the input typeIn is not a valid type for
       filtering.
   :raises TypeError: If the input typeIn is not a valid type for
       filtering.
   :raises Exception: If the filtering operation fails for any reason.


.. py:function:: read_axis_options(axis_options_txt)

   Return axis customization options.

   :param axis_options_txt: a copy of the last line ``Axis Options``
       in ``Custom.in`` templates
   :type  axis_options_txt: str
   :return: X-axis bounds as a numpy array or ``None`` if undedefined
   :rtype:  array or None
   :return: Y-axis bounds as a numpy array or ``None`` if undedefined
   :rtype:  array or None
   :return: colormap (e.g., ``jet``, ``nipy_spectral``) or line
       options (e.g., ``--r`` for dashed red)
   :rtype:  str
   :return: linear (``lin``) or logarithmic (``log``) color scale
   :rtype:  str
   :return: projection (e.g., ``ortho -125,45``)
   :rtype:  str
   :raises ValueError: If the input axis_options_txt is not a valid
       type for axis options.


.. py:function:: remove_whitespace(raw_input)

   Remove whitespace inside an expression.

   This is different from the ``.strip()`` method, which only removes
   whitespaces at the edges of a string.

   :param raw_input: user input for variable, (e.g.,
       ``[atmos_average.temp] + 2)``
   :type  raw_input: str
   :return: raw_input without whitespaces (e.g.,
       ``[atmos_average.temp]+2)``
   :rtype:  str
   :raises ValueError: If the input raw_input is not a valid type for
       whitespace removal.


.. py:function:: select_range(Ncdf_num, bound)

   Return the prefix numbers for the netCDF files in the directory
   within the user-defined range.

   :param Ncdf_num: a sorted array of sols
   :type  Ncdf_num: array
   :param bound: a sol (e.g., 0350) or range of sols ``[min max]``
   :type  bound: int or array
   :return: a sorted array of sols within the bounds
   :rtype:  array
   :raises ValueError: If the input Ncdf_num is not a valid type for
       file name.
   :raises TypeError: If the input bound is not a valid type for
       file name.
   :raises Exception: If the range selection fails for any reason.


.. py:function:: shift_data(lon, data)

   Shifts the longitude data from 0-360 to -180/180 and vice versa.

   :param lon: 1D array of longitude
   :type  lon: array [lon]
   :param data: 2D array with last dimension = longitude
   :type  data: array [1,lon]
   :return: 1D array of longitude in -180/180 or 0-360
   :rtype:  array [lon]
   :return: 2D array with last dimension = longitude
   :rtype:  array [1,lon]
   :raises ValueError: If the longitude coordinate type is invalid.
   :raises TypeError: If the input data is not a valid type for
       shifting.
   :raises Exception: If the shifting operation fails for any reason.


   .. note::
       Use ``np.ma.hstack`` instead of ``np.hstack`` to keep the
       masked array properties


.. py:function:: split_varfull(varfull)

   Split ``varfull`` object into its component parts

   :param varfull: a ``varfull`` object (e.g,
       ``atmos_average@2.zsurf``, ``02400.atmos_average@2.zsurf``)
   :type  varfull: str
   :return: (sol_array) a sol number or ``None`` (if none provided)
   :rtype:  int or None
   :return: (filetype) file type (e.g, ``atmos_average``)
   :rtype:  str
   :return: (var) variable of interest (e.g, ``zsurf``)
   :rtype:  str
   :return: (``simuID``) simulation ID (Python indexing starts at 0)
   :rtype:  int
   :raises ValueError: If the input varfull is not a valid type for
       splitting.


.. py:data:: add_sol_time_axis

   

.. py:data:: args

   

.. py:data:: current_version
   :value: 3.5

   

.. py:data:: debug

   

.. py:data:: degr
   :value: '°'

   

.. py:data:: exit_code

   

.. py:data:: include_NaNs

   

.. py:data:: lon_coord_type

   

.. py:data:: namespace

   

.. py:data:: parser

   

