# !/usr/bin/env python3
"""
Ncdf_wrapper archives data into netCDF format. It serves as a wrapper
for creating netCDF files.

Third-party Requirements:

    * ``numpy``
    * ``amescap.FV3_utils``
    * ``scipy.io``
    * ``netCDF4``
    * ``os``
    * ``datetime``
"""

# Load generic Python modules
import numpy as np
from netCDF4 import Dataset, MFDataset
from scipy.io import FortranFile
from amescap.FV3_utils import daily_to_average, daily_to_diurn
import os
import datetime

# ======================================================================
#                           DEFINITIONS
# ======================================================================

class Ncdf(object):
    """
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
    """

    def __init__(self, filename=None, description_txt="", action="w",
                 ncformat="NETCDF4_CLASSIC"):
        if filename:
            if filename[-3:] != ".nc":
            # Assume that only path is provided. Create file name
                now = datetime.datetime.now()
                filename = (f"{filename}/run_{now.day:02}-{now.month:02}-"
                            f"{now.year:04}_{now.hour}-{now.minute}-"
                            f"{now.second}.nc")
        else:
            # Create a default file name if path/filename not provided
            pathname = f"{os.getcwd()}/"
            now = datetime.datetime.now()
            filename = (f"{pathname}run_{now.day:02}-{now.month:02}-"
                      f"{now.year:04}_{now.hour}-{now.minute}-{now.second}.nc")
        self.filename = filename
        if action=="w":
            self.f_Ncdf = Dataset(filename, "w", format = ncformat)
            self.f_Ncdf.description = description_txt
        elif action=="a":
            # Append to file
            self.f_Ncdf = Dataset(filename, "a", format = ncformat)

        # Create dictionaries to hold dimensions and variables
        self.dim_dict = dict()
        self.var_dict = dict()
        # print(f"{filename} was created")


    def close(self):
        self.f_Ncdf.close()
        print(f"{self.filename} was created")


    def add_dimension(self, dimension_name, length):
        self.dim_dict[dimension_name] = (
            self.f_Ncdf.createDimension(dimension_name, length))


    def print_dimensions(self):
        print(self.dim_dict.items())


    def print_variables(self):
        print(self.var_dict.keys())


    def add_constant(self, variable_name, value, longname_txt="",
                     units_txt=""):
        if "constant" not in self.dim_dict.keys():
            self.add_dimension("constant", 1)

        # Add the value to the longname
        longname_txt = f"{longname_txt} ({value})"
        self._def_variable(variable_name, ("constant"), longname_txt,
                           units_txt)
        self.var_dict[variable_name][:] = value


    # Private Definitions
    def _def_variable(self, variable_name, dim_array, longname_txt="",
                      units_txt="",datatype="f4"):
        self.var_dict[variable_name] = self.f_Ncdf.createVariable(variable_name,
                                                                  datatype,
                                                                  dim_array)
        self.var_dict[variable_name].units = units_txt
        self.var_dict[variable_name].long_name = longname_txt
        self.var_dict[variable_name].dim_name = str(dim_array)


    def _def_axis1D(self, variable_name, dim_array, longname_txt="",
                    units_txt="", cart_txt="",datatype="f8"):
        self.var_dict[variable_name] = self.f_Ncdf.createVariable(variable_name,
                                                                 datatype,
                                                                 dim_array)
        self.var_dict[variable_name].units = units_txt
        self.var_dict[variable_name].long_name = longname_txt
        self.var_dict[variable_name].cartesian_axis = cart_txt


    def _test_var_dimensions(self, Ncvar):
        all_dim_OK = True
        for s in Ncvar.dimensions:
            if s not in self.dim_dict.keys():
                print(f"***Warning***, dimension '{Ncvar._name}' not yet "
                      f"defined, skipping variable")
                all_dim_OK = False
        return all_dim_OK


    def _is_cart_axis(self, Ncvar):
        """
        The cartesian axis attribute is replicated if ``cartesian_axis``
        is in the original variable and if the dimension = 1. This will
        exclude FV3 T-cell latitudes ``grid_xt_bnds`` and
        ``grid_yt_bnds``, which are of size ``(lon, bnds)`` &
        ``(lat, bnds)`` with dimension = 2.
        """

        cart_axis = False
        tmp_cart = getattr(Ncvar, 'cartesian_axis', False)
        tmp_size = getattr(Ncvar, 'dimensions')

        if tmp_cart and len(tmp_size)==1:
            cart_axis = True
        return cart_axis


    def log_variable(self, variable_name, DATAin, dim_array, longname_txt="",
                     units_txt="",datatype="f4"):
        """
        EX::

            Log.log_variable("sfcT", sfcT, ("time", "Nx"),
                             "soil temperature", "K")
        """
        if variable_name not in self.var_dict.keys():
            self._def_variable(variable_name, dim_array, longname_txt,
                               units_txt,datatype)
        self.var_dict[variable_name].long_name = longname_txt
        self.var_dict[variable_name].dim_name = str(dim_array)
        self.var_dict[variable_name].units = units_txt
        self.var_dict[variable_name][:] = DATAin


    def log_axis1D(self, variable_name, DATAin, dim_name, longname_txt="",
                   units_txt="", cart_txt=""):
        """
        EX::

            Log.log_axis1D("areo", areo, "time", "degree", "T")
        """
        if variable_name not in self.var_dict.keys():
            self._def_axis1D(variable_name, dim_name, longname_txt, units_txt,
                             cart_txt)
        self.var_dict[variable_name].long_name = longname_txt
        self.var_dict[variable_name].units = units_txt
        self.var_dict[variable_name].cartesian_axis = cart_txt
        self.var_dict[variable_name][:] = DATAin


    def add_dim_with_content(self, dimension_name, DATAin, longname_txt="",
                             units_txt="", cart_txt=''):
        """
        Function to define a dimension and add a variable at the
        same time. Equivalent to ``add_dimension()`` followed by
        ``log_axis1D()``::

            lon_array = np.linspace(0, 360)

        EX::

            Log.add_dim_with_content("lon", lon_array, "longitudes",
                                     "degree", "X")
        """

        if dimension_name not in self.dim_dict.keys():
            self.add_dimension(dimension_name, len(DATAin))

        if longname_txt == "":
            # If no longname provided, use ``dimension_name``
            longname_txt = dimension_name

        if dimension_name not in self.var_dict.keys():
            self._def_axis1D(dimension_name, dimension_name, longname_txt,
                             units_txt, cart_txt)
        self.var_dict[dimension_name].long_name = longname_txt
        self.var_dict[dimension_name].units = units_txt
        self.var_dict[dimension_name].cartesian_axis = cart_txt
        self.var_dict[dimension_name][:] = DATAin

    # .. note:: The attribute ``name``  was replaced by ``_name`` for
    # compatibility with MFDataset:
    # When using ``f=MFDataset(fname,"r")``, ``f.variables[var]`` does
    # not have a ``name`` attribute but does have ``_name``


    def copy_Ncaxis_with_content(self, Ncdim_var):
        """
        Copy a netCDF DIMENSION variable (e.g.,
        ``Ncdim = f.variables["lon"]``). If the dimension does not exist
        yet, it will be created
        """

        longname_txt = getattr(Ncdim_var, "long_name", Ncdim_var._name)
        units_txt = getattr(Ncdim_var, "units", "")
        cart_txt = getattr(Ncdim_var, "cartesian_axis", "")
        self.add_dim_with_content(Ncdim_var._name, Ncdim_var[:], longname_txt,
                                  units_txt, cart_txt)


    def copy_Ncvar(self, Ncvar, swap_array=None):
        """
        Copy a netCDF variable from another file (e.g.,
        ``Ncvar = f.variables["ucomp"]``). All dimensions must already
        exist. If ``swap_array`` is provided, the original values are
        swapped with this array.
        """

        if Ncvar._name not in self.var_dict.keys():
            dim_array = Ncvar.dimensions
            longname_txt = getattr(Ncvar, "long_name", Ncvar._name)
            units_txt = getattr(Ncvar, "units", "")
            self._def_variable(Ncvar._name, Ncvar.dimensions, longname_txt,
                               units_txt,Ncvar.dtype)
            if np.any(swap_array):
                self.log_variable(Ncvar._name, swap_array[:], Ncvar.dimensions,
                                  longname_txt, units_txt)
            else:
                self.log_variable(Ncvar._name, Ncvar[:], Ncvar.dimensions,
                                  longname_txt, units_txt)
        else:
            print(f"***Warning***, '{Ncvar._name}' is already defined, "
                  f"skipping it")


    def copy_all_dims_from_Ncfile(self, Ncfile_in, exclude_dim=[],
                                  time_unlimited=True):
        """
        Copy all variables, dimensions, and attributes from another
        netCDF file
        """

        # First include dimensions
        all_dims = Ncfile_in.dimensions.keys()
        for idim in all_dims:
            if idim not in exclude_dim:
                if idim == 'time' and time_unlimited:
                    self.add_dimension(Ncfile_in.dimensions[idim]._name, None)
                else:
                    self.add_dimension(Ncfile_in.dimensions[idim]._name,
                                       Ncfile_in.dimensions[idim].size)


    def copy_all_vars_from_Ncfile(self, Ncfile_in, exclude_var=[]):
        # First include variables
        all_vars = Ncfile_in.variables.keys()
        for ivar in all_vars:
            if ivar not in exclude_var:
                if self._test_var_dimensions(Ncfile_in.variables[ivar]):
                    # Skip variables if not all dimensions are available
                    if self._is_cart_axis(Ncfile_in.variables[ivar]):
                        self.copy_Ncaxis_with_content(
                            Ncfile_in.variables[ivar])
                    else:
                        self.copy_Ncvar(Ncfile_in.variables[ivar])


    def merge_files_from_list(self, Ncfilename_list, exclude_var=[]):
        Mf_IN = MFDataset(Ncfilename_list, "r")
        self.copy_all_dims_from_Ncfile(Mf_IN)
        self.copy_all_vars_from_Ncfile(Mf_IN, exclude_var = exclude_var)
        Mf_IN.close()


# ======================================================================
#       Wrapper for creating netCDF-like objects from Legacy GCM
#                         Fortran binary files
# ======================================================================

class Fort(object):
    """
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
    """

    class Fort_var(np.ndarray):
        """
        Sub-class that emulates a netCDF-like variable by adding the
        ``name``, ``long_name``, ``units``, and ``dimensions``
        attributes to a numpy array. Inner class for
        ``fortran_variables`` (Fort_var) that comprise the Fort file.
        Alex Kling

        A useful resource on subclassing is available at:
        https://numpy.org/devdocs/reference/arrays.classes.html

        .. note::
            Because we use an existing ``numpy.ndarray`` to define
            the object, we do not call ``__array_finalize__(self, obj)``

        :param np.ndarray: _description_
        :type np.ndarray: _type_
        :return: _description_
        :rtype: _type_
        """


        def __new__(cls, input_vals, *args, **kwargs):
            return np.asarray(input_vals).view(cls)


        def __init__(self, input_vals, name_txt, long_name_txt, units_txt,
                     dimensions_tuple):
            self.name = name_txt
            self.long_name = long_name_txt
            self.units = units_txt
            self.dimensions = dimensions_tuple
    # End inner class


    def __init__(self, filename=None, description_txt=""):
        self.filename = filename
        self.path, self.name = os.path.split(filename)
        print(f"Reading {filename}...")
        self.f = FortranFile(filename)

        if len(self.name)==12:
            # Get output number (e.g. 11 for fort.11_0070)
            self.fort_type = filename[-7:-5]
        else:
            # Case if file is simply named ``fort.11``, which is true
            # for the first file of a cold start
            self.fort_type = filename[-2:]

        self.nperday = 16 # TODO Hard-coded: 16 outputs per day
        self.nsolfile = 10 # TODO Hard-coded: 10 sols per output
        # Add time of day dimensions
        self.tod_name = tod_name = f"time_of_day_{self.nperday:02}"
        # This samples every 1.5 hours, centered at 1/2 timesteps (0.75)
        # i.e., ``np.arange(0.75, 24, 1.5)``
        self.tod = np.arange(0.5*24/self.nperday, 24, 24/self.nperday)

        # Initialize dictionaries
        self.dimensions = {}
        self.variables = {}

        if self.fort_type == "11":
            self._read_Fort11_header()
            self._read_Fort11_constants()
            self._read_Fort11_static()
            self._create_dims()
            self._read_Fort11_dynamic()
            self._add_axis_as_variables()
            # TODO monotically increasing MY: get date in MGCM date
            # format (00000)
            # -1 round to nearest 10
            self.fdate = ("%05i"%np.round(self.variables["time"][0], -1))


    # Public methods
    def write_to_fixed(self):
        """
        Create ``fixed`` file (all static variables)
        """

        Log = Ncdf(f"{self.path}/{self.fdate}.fixed.nc")

        # Define dimensions
        for ivar in ["lat", "lon", "pfull", "phalf", "zgrid"]:
            if ivar =="lon":
                cart_ax="X"
            if ivar =="lat":
                cart_ax="Y"
            if ivar in ["pfull", "phalf", "zgrid"]:
                cart_ax="Z"
            fort_var = self.variables[ivar]
            Log.add_dim_with_content(dimension_name = ivar,
                                     DATAin = fort_var,
                                     longname_txt = fort_var.long_name,
                                     units_txt = fort_var.units,
                                     cart_txt = cart_ax)
        # Log static variables
        for ivar in self.variables.keys():
            if "time" not in self.variables[ivar].dimensions:
                fort_var = self.variables[ivar]
                Log.log_variable(variable_name = ivar,
                                 DATAin = fort_var,
                                 dim_array = fort_var.dimensions,
                                 longname_txt = fort_var.long_name,
                                 units_txt = fort_var.units)
        Log.close()


    def write_to_daily(self):
        """
        Create daily file (continuous time series)
        """

        Log = Ncdf(f"{self.path}/{self.fdate}.atmos_daily.nc")

        # Define dimensions
        for ivar in ["lat", "lon", "pfull", "phalf", "zgrid"]:
            if ivar =="lon":
                cart_ax="X"
            if ivar =="lat":
                cart_ax="Y"
            if ivar in ["pfull", "phalf", "zgrid"]:
                cart_ax="Z"
            fort_var = self.variables[ivar]
            Log.add_dim_with_content(dimension_name = ivar,
                                     DATAin = fort_var,
                                     longname_txt = fort_var.long_name,
                                     units_txt = fort_var.units,
                                     cart_txt = cart_ax)

        # Add ``scalar_axis`` dimension (size 1, only used with areo)
        Log.add_dimension("scalar_axis", 1)

        # Add aggregation dimension (None size for unlimited)
        Log.add_dimension("time", None)
        fort_var = self.variables["time"]
        Log.log_axis1D(variable_name = "time", DATAin = fort_var,
                       dim_name = "time", longname_txt = fort_var.long_name,
                       units_txt = fort_var.units, cart_txt = "T")

        # Special case for the solar longitude (``areo``): needs to be
        # interpolated linearly every 16 timesteps
        ivar = "areo"
        fort_var = self.variables[ivar]
        # ``areo`` is reshaped as ``[time, scalar_axis] = [160, 1]``
        var_out = self._linInterpLs(
            np.squeeze(fort_var[:]), 16).reshape([len(fort_var), 1])
        Log.log_variable(variable_name = ivar, DATAin = var_out,
                         dim_array = fort_var.dimensions,
                         longname_txt = fort_var.long_name,
                         units_txt = fort_var.units)

        # Log dynamic variables as well as ak (pk), bk
        for ivar in self.variables.keys():
            if ("time" in self.variables[ivar].dimensions and
                ivar != "areo" or
                ivar in ["pk", "bk"]):
                fort_var = self.variables[ivar]
                Log.log_variable(variable_name = ivar, DATAin = fort_var,
                                 dim_array = fort_var.dimensions,
                                 longname_txt = fort_var.long_name,
                                 units_txt = fort_var.units)
        Log.close()


    def write_to_average(self, day_average=5):
        """
        Create average file (e.g., N-day averages [N=5 usually])
        """

        Log = Ncdf(f"{self.path}/{self.fdate}.atmos_average.nc")

        # Define dimensions
        for ivar in ["lat", "lon", "pfull", "phalf", "zgrid"]:
            if ivar == "lon":
                cart_ax = "X"
            if ivar == "lat":
                cart_ax = "Y"
            if ivar in ["pfull", "phalf", "zgrid"]:
                cart_ax = "Z"
            fort_var = self.variables[ivar]
            Log.add_dim_with_content(dimension_name = ivar,
                                     DATAin = fort_var,
                                     longname_txt = fort_var.long_name,
                                     units_txt = fort_var.units,
                                     cart_txt = cart_ax)

        # Add ``scalar_axis`` dimension (size 1, only used with areo)
        Log.add_dimension("scalar_axis", 1)

        # Add aggregation dimension (None size for unlimited)
        Log.add_dimension("time", None)

        # Perform day average and log new time axis
        time_in = self.variables["time"]
        time_out = daily_to_average(varIN = fort_var,
                                    dt_in = (time_in[1]-time_in[0]),
                                    nday = day_average,
                                    trim = True)
        Log.log_axis1D(variable_name = "time",
                       DATAin = time_out,
                       dim_name = "time",
                       longname_txt = time_in.long_name,
                       units_txt = time_in.units,
                       cart_txt = "T")

        # Log static variables
        for ivar in ["pk", "bk"]:
            fort_var = self.variables[ivar]
            Log.log_variable(variable_name = ivar,
                             DATAin = fort_var,
                             dim_array = fort_var.dimensions,
                             longname_txt = fort_var.long_name,
                             units_txt = fort_var.units)

        #Log dynamic variables
        for ivar in self.variables.keys():
            if "time" in self.variables[ivar].dimensions:
                fort_var = self.variables[ivar]
                var_out = daily_to_average(varIN = fort_var,
                                           dt_in = (time_in[1]-time_in[0]),
                                           nday = day_average,
                                           trim = True)
                Log.log_variable(variable_name = ivar,
                                 DATAin = var_out,
                                 dim_array = fort_var.dimensions,
                                 longname_txt = fort_var.long_name,
                                 units_txt = fort_var.units)
        Log.close()


    def write_to_diurn(self, day_average=5):
        """
        Create diurn file (variables organized by time of day & binned
        (typically 5-day bins)
        """

        Log = Ncdf(f"{self.path}/{self.fdate}.atmos_diurn.nc")

        # Define dimensions
        for ivar in ["lat", "lon", "pfull", "phalf", "zgrid"]:
            if ivar =="lon":
                cart_ax="X"
            if ivar =="lat":
                cart_ax="Y"
            if ivar in ["pfull" , "phalf", "zgrid"]:
                cart_ax="Z"
            fort_var=self.variables[ivar]
            Log.add_dim_with_content(dimension_name = ivar,
                                     DATAin = fort_var,
                                     longname_txt = fort_var.long_name,
                                     units_txt = fort_var.units,
                                     cart_txt = cart_ax)

        # Add scalar_axis dimension (size 1, only used with areo)
        Log.add_dimension("scalar_axis", 1)

        # Add time_of_day dimensions
        Log.add_dim_with_content(dimension_name = self.tod_name,
                                 DATAin = self.tod,
                                 longname_txt = "time of day",
                                 units_txt = "hours since 0000-00-00 00:00:00",
                                 cart_txt = "N")

        # Add aggregation dimension (None size for unlimited)
        Log.add_dimension("time", None)

        # Perform day average and log new time axis
        time_in = self.variables["time"]
        time_out = daily_to_average(varIN = time_in,
                                    dt_in = (time_in[1]-time_in[0]),
                                    nday = day_average,trim = True)
        Log.log_axis1D(variable_name = "time", DATAin = time_out,
                       dim_name = "time", longname_txt = time_in.long_name,
                       units_txt = time_in.units, cart_txt = "T")

        # Log static variables
        for ivar in ["pk", "bk"]:
            fort_var = self.variables[ivar]
            Log.log_variable(variable_name = ivar, DATAin = fort_var,
                             dim_array = fort_var.dimensions,
                             longname_txt = fort_var.long_name,
                             units_txt = fort_var.units)

        # Loop over all variables in file
        for ivar in self.variables.keys():
            if "time" in self.variables[ivar].dimensions:
                fort_var = self.variables[ivar]
                if "time" in fort_var.dimensions and ivar != "time":
                    # If time is the dimension (& not just a time array)
                    dims_in = fort_var.dimensions
                    if type(dims_in) == str:
                        # If dimension = "time" only, it is a string
                        dims_out = (dims_in,) + (self.tod_name,)
                    else:
                        # If dimensions = tuple
                        # (e.g., ``[time,lat,lon]``)
                        dims_out = ((dims_in[0],)
                                    + (self.tod_name,)
                                    + dims_in[1:])

                    var_out = daily_to_diurn(fort_var[:],
                                             time_in[0:self.nperday])
                    if day_average != 1:
                        # dt = 1 sol between two diurn timesteps
                        var_out = daily_to_average(var_out, 1., day_average)
                    Log.log_variable(ivar, var_out, dims_out,
                                     fort_var.long_name, fort_var.units)
        Log.close()


    # Public method
    def close(self):
        self.f.close()
        print(f"{self.filename} was closed")


    # Private methods
    def _read_Fort11_header(self):
        """
        Return values from ``fort.11`` header.
        ``f`` is an open ``scipy.io.FortranFile`` object

        :return: ``RUNNUM``, ``JM``, ``IM``, ``LM``, ``NL``, ``ntrace``,
            ``version``, and ``SM``

        .. note::
            In ``myhist.f``:

            write(11) RUNNUM (float), JM, IM, LAYERS, NL, NTRACE (ints),
            version (char= 7)

        These are saved as attributes (e.g., uses ``f.LAYERS`` to
        access the number of layers).
        """

        Rec = self.f.read_record("f4", "(1, 5)i4", "S7")
        self.RUNNUM = Rec[0][0]
        self.JM = Rec[1][0, 0]
        self.IM = Rec[1][0, 1]
        self.LM = Rec[1][0, 2]
        self.NL = Rec[1][0, 3]
        self.ntrace = Rec[1][0, 4]
        self.version = Rec[2][0]

        # Also compute subsurface grid (boundaries)
        self.SM = 2*self.NL + 1

    def _read_Fort11_constants(self):
        """
        Return run constants from ``fort.11`` header.
        ``f`` is an open ``scipy.io.FortranFile`` object

        .. note::
            In ``myhist.f``:

            write(11) DSIG, DXYP, GRAV, RGAS, cp, stbo, xlhtc, kapa,
            *          cmk, decmax, eccn, orbinc, vinc, sdepth, alicen,
            *          alices, egoco2n, egoco2s, npcwikg, gidn, gids

        These are saved as attributes (e.g., uses ``f.rgas`` to access
        the gas constant for the simulation.
        """

        Rec = self.f.read_record(f"(1, {self.LM})f4",
                                 f"(1, {self.JM})f4",
                                 "f4", "f4", "f4", "f4", "f4", "f4", "f4",
                                 "f4", "f4", "f4", "f4",
                                 f"(1, {self.SM})f4", "f4", "f4", "f4", "f4",
                                 "f4", "f4", "f4")
        self.dsig = np.array(Rec[0][0, :])
        self.dxyp = np.array(Rec[1][0, :])
        self.grav = Rec[2][0]
        self.rgas = Rec[3][0]
        self.cp = Rec[4][0]
        self.stbo = Rec[5][0]
        self.xlhtc = Rec[6][0]
        self.kapa = Rec[7][0]
        self.cmk = Rec[8][0]
        self.decmax = Rec[9][0]
        self.eccn = Rec[10][0]
        self.orbinc = Rec[11][0]
        self.vinc = Rec[12][0]
        self.sdepth = np.array(Rec[13][0, :])
        self.alicen = Rec[14][0]
        self.alices = Rec[15][0]
        self.egoco2n = Rec[16][0]
        self.egoco2s = Rec[17][0]


    def _read_Fort11_static(self):
        """
        Return values from ``fort.11`` header.
        ``f`` is an open ``scipy.io.FortranFile`` object

        .. note::
            In ``myhist.f``:

            write(11) TOPOG, ALSP, ZIN, NPCFLAG

        These are saved as variables (e.g., uses
        ``f.variables["zsurf"]`` to access the topography.
        """

        Rec = self.f.read_record(f"({self.IM},{self.JM})f4",
                                 f"({self.IM},{self.JM})f4",
                                 f"({self.NL},{self.IM},{self.JM})f4",
                                 f"({self.IM},{self.JM})f4")

        # Add static variables to the variables dictionary.
        self.variables["zsurf"] = self.Fort_var(-np.array(Rec[0].T/self.grav),
                                                "zsurf", "surface height",
                                                "m", ("lat", "lon"))
        self.variables["alb"] = self.Fort_var(np.array(Rec[1].T), "alb",
                                              "Surface Albedo", "mks",
                                              ("lat", "lon"))
        self.variables["thin"] = self.Fort_var(
            np.array(Rec[2].transpose([0, 2, 1])), "thin",
            "Surface Thermal Inertia", "J/m2/K/s1/2", ("zgrid", "lat", "lon"))
        self.variables["npcflag"] = self.Fort_var(np.array(Rec[3].T),
                                                  "npcflag", "Polar ice flag",
                                                  "none", ("lat", "lon"))


    def _create_dims(self):
        """
        Create dimension axis from ``IM``, ``JM`` after reading the
        header. Also compute a vertical grid structure that includes
        sigma values at the layer boundaries AND midpoints for the
        radiation code. Total size is ``2*LM+2``
        """

        JM = self.JM # JM = 36
        IM = self.IM # IM = 60
        LM = self.LM
        NL = self.NL
        self.lat = -90.0 + (180.0/JM)*np.arange(1, JM+1)
        self.lon = -180. + (360./IM)*np.arange(0, IM)

        # Compute sigma layer. Temporary arrays:
        sigK = np.zeros(2*LM + 3) # Layer midpoints + boundaries
        sigL = np.zeros(LM) # Layer midpoints only

        # These are edges and midpoints for output
        self.sigm = np.zeros(2*LM + 1)

        sigK[0:3] = 0.

        for l in range(0, LM):
            k = (2*(l) + 3)
            sigK[k] = (sigK[k-2] + self.dsig[l])
        for k in range(4, (2*LM + 3 - 1), 2):
            sigK[k] = (0.5*(sigK[k+1] + sigK[k-1]))
        for l in range(0, LM):
            sigL[l] = sigK[l*2 + 1]

        sigK[2*LM + 2] = 1.0
        self.sigm[:] = sigK[2:]

        # Subsurface layer
        # Assume this is midpoint bound so we take every other point
        # starting with the 2nd
        self.zgrid = self.sdepth[1::2] # TODO check


    def _ra_1D(self, new_array, name_txt):
        """
        ``_ra`` stands for "Return array": Append single timesteps
        along the first (``time``) dimensions
        """

        if type(new_array) != np.ndarray:
            new_array = np.array([new_array])

        # Add ``time`` axis to new data (e.g. [lat, lon]
        # -> [1, lat, lon])
        new_shape = np.append([1], new_array.shape)
        if name_txt not in self.variables.keys():
            # First time that the variable is encountered
            return new_array
        else:
            # Get values from existing array and append to it. Note
            # that ``np.concatenate((x,y))`` takes a tuple as argument.
            return np.append(self.variables[name_txt], new_array)


    def _log_var(self,name_txt, long_name, unit_txt, dimensions, Rec=None,
                 scaling=None):
        if Rec is None:
            # If no record is provided, read directly from file. Note
            # that this is reading only one timestep at the time!
            if dimensions == ("time", "lat", "lon"):
                Rec = self.f.read_reals("f4").reshape(self.JM, self.IM,
                                                      order = "F")
            if dimensions == ("time", "pfull", "lat", "lon"):
                Rec = self.f.read_reals("f4").reshape(self.JM, self.IM,
                                                      self.LM, order = "F")
            if dimensions == ("time", "zgrid", "lat", "lon"):
                Rec = self.f.read_reals("f4").reshape(self.JM, self.IM,
                                                      self.NL, order = "F")
            # If scaling, scale it!
            if scaling:
                Rec *= scaling

        # Reorganize 2D and 3D vars from
        # ``[lat, lon, lev]`` -> ``[lev, lat, lon]``
        if (dimensions == ("time", "pfull", "lat", "lon") or
            dimensions == ("time", "zgrid", "lat", "lon")):
            Rec = Rec.transpose([2, 0, 1])

        # Set to pole point to value at N-1
        Rec[..., -1, :] = Rec[..., -2, :]

        # Add time axis to new data (e.g. [lat, lon] -> [1, lat, lon])
        new_shape = np.append([1], Rec.shape)
        if name_txt not in self.variables.keys():
            # First time that the variable is encountered
            Rec = Rec.reshape(new_shape)
        else:
            # Get values from existing array and append to it. Note
            # that ``np.concatenate((x, y))`` takes a tuple as argument.
            Rec = np.concatenate((self.variables[name_txt],
                                  Rec.reshape(new_shape)))
        # Log the variable
        self.variables[name_txt] = self.Fort_var(Rec, name_txt, long_name,
                                                 unit_txt, dimensions)


    def _read_Fort11_dynamic(self):
        """
        Read variables from ``fort.11`` files that change with each
        timestep.

        In ``mhistv.f``::

            WRITE(11) TAU, VPOUT, RSDIST, TOFDAY, PSF, PTROP, TAUTOT,
                *          RPTAU, SIND, GASP
                WRITE(11) NC3, NCYCLE

                WRITE(11) P
                WRITE(11) T
                WRITE(11) U
                WRITE(11) V
                WRITE(11) GT
                WRITE(11) CO2ICE
                WRITE(11) STRESSX
                WRITE(11) STRESSY
                WRITE(11) TSTRAT
                WRITE(11) TAUSURF
                WRITE(11) SSUN
                WRITE(11) QTRACE
                WRITE(11) QCOND
                write(11) STEMP
                write(11) fuptopv, fdntopv, fupsurfv, fdnsurfv
                write(11) fuptopir, fupsurfir, fdnsurfir
                write(11) surfalb
                write(11) dheat
                write(11) geot
        """

        # Typically ``nsteps = 16 x 10 = 160``
        nsteps = self.nperday * self.nsolfile
        append = False
        for iwsol in range(0, nsteps):
            Rec = self.f.read_record("f4")
            # TAU = Rec[0]
            # VPOUT = Rec[1]
            # RSDIST = Rec[2]
            # TOFDAY = Rec[3]
            # PSF = Rec[4]
            # PTROP = Rec[5]
            # TAUTOT = Rec[6]
            # RPTAU = Rec[7]
            # SIND = Rec[8]
            # GASP2 = Rec[9]

            self.variables["time"] = self.Fort_var(
                self._ra_1D(Rec[0]/24, "time"), "time",
                "elapsed time from the start of the run",
                "days since 0000-00-00 00:00:00", ("time"))

            self.variables["areo"] = self.Fort_var(
                self._ra_1D(Rec[1].reshape([1, 1]), "areo"), "areo",
                "solar longitude", "degree", ("time", "scalar_axis"))
            # TODO "areo" monotically increasing?

            self.variables["rdist"] = self.Fort_var(
                self._ra_1D(Rec[2], "rdist"), "rdist",
                "square of the Sun-Mars distance", "(AU)**2", ("time"))

            self.variables["tofday"]=self.Fort_var(
                self._ra_1D(Rec[3], "tofday"), "npcflag", "time of day",
                "hours since 0000-00-00 00:00:00", ("time"))
            # TODO "tofday" edge or center?

            self.variables["psf"] =  self.Fort_var(
                self._ra_1D(Rec[4]*100, "psf"), "psf",
                "Initial global surface pressure", "Pa", ("time"))

            self.variables["ptrop"] = self.Fort_var(
                self._ra_1D(Rec[5], "ptrop"), "ptrop",
                "pressure at the tropopause", "Pa", ("time"))

            self.variables["tautot"]=self.Fort_var(
                self._ra_1D(Rec[6], "tautot"), "tautot",
                "Input (global) dust optical depth at the reference pressure",
                "none", ("time"))

            self.variables["rptau"] = self.Fort_var(
                self._ra_1D(Rec[7]*100, "rptau"), "rptau",
                "reference pressure for dust optical depth", "Pa", ("time"))

            self.variables["sind"] = self.Fort_var(
                self._ra_1D(Rec[8], "sind"), "sind",
                "sine of the sub-solar latitude", "none", ("time"))

            self.variables["gasp"] = self.Fort_var(
                self._ra_1D(Rec[9]*100, "gasp"), "gasp",
                "global average surface pressure", "Pa", ("time"))

            Rec = self.f.read_record("i4")
            # NC3 = Rec[0]
            # NCYCLE = Rec[1]

            self.variables["nc3"] = self.Fort_var(
                self._ra_1D(Rec[0], "nc3"), "nc3",
                "full COMP3 is done every nc3 time steps.", "None", ("time"))

            self.variables["ncycle"] = self.Fort_var(
                self._ra_1D(Rec[1], "ncycle"), "ncycle", "ncycle", "none",
                ("time"))

            self._log_var("ps", "surface pressure", "Pa",
                          ("time", "lat", "lon"), scaling = 100)

            self._log_var("temp", "temperature", "K",
                          ("time", "pfull", "lat", "lon"))

            self._log_var("ucomp", "zonal wind", "m/sec",
                          ("time", "pfull", "lat", "lon"))

            self._log_var("vcomp", "meridional wind", "m/s",
                          ("time", "pfull", "lat", "lon"))

            self._log_var("ts", "surface temperature", "K",
                          ("time", "lat", "lon"))

            self._log_var("snow", "surface amount of CO2 ice on the ground",
                          "kg/m2", ("time", "lat", "lon"))

            self._log_var("stressx", "zonal component of surface stress",
                          "kg/m2", ("time", "lat", "lon"))

            self._log_var("stressy", "merdional component of surface stress",
                          "kg/m2", ("time", "lat", "lon"))

            self._log_var("tstrat", "stratosphere temperature", "K",
                          ("time", "lat", "lon"))

            self._log_var("tausurf",
                          "visible dust optical depth at the surface.", "none",
                          ("time", "lat", "lon"))

            self._log_var("ssun", "solar energy absorbed by the atmosphere",
                          "W/m2", ("time", "lat", "lon"))

            # Write(11) QTRACE # dust mass:1, dust number 2||
            # water ice mass: 3 and water ice number 4||
            # dust core mass:5|| water vapor mass: 6
            Rec = self.f.read_reals("f4").reshape(self.JM, self.IM, self.LM,
                                                  self.ntrace, order = "F")

            self._log_var("dst_mass", "dust aerosol mass mixing ratio",
                          "kg/kg", ("time", "pfull", "lat", "lon"),
                          Rec = Rec[..., 0])

            self._log_var("dst_num", "dust aerosol number", "number/kg",
                          ("time", "pfull", "lat", "lon"), Rec = Rec[..., 1])

            self._log_var("ice_mass", "water ice aerosol mass mixing ratio",
                          "kg/kg", ("time", "pfull", "lat", "lon"),
                          Rec = Rec[..., 2])

            self._log_var("ice_num", "water ice  aerosol number", "number/kg",
                          ("time", "pfull", "lat", "lon"), Rec = Rec[..., 3])

            self._log_var("cor_mass",
                          "dust core mass mixing ratio for water ice", "kg/kg",
                          ("time", "pfull", "lat", "lon"), Rec = Rec[..., 4])

            self._log_var("vap_mass", "water vapor mass mixing ratio", "kg/kg",
                          ("time", "pfull", "lat", "lon"), Rec = Rec[..., 5])

            # write(11) QCOND   dust mass:1, dust number 2||
            # water ice mass: 3 and water ice number 4||
            # dust core mass:5|| water vapor mass: 6
            Rec = self.f.read_reals("f4").reshape(self.JM, self.IM,
                                                  self.ntrace, order = "F")

            self._log_var("dst_mass_sfc", "dust aerosol mass on the surface",
                          "kg/m2", ("time", "lat", "lon"), Rec = Rec[..., 0])

            self._log_var("dst_num_sfc", "dust aerosol number on the surface",
                          "number/m2", ("time", "lat", "lon"),
                          Rec = Rec[..., 1])

            self._log_var("ice_mass_sfc",
                          "water ice aerosol mass on the surface", "kg/m2",
                          ("time", "lat", "lon"), Rec = Rec[..., 2])

            self._log_var("ice_num_sfc",
                          "water ice  aerosol number on the surface",
                          "number/m2", ("time", "lat", "lon"),
                          Rec = Rec[..., 3])

            self._log_var("cor_mass_sfc",
                          "dust core mass for water ice on the surface",
                          "kg/m2", ("time", "lat", "lon"), Rec = Rec[..., 4])

            self._log_var("vap_mass_sfc", "water vapor mass on the surface",
                          "kg/m2", ("time", "lat", "lon"), Rec = Rec[..., 5])

            # write(11) stemp
            Rec = self.f.read_reals("f4").reshape(self.JM, self.IM, self.NL,
                                                  order = "F")

            self._log_var("soil_temp", "sub-surface soil temperature", "K",
                          ("time", "zgrid", "lat", "lon"), Rec = Rec)

            # write(11) fuptopv,  fdntopv,  fupsurfv,  fdnsurfv
            #.. NOTE: the following are read in Fortran order:
            # (IM, JM) > (60, 36) and not (JM, IM) > (36, 60) since we
            # are not using the ``order = "F"`` flag. These need to be
            # transposed.
            Rec = self.f.read_record(f"({self.IM}, {self.JM})f4",
                                     f"({self.IM}, {self.JM})f4",
                                     f"({self.IM}, {self.JM})f4",
                                     f"({self.IM}, {self.JM})f4")

            self._log_var("fuptopv",
                          "upward visible flux at the top of the atmosphere",
                          "W/m2", ("time", "lat", "lon"), Rec = Rec[0].T)

            self._log_var("fdntopv",
                          "downward visible flux at the top of the atmosphere",
                          "W/m2", ("time", "lat", "lon"), Rec = Rec[1].T)

            self._log_var("fupsurfv", "upward visible flux at the surface",
                          "W/m2", ("time", "lat", "lon"), Rec = Rec[2].T)

            self._log_var("fdnsurfv", "downward visible flux at the surface",
                          "W/m2", ("time", "lat", "lon"), Rec = Rec[3].T)

            # write(11) fuptopir, fupsurfir, fdnsurfir

            # the following are read in fortran order:
            # (IM, JM) > (60, 36) and not (JM, IM) > (36, 60) since we
            # are not using the ``order = "F"`` flag. These need to be
            # transposed.
            Rec = self.f.read_record(f"({self.IM}, {self.JM})f4",
                                     f"({self.IM}, {self.JM})f4",
                                     f"({self.IM}, {self.JM})f4")

            self._log_var("fuptopir",
                          "upward IR flux at the top of the atmosphere",
                          "W/m2", ("time", "lat", "lon"), Rec = Rec[0].T)

            self._log_var("fupsurfir",
                          "upward IR flux at the surface", "W/m2",
                          ("time", "lat", "lon"), Rec = Rec[1].T)

            self._log_var("fdnsurfir",
                          "downward IR flux at the surface", "W/m2",
                          ("time", "lat", "lon"), Rec = Rec[2].T)

            # write(11) surfalb
            self._log_var("surfalb",
                          ("surface albedo in the visible, soil or H2O, CO2 "
                          "ices if present"), "none", ("time", "lat", "lon"))

            # write(11) dheat
            # write(11) geot
            self._log_var("dheat", "diabatic heating rate", "K/sol",
                          ("time", "pfull", "lat", "lon"))

            self._log_var("geot", "geopotential", "m2/s2",
                          ("time", "pfull", "lat", "lon"))


    def _add_axis_as_variables(self):
        """
        Add dimensions to the file as variables
        """

        self.variables["lat"] = self.Fort_var(self.lat, "lat", "latitude",
                                              "degrees_N", ("lat"))
        self.variables["lon"] = self.Fort_var(self.lon, "lon", "longitude",
                                              "degrees_E", ("lon"))
        sgm = self.sigm
        pref = self.variables['psf'][0] # 7.01*100 Pa
        pk = np.zeros(self.LM + 1)
        bk = np.zeros(self.LM + 1)
        pfull = np.zeros(self.LM)
        phalf = np.zeros(self.LM + 1)

        pk[0] = 0.08/2
        # TODO [AK] change pk[0]=.08 to pk[0]=.08/2, otherwise
        # phalf[0] > phalf[1]

        for iz in range(self.LM):
            bk[iz+1] = sgm[2*iz + 2]

        # Output in Pa
        phalf[:] = (pk[:] + pref*bk[:])

        # First layer
        if pk[0] == 0 and bk[0] == 0:
            pfull[0] = 0.5*(phalf[0] + phalf[1])
        else:
            pfull[0] = ((phalf[1] - phalf[0])
                        /(np.log(phalf[1]) - np.log(phalf[0])))
        # Rest of layers:
        pfull[1:] = ((phalf[2:] - phalf[1:-1])
                     /(np.log(phalf[2:]) - np.log(phalf[1:-1])))

        self.variables["phalf"] = self.Fort_var(
            phalf, "phalf","ref half pressure level", "Pa", ("phalf"))
        self.variables["pfull"] = self.Fort_var(
            pfull, "pfull", "ref full pressure level", "Pa", ("pfull"))
        self.variables["bk"] = self.Fort_var(
            bk, "bk", "vertical coordinate sigma value", "none", ("phalf"))
        self.variables["pk"] = self.Fort_var(
            pk, "pk", "pressure part of the hybrid coordinate", "Pa",
            ("phalf"))
        self.variables["zgrid"] = self.Fort_var(
            self.zgrid, "zgrid", "depth at the mid-point of each soil layer",
            "m", ("zgrid"))


    def _linInterpLs(self, Ls, stride=16):
        """
        Linearly interpolate a step-wise 1D array

        :param Ls: input solar longitude
        :type Ls: float
        :param stride: default stride
        :type stride: int
        :return Ls_out: solar longitude (Ls) [float]

        ..note::
            In the Legacy GCM fortran binaries, the solar
            longitude is only updated once per day, implying that 16
            successive timesteps would have the same ls value. This routine
            linearly interpolates the Ls between those successive values.
        """

        Ls = np.array(Ls)
        Ls_out = np.zeros_like(Ls)
        Lsdi = Ls[::stride]
        # Add an end point using the last Delta Ls:
        Lsdi = np.append(Lsdi, (2*Lsdi[-1] - Lsdi[-2]))

        for i in range(len(Ls)//stride):
            Ls_out[(i*stride):((i+1)*stride)] = (np.arange(0, stride)
                                                 / np.float32(stride)
                                                 * (Lsdi[i+1] - Lsdi[i])
                                                 + Lsdi[i])
        return Ls_out
