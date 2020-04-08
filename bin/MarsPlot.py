#!/usr/bin/env python3

#Load generic Python Modules
import argparse #parse arguments
import os       #access operating systems function
import subprocess #run command
import sys       #system command


#TODO remove this block to use package instead
#==============

#sys.path.append('/Users/akling/amesgcm/amesgcm/')
#from Script_utils import          check_file_tape,prYellow,prRed,prCyan,prGreen,prPurple, print_fileContent
from amesgcm.Script_utils import check_file_tape,prYellow,prRed,prCyan,prGreen,prPurple, print_fileContent
from amesgcm.FV3_utils import lon360_to_180,lon180_to_360
#=====Attempt to import specific scientic modules one may not find in the default python on NAS ====
try:
    import matplotlib
    matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.cm import get_cmap
    from matplotlib.ticker import MultipleLocator, FuncFormatter  #format ticks
    from netCDF4 import Dataset, MFDataset
    from numpy import sqrt, exp, max, mean, min, log, log10

except ImportError as error_msg:
    prYellow("Error while importing modules")
    prYellow('Your are using python '+str(sys.version_info[0:3]))
    prYellow('Please, source your virtual environment');prCyan('    source envPython3.7/bin/activate.csh \n')
    print("Error was: "+ error_msg.message)
    exit()
except Exception as exception:
    # Output unexpected Exceptions.
    print(exception, False)
    print(exception.__class__.__name__ + ": " + exception.message)
    exit()


#======================================================
#                  ARGUMENTS PARSER
#======================================================

global current_version;current_version=1.4
parser = argparse.ArgumentParser(description="""\033[93mAnalysis Toolkit for the Ames GCM, V%s\033[00m """%(current_version),
                                formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('custom_file', nargs='?',type=argparse.FileType('r'),default=None, #sys.stdin
                             help='Use optional input file Custom.in to plot the graphs \n'
                                  '> Usage: MarsPlot Custom.in  [other options]')

parser.add_argument('-i', '--inspect_file', default=None,
                 help='Inspect Netcdf file content. Variables are sorted by dimensions \n'
                      '> Usage: MarsPlot -i 00000.atmos_davg.nc')

parser.add_argument('-d','--date', nargs='+',default=None,
                 help='Specify the range of files to use, default is the last file  \n'
                      '> Usage: MarsPlot Custom.in -d 700     (one file) \n'
                      '         MarsPlot Custom.in -d 350 700 (start file end file)')

parser.add_argument('--template','-template', action='store_true',
                        help="""Generate a template Custom.in for customization of the plots.\n """
                             """(Use '--temp' for a skinned version of the template)\n""")
parser.add_argument('-temp','--temp', action='store_true',help=argparse.SUPPRESS) #same as --template but without the instructions

parser.add_argument('-do','--do', nargs=1,type=str,default=None, #sys.stdin
                             help='(Re)-use a template file my_custom.in. First search in ~/env/mars_templates/,\n'
                                 '                                                then in /u/mkahre/MCMC/analysis/working/templates/ \n'
                                  '> Usage: MarsPlot -do my_custom [other options]')

parser.add_argument("-o", "--output",default="pdf",
                 choices=['pdf','eps','png'],
                 help='Output file format\n'
                       'Default is pdf if ghostscript (gs) is available and png otherwise\n'
                        '> Usage: MarsPlot Custom.in -o png \n'
                        '       : MarsPlot Custom.in -o png -pw 500 (set pixel width to 500, default is 2000)\n')

parser.add_argument("-pw", "--pwidth",default=2000,type=float,
                 help=argparse.SUPPRESS)


parser.add_argument('-dir', '--directory', default=os.getcwd(),
                 help='Target directory if input files are not present in current directory \n'
                      '> Usage: MarsPlot Custom.in [other options] -dir /u/akling/FV3/verona/c192L28_dliftA/history')

parser.add_argument('--debug',  action='store_true', help='Debug flag: do not by-pass errors on a particular figure')
#======================================================
#                  MAIN PROGRAM
#======================================================
def main():

    global output_path ; output_path = os.getcwd()
    global input_paths  ; input_paths=[];input_paths.append(parser.parse_args().directory)
    global out_format  ; out_format=parser.parse_args().output
    global debug       ;debug =parser.parse_args().debug
    global Ncdf_num         #host the simulation timestamps
    global objectList      #contains all figure object
    global customFileIN    #template name
    global simulation_name #ref simulation
    global levels;levels=21 #number of contour for 2D plots
    global label_size;label_size=12 #Label size for title, xalbel, ylabel
    global my_dpi;my_dpi=96.        #pixel per inch for figure output
    global pixel_width;pixel_width=parser.parse_args().pwidth #pixel width for saving figure



    objectList=[Fig_2D_lon_lat('fixed.zsurf',True),\
                Fig_2D_lat_press('atmos_average.ucomp',True),\
                Fig_2D_time_lat('atmos_average.taudust_IR',False),\
                Fig_2D_lon_press('atmos_average_pstd.temp',False),\
                Fig_2D_time_press('atmos_average_pstd.temp',False),\
                Fig_2D_lon_time('atmos_average.temp',False),\
                Fig_1D('atmos_average.temp',False)]
        #=============================
    #----------Group together the 1st two figures----
    objectList[0].subID=1;objectList[0].nPan=2 #1st object of a 2 pannel figure
    objectList[1].subID=2;objectList[1].nPan=2 #2nd object of a 2 pannel figure

    # Begin main loop:


    # ----- Option 1 :Inspect content of a Netcdf file ----
    if parser.parse_args().inspect_file:
        check_file_tape(parser.parse_args().inspect_file,abort=False) #NAS-specific, check if the file is on tape
        print_fileContent(parser.parse_args().inspect_file)

        # ----- Option 2: Generate a template file ----
    elif parser.parse_args().template or parser.parse_args().temp:
        make_template()

    # --- Gather simulation information from template or inline argument
    else:

        # --- Option 2, case A:   Use Custom.in  for everything ----
        if parser.parse_args().custom_file:
           print('Reading '+parser.parse_args().custom_file.name)
           namelist_parser(parser.parse_args().custom_file.name)


        # --- Option 2, case B:   Use Custom.in in ~/FV3/templates for everything ----
        if parser.parse_args().do:
           print('Reading '+path_to_template(parser.parse_args().do))
           namelist_parser(path_to_template(parser.parse_args().do))


        # set bounds  (e.g. starting file, end file)
        if parser.parse_args().date: #a date single date or a range is provided
            # first check if the value provided is the right type
            try:
                bound=np.asarray(parser.parse_args().date).astype(np.float)
            except Exception as e:
                prRed('*** Syntax Error***')
                prRed("""Please use:   'MarsPlot Custom.in -d XXXX [YYYY] -o out' """)
                exit()

        else: # no date is provided, default is last file only
            bound=np.asarray(get_Ncdf_num()[-1])

        #-----

        Ncdf_num=get_Ncdf_num() #Get all timestamps in directory
        Ncdf_num=select_range(Ncdf_num,bound)  # Apply bounds to the desired dates
        nfiles=len(Ncdf_num)                   #number of timestamps


        #print('MarsPlot is running...')
        simulation_name=input_paths[0].split('/')[-2]
        #Make a ./plots folder in the current directory if it does not exist
        dir_plot_present=os.path.exists(output_path+'/'+'plots')
        if not dir_plot_present:
            os.makedirs(output_path+'/'+'plots')

        fig_list=list()#list of figures

        #============Do plots==================
        global i_list;
        for i_list in range(0,len(objectList)):

            status=objectList[i_list].plot_type+' :'+objectList[i_list].varfull
            progress(i_list,len(objectList),status,None) #display the figure in progress

            objectList[i_list].do_plot()

            if objectList[i_list].success and out_format=='pdf' and not debug : sys.stdout.write("\033[F");sys.stdout.write("\033[K") #if success,flush the previous output

            status=objectList[i_list].plot_type+' :'+objectList[i_list].varfull+objectList[i_list].fdim_txt
            progress(i_list,len(objectList),status,objectList[i_list].success)
            # Add the figure to the list of figures
            if objectList[i_list].subID==objectList[i_list].nPan: #only for the last pannel of a subplot
                if i_list< len(objectList)-1 and not objectList[i_list+1].addLine:
                    fig_list.append(objectList[i_list].fig_name)
                #Last subplot
                if i_list== len(objectList)-1 :fig_list.append(objectList[i_list].fig_name)

        progress(100,100,'Done')# 100% completed


        #========Making multipage pdf=============
        if out_format=="pdf" and len(fig_list)>0:
            print('Merging figures...')
            #print("Plotting figures:",fig_list)
            debug_filename=output_path+'/.debug_MCMC_plots.txt' #debug file (masked), use to redirect the outputs from ghost script
            fdump = open(debug_filename, 'w') #
            #Construct list of figures----

            all_fig=' '
            for figID in fig_list:
                all_fig+=figID+' '

            #Output name for the pdf
            try:
                input_file=output_path+'/'+parser.parse_args().custom_file.name
                basename=input_file.split('/')[-1].split('.')[0].strip() #get the input file name, e.g "Custom_01" or
            except: #Special case where no Custom.in is provided
                basename='Custom'


            #default name is Custom.in, output Diagnostics.pdf
            if basename=='Custom':
                output_pdf=fig_name=output_path+'/'+'Diagnostics.pdf'
            #default name is Custom_XX.in, output Diagnostics_XX.pdf
            elif  basename[0:7]=="Custom_":
                output_pdf=fig_name=output_path+'/Diagnostics_'+basename[7:9]+'.pdf' #same name as input file
            #name is different use it
            else:
                output_pdf=fig_name=output_path+'/'+basename+'.pdf' #same name as input file

            #command to make a multipage pdf out of the the individual figures using ghost scritp.
            # Also remove the temporary files when done
            cmd_txt='gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dEPSCrop -sOutputFile='+output_pdf+' '+all_fig
            try:
                #Test the ghost scrit and remove command, exit otherwise--
                subprocess.check_call(cmd_txt,shell=True, stdout=fdump, stderr=fdump)
                #Execute the commands now
                subprocess.call(cmd_txt,shell=True, stdout=fdump, stderr=fdump) #run ghostscript to merge the pdf
                subprocess.call('rm -f '+all_fig,shell=True, stdout=fdump, stderr=fdump)#remove temporary pdf figs
                subprocess.call('rm -f '+debug_filename,shell=True)#remove debug file

                #If the plot directory was not present initially, remove it
                if not dir_plot_present:
                    subprocess.call('rm -r '+output_path+'/plots',shell=True)
                give_permission(output_pdf)
                print(output_pdf + ' was generated')

            except subprocess.CalledProcessError:
                print("ERROR with ghostscript when merging pdf, please try alternative formats")
                if debug:raise



#======================================================
#                  DATA OPERATION UTILITIES
#======================================================

def shift_data(lon,data):
    '''
    This function shift the longitude and data from a 0->360 to a -180/+180 grid.
    Args:
        lon: 1D array of longitude 0->360
        data: 2D array with last dimension being the longitude
    Returns:
        lon: 1D array of longitude -180/+180
        data: shifted data
    Note: Use np.ma.hstack instead of np.hstack to keep the masked array properties
    '''
    lon_180=lon.copy()
    nlon=len(lon_180)
    # for 1D plots: if 1D, reshape array
    if len(data.shape) <=1:
        data=data.reshape(1,nlon)
    #===
    lon_180[lon_180>180]-=360.
    data=np.hstack((data[:,lon_180<0],data[:,lon_180>=0]))
    lon_180=np.append(lon_180[lon_180<0],lon_180[lon_180>=0])
    # If 1D plot, squeeze array
    if data.shape[0]==1:
        data=np.squeeze(data)
    return lon_180,data


def MY_func(Ls_cont):
    '''
    This function return the Mars Year
    Args:
        Ls_cont: solar longitude, contineuous
    Returns:
        MY : int the Mars year
    '''
    return (Ls_cont)//(360.)+1


def get_topo_2D(simuID,sol_array):

    '''
    This function shift the topography from a file
    Returns:
        zsurf: the topography
    '''
    global input_paths
    global Ncdf_num
    if sol_array:
        Sol_num_current=sol_array
    else : #no sol requested, use default as provided by MarsPlot Custom.in -d sol
        Sol_num_current =Ncdf_num
    file_list = input_paths[simuID]+'/%05d.'%( Sol_num_current[0])+'fixed.nc' #TODO file list multiple simulations
    f=Dataset(file_list, 'r', format='NETCDF4_classic')
    zsurf=f.variables['zsurf'][:]
    f.close()
    return zsurf


def get_lon_index(lon_query_180,lons):
    '''
    Given a range of requested longitudes, return the indexes to extract data from the netcdf file
    Args:
        lon_query_180: requested longitudes in -180/+180 units: value, [min, max] or None
        lons:          1D array of longitude in the native coordinates (0->360)
    Returns:
        loni: 1D array of file indexes
        txt_lon: text descriptor for the extracted longitudes
    *** Note that the keyword 'all' is passed as -99999 by the rT() functions
    '''
    Nlon=len(lons)
    lon_query_180=np.array(lon_query_180)

    #If None, set to default, i.e 'all' for a zonal average
    if lon_query_180.any()==None: lon_query_180=np.array(-99999)
    #one longitude is provided
    if lon_query_180.size==1:
        #request zonal average
        if lon_query_180==-99999:
            loni=np.arange(0,Nlon)
            txt_lon=', zonal avg'
        else:
            #get closet value
            lon_query_360=lon180_to_360(lon_query_180)
            loni=np.argmin(np.abs(lon_query_360-lons))
            txt_lon=', lon=%.1f'%(lon360_to_180(lons[loni]))
    # a range is requested
    elif lon_query_180.size==2:
        lon_query_360=lon180_to_360(lon_query_180)
        loni_bounds=np.array([np.argmin(np.abs(lon_query_360[0]-lons)),np.argmin(np.abs(lon_query_360[1]-lons))])
        #if loni_bounds[0]>loni_bounds[1]:loni_bounds=np.flipud(loni_bounds) #lon should be increasing for extraction #TODO
        loni=np.arange(loni_bounds[0],loni_bounds[1])
        lon_bounds_180=lon360_to_180([lons[loni[0]],lons[loni[-1]]])
        if lon_bounds_180[0]>lon_bounds_180[1]:lon_bounds_180=np.flipud(lon_bounds_180) #lon should be also increasing for display
        txt_lon=', lon=avg[%.1f<->%.1f]'%(lon_bounds_180[0],lon_bounds_180[1])
    return loni,txt_lon

def get_lat_index(lat_query,lats):
    '''
    Given a range of requested latitudes, return the indexes to extract data from the netcdf file
    Args:
        lat_query: requested latitudes -90/+90
        lats:      1D array of latitudes in the native coordinates
    Returns:
        lati: 1D array of file indexes
        txt_lat: text descriptor for the extracted latitudes
    *** Note that the keyword 'all' is passed as -99999 by the rT() functions
    '''
    Nlat=len(lats)
    lat_query=np.array(lat_query)
    #If None, set to default, i.e equator
    if lat_query.any()==None: lat_query=np.array(0.)
    #one latitude is provided
    if lat_query.size==1:
        #request meridional average
        if lat_query==-99999:
            lati=np.arange(0,Nlat)
            txt_lat=', merid. avg'
        else:
            #get closet value
            lati=np.argmin(np.abs(lat_query-lats))
            txt_lat=', lat=%g'%(lats[lati])
    # a range is requested
    elif lat_query.size==2:
        lat_bounds=np.array([np.argmin(np.abs(lat_query[0]-lats)),np.argmin(np.abs(lat_query[1]-lats))])
        if lat_bounds[0]>lat_bounds[1]:lat_bounds=np.flipud(lat_bounds) #lat should be incresing for extraction
        lati=np.arange(lat_bounds[0],lat_bounds[1]+1)
        txt_lat=', lat=avg[%g<->%g]'%(lats[lati[0]],lats[lati[-1]])
    return lati,txt_lat


def get_level_index(level_query,levs):
    '''
    Given a range of requested pressures (resp. depth for 'zgrid'), return the indexes to extract data from the netcdf file
    Args:
        level_query: requested  pressure in [Pa] (resp. depth in [m])
        levs:         1D array of levels in the native coordinates [Pa] (resp. [m])
    Returns:
        levi: 1D array of file indexes
        txt_lev: text descriptor for the extracted pressure (resp. depth)
    *** Note that the keyword 'all' is passed as -99999 by the rT() functions
    '''
    level_query=np.array(level_query)
    Nz=len(levs)
    #If None, set to default, i.e  surface
    if level_query.any()== None: level_query=np.array(2*10**7) #a big number > Psfc (even for a 10bar Early Mars simulation)

    #one level is provided
    if level_query.size==1:
        #average
        if level_query==-99999:
            levi=np.arange(0,Nz)
            txt_level=', column avg'
        #specific level
        else:
            levi=np.argmin(np.abs(level_query-levs))

        # provide smart labelling
            if level_query>10.**7: #None, i.e sfc was requested
                txt_level=', at sfc'
            else:
                #txt_level=', lev=%g Pa'%(levs[levi])
                txt_level=', lev={0:1.2e} Pa'.format(levs[levi])

    elif level_query.size==2: #bounds are provided
        levi_bounds=np.array([np.argmin(np.abs(level_query[0]-levs)),np.argmin(np.abs(level_query[1]-levs))])
        if levi_bounds[0]>levi_bounds[1]:levi_bounds=np.flipud(levi_bounds) #level should be incresing for extraction
        levi=np.arange(levi_bounds[0],levi_bounds[1]+1)
        lev_bounds=[levs[levi[0]],levs[levi[-1]]] #this is for display
        if lev_bounds[0]<lev_bounds[1]:lev_bounds=np.flipud(lev_bounds) #lev should be also decreasing for display
        txt_level=', lev=avg[{0:1.2e}<->{1:1.2e}] Pa'.format(lev_bounds[0],lev_bounds[1])

    return levi,txt_level

def get_time_index(Ls_query_360,Ls,t_day):
    '''
    Given a range of requested solar longitude [0-360], return the indexes to extract data from the netcdf file.
    First try the Mars year of the last timestep, then try the year before then picks whichever Ls period is closest to the requested date.

    Args:
        Ls_query_360: requested  solar longitudes
        Ls_c:         1D array of continueous solar longitudes
        t_day:        time array in days
    Returns:
        ti: 1D array of file indexes
        txt_time: text descriptor for the extracted solar longitudes
    *** Note that the keyword 'all' is passed as -99999 by the rT() functions
    '''
    #TODO t_day is not used at the moment
    Nt=len(Ls)
    Ls_query_360=np.array(Ls_query_360)


    #If None, set to default, i.e last time step
    if Ls_query_360.any()==None: Ls_query_360=np.mod(Ls[-1],360.)

    #one time is provided
    if Ls_query_360.size==1:
        #time average average requested
        if Ls_query_360==-99999:
            ti=np.arange(0,Nt)
            txt_time=', time avg'
        else:
            #get the Mars year of the last timestep in the file
            MY_end=MY_func(Ls[-1]) #number of Mars year at the end of the file.
            if MY_end >=1:
            #check if the desired Ls is available for this Mars Year
                Ls_query=Ls_query_360+(MY_end-1)*360. #(MY starts at 1, not zero)
            else: 
                Ls_query=Ls_query_360    
            #If this time is greater that the last Ls, look one year back
            if Ls_query>Ls[-1] and MY_end>1:
                MY_end-=1 #one year back
                Ls_query=Ls_query_360+(MY_end-1)*360.
            ti=np.argmin(np.abs(Ls_query-Ls))
            txt_time=', Ls= (MY%2i) %.2f'%(MY_end,np.mod(Ls[ti],360.))

    # a range is requested
    elif Ls_query_360.size==2:

        #get the Mars year of the last timestep in the file
        MY_last=MY_func(Ls[-1]) #number of Mars year at the end of the file.
        if MY_last >=1:
        #try the mars year of the last time step
            Ls_query_last=Ls_query_360[1]+(MY_last-1)*360.
        else:
            Ls_query_last=Ls_query_360[1]
        #First consider the further end of the desired range
        #This time is greater that the last Ls, look one year back
        if Ls_query_last>Ls[-1] and  MY_last>1:
            MY_last-=1
            Ls_query_last=Ls_query_360[1]+(MY_last-1)*360. #(MY starts at 1, not zero)
        ti_last=np.argmin(np.abs(Ls_query_last-Ls))
        #then get the first value, for that Mars year
        MY_beg=MY_last.copy()
        #try the mars year of the last time step
        Ls_query_beg=Ls_query_360[0]+(MY_beg-1)*360.
        ti_beg=np.argmin(np.abs(Ls_query_beg-Ls))

        #if the begining value is higher, search in the year before for ti_beg
        if ti_beg>=ti_last:
            MY_beg-=1
            Ls_query_beg=Ls_query_360[0]+(MY_beg-1)*360.
            ti_beg=np.argmin(np.abs(Ls_query_beg-Ls))
            
   
        ti=np.arange(ti_beg,ti_last+1)
        
        Ls_bounds=[Ls[ti[0]],Ls[ti[-1]]] #this is for display
        txt_time=', Ls= avg [(MY%2i) %.2f <-> (MY%2i) %.2f]'%(MY_beg,np.mod(Ls_bounds[0],360.),MY_last,np.mod(Ls_bounds[1],360.))


    return ti,txt_time

#======================================================
#                  TEMPLATE UTILITIES
#======================================================

def filter_input(txt,typeIn='char'):
    '''
    Read Template for the type of data expected
    Args:
        txt: a string, typical the right-hand sign of an equal sign '3', '3,4', or 'all'
        typeIn: type of data expected: 'char', 'float', 'int', 'bool'
    Returns:
        out: float or 1D array [val1,val2] in the expected format

    '''
    if txt =='None' or not txt: #None or empty string
        return None

    if "," in txt: #two values are provided
        answ = []
        for i in range(0,len(txt.split(','))):
            #== For a 'char', read all text as one
            #if typeIn=='char': answ.append(txt.split(',')[i].strip())
            if typeIn=='char': answ= txt
              #====
            if typeIn=='float':answ.append(np.float(txt.split(',')[i].strip()))
            if typeIn=='int':  answ.append(np.int(txt.split(',')[i].strip()))
            if typeIn=='bool': answ.append(txt.split(',')[i].strip()=='True')
        return answ
    else:
        if typeIn=='char':
            answ= txt
        if typeIn=='bool':
            answ=  ('True'==txt)
        #for float and int type, pass the 'all' key word as -99999
        if typeIn=='float':
            if txt=='all':
                answ= -99999.
            elif txt=='AXIS':
                answ= -88888.
            else:
                answ= np.float(txt)
        if typeIn=='int':
            if txt=='all':
                answ= -99999
            else:
                answ=  np.int(txt)
  #would be True is text matches
        return answ

def rT(typeIn='char'):
    '''
    Read Template for the type of data expected
    Args:
        typeIn: type of data expected: 'char', 'float', 'int', 'bool'
    Returns:
        out: float or 1D array [val1,val2] in the expected format

    '''
    global customFileIN
    raw_input=customFileIN.readline()


    #get text on the right side of the equal sign if there is only one  equal '=' sign
    if len(raw_input.split('='))==2:
        txt=raw_input.split('=')[1].strip()

    #---read the string manually if there is more than one'=' signs: e.g '02400.atmos_average2.{lat =20}'
    elif len(raw_input.split('='))>2:
        current_varfull='';record=False
        for i in range(0,len(raw_input)):
            if record: current_varfull+=raw_input[i]
            if raw_input[i]=='=': record=True
        txt=current_varfull.strip()

    return  filter_input(txt,typeIn)




def read_axis_options(axis_options_txt):
    '''
    Return axis customization options
    Args:
        axis_options_txt: One liner string: 'Axis Options  : lon = [5,8] | lat = [None,None] | cmap = jet'
    Returns:
        Xaxis: X-axis bounds as a numpy array or None if undedefined
        Yaxis: Y-axis bounds as a numpy array or None if undedefined
        custom_line: string, i.e colormap ('jet', 'spectral') or line options, e.g '--r' for dashed red

    '''
    list_txt=axis_options_txt.split(':')[1].split('|')
    #Xaxis: get bound
    txt=list_txt[0].split('=')[1].replace('[','').replace(']','')
    Xaxis=[]
    for i in range(0,len(txt.split(','))):
        if txt.split(',')[i].strip()=='None':
            Xaxis=None
            break
        else:
            Xaxis.append(np.float(txt.split(',')[i].strip()))
    #Yaxis: get bound
    txt=list_txt[1].split('=')[1].replace('[','').replace(']','')
    Yaxis=[]
    for i in range(0,len(txt.split(','))):
        if txt.split(',')[i].strip()=='None':
            Yaxis=None
            break
        else:
            Yaxis.append(np.float(txt.split(',')[i].strip()))
    #Line or colormap
    custom_line=list_txt[2].split('=')[1].strip()
    return Xaxis, Yaxis,custom_line

def split_varfull(varfull):
    '''
    Split  the varfull object into its different components.
    Args:
        varfull: a varfull object, for example 'atmos_average2.zsurf',
                                               '02400.atmos_average2.zsurf'
    Returns:
        sol_array: a sol number e.g 2400 or None if none is provided
        filetype:  file type, i.e 'atmos_average'
        var:       variable of interest, i.e 'zsurf'
        simuID:    int, simulation ID = 2-1= 1 as Python indexes start at zero

    '''

    #---Default case: no sols number is provided, e.g 'atmos_average2.zsurf'--
    #extract variables and file from varfull

    if len(varfull.split('.'))==1+1: # atmos_average2.zsurf'
        sol_array=np.array([None])
        filetypeID=varfull.split('.')[0].strip() #file and ID
        var=varfull.split('.')[1].strip()        #variable name
    #---A sol number is profided, e.g '02400.atmos_average2.zsurf'
    elif len(varfull.split('.'))==2+1:
        sol_array=np.array([int(varfull.split('.')[0].strip())])   #sol number
        filetypeID=varfull.split('.')[1].strip() #file and ID
        var=varfull.split('.')[2].strip()        #variable name
    # in case more than 9 simulations are requested
    if filetypeID[-2:].isdigit():
        simuID=int(filetypeID[-2:])-1
        filetype=filetypeID[:-2]

    elif filetypeID[-1:].isdigit(): #only last digit
        simuID=int(filetypeID[-1])-1
        filetype=filetypeID[:-1]
    else: #no digit, i.e reference simulation
        simuID=0
        filetype=filetypeID
    if simuID<0:
        prRed('*** Error ***, reading %s: only simulations # 1-99 are supported'%(varfull))
    return sol_array,filetype,var,simuID


def remove_whitespace(raw_input):
    '''
    Remove the white space inside an expression. This is different from the '.strip()' method that only remove white spaces at the edges of the string
    Args:
        raw_input: a string, e.g '[atmos_average.temp] +  2'
    Returns:
        processed_input  the string without white spaces, e.g [atmos_average.temp]+2'

    '''
    processed_input=''
    for i in range(0,len(raw_input)):
        if raw_input[i]!=' ': processed_input+=raw_input[i]
    return processed_input

def clean_comma_whitespace(raw_input):
    '''
    Remove the commas and white spaces inside an expression.
    Args:
        raw_input: a string, e.g 'lat=3. ,'
    Returns:
        processed_input  the string without white spaces or commas e.g 'lat=3.lon=2lev=10.'

    '''
    processed_input=''
    for i in range(0,len(raw_input)):
        if raw_input[i]!=',': processed_input+=raw_input[i]
    return remove_whitespace(processed_input)


def get_list_varfull(raw_input):
    '''
    Given an expression object with '[]' return the different variable needed
    Args:
        raw_input: a complex varfull object, for example '2*[atmos_average.temp]+[atmos_average2.ucomp]*1000'
    Returns:
        var_list  a list of variable to load, e.g ['atmos_average.temp', 'atmos_average2.ucomp']

    '''
    var_list=[]
    record = False
    current_name=''
    for i in range(0,len(raw_input)):
        if raw_input[i]==']':
            record=False
            var_list.append(current_name.strip())
            current_name=''
        if record: current_name+=raw_input[i]
        if raw_input[i]=='[': record=True
    return var_list

def get_overwrite_dim_2D(varfull_bracket,plot_type,fdim1,fdim2):
    '''
    Given a single varfull object with '{}' return the new dimensions to overwrite the default dimensions
    Args:
        varfull_bracket: a  varfull object with any of the following atmos_average.temp{lev=10;time=350;lon=155;lat=25} (brackets and semi-colons separated)
        plot_type: the type of plot

    Returns:
        varfull the varfull without brackets: e.g 'atmos_average.temp'
        fdim_out1,fdim_out1: the dimensions to update
    NOTE:
    2D_lon_lat:   fdim1=time
                  fdim2=lev

    2D_lat_press: fdim1=time
                  fdim2=lon

    2D_time_lat: fdim1=lon
                  fdim2=lev

    2D_lon_press: fdim1=time
                  fdim2=lat

    2D_time_press:fdim1=lat
                  fdim2=lon

    2D_lon_time:  fdim1=lat
                  fdim2=lev
    '''
    #Initialization: use the dimension provided in the template
    fdim_out1=fdim1; fdim_out2=fdim2
    varfull_no_bracket=varfull_bracket.split('{')[0].strip()    #left of the '{' character
    overwrite_txt=remove_whitespace(varfull_bracket.split('{')[1][:-1])  #right of the'{' character, with the last '}' removed
    ndim_update=overwrite_txt.count('=') #count the number of '=' in the string
    split_dim=overwrite_txt.split(';');  #split to different blocs e.g 'lat =3.' and 'lon=20'
    if overwrite_txt.count(';')<overwrite_txt.count('=')-1: prYellow("""*** Error:, use semicolon ';' to separate dimensions '{}'""")
    for i in range(0,ndim_update):
        #Check if the requested dimension exists:
        if split_dim[i].split('=')[0] not in ['time','lev','lon','lat']:
            prYellow("""*** Warning***, ignoring dimension: '"""+split_dim[i].split('=')[0]+"""' not recognized: must be 'time','lev','lon' or 'lat'""")

        if plot_type=='2D_lon_lat':
            if split_dim[i].split('=')[0]=='time':fdim_out1=filter_input(split_dim[i].split('=')[1],'float')
            if split_dim[i].split('=')[0]=='lev': fdim_out2=filter_input(split_dim[i].split('=')[1],'float')
        if plot_type=='2D_lat_press':
            if split_dim[i].split('=')[0]=='time':fdim_out1=filter_input(split_dim[i].split('=')[1],'float')
            if split_dim[i].split('=')[0]=='lon': fdim_out2=filter_input(split_dim[i].split('=')[1],'float')
        if plot_type=='2D_time_lat':
            if split_dim[i].split('=')[0]=='lon': fdim_out1=filter_input(split_dim[i].split('=')[1],'float')
            if split_dim[i].split('=')[0]=='lev': fdim_out2=filter_input(split_dim[i].split('=')[1],'float')
        if plot_type=='2D_lon_press':
            if split_dim[i].split('=')[0]=='time':fdim_out1=filter_input(split_dim[i].split('=')[1],'float')
            if split_dim[i].split('=')[0]=='lat': fdim_out2=filter_input(split_dim[i].split('=')[1],'float')
        if plot_type=='2D_time_press':
            if split_dim[i].split('=')[0]=='lat': fdim_out1=filter_input(split_dim[i].split('=')[1],'float')
            if split_dim[i].split('=')[0]=='lon': fdim_out2=filter_input(split_dim[i].split('=')[1],'float')
        if plot_type=='2D_lon_time':
            if split_dim[i].split('=')[0]=='lat': fdim_out1=filter_input(split_dim[i].split('=')[1],'float')
            if split_dim[i].split('=')[0]=='lev': fdim_out2=filter_input(split_dim[i].split('=')[1],'float')

    # NOTE: filter_input() convert the text '3' or '4,5' to real variable, e.g numpy.array([3.]) numpy.array([4.,5.])
    return varfull_no_bracket, fdim_out1, fdim_out2

def get_overwrite_dim_1D(varfull_bracket,t_in,lat_in,lon_in,lev_in):
    '''
    Given a single varfull object with '{}' return the new dimensions to overwrite the default dimensions
    Args:
        varfull_bracket: a  varfull object with any of the following atmos_average.temp{lev=10;time=350;lon=155;lat=25}
        t_in,lat_in,lon_in,lev_in: the variables as defined by self.t ,self.lat,self.lon,self.lev

    Returns:
        varfull the varfull without brackets: e.g 'atmos_average.temp'
        t_out,lat_out,lon_out,lev_out: the dimensions to update
    NOTE:

    '''
    #Initialization: use the dimension provided in the template
    t_out=t_in; lat_out=lat_in; lon_out=lon_in;lev_out=lev_in
    varfull_no_bracket=varfull_bracket.split('{')[0].strip()    #left of the '{' character
    overwrite_txt=remove_whitespace(varfull_bracket.split('{')[1][:-1])  #right of the'{' character, with the last '}' removed
    ndim_update=overwrite_txt.count('=') #count the number of '=' in the string
    split_dim=overwrite_txt.split(';');  #split to different blocs e.g 'lat =3.' and 'lon=20'
    for i in range(0,ndim_update):
        #Check if the requested dimension exists:
        if split_dim[i].split('=')[0] not in ['time','lev','lon','lat']:
            prYellow("""*** Warning***, ignoring dimension: '"""+split_dim[i].split('=')[0]+"""' not recognized: must be 'time','lev','lon' or 'lat'""")


        if split_dim[i].split('=')[0]=='time':t_out=  filter_input(split_dim[i].split('=')[1],'float')
        if split_dim[i].split('=')[0]=='lat': lat_out=filter_input(split_dim[i].split('=')[1],'float')
        if split_dim[i].split('=')[0]=='lon': lon_out=filter_input(split_dim[i].split('=')[1],'float')
        if split_dim[i].split('=')[0]=='lev': lev_out=filter_input(split_dim[i].split('=')[1],'float')

    # NOTE: filter_input() convert the text '3' or '4,5' to real variable, e.g numpy.array([3.]) numpy.array([4.,5.])
    return varfull_no_bracket, t_out,lat_out,lon_out,lev_out


def create_exec(raw_input,varfull_list):
    expression_exec=raw_input
    for i in range(0,len(varfull_list)):
        swap_txt='['+varfull_list[i]+']'
        expression_exec=expression_exec.replace(swap_txt,'VAR[%i]'%(i))
    return expression_exec

def fig_layout(subID,nPan):
    '''
    Return figure layout
    Args:
        subID:    integer, current subplot number
        nPan : integer, number of pannels desired on the figure up 36 (6x6 pannel)
    Returns:
        out: tuple with approriate layout: plt.subplot(nrows=out[0],ncols=out[1],plot_number=out[2])
    '''
    out=list((0,0,0)) #initialization

    if nPan==1:out[0:2]=(1,1) #nrow,ncol
    if nPan==2:out[0:2]=(1,2)
    if nPan==3 or nPan==4 :out[0:2]=(2,2)
    if nPan==5 or nPan==6 :out[0:2]=(2,3)
    if nPan==7 or nPan==8 :out[0:2]=(2,4)
    if nPan==9:            out[0:2]=(3,3)
    if 10<=nPan<=12:out[0:2]=(3,4)
    if 13<=nPan<=16:out[0:2]=(4,4)
    if 17<=nPan<=20:out[0:2]=(4,5)
    if 21<=nPan<=25:out[0:2]=(5,5)
    if 26<=nPan<=30:out[0:2]=(5,6)
    if 30<=nPan<=36:out[0:2]=(6,6)

    out[2]=subID #finally the current plot

    return out

def make_template():
    global customFileIN # (will be modified)
    global current_version
    newname=output_path+'/Custom.in'
    newname= create_name(newname)

    customFileIN=open(newname,'w')

    lh="""# """ #Add a line header. This is primary use to change the coloring of the text when using vim
    #==============Create header with instructions, and add the version number to the title====
    customFileIN.write("===================== |MarsPlot V%s|===================\n"%(current_version))
    if parser.parse_args().template: #Additional instructions if requested
        customFileIN.write(lh+"""INSTRUCTIONS:\n""")
        customFileIN.write(lh+"""> Find the matching  template for the desired plot type. Do not edit any labels left of any '=' sign \n""")
        customFileIN.write(lh+"""> Duplicate/remove any of the <<<< blocks>>>>, skip by setting <<<< block = False >>>> \n""")
        customFileIN.write(lh+"""> 'True', 'False' and 'None' are capitalized. Do not use quotes '' anywhere in this file \n""")
        customFileIN.write(lh+"""> Cmin, Cmax define the colorbar range. Scientific notation (e.g. 1e-6, 2e3) is supported \n""")
        customFileIN.write(lh+"""> 'Level' refers to 'level'[Pa], 'pfull'[Pa] or 'zgrid' [m] depending on the type of *.nc file \n""")
        customFileIN.write(lh+"""FREE DIMENSIONS:\n""")
        customFileIN.write(lh+"""> Use 'Dimension = 55.' to set to the closest value\n""")
        customFileIN.write(lh+"""> Use 'Dimension = all' to average over all values\n""")
        customFileIN.write(lh+"""> Use 'Dimension = -55.,55.' to get the average between -55. and 55. \n""")
        customFileIN.write(lh+"""> 'None' refers to the default setting for that Dimension: \n""")
        customFileIN.write(lh+"""    -A) time  = instant time step at Nt (i.e last timestep) \n""")
        customFileIN.write(lh+"""    -B) lev = sfc (i.e, Nz for *.nc files and 0 for *_pstd.nc files) \n""")
        customFileIN.write(lh+"""    -C) lat   = equator slice \n""")
        customFileIN.write(lh+"""    -D) lon   = 'all', i.e zonal average over all longitudes\n""")
        customFileIN.write(lh+"""> Overwrite the dimension using atmos_average.temp{time = 100 ; lev= 5.; lon= all ; lat=45} Use brackets '{}' and SEMI-COLONS ';'\n""")
        customFileIN.write(lh+""">    Units must be the same as the free dimension block, i.e time [Ls], lev [Pa], lon [deg], and lat [deg]   \n""")
        customFileIN.write(lh+"""TIME SERIES AND 1D PLOTS:\n""")
        customFileIN.write(lh+"""> Use 'Dimension = AXIS' to set the varying axis\n""")
        customFileIN.write(lh+"""> The other free dimensions accept value, 'all' or valmin,valmax as above\n""")
        customFileIN.write(lh+"""AXIS OPTIONS:\n""")
        customFileIN.write(lh+"""Set the x-axis and y-axis limits in the figure units. All Matplolib styles are supported:\n""")
        customFileIN.write(lh+"""> 'cmap' changes the colormap: 'jet' (winds), 'spectral' (temperature), 'bwr' (diff plot)\n""")
        customFileIN.write(lh+"""> 'line' sets the line style:  '-r' (solid red), '--g' (dashed green), '-ob' (solid & blue markers)\n""")
        customFileIN.write(lh+"""KEYWORDS:\n""")
        customFileIN.write(lh+"""> 'HOLD ON' [blocks of figures] 'HOLD OFF' groups the figures as a multi-pannel page\n""")
        customFileIN.write(lh+"""> [line plot 1] 'ADD LINE' [line plot 2] adds similar 1D-plots on the same figure)\n""")
        customFileIN.write(lh+"""> 'START' and (optionally) 'STOP' can be used to conveniently skip plots below. Use '#' to add comments. \n""")
        customFileIN.write(lh+"""ALGEBRA AND CROSS-SIMULATIONS PLOTS:\n""")
        customFileIN.write(lh+"""Use 'N>' to add a Nth simulation with matching timesteps to the <<< Simulations >>> block  \n""")
        customFileIN.write(lh+"""Use full path: '2> /u/akling/FV3/verona/simu2/history' Empty fields are ignored, comment out with '#' \n""")
        customFileIN.write(lh+"""A variable 'var' in a 'XXXXX.file.nc' from this simulation is accessed using 'XXXXX.fileN.var' syntax \n""")
        customFileIN.write(lh+"""Encompass raw outputs with square brackets '[]' for element-wise operations, e.g: \n""")
        customFileIN.write(lh+"""> '[fixed.zurf]/(10.**3)'                               (convert topography from [m] to [km])\n""")
        customFileIN.write(lh+"""> '[atmos_average.taudust_IR]/[atmos_average.ps]*(610)' (normalize the dust opacity)     \n""")
        customFileIN.write(lh+"""> '[atmos_average.temp]-[atmos_average2.temp]'    (temp. difference between ref simu and simu 2)\n""")
        customFileIN.write(lh+"""> '[atmos_average.temp]-[atmos_average.temp{lev=10}]'    (temp. difference between the default (near surface) and the 10 Pa level\n""")

        customFileIN.write(lh+"""        Supported expressions are: sqrt, log, exp, min, max, mean\n""")
    customFileIN.write("<<<<<<<<<<<<<<<<<<<<<< Simulations >>>>>>>>>>>>>>>>>>>>>\n")
    customFileIN.write("ref> None\n")
    customFileIN.write("2>\n")
    customFileIN.write("3>\n")
    customFileIN.write("=======================================================\n")
    customFileIN.write("START\n")
    customFileIN.write("\n") #new line
    #===============================================================
    #For the default list of figures in main(), create a  template.
    for i in range(0,len(objectList)):
        if objectList[i].subID==1 and objectList[i].nPan>1: customFileIN.write('HOLD ON\n')
        objectList[i].make_template()
        customFileIN.write('\n')
        if objectList[i].nPan>1 and objectList[i].subID==objectList[i].nPan: customFileIN.write('HOLD OFF\n')

        #Separate the empty templates
        if  i==1:
            customFileIN.write("""#=========================================================================\n""")
            customFileIN.write("""#================== Empty Templates (set to False)========================\n""")
            customFileIN.write("""#========================================================================= \n""")
            customFileIN.write(""" \n""")



    customFileIN.close()

    # NAS system only: set group permission to the file and print completion message
    give_permission(newname)
    print(newname +' was created ')
    #---

def give_permission(filename):
    # NAS system only: set group permission to the file
    try:
        subprocess.check_call(['setfacl -v'],shell=True,stdout=open(os.devnull, "w"),stderr=open(os.devnull, "w")) #catch error and standard output
        cmd_txt='setfacl -R -m g:s0846:r '+filename
        subprocess.call(cmd_txt,shell=True)
    except subprocess.CalledProcessError:
        pass

def namelist_parser(Custom_file):
    '''
    Parse a template
    Args:
        Custom_file: full path to Custom.in file
    Actions:
        Update  global variableFigLayout, objectList
    '''
    global objectList
    global customFileIN
    global input_paths
    # A Custom file is provided, flush the default figures defined in main()
    #---
    objectList=[] #all individual plots

    pannelList=[] #list of pannels
    subplotList=[] #layout of figures
    addLineList=[] #add several line plot on the same graphs
    nobj=0        #number for the object: e.g 1,[2,3],4... with 2 & 3 plotted as a two pannels plot
    nPannel=1        #number of pannels ploted along this object, e.g: '1' for object #1 and '2' for the objects #2 and #3
    subplotID=1  #subplot ID per object: e.g '1' for object #1, '1' for object #2 and '2' for object #3
    holding=False
    addLine=False
    addedLines=0  #line plots
    npage=0       #plot number at the begining of a new page (e.g 'HOLD ON')

    customFileIN=open(Custom_file,'r')
    #===Get version in the header====
    version=np.float(customFileIN.readline().split('|')[1].strip().split('V')[1].strip())
    # Check if the main versions are compatible,  (1.1 and 1.2 are OK but not 1.0 and 2.0)
    if np.int(version)!=np.int(current_version):
         prYellow('*** Warning ***')
         prYellow('Using MarsPlot V%s but Custom.in template is depreciated (using V%s)'%(current_version,version))
         prYellow('***************')

    #==========Skip the header======
    while (customFileIN.readline()[0]!='<'):
        pass
    #==========Read simulations in <<<<<<<<<< Simulations >>>>>> ======
    while True:
        line=customFileIN.readline()
        if line[0]=='#': #skip comment
            pass
        else:
            if line[0]=='=': break #finished reading
            # Special case reference simulation
            if line.split('>')[0]=='ref':
                # if is different from default, overwrite it
                if line.split('>')[1].strip()!='None':
                    input_paths[0]=line.split('>')[1].strip()
            else:
                if '>' in line: #line contains '>' symbol
                    if line.split('>')[1].strip(): #line exist and is not blank
                        input_paths.append(line.split('>')[1].strip())

    #===========skip lines until the kweywor 'START' is found================
    nsafe=0 #initialize counter for safety
    while True and nsafe<1000:
        line=customFileIN.readline()
        if line.strip()=='START':break
        nsafe+=1
    if nsafe==1000:prRed(""" Custom.in is missing a 'START' keyword after the '=====' simulation block""")

    #=============Start reading the figures=================
    while True:
        line=customFileIN.readline()

        if not line or line.strip()=='STOP':
            break #Reached End Of File

        if line.strip()=='HOLD ON':
            holding=True
            subplotID=1
        #adding a 1D plot to an existing line plot
        if line.strip()=='ADD LINE':
            addLine=True

        if line[0]== '<': #If new figure
            figtype,boolPlot=get_figure_header(line)
            if boolPlot : #only if we want to plot the field
            #Add object to the list
                if figtype =='Plot 2D lon X lat'   :objectList.append(Fig_2D_lon_lat())
                if figtype =='Plot 2D time X lat'  :objectList.append(Fig_2D_time_lat())
                if figtype =='Plot 2D lat X press' :objectList.append(Fig_2D_lat_press())
                if figtype =='Plot 2D lon X press' :objectList.append(Fig_2D_lon_press())
                if figtype =='Plot 2D time X press':objectList.append(Fig_2D_time_press())
                if figtype =='Plot 2D lon X time'  :objectList.append(Fig_2D_lon_time())
                if figtype =='Plot 1D'             :objectList.append(Fig_1D())
                objectList[nobj].read_template()
                nobj+=1
                #====debug only===========
                #print('------nobj=',nobj,' npage=',npage,'-------------------')


                #===================
                if holding and not addLine:
                    subplotList.append(subplotID)
                    pannelList.append(subplotID)
                    subplotID+=1
                    #Add +1 pannel to all plot in current page
                    for iobj in range(npage,nobj-1):
                        pannelList[iobj]+=1

                elif holding and addLine:
                    #Do not update  subplotID if we are adding lines
                    subplotList.append(subplotID-1)
                    pannelList.append(subplotID-1)



                else :
                    #We are not holding: there is a single pannel per page and we reset the page counter
                    pannelList.append(1)
                    subplotList.append(1)
                    npage=nobj

                #====================

                if addLine:
                    addedLines+=1
                    addLineList.append(addedLines)
                else:
                     addLineList.append(0) #no added lines
                     addedLines=0 #reset line counter

                #====debug only====
                #for ii in range(0,len(   subplotList)):
                #    prCyan('[X,%i,%i,%i]'%(subplotList[ii],pannelList[ii],addLineList[ii]))
                #=================


                #============Depreciated=(old way to attribute the plot numbers without using npage)=============
                # if holding:
                #     subplotList.append(subplotID-addedLines)
                #     pannelList.append(subplotID-addedLines)
                #     if not addLine:
                #         # add +1 to the number of pannels for the previous plots
                #         n=1
                #         while n<=subplotID-1:
                #             pannelList[nobj-n-1]+=1 #print('editing %i pannels, now %i'%(subplotID-1,nobj-n-1))
                #             n+=1
                #     subplotID+=1
                # else :
                #     pannelList.append(1)
                #     subplotList.append(1)
                #========================================================


            addLine=False #reset after reading each block
        if line.strip()=='HOLD OFF':
            holding=False
            subplotID=1
            npage=nobj


    #Make sure we are not still holding figures
    if holding:
        prRed('*** Error ***')
        prRed("""Missing 'HOLD OFF' statement in """+Custom_file)
        exit()
    #Make sure we are not still holding figures
    if addLine:
        prRed('*** Error ***')
        prRed("""Cannot have 'ADD LINE' after the last figure in """+Custom_file)
        exit()
    #Finished reading the file, attribute the right number of figure and pannels for each plot
    #print('=======Summary=========')
    for i in range(0,nobj):
        objectList[i].subID=subplotList[i]
        objectList[i].nPan=pannelList[i]
        objectList[i].addLine=addLineList[i]
        #==debug only====
        #prPurple('%i:[%i,%i,%i]'%(i,objectList[i].subID,objectList[i].nPan,objectList[i].addLine))
    customFileIN.close()


def get_figure_header(line_txt):
    '''
    This function return the type of a figure and tells us
    Args:
        line_txt: string, gigure header from Custom.in, i.e '<<<<<<<<<<<<<<| Plot 2D lon X lat = True |>>>>>>>>>>>>>'
    Returns:
        figtype : string, figure type, i.e:  Plot 2D lon X lat
        boolPlot: boolean, is the plot wanted?
    '''
    line_cmd=line_txt.split('|')[1].strip() #Plot 2D lon X lat = True
    figtype=line_cmd.split('=')[0].strip()  #Plot 2D lon X lat
    boolPlot=line_cmd.split('=')[1].strip()=='True' # Return True
    return figtype, boolPlot

#======================================================
#                  FILE SYSTEM UTILITIES
#======================================================



def get_Ncdf_num():
    '''
    Get the sol numbers of all the netcdf files in directory
    This test is based on the existence of a least one  XXXXX.fixed.nc in the current directory.
    Args:
        None
    Returns:
        Ncdf_num: a sorted array of sols
    '''
    list_dir=os.listdir(input_paths[0])
    avail_fixed = [k for k in list_dir if '.fixed.nc' in k] #e.g. '00350.fixed.nc', '00000.fixed.nc'
    list_num = [item[0:5] for item in avail_fixed]          #remove .fixed.nc, e.g. '00350', '00000'
    Ncdf_num=np.sort(np.asarray(list_num).astype(np.float)) # transform to array, e.g. [0, 350]
    if Ncdf_num.size==0:
        print("No XXXXX.fixed.nc detected in "+input_paths[0])
        raise SystemExit #Exit cleanly
    return Ncdf_num

def select_range(Ncdf_num,bound):
    '''
    Args:
        Ncdf_num:  a sorted array of sols
        bound: a integer representing a date (e.g. 0350) or an array containing the sol bounds (e.g [min max])
    Returns:
        Ncdf_num: a sorted array of sols within the prescribed bounds
    '''
    bound=np.array(bound)
    if bound.size==1:
        Ncdf_num=Ncdf_num[Ncdf_num==bound]
        if Ncdf_num.size==0:
            prRed('*** Error ***')
            prRed("File %05d.fixed.nc not detected"%(bound))
            exit()
    elif bound.size==2:
        Ncdf_num=Ncdf_num[Ncdf_num>=bound[0]]
        Ncdf_num=Ncdf_num[Ncdf_num<=bound[1]]
        if Ncdf_num.size==0:
            prRed('*** Error ***')
            prRed("No XXXXX.fixed.nc detected between sols [%05d-%05d] please check date range"%(bound[0],bound[1]))
            exit()
    return Ncdf_num

def create_name(root_name):
    '''
    Create a file name based on its existence in the current directory.
    Args:
        root_name: desired name for the file: "/path/custom.in" or "/path/figure.png"
    Returns:
        new_name: new name if the file already exists: "/path/custom_01.in" or "/path/figure_01.png"
    '''
    n=1
    len_ext=len(root_name.split('.')[-1]) #get extension lenght (e.g 2 for *.nc, 3 for *.png)
    ext=root_name[-len_ext:]              #get extension
    new_name=root_name #initialization
    #if example.png already exist, create example_01.png
    if os.path.isfile(new_name):
        new_name=root_name[0:-(len_ext+1)]+'_%02d'%(n)+'.'+ext
    #if example_01.png already exist, create example_02.png etc...
    while os.path.isfile(root_name[0:-(len_ext+1)]+'_%02d'%(n)+'.'+ext):
        n=n+1
        new_name=root_name[0:-(len_ext+1)]+'_%02d'%(n)+'.'+ext
    return new_name


def path_to_template(custom_name):
    '''
    Create a file name based on its existence in the current directory.
    Args:
        custom_name: custom file name, accepted formats are my_custom or my_custom.in
    Returns:
        full_path: full_path to /u/user/FV3/templates/my_custom.in

         If file not found,try:/lou/s2n/mkahre/MCMC/analysis/working/templates/my_custom.in
    '''

    #local_dir=os.path.expanduser("~") +'/FV3/templates'
    local_dir=sys.prefix+'/mars_templates'
    shared_dir='/lou/s2n/mkahre/MCMC/analysis/working/templates'

    #---
    custom_name=custom_name[0] #convert the 1-element list to a string
    if custom_name[-3:]!='.in':  custom_name=custom_name+'.in'#add extension if not provided
    #first look in  '~/FV3/templates'
    if not os.path.isfile(local_dir+'/'+custom_name):
        #then look in  '/lou/s2n/mkahre/MCMC/analysis/working/templates'
        if not os.path.isfile(shared_dir+'/'+custom_name):
            prRed('*** Error ***')
            prRed('File '+custom_name+' not found in '+local_dir+' ... nor in : \n                          '+shared_dir)
            # if a local ~/FV3/templates path does not exist, suggest to create it
            if not os.path.exists(local_dir):
                prYellow('Note: directory: ~/FV3/templates'+' does not exist, create it with:')
                prCyan('mkdir '+local_dir)
            exit()
        else:
            return shared_dir+'/'+custom_name
    else:
        return local_dir+'/'+custom_name


def progress(k,Nmax,txt='',success=True):
    """
    Display a progress bar to monitor heavy calculations.
    Args:
        k: current iteration of the outer loop
        Nmax: max iteration of the outer loop
    Returns:
        Running... [#---------] 10.64 %
    """
    import sys
    progress=float(k)/Nmax
    barLength = 10 # Modify this to change the length of the progress bar
    block = int(round(barLength*progress))
    bar = "[{0}]".format( "#"*block + "-"*(barLength-block))
    #bar = "Running... [\033[96m{0}\033[00m]".format( "#"*block + "-"*(barLength-block))  #add color
    if success==True:
        #status="%i %% (%s)"%(100*progress,txt) #no color
        status="%3i %% \033[92m(%s)\033[00m"%(100*progress,txt)  #green
    elif success==False:
        status="%3i %% \033[91m(%s)\033[00m"%(100*progress,txt) #red
    elif success==None:
        status="%3i %% (%s)"%(100*progress,txt) #red
    text='\r'+bar+status+'\n'
    sys.stdout.write(text)
    if not debug: sys.stdout.flush()

#======================================================
#                  FIGURE DEFINITIONS
#======================================================
class Fig_2D(object):
    # Parent class for 2D figures
    def __init__(self,varfull='fileYYY.XXX',doPlot=False,varfull2=None):

        self.title=None
        self.varfull=varfull
        self.range=None
        self.fdim1=None
        self.fdim2=None
        self.varfull2=varfull2
        # Logic
        self.doPlot=doPlot
        self.plot_type=self.__class__.__name__[4:]

        #Extract filetype, variable, and simulation ID (initialization only)
        self.sol_array,self.filetype,self.var,self.simuID=split_varfull(self.varfull)
        if self.varfull2: self.sol_array2,self.filetype2,self.var2,self.simuID2=split_varfull(self.varfull2)

        #Multi pannel
        self.nPan=1
        self.subID=1
        #Annotation for free dimensions
        self.fdim_txt=''
        self.success=False
        self.addLine=False
        #Axis options

        self.Xlim=None
        self.Ylim=None
        self.axis_opts='jet'


    def make_template(self,plot_txt,fdim1_txt,fdim2_txt,Xaxis_txt,Yaxis_txt):
        customFileIN.write("<<<<<<<<<<<<<<| {0:<15} = {1} |>>>>>>>>>>>>>\n".format(plot_txt,self.doPlot))
        customFileIN.write("Title          = %s\n"%(self.title))             #1
        customFileIN.write("Main Variable  = %s\n"%(self.varfull))           #2
        customFileIN.write("Cmin, Cmax     = %s\n"%(self.range))             #3
        customFileIN.write("{0:<15}= {1}\n".format(fdim1_txt,self.fdim1))    #4
        customFileIN.write("{0:<15}= {1}\n".format(fdim2_txt,self.fdim2))    #4
        customFileIN.write("2nd Variable   = %s\n"%(self.varfull2))          #6
        customFileIN.write("Axis Options  : {0} = [None,None] | {1} = [None,None] | cmap = jet \n".format(Xaxis_txt,Yaxis_txt))

    def read_template(self):
        self.title= rT('char')                   #1
        self.varfull=rT('char')                  #2
        self.range=rT('float')                   #3
        self.fdim1=rT('float')                   #4
        self.fdim2=rT('float')                   #5
        self.varfull2=rT('char')                 #6
        self.Xlim,self.Ylim,self.axis_opts=read_axis_options(customFileIN.readline())     #7

        #Various sanity checks
        if self.range and len(np.atleast_1d(self.range))==1:
            prYellow('*** Warning ***, In plot %s, Cmin, Cmax must be two values, resetting to default'%(self.varfull))
            self.range=None

        #Update the variable after reading template
        #self.sol_array,self.filetype,self.var,self.simuID=split_varfull(self.varfull)
        #if self.varfull2: self.sol_array2,self.filetype2,self.var2,self.simuID2=split_varfull(self.varfull2)


    def prep_file(self,var_name,file_type,simuID,sol_array):
        global input_paths
        global Ncdf_num #global variable that holds  the different ols number, e.g [1500,2400]
        # A specific sol was requested, e.g [2400]
        if sol_array:
            Sol_num_current=sol_array
        else : #no sol requested, use default as provided by MarsPlot Custom.in -d sol
            Sol_num_current =Ncdf_num

        nfiles=len(Sol_num_current)

        file_list = [None]*nfiles #initialize the list


        #Loop over the requested time steps
        for i in range(0,nfiles):
            file_list[i] = input_paths[simuID]+'/%05d.'%(Sol_num_current[i])+file_type+'.nc'
            check_file_tape(file_list[i],abort=False)

        if file_type=='fixed': # XXXX.fixed.nc does not have an aggregation dimension so we use Dataset
            f=Dataset(file_list[0], 'r')
        else:
            f=MFDataset(file_list, 'r') #use MFDataset instead

        var_info=f.variables[var_name].long_name+' ['+ f.variables[var_name].units+']'
        dim_info=f.variables[var_name].dimensions
        dims=f.variables[var_name].shape
        return f, var_info,dim_info, dims

    def data_loader_2D(self,varfull,plot_type):

        #Simply plot one of the variable in the file
        if not '[' in varfull:
            #---If overwriting dimensions, get the new dimensions and trim varfull from the '{lev=5.}' part
            if '{' in varfull :
                varfull,fdim1_extract,fdim2_extract=get_overwrite_dim_2D(varfull,plot_type,self.fdim1,self.fdim2)
                # fdim1_extract,fdim2_extract constains the dimensions to overwrite is '{}' are provided of the default self.fdim1, self.fdim2  otherwise
            else: # no '{ }' use to overwrite the dimensions, copy the plots' defaults
                fdim1_extract,fdim2_extract=self.fdim1, self.fdim2
            sol_array,filetype,var,simuID=split_varfull(varfull)
            xdata,ydata,var,var_info=self.read_NCDF_2D(var,filetype,simuID,sol_array,plot_type,fdim1_extract,fdim2_extract)
        #Realize a operation on the variables
        else:
            VAR=[]
            # Extract individual variables and prepare for execution
            varfull=remove_whitespace(varfull)
            varfull_list=get_list_varfull(varfull)
            #Initialize list of requested dimensions;
            fdim1_list=[None]*len(varfull_list)
            fdim2_list=[None]*len(varfull_list)
            expression_exec=create_exec(varfull,varfull_list)



            for i in range(0,len(varfull_list)):
                #---If overwriting dimensions, get the new dimensions and trim varfull from the '{lev=5.}' part
                if '{' in varfull_list[i] :
                    varfull_list[i],fdim1_list[i],fdim2_list[i]=get_overwrite_dim_2D(varfull_list[i],plot_type,self.fdim1,self.fdim2)
                else: # no '{ }' use to overwrite the dimensions, copy the plots' defaults
                    fdim1_list[i],fdim2_list[i]=self.fdim1, self.fdim2

                sol_array,filetype,var,simuID=split_varfull(varfull_list[i])
                xdata,ydata,temp,var_info=self.read_NCDF_2D(var,filetype,simuID,sol_array,plot_type,fdim1_list[i],fdim2_list[i])
                VAR.append(temp)
            var_info=varfull
            var=eval(expression_exec)

        return xdata,ydata,var,var_info

    def read_NCDF_2D(self,var_name,file_type,simuID,sol_array,plot_type,fdim1,fdim2):
        f, var_info,dim_info, dims=self.prep_file(var_name,file_type,simuID,sol_array)

        #If self.fdim is empty, add the variable (do only once)
        add_fdim=False
        if not self.fdim_txt.strip():add_fdim=True

        #Initialize dimensions (These are in all the .nc files)

        lat=f.variables['lat'][:];lati=np.arange(0,len(lat))
        lon=f.variables['lon'][:];loni=np.arange(0,len(lon))

        #Load variable depending on the requested free dimensions

        #======static======= , ignore level and time dimension
        if dim_info==(u'lat', u'lon'):
            var=f.variables[var_name][lati,loni]
            f.close()
            return lon,lat,var,var_info
        #if plot_type=='2D_lon_lat':
         #   return lon,lat,var,var_info

        #======time,lat,lon=======
        if f.variables[var_name].dimensions==(u'time', u'lat', u'lon'):
        #Initialize dimension
            t=f.variables['time'][:];Ls=np.squeeze(f.variables['areo'][:]);ti=np.arange(0,len(t))
            t_stack=np.vstack((t,Ls)) #stack the time and ls array as one variable

            if plot_type=='2D_lon_lat': ti,temp_txt =get_time_index(fdim1,Ls,t)
            if plot_type=='2D_time_lat':loni,temp_txt =get_lon_index(fdim1,lon)
            if plot_type=='2D_lon_time':lati,temp_txt =get_lat_index(fdim1,lat)

            if add_fdim:self.fdim_txt+=temp_txt
            #Extract data and close file
            var=f.variables[var_name][ti,lati,loni].reshape(len(np.atleast_1d(ti)),len(np.atleast_1d(lati)),len(np.atleast_1d(loni)))
            f.close()
            #Return data
            if plot_type=='2D_lon_lat': return lon,lat,np.mean(var,axis=0),var_info #time average
            if plot_type=='2D_time_lat':return t_stack,lat,np.mean(var,axis=2).T,var_info #transpose, Xdim must be in last column of var
            if plot_type=='2D_lon_time':return lon,t_stack,np.mean(var,axis=1),var_info


        #======time,level,lat,lon=======
        if (dim_info==(u'time', u'pfull', u'lat', u'lon')
           or dim_info==(u'time', u'pstd', u'lat', u'lon')
           or dim_info==(u'time', u'zgrid', u'lat', u'lon')):

            #Initialize dimensions
            if dim_info[1]=='pfull': levs=f.variables[dim_info[1]][:] 
            if dim_info[1]=='pstd': levs=f.variables[dim_info[1]][:] 
            if dim_info[1]=='zgrid': levs=     f.variables[dim_info[1]][:] # meters
            zi=np.arange(0,len(levs))
            t=f.variables['time'][:];Ls=np.squeeze(f.variables['areo'][:]);ti=np.arange(0,len(t))
            t_stack=np.vstack((t,Ls)) #stack the time and ls array as one variable

            if plot_type=='2D_lon_lat':
                ti,temp_txt =get_time_index(fdim1,Ls,t)
                if add_fdim:self.fdim_txt+=temp_txt
                zi,temp_txt =get_level_index(fdim2,levs)
                if add_fdim:self.fdim_txt+=temp_txt

            if plot_type=='2D_time_lat':
                loni,temp_txt =get_lon_index(fdim1,lon)
                if add_fdim:self.fdim_txt+=temp_txt
                zi,temp_txt =get_level_index(fdim2,levs)
                if add_fdim:self.fdim_txt+=temp_txt

            if plot_type=='2D_lat_press':
                ti,temp_txt =get_time_index(fdim1,Ls,t)
                if add_fdim:self.fdim_txt+=temp_txt
                loni,temp_txt =get_lon_index(fdim2,lon)
                if add_fdim:self.fdim_txt+=temp_txt

            if plot_type=='2D_lon_press':
                ti,temp_txt =get_time_index(fdim1,Ls,t)
                if add_fdim:self.fdim_txt+=temp_txt
                lati,temp_txt =get_lat_index(fdim2,lat)
                if add_fdim:self.fdim_txt+=temp_txt


            if plot_type=='2D_time_press':
                lati,temp_txt =get_lat_index(fdim1,lat)
                if add_fdim:self.fdim_txt+=temp_txt
                loni,temp_txt =get_lon_index(fdim2,lon)
                if add_fdim:self.fdim_txt+=temp_txt

            if plot_type=='2D_lon_time':
                lati,temp_txt =get_lat_index(fdim1,lat)
                if add_fdim:self.fdim_txt+=temp_txt
                zi,temp_txt =get_level_index(fdim2,levs)
                if add_fdim:self.fdim_txt+=temp_txt


            var=f.variables[var_name][ti,zi,lati,loni].reshape(len(np.atleast_1d(ti)),\
                                                               len(np.atleast_1d(zi)),\
                                                               len(np.atleast_1d(lati)),\
                                                               len(np.atleast_1d(loni)))
            f.close()
            #(u'time', u'pfull', u'lat', u'lon')
            if plot_type=='2D_lon_lat': return  lon,   lat,   np.mean(np.mean(var,axis=1),axis=0),var_info
            if plot_type=='2D_time_lat':return t_stack,lat,   np.mean(np.mean(var,axis=1),axis=2).T,var_info #transpose
            if plot_type=='2D_lat_press':return  lat, levs,   np.mean(np.mean(var,axis=3),axis=0),var_info
            if plot_type=='2D_lon_press':return  lon, levs,   np.mean(np.mean(var,axis=2),axis=0),var_info
            if plot_type=='2D_time_press':return t_stack,levs,np.mean(np.mean(var,axis=3),axis=2).T,var_info #transpose
            if plot_type=='2D_lon_time':  return lon,t_stack, np.mean(np.mean(var,axis=2),axis=1),var_info



        #stack the time and ls array as one variable
        #t_stack=np.vstack((t,Ls))

        #return lon,lat,var,var_info


    def make_title(self,var_info,xlabel,ylabel):
        if self.title:
            plt.title(self.title,fontsize=label_size-self.nPan//2)
        else:
            plt.title(var_info+'\n'+self.fdim_txt[1:],fontsize=12-self.nPan//2) #we remove the first coma ',' of fdim_txt to print to the new line
        plt.xlabel(xlabel,fontsize=label_size-self.nPan//2)
        plt.ylabel(ylabel,fontsize=label_size-self.nPan//2)


    def exception_handler(self,e,ax):
        if debug:raise
        sys.stdout.write("\033[F");sys.stdout.write("\033[K")#cursor up one line, then clear the whole line previous output
        prYellow('*** Warning *** %s'%(e))
        ax.text(0.5, 0.5, 'ERROR:'+str(e),horizontalalignment='center',verticalalignment='center', \
            bbox=dict(boxstyle="round",ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8),),\
            transform=ax.transAxes,wrap=True,fontsize=16)

    def fig_init(self):
        #create figure
        out=fig_layout(self.subID,self.nPan)
        if self.subID==1:
            fig= plt.figure(facecolor='white',figsize=(pixel_width/my_dpi, pixel_width/1.4/my_dpi)) #create figure if 1st pannel, 1.4 is ratio (16:9 screen would be 1.77)
            #plt.suptitle(simulation_name)

        ax = plt.subplot(out[0],out[1],out[2]) #nrow,ncol,subID
        ax.patch.set_color('.1') #Nan are grey
        return ax

    def fig_save(self):
        #save the figure
        if  self.subID==self.nPan: #Last subplot
            if  self.subID==1: #1 plot
                if not '[' in self.varfull:
                    sensitive_name=self.varfull.split('{')[0].strip()  #add split '{' in case varfull contains layer, does not do anything otherwise
                    # varfull is a complex expression
                else:
                    sensitive_name='expression_'+get_list_varfull(self.varfull)[0].split('{')[0].strip()
            else: #multi pannel
                sensitive_name='multi_pannel'
            plt.tight_layout()
            self.fig_name=output_path+'/plots/'+sensitive_name+'.'+out_format
            self.fig_name=create_name(self.fig_name)
            plt.savefig(self.fig_name,dpi=my_dpi )
            if out_format!="pdf":print("Saved:" +self.fig_name)


    def filled_contour(self,xdata,ydata,var):
        if self.range:
            plt.contourf(xdata, ydata,var,np.linspace(self.range[0],self.range[1],levels),extend='both',cmap=self.axis_opts)
        else:
            plt.contourf(xdata, ydata,var,levels,cmap=self.axis_opts)
        cbar=plt.colorbar(orientation='horizontal',aspect=50)
        cbar.ax.tick_params(labelsize=label_size-self.nPan//2) #shrink the colorbar label as the number of subplot increase

        #self.nPan

    def solid_contour(self,xdata,ydata,var):
       np.seterr(divide='ignore', invalid='ignore') #prevent error message when making contour
       CS=plt.contour(xdata, ydata,var,8,colors='k',linewidths=3)
       plt.clabel(CS, inline=1, fontsize=14,fmt='%g')


#===============================

class Fig_2D_lon_lat(Fig_2D):

    #make_template is calling method from the parent class
    def make_template(self):
        super(Fig_2D_lon_lat, self).make_template('Plot 2D lon X lat','Ls 0-360','Level [Pa]','lon','lat')

    def do_plot(self):
        #create figure
        ax=super(Fig_2D_lon_lat, self).fig_init()
        try:    #try to do the figure, will return the error otherwise
            lon,lat,var,var_info=super(Fig_2D_lon_lat, self).data_loader_2D(self.varfull,self.plot_type)
            lon180,var=shift_data(lon,var)

            super(Fig_2D_lon_lat, self).filled_contour(lon180, lat,var)

            #---Add topo contour---
            zsurf=get_topo_2D(self.simuID,self.sol_array) #get topo

            lon180,zsurf=shift_data(lon,zsurf)
            plt.contour(lon180, lat,zsurf,10,colors='k',linewidths=0.5,linestyles='solid')   #topo
            #----

            if self.varfull2:
                _,_,var2,var_info2=super(Fig_2D_lon_lat, self).data_loader_2D(self.varfull2,self.plot_type)
                lon180,var2=shift_data(lon,var2)
                super(Fig_2D_lon_lat, self).solid_contour(lon180, lat,var2)
                var_info+=" (& "+var_info2+")"


            if self.Xlim:plt.xlim(self.Xlim)
            if self.Ylim:plt.ylim(self.Ylim)

            super(Fig_2D_lon_lat, self).make_title(var_info,'Longitude','Latitude')
             #--- Annotation---

            ax.xaxis.set_major_locator(MultipleLocator(30))
            ax.xaxis.set_minor_locator(MultipleLocator(10))
            ax.yaxis.set_major_locator(MultipleLocator(15))
            ax.yaxis.set_minor_locator(MultipleLocator(5))
            plt.xticks(fontsize=label_size-self.nPan//2, rotation=0)
            plt.yticks(fontsize=label_size-self.nPan//2, rotation=0)

            self.success=True

        except Exception as e: #Return the error
            super(Fig_2D_lon_lat, self).exception_handler(e,ax)
        super(Fig_2D_lon_lat, self).fig_save()

class Fig_2D_time_lat(Fig_2D):

    def make_template(self):
        #make_template is calling method from the parent class
        super(Fig_2D_time_lat, self).make_template('Plot 2D time X lat','Lon +/-180','Level [Pa]','sols','lat')
                                                                        #self.fdim1,  self.fdim2, self.Xlim,self.Ylim

    def do_plot(self):
        #create figure
        ax=super(Fig_2D_time_lat, self).fig_init()
        try:    #try to do the figure, will return the error otherwise

            t_stack,lat,var,var_info=super(Fig_2D_time_lat, self).data_loader_2D(self.varfull,self.plot_type)
            tim=t_stack[0,:];Ls=t_stack[1,:]

            super(Fig_2D_time_lat, self).filled_contour(Ls, lat,var)

            if self.varfull2:
                _,_,var2,var_info2=super(Fig_2D_time_lat, self).data_loader_2D(self.varfull2,self.plot_type)
                super(Fig_2D_time_lat, self).solid_contour(Ls, lat,var2)
                var_info+=" (& "+var_info2+")"


            #Axis formatting
            if self.Xlim:
                plt.xlim(self.Xlim)
            if self.Ylim:plt.ylim(self.Ylim)


            Ls_ticks = [item for item in ax.get_xticks()]
            labels = [item for item in ax.get_xticklabels()]


            for i in range(0,len(Ls_ticks)):
                id=np.argmin(np.abs(Ls-Ls_ticks[i])) #find tmstep closest to this tick
                labels[i]='Ls %g\nsol %i'%(np.mod(Ls_ticks[i],360.),tim[id])


            ax.set_xticklabels(labels)

            super(Fig_2D_time_lat, self).make_title(var_info,'','Latitude') #no 'Time' label as it is obvious

            ax.yaxis.set_major_locator(MultipleLocator(15))
            ax.yaxis.set_minor_locator(MultipleLocator(5))
            plt.xticks(fontsize=label_size-self.nPan//2, rotation=0)
            plt.yticks(fontsize=label_size-self.nPan//2, rotation=0)

            self.success=True

        except Exception as e: #Return the error
            super(Fig_2D_time_lat, self).exception_handler(e,ax)
        super(Fig_2D_time_lat, self).fig_save()

class Fig_2D_lat_press(Fig_2D):

    def make_template(self):
        #make_template is calling method from the parent class
        super(Fig_2D_lat_press, self).make_template('Plot 2D lat X press','Ls 0-360 ','Lon +/-180','Lat','level[Pa]')
                                                                          #self.fdim1,  self.fdim2, self.Xlim,self.Ylim
    def do_plot(self):
        #create figure
        ax=super(Fig_2D_lat_press, self).fig_init()
        try:    #try to do the figure, will return the error otherwise

            lat,pfull,var,var_info=super(Fig_2D_lat_press, self).data_loader_2D(self.varfull,self.plot_type)
            super(Fig_2D_lat_press, self).filled_contour(lat,pfull,var)

            if self.varfull2:
                _,_,var2,var_info2=super(Fig_2D_lat_press, self).data_loader_2D(self.varfull2,self.plot_type)
                super(Fig_2D_lat_press, self).solid_contour(lat, pfull,var2)
                var_info+=" (& "+var_info2+")"


            ax.set_yscale("log")
            ax.invert_yaxis()

            if self.Xlim:plt.xlim(self.Xlim)
            if self.Ylim:plt.ylim(self.Ylim)

            super(Fig_2D_lat_press, self).make_title(var_info,'Latitude','Pressure [Pa]')


            ax.xaxis.set_major_locator(MultipleLocator(15))
            ax.xaxis.set_minor_locator(MultipleLocator(5))
            plt.xticks(fontsize=label_size-self.nPan//2, rotation=0)
            plt.yticks(fontsize=label_size-self.nPan//2, rotation=0)

            def alt_KM(press,scale_height_KM=10.,reference_press=1000.):
                return -scale_height_KM*np.log(press/reference_press) # p to altitude in km
        # bound_P= np.sort([pfull[0],pfull[-1]]) #depending if *_plevs.nc or *.nc, the axis may be inverted
        #
        # def format_axis():
        #     ax2 = ax.twinx()
        #     plt.ylim([alt_KM(bound_P[-1],bound_P[-1]),alt_KM(bound_P[0],bound_P[-1])])
        #     plt.ylabel('Z [km]',fontsize=8)
        # format_axis()
            self.success=True
        except Exception as e: #Return the error
            super(Fig_2D_lat_press, self).exception_handler(e,ax)
        super(Fig_2D_lat_press, self).fig_save()

class Fig_2D_lon_press(Fig_2D):

    def make_template(self):
        #make_template is calling method from the parent class
        super(Fig_2D_lon_press, self).make_template('Plot 2D lon X press','Ls 0-360 ','Latitude','Lon +/-180','level[Pa]')

    def do_plot(self):
        #create figure
        ax=super(Fig_2D_lon_press, self).fig_init()
        try:    #try to do the figure, will return the error otherwise

            lon,pfull,var,var_info=super(Fig_2D_lon_press, self).data_loader_2D(self.varfull,self.plot_type)
            lon180,var=shift_data(lon,var)

            super(Fig_2D_lon_press, self).filled_contour(lon180,pfull,var)

            if self.varfull2:
                _,_,var2,var_info2=super(Fig_2D_lon_press, self).data_loader_2D(self.varfull2,self.plot_type)
                _,var2=shift_data(lon,var2)
                super(Fig_2D_lon_press, self).solid_contour(lon180, pfull,var2)
                var_info+=" (& "+var_info2+")"


            ax.set_yscale("log")
            ax.invert_yaxis()

            if self.Xlim:plt.xlim(self.Xlim)
            if self.Ylim:plt.ylim(self.Ylim)

            super(Fig_2D_lon_press, self).make_title(var_info,'Longitude','Pressure [Pa]')

            ax.xaxis.set_major_locator(MultipleLocator(30))
            ax.xaxis.set_minor_locator(MultipleLocator(10))
            plt.xticks(fontsize=label_size-self.nPan//2, rotation=0)
            plt.yticks(fontsize=label_size-self.nPan//2, rotation=0)

            self.success=True
        except Exception as e: #Return the error
            super(Fig_2D_lon_press, self).exception_handler(e,ax)
        super(Fig_2D_lon_press, self).fig_save()

class Fig_2D_time_press(Fig_2D):

    def make_template(self):
        #make_template is calling method from the parent class
        super(Fig_2D_time_press, self).make_template('Plot 2D time X press','Latitude','Lon +/-180','sols','level[Pa]')

    def do_plot(self):
        #create figure
        ax=super(Fig_2D_time_press, self).fig_init()
        try:    #try to do the figure, will return the error otherwise

            t_stack,pfull,var,var_info=super(Fig_2D_time_press, self).data_loader_2D(self.varfull,self.plot_type)
            tim=t_stack[0,:];Ls=t_stack[1,:]
            super(Fig_2D_time_press, self).filled_contour(Ls,pfull,var)

            if self.varfull2:
                _,_,var2,var_info2=super(Fig_2D_time_press, self).data_loader_2D(self.varfull2,self.plot_type)
                super(Fig_2D_time_press, self).solid_contour(Ls, pfull,var2)
                var_info+=" (& "+var_info2+")"


            #Axis formatting
            if self.Xlim:plt.xlim(self.Xlim)
            if self.Ylim:plt.ylim(self.Ylim)

            Ls_ticks = [item for item in ax.get_xticks()]
            labels = [item for item in ax.get_xticklabels()]


            for i in range(0,len(Ls_ticks)):
                id=np.argmin(np.abs(Ls-Ls_ticks[i])) #find tmstep closest to this tick
                labels[i]='Ls %g\nsol %i'%(np.mod(Ls_ticks[i],360.),tim[id])


            ax.set_xticklabels(labels)
            plt.xticks(fontsize=label_size-self.nPan//2, rotation=0)
            plt.yticks(fontsize=label_size-self.nPan//2, rotation=0)

            ax.set_yscale("log")
            ax.invert_yaxis()

            super(Fig_2D_time_press, self).make_title(var_info,'','Pressure [Pa]')

            self.success=True
        except Exception as e: #Return the error
            super(Fig_2D_time_press, self).exception_handler(e,ax)
        super(Fig_2D_time_press, self).fig_save()

class Fig_2D_lon_time(Fig_2D):

    def make_template(self):
        #make_template is calling method from the parent class
        super(Fig_2D_lon_time, self).make_template('Plot 2D lon X time','Latitude','Level [Pa]','Lon +/-180','sols')

    def do_plot(self):
        #create figure
        ax=super(Fig_2D_lon_time, self).fig_init()
        try:    #try to do the figure, will return the error otherwise

            lon,t_stack,var,var_info=super(Fig_2D_lon_time, self).data_loader_2D(self.varfull,self.plot_type)
            lon180,var=shift_data(lon,var)
            tim=t_stack[0,:];Ls=t_stack[1,:]
            super(Fig_2D_lon_time, self).filled_contour(lon180,Ls,var)

            if self.varfull2:
                _,_,var2,var_info2=super(Fig_2D_lon_time, self).data_loader_2D(self.varfull2,self.plot_type)
                _,var2=shift_data(lon,var2)
                super(Fig_2D_lon_time, self).solid_contour(lon180,Ls,var2)
                var_info+=" (& "+var_info2+")"


            #Axis formatting
            if self.Xlim:plt.xlim(self.Xlim)
            if self.Ylim:plt.ylim(self.Ylim)


            Ls_ticks = [item for item in ax.get_yticks()]
            labels = [item for item in ax.get_yticklabels()]


            for i in range(0,len(Ls_ticks)):
                id=np.argmin(np.abs(Ls-Ls_ticks[i])) #find tmstep closest to this tick
                labels[i]='Ls %g\nsol %i'%(np.mod(Ls_ticks[i],360.),tim[id])

            ax.set_yticklabels(labels)

            ax.xaxis.set_major_locator(MultipleLocator(30))
            ax.xaxis.set_minor_locator(MultipleLocator(10))

            super(Fig_2D_lon_time, self).make_title(var_info,'Longitude','')
            plt.xticks(fontsize=label_size-self.nPan//2, rotation=0)
            plt.yticks(fontsize=label_size-self.nPan//2, rotation=0)

            self.success=True
        except Exception as e: #Return the error
            super(Fig_2D_lon_time, self).exception_handler(e,ax)
        super(Fig_2D_lon_time, self).fig_save()

class Fig_1D(object):
    # Parent class for 1D figures
    def __init__(self,varfull='atmos_average.ts',doPlot=True):

        self.legend=None
        self.varfull=varfull
        self.t='AXIS' #default value for AXIS
        self.lat=None
        self.lon=None
        self.lev=None
        # Logic
        self.doPlot=doPlot
        self.plot_type='1D_time'

        #Extract filetype, variable, and simulation ID (initialization only)
        self.sol_array,self.filetype,self.var,self.simuID=split_varfull(self.varfull)

        #Multi pannel
        self.nPan=1
        self.subID=1
        self.addLine=False
        #Annotation for free dimensions
        self.fdim_txt=''
        self.success=False
        #Axis options

        self.Dlim=None #Dimension limit
        self.Vlim=None #variable limits
        self.axis_opts='-'



    def make_template(self):
        customFileIN.write("<<<<<<<<<<<<<<| Plot 1D = {0} |>>>>>>>>>>>>>\n".format(self.doPlot))
        customFileIN.write("Legend         = %s\n"%(self.legend))             #1
        customFileIN.write("Main Variable  = %s\n"%(self.varfull))            #2
        customFileIN.write("Ls 0-360       = {0}\n".format(self.t))           #3
        customFileIN.write("Latitude       = {0}\n".format(self.lat))         #4
        customFileIN.write("Lon +/-180     = {0}\n".format(self.lon))         #5
        customFileIN.write("Level [Pa]     = {0}\n".format(self.lev))         #6
        customFileIN.write("Axis Options  : lat,lon+/-180,[Pa],sols = [None,None] | var = [None,None] | linestyle = - \n")#7

    def read_template(self):
        self.legend= rT('char')             #1
        self.varfull=rT('char')             #2
        self.t=rT('float')                  #3
        self.lat=rT('float')                #4
        self.lon=rT('float')                #5
        self.lev=rT('float')                #6
        self.Dlim,self.Vlim,self.axis_opts=read_axis_options(customFileIN.readline())     #7

        self.plot_type=self.get_plot_type()


    def get_plot_type(self):
        '''
        Note that the or self.t =='AXIS' test  and the  self.t  =-88888 assignment are only used when MarsPlot is used without a template
        '''
        ncheck=0
        graph_type='Error'
        if self.t  ==-88888 or self.t    =='AXIS': self.t  =-88888;graph_type='1D_time';ncheck+=1
        if self.lat==-88888 or self.lat  =='AXIS': self.lat=-88888;graph_type='1D_lat' ;ncheck+=1
        if self.lon==-88888 or self.lon  =='AXIS': self.lon=-88888;graph_type='1D_lon' ;ncheck+=1
        if self.lev==-88888 or self.lev  =='AXIS': self.lev=-88888;graph_type='1D_lev' ;ncheck+=1
        if ncheck==0:
            prYellow('''*** Warning *** In 1D plot, %s: use 'AXIS' to set varying dimension '''%(self.varfull))
        if ncheck>1:
            prYellow('''*** Warning *** In 1D plot, %s: 'AXIS' keyword may only be used once '''%(self.varfull))
        return graph_type


    def prep_file(self,var_name,file_type,simuID,sol_array):
        global input_paths
        global Ncdf_num
        #first file in the list
        # A specific sol was requested, e.g [2400]
        if sol_array:
            Sol_num_current=sol_array
        else : #no sol requested, use default as provided by MarsPlot Custom.in -d sol
            Sol_num_current =Ncdf_num

        nfiles=len(Sol_num_current)


        file_list = [None]*nfiles #initialize the list

        #Loop over the requested time steps
        for i in range(0,nfiles):
            file_list[i] = input_paths[simuID]+'/%05d.'%(Sol_num_current[i])+file_type+'.nc'
            check_file_tape(file_list[i],abort=False)

        if file_type=='fixed': # XXXX.fixed.nc does not have an aggregation dimension so we use Dataset
            f=Dataset(file_list[0], 'r')
        else:
            f=MFDataset(file_list, 'r') #use MFDataset instead

        var_info=f.variables[var_name].long_name+' ['+ f.variables[var_name].units+']'
        dim_info=f.variables[var_name].dimensions
        dims=f.variables[var_name].shape
        return f, var_info,dim_info, dims

    def data_loader_1D(self,varfull,plot_type):


        if not '[' in varfull:
            if '{' in varfull :
                varfull,t_req,lat_req,lon_req,lev_req=get_overwrite_dim_1D(varfull,self.t,self.lat,self.lon,self.lev)
                # t_req,lat_req,lon_req,lev_req constain the dimensions to overwrite is '{}' are provided of the defaultself.t,self.lat,self.lon,self.lev otherwise
            else: # no '{ }' use to overwrite the dimensions, copy the plots' defaults
                t_req,lat_req,lon_req,lev_req= self.t,self.lat,self.lon,self.lev
            sol_array,filetype,var,simuID=split_varfull(varfull)
            xdata,var,var_info=self.read_NCDF_1D(var,filetype,simuID,sol_array,plot_type,t_req,lat_req,lon_req,lev_req)

        else:
            VAR=[]
            # Extract individual variables and prepare for execution
            varfull=remove_whitespace(varfull)
            varfull_list=get_list_varfull(varfull)
            expression_exec=create_exec(varfull,varfull_list)

            #Initialize list of requested dimensions;
            t_list=[None]*len(varfull_list)
            lat_list=[None]*len(varfull_list)
            lon_list=[None]*len(varfull_list)
            lev_list=[None]*len(varfull_list)
            expression_exec=create_exec(varfull,varfull_list)

            for i in range(0,len(varfull_list)):
                #---If overwriting dimensions, get the new dimensions and trim varfull from the '{lev=5.}' part
                if '{' in varfull_list[i] :
                    varfull_list[i],t_list[i],lat_list[i],lon_list[i],lev_list[i]=get_overwrite_dim_1D(varfull_list[i],self.t,self.lat,self.lon,self.lev)
                else: # no '{ }' use to overwrite the dimensions, copy the plots' defaults
                    t_list[i],lat_list[i],lon_list[i],lev_list[i]=self.t,self.lat,self.lon,self.lev
                sol_array,filetype,var,simuID=split_varfull(varfull_list[i])
                xdata,temp,var_info=self.read_NCDF_1D(var,filetype,simuID,sol_array,plot_type,t_list[i],lat_list[i],lon_list[i],lev_list[i])
                VAR.append(temp)
            var_info=varfull
            var=eval(expression_exec)

        return xdata,var,var_info

    def read_NCDF_1D(self,var_name,file_type,simuID,sol_array,plot_type,t_req,lat_req,lon_req,lev_req):
        '''
        Given an expression object with '[]' return the different variable needed
        Args:
            var_name: variable name, e.g 'temp'
            file_type: 'fixed' or 'atmos_average'
            sol_array: sol if different from default e.g '02400'
            plot_type: e.g '1D_lon', '1D_lat'
            t_req,lat_req,lon_req,lev_req: the Ls, lat, lon and level [Pa] requested
        Returns:
            dim_array: the axis, e.g one array of longitudes
            var_array: the variable extracted

        '''

        f, var_info,dim_info, dims=self.prep_file(var_name,file_type,simuID,sol_array)

        #If self.fdim is empty, add the variable (do only once)
        add_fdim=False
        if not self.fdim_txt.strip():add_fdim=True

        #Initialize dimensions (These are in all the .nc files)

        lat=f.variables['lat'][:];lati=np.arange(0,len(lat))
        lon=f.variables['lon'][:];loni=np.arange(0,len(lon))

        #Load variable depending on the requested free dimensions

        #======static======= , ignore level and time dimension
        if dim_info==(u'lat', u'lon'):
            if plot_type=='1D_lat':
                loni,temp_txt =get_lon_index(lon_req,lon)
            elif plot_type=='1D_lon':
                lati,temp_txt =get_lat_index(lat_req,lat)

            if add_fdim:self.fdim_txt+=temp_txt
            var=f.variables[var_name][lati,loni].reshape(len(np.atleast_1d(lati)),len(np.atleast_1d(loni)))
            f.close()

            if plot_type=='1D_lat': return  lat, np.mean(var,axis=1),var_info
            if plot_type=='1D_lon': return  lon, np.mean(var,axis=0),var_info

        #======time,lat,lon=======
        if f.variables[var_name].dimensions==(u'time', u'lat', u'lon'):
        #Initialize dimension
            t=f.variables['time'][:];Ls=np.squeeze(f.variables['areo'][:]);ti=np.arange(0,len(t))
            t_stack=np.vstack((t,Ls)) #stack the time and ls array as one variable

            if plot_type=='1D_lat':
                ti,temp_txt =get_time_index(t_req,Ls,t)
                if add_fdim:self.fdim_txt+=temp_txt
                loni,temp_txt =get_lon_index(lon_req,lon)
                if add_fdim:self.fdim_txt+=temp_txt
            if plot_type=='1D_lon':
                lati,temp_txt =get_lat_index(lat_req,lat)
                if add_fdim:self.fdim_txt+=temp_txt
                ti,temp_txt =get_time_index(t_req,Ls,t)
                if add_fdim:self.fdim_txt+=temp_txt
            if plot_type=='1D_time':
                loni,temp_txt =get_lon_index(lon_req,lon)
                if add_fdim:self.fdim_txt+=temp_txt
                lati,temp_txt =get_lat_index(lat_req,lat)
                if add_fdim:self.fdim_txt+=temp_txt

            #Extract data and close file
            var=f.variables[var_name][ti,lati,loni].reshape(len(np.atleast_1d(ti)),len(np.atleast_1d(lati)),len(np.atleast_1d(loni)))
            f.close()
            #Return data
            if plot_type=='1D_lat': return lat,    np.mean(np.mean(var,axis=2),axis=0),var_info
            if plot_type=='1D_lon': return lon,    np.mean(np.mean(var,axis=1),axis=0),var_info #transpose, Xdim must be in last column of var
            if plot_type=='1D_time':return t_stack,np.mean(np.mean(var,axis=2),axis=1),var_info


        #======time,level,lat,lon=======
        if (dim_info==(u'time', u'pfull', u'lat', u'lon')
           or dim_info==(u'time', u'pstd', u'lat', u'lon')
           or dim_info==(u'time', u'zgrid', u'lat', u'lon')):

            #Initialize dimensions
            if dim_info[1]=='pfull': levs=f.variables[dim_info[1]][:] 
            if dim_info[1]=='pstd': levs=f.variables[dim_info[1]][:] 
            if dim_info[1]=='zgrid': levs=     f.variables[dim_info[1]][:] # meters
            zi=np.arange(0,len(levs))
            t=f.variables['time'][:];Ls=np.squeeze(f.variables['areo'][:]);ti=np.arange(0,len(t))
            t_stack=np.vstack((t,Ls)) #stack the time and ls array as one variable

            if plot_type=='1D_lat':
                ti,temp_txt =get_time_index(t_req,Ls,t)
                if add_fdim:self.fdim_txt+=temp_txt
                loni,temp_txt =get_lon_index(lon_req,lon)
                if add_fdim:self.fdim_txt+=temp_txt
                zi,temp_txt =get_level_index(lev_req,levs)
                if add_fdim:self.fdim_txt+=temp_txt

            if plot_type=='1D_lon':
                lati,temp_txt =get_lat_index(lat_req,lat)
                if add_fdim:self.fdim_txt+=temp_txt
                ti,temp_txt =get_time_index(t_req,Ls,t)
                if add_fdim:self.fdim_txt+=temp_txt
                zi,temp_txt =get_level_index(lev_req,levs)
                if add_fdim:self.fdim_txt+=temp_txt

            if plot_type=='1D_time':
                loni,temp_txt =get_lon_index(lon_req,lon)
                if add_fdim:self.fdim_txt+=temp_txt
                lati,temp_txt =get_lat_index(lat_req,lat)
                if add_fdim:self.fdim_txt+=temp_txt
                zi,temp_txt =get_level_index(lev_req,levs)
                if add_fdim:self.fdim_txt+=temp_txt

            if plot_type=='1D_lev':
                ti,temp_txt =get_time_index(t_req,Ls,t)
                if add_fdim:self.fdim_txt+=temp_txt
                lati,temp_txt =get_lat_index(lat_req,lat)
                if add_fdim:self.fdim_txt+=temp_txt
                loni,temp_txt =get_lon_index(lon_req,lon)
                if add_fdim:self.fdim_txt+=temp_txt


            var=f.variables[var_name][ti,zi,lati,loni].reshape(len(np.atleast_1d(ti)),\
                                                               len(np.atleast_1d(zi)),\
                                                               len(np.atleast_1d(lati)),\
                                                               len(np.atleast_1d(loni)))
            f.close()
            #(u'time', u'pfull', u'lat', u'lon')
            if plot_type=='1D_lat': return lat,    np.mean(np.mean(np.mean(var,axis=3),axis=1),axis=0),var_info
            if plot_type=='1D_lon': return lon,    np.mean(np.mean(np.mean(var,axis=1),axis=1),axis=0),var_info
            if plot_type=='1D_time':return t_stack,np.mean(np.mean(np.mean(var,axis=1),axis=2),axis=1),var_info
            if plot_type=='1D_lev': return levs,   np.mean(np.mean(np.mean(var,axis=3),axis=2),axis=0),var_info



    def exception_handler(self,e,ax):
        if debug:raise
        sys.stdout.write("\033[F");sys.stdout.write("\033[K")
        prYellow('*** Warning *** Attempting %s profile for %s: %s'%(self.plot_type,self.varfull,str(e)))
        ax.text(0.5, 0.5, 'ERROR:'+str(e),horizontalalignment='center',verticalalignment='center', \
            bbox=dict(boxstyle="round",ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8),),\
            transform=ax.transAxes,wrap=True,fontsize=16)

    def fig_init(self):
        #create figure
        out=fig_layout(self.subID,self.nPan)
        if self.subID==1 and not self.addLine:
            fig= plt.figure(facecolor='white',figsize=(pixel_width/my_dpi, pixel_width/1.4/my_dpi)) #create figure if 1st pannel
            #plt.suptitle(simulation_name) #TODO remove 
        if not self.addLine:
            ax = plt.subplot(out[0],out[1],out[2]) #nrow,ncol,subID
        else:

            ax=plt.gca()

        return ax

    def fig_save(self):

        #save the figure
        if  self.subID==self.nPan : #Last subplot
            if  self.subID==1: #1 plot
                if not '[' in self.varfull:
                    sensitive_name=self.varfull.split('{')[0].strip()  #add split '{' in case varfull contains layer, does not do anything otherwise
                else:
                    sensitive_name='expression_'+get_list_varfull(self.varfull)[0].split('{')[0].strip()
            else: #multi pannel
                sensitive_name='multi_pannel'

            self.fig_name=output_path+'/plots/'+sensitive_name+'.'+out_format
            self.fig_name=create_name(self.fig_name)

            if i_list< len(objectList)-1 and not objectList[i_list+1].addLine:
                plt.savefig(self.fig_name,dpi=my_dpi )
                if out_format!="pdf":print("Saved:" +self.fig_name)
            #Last subplot
            if i_list== len(objectList)-1 :
                plt.savefig(self.fig_name,dpi=my_dpi )
                if out_format!="pdf":print("Saved:" +self.fig_name)


    def do_plot(self):
        #create figure
        ax=self.fig_init()


        try:    #try to do the figure, will return the error otherwise

            xdata,var,var_info=self.data_loader_1D(self.varfull,self.plot_type)

            if self.legend:
                txt_label=self.legend
            else:
                txt_label=var_info+'\n'+self.fdim_txt[1:]#we remove the first coma ',' of fdim_txt to print to the new line


            if self.plot_type=='1D_lat':

                plt.plot(var,xdata,self.axis_opts,lw=2,label=txt_label)
                plt.xlabel(var_info)
                plt.ylabel('Latitude')

                ax.yaxis.set_major_locator(MultipleLocator(15))
                ax.yaxis.set_minor_locator(MultipleLocator(5))
                if self.Dlim:plt.ylim(self.Dlim)
                if self.Vlim:plt.xlim(self.Vlim)

            if self.plot_type=='1D_lon':
                lon180,var=shift_data(xdata,var)

                plt.plot(lon180,var,self.axis_opts,lw=2,label=txt_label)
                plt.xlabel('Longitude')
                plt.ylabel(var_info)

                ax.xaxis.set_major_locator(MultipleLocator(30))
                ax.xaxis.set_minor_locator(MultipleLocator(10))
                if self.Dlim:plt.xlim(self.Dlim)
                if self.Vlim:plt.ylim(self.Vlim)


            if self.plot_type=='1D_time':
                tim=xdata[0,:];Ls=xdata[1,:]

                plt.plot(Ls,var,self.axis_opts,lw=2,label=txt_label)
                plt.ylabel(var_info)
                #Axis formatting
                if self.Vlim:plt.ylim(self.Vlim)

                if self.Dlim:
                    plt.xlim(self.Dlim) #TODO

                Ls_ticks = [item for item in ax.get_xticks()]
                labels = [item for item in ax.get_xticklabels()]

                for i in range(0,len(Ls_ticks)):
                    id=np.argmin(np.abs(Ls-Ls_ticks[i])) #find tmstep closest to this tick
                    labels[i]='Ls %g\nsol %i'%(np.mod(Ls_ticks[i],360.),tim[id])

                ax.set_xticklabels(labels)

            if self.plot_type=='1D_lev':

                plt.plot(var,xdata,self.axis_opts,lw=2,label=txt_label)
                plt.xlabel(var_info)
                plt.ylabel('Pressure [Pa]')

                ax.set_yscale("log")
                ax.invert_yaxis()

                if self.Dlim:plt.ylim(self.Dlim)
                if self.Vlim:plt.xlim(self.Vlim)


            #====comon labelling====
            plt.xticks(fontsize=label_size-self.nPan//2, rotation=0)
            plt.yticks(fontsize=label_size-self.nPan//2, rotation=0)
            plt.legend(fontsize=label_size-self.nPan//2)
            plt.grid(True)


            self.success=True
        except Exception as e: #Return the error
            self.exception_handler(e,ax)
        self.fig_save()

#======================================================
#                  END OF PROGRAM
#======================================================

if __name__ == '__main__':
    main()
