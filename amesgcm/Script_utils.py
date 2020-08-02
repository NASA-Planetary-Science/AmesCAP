import os
import sys
import subprocess
from netCDF4 import Dataset, MFDataset
import numpy as np
import re
#=========================================================================   
#=========================Scripts utilities===============================
#=========================================================================


# The functions below allow to print in different color
def prRed(skk): print("\033[91m{}\033[00m".format(skk)) 
def prGreen(skk): print("\033[92m{}\033[00m".format(skk)) 
def prCyan(skk): print("\033[96m{}\033[00m".format(skk)) 
def prYellow(skk): print("\033[93m{}\033[00m".format(skk)) 
def prPurple(skk): print("\033[95m{}\033[00m".format(skk)) 
def prLightPurple(skk): print("\033[94m{}\033[00m".format(skk)) 

def MY_func(Ls_cont):
    '''
    This function return the Mars Year
    Args:
        Ls_cont: solar longitude, contineuous
    Returns:
        MY : int the Mars year
    '''
    return (Ls_cont)//(360.)+1
 
def find_tod_in_diurn(fNcdf): 
    '''
    Return the variable for the local time axis in diurn files. 
    Original implementation by Victoria H.
    Args:
        fNcdf: an (open) Netcdf file object 
    Return:
        tod (string): 'time_of_day_16'or 'time_of_day_24'
    '''
    regex=re.compile('time_of_day.')
    varset=fNcdf.variables.keys()
    return [string for string in varset if re.match(regex, string)][0] #Exctract the 1st element of the list 
 


 
    
def print_fileContent(fileNcdf):
    '''
    Print the content of a Netcdf file in a compact format. Variables are sorted by dimensions.
    Args:
        fileNcdf: full path to netcdf file
    Returns: 
        None (print in the terminal)
    '''    
    #Define Colors for printing
    def Green(skk): return"\033[92m{}\033[00m".format(skk)
    def Cyan(skk): return "\033[96m{}\033[00m".format(skk)
    def Yellow(skk):return"\033[93m{}\033[00m".format(skk)
    def Purple(skk):return"\033[95m{}\033[00m".format(skk)
    if not os.path.isfile(fileNcdf):
        print(fileNcdf+' not found')
    else:    
        f=Dataset(fileNcdf, 'r')
        print("===================DIMENSIONS==========================")
        print(list(f.dimensions.keys()))
        print(str(f.dimensions))
        print("====================CONTENT==========================")
        all_var=f.variables.keys() #get all variables
        all_dims=list() #initialize empty list
        for ivar in all_var:
            all_dims.append(f.variables[ivar].dimensions) #get all the variables dimensions
        all_dims=set(all_dims) #filter duplicates (an object of type set() is an unordered collections of distinct objects 
        all_dims=sorted(all_dims,key=len) #sort dimensions by lenght, e.g ('lat') will come before ('lat','lon')
        var_done=list()
        for idim in all_dims:
            for ivar in all_var:
                if f.variables[ivar].dimensions==idim :
                    txt_dim=getattr(f.variables[ivar],'dimensions','') 
                    txt_shape=getattr(f.variables[ivar],'shape','')
                    txt_long_name=getattr(f.variables[ivar],'long_name','')
                    txt_units=getattr(f.variables[ivar],'units','')
                    print(Green(ivar.ljust(15))+': '+Purple(txt_dim)+'= '+Cyan(txt_shape)+', '+Yellow(txt_long_name)+\
                    '  ['+txt_units+']')

        try: #This part will be skipped if  the netcdf file does not contains a 'time' variable
            t_ini=f.variables['time'][0];t_end=f.variables['time'][-1]
            Ls_ini=np.squeeze(f.variables['areo'])[0];Ls_end=np.squeeze(f.variables['areo'])[-1]
            MY_ini=MY_func(Ls_ini);MY_end=MY_func(Ls_end)
            print('')
            print('Ls ranging from %6.2f to %6.2f: %.2f days'%(np.mod(Ls_ini,360.),np.mod(Ls_end,360.),t_end-t_ini))
            print('               (MY %02i)   (MY %02i)'%(MY_ini,MY_end))
        except:
            pass
        f.close()
        print("=====================================================") 

def print_varContent(fileNcdf,list_varfull,print_stat=False):
    '''
    Print the content of a variable inside a Netcdf file
    This test is based on the existence of a least one  00XXX.fixed.nc in the current directory.
    Args:
        fileNcdf:      full path to netcdf file
        list_varfull:  list of variable names and optional slices, e.g ['lon' ,'ps[:,10,20]']
        print_stat:  if true, print min, mean and max instead of values
    Returns: 
        None (print in the terminal)
    '''    
    #Define Colors for printing
    def Cyan(skk): return "\033[96m{}\033[00m".format(skk)
    def Red(skk):  return "\033[91m{}\033[00m".format(skk)
    if not os.path.isfile(fileNcdf):
        print(fileNcdf+' not found')
    else: 
    
        if print_stat:
            print(Cyan('__________________________________________________________________________'))
            print(Cyan('           VAR            |      MIN      |      MEAN     |      MAX      |'))
            print(Cyan('__________________________|_______________|_______________|_______________|'))
        for varfull in list_varfull:
            try:
                slice='[:]'
                if '[' in varfull:
                    varname,slice=varfull.strip().split('[');slice='['+slice
                else:
                    varname=varfull.strip()
                cmd_txt="""f.variables['"""+varname+"""']"""+slice
                f=Dataset(fileNcdf, 'r')
                var=eval(cmd_txt)
                
                if print_stat:
                    Min=np.nanmin(var)   
                    Mean=np.nanmean(var)
                    Max=np.nanmax(var)
                    print(Cyan('%26s|%15g|%15g|%15g|'%(varfull,Min,Mean,Max)))
                else:
                    print(Cyan(varfull+'= '))
                    print(Cyan(var))
                    print(Cyan('______________________________________________________________________')) 
            except: 
                if print_stat:
                    print(Red('%26s|%15s|%15s|%15s|'%(varfull,'','','')))   
                else:
                    print(Red(varfull))
        #Last line for the table
        if print_stat:
            print(Cyan('__________________________|_______________|_______________|_______________|'))           
        f.close()
        
         



def give_permission(filename): 
    '''
    # NAS system only: set group permission to the file
    '''
    import subprocess
    import os
    
    try:
        subprocess.check_call(['setfacl -v'],shell=True,stdout=open(os.devnull, "w"),stderr=open(os.devnull, "w")) #catch error and standard output
        cmd_txt='setfacl -R -m g:s0846:r '+filename
        subprocess.call(cmd_txt,shell=True)
    except subprocess.CalledProcessError: 
        pass

def check_file_tape(fileNcdf,abort=False):
    '''
    Relevant for use on the NASA Advanced Supercomputing (NAS) environnment only
    Check if a file is present on the disk by running the NAS dmls -l data migration command.
    This avoid the program to stall if the files need to be migrated from the disk to the tape
    Args:
        fileNcdf: full path to netcdf file
        exit: boolean. If True, exit the program (avoid stalling the program if file is not on disk)
    Returns:
        None (print status and abort program)
    '''
    # If the filename provided is not a netcdf file, exit program right away
    if fileNcdf[-3:]!='.nc':
        prRed('*** Error ***')
        prRed(fileNcdf + ' is not a netcdf file \n' )
        exit()
    #== Then check if the file actually exists on the system,  exit otherwise.


    try:
        #== NAS system only: file exists, check if it is active on disk or needs to be migrated from Lou
        subprocess.check_call(["dmls"],shell=True,stdout=open(os.devnull, "w"),stderr=open(os.devnull, "w")) #check if dmls command is available (NAS systems only)
        cmd_txt='dmls -l '+fileNcdf+"""| awk '{print $8,$9}'""" #get the last columns of the ls command with filename and status
        dmls_out=subprocess.check_output(cmd_txt,shell=True).decode('utf-8')  # get 3 letter identifier from dmls -l command, convert byte to string for Python 3
        if dmls_out[1:4] not in ['DUL','REG','MIG']: #file is OFFLINE, UNMIGRATING etc...
            if abort :
                prRed('*** Error ***')
                print(dmls_out)
                prRed(dmls_out[6:-1]+ ' is not available on disk, status is: '+dmls_out[0:5])
                prRed('CHECK file status with  dmls -l *.nc and run  dmget *.nc to migrate the files')
                prRed('Exiting now... \n')
                exit()
            else:
                prYellow('*** Warning ***')
                prYellow(dmls_out[6:-1]+ ' is not available on disk, status is: '+dmls_out[0:5])
                prYellow('Consider checking file status with  dmls -l *.nc and run  dmget *.nc to migrate the files')
                prYellow('Waiting for file to be migrated to disk, this may take a while...')
    except subprocess.CalledProcessError: #subprocess.check_call return an eror message
        if abort :
             exit()
        else:
            pass


def get_Ncdf_path(fNcdf):
    '''
    Return the full path of a Netcdf object.
    Note that 'Dataset' and multi-files dataset (i.e. 'MFDataset') have different
    attributes for the path, hence the need for this function.
    Args:
        fNcdf : Dataset or  MFDataset object          
    Returns :
        Path: string (list) for Dataset (MFDataset) 
    '''
    fname_out=getattr(fNcdf,'_files',False) #Only MFDataset has the_files attribute
    if not fname_out: fname_out=getattr(fNcdf,'filepath')() #Regular Dataset
    return fname_out


def FV3_file_type(fNcdf): 
    '''
    Return the type of output files:
    Args:
        fNcdf: an (open) Netcdf file object 
    Return:
       f_type (string): 'fixed', 'average', 'daily', or 'diurn'
       interp_type (string): 'pfull','pstd','zstd','zagl'
    '''
    #Get the full path from the file
    fullpath=get_Ncdf_path(fNcdf)
    #If MFDataset, get the 1st file in the list
    if type(fullpath)==list:fullpath=fullpath[0]
    
    #Get the filename without the path
    _,filename=os.path.split(fullpath)
    
    #Initialize
    f_type='other'
    interp_type='unknown'
    
    if 'fixed'         in filename:f_type='fixed'
    if 'atmos_average' in filename:f_type='average'
    if 'atmos_daily'   in filename:f_type='daily'
    if 'atmos_diurn'   in filename:f_type='diurn'
    
    dims=fNcdf.dimensions.keys()
    if 'pfull' in dims: interp_type='pfull'
    if 'pstd'  in dims: interp_type='pstd'
    if 'zstd'  in dims: interp_type='zstd'
    if 'zagl'  in dims: interp_type='zagl'
    return f_type,interp_type



def alt_FV3path(fullpaths,alt,test_exist=True):
    '''
    Internal function. given an interpolated daily, diurn or average file
    return the raw or fixed file. Accept string or list as input
    Args:
        fullpaths : e.g '/u/path/00010.atmos_average_pstd.nc' or LIST
            alt: alternative type 'raw' or 'fixed' 
            test_exist=True test file existence on disk          
    Returns :
        Alternative path to raw or fixed file, e.g.
                    '/u/path/00010.atmos_average.nc'
                    '/u/path/00010.fixed.nc' 
    '''
    out_list=[] 
    one_element=False
    #Convert as list for generality
    if type(fullpaths)==str:
        one_element=True
        fullpaths=[fullpaths]
        
    for fullpath in fullpaths:    
        path, filename = os.path.split(fullpath)
        DDDDD=filename.split('.')[0] #Get the date
        ext=filename[-8:] #Get the extension
        #This is an interpolated file
        if alt=='raw':
            if ext in ['_pstd.nc','_zstd.nc','_zagl.nc','plevs.nc']:
                if ext =='plevs.nc':
                    file_raw=filename[0:-9]+'.nc'
                else:   
                    file_raw=filename[0:-8]+'.nc'
            else:
                raise ValueError('In alt_FV3path(), FV3 file %s not recognized'%(filename))
            new_full_path=path+'/'+file_raw
        if alt=='fixed':
            new_full_path=path+'/'+DDDDD+'.fixed.nc'
        if test_exist and not (os.path.exists(new_full_path)):
            raise ValueError('In alt_FV3path() %s does not exist '%(new_full_path))
               
        out_list.append(new_full_path)  
    if one_element:out_list=out_list[0]  
    return  out_list
    
def smart_reader(fNcdf,var_list,suppress_warning=False):
    """
    Smarter alternative to using var=fNcdf.variables['var'][:] when handling PROCESSED files that also check 
    matching XXXXX.atmos_average.nc (or daily...) and XXXXX.fixed.nc files
    
    Args:
        fNcdf: Netcdf file object (i.e. already opened with Dataset or MFDataset)
        var_list: variable or list of variables, e.g 'areo' or ['pk','bk','areo']
        suppress_warning: Suppress debug statement, useful if variable is not expected to be found in the file anyway  
    Returns:
        out_list: variables content as singleton or values to unpack
        
    -------    
    Example: 
    
    from netCDF4 import Dataset
    
    fNcdf=Dataset('/u/akling/FV3/00668.atmos_average_pstd.nc','r') 
     
    ucomp= fNcdf.variables['ucomp'][:]   # << this is the regular way
    vcomp= smart_reader(fNcdf,'vcomp')   # << this is exacly equivalent  
    pk,bk,areo= smart_reader(fNcdf,['pk','bk','areo'])  # this will get 'areo' from 00668.atmos_average.nc is not available in the original _pstd.nc file
                                                        # if pk and bk are absent from 0668.atmos_average.nc, it will also check 00668.fixed.n
    *** NOTE ***
        -Only the variables' content is returned, not the attributes
    """  
    
    #This out_list is for the variable
    out_list=[]
    one_element=False
    file_is_MF=False
    
    Ncdf_path= get_Ncdf_path(fNcdf) #Return string (Dataset) or list (MFDataset)
    if type(Ncdf_path)==list:file_is_MF=True
    
    #For generality convert to list if only one variable is provided, e.g 'areo'>['areo']
    if type(var_list)==str:
        one_element=True
        var_list=[var_list]

    for ivar in var_list:
    #First try to read in the original file
        if ivar in fNcdf.variables.keys():
            out_list.append(fNcdf.variables[ivar][:])
        else:
            full_path_try=alt_FV3path(Ncdf_path,alt='raw',test_exist=True) 
            if file_is_MF:
                f_tmp=MFDataset(full_path_try,'r')
            else:    
                f_tmp=Dataset(full_path_try,'r')
                
            if ivar in f_tmp.variables.keys():
                out_list.append(f_tmp.variables[ivar][:])
                if not suppress_warning: print('**Warning*** Using variable %s in %s instead of original file(s)'%(ivar,full_path_try))
                f_tmp.close()
            else:    
                f_tmp.close()
                full_path_try=alt_FV3path(Ncdf_path,alt='fixed',test_exist=True) 
                if file_is_MF:full_path_try=full_path_try[0]
                
                f_tmp=Dataset(full_path_try,'r')
                if ivar in f_tmp.variables.keys():
                    out_list.append(f_tmp.variables[ivar][:])
                    f_tmp.close()
                    if not suppress_warning: print('**Warning*** Using variable %s in %s instead of original file(s)'%(ivar,full_path_try))
                else: 
                    print('***ERROR*** Variable %s not found in %s, NOR in raw output or fixed file'%(ivar,full_path_try))
                    print('            >>> Assigning  %s  to NaN'%(ivar))
                    f_tmp.close()
                    out_list.append(np.NaN)
    if one_element:out_list=out_list[0]
    return out_list

            
def progress(k,Nmax):
    """
    Display a progress bar to monitor heavy calculations. 
    Args:
        k: current iteration of the outer loop
        Nmax: max iteration of the outer loop
    Returns:
        Running... [#---------] 10.64 %
    """  
    import sys
    from math import ceil #round yo the 2nd digit
    progress=float(k)/Nmax
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rRunning... [{0}] {1} {2}%".format( "#"*block + "-"*(barLength-block), ceil(progress*100*100)/100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


def section_content_amesgcm_profile(section_ID):
    '''
    Execude code section in /home/user/.amesgcm_profile
    Args:
        section_ID: string defining the section to loa, e.g 'Pressure definitions for pstd'
    Returns
        return line in that section as python code
    '''
    import os
    import numpy as np
    input_file=os.environ['HOME']+'/.amesgcm_profile'
    try:
        f=open(input_file, "r")
        contents='' 
        rec=False
        while True:
            line = f.readline()
            if not line :break #End of File
            if line[0]== '<':
                rec=False
                if line.split('|')[1].strip() ==section_ID:
                    rec=True
                    line = f.readline()
            if rec : contents+=line                
        f.close()     
        if contents=='':
            prRed("No content found for <<< %s >>> block"%(section_ID))     
        return contents

    except FileNotFoundError:
        prRed("Error: %s config file not found "%(input_file)) 
        prYellow("To use this feature, create a hidden config file from the template in your home directory with:") 
        prCyan("    cp amesGCM3/mars_templates/amesgcm_profile  ~/.amesgcm_profile")      
        exit()
    except Exception as exception: #Return the error
        prRed('Error')
        print(exception)  
        
def wbr_cmap():
    '''
    Returns a color map that goes from white>blue>green>yellow>red
    '''
    from matplotlib.colors import ListedColormap
    tmp_cmap = np.zeros((254,4))
    tmp_cmap [:,3]=1. #set alpha

    tmp_cmap[:,0:3]=np.array([[255,255,255],
    [252,254,255], [250,253,255], [247,252,254], [244,251,254], [242,250,254], 
    [239,249,254], [236,248,253], [234,247,253], [231,246,253], [229,245,253], 
    [226,244,253], [223,243,252], [221,242,252], [218,241,252], [215,240,252], 
    [213,239,252], [210,238,251], [207,237,251], [205,236,251], [202,235,251], 
    [199,234,250], [197,233,250], [194,232,250], [191,231,250], [189,230,250], 
    [186,229,249], [183,228,249], [181,227,249], [178,226,249], [176,225,249], 
    [173,224,248], [170,223,248], [168,222,248], [165,221,248], [162,220,247], 
    [157,218,247], [155,216,246], [152,214,245], [150,212,243], [148,210,242], 
    [146,208,241], [143,206,240], [141,204,238], [139,202,237], [136,200,236], 
    [134,197,235], [132,195,234], [129,193,232], [127,191,231], [125,189,230], 
    [123,187,229], [120,185,228], [118,183,226], [116,181,225], [113,179,224], 
    [111,177,223], [109,175,221], [106,173,220], [104,171,219], [102,169,218], 
    [100,167,217], [ 97,165,215], [ 95,163,214], [ 93,160,213], [ 90,158,212], 
    [ 88,156,211], [ 86,154,209], [ 83,152,208], [ 81,150,207], [ 79,148,206], 
    [ 77,146,204], [ 72,142,202], [ 72,143,198], [ 72,144,195], [ 72,145,191], 
    [ 72,146,188], [ 72,147,184], [ 72,148,181], [ 72,149,177], [ 72,150,173], 
    [ 72,151,170], [ 72,153,166], [ 72,154,163], [ 72,155,159], [ 72,156,156], 
    [ 72,157,152], [ 72,158,148], [ 72,159,145], [ 72,160,141], [ 72,161,138], 
    [ 73,162,134], [ 73,163,131], [ 73,164,127], [ 73,165,124], [ 73,166,120], 
    [ 73,167,116], [ 73,168,113], [ 73,169,109], [ 73,170,106], [ 73,172,102], 
    [ 73,173, 99], [ 73,174, 95], [ 73,175, 91], [ 73,176, 88], [ 73,177, 84], 
    [ 73,178, 81], [ 73,179, 77], [ 73,181, 70], [ 78,182, 71], [ 83,184, 71], 
    [ 87,185, 72], [ 92,187, 72], [ 97,188, 73], [102,189, 74], [106,191, 74], 
    [111,192, 75], [116,193, 75], [121,195, 76], [126,196, 77], [130,198, 77], 
    [135,199, 78], [140,200, 78], [145,202, 79], [150,203, 80], [154,204, 80], 
    [159,206, 81], [164,207, 81], [169,209, 82], [173,210, 82], [178,211, 83], 
    [183,213, 84], [188,214, 84], [193,215, 85], [197,217, 85], [202,218, 86], 
    [207,220, 87], [212,221, 87], [217,222, 88], [221,224, 88], [226,225, 89], 
    [231,226, 90], [236,228, 90], [240,229, 91], [245,231, 91], [250,232, 92], 
    [250,229, 91], [250,225, 89], [250,222, 88], [249,218, 86], [249,215, 85], 
    [249,212, 84], [249,208, 82], [249,205, 81], [249,201, 80], [249,198, 78], 
    [249,195, 77], [248,191, 75], [248,188, 74], [248,184, 73], [248,181, 71], 
    [248,178, 70], [248,174, 69], [248,171, 67], [247,167, 66], [247,164, 64], 
    [247,160, 63], [247,157, 62], [247,154, 60], [247,150, 59], [247,147, 58], 
    [246,143, 56], [246,140, 55], [246,137, 53], [246,133, 52], [246,130, 51], 
    [246,126, 49], [246,123, 48], [246,120, 47], [245,116, 45], [245,113, 44], 
    [245,106, 41], [244,104, 41], [243,102, 41], [242,100, 41], [241, 98, 41], 
    [240, 96, 41], [239, 94, 41], [239, 92, 41], [238, 90, 41], [237, 88, 41], 
    [236, 86, 41], [235, 84, 41], [234, 82, 41], [233, 80, 41], [232, 78, 41], 
    [231, 76, 41], [230, 74, 41], [229, 72, 41], [228, 70, 41], [228, 67, 40], 
    [227, 65, 40], [226, 63, 40], [225, 61, 40], [224, 59, 40], [223, 57, 40], 
    [222, 55, 40], [221, 53, 40], [220, 51, 40], [219, 49, 40], [218, 47, 40], 
    [217, 45, 40], [217, 43, 40], [216, 41, 40], [215, 39, 40], [214, 37, 40], 
    [213, 35, 40], [211, 31, 40], [209, 31, 40], [207, 30, 39], [206, 30, 39], 
    [204, 30, 38], [202, 30, 38], [200, 29, 38], [199, 29, 37], [197, 29, 37], 
    [195, 29, 36], [193, 28, 36], [192, 28, 36], [190, 28, 35], [188, 27, 35], 
    [186, 27, 34], [185, 27, 34], [183, 27, 34], [181, 26, 33], [179, 26, 33], 
    [178, 26, 32], [176, 26, 32], [174, 25, 31], [172, 25, 31], [171, 25, 31], 
    [169, 25, 30], [167, 24, 30], [165, 24, 29], [164, 24, 29], [162, 23, 29], 
    [160, 23, 28], [158, 23, 28], [157, 23, 27], [155, 22, 27], [153, 22, 27], 
    [151, 22, 26], [150, 22, 26], [146, 21, 25]])/255.
    
    return ListedColormap(tmp_cmap)  

def rjw_cmap():
    '''
    Returns a color map that goes from red<jade<wisteria or 'rjw'
    '''
    from matplotlib.colors import ListedColormap
    tmp_cmap = np.zeros((55,4))
    tmp_cmap [:,3]=1. #set alpha

    tmp_cmap[:,0:3]=np.array([[255,   0, 244],
       [248,  40, 244],[241,  79, 244],[234, 119, 244],[228, 158, 244],
       [221, 197, 245],[214, 190, 245],[208, 182, 245],[201, 175, 245],
       [194, 167, 245],[188, 160, 245],[181, 152, 246],[175, 145, 246],
       [140, 140, 247],[105, 134, 249],[ 70, 129, 251],[ 35, 124, 253],
       [  0, 119, 255],[  0, 146, 250],[  0, 173, 245],[  0, 200, 241],
       [  0, 227, 236],[  0, 255, 231],[  0, 255, 185],[  0, 255, 139],
       [  0, 255,  92],[  0, 255,  46],[  0, 255,   0],[ 63, 247,  43],
       [127, 240,  87],[191, 232, 130],[255, 225, 174],[255, 231, 139],
       [255, 237, 104],[255, 243,  69],[255, 249,  34],[255, 255,   0],
       [255, 241,  11],[255, 227,  23],[255, 213,  35],[255, 199,  47],
       [255, 186,  59],[255, 172,  71],[255, 160,  69],[255, 148,  67],
       [255, 136,  64],[255, 124,  62],[255, 112,  60],[255, 100,  58],
       [255,  80,  46],[255,  60,  34],[255,  40,  23],[255,  20,  11],
       [255,   0,   0],[237,  17,   0]])/255.
    return ListedColormap(tmp_cmap)    
                 
          
#=========================================================================   
#=========================================================================
#=========================================================================
