import os
import sys
import subprocess
from netCDF4 import Dataset
import numpy as np
#=========================================================================   
#=========================Scripts utilities===============================
#=========================================================================


# The functions below allow to print in different color
def prRed(skk): print("\033[91m{}\033[00m" .format(skk)) 
def prGreen(skk): print("\033[92m{}\033[00m" .format(skk)) 
def prCyan(skk): print("\033[96m{}\033[00m" .format(skk)) 
def prYellow(skk): print("\033[93m{}\033[00m" .format(skk)) 
def prPurple(skk): print("\033[95m{}\033[00m" .format(skk)) 
def prLightPurple(skk): print("\033[94m{}\033[00m" .format(skk)) 

def MY_func(Ls_cont):
    '''
    This function return the Mars Year
    Args:
        Ls_cont: solar longitude, contineuous
    Returns:
        MY : int the Mars year
    '''
    return (Ls_cont)//(360.)+1
    
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
        prYellow("To use this feature, create a hidden config file in your home directory with:") 
        prCyan("    cp amesGCM3/mars_templates/amesgcm_profile  ~/.amesgcm_profile")      
        
    except Exception as exception: #Return the error
        prRed('Error')
        print(exception)    
#=========================================================================   
#=========================================================================
#=========================================================================
