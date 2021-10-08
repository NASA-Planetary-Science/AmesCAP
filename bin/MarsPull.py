#!/usr/bin/env python

# The functions below allow to print in different color
def prCyan(skk):        print("\033[96m{}\033[00m".format(skk))
def prYellow(skk):      print("\033[93m{}\033[00m".format(skk))

import sys
import os

#=====Attempt to import specific scientic modules one may not find in the default python on NAS ====
try:
    import numpy as np
    import argparse #parse arguments
    import requests

except ImportError as error_msg:
    prYellow("Error while importing modules")
    prYellow('Your are using python '+str(sys.version_info[0:3])+', recommend version is (2, 7, 12)' )
    prYellow('Please load recommended version with:')
    prCyan('        module load python/2.7.12\n')
    print("Error was: "+ error_msg.message)
    exit()

except Exception as exception:
    # Output unexpected Exceptions.
    print(exception, False)
    print(exception.__class__.__name__ + ": " + exception.message)
    exit()



parser = argparse.ArgumentParser(description="""\033[93mUilities for accessing files on the MCMC NAS repository\033[00m """,       formatter_class=argparse.RawTextHelpFormatter)


parser.add_argument('-id','--id',type=str,
                             help="Simulation identifier corresponding to a sub-directory of /legacygcmdata on the MCMC Data portal: \n"
                                "\033[96mhttps://data.nas.nasa.gov/legacygcm/data_legacygcm.php?dir=/legacygcmdata\033[00m \n"
                                "\n"
                                  "Current options are : '\033[93mACTIVECLDS\033[00m', '\033[93mINERTCLDS\033[00m' and  '\033[93mACTIVECLDS_NCDF\033[00m' "
                                  "Usage: ./MarsPull.py -id  INERTCLDS \n")

parser.add_argument('-ls','--ls', nargs='+',type=float,
                             help='Query data by Solar Longitude \n'
                                  'Usage: ./MarsPull.py -ls  90. \n'
                                  '       ./MarsPull.py -ls  start   stop \n')

parser.add_argument('-f', '--filename',nargs='+',type=str,
                 help=' Query specific file(s) \n'
                      '> Usage:  ./MarsPull -id  ACTIVECLDS_NCDF -f LegacyGCM_fixed.nc  LegacyGCM_Ls000_Ls004.nc')
#===Define variables===


saveDir=os.getcwd()+'/'
#============================ Files available in the database ======================
Ls_ini=np.array([0,5,10,15,19,24,29,34,38,43,48,52,57,61,66,70,75,79,84,88,93,97,102,\
106,111,116,121,125,130,135,140,146,151,156,162,167,173,179,184,190,196,202,\
209,215,221,228,234,241,247,254,260,266,273,279,286,292,298,304,310,316,322,\
328,333,339,344,350,355])

Ls_end=np.array([4,9,14,19,24,29,33,38,42,47,52,56,61,65,70,74,79,83,88,92,97,\
101,106,111,115,120,125,130,135,140,145,150,156,161,167,172,178,184,190,196,\
202,208,214,221,227,233,240,246,253,259,266,272,279,285,291,297,304,310,316,\
321,327,333,338,344,349,354,0])
#====================================================================================
    
def download(url, filename):
    '''
    Save a file from a URL locally
    Args:
        URL: The URL to download, e.g   'https://data.nas.nasa.gov/legacygcm/download_data.php?file=/legacygcmdata/LegacyGCM_Ls000_Ls004.nc'
        filename: The local filename e.g  '/lou/la4/akling/Data/LegacyGCM_Ls000_Ls004.nc'
    '''
    
    _ , fname=os.path.split(filename)
    response = requests.get(url, stream=True)
    total = response.headers.get('content-length')

    if response.status_code == 404:
        print('Error during download, error code is: ',response.status_code)
    else:
        
        #If we have access to the size of the file, return progress bar
        if total is not None:
            with open(filename, 'wb') as f:
                downloaded = 0
                if total :total = int(total)
                for data in response.iter_content(chunk_size=max(int(total/1000), 1024*1024)):
                    downloaded += len(data)
                    f.write(data)
                    done = int(50*downloaded/total)
                    sys.stdout.write('\r[{}{}]'.format('#' * done, '.' * (50-done)))
                    sys.stdout.flush()
            sys.stdout.write('\n') 
        else:
            
            #Header is unknown yet, skip the use of a progress bar
            print('Downloading %s ...'%(fname))
            with open(local_file, 'wb')as f:
                f.write(data.content)
            print('%s Done'%(fname))           
                

#======================================================
#                  MAIN PROGRAM
#======================================================
def main():
    
    #Original
    #URLbase="https://data.nas.nasa.gov/legacygcm/download_data_legacygcm.php?file=/legacygcmdata/"
    simu_ID=parser.parse_args().id
    
    if simu_ID is None :
        prYellow("***Error*** simulation ID [-id --id] is required. See 'MarsPull.py -h' for help")
        exit()
    
    URLbase='https://data.nas.nasa.gov/legacygcm/download_data_legacygcm.php?file=/legacygcmdata/'+simu_ID+'/'
    
    
    if parser.parse_args().ls :
        data_input=np.asarray(parser.parse_args().ls)
        if len(data_input)==1: #query only  the file that contains this Ls
            i_start=np.argmin(np.abs(Ls_ini-data_input))
            if data_input<Ls_ini[i_start]:i_start-=1
            i_request=np.arange(i_start,i_start+1)
        
        elif len(data_input)==2: #start stop  is provided
            i_start=np.argmin(np.abs(Ls_ini-data_input[0]))
            if data_input[0]<Ls_ini[i_start]:i_start-=1

            i_end=np.argmin(np.abs(Ls_end-data_input[1]))
            if data_input[1]>Ls_end[i_end]:i_end+=1
         
            i_request=np.arange(i_start,i_end+1) 

        print('Saving  %i files in %s '%(len(i_request),saveDir))
        for ii in i_request:
            #Legacy .nc files  
            if simu_ID=='ACTIVECLDS_NCDF':
                fName='LegacyGCM_Ls%03d_Ls%03d.nc'%(Ls_ini[ii],Ls_end[ii])
            #fort.11 files    
            else:
                fName='fort.11_%04d'%(670+ii)
                
            url = URLbase+fName
            filename=saveDir+fName
            #print('Downloading '+ fName+ '...')
            print('Downloading '+ url+ '...')
            download(url,filename)

    elif parser.parse_args().filename:
        f_input=np.asarray(parser.parse_args().filename)
        for ff in f_input :
            url = URLbase+ff
            filename=saveDir+ff
            print('Downloading '+ url+ '...')#ff
            download(url,filename)
    else:
        prYellow('No data requested, use -ls or -f options')
#======================================================
#                  END OF PROGRAM
#======================================================

if __name__ == '__main__':
    main()

import wget
def bar_custom(current, total, width=80):
    print("Downloading: %d%% [%d / %d] bytes" % (current / total * 100, current, total))
wget.download('https://data.nas.nasa.gov/legacygcm/download_data_legacygcm.php?file=/legacygcmdata/INERTCLDS/fort.11_0674', bar=bar_custom)