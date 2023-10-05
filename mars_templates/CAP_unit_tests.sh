#!/bin/bash

# ===============================================================================================
# Updated April 6, 2023
#
# PURPOSE:            *** FOR INTERNAL MCMC USE ONLY ***
#            This script tests every function in CAP. It is meant for diagnosing problems
#            with development versions of CAP and is not intended for use outside of the MCMC.
#
# NOTES:     This script requires:
#               1. The folder test_files/ and its contents be in your current working directory
#               2. The *.in and *.nc files in test_files/ be copied to your current working
#                  directory
#               3. CAP_Test_All_Plots.in be copied to your amesCAP/mars_templates directory
#
#            Compatible with all shells: bash, csh, tcsh, or zsh.
#            Windows Cygwin users should be able to run this too.
#
# ===============================================================================================

# Colorized print statements
r='\033[1;31m'  # bold red
nc='\033[1;0m'  # bold no color
bb='\033[1;34m' # bold blue
g='\033[1;32m'  # bold green
y='\033[1;33m'  # bold yellow

# Seperator
sep='========================================================================================='
frame='***********************'
cmark="\xE2\x9C\x94"
redX="\xE2\x9D\x8C"

# ===============================================================================================
#                                          Intro
# ===============================================================================================
echo -e "\n\n ${bb}Created by Courtney Batterson | courtney.m.batterson@nasa.gov"
echo -e "NASA Ames Mars Climate Modeling Center"
echo -e "Last updated April 6, 2023"

# ===============================================================================================
#                          Set path (PTH) to netCDF output files 
# ===============================================================================================

echo -e "\n\n ${bb}${sep}${nc}"
echo -e "\n ${bb}>>> Running ${y}CAP_Diagnostics.sh${bb} in the directory: ${g}${PWD}${bb}"
echo -e "    If this is ${r}NOT${bb} where your output files are, provide the ${g}<path>${bb}"
echo -e "    to the correct directory. Otherwise, type ${g}Enter${bb} to continue:"

read PTH
if [ -z "$PTH" ]; then
    PTH=${PWD}
    echo -e "${g}>>> Confirm path = ${PTH}"
else
    echo -e "${g}>>> Confirm path = ${PTH}"
fi
cd $PTH

echo -e "\n ${r}      =============================================================${nc}"
echo -e "       ${r}PLEASE NOTE:"
echo -e "           You can stop this program anytime by typing ${y}Control-C${r}"
echo -e "           The program should complete in ${y}~15 minutes${nc}"
echo -e "${r}       =============================================================${nc}"

echo -e "\n\n ${bb}${sep}${nc}"
echo -e "\n ${bb}>>> Please ${g}confirm${bb} you are compliant with the following requirements:\n"
echo -e "    1. The folder ${y}test_files/${bb} and its contents are in your ${y}current directory${bb}"
echo -e "    2. You have copied the ${y}*.nc${bb} files in ${y}test_files/${bb} to your ${y}current directory${bb}"
echo -e "    3. If you have CAP installed in your ${y}HOME${bb} directory, then ${y}test_files/CAP_Test_All_Plots.in${bb} "
echo -e "       will be copied to ${y}~amesCAP/mars_templates${bb} automatically."
echo -e "       If you have CAP installed elsewhere, copy ${y}CAP_Test_All_Plots.in${bb} to "
echo -e "       ${y}mars_templates/${bb} manually and comment out line 451 in ${y}CAP_Diagnostics.sh${nc}"
echo -e "\n ${bb}>>> If these requirements are satisfied, type ${g}Enter${bb}"

read CONFIRMATION
if [ -z "$CONFIRMATION" ]; then
    echo -e "${g}>>> Confirmed. Let the program begin!"
else
    exit 1
fi

# ===============================================================================================
#                                    File Manipulations 
# ===============================================================================================

# =============================== Function file_exists =========================================
file_exists () {
    # USE: file_exists "filename"
    # $1=FILE 
    if ls ${PTH}/$1 1> /dev/null 2>&1; then
        echo -e "${g} ${cmark}  ${bb}$PTH/${g}$1${nc}"
        MarsPlot.py -i $1 > inspect.out
        true
    else
        if $DERIVED_FILE; then
            # File needs to be created from CAP functions
            echo -e "${y} Creating $1${nc}"
            false
        else
            # File doesn't exist in the folder at all, cannot be derived from CAP functions
            echo -e "${r} ${redX} $PTH/$1 MISSING${nc}"
            false
        fi
    fi
}

# ================================= Function confirm_var ======================================
confirm_var () {
    # USE: confirm_var "varname" "filename"
    # $1=VNAME, $2=FILE 
    #echo -e "${g} ${cmark}  ${y}$1${g} is in ${y}$2${nc}"
    echo -e "${bb}$2 ${y}$1${g} ${cmark}"
}

# ============================= Function report_missing_var ==================================
report_missing_var () {
    # USE: confirm_var "varname" "filename"
    # $1=VNAME, $2=FILE 
    #echo -e "${r} ${redX}  MISSING ${y}$1${r} in ${y}$2${nc}"
    echo -e "${r}$2 $1 ${redX}"
}

# ============================== Function check_for_var =====================================
check_for_var () {
    # USE: check_for_var "varname" "filename"
    # $1=VNAME, $2=FILE

    # Dump variables to inspect.out
    MarsPlot.py -i $2 > inspect.out
    if grep -q "$1" inspect.out; then
        confirm_var "$1" "$2"
        true
    else
        report_missing_var "$1" "$2"
        false
    fi
}

# ===============================================================================================
#                                           Main Program 
# ===============================================================================================
echo -e "\n\n ${bb}${sep}${nc}"
echo -e "${bb}      ${frame} CHECKING FOR FILES IN DIRECTORY ${frame} \n"

DAILY="01336.atmos_daily.nc"
DIURN="01336.atmos_diurn.nc"
AVERAGE="01336.atmos_average.nc"
PAVERAGE="01336.atmos_average_pstd.nc" # created in this script
AVERAGE2="00668.atmos_average.nc"
PDIURN="01336.atmos_diurn_pstd.nc"
PAVERAGEc48="00669.atmos_average_pstd_c48.nc"
FIXED="01336.fixed.nc"
COMBINE1="01336.atmos_average_for_combine1.nc"
COMBINE2="01336.atmos_average_for_combine2.nc"

# Check that all required FILES are present in the current directory
DERIVED_FILE=false
if file_exists $DAILY; then DAILY_STATUS=true; else DAILY_STATUS=false; fi
if file_exists $DIURN; then DIURN_STATUS=true; else DIURN_STATUS=false; fi
if file_exists $AVERAGE; then AVERAGE_STATUS=true; else AVERAGE_STATUS=false; fi
if file_exists $AVERAGE2; then AVERAGE2_STATUS=true; else AVERAGE2_STATUS=false; fi
if file_exists $PDIURN; then PDIURN_STATUS=true; else PDIURN_STATUS=false; fi
if file_exists $PAVERAGEc48; then PAVERAGEc48_STATUS=true; else PAVERAGEc48_STATUS=false; fi
if file_exists $FIXED; then FIXED_STATUS=true; else FIXED_STATUS=false; fi
if file_exists $COMBINE1; then COMBINE1_STATUS=true; else COMBINE1_STATUS=false; fi
if file_exists $COMBINE2; then COMBINE2_STATUS=true; else COMBINE2_STATUS=false; fi

echo -e "\n\n ${bb}${sep}${nc}"
echo -e "${bb}      ${frame} CHECKING FOR VARIABLES IN FILES ${frame} \n"

# Check that all required VARIABLES are present in above FILES
AVG_VARS=(pk bk temp ps ts dst_mass_micro vcomp ucomp r cldcol taudust_IR areo cld)
for VAR in ${AVG_VARS[@]}; do
    if $AVERAGE_STATUS; then check_for_var $VAR $AVERAGE; fi
done
AVG_VARS=(temp taudust_IR areo)
for VAR in ${AVG_VARS[@]}; do
    if $AVERAGE2_STATUS; then check_for_var $VAR $AVERAGE2; fi
done

DAY_VARS=(temp ps areo)
for VAR in ${DAY_VARS[@]}; do
    if $DAILY_STATUS; then check_for_var $VAR $DAILY; fi
    if $DIURN_STATUS; then check_for_var $VAR $DIURN; fi
done

DAY_VARS=(ps areo)
for VAR in ${DAY_VARS[@]}; do
    if $PDIURN_STATUS; then check_for_var $VAR $PDIURN; fi
done

FIX_VARS=(zsurf)
for VAR in ${FIX_VARS[@]}; do
    if $FIXED_STATUS; then check_for_var $VAR $FIXED; fi
done

# From now on, any call to file_exists refers to a file that CAN be created with CAP
DERIVED_FILE=true
# ================================ Interpolation ====================================
echo -e "\n\n ${bb}${sep}${nc}"
echo -e "${bb}      ${frame} TESTING MARSINTERP FUNCTIONS ${frame} \n"

if $AVERAGE_STATUS; then
    TYPES=(pstd zstd zagl)
    for TYPE in ${TYPES[@]}; do
        if file_exists ${AVERAGE%.*}_${TYPE}.nc; then
            echo -e "\n ${bb}>>> ${y}${AVERAGE%.*}_${TYPE}.nc${bb} exists already, skipping the pressure interpolation."
        else
            echo -e "\n ${bb}>>> Pressure-interpolating ${y}$AVERAGE${bb} to ${y}${TYPE}${bb}. \n    Only interpolating ${y}temp${bb} and ${y}r${bb} variables to save time.${nc} \n"
            MarsInterp.py $AVERAGE -t ${TYPE} --include temp r
            echo -e "\n ${g} ${cmark}  Interpolation COMPLETE: ${y}$AVERAGE${bb} to ${y}${TYPE}${nc} \n"
        fi
    done

    echo -e "\n\n ${bb}${sep}${nc}"

    if file_exists ${AVERAGE%.*}_pstd_phalf_mb.nc; then
        echo -e "\n ${bb}>>> ${y}${AVERAGE%.*}_pstd_phalf_mb.nc${bb} exists already, skipping the pressure interpolation."
    else
        echo -e "\n ${bb}>>> Pressure-interpolating ${y}$AVERAGE${bb} to ${y}pstd (phalf_mb)${bb} including only ${y}temp${bb} and creating a file called ${y}${AVERAGE%.*}_pstd_phalf_mb.nc${bb}...${nc} \n"
        MarsInterp.py $AVERAGE -t pstd -l phalf_mb --include temp -e phalf_mb
        echo -e "\n ${g} ${cmark}  Interpolation COMPLETE: ${y}pstd phalf_mb${g} created ${y}${AVERAGE%.*}_pstd_phalf_mb.nc${nc} \n"
    fi

    echo -e "\n\n ${bb}${sep}${nc}"
    echo -e "\n ${bb}>>> Outputting grid information for ${y}p44${bb}, not actually doing interpolation"
    echo -e "    ${y}(i.e. MarsInterp.py $AVERAGE -t pstd -l p44 -g)${nc} \n"
    MarsInterp.py $AVERAGE -t pstd -l p44 -g
    echo -e "\n ${g} ${cmark}  COMPLETE${nc} \n"
fi
# ===============================================================================================


# ====================================== Add Variables =========================================
echo -e "\n\n ${bb}${sep}${nc}"
echo -e "${bb}      ${frame} TESTING MARSVARS FUNCTIONS ${frame} \n"

if $AVERAGE_STATUS; then 
    echo -e "${bb}>>> Removing ${y}r${bb} to ${y}$AVERAGE${bb}...${nc} \n"
    MarsVars.py $AVERAGE -rm cldcol
    echo -e "\n ${g} ${cmark}  Removal COMPLETE: ${y}cldcol${g} from ${y}$AVERAGE${nc} \n"
    echo -e "\n\n ${bb}${sep}${nc}"

    echo -e "${bb}>>> Adding ${y}DP & DZ${bb} to ${y}$AVERAGE${bb}...${nc} \n"
    MarsVars.py $AVERAGE -add DP DZ
    echo -e "\n ${g} ${cmark}  Add COMPLETE: ${y}DP${g} & ${y}DZ${g} to ${y}$AVERAGE${nc} \n"
    echo -e "\n\n ${bb}${sep}${nc}"
    
    echo -e "${bb}>>> Column-Integrate ${y}cld${bb} in ${y}$AVERAGE${bb}...${nc} \n"
    MarsVars.py $AVERAGE -col cld
    echo -e "\n ${g} ${cmark}  Column integration COMPLETE: ${y}cld${g} in ${y}$AVERAGE${nc} \n"
    echo -e "\n\n ${bb}${sep}${nc}"

    echo -e "${bb}>>> Do a ${y}zonal detrend${bb} on ${y}ucomp${bb}...${nc} \n"
    MarsVars.py $AVERAGE -zd ucomp
    echo -e "\n ${g} ${cmark}  Zonal detrending COMPLETE: ${y}ucomp_p${g} added to ${y}$AVERAGE${nc} \n"
    echo -e "\n\n ${bb}${sep}${nc}"

    echo -e "${bb}>>> Editing variable attributes."
    echo -e "       - Multiply ${y}ps${bb} by ${y}0.01${bb} to convert Pa -> mb;"
    echo -e "       - Change the variable name to ${y}Pmb${bb} and longname to ${y}Pressure mb${bb};"
    echo -e "       - Change the unit from ${y}Pa${bb} to ${y}mb${nc}"
    MarsVars.py $AVERAGE --edit ps -multiply 0.01 -rename Pmb -longname 'Pressure mb' -unit 'mbar'
    echo -e "\n ${g} ${cmark}  Attributes for ${y}ps${g} changed in ${y}$AVERAGE${g}. See the output below:${nc} \n"
    
    echo -e "${bb}>>> Testing ${y}inspect${bb} from MarsPlot...${nc} \n"
    MarsPlot.py -i $AVERAGE
    echo -e "\n ${g} ${cmark}  MarsPlot ${y}inspect${g} function works${nc} \n"

    echo -e "\n\n ${bb}${sep}${nc}"
    echo -e "\n ${bb}>>> Revert those changes back...${nc} \n"
    MarsVars.py $AVERAGE --edit Pmb -multiply 100 -rename ps -longname 'surface pressure' -unit 'Pa'
    echo -e "\n\n ${bb}${sep}${nc}"

    echo -e "${bb}>>> Differentiate ${y}temp${bb} with respect to the Z axis in ${y}$AVERAGE${bb}...${nc} \n"
    MarsVars.py $AVERAGE -zdiff temp
    echo -e "\n ${g} ${cmark}  Differentiation COMPLETE: ${y}(temp in $AVERAGE)${nc}"
    echo -e "\n\n ${bb}${sep}${nc}"

    echo -e "${bb}>>> Convert aerosol opacities between ${y}op/m${bb} and ${y}op/Pa${bb} in ${y}$AVERAGE${bb}...${nc} \n"
    MarsVars.py $AVERAGE -dp_to_dz taudust_IR # op/Pa -> op/m; fi
    echo -e "\n ${g} ${cmark}  1/2 conversion to taudust_IR_dp_to_dz COMPLETE: ${y}op/Pa -> op/m${nc} \n"

    MarsVars.py $AVERAGE -dz_to_dp taudust_IR # op/m  -> op/Pa; fi
    echo -e "\n ${g} ${cmark}  2/2 conversion to taudust_IR_dz_to_dp COMPLETE: ${y}op/m -> op/Pa${nc}"
    echo -e "\n\n ${bb}${sep}${nc}"

    echo -e "${bb}>>> Extracting ${y}r${bb} from ${y}$PAVERAGE${bb} and saving them to ${y}${PAVERAGE%.*}_extract.nc${bb}...${nc} \n"
    MarsVars.py $PAVERAGE -extract r
    echo -e "${g} ${cmark}  EXTRACT COMPLETE: moved ${y}r${g} to ${y}${PAVERAGE%.*}_extract.nc${nc} \n"

    if file_exists ${COMBINE1}; then
        if file_exists ${COMBINE2}; then
            echo -e "${bb}>>> Combining ${y}${COMBINE1}${bb} and ${y}${COMBINE2}${bb} into ${y}${COMBINE1}${bb}...${nc}"
            echo -e "${bb}    Note ${y}${COMBINE1}${bb} and ${y}${COMBINE2}${bb} are 250 days long:${nc}\n"
            MarsPlot.py -i ${COMBINE1}
            MarsFiles.py *combine*.nc -c
            echo -e "${g} ${cmark}  Combine COMPLETE: ${y}${COMBINE1}${bb} is now 500 days long:${nc} \n"
            MarsPlot.py -i ${COMBINE1}
        else
            echo -e "${r}>>> ${y}${COMBINE2}${bb} does not exist, cannot test ${y}MarsFiles.py --combine${nc}"
        fi
    else
        echo -e "${r}>>> ${y}${COMBINE1}${bb} does not exist, cannot test ${y}MarsFiles.py --combine${nc}"
    fi
fi
# ===============================================================================================


# ===================================== Try to Query Data =======================================
# echo -e "\n\n ${bb}${sep}${nc}"
# echo -e "${bb}      ${frame} TESTING MARSPULL FUNCTIONS ${frame} \n"

# echo -e "${bb}>>> Querying Legacy GCM data from the NAS Data Portal and porting it into Legacy_Data_By_Ls...${nc} \n"
# mkdir Legacy_Data_By_Ls && cd Legacy_Data_By_Ls
# echo -e "${bb}>>> Pulling data by Ls...${nc} \n"
# MarsPull.py -id INERTCLDS -ls 90. 
# echo -e "\n ${g} ${cmark}  COMPLETE: Data saved in ${y}Legacy_Data_By_Ls/${nc}"
# cd ..
# echo -e "\n\n ${bb}${sep}${nc}"

# echo -e "${bb}>>> Pulling data by ID and porting it into Legacy_Data_By_ID...${nc} \n"
# mkdir Legacy_Data_By_ID && cd Legacy_Data_By_ID
# MarsPull.py -id INERTCLDS -f LegacyGCM_Ls000_Ls004.nc LegacyGCM_Ls005_Ls009.nc
# echo -e "\n ${g} ${cmark}  COMPLETE: Data saved in ${y}Legacy_Data_By_ID/${nc}"
# ===============================================================================================


# ================================= Legacy File Manipulations ===================================
# echo -e "\n\n ${bb}${sep}${nc}"
# echo -e "${bb}      ${frame} TESTING MARSFILES FUNCTIONS RELATED TO MARSPULL ${frame} \n"

# echo -e "${bb}>>> Creating ${y}fixed${bb}, ${y}average${bb}, ${y}daily${bb}, and ${y}diurn${bb} files from Legacy GCM data in Legacy_Data_By_ID...${nc} \n"
# MarsFiles.py LegacyGCM*.nc -fv3 fixed average daily diurn
# echo -e "\n ${g} ${cmark}  ${y}fixed${g}, ${y}average${g}, ${y}daily${g}, and ${y}diurn${g} COMPLETE ${nc}"
# echo -e "\n\n ${bb}${sep}${nc}"

# echo -e "${bb}>>> Combining those files...${nc} \n"
# MarsFiles.py $AVERAGE --combine
# echo -e "\n ${g} ${cmark}  Combine COMPLETE: ${nc}"
# cd ..

# ===============================================================================================



# ===================================== File Manipulations ======================================
echo -e "\n\n ${bb}${sep}${nc}"
echo -e "${bb}      ${frame} TESTING MARSFILES FUNCTIONS ${frame} \n"

echo -e "${bb}>>> Time-shifting MGCM data to ${y}3 AM/3 PM${bb} local time...${nc} \n"
if $DIURN_STATUS; then MarsFiles.py $DIURN -t '3. 15.'; fi # target local times only
echo -e "\n ${g} ${cmark}  Time-shift COMPLETE ${nc} \n"
echo -e "\n\n ${bb}${sep}${nc}"

echo -e "${bb}>>> Binning ${y}daily${bb} files to create ${y}average${bb} files (5-sol averages)...${nc} \n"
if $DAILY_STATUS; then MarsFiles.py $DAILY -bd -ba 5; fi # 5-day bin avereage file
echo -e "\n ${g} ${cmark}  Binning COMPLETE: ${y}Daily${g} -> ${y}average${nc} \n"
echo -e "\n\n ${bb}${sep}${nc}"

if $PAVERAGEc48_STATUS; then
    if $PAVERAGE_STATUS; then 
        echo -e "${bb}>>> Regridding ${y}c24${bb} average file to ${y}c48${bb}...${nc} \n"
        MarsFiles.py ${PAVERAGE%.*}_extract.nc -rs $PAVERAGEc48
        echo -e "\n ${g} ${cmark}  Regridding COMPLETE: ${y}c24${g} -> ${y}c48${nc} \n"
        echo -e "\n\n ${bb}${sep}${nc}"
    fi
fi

echo -e "${bb}>>> Computing zonal average file from ${y}$AVERAGE${bb}, saving to ${y}${AVERAGE%.*}_zonal_avg_CAP_diag.nc${bb}...${nc} \n"
if $AVERAGE_STATUS; then
    MarsFiles.py $AVERAGE -za -e CAP_diag
fi
echo -e "\n ${g} ${cmark}  Zonal averaging COMPLETE: ${y}${AVERAGE%.*}_zonal_avg_CAP_diag.nc${nc} \n"
# ===============================================================================================


# ========================================== Filtering ==========================================
echo -e "\n\n ${bb}${sep}${nc}"
echo -e "${bb}      ${frame} TESTING MARSFILES FILTERING FUNCTIONS ${frame} \n"

if $DAILY_STATUS; then 
    if file_exists ${DAILY%.*}_bpf.nc; then
        echo -e "\n ${bb}>>> ${y}${AVERAGE%.*}_bpf.nc${bb} exists already, skipping the filter. \n"
        echo -e "\n\n ${bb}${sep}${nc}"
    else
        echo -e "${bb}>>> Testing ${y}band-pass${bb} filter on ${y}daily${bb} file. ${y}sol_min=1${bb}, ${y}sol_max=10${bb}...${nc} \n"
        MarsFiles.py $DAILY -bpf 1 10 # sol min, sol max
        echo -e "\n ${g} ${cmark}  ${y}Band-pass${g} filter COMPLETE ${nc} \n"
        echo -e "\n\n ${bb}${sep}${nc}"
    fi

    if file_exists ${DAILY%.*}_hpf.nc; then
        echo -e "\n ${bb}>>> ${y}${AVERAGE%.*}_hpf.nc${bb} exists already, skipping the filter. \n"
        echo -e "\n\n ${bb}${sep}${nc}"
    else
        echo -e "${bb}>>> Testing ${y}high-pass${bb} filter on ${y}daily${bb} file. ${y}sol_min=1${bb}...${nc} \n"
        MarsFiles.py $DAILY -hpf 1 --no_trend # sol min, amplitudes only
        echo -e "\n ${g} ${cmark}  ${y}High-pass${g} filter COMPLETE ${nc} \n"
        echo -e "\n\n ${bb}${sep}${nc}"
    fi

    if file_exists ${DAILY%.*}_lpf.nc; then
        echo -e "\n ${bb}>>> ${y}${AVERAGE%.*}_lpf.nc${bb} exists already, skipping the filter. \n"
        echo -e "\n\n ${bb}${sep}${nc}"
    else
        echo -e "${bb}>>> Testing ${y}low-pass${bb} filter on ${y}daily${bb} file. ${y}sol_min=10${bb}...${nc} \n"
        MarsFiles.py $DAILY -lpf 10 # sol max
        echo -e "\n ${g} ${cmark}  ${y}Low-pass${g} filter COMPLETE ${nc} \n"
        echo -e "\n\n ${bb}${sep}${nc}"
    fi
fi
# ===============================================================================================


# ======================================== Tide analysis ========================================
echo -e "\n\n ${bb}${sep}${nc}"
echo -e "${bb}      ${frame} TESTING MARSFILES TIDAL FUNCTIONS ${frame} \n"

# tide analysis (1=diurn, 2=semi-diurn, 3=ter-diurn...)
if $PDIURN_STATUS; then 
    echo -e "${bb}>>> Tidal analysis on ${y}diurn${bb} file: 4 harmonics, ${y}ps${bb} only...${nc} \n"
    MarsFiles.py $PDIURN -tidal 4 --include ps # 4 harmonics, only do for ps
    echo -e "\n ${g} ${cmark}  COMPLETE: 4 harmonics, ${y}ps${g} only ${nc}"
    echo -e "\n\n ${bb}${sep}${nc}"

    echo -e "${bb}>>> Tidal analysis on ${y}diurn${bb} file: reconstruct first 6 harmonics, ${y}ps${bb} only...${nc} \n"
    MarsFiles.py $PDIURN -tidal 6 --reconstruct --include ps # reconstruct first 6 harmonics
    echo -e "\n ${g} ${cmark}  COMPLETE: reconstructed 6 harmonics, ${y}ps${g} only ${nc}"
    echo -e "\n\n ${bb}${sep}${nc}"

    echo -e "${bb}>>> Tidal analysis on ${y}diurn${bb} file: normalize 6 harmonics, ${y}ps${bb} only...${nc} \n"
    MarsFiles.py $PDIURN -tidal 6 --normalize --include ps # amplitude=%
    echo -e "\n ${g} ${cmark}  COMPLETE: normalized 6 harmonics, ${y}ps${g} only ${nc}"
    echo -e "\n\n ${bb}${sep}${nc}"
fi
# ===============================================================================================


# ====================================== MarsPlot Functions =====================================
echo -e "\n\n ${bb}${sep}${nc}"
echo -e "${bb}      ${frame} TESTING MARSPLOT FUNCTIONS ${frame} \n"

if $AVERAGE_STATUS; then
    echo -e "${bb}>>> Testing ${y}dump${bb} on ${y}average${bb} file...${nc} \n"
    MarsPlot.py -i $AVERAGE -dump ts[0,:2,:2]
    echo -e "\n ${g} ${cmark}  Dump COMPLETE${nc}"
    echo -e "\n\n ${bb}${sep}${nc}"

    echo -e "${bb}>>> Testing ${y}stat${bb} on ${y}average${bb} file...${nc} \n"
    MarsPlot.py -i $AVERAGE -stat ts
    echo -e "\n ${g} ${cmark}  Stat COMPLETE${nc}"
fi
# ===============================================================================================


# ==================================== Custom.in Functions ======================================
echo -e "\n\n ${bb}${sep}${nc}"
echo -e "${bb}      ${frame} TESTING PLOTTING FUNCTIONS ${frame} \n"

rm Custom*.in

echo -e "${bb}>>> Creating ${y}Custom.in${bb}...${nc} \n"
MarsPlot.py --template
if file_exists Custom.in; then
    echo -e "\n ${g} ${cmark}  Templates created${nc}"
else
    echo -e "${r}>>> ${y}Custom.in${bb} was not created ${nc}"
fi
echo -e "\n\n ${bb}${sep}${nc}"

echo -e "${bb}>>> Running ${y}Custom.in${bb} in debug mode...${nc} \n"
MarsPlot.py Custom.in -d 01336 --debug # outputs errors
echo -e "\n ${g} ${cmark}  Ran ${y}Custom.in${g} in ${y}debug${g} mode${nc}"
echo -e "\n\n ${bb}${sep}${nc}"

if file_exists test_files/StackYear.in; then
    echo -e "${bb}>>> Running ${y}StackYear.in${bb} with ${y}--stack_year${bb} on multiple files ${y}(00668-01336)${bb}...${nc} \n"
    cp test_files/StackYear.in .
    MarsPlot.py StackYear.in -d 668 1336 -sy
    echo -e "\n ${g} ${cmark}  Ran ${y}StackYear.in${g} with ${y}--stack_year${bb} on multiple files ${y}(00668-01336)${bb}${nc}"
else
    echo -e "${r}>>> ${y}test_files/StackYear.in${bb} does not exist ${nc}"
fi
echo -e "\n\n ${bb}${sep}${nc}"

cp test_files/CAP_Test_All_Plots.in ~amesCAP-dev/mars_templates/.
if ls ~amesCAP-dev/mars_templates/CAP_Test_All_Plots.in 1> /dev/null 2>&1; then
    echo -e "${bb}>>> Running ${y}CAP_Test_All_Plots.in${bb} using ${y}-do${bb}...${nc} \n"
    MarsPlot.py -do CAP_Test_All_Plots.in -d 01336 # searches ~/amesCAP/mars_templates, then /u/mkahre...shared_templates/
    echo -e "\n ${g} ${cmark}  Ran ${y}CAP_Test_All_Plots.in${nc}"
else
    echo -e "${r}>>> ${y}~amesCAP-dev/mars_templates/CAP_Test_All_Plots.in${bb} does not exist. You probably do NOT have CAP installed in your home directory and you need to move CAP_Test_All_Plots.in to your mars_templates directory so that CAP can find it.${nc}"
fi
echo -e "\n\n ${bb}${sep}${nc}"

echo -e "${bb}>>> Running ${y}Custom.in${bb} to output ${y}PNG${bb}...${nc} \n"
MarsPlot.py Custom.in -d 01336 -o png -pw 500 # sets pixel width=500, default 2000
echo -e "\n ${g} ${cmark}  Ran ${y}Custom.in${g} and output a ${y}PNG${nc}"
echo -e "\n\n ${bb}${sep}${nc}"

echo -e "${bb}>>> Running ${y}Custom.in${bb} in ${y}portrait (-vert)${bb} mode...${nc} \n"
MarsPlot.py Custom.in -d 01336 -vert # portrait mode
echo -e "\n ${g} ${cmark}  Ran ${y}Custom.in${g} in ${y}portrait${g} mode${nc}"
echo -e "\n\n ${bb}${sep}${nc}"

echo -e "${bb}>>> Pointing ${y}Custom.in${bb} to files not in current directory...${nc} \n"
MarsPlot.py Custom.in -d 01336 -dir $PTH/test_files # use files not in pwd
echo -e "\n ${g} ${cmark}  Ran ${y}Custom.in${g} pointing to files not in current directory${nc}"
echo -e "\n\n ${bb}${sep}${nc}"

# ===============================================================================================

rm inspect.out

# ===============================================================================================
echo -e "\n ${g}                ${frame} CAP_Diagnostics.sh COMPLETE! ${frame} "
echo -e "\n ${bb}${sep}${nc}"
# ===============================================================================================


