echo " 
**********************************************************************
This is the documentation of the command-line MCMC tools available in:
/u/mkahre/MCMC/analysis/bin/
**********************************************************************

MarsDocumentation
		--> This file, quick documentation of the tools available.

MarsTips        
                --> Summary of useful command line functions, everybody is welcomed to contribute to this file
                   
MarsDir
		--> simple alias to the directory /u/mkahre/MCMC/analysis/

MarsCalendar
		--> Tanguy's SOL to Ls converter 
                USAGE: MarsCalendar 750.
                       MarsCalendar start end step
                       MarsCalendar --help
MarsPlot
		--> Simple vizualization toolkit for FV3
                USAGE: MarsPlot --help (in the FV3/verona/simuID/history directory)
                   
MarsViewer                       
	       --> Display a pdf, png or eps preview through the terminal as a pixelated image (ugly but does not require X11)
               If you terminal is 100 character large, you will get a 100x50 image. 
               (Consider opening a dedicated terminal with a small font (e.g fontsize =5)
               USAGE: MarsViewer image.png
                      MarsViewer Diagnostics.pdf

MarsLayers     --> Create a new set of vertical layers for use in fv_eta.F90
               USAGE: MarsLayers 38 (the number of layers)
                      MarsLayers 38 --save (save the ak and bk as a Netcdf file in the home directory)
                      MarsLayers 38 --help

MarsVars     --> Add variable (like the density), remove variable or differentiate variable from the diagnostic files. 
               USAGE: MarsVars 02400.atmos_average.nc -add rho
                      MarsVars --help
                      (IN PRORESS)


                  
"
