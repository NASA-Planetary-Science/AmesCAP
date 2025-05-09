---  
title: 'Community Analysis Pipeline: A Python package for processing Mars climate model data'  
tags:  
  - Python  
  - astronomy  
  - Mars global climate model  
  - data processing  
  - data visualization  
authors:  
  - name: Alexandre M. Kling  
    orcid: 0000-0002-2980-7743
    equal-contrib: true  
    affiliation: 1  
    corresponding: true \# (This is how to denote the corresponding author)  
  - name: Courtney M. L. Batterson  
    orcid: 0000-0001-5894-095X  
    equal-contrib: true  
    affiliation: 1  
  - name: Richard A. Urata  
    orcid: 0000-0001-8497-5718  
    equal-contrib: true  
    affiliation: 1  
  - name: Victoria L. Hartwick  
    orcid: 0000-0002-2082-8986  
    equal-contrib: true  
    affiliation: 3  
  - name: Melinda A. Kahre  
    orcid: 0000-0002-0935-5532
    equal-contrib: true  
    affiliation: 2  
affiliations:  
 - name: Bay Area Environmental Research Institute, United States  
   index: 1  
 - name: NASA Ames Research Center, United States  
   index: 2  
 - name: Southwest Research Institute, United States  
   index: 3  
date: 22 April 2025  
bibliography: paper.bib

---

# Summary

The Community Analysis Pipeline (CAP) is a Python package designed to streamline and simplify the complex process of analyzing large datasets created by global climate models (GCMs). CAP consists of a suite of tools that manipulate NetCDF files in order to produce secondary datasets and figures useful for science and engineering applications. CAP also facilitates inter-model and model-observation comparisons, and it is the first software of its kind to standardize these comparisons. The goal is to enable users with varying levels of programming experience to work with complex data products from a variety of GCMs and thereby lower the barrier to entry for planetary science research. 

# Statement of need

GCMs perform numerical simulations that describe the evolution of climate systems on planetary bodies. GCMs simulate physical processes within the atmosphere (and, if applicable, within the surface of the planet, ocean, and any interactions therein), calculate radiative transfer within those mediums, and use a computational fluid dynamics (CFD) solver (the “dynamical core”) to predict the transport of heat and momentum within the atmosphere. Typical GCM products include surface and atmospheric variables such as wind, temperature, and aerosol concentrations. While GCMs have been applied to planetary bodies in our Solar System (e.g. Earth, Venus, Pluto) and in other stellar systems (e.g. [@Hartwick:2023]), CAP is currently compatible with Mars GCMs (MGCMs). Several MGCMs are actively in use and under development in the Mars community, including the NASA Ames MGCM (Legacy and FV3-based versions), NASA Goddard ROCKE-3D, the Laboratoire de Météorologie Dynamique (LMD) Mars Planetary Climate Model (PCM), the Open University OpenMars, NCAR MarsWRF, NCAR MarsCAM, GFDL Mars GCM, Harvard DRAMATIC Mars GCM, Max Planck Institute Mars GCM, and GEM-Mars. Of these, CAP is compatible with four models so far: the NASA Ames MGCM, PCM, OpenMars, and MarsWRF. 

MGCM output is complex in both size and structure. Analyzing the output requires GCM-specific domain knowledge. We identify the following major challenges for working with MGCM output:

\* Files tend to be fairly complex in structure, with output fields represented by multiple variables (e.g. air vs surface temperature), varying units (e.g. Kelvin), complex dimensional structures (e.g. 2–5 dimensions), and a variety of sampling frequencies (e.g. temporally averaged or instantaneous) on different horizontal and vertical grids.  
\* File sizes typically range from \~10 Gb–10 Tb for simulations describing the Martian climate over a full orbit around the Sun (depending on the number of atmospheric fields being analyzed, time sampling, and the horizontal and vertical resolutions of the run). Large files require curated processing pipelines in order to manage memory storage. This can be particularly challenging for users that do not have access to academic or enterprise clusters or supercomputers for their analyses.  
\* Domain-specific knowledge is required to derive secondary variables, manipulate complex data structures, and visualize results. Working with MGCM data is especially difficult for users unfamiliar with the fields commonly output by MGCMs or the mathematical methods used in climate science.

CAP offers a streamlined workflow for processing and analyzing MGCM data products by providing a set of libraries and executables that facilitate file manipulation and data visualization from the command-line. This benefits existing modelers by automating both routine and sophisticated post-processing tasks. It also expands access to MGCM products by removing some of the technical roadblocks associated with processing these complex data products.

CAP is the first software package to provide data visualization, file manipulation, variable derivation, and inter-model or model-observation comparison features in one software suite. Existing tools perform a subset of the functions that CAP offers, but none of them provide both complex data analysis tools and visualization, nor are they specifically designed for climate modeling. Some of the more popular tools include Panoply [@Schmunk:2024], Ncview [@Pierce:2024], Grid Analysis and Display System (GrADS; @GMU), and Paraview [@Kitware:2023]. Each tool offers simple solutions for visualizing NetCDF data and some provide minimal flexibility for user-defined computations. However, CAP is the only software package with an open-source Python framework for analyzing and plotting climate data and performing inter-model and model-observation comparisons.

CAP has been used in multiple research projects that have been published and/or presented at conferences worldwide (e.g., [@Urata:2025, @Batterson:2023, @Hartwick:2022a, @Hartwick:2022b, @Steakley:2023, @Steakley:2024, @Kahre:2022, @Kahre:2023, @Nagata:2025, @Hartwick:2024]).

# Functionality

CAP consists of six command-line executables that can be used sequentially or individually to derive secondary data products, thus offering a high level of flexibility. A configuration text file is provided so that users can define the input file structure (e.g., variable names, longitudinal structure, and interpolation levels) and preferred plotting style (e.g., time axis units) for their analysis. The six executables in CAP are described below: 

## MarsPull

MarsPull is a data pipeline utility for downloading MGCM data products from the NAS Data Portal ([https://data.nas.nasa.gov/](https://data.nas.nasa.gov/)) . Recognizing that each member within the science and engineering community has their own requirements for hosting proprietary Mars climate datasets (e.g. institutional servers, Zenodo, GitHub, etc.), MarsPull is intended to be a mechanism for interfacing those datasets. MarsPull enables users to query data meeting a specific criteria, such as a date range (e.g., solar longitude), which allows users to parse repositories first and download only the necessary data, thus avoiding downloading entire repositories which can be large (\>\>15Gb). A typical application of MarsPull is:

`> MarsPull directory_name -f MGCM_file1.nc MGCM_file2.nc` 

## MarsFormat

MarsFormat is a utility for converting non-NASA Ames MGCM products into NASA Ames-like MGCM products for compatibility with CAP. MarsFormat reorders dimensions, adds standardized coordinates that are expected by other executables for various computations (e.g., pressure interpolation), converts variable units to conform to the International System of Units (e.g., Pa for pressure), and reorganizes coordinate values as needed (e.g., reversing the vertical pressure array for plotting). Additional, model-specific operations are performed as necessary. For example, MarsWRF data requires un-staggering latitude-longitude grids and calculating absolute fields from perturbation fields. A typical application of MarsFormat is:

`> MarsFormat MGCM_file.nc -gcm model_name`

## MarsFiles

MarsFiles provides several tools for file manipulation such as file size reduction, temporal and spatial filtering, and splitting or concatenating data along specified dimensions. Operations performed by MarsFiles are applied to entire NetCDF files producing new data structures with amended file names. A typical application of MarsFiles is:

`> MarsFiles MGCM_file.nc -flags`

## MarsVars

MarsVars performs variable operations such as adding, removing, and editing variables and computing column integrations. It is standard practice within the modeling community to avoid outputting variables that can be derived outside of the MGCM in order to minimize file size. For example, atmospheric density (rho) is easily derived from temperature and pressure and therefore typically not included in output files. MarsVars derives rho from temperature and pressure and adds it to the file with a single command line argument. A typical application of MarsVars is:

`> MarsVars MGCM_file.nc –add rho`

## MarsInterp

MarsInterp interpolates the vertical coordinate to a standard grid: pressure, altitude, or altitude above ground level. Vertical grids vary considerably from model to model. Most MGCMs use a pressure or hybrid pressure vertical coordinate (e.g. terrain-following, pure pressure levels, or sigma levels) in which the geometric heights and mid-layer pressures of the atmospheric layers vary in latitude and longitude. It is therefore necessary to interpolate to a standard vertical grid in order to do any rigorous spatial averaging or inter-model or observation-to-model comparisons. A typical application of MarsInterp is:

`> MarsInterp MGCM_file.nc -t pstd`

## MarsPlot

MarsPlot is the plotting utility for CAP. It accepts a modifiable text template containing a list of plots to generate (Custom.in) as input and outputs graphics to PDF or PNG. It supports multiple types of 1-D or 2-D plots, color schemes, map projections, and can customize axes range, plot titles, or contour intervals. It also supports some simple math functions to derive secondary fields not supported by MarsVars. A typical application of MarsPlot is:

`> MarsPlot Custom.in`

# Acknowledgements

This work is supported by the Planetary Science Division of the National Aeronautics and Space Administration as part of the Mars Climate Modeling Center funded by the Internal Scientist Funding Model.

# References

