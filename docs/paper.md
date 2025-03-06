---
title: 'Community Analysis Pipeline: A Python package for processing Mars climate model data'
tags:
  - Python
  - astronomy
  - Mars global climate model
  - data processessing
  - data visualization
authors:
  - name: Alexandre M. Kling
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: 1 
    corresponding: true # (This is how to denote the corresponding author)
  - name: Courtney M. Batterson
    orcid: 0000-0000-0000-0000
    equal-contrib: true 
    affiliation: 1
  - name: Richard A. Urata
    orcid: 0000-0001-8497-5718
    equal-contrib: true 
    affiliation: 1
  - name: Melinda A. Kahre
    orcid: 0000-0000-0000-0000
    equal-contrib: true 
    affiliation: 2
affiliations:
 - name: Bay Area Environmental Research Institute, United States
   index: 1
 - name: NASA Ames Research Center, United States
   index: 2
date: 4 March 2025
bibliography: paper.bib

---

# Summary
The Community Analysis Pipeline (CAP) is a Python package designed to streamline and simplify the complex process of analyzing large datasets that are created by global climate models (GCMs). The package consists of a suite of tools that easily manipulate NetCDF climate data with the goal of producing desired figures for validation and comparison for scientists with any level of prior programming experience. CAP currently consists of six tools. *MarsPull* accesses and downloads publicly available Mars GCM data so independently created data is not required. *MarsFormat* converts a selection of Mars GCM output to be compatible with CAP functions. *MarsFiles* performs various file operations such as splitting and concatenation. *MarsVars* processes available data to derive secondary variables such as column integration, potential temperature, relative vorticity, etc. *MarsInterp* performs vertical interpolation on the output to convert the model native grid to a standard pressure or altitude grid. *MarsPlot* visualizes the output as standard 2-dimensional or 1-dimensional plots, with flexibility for different plot types, as well as averaging options.


# Statement of need
Global climate models perform computational simulations to describe the evolution of climate systems on planetary bodies.  GCMs  simulate physical processes within the atmosphere (and if applicable, within the surface of the planet, ocean and any interactions therein), calculate radiative transfer within those medium, and use a computational fluid dynamics solver (the “dynamical core”) to calculate the transport of heat and momentum within the atmosphere. Typical products of the GCM describe the evolutions of surface and atmosphere variables such as winds, temperature and aerosols. While they have been applied to many objects in our Solar System (e.g. Earth, Venus, Pluto) and in other stellar systems (e.g. [REF]), this tool specifically focuses on Mars Global Climate Models (MGCM).
Mars Global Climate Model (MGCM) outputs can be complex to analyze both in terms of sizes, structures, and also require knowledge-based processing tools for certain applications. We identify the following major challenges for new users:
* The outputs tend to be fairly complex in structure, with output fields of various meanings (e.g air vs surface temperature), units (e.g. Kelvin),  dimensions (e.g. 2-dimensional to 5-dimensional), sampling (e.g. temporally averaged) on various modeling grids.
* Typical disk file sizes range between ~10 Gb-10Tb  for simulations describing the Martian climate over a full orbit around the Sun (depending on the choices of atmospheric fields being analyzed, time sampling, horizontal and vertical resolution of the run) which requires intentional handling of the files to prevent memory issues. This can be particularly challenging for users that do not readily have access to academic or entreprise clusters/supercomputer for their analysis.
* Knowledge-based tools are required to adapt those files when a user requires secondary variables or specific data structures not immediately available through the native (raw) MGCM files. 

The objective of CAP is to provide a streamlined workflow for processing and analyzing Mars Global Climate Model (MGCM) products. This benefits existing modelers by automating both routine and  more sophisticated post-processing tasks. It also allows to expand access to MGCM products to new users by removing some of the roadblocks traditionally associated with handling those complex products, such as interpreting large datasets, deriving secondary fields, and smoothly visualizing the data. CAP accomplishes this by providing a set of libraries and executables that facilitate the processing of the MGCM files from the command-line, as well as providing a visualization tool that can be used to generate scientific plots. 

Compared to the traditional “one-step” scripted programming approach that requires extended  modifications to execute any particular analysis, the executables and the visualization tool in CAP do not require extensive programming knowledge. This makes CAP accessible to users of many backgrounds and various skill levels. A feature of CAP is its extensive input/output (I/O) capabilities which allow executables to generate or update intermediary files at each processing step. This offers a high degree of flexibility regarding the choice and extent of the post-processing applied to the raw MGCM outputs. This enables users of various backgrounds to process MGCM products to the degree they need for their particular application from the command-line. In addition, users have the possibility either to directly import CAP libraries within their own-python based applications (without using the executables), or generate files for use in external applications in other programming languages. 

We conclude this summary of CAP by listing some of its design attributes: 

* CAP is written in Python, an open-source programming language with extensive scientific libraries available.  
* The files are written in the netCDF4 data format developed by UCAR/Unidata [@Unidata2024]. This format is widely used in the climate modeling community and is self-descriptive (meaning that a file contains metadata about its content in terms of variables names, units, etc.)
* CAP uses a convention for output formatting inherited from the Geophysical Fluid Dynamics Laboratory (GFDL) Finite-Volume Cubed-Sphere Dynamical Core, referred here as "FV3 format" [@Zhao2018]. Specifically, outputs may be binned and averaged in time in various ways for analysis.
* CAP offers multi-model support. At the time of the writing, actively supported models are the NASA Ames Legacy GCM [@Haberle2019], and the NASA Ames GCM [@Bertrand2020] are supported. CAP also offers compatibility with the openMars [@Holmes2020], MarsWRF [@Newman2019], PCM [REF], EMARS [REF] data products.

#Optional functionality discussion



# Acknowledgements
This work is supported by the Planetary Science Division of the National Aeronautics and Space Administration as part of the Mars Climate Modeling Center funded by the Internal Scientist Funding Model.