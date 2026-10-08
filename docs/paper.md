---
title: 'Community Analysis Pipeline: A Python package for processing Mars climate model data'
tags:
  - Python
  - astronomy
  - Mars global climate model
  - data processing
  - data visualization
authors:
  - name: Courtney M. L. Batterson
    orcid: 0000-0001-5894-095X
    equal-contrib: true
    affiliation: 1
  - name: Richard A. Urata
    orcid: 0000-0001-8497-5718
    equal-contrib: true
    affiliation: 1
    corresponding: true
  - name: Victoria L. Hartwick
    orcid: 0000-0002-2082-8986
    equal-contrib: true
    affiliation: 3
  - name: Alexandre M. Kling
    orcid: 0000-0002-2980-7743
    equal-contrib: true
    affiliation: "1, 4"
  - name: Melinda A. Kahre
    orcid: 0000-0002-0935-5532
    equal-contrib: true
    affiliation: 2
affiliations:
 - name: Bay Area Environmental Research Institute, United States
   index: 1
   ror: 024tt5x58
 - name: NASA Ames Research Center, United States
   index: 2
   ror: 02acart68
 - name: Southwest Research Institute, United States
   index: 3
   ror: 03tghng59
 - name: Astera Institute, United States
   index: 4
   ror: 00ydx1s47
date: 2 October 2026
bibliography: paper.bib

---

# Summary

The Community Analysis Pipeline (CAP) is a Python package designed to streamline and simplify the complex process of analyzing the large datasets output by global climate models (GCMs). CAP consists of a suite of command-line tools that manipulate netCDF files to produce secondary datasets and figures useful for science and engineering. By converting output from several Mars GCMs and reanalyses to a common format, CAP facilitates inter-model and model-to-observation comparison for Mars. The goal is to enable users of varying levels of programming experience to work more easily with complex data products from GCMs, thereby lowering the barrier to entry into planetary science research.

# Statement of need

GCMs perform numerical simulations that describe the evolution of climate systems on planetary bodies. GCM data products routinely include surface and atmospheric variables such as wind, temperature, and aerosol concentrations. While GCMs have been applied to planetary bodies in our Solar System (e.g., Earth, Venus, Pluto) and in other stellar systems (e.g., @Hartwick2023), CAP focuses on Mars GCMs (MGCMs). Several MGCMs are actively in use and under development in the Mars community, including the NASA Ames MGCM (Legacy and FV3-based versions, the latter built on the GFDL AM4 framework [@Zhao2018]), NASA Goddard ROCKE-3D [@Way2017], the Laboratoire de Météorologie Dynamique (LMD) Mars Planetary Climate Model (PCM) [@Forget1999], MarsWRF [@Richardson2007; @Newman2019], the Tohoku University DRAMATIC Mars GCM [@Kuroda2005; @Kobayashi2026], the Max Planck Institute Mars GCM [@Hartogh2005], and GEM-Mars [@Neary2018], as well as the OpenMARS [@Holmes2020] and EMARS [@Greybush2019] reanalyses. CAP was originally developed as a suite of tools specific to the NASA Ames MGCM [@Haberle2019; @Bertrand2020; @AmesGCM2023], covering everything from plotting in Mars conventions (e.g., solar longitude (L$_s$) rather than an Earth calendar) to vertical interpolation and derived variables that depend on Mars constants such as gravity and the specific heat of the atmosphere ($C_p$). Because we viewed such tools as a need shared across the Mars climate modeling community, we later added compatibility with other MGCMs: CAP reads output from the NASA Ames MGCM directly and, through its `MarsFormat` tool, converts output from PCM, MarsWRF, OpenMARS, and EMARS (see the MarsFormat documentation at <https://amescap.readthedocs.io>). Model-specific variable and dimension names are defined in a user-editable variable dictionary, so additional models can be supported without changes to the rest of the pipeline.

GCM output is complex in both size and structure, and analyzing the data requires GCM-specific domain knowledge. We highlight the following major challenges for working with MGCM output:

- Files are complex in structure with output fields represented by multiple variables (e.g., air and surface temperature) with varying units (e.g., Kelvin, Celsius) in multi-dimensional structures (e.g., 2–5 dimensions) at a variety of sampling frequencies (e.g., temporally averaged, instantaneous) and on custom horizontal and vertical grids.
- Output from a simulation of one Mars year ranges from \~10 GB–10 TB, depending on the fields saved, time sampling, and resolution. Such files require memory-aware processing, which is particularly challenging for users without access to clusters or supercomputers.
- Domain-specific knowledge is required to derive secondary variables, manipulate complex data structures, and visualize results. Such information is not always publicly available, and many MGCMs do not output data in self-describing formats like netCDF.

CAP provides libraries and command-line executables for file manipulation, analysis, and visualization. It automates routine and sophisticated post-processing for experienced modelers and removes technical roadblocks for new users of MGCM data.

# State of the Field

General-purpose netCDF viewers such as Panoply [@Schmunk2024], Ncview [@Pierce2024], the Grid Analysis and Display System (GrADS; @GrADS), and ParaView [@Kitware2023] offer quick visualization of gridded data, and some (e.g., GrADS expressions and ParaView's Python interface) support user-defined computations. For Earth science, Python packages such as xarray [@Hoyer2017], Iris [@Iris], GeoCAT [@GeoCAT2020], and ESMValTool [@Righi2020] provide labeled-array analysis, model evaluation, and plotting, while planetary tools such as planetoplot [@planetoplot], aeolus [@Sergeev2021], and wrf-python [@Ladwig2017] address plotting or diagnostics for particular models. CAP builds on this ecosystem: it uses NumPy [@Harris2020], Matplotlib [@Hunter2007], netCDF4 [@Unidata2024], and xarray rather than reimplementing array storage, file input/output, or plotting, and it adapts the wrf-python approach to destaggering MarsWRF winds.

What CAP adds is a Mars-specific layer that these packages do not provide as a single workflow: conversion of several MGCM and reanalysis formats to one convention; Mars time handling (sols, L$_s$, and local true solar time); derived variables that need MGCM vertical coordinates and Mars constants (e.g., mass streamfunction, column-integrated and vertically differentiated quantities); vertical interpolation from hybrid sigma-pressure levels to standard pressure, altitude, or height above ground; and template-driven plotting, all behind one command-line interface. We built CAP as a separate package because we needed a complete toolkit for the NASA Ames MGCM, from Mars-specific plotting and vertical interpolation to derived variables computed with Mars constants, and the Earth-oriented packages assume Earth calendars, constants, and model conventions (e.g., CF-compliant time coordinates and Earth model output standards). Support for other MGCMs was added afterward, once we saw a broader need for such tools in the Mars climate modeling community.

# Software Design

CAP is written in Python for its approachability (for developers and users), open-source-friendly design, and integration of third-party packages (NumPy, netCDF4, xarray, Matplotlib). CAP is portable to Linux, macOS, and Windows, which matters because model data often moves between systems, whether from a supercomputer to local storage or between collaborators. A compiled language such as C++ could also combine analysis and visualization in one package (e.g., with VTK for visualization and Eigen or xtensor for array computation). We chose Python because the intended users are planetary scientists, many of whom analyze data in Python. Python lets them read, modify, and extend CAP's functions and plotting templates without a compilation step, and the scientific Python stack already provides netCDF input/output and publication-quality plotting. The cost is lower computational efficiency, which CAP mitigates by vectorizing calculations with NumPy and by processing large files in chunks (e.g., `MarsInterp` interpolates in time chunks and `MarsFormat` can use dask for lazy, chunked processing).

Several design choices follow from the structure of Mars model output:

- **Time.** The netCDF Climate and Forecast (CF) conventions define time relative to Earth calendars. Mars days (sols) are about 2.7% longer than Earth days, and a Mars year has about 668.6 sols, so no CF calendar represents Mars time. CAP therefore stores time as sols since the start of the simulation and keeps L$_s$ in a separate `areo` variable. `MarsCalendar` converts between sols and L$_s$, and `MarsFormat` computes L$_s$ for converted files that lack it. Diurnal ("diurn") files store each hour of the day in universal time, so the local time differs with longitude; `time_shift_calc` in `FV3_utils.py` shifts these fields to uniform local solar time, including an equation-of-time correction when the model provides it, which `MarsFiles` exposes through its `-t` option.
- **Common format.** All tools operate on files in the NASA Ames FV3-based MGCM layout. `MarsFormat` maps other models onto this layout, so that `MarsFiles`, `MarsVars`, `MarsInterp`, and `MarsPlot` need not handle each model separately.
- **Executables and libraries.** Each executable in `bin/` (`MarsPull`, `MarsFiles`, `MarsVars`, `MarsInterp`, `MarsPlot`, `MarsFormat`, `MarsCalendar`, `MarsNest`) is a script that parses its arguments, processes one file at a time, and writes a new file, while shared numerical routines (e.g., vertical coordinates, interpolation, and time shifting in `FV3_utils.py`, and spectral filtering in `Spectral_utils.py`, which optionally uses SHTools [@Wieczorek2018]) are plain functions that can be imported into user scripts. Executables use module-level constants and settings, for example the planetary constants in `MarsVars`, which are reset for each input file from the file's metadata. This keeps each executable readable as a linear script, at the cost of making those routines harder to reuse than the library functions.
- **Writing files.** Reading is functional, using netCDF4 and xarray directly, while writing is object-oriented: the `Ncdf` class in `Ncdf_wrapper.py` wraps a netCDF4 dataset and provides methods to copy dimensions and variables, log new variables with their metadata, and copy large variables in chunks. Tools such as `MarsFiles` (e.g., when concatenating or splitting files) and `MarsInterp` build their output through `Ncdf`, so metadata handling is shared, while the processing itself remains in functions. The related `Fort` class converts binary output from the Legacy NASA Ames MGCM to the same netCDF layout.

CAP is tested with unit tests of numerical routines and integration tests of each executable on synthetic files that mimic each supported model, run on Linux, macOS, and Windows with GitHub Actions.

# Research Impact Statement

CAP is developed by the NASA Ames Mars Climate Modeling Center alongside the NASA Ames MGCM. It was the analysis tool taught at the NASA Ames MGCM and CAP Tutorial held virtually on 13–15 November 2023, whose exercises are part of CAP's documentation and can now be completed either with the publicly released NASA Ames MGCM simulations, which `MarsPull` downloads from the NASA Advanced Supercomputing Data Portal, or with synthetic files generated by CAP's test suite.

CAP has begun to be adopted outside the development team. In 2026, Mark I. Richardson, first author of the PlanetWRF model description [@Richardson2007], contributed support for current PlanetWRF output, native staggered winds, true local solar time, and faster chunked interpolation and conversion.

Within the Mars Climate Modeling Center, CAP has been used to analyze NASA Ames MGCM simulations in studies including uniform local time shifting, vertical interpolation to a standard altitude grid, and migrating tide analysis [@Urata2025], adding streamfunctions[@Batterson2023], preliminary analysis figures [@Hartwick2022a; @Hartwick2022b; @Batterson2023; @Urata2025], figures for conference presentations [@Kahre2022; @Kahre2023; @Steakley2023; @Steakley2024; @Hartwick2024].

Because CAP has not previously had a DOI or a peer-reviewed description, studies that used it could not cite it formally, and it does not appear in their reference lists. This paper, together with an archived, versioned release, gives CAP a citable reference for future work.

# AI Use Disclosure

Claude (Sonnet 4.5) was used to assist with software execution error analysis and for the development of the automated testing suite via GitHub Actions. During the review of this paper, Claude (Claude Code with Opus 5.5) was used to help implement reviewer-requested changes, including the MarsPull download fixes, compact test fixtures, regression tests, and documentation and manuscript revisions. Quality and correctness of AI-generated content were verified by comparing results to expected output, by running the automated test suite, and by review from multiple team members.

# Acknowledgements

This work is supported by the Planetary Science Division of the National Aeronautics and Space Administration as part of the Mars Climate Modeling Center funded by the Internal Scientist Funding Model.

# References
