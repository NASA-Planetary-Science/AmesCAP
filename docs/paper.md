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
    affiliation: 1
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
date: 9 May 2025
bibliography: paper.bib

---

# Summary

The Community Analysis Pipeline (CAP) is a Python package designed to streamline and simplify the complex process of analyzing the large datasets output by global climate models (GCMs). CAP consists of a suite of tools that manipulate NetCDF files to produce secondary datasets and figures useful for science and engineering. CAP facilitates inter-model and model-to-observation comparison, and it is the first software of its kind to standardize these comparisons. The goal is to enable users of varying levels of programming experience to work more easily with complex data products from GCMs, thereby lowering the barrier to entry into planetary science research.

# Statement of need

GCMs perform numerical simulations that describe the evolution of climate systems on planetary bodies. GCMs simulate physical processes within the atmosphere (and, if applicable, within the surface of the planet, ocean, and any interactions between them), calculate radiative transfer within those mediums, and predict the transport of heat and momentum within the atmosphere using a computational fluid dynamics (CFD) solver (the “dynamical core”). GCM data products routinely include surface and atmospheric variables such as wind, temperature, and aerosol concentrations. While GCMs have been applied to planetary bodies in our Solar System (e.g., Earth, Venus, Pluto) and in other stellar systems (e.g., @Hartwick2023), CAP is extends compatibility to Mars GCMs (MGCMs). Several MGCMs are actively in use and under development in the Mars community, including the NASA Ames MGCM (Legacy and FV3-based versions), NASA Goddard ROCKE-3D, the Laboratoire de Météorologie Dynamique (LMD) Mars Planetary Climate Model (PCM), the Open University OpenMars, NCAR MarsWRF, NCAR MarsCAM, GFDL Mars GCM, Harvard DRAMATIC Mars GCM, Max Planck Institute Mars GCM, and GEM-Mars. To date, CAP is compatible with four of these models (the NASA Ames MGCM, PCM, OpenMars, and MarsWRF) and undergoing development to further extend compatibility to other models.

GCM output is complex in both size and structure, and analyzing the data requires GCM-specific domain knowledge. We highlight the following major challenges for working with MGCM output:

- Files are complex in structure with output fields represented by multiple variables (e.g., air and surface temperature) with varying units (e.g., Kelvin, Celcius) in multi-dimensional structures (e.g., 2–5 dimensions) at a variety of sampling frequencies (e.g., temporally averaged, instantaneous) and on custom horizontal and vertical grids.
- File sizes range from \~10 Gb–10 Tb for simulations describing the Martian climate over a full orbit around the Sun, and vary significantly depending on the number of atmospheric fields being analyzed, time sampling, and horizontal and vertical resolutions. Large files require curated processing pipelines to manage memory storage, which can be particularly challenging for users that do not have access to academic or enterprise clusters or supercomputers for their analyses.
- Domain-specific knowledge is required to derive secondary variables, manipulate complex data structures, and visualize results. Such information is not always publicly available online and many MGCMs do not output data into self-describing formats like netCDF. Working with MGCM data is especially difficult for users unfamiliar with the fields commonly output by MGCMs or the mathematical methods used in climate science.

CAP offers a streamlined workflow for processing and analyzing MGCM data products by providing a set of libraries and executables that facilitate file manipulation and data visualization from the command line. This benefits existing modelers by automating both routine and sophisticated post-processing tasks. It also expands access to MGCM products by removing some of the technical roadblocks associated with processing these complex data products.



# State of the Field
<!-- A description of how this software compares to other commonly-used packages in the research area. If related tools exist, provide a clear “build vs. contribute” justification explaining your unique scholarly contribution and why existing alternatives are insufficient. -->
CAP is the first software package to provide data visualization, file manipulation, variable derivation, and inter-model or model-to-observation comparison tools for MGCM output in one software suite. Some existing tools perform a subset of these functions, but none of them provide both complex data analysis and visualization tools and they are not designed specifically for climate modeling. Some of the more popular tools include Panoply [@Schmunk2024], Ncview [@Pierce2024], Grid Analysis and Display System (GrADS; @GMU), and Paraview [@Kitware2023]. Each tool offers simple solutions for visualizing NetCDF data. Some provide minimal flexibility for user-defined computations. However, existing tools are either Java or C based, which preclude the use of powerful analysis and visualization packages available in Python such as NumPy and Matplotlib. CAP is the only software package with an open-source Python framework for analyzing and plotting climate model data and performing inter-model or model-to-observation comparison. 

# Software Design
<!-- An explanation of the trade-offs you weighed, the design/architecture you chose, and why it matters for your research application. This should demonstrate meaningful design thinking beyond a superficial code structure description. -->
CAP was originally designed to be an accessible toolkit for manipulating, analyzing, and visualizing data from the NASA Ames MGCM. It is written in Python for its approachability (for developers and users), open-source-friendly design, and integration of third party packages (i.e., NumPy, NetCDF4, Matplotlib). CAP is easily portable to multiple operating systems, which is a key feature of the software as the model data often moves between systems, whether from a cluster to local storage or between users who have differing operating systems. The use of Python over other languages (e.g., C/C++, or Fortran) means CAP sacrifices some computational efficiency for the many benefits of having all of the analysis and visualization tools in a single, open-source package.

CAP tools are grouped by functionality for an intuitive experience. For example, file manipulation functions (like splitting a dataset by date or slicing by latitude) are stored in MarsFiles, functions that add diagnostic variables are in MarsVars, and plotting is handled by MarsPlot. The main executables (MarsFiles, MarsVars, MarsPlot, MarsInterp, MarsFormat, MarsPull, MarsCalendar) contain readable instructions and function definitions, many of which are accessible via the command line with the `-h` flag, while the more complicated backend functions are stored in separate files (e.g., FV3_utils.py, Ncdf_wrapper.py). 

# Research Impact Statement
<!-- Evidence of realized impact (publications, external use, integrations) or credible near-term significance (benchmarks, reproducible materials, community-readiness signals). The evidence should be compelling and specific, not aspirational. -->
CAP has been used in multiple research projects that have been published and/or presented at conferences worldwide (e.g., @Urata2025; @Batterson2023; @Hartwick2022a; @Hartwick2022b; @Steakley2023; @Steakley2024; @Kahre2022; @Kahre2023; @Nagata2025; @Hartwick2024).

# AI Use Disclosure
<!-- Transparent disclosure of any use of generative AI in the software creation, documentation, or paper authoring. If no AI tools were used, state this explicitly. If AI tools were used, describe how they were used and how the quality and correctness of AI-generated content was verified. -->
Claude (Sonnet 4.5) was used to assist with software execution error analysis and for the development of the automated testing suite via GitHub Actions. Quality and correctness of AI-generated content was verified by comparing results to expected output, and was tested by multiple team members.

# Acknowledgements

This work is supported by the Planetary Science Division of the National Aeronautics and Space Administration as part of the Mars Climate Modeling Center funded by the Internal Scientist Funding Model.

# References
