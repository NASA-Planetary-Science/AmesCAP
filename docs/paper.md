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
    corresponding: true
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

The Community Analysis Pipeline (CAP) is a Python package designed to streamline and simplify the complex process of analyzing large datasets created by global climate models (GCMs). CAP consists of a suite of tools that manipulate NetCDF files in order to produce secondary datasets and figures useful for science and engineering applications. CAP also facilitates inter-model and model-observation comparisons, and it is the first software of its kind to standardize these comparisons. The goal is to enable users with varying levels of programming experience to work with complex data products from a variety of GCMs and thereby lower the barrier to entry for planetary science research.

# Statement of need

GCMs perform numerical simulations that describe the evolution of climate systems on planetary bodies. GCMs simulate physical processes within the atmosphere (and, if applicable, within the surface of the planet, ocean, and any interactions therein), calculate radiative transfer within those mediums, and use a computational fluid dynamics (CFD) solver (the “dynamical core”) to predict the transport of heat and momentum within the atmosphere. Typical GCM products include surface and atmospheric variables such as wind, temperature, and aerosol concentrations. While GCMs have been applied to planetary bodies in our Solar System (e.g., Earth, Venus, Pluto) and in other stellar systems (e.g., @Hartwick2023), CAP is currently compatible with Mars GCMs (MGCMs). Several MGCMs are actively in use and under development in the Mars community, including the NASA Ames MGCM (Legacy and FV3-based versions), NASA Goddard ROCKE-3D, the Laboratoire de Météorologie Dynamique (LMD) Mars Planetary Climate Model (PCM), the Open University OpenMars, NCAR MarsWRF, NCAR MarsCAM, GFDL Mars GCM, Harvard DRAMATIC Mars GCM, Max Planck Institute Mars GCM, and GEM-Mars. Of these, CAP is compatible with four models so far: the NASA Ames MGCM, PCM, OpenMars, and MarsWRF.

MGCM output is complex in both size and structure. Analyzing the output requires GCM-specific domain knowledge. We identify the following major challenges for working with MGCM output:

- Files tend to be fairly complex in structure, with output fields represented by multiple variables (e.g., air vs surface temperature), varying units (e.g., Kelvin), complex dimensional structures (e.g., 2–5 dimensions), and a variety of sampling frequencies (e.g., temporally averaged or instantaneous) on different horizontal and vertical grids.
- File sizes typically range from \~10 Gb–10 Tb for simulations describing the Martian climate over a full orbit around the Sun (depending on the number of atmospheric fields being analyzed, time sampling, and the horizontal and vertical resolutions of the run). Large files require curated processing pipelines in order to manage memory storage. This can be particularly challenging for users that do not have access to academic or enterprise clusters or supercomputers for their analyses.
- Domain-specific knowledge is required to derive secondary variables, manipulate complex data structures, and visualize results. Working with MGCM data is especially difficult for users unfamiliar with the fields commonly output by MGCMs or the mathematical methods used in climate science.

CAP offers a streamlined workflow for processing and analyzing MGCM data products by providing a set of libraries and executables that facilitate file manipulation and data visualization from the command-line. This benefits existing modelers by automating both routine and sophisticated post-processing tasks. It also expands access to MGCM products by removing some of the technical roadblocks associated with processing these complex data products.



# State of the Field
<!-- A description of how this software compares to other commonly-used packages in the research area. If related tools exist, provide a clear “build vs. contribute” justification explaining your unique scholarly contribution and why existing alternatives are insufficient. -->
CAP is the first software package to provide data visualization, file manipulation, variable derivation, and inter-model or model-observation comparison features in one software suite. Existing tools perform a subset of the functions that CAP offers, but none of them provide both complex data analysis tools and visualization, nor are they specifically designed for climate modeling. Some of the more popular tools include Panoply [@Schmunk2024], Ncview [@Pierce2024], Grid Analysis and Display System (GrADS; @GMU), and Paraview [@Kitware2023]. Each tool offers simple solutions for visualizing NetCDF data and some provide minimal flexibility for user-defined computations. However, the listed existing tools are either Java or C based, which precludes the use of powerful analysis and visualization packages available in Python such as NumPy and Matplotlib. CAP is the only software package with an open-source Python framework for analyzing and plotting climate model data and performing inter-model or model-observation comparisons. 

# Software Design
<!-- An explanation of the trade-offs you weighed, the design/architecture you chose, and why it matters for your research application. This should demonstrate meaningful design thinking beyond a superficial code structure description. -->
CAP was designed to provide an easily accessible toolkit for the manipulation, analysis, and visualization of data originating from the NASA Ames MGCM. Python was chosen as the basis for the software toolkit for its ease-of-use (for developers and users) and the availability of useful third party packages (i.e., NumPy, NetCDF4, Matplotlib). The software being easily portable to multiple systems is a key feature for this toolkit, as the model data can quickly move between different systems, or be shared between users with different computing environments. While this decision to use Python over other languages (e.g., C/C++, or Fortran) sacrifices some computational efficiency with some of the more intense analyses, the benefits of having the analysis tools and visualization tools all in one easy-to-install package outweigh the negatives.

The tools in the software package are grouped by general functionality to make it easier to use. For example, file manipulation functions such as splitting a dataset by date or slicing by latitude are in MarsFiles, functions to add diagnostic variables are in MarsVars, and plotting is handled by MarsPlot. We try to maintain readability of the main modules by moving backend functions to separate files (e.g., FV3_utils.py, or Ncdf_wrapper.py). 

# Research Impact Statement
<!-- Evidence of realized impact (publications, external use, integrations) or credible near-term significance (benchmarks, reproducible materials, community-readiness signals). The evidence should be compelling and specific, not aspirational. -->
CAP has been used in multiple research projects that have been published and/or presented at conferences worldwide (e.g., @Urata2025; @Batterson2023; @Hartwick2022a; @Hartwick2022b; @Steakley2023; @Steakley2024; @Kahre2022; @Kahre2023; @Nagata2025; @Hartwick2024).

# AI Use Disclosure
<!-- Transparent disclosure of any use of generative AI in the software creation, documentation, or paper authoring. If no AI tools were used, state this explicitly. If AI tools were used, describe how they were used and how the quality and correctness of AI-generated content was verified. -->
Claude (Sonnet 4.5) was used to assist with software execution error analysis and for the development of the automated testing suite via GitHub Actions. Quality and correctness of AI-generated content was verified by comparing results to expected output, and was tested by multiple team members.

# Acknowledgements

This work is supported by the Planetary Science Division of the National Aeronautics and Space Administration as part of the Mars Climate Modeling Center funded by the Internal Scientist Funding Model.

# References

