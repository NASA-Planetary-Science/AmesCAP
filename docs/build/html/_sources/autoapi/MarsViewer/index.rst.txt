:py:mod:`MarsViewer`
====================

.. py:module:: MarsViewer

.. autoapi-nested-parse::

   The MarsViewer executable is for ...

   Alex 11/29/18
   This is a pure-Python implementation of a primitive image viewer.
   It is designed to display an image preview in the terminal where accessing NAS through 
   X11 is not possible or is slow. Supported formats are PNG, EPS, and PDF.

   Adapted from:
   https://github.com/Belval/pdf2image

   The executable requires two arguments:
       * [-id --id]      the simulation identifier, AND
       * [-ls --ls]      the desired solar longitude(s), OR
       * [-f -filename]  the name(s) of the desired file(s)

   Third-party Requirements:
       * numpy
       * argparse
       * requests

   List of Functions:
       * download - queries the requested file from the NAS Data Portal



