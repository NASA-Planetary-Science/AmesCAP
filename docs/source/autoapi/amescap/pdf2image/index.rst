:py:mod:`amescap.pdf2image`
===========================

.. py:module:: amescap.pdf2image

.. autoapi-nested-parse::

       pdf2image is a light wrapper for the poppler-utils tools that can convert your
       PDFs into Pillow images.
   Reference:
   https://github.com/Belval/pdf2image



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   amescap.pdf2image.convert_from_path
   amescap.pdf2image.convert_from_bytes



.. py:function:: convert_from_path(pdf_path, dpi=200, output_folder=None, first_page=None, last_page=None, fmt='ppm', thread_count=1, userpw=None, use_cropbox=False)

   Description: Convert PDF to Image will throw whenever one of the condition is reached
   Parameters:
       pdf_path -> Path to the PDF that you want to convert
       dpi -> Image quality in DPI (default 200)
       output_folder -> Write the resulting images to a folder (instead of directly in memory)
       first_page -> First page to process
       last_page -> Last page to process before stopping
       fmt -> Output image format
       thread_count -> How many threads we are allowed to spawn for processing
       userpw -> PDF's password
       use_cropbox -> Use cropbox instead of mediabox 


.. py:function:: convert_from_bytes(pdf_file, dpi=200, output_folder=None, first_page=None, last_page=None, fmt='ppm', thread_count=1, userpw=None, use_cropbox=False)

   Description: Convert PDF to Image will throw whenever one of the condition is reached
   Parameters:
       pdf_file -> Bytes representing the PDF file
       dpi -> Image quality in DPI
       output_folder -> Write the resulting images to a folder (instead of directly in memory)
       first_page -> First page to process
       last_page -> Last page to process before stopping
       fmt -> Output image format
       thread_count -> How many threads we are allowed to spawn for processing
       userpw -> PDF's password
       use_cropbox -> Use cropbox instead of mediabox


