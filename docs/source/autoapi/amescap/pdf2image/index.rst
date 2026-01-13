:py:mod:`amescap.pdf2image`
===========================

.. py:module:: amescap.pdf2image

.. autoapi-nested-parse::

   pdf2image is a light wrapper for the poppler-utils tools that can
   convert PDFs into Pillow images.

   Reference: https://github.com/Belval/pdf2image

   Third-party Requirements:

       * ``io``
       * ``tempfile``
       * ``re``
       * ``os``
       * ``subprocess``
       * ``PIL``
       * ``uuid``
       * ``Pillow``



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   amescap.pdf2image.convert_from_bytes
   amescap.pdf2image.convert_from_path



.. py:function:: convert_from_bytes(pdf_file, dpi=200, output_folder=None, first_page=None, last_page=None, fmt='ppm', thread_count=1, userpw=None, use_cropbox=False)

   Convert PDF to Image will throw an error whenever one of the
   condition is reached

   :param pdf_file: Bytes representing the PDF file
   :type pdf_file: float
   :param dpi: image quality in DPI (default 200)
   :type dpi: int
   :param output_folder: folder to write the images to (instead of
       directly in memory)
   :type output_folder: str
   :param first_page: first page to process
   :type first_page: int
   :param last_page: last page to process before stopping
   :type last_page: int
   :param fmt: output image format
   :type fmt: str
   :param thread_count: how many threads to spawn for processing
   :type thread_count: int
   :param userpw: PDF password
   :type userpw: str
   :param use_cropbox: use cropbox instead of mediabox
   :type use_cropbox: bool


.. py:function:: convert_from_path(pdf_path, dpi=200, output_folder=None, first_page=None, last_page=None, fmt='ppm', thread_count=1, userpw=None, use_cropbox=False)

   Convert PDF to Image will throw an error whenever one of the
   conditions is reached.

   :param pdf_path: path to the PDF that you want to convert
   :type pdf_path: str
   :param dpi: image quality in DPI (default 200)
   :type dpi: int
   :param output_folder: folder to write the images to (instead of
       directly in memory)
   :type output_folder: str
   :param first_page: first page to process
   :type first_page: int
   :param last_page: last page to process before stopping
   :type last_page: int
   :param fmt: output image format
   :type fmt: str
   :param thread_count: how many threads to spawn for processing
   :type thread_count: int
   :param userpw: PDF password
   :type userpw: str
   :param use_cropbox: use cropbox instead of mediabox
   :type use_cropbox: bool


