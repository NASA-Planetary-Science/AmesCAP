# !/usr/bin/env python3
"""
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
"""

# Load generic Python modules
import os
import re
import tempfile
import uuid
from io import BytesIO
from subprocess import Popen, PIPE
from PIL import Image

def convert_from_path(pdf_path, dpi=200, output_folder=None, first_page=None,
                      last_page=None, fmt="ppm", thread_count=1, userpw=None,
                      use_cropbox=False):
    """
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
    """
    page_count = __page_count(pdf_path, userpw)

    if thread_count < 1:
        thread_count = 1

    if first_page is None:
        first_page = 1

    if (last_page is None or
        last_page > page_count):
        last_page = page_count

    # Recalculate page count based on first and last page
    page_count = last_page - first_page + 1

    if thread_count > page_count:
        thread_count = page_count

    reminder = page_count % thread_count
    current_page = first_page
    processes = []
    for _ in range(thread_count):
        # A unique identifier for our files if the directory is not empty
        uid = str(uuid.uuid4())
        # Get the number of pages the thread will be processing
        thread_page_count = page_count // thread_count + int(reminder > 0)
        # Build the command accordingly
        args, parse_buffer_func = __build_command(
            ["pdftoppm", "-r", str(dpi), pdf_path], output_folder,
            current_page, (current_page + thread_page_count - 1), fmt, uid,
            userpw, use_cropbox)
        # Update page values
        current_page = current_page + thread_page_count
        reminder -= int(reminder > 0)
        # Spawn the process and save its uuid
        processes.append((uid, Popen(args, stdout = PIPE, stderr = PIPE)))

    images = []
    for uid, proc in processes:
        data, _ = proc.communicate()
        if output_folder is not None:
            images += __load_from_output_folder(output_folder, uid)
        else:
            images += parse_buffer_func(data)
    return images

def convert_from_bytes(pdf_file, dpi=200, output_folder=None, first_page=None,
                       last_page=None, fmt='ppm', thread_count=1, userpw=None,
                       use_cropbox=False):
    """
    Convert PDF to Image will throw an error whenever one of the condition is reached

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
    """
    with tempfile.NamedTemporaryFile('wb') as f:
        f.write(pdf_file)
        f.flush()
        return convert_from_path(f.name,
                                 dpi = dpi,
                                 output_folder = output_folder,
                                 first_page = first_page,
                                 last_page = last_page,
                                 fmt = fmt,
                                 thread_count = thread_count,
                                 userpw = userpw,
                                 use_cropbox = use_cropbox)

def __build_command(args, output_folder, first_page, last_page, fmt, uid,
                    userpw, use_cropbox):
    if use_cropbox:
        args.append("-cropbox")
    if first_page is not None:
        args.extend(["-f", str(first_page)])
    if last_page is not None:
        args.extend(["-l", str(last_page)])
    parsed_format, parse_buffer_func = __parse_format(fmt)
    if parsed_format != "ppm":
        args.append("-" + parsed_format)
    if output_folder is not None:
        args.append(os.path.join(output_folder, uid))
    if userpw is not None:
        args.extend(["-upw", userpw])
    return args, parse_buffer_func

def __parse_format(fmt):
    if fmt[0] == ".":
        fmt = fmt[1:]
    if (fmt == "jpeg" or
        fmt == "jpg"):
        return "jpeg", __parse_buffer_to_jpeg
    if fmt == "png":
        return "png", __parse_buffer_to_png

    # Unable to parse the format so use the default
    return "ppm", __parse_buffer_to_ppm

def __parse_buffer_to_ppm(data):
    images = []
    index = 0
    while index < len(data):
        code, size, rgb = tuple(data[index:index+40].split(b"\n")[0:3])
        size_x, size_y = tuple(size.split(b" "))
        file_size = (len(code) + len(size) + len(rgb) + 3
                     + int(size_x) * int(size_y) * 3)
        images.append(Image.open(BytesIO(data[index:index + file_size])))
        index += file_size
    return images

def __parse_buffer_to_jpeg(data):
    # Last element is obviously empty
    return [Image.open(BytesIO(image_data + b"\xff\xd9")) for image_data
            in data.split(b"\xff\xd9")[:-1]]

def __parse_buffer_to_png(data):
    images = []
    index = 0
    while index < len(data):
        # 4 bytes for IEND + 4 bytes for CRC
        file_size = data[index:].index(b"IEND") + 8
        images.append(Image.open(BytesIO(data[index:index+file_size])))
        index += file_size
    return images

def __page_count(pdf_path, userpw=None):
    try:
        if userpw is not None:
            proc = Popen(["pdfinfo", pdf_path, "-upw", userpw], stdout = PIPE,
                         stderr = PIPE)
        else:
            proc = Popen(["pdfinfo", pdf_path], stdout = PIPE, stderr = PIPE)
        out, err = proc.communicate()
    except:
        raise Exception("Unable to get page count. If not on a Linux system, "
                        "please install the Poppler PDF rendering library")

    try:
        # This will throw an error if we are unable to get page count
        return int(re.search(r"Pages:\s+(\d+)",
                             out.decode("utf8", "ignore")).group(1))
    except:
        raise Exception(f"Unable to get page count. "
                        f"{err.decode('utf8', 'ignore')}")

def __load_from_output_folder(output_folder, uid):
    return [Image.open(os.path.join(output_folder, f)) for f
            in sorted(os.listdir(output_folder)) if uid in f]
