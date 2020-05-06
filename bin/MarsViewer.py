#!/usr/bin/env python3
'''
Alex 11/29/18
This is a pure Python implementation of a primitive image viewer. 
It is designed to display an image preview in the terminal where accessing NAS through X11 is not possible or done over a slow internet connection
Supported format are png, eps and pdf.

Adapted from:
https://github.com/Belval/pdf2image
'''
import sys
import numpy as np
import os

def prCyan(skk): print("\033[96m{}\033[00m" .format(skk)) 
def prYellow(skk): print("\033[93m{}\033[00m" .format(skk))

try:
	from PIL import Image
except ImportError as error_msg:
	prYellow("Error, this function requires the Python Image Library (Pillow)")
	prYellow('Please install Pillow with:') 
	prCyan('        pip install Pillow\n')
	print("Error was: "+ error_msg.message)
	exit()

from amesgcm.pdf2image import convert_from_bytes

def get_ansi_color_code(r, g, b):
	if r == g and g == b:
		if r < 8:
			return 16
		if r > 248:
			return 231
		return round(((r - 8) / 247) * 24) + 232
	return 16 + (36 * round(r / 255 * 5)) + (6 * round(g / 255 * 5)) + round(b / 255 * 5)


def get_color(r, g, b):
	return "\x1b[48;5;{}m \x1b[0m".format(int(get_ansi_color_code(r,g,b)))


def show_image(img_path):
	try:
		# If pdf, merge image together
		if img_path[-3:]=='pdf' :
			print('Extracting pdf pages...')	
			images = convert_from_bytes(open(img_path, 'rb').read())
			widths, heights = zip(*(i.size for i in images))
			total_height = sum(heights)
			max_width = max(widths)
			
			img = Image.new('RGB', (max_width, total_height))
			
			y_offset = 0
			for im in images:
				img.paste(im, (0,y_offset))
				y_offset += im.size[1]
		else: 		
			img = Image.open(img_path)
			
	except Exception as e:
		exit(e)
	rows, columns = os.popen('stty size', 'r').read().split()
	mywidth = int(columns)
	scale_factor=0.5 #account for height of line in terminal
	
	wpercent = (mywidth/float(img.size[0]))
	hsize = int(scale_factor*(float(img.size[1])*float(wpercent)))
	img = img.resize((mywidth,hsize), Image.ANTIALIAS)
	

	img_arr = np.asarray(img)
	h,w,c = img_arr.shape

	for x in range(h):
		for y in range(w):
			pix = 1.*img_arr[x][y]
			txt_col=get_color(pix[0], pix[1], pix[2])
			sys.stdout.write(txt_col)
		print('')


if __name__ == '__main__':
	if len(sys.argv) > 1:
		img_path = sys.argv[1]
		show_image(img_path)
