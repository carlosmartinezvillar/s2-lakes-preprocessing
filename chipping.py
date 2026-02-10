'''
This script chips 10980x10980 Sentinel-2 images into 256x256 (or 512x512) images.
'''

import os
import xml.etree.ElementTree as ET
import rasterio as rio
from rasterio.windows import Window
import numpy as np
import matplotlib.pyplot as plt
import glob
import math
import multiprocessing as mp
import time
from PIL import Image
import sys
import argparse

# Typing
from typing import Tuple, List
ndarray = np.ndarray #quick fix...

# XML namespace in METADATA.xml
ns = {
	'n1':"https://psd-14.sentinel2.eo.esa.int/PSD/User_Product_Level-2A.xsd",
	'other':"http://www.w3.org/2001/XMLSchema-instance",
	'another':"https://psd-14.sentinel2.eo.esa.int/PSD/User_Product_Level-2A.xsd"
	}

DATA_DIR = None
LABEL_DIR = None
CHIP_DIR = None
#SET DIRS HERE BECAUSE THREAD ACCESS
# DATA_DIR  = args.data_dir
# LABEL_DIR = DATA_DIR+'/dynamicworld' 
# CHIP_DIR  = args.chip_dir
# if CHIP_DIR is None:
# 	CHIP_DIR  = DATA_DIR+'/chips' #Subdir in same place as .SAFE folders

# PIXEL LIMITS -- I suppose these are good here too bc threads/processes
CHIP_SIZE = 512
# WATER_MIN = 128*128 #1/4 of the image
WATER_MIN = CHIP_SIZE*CHIP_SIZE // 4
WATER_MAX = CHIP_SIZE*CHIP_SIZE-WATER_MIN #balanced for 1/8 land
# BAD_PX    = 3276 #unused
STRIDE    = 256
N_PROC    = 32

####################################################################################################
# CLASSES
####################################################################################################
class EmptyLabelError(Exception): #fancy.
	pass

class IncompleteDirError(Exception):
	pass

class Product():
	'''
	An object referencing a single Sentinel-2 product in the ESA database.

	Parameters
	----------
	id:
	tile:
	date:
	orbit:
	s2_fnames:
	s2_readers:
	gee_id:
	dw_path:
	dw_reader:
	s2_borders:
	dw_borders:
	base_chip_id:

	Methods
	-------
	get_band_filenames()
	get_gee_id()

	'''
	def __init__(self,safe_id):
		self.id    = safe_id
		self.tile  = self.id[38:44]
		self.date  = self.id[11:26]
		self.orbit = self.id[33:37]

		#1.1 ID -> BAND READERS
		self.s2_fnames  = self.get_band_filenames() #sorted
		self.s2_readers = []
		for f in self.s2_fnames:
			band_path = f'{DATA_DIR}/{safe_id}/{f}'
			if not os.path.isfile(band_path):
				raise IncompleteDirError(f"Missing band file {f}")
			self.s2_readers += [rio.open(band_path,'r',tiled=True)]

		#1.2 ID -> XML PATH
		#2.XML -> DW PATH
		#3.DW PATH -> DW READER
		self.gee_id    = self.get_gee_id()
		self.dw_path   = f'{LABEL_DIR}/{self.gee_id}.tif'		
		self.dw_reader = rio.open(self.dw_path,'r',tiled=True)

		#Check label
		if self.dw_reader.statistics(1).max == 0:
			raise EmptyLabelError("Label is zero everywhere.")

		#4.DW READER -> BOUNDS DW
		#5.DW READER+BAND2 READER -> BOUNDS S2 & BOUNDS DW
		self.s2_borders,self.dw_borders = align(self.s2_readers[0],self.dw_reader)
	
		#format: DATE_DSTRIP_TILE_ROTATION_WINROW_WINCOL_B0*.tif
		#format: DATE_DSTRIP_TILE_ROTATION_WINROW_WINCOL_LBL.tif	
		self.base_chip_id = self.gee_id + '_' + self.orbit


	def parse_xml(self):
		assert os.path.isfile(path), "No file found in path %s" % path

		# get datastrip
		root      = ET.parse(path).getroot()
		prod_info = root.find('n1:General_Info',namespaces=ns).find('Product_Info')
		granule   = prod_info.find('Product_Organisation').find('Granule_List').find('Granule')
		datastrip = granule.attrib['datastripIdentifier'].split('_')[-2][1:]

	return datastrip

	def get_band_filenames(self):
		return [f'{self.tile}_{self.date}_{b}_10m.jp2' for b in ['B02','B03','B04','B08']]

	def get_gee_id(self):
		xml_path  = glob.glob(f'{DATA_DIR}/{self.id}/*.xml')[0]
		datastrip = parse_xml(xml_path)
		return '_'.join([self.date,datastrip,self.tile])


####################################################################################################
# STRINGS+PARSING
####################################################################################################
def parse_xml(path):

	# check path
	assert os.path.isfile(path), "No file found in path %s" % path

	# get datastrip
	root      = ET.parse(path).getroot()
	prod_info = root.find('n1:General_Info',namespaces=ns).find('Product_Info')
	granule   = prod_info.find('Product_Organisation').find('Granule_List').find('Granule')
	datastrip = granule.attrib['datastripIdentifier'].split('_')[-2][1:]

	return datastrip


def get_gee_id(s2_id: str) -> str:
	# xml_name    = [f for f in os.listdir(DATA_DIR+'/'+s2_id) if f[-4:]=='.xml'][0]
	# xml_path    = DATA_DIR + '/' + '/'.join([s2_id,xml_name])	
	xml_path  = glob.glob(DATA_DIR + '/' + s2_id + '/*.xml')[0]
	datastrip = parse_xml(xml_path)
	date,tile = s2_id.split('_')[2:6:3]
	gee_id    = '_'.join([date,datastrip,tile])
	return gee_id


####################################################################################################
# RASTER PROCESSING.
####################################################################################################
def remove_label_borders(src: rio.DatasetReader) -> dict:
	'''
	Take a rasterio DatasetReader for a dynamicworld image and get the indices 
	where non-zeros begin at the top, bottom, left, and right.

	Parameters
	----------
	src: rasterio.DatasetReader
		Dataset reader for a dynamic world array (which has zeroes where S2
		still has data, making it redundant to check for zeroes in the S2 array).

	Returns
	-------
	dict
		dictionary with indices of first non-zero values at top, left, right, 
		bottom

	'''
	top    = 0
	bottom = src.height-1
	left   = 0
	right  = src.width-1

	while(True):
		row = src.read(1,window=rio.windows.Window(0,top,src.width,1))
		if row.sum() == 0:
			top += 1
		else:
			break

	while(True):
		row = src.read(1,window=rio.windows.Window(0,bottom,src.width,1))
		if row.sum() == 0:
			bottom -= 1
		else:
			break

	while(True):
		col = src.read(1,window=rio.windows.Window(left,0,1,src.height))
		if col.sum() == 0:
			left += 1
		else:
			break

	while(True):
		col = src.read(1,window=rio.windows.Window(right,0,1,src.height))
		if col.sum() == 0:
			right -= 1
		else:
			break

	return {'top':top, 'bottom':bottom, 'left':left, 'right':right}


def align(s2_src: rio.DatasetReader,dw_src: rio.DatasetReader) -> Tuple:
	'''
	Do everything: match indices and remove borders.
	'''
	# 1. REMOVE DW NO-DATA BORDERS(~1-2px each side)
	dw_ij = remove_label_borders(dw_src) # <---- THIS CAN BE COMBINED

	# 2. MATCH DW to S2 (DW has ~20px less on each side) 
	# DW ij's (px index) -> DW xy's (coords)
	dw_xy_ul = dw_src.xy(dw_ij['top'],dw_ij['left'],offset='center')
	dw_xy_lr = dw_src.xy(dw_ij['bottom'],dw_ij['right'],offset='center')
	# DW xy's (coords) -> S2 ij's (px index)
	s2_ij = {}
	s2_ij['top'],s2_ij['left']     = s2_src.index(dw_xy_ul[0],dw_xy_ul[1],op=math.floor)
	s2_ij['bottom'],s2_ij['right'] = s2_src.index(dw_xy_lr[0],dw_xy_lr[1],op=math.floor)

	# 3. TRIM S2 -- REMOVE S2 TILE OVERLAP & ADJUST DW
	if s2_ij['top'] < 492: #shift top down
		delta        = 492 - s2_ij['top']
		s2_ij['top'] = 492
		dw_ij['top'] = dw_ij['top'] + delta

	if s2_ij['bottom'] > 10487: #shift bottom up
		delta           = s2_ij['bottom'] - 10487
		s2_ij['bottom'] = 10487	
		dw_ij['bottom'] = dw_ij['bottom'] - delta

	if s2_ij['left'] < 492: #shift left right
		delta         = 492 - s2_ij['left']
		s2_ij['left'] = 492	
		dw_ij['left'] = dw_ij['left'] + delta

	if s2_ij['right'] > 10487: #shift right left
		delta          = s2_ij['right'] - 10487
		s2_ij['right'] = 10487		
		dw_ij['right'] = dw_ij['right'] - delta

	return s2_ij,dw_ij	


def get_windows_strided(borders: dict, chip_size: int, stride: int) -> [Tuple]:
	# number of rows and cols takin' the boundaries into acct
	n_px_rows = borders['bottom'] + 1 - borders['top']
	n_px_cols = borders['right'] + 1 - borders['left']

	#nr of overlapping (or not) blocks in each direction
	block_rows = (n_px_rows - chip_size) // stride + 1
	block_cols = (n_px_cols - chip_size) // stride + 1

	#total blocks
	N = block_rows * block_cols

	windows = []

	for k in range(N):
		i = k // block_cols
		j = k % block_cols
		row_start = i * stride + borders['top']
		col_start = j * stride + borders['left']
		W = Window(col_start,row_start,chip_size,chip_image)
		windows += [[(str(i),str(j)),W]]

	return windows


def get_windows(borders: dict) -> [Tuple]:
	'''
	Given a dicts of boundaries, returns an array list with tuples (i,j) for block indices i,j and 
	window objects corresponding to the block i,j while considering only the area of the raster
	within the boundaries defined by the indices in the dict. For example, if the array had two rows
	and a column of no data (top and left) the blocks are offseted and defined as:

			    left    256      512
				| 0 0 ..  	      |
				| 0 0... |		  | 
	    top ----+--------+--------+----
		    0 0 |        |        |
		    0 0 | (0, 0) | (0, 1) |
		     .  |        |        |
		     .  +--------+--------+
		     .  |        |        |
		        | (1, 0) | (1, 1) |
		        |        |        |
		512 ----+--------+--------+---
				|                 |

	Parameters
	----------
	borders: dict
		The dictionary containing the first and last indices of usable data in
		both directions.
	'''

	# number of rows and cols takin' the boundaries into acct
	n_px_rows = borders['bottom'] + 1 - borders['top']
	n_px_cols = borders['right'] + 1 - borders['left']

	#nr of blocks in each direction
	block_rows = n_px_rows // CHIP_SIZE
	block_cols = n_px_cols // CHIP_SIZE

	#total blocks
	N = block_rows * block_cols

	windows = []

	for k in range(N):
		i = k // block_cols
		j = k % block_cols
		row_start = i * CHIP_SIZE + borders['top']
		col_start = j * CHIP_SIZE + borders['left']
		W = Window(col_start,row_start,CHIP_SIZE,CHIP_SIZE)
		windows += [[(str(i),str(j)),W]]

	return windows


def chip_image(product,index,N):
	print(f'[{index}/{N-1}] PROCESSING {product.id} ')
	start_time = time.time()

	# LOAD ARRAYS AND NORMALIZE BANDS
	rgbn = []
	for reader in product.s2_readers:
		# print(f'Loading {reader.name[-34:]}')
		band_array = reader.read(1)
		zero_mask  = band_array == 0
		cutoff     = int(np.percentile(band_array[~zero_mask],99)) # DO NOT PASS FLOAT TO CLIP HERE!!!
		# cutoff     = np.percentile(band_array[~zero_mask],99) #This might have to be lower?
		band_array = np.clip(band_array,0,cutoff)
		band_array = (band_array / cutoff * 255).astype(np.uint8)
		rgbn.append(band_array)

	#SPLIT WINDOWS
	# s2_windows = get_windows(product.s2_borders)
	# dw_windows = get_windows(product.dw_borders)	
	s2_windows = get_windows_strided(product.s2_borders,CHIP_SIZE,STRIDE)
	dw_windows = get_windows_strided(product.dw_borders,CHIP_SIZE,STRIDE)
	share    = len(s2_windows) // N_PROC
	leftover = len(s2_windows) % N_PROC
	start    = [i*share for i in range(N_PROC)]
	stop     = [i*share+share for i in range(N_PROC)]
	stop[-1] += leftover
	s2_chunks = [s2_windows[s0:s1] for s0,s1 in zip(start,stop)]
	dw_chunks = [dw_windows[s0:s1] for s0,s1 in zip(start,stop)]	

	lock = mp.Lock()

	#THROW WORKERS AT ARRAYS
	processes = []
	for i in range(N_PROC):
		p = mp.Process(
			target=chip_image_worker,
			args=(rgbn,product.dw_path,s2_chunks[i],dw_chunks[i],product.base_chip_id,lock)
			)
		p.start()
		processes.append(p)

	for p in processes:
		p.join(timeout=60)
	print("All workers done. ",end='')
	exec_time = time.time() - start_time
	print(f"({exec_time:.3f} seconds).")


def chip_image_worker(rgbn,dw_path,s2_windows,dw_windows,base_id,lock):

	stats = []
	lbl_rdr = rio.open(dw_path,'r',tiled=True)

	for k,(rowcol,w) in enumerate(s2_windows):

		lbl_arr = lbl_rdr.read(1,window=dw_windows[k][1])

		# CHECK LABEL NO DATA
		if (lbl_arr == 0).any():
			continue

		# CHECK WATER/LAND RATIO
		n_water = (lbl_arr==1).sum()
		if n_water < WATER_MIN or n_water > WATER_MAX:
			continue

		# ALL GOOD -- SAVE BANDS
		row = rowcol[0]
		col = rowcol[1]

		# SAVE BANDS IN SINGLE [R,G,B,NIR] FILE (NIR stored in alpha)
		outfile = f'{CHIP_DIR}/{base_id}_{row:02}_{col:02}_B0X.tif'
		r = Image.fromarray(rgbn[0][w.row_off:w.row_off+CHIP_SIZE, w.col_off:w.col_off+CHIP_SIZE])
		g = Image.fromarray(rgbn[1][w.row_off:w.row_off+CHIP_SIZE, w.col_off:w.col_off+CHIP_SIZE])
		b = Image.fromarray(rgbn[2][w.row_off:w.row_off+CHIP_SIZE, w.col_off:w.col_off+CHIP_SIZE])
		n = Image.fromarray(rgbn[3][w.row_off:w.row_off+CHIP_SIZE, w.col_off:w.col_off+CHIP_SIZE])
		img = Image.merge('RGBA',(r,g,b,n))
		img.save(outfile)

		# ALL GOOD -- SAVE LABEL
		outfile = f'{CHIP_DIR}/{base_id}_{row:02}_{col:02}_LBL.tif'
		lbl_arr[lbl_arr!=1] = 0 #everything else (already checked for zeroes above)
		lbl_arr[lbl_arr==1] = 255 #water
		img = Image.fromarray(lbl_arr)
		img.save(outfile)

		stats.append(f'{outfile.split("/")[-1]}\t{n_water}\n')

	# LOG
	lock.acquire()
	# print(f'Worker {mp.current_process()} done.')	
	with open(f'{CHIP_DIR}/stats.txt','a') as fp:
		for line in stats:
			fp.write(line)
	lock.release()


if __name__ == '__main__':

	########## ARGV CONFIG ##########
	parser = argparse.ArgumentParser(
		prog="chipping.py",
		description="Preprocessing (chipping) of Sentinel-2 and DynamicWorld V1 images.")

	# PATHS
	parser.add_argument('--data-dir',default='./dat',
		help="Dataset directory")
	parser.add_argument('--chip-dir',
		help="Chip directory")

	# LOAD
	args = parser.parse_args()

	########## SET ARGS ##########
	DATA_DIR  = args.data_dir
	LABEL_DIR = DATA_DIR+'/dynamicworld' #<-- fix this at some point...
	CHIP_DIR  = args.chip_dir

	if not os.path.isdir(DATA_DIR):
		print("DATA_DIR not found. EXITING.")
		sys.exit()
	print(f"DATA_DIR:  {DATA_DIR}")	

	if len(glob.glob('*.SAFE',root_dir=DATA_DIR)) == 0 :
		print("EMPTY DATA_DIR")
		sys.exit()

	print(f"LABEL_DIR set to: {LABEL_DIR}")
	print(f"CHIP_DIR set to:  {CHIP_DIR}")


	#.SAFE folders in data directory
	folders = glob.glob('*.SAFE',root_dir=DATA_DIR)
	paths   = glob.glob(DATA_DIR+'/*.SAFE')

	# Check everything is there
	if not os.path.isdir(LABEL_DIR):
		print("LABEL_DIR not found. EXITING.")
		sys.exit()

	#make chip dir if not already there
	if not os.path.isdir(CHIP_DIR):
		os.mkdir(CHIP_DIR) 

	# clean log file
	if os.path.isfile(f"{CHIP_DIR}/stats.txt"):
		os.remove(CHIP_DIR+'/stats.txt')


	########## PROCESS .SAFE FOLDERS ##########
	N = len(folders)
	for i,f in enumerate(folders):
		try:
			product = Product(f) #load metadata
		except (EmptyLabelError,IncompleteDirError) as e:
			print(f'ERROR: {e}')
			print(f'---> SKIPPING {f}')
			with open(f'{CHIP_DIR}/errored.txt','a') as fp:
				fp.write(f'{f}\n')
			continue

		# <----- CHIP ----->
		chip_image(product,i,N)

	print("DONE.")