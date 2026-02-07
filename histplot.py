'''
THIS script generates some images and plots to evaluate the quality of the products, labels, and
chips produced.
'''
import os
import zipfile
import rasterio as rio
from rasterio.windows import Window
import matplotlib.pyplot as plt
import glob
import time
import argparse
import numpy as np
import xml.etree.ElementTree as ET
from lxml import etree as LT
import geopandas as gpd
from fiona.drvsupport import supported_drivers
import sys
import pyproj
pyproj.network.set_network_enabled(False) #weird use where to_crs()

import chipping

plt.style.use('fast')
plt.rcParams['font.family'] = 'Courier'
plt.rcParams['font.size'] = 10

####################################################################################################
# WHOLE RASTER PLOTS
####################################################################################################
def plot_raster_lbl_binary(out_path,P,chip_size=512):
	'''
	Plot workable area overlapping windows. Black and white classes.

	Parameters
	----------
	out_path:
	P: chipping.Product() class.
	chip_size:

	'''
	h = P.dw_borders['bottom'] + 1 - P.dw_borders['top']
	w = P.dw_borders['right'] + 1 - P.dw_borders['left']
	windows_height = h - (h % chip_size) + 1
	windows_width  = w - (w % chip_size) + 1

	kwargs = P.dw_reader.meta.copy()
	kwargs.update({'height':windows_height,'width':windows_width,'count':3,'compress':'lzw'})
	out_ptr = rio.open(out_path,'w',**kwargs)

	dw_windows = chipping.get_windows_strided(P.dw_borders,chip_size,stride=chip_size)

	for _,w in dw_windows:
		#read
		arr        = P.dw_reader.read(1,window=w)
		white_mask = arr == 1 #water
		gray_mask  = arr == 0 #nodata
		black_mask = ~(white_mask | gray_mask) #land

		#change
		arr[white_mask] = 255
		arr[black_mask] = 0
		arr[gray_mask]  = 128
		arr_3d = np.repeat(arr[np.newaxis,:,:],repeats=3,axis=0)

		#write
		out_win = Window(w.col_off-P.dw_borders['left'],w.row_off-dw_borders['top'],chip_size,chip_size)
		out_ptr.write(arr_3d,window=out_win)

	out_ptr.close()
	print("IMAGE WRITTEN TO: %s" % out_path)


def plot_raster_lbl_colors():
	DW_PALETTE_10 = ['000000','419bdf','397d49','88b053','7a87c6','e49635', 
	    'dfc35a','c4281b','a59b8f','b39fe1'];

	pass


def plot_raster_lbl_binary_windows(out_path,P,chip_size=512,stride=512):
	'''
	Plot raster with windows marked. Black and white classes.

	Parameters
	----------
	out_path:
	P: chipping.Product() class.
	chip_size:

	'''
	H = P.dw_reader.height #original raster
	W = P.dw_reader.width
	h = P.dw_borders['bottom'] + 1 - P.dw_borders['top'] #trimmed raster (workable area)
	w = P.dw_borders['right'] + 1 - P.dw_borders['left']
	w_h = h - (h % chip_size) + 1 #raster area overlapping windows
	w_w = w - (w % chip_size) + 1


	#Get lines -- red FF0000, yellow FFFF00, green 00FF00
	# yellow_line    = np.ones((3,CHIP_SIZE))
	# yellow_line[0] = 65535
	# yellow_line[1] = 65535
	# yellow_line[2] = 0
	green_line     = np.ones((3,chip_size))
	green_line[0]  = 0
	green_line[1]  = 65535
	green_line[2]  = 0	
	red_line       = np.ones((3,chip_size))
	red_line[0]    = 65535
	red_line[1]    = 0
	red_line[2]    = 0

	kwargs = P.dw_reader.meta.copy()
	kwargs.update({'height':H,'width':W,'count':3,'compress':'lzw'})
	out_ptr = rio.open(out_path,'w',**kwargs)

	dw_windows = chipping.get_windows_strided(P.dw_borders,chip_size,stride=chip_size)

	for _,w in dw_windows:
		#read
		arr        = P.dw_reader.read(1,window=w)
		white_mask = arr == 1 #water
		gray_mask  = arr == 0 #nodata
		black_mask = ~(white_mask | gray_mask) #land

		#change
		arr[white_mask] = 255
		arr[black_mask] = 0
		arr[gray_mask]  = 128
		arr_3d = np.repeat(arr[np.newaxis,:,:],repeats=3,axis=0)

		#write
		out_win = Window(w.col_off-P.dw_borders['left'],w.row_off-dw_borders['top'],chip_size,chip_size)
		out_ptr.write(arr_3d,window=out_win)

	out_ptr.close()
	print("IMAGE WRITTEN TO: %s" % out_path)


def plot_raster_rgb_windows():
	pass


def plot_raster_rgb_bounded():
	pass

####################################################################################################
# POLYGONS
####################################################################################################
def filter_tile_kml(dropped=False):
	'''
	Filter the KML containing the complete list of MGRS tiles used by Sentinel-2.
	The processed result keeps only the tiles present in the data directory (using DynamicWorld
	labels).
	'''

	# CHECK FILE OR .ZIP OF EXISTS 
	kml_path = './kml/S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000_21000101T000000_B00.kml'

	if not os.path.isfile(kml_path):
		#unzip
		zip_path = kml_path[0:-4]+'.zip'
		if not os.path.isfile(zip_path):
			print("No KML or ZIP file found for plotting tiles.")
			return
		else:
			try:
				print(f"Extracting {zip_path}")
				with zipfile.ZipFile(zip_path,'r') as zp:
					zp.extractall('./')
			except:
				print("Could not extract KML from .zip file.")
				return

	#LOAD LIST OF TILES IN DATASET
	# products = glob.glob('*.SAFE',root_dir=DATA_DIR) #Using .SAFE folders
	# tiles = [p.split('-'[5][1:] for p in products)]

	products = glob.glob('*.tif',root_dir=f'{DATA_DIR}/dynamicworld') #Using dynamicworld
	products = [p.split('_')[2][1:6] for p in products]
	tiles_unique,counts = np.unique(products,return_counts=True)
	if dropped:
		tiles = ['11TKE','11SKD']
		out_file_name = 'filtered_tiles_overlapping'
	else:
		tiles = list(tiles_unique)
		out_file_name = 'filtered_tiles_nonoverlapping'		

	print("RASTERS PER TILE")
	print("-"*80)
	for t,c in zip(tiles_unique,counts):
		print(f"{t} | {c} | " + "*"*(c//2))

	kml_ns = {
		None:"http://www.opengis.net/kml/2.2",
		'gx':"http://www.opengis.net/kml/ext/2.2",
		'kml':"http://www.opengis.net/kml/2.2",
		'atom':"http://www.w3.org/2005/Atom"
	}

	#PARSE ORIGINAL SENTINEL-2 KML
	source_root = LT.parse(kml_path).getroot() #<kml>, source_root[0] is <Document>

	#NAMESPACES GOT SILLY JUST REMOVE THEM
	for e in source_root.iter():
		ns,_ = e.tag.split('}')
		e.tag = _

	source_folder = source_root[0][5] # this is <Folder>
	# source_root[0][6] second folder with info at the end.
	document_keep = [e for e in source_root[0][0:5]] + [source_root[0][6]] #first 4 elems and folder at the end

	# START APPENDING STUFF TO NEW KML
	# Append stuff to <FOLDER>
	target_folder = LT.Element("Folder")
	target_folder.text = '\n\t\t'
	target_folder.append(source_folder[0]) #<name>
	target_folder.append(source_folder[0]) #<Snippet maxLines="2">Legend: Single Symbol</Snippet>

	# Append <placemark>'s' to new <Folder>
	for placemark in source_folder:
		if placemark[0].text in tiles:
			target_folder.append(placemark) # Append <Placemark> to new <Folder>
			target_folder[-1][2].text = LT.CDATA(placemark[2].text) # <Add ![CDATA[]] tah in <description>

	# Append new <Folder> to new <DOCUMENT>
	target_document = LT.Element("Document")
	for e in document_keep[0:5]: #first 4 elements before <folder>
		target_document.append(e)
	target_document.text = '\n\t'
	target_document[0].text = 'filtered'
	target_document.append(target_folder) #append <folder> with placemarks

	# Append document to <XML>
	target_root = LT.Element("kml",nsmap=kml_ns)
	target_root.append(target_document)
	target_tree = LT.ElementTree(target_root)	
	target_tree.write(f"./kml/{out_file_name}.kml",encoding='UTF-8',xml_declaration=True,pretty_print=True)
	print(f"KML file written to {out_file_name}.kml")


def plot_map_tile_polygons():
	'''
	Plot a polygons of the United States along with the polygons of the 
	tiles in the dataset.
	'''

	# file paths
	US_SHP_PATH   = "./kml/cb_2024_us_state_500k/cb_2024_us_state_500k.shp"
	S2_KML_PATH_1 = "./kml/filtered_tiles_nonoverlapping.kml"
	S2_KML_PATH_2 = "./kml/filtered_tiles_overlapping.kml"
	WATERB_PATH   = "./kml/sites/sites.shp"
	OUT_PATH      = "./fig/tiles.png"
	common_crs = "EPSG:5070"

	# 1. LOAD FILES
	#---------------
	# 1.1 LOAD US STATES
	states      = gpd.read_file(US_SHP_PATH)
	territories = ['PR','AS','VI','MP','GU','AK','HI']
	contiguous  = states[~states['STUSPS'].isin(territories)]

	# 1.2 LOAD SENTINEL TILES USED
	supported_drivers['LIBKML'] = 'rw' 	# Enable KML driver in fiona
	tiles_gut = gpd.read_file(S2_KML_PATH_1, driver='KML') # read
	tiles_gut['geometry'] = tiles_gut.geometry.apply(lambda x: x.geoms[0]) #POLYGON in GEOMETRYCOLLECTION
	tiles_bad = gpd.read_file(S2_KML_PATH_2, driver='KML') # read
	tiles_bad['geometry'] = tiles_bad.geometry.apply(lambda x: x.geoms[0]) #POLYGON in GEOMETRYCOLLECTION

	# 1.3 LOAD WATERBODY POLYGONS
	water = gpd.read_file(WATERB_PATH)

	# 2. PROJECT TO COMMON CRS
	#-------------------------
	contiguous = contiguous.to_crs(common_crs) #<---- break
	tiles_gut  = tiles_gut.to_crs(common_crs)
	tiles_bad  = tiles_bad.to_crs(common_crs)
	water      = water.to_crs(common_crs)

	# 3. PLOT LAYERS
	# --------------
	fig, ax = plt.subplots(1,1,figsize=(24,20))

	# contiguous.plot(ax=ax,color='white',alpha=1.0,edgecolor='black',linewidth=0.2)
	# water.plot(ax=ax,color='#88D4E9',alpha=1.0,edgecolor='blue',linewidth=0.05)
	# tiles_gut.plot(ax=ax,facecolor='none',alpha=1.0,edgecolor='red',linewidth=1.0)
	# tiles_bad.plot(ax=ax,color='none',alpha=1.0,edgecolor='blue',linewidth=1.0)

	contiguous.plot(ax=ax,color='white',alpha=1.0,edgecolor='black',linewidth=0.2)
	water.plot(ax=ax,color='#88D4E9',alpha=1.0,edgecolor='blue',linewidth=0.1)
	tiles_gut.plot(ax=ax,color='red',alpha=0.3,edgecolor='red',linewidth=1.0)
	tiles_bad.plot(ax=ax,color='blue',alpha=0.3,edgecolor='blue',linewidth=1.5)

	# zoom in
	xmin, ymin, xmax, ymax = contiguous.total_bounds
	print(f"X:{xmin}--{xmax} | Y: {ymin}--{ymax}")

	x_range = xmax - xmin
	y_range = ymax - ymin
	xmax = x_range*0.35 + xmin #~1/3 of US in plot
	xmin = xmin - x_range*0.02
	ymax = ymin + y_range*0.85
	ymin = ymin + y_range*0.30

	ax.set_xlim(xmin,xmax)
	ax.set_ylim(ymin,ymax)

	ax.set_title("Sentinel-2 (MGRS) Tiles",fontsize=24)
	ax.set_axis_off()
	plt.tight_layout()
	plt.savefig(OUT_PATH)

####################################################################################################
# TILE HISTOGRAMS
####################################################################################################
def tile_date_distribution():
	chip_files_tr = glob("*_B0X.tif",root_dir=f"{chip_dir}/training")
	chip_files_va = glob("*_B0X.tif",root_dir=f"{chip_dir}/validation")
	chip_files_te = glob("*_B0X.tif",root_dir=f"{chip_dir}/testing")	

	chip_files = chip_files_tr + chip_files_va + chip_files_te

	pass


def tile_month_distribution():
	pass


####################################################################################################
# OTHER STATISTICS
####################################################################################################
def calculate_band_normal_parameters():
	'''
	Get and store the bands' mean and standard deviation of the pixel values for the entire dataset
	of chips (a la imagenet).
	'''
	training_chips = glob.glob("*_B0X.tif",root_dir=f"{DATA_DIR}/training")

	r_avg_sum = 0.0
	g_avg_sum = 0.0
	b_avg_sum = 0.0
	n_avg_sum = 0.0

	r_std_sum = 0.0
	g_std_sum = 0.0
	b_std_sum = 0.0
	n_std_sum = 0.0

	total_pixels = 0.0

	for chip in training_chips:
		reader = rio.open(chip,'r')

		r_min_i,r_max_i,r_avg_i,r_std_i = reader.statistics(bidx=1,approx=False,clear_cache=False)
		g_min_i,g_max_i,g_avg_i,g_std_i = reader.statistics(bidx=2,approx=False,clear_cache=False)
		b_min_i,b_max_i,b_avg_i,b_std_i = reader.statistics(bidx=3,approx=False,clear_cache=False)
		n_min_i,n_max_i,n_avg_i,n_std_i = reader.statistics(bidx=4,approx=False,clear_cache=False)

		chip_size = reader.height * reader.width
		total_pixels += chip_size

		r_avg_sum += r_avg_i * chip_size
		g_avg_sum += g_avg_i * chip_size
		b_avg_sum += b_avg_i * chip_size
		n_avg_sum += n_avg_i * chip_size

		r_std_sum += r_std_i**2 * chip_size #store sum of squared-differences
		g_std_sum += g_std_i**2 * chip_size
		b_std_sum += b_std_i**2 * chip_size
		n_std_sum += n_std_i**2 * chip_size

	r_avg = r_avg_sum / total_pixels
	g_avg = g_avg_sum / total_pixels
	b_avg = b_avg_sum / total_pixels
	n_avg = n_avg_sum / total_pixels

	r_std = np.sqrt(r_std_sum / total_pixels)
	g_std = np.sqrt(g_std_sum / total_pixels)
	b_std = np.sqrt(b_std_sum / total_pixels)
	n_std = np.sqrt(n_std_sum / total_pixels)

	return ([r_avg,g_avg,b_avg,n_avg],[r_std,g_std,b_std,n_std])


def calculate_band_medians():
	'''
	Same as above for medians.
	'''
	pass


def calculate_class_mean_rgb():
	pass

	
def count_classes():
	pass


####################################################################################################
# MAIN
####################################################################################################
def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--data-dir',default='./dat',
		help="Dataset directory")
	parser.add_argument('--chip-dir',  #<--------------------------------- FIX | separate TILES v chips
		help="Chip directory")
	parser.add_argument('--kml',action='store_true',default=False,
		help="Read tile kml and filter.")
	parser.add_argument('--plots',action='store_true',default=False,
		help='Plot a sample of inputs and labels.')
	args = parser.parse_args()
	#check something here..
	return args


if __name__ == '__main__':

	args = parse_args()

	if args.plots:
		# plot_raster_lbl_binary('./fig/raster_lbl_binary.png',sample_product)
		plot_map_tile_polygons()
		sys.exit(0)

	if args.data_dir is None:
		print("No data directory given.")
		sys.exit(1)

	DATA_DIR = args.data_dir

	if not os.path.isdir(DATA_DIR):
		print("Data directory path not found.")
		sys.exit(1)





