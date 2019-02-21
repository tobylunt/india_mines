# This file takes a raw KML converted to .txt file from Elise's work
# identifying visible mineral deposits in Google Earth. We parse the
# locations and any other characteristics into a tabular format,
# output CSV, and then scrape image tiles from the location data.

# NOTE! This needs to run on Python 2.7. This is due to peculiarities
# of getting GDAL and osgeo to work properly on CENTOS. Build a conda
# environment for best results.


################
# (1) Preamble #
################

# necessary imports
import pandas as pd
from datetime import date
from bs4 import BeautifulSoup as Soup
from sentinelsat import SentinelAPI, read_geojson, geojson_to_wkt
from pathlib import Path
import os
import shutil
import subprocess
import glob
import zipfile
import re
from joblib import Parallel, delayed
import multiprocessing
import ogr, osr # for spatialreference
from osgeo import gdal # we will use the gdal python library (python 2.7 only)

# prompt user for copernicus credentials
user = raw_input("Please enter copernicus username: ")
password = raw_input("Please enter password: ")

# connect to the API
api = SentinelAPI(user, password, 'https://scihub.copernicus.eu/dhus')

# get user name
myname = subprocess.check_output('whoami')[:-1]

# get scratch into an object
tmppath = os.path.expanduser('/scratch/' + myname)

# define working directory for contained unzipping
unzippath = tmppath + '/unzip_tmp'

# set imagery storage path for master full-size scenes (not chips)
imgpath = tmppath + '/raw_sentinel_imgs'

# make sure we have the necessary subfolders - tif includes R, G, B,
# and IR layers in .tif chips; TC for true color .jpgs. both will be
# named by the mine_id. establish paths for these folders
tifchippath = tmppath + "/chips/tif"
jpgchippath = tmppath + "/chips/jpg"

# define a function for checking directory existence, and clearing if desired
def checkdir(dirstring):
    if os.path.isdir(dirstring):
        files = glob.glob(dirstring + '/*')
        for f in files:
            os.system('rm -rf '+ f)
        print "original contents removed"
    else:
        print "creating directory " + dirstring
        os.makedirs(dirstring)

# run for the two required folders
checkdir(tifchippath)
checkdir(jpgchippath)
checkdir(unzippath)

# don't want to wipe old raw images, unlike the other dirs. so make
# this manually.
if os.path.isdir(imgpath):
    print imgpath + " already exists. continuing"
else:
    os.makedirs(imgpath)

# define a function for translating points across CRSs
def transform_coords(infile, x_in, y_in):
    
    # get WKT projection from input 4-layer tif
    stacked_full = gdal.Open(infile).GetProjection()
    
    # retrieve the input CRS - this will not be constant across India! 
    src = osr.SpatialReference()
    src.ImportFromWkt(stacked_full)
    
    # create a geometry from input coordinates
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(x_in, y_in)
    
    # create coordinate transformation. CRS 4326 is WGS84, which is the
    # CRS we have our coordinates in from google earth. need to convert
    # these into the georeference used by the Sentinel scenes.
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(4326)
    
    # the output spatial reference is pulled from the 4-layer tif
    coordTransform = osr.CoordinateTransformation(inSpatialRef, src)
    
    # transform the point
    point.Transform(coordTransform)
    
    # return new point in correct georeference
    return [point.GetX(), point.GetY()]

# RUN A TEST. This is from a corner of a sample TIF. 
# degrees: 76d51'59.81"E, 15d41'14.02"N
# decimal degrees: 76.866614, 15.687228
# georeferenced coords: 700050.000, 1735220.000
# transform_coords(outfile, 76.866614, 15.687228) - WORKS


############################
# (2) Clean KML to tabular #
############################

# initialize our output dataframe
output = pd.DataFrame()

# get KML into beautiful soup, parsing as XML
filename = os.path.expanduser('~/iec1/minec_new/elise_mines.kml')
with open(filename) as data:
    kml_soup = Soup(data) 

# now we go into the KML and tabulate the HTML bits
placemarks = kml_soup.find_all('placemark')

# iterate all places marked in the KML
for place in placemarks:
    
    # get the single variables for this place that we want to record
    name = place.find_all('name')[0].text
    desc_raw = place.find_all('description')[0].text
    
    # process the description into key-value pairs.
    try:
        desc_raw = desc_raw.rstrip('\n')
        new_row = dict(item.split(':') for item in desc_raw.split('\n'))
    except:
        new_row = {'notes' : desc_raw}
        
    # get GE's lat/lon from the LookAt tag - this is the lat/lon for
    # the VIEW, not the feature itself.
    lookat = place.find_all('lookat')
    
    # there should only be one LookAt per placemark. assert this
    assert len(lookat) == 1
    
    # get lat and lon from lookat
    view_lat = lookat[0].find_all('latitude')[0].text
    view_lon = lookat[0].find_all('longitude')[0].text
    
    # now extract the lat/lon for the pins. these are in nested <point> and <coordinates> tags
    coords = place.find('point').find_all('coordinates')[0].text
    
    # split into a list so we can separate lat and lon. 
    coords = coords.split(',')
    lon = coords[0]
    lat = coords[1]
    
    # add the other variables to our dict for this placemark
    new_row.update({'view_lat': view_lat, 'view_lon': view_lon, 'name': name, 'lat' : lat, 'lon' : lon})
    
    # append this new row to our output dataframe
    output = output.append(new_row, ignore_index=True)

# clean up Type and type columns
output['type'] = output['type'].fillna(output['Type'])
output = output.drop('Type', 1)

# remove any slash characters in the type note
output['type'] = output['type'].str.replace('/','-')

# remove stock google pins (from data entry)
target_strings = ['Eiffel', 'Christ', 'Canyon', 'Sydney', 'Basilica', 'London', 'Titanic', 'Forbidden', 'Fuji', 'Google']
for string in target_strings: 
    output = output[~output.name.str.contains(string)]

# drop mines where no features were found
output = output[pd.notnull(output['type'])]
output = output[~output.type.str.contains('none')]

# drop known small mines
output = output[~output.w.astype(str).str.contains('<')]
output = output[~output.w.astype(str).str.contains('.5')]

# create a unique index for each observation, starting at 1 (string fmt)
output = output.reset_index(drop=True)
output['scrape_id'] = output.index + 1
output['scrape_id'] = output['scrape_id'].astype(str)

# coerce lat and lon datatypes
output['lat'] = output['lat'].astype(float)
output['lon'] = output['lon'].astype(float)

# quick check - note that there are likely duplicates in view_lat/view_lon, but not lat/lon.
#output.info()
#output.describe(include = 'all').transpose()

# write out to CSV
output.to_csv('~/iec1/minec_new/elise_mines.csv', index=False)

# for testing:
#output = pd.DataFrame({'scrape_id':['15'], 'lat':[22.10997928], 'lon':[85.43066325], 'name':['Iron Ore 53 - feature 1'], 'notes':['pretty clear mining activity'],'view_lat':[22.103936828],'view_lon':[5.43165325],'type':['terraced land'],'w':[2.5]})


###########################################
# (3) Assign sentinel scenes to each mine #
###########################################

# initialize a dictionary for matching sentinel product ids to mine ids
mine_sent_ids = {}

# define a function for retrieving the most recent cloud-free sentinel
# image covering a particular point. 
def get_best_scene_pid(row, in_dict):
    
    # extract mine ID
    mine_id = row['scrape_id']

    # get lat and lon into objects
    lon = row['lon']
    lat = row['lat']

    # build your footprint from your point. needs to be in the following format:
    # 'POLYGON((ul, ur, lr, ll, ul))'. this is just a rectangle expanded
    # around the lat/lon point by a fraction of a degree. the sentinelsat API
    # will return the whole scene (raster image) that covers this
    # rectangle (but need a rectangle, can't just use a point)
    # python3 version: footprint = "POLYGON((" + str(round(lon - .05, 8)) + " " + str(round(lat + .05, 8)) + "," + str(round(lon + .1, 8)) + " " + str(round(lat + .1, 8)) + "," + str(round(lon + .1, 8)) + " " + str(round(lat - .1, 8)) + "," + str(round(lon - .1, 8)) + " " + str(round(lat - .1, 8)) + "," + str(round(lon - .05, 8)) + " " + str(round(lat + .05, 8)) + "))"
    footprint = "POLYGON((" + "%.8f" % round(lon - .05, 8) + " " + "%.8f" % round(lat + .05, 8) + "," + "%.8f" % round(lon + .1, 8) + " " + "%.8f" % round(lat + .1, 8) + "," + "%.8f" % round(lon + .1, 8) + " " + "%.8f" % round(lat - .1, 8) + "," + "%.8f" % round(lon - .1, 8) + " " + "%.8f" % round(lat - .1, 8) + "," + "%.8f" % round(lon - .05, 8) + " " + "%.8f" % round(lat + .05, 8) + "))"

    # search by polygon, time, and SciHub query keywords
    products = api.query(footprint,
                         date=('20180101',
                               date(2018, 12, 29)),
                         platformname='Sentinel-2',
                         cloudcoverpercentage=(0, 1))

    # GeoPandas GeoDataFrame with the metadata of the scenes and the footprints as geometries
    scenes = api.to_geodataframe(products)

    # extract date from the summary variable - it's the first thing before
    # the column. clip the leading "date: " and convert to datetime
    scenes['sdate'] = pd.to_datetime(scenes['summary'].str.split(', ', expand=True)[0].str.slice(6,))

    # sort by most recent date and lowest cloud cover percentage
    scenes = scenes.sort_values(['sdate', 'cloudcoverpercentage'], ascending=[False, True])

    # get the product ID for this most desirable scene
    pid = scenes.index[0]
    
    # append the id pair to the input dict
    in_dict[mine_id] = pid

# get rows of output dataframe into a list for passing into the above function
rowlist = []
for index, row in output.iterrows():
    rowlist.append(row)

# execute the function to expand the mine_sent_ids dict to contain all ids
map(lambda x:get_best_scene_pid(x, mine_sent_ids), rowlist)


########################################################
# (4) Get sentinel scenes using sentinelsat python API #
########################################################

# initialize a dictionary for matching sentinel filename stubs to mine ids
sent_mine_keys = {}

# define a fucntion for downloading scenes with basic processing
def download_scene(pid):
    
    # get the filename stub and the TCI fn into objects
    metadata = api.get_product_odata(pid)
    downtitle = metadata['title']
    tci_fn = imgpath + '/' + pid + '_tci.jp2'

    # set the unzipped root directory into an object for ease of use
    unzipped_root = unzippath + '/' + downtitle + '.SAFE'

    # check to see if we've got this Sentinel scene already. put it in a
    # try block for safety. look for the true color image
    try:
        my_abs_path = Path(tci_fn).resolve()
        
    # if we don't have the processed images ready, we need to download and unzip them.
    except OSError:

        print "file does not exist - downloading pid: " + pid
        downinfo = api.download(pid, directory_path = tmppath, checksum=False)
        downpath = downinfo['path']

        # unzip the download to a new tmp dir
        print "unzipping from: " + downpath
        with zipfile.ZipFile(downpath,"r") as zip_ref:
            for info in zip_ref.infolist():
                print info.filename, info.date_time, info.file_size, info.compress_size
            zip_ref.extractall(unzippath)
        print "successfully unzipped the download from " + downpath

        # need to glob for the image path, as there's a variant subfolder that
        # contains images but consistently begins with 'L1C'
        unzipped_img = glob.glob(unzipped_root + '/GRANULE/' + 'L1C*')[0] + '/IMG_DATA'

        # find our target layers using wildcards
        jp2_ir_dl = {'ir' : glob.glob(unzipped_img + '/' + '*_B08.jp2')[0]}
        jp2_r_dl  = {'r'  : glob.glob(unzipped_img + '/' + '*_B04.jp2')[0]}
        jp2_g_dl  = {'g'  : glob.glob(unzipped_img + '/' + '*_B03.jp2')[0]}
        jp2_b_dl  = {'b'  : glob.glob(unzipped_img + '/' + '*_B02.jp2')[0]}
        jp2_tc_dl = {'tci' : glob.glob(unzipped_img + '/' + '*_TCI.jp2')[0]}

        # move these to the temporary unzip location to the raw image path
        for jp2dict in [jp2_ir_dl, jp2_r_dl, jp2_g_dl, jp2_b_dl, jp2_tc_dl]:
            shutil.copyfile(jp2dict.values()[0], imgpath + '/' + pid + '_' + jp2dict.keys()[0] + '.jp2')

        # remove the raw zip
        os.remove(downpath)
        
        # remove the full unzipped path
        os.system('rm -rf ' + unzipped_root)

    # if we do have the processed images, we're good.
    else:
        # inform user that file exists
        print "file exists - continuing"

# execute the fuction for all mines. copernicus data hub only allows
# two concurrent downloads!
pids = mine_sent_ids.values()
num_cores = 2

# use a try/catch in case of unexpected errors 
def dl_recursive():
    try:
        Parallel(n_jobs=num_cores)(delayed(download_scene)(pid) for pid in pids)
    except Exception as e:
        print e
        print "error - trying again"
        dl_recursive()
dl_recursive()

# for testing:
#sent_mine_keys = {15: u'/scratch/lunt/unzip_tmp/31c86ad0-bdbb-440e-bb66-afa566890db1'}


#####################################
# (5) Clip, jitter, and merge bands #
#####################################

# define a function for processing images
def sent_postprocess(mine_id, fn_path):
    
    # extract pid from fn_path
    pid = fn_path.split('/')[-1]
    
    # set file names into objects
    jp2_ir   = fn_path + '_ir.jp2'
    jp2_r    = fn_path + '_r.jp2'
    jp2_b    = fn_path + '_b.jp2'
    jp2_g    = fn_path + '_g.jp2'
    jp2_tci  = fn_path + '_tci.jp2'
    
    # specify the output filename for the 4-layer output tif for the
    # sentinel scene, and the chipped version for this mine
    stacked_out = imgpath + '/' + pid + '.tif'
    stacked_chip = tifchippath + '/' + mine_id + '.tif'
    
    # output fn for TCI jpg
    tci_chunk = imgpath + '/' + mine_id + '_tci.tif'
    
    # consolidate ir, r, g, and b into a 4-layer tif
    # NOTE: don't have a way to get the fourth alpha band in there yet.
    # os.system('/users/tobiaslunt/anaconda3/bin/gdal_merge.py -separate -co PHOTOMETRIC=RGB -o ' + stacked_out + ' ' + jp2_r + ' ' + jp2_b + ' ' + jp2_g)
    os.system('gdal_merge.py -separate -co PHOTOMETRIC=RGB -o ' + stacked_out + ' ' + jp2_r + ' ' + jp2_b + ' ' + jp2_g)
    
    # retrive lat and lon from the dataframe for this mine - mine_id
    # is just the row index plus one to be indexed starting at 1
    # rather than zero.
    lat = output.loc[output['scrape_id'] == mine_id, 'lat'][0]
    lon = output.loc[output['scrape_id'] == mine_id, 'lon'][0]
    
    # chip the 4-layer tif. projwin is two points: upper left; lower right.
    ulx = lon - .05
    uly = lat + .05
    lrx = lon + .05
    lry = lat - .05
    [ulx_trans, uly_trans] = transform_coords(stacked_out, ulx, uly)
    [lrx_trans, lry_trans] = transform_coords(stacked_out, lrx, lry)
    
    # chip the 4-layer tif
    gdal.Translate(stacked_chip, stacked_out, projWin = [ulx_trans, uly_trans, lrx_trans, lry_trans])
    
    # remove the georeferencing info from the chip - just pixels (to
    # match resnet input data format). mogrify is an imagemagick command.
    os.system('mogrify -strip ' + stacked_chip)
    
    # chip the jpg. step 1 is jp2 to tif, then we'll chip the tif with jittering below
    gdal.Translate(tci_chunk, jp2_tci, projWin = [ulx_trans, uly_trans, lrx_trans, lry_trans])
    os.system('mogrify -strip ' + tci_chunk)
    
    # store the center pixel output string from gdalinfo
    centerpoints = os.popen('gdalinfo ' + stacked_chip + ' 2>/dev/null | grep Center').read()
    
    # extract the center x and y pixels. start with x.
    center_x = re.match("^Center\s+\((?P<match>[^,]+)", centerpoints)
    center_x = center_x.group('match')
    center_x = int(float(center_x))
    
    # now y
    center_y = centerpoints.split('(')[1].split(',')[1].split(')')[0]
    center_y = int(float(center_y))
    
    # chips need to be 256*256 centered on the point. we know the
    # center is at the correct location. we will jitter by 100 pixels
    # in each combination of the cardinal directions, for 9 total
    # images. first, create an array of x,y offsets and suffixes
    offsets = {'1' : [-50, -50], '2' : [0, -50], '3' : [50, -50], '4' : [-50, 0], '5' : [0, 0], '6' : [50, 0], '7' : [-50, 50], '8' : [0, 50], '9' : [50, 50]}
    
    # now create the jittered tifs and jpgs
    for key, value in offsets.iteritems():
        
        # set center pixel values
        jittered_x = center_x + value[0]
        jittered_y = center_y + value[1]
        
        # set corner pixel values
        ulx = "%.0f" % (jittered_x - 128)
        uly = "%.0f" % (jittered_y - 128)
        lrx = "%.0f" % (jittered_x + 128)
        lry = "%.0f" % (jittered_y + 128)
        
        # chip the true color tif to jpg
        # NOTE: NOT SURE WHETHER TO -SCALE OR NOT. Choosing to use original scale.
        os.system('gdal_translate -of JPEG --config GDAL_PAM_ENABLED NO -projwin ' + ulx + ' ' + uly + ' ' + lrx + ' ' + lry + ' ' + tci_chunk + ' ' + jpgchippath + '/' + mine_id + '_' + key + '.jpg')
        
        # chip the stacked tif
        os.system('gdal_translate -of Gtiff --config GDAL_PAM_ENABLED NO -projwin ' + ulx + ' ' + uly + ' ' + lrx + ' ' + lry + ' ' + stacked_chip + ' ' + tifchippath + '/' + mine_id + '_' + key + '.tif')
    
    # remove the individual bands now that we have the stacked scene.
    for file in [jp2_r, jp2_ir, jp2_b, jp2_g]:
        os.remove(file)
    
    # remove the master stacked chip now that we've created our jitters, and the master TCI chunk
    os.remove(stacked_chip)
    os.remove(tci_chunk)


# add the associated sentinel filename stub to a dictionary keyed on the mine_id
sent_mine_keys = {}
for mine_id in sorted(mine_sent_ids):
     sent_mine_keys[mine_id] = imgpath + '/' + mine_sent_ids[mine_id]

# 2.7 python iteration syntax - 3 is for key, value in d.items():
for key, value in sent_mine_keys.iteritems():
    sent_postprocess(key, value)



##############################################################
# (6) publish the raw color images to HTML for manual review #
##############################################################

# set the public jpg path into an object
jpg_public = os.path.expanduser('~/public_html/png/mine_chips')

# make a directory to store jpg files in the public html folder
checkdir(jpg_public)

# copy all the jpg chips to the public html directory for assessment
os.system('cp ' + jpgchippath + '/* ' + jpg_public + '/')

# define header and footer html blocks
header = '''
<html>
  <body>
'''

footer = '''
  </body>
</html>
'''

# this defines the block of a 3x3 minigrid for each mine_id -
# centered, plus 8 jitters around it.
s = '''  <h1>Mine id: $0</h1>
	 <table>
           <tr>
          <td align=center>
	  fn: $1.jpg<br>
	  <a href="$11">
	    <img src="$11" width=256>
	  </a>
	</td>
	<td align=center>
	  fn: $2.jpg<br>
	  <a href="$22">
	    <img src="$22" width=256>
	  </a>
	</td>
	<td align=center>
	  fn: $3.jpg<br>
	  <a href="$33">
	    <img src="$33" width=256>
	  </a>
	</td>
      </tr>
       <tr>
	<td align=center>
	  fn: $4.jpg<br>
	  <a href="$44">
	    <img src="$44" width=256>
	  </a>
	</td>
	<td align=center>
	  fn: $5.jpg<br>
	  <a href="$55">
	    <img src="$55" width=256>
	  </a>
	</td>
	<td align=center>
	  fn: $6.jpg<br>
	  <a href="$66">
	    <img src="$66" width=256>
	  </a>
	</td>
      </tr>
       <tr>
	<td align=center>
	  fn: $7.jpg<br>
	  <a href="$77">
	    <img src="$77" width=256>
	  </a>
	</td>
	<td align=center>
	  fn: $8.jpg<br>
	  <a href="$88">
	    <img src="$88" width=256>
	  </a>
	</td>
	<td align=center>
	  fn: $9.jpg<br>
	  <a href="$99">
	    <img src="$99" width=256>
	  </a>
	</td>
      </tr>
    </table>
    <hr>
'''

# define a function for adding images to our html body in groups of
# nine (jitter sets for each mine)
def replace_img(mine_idlist, header):
    stub = 'http://caligari.dartmouth.edu/~' + myname + '/png/mine_chips'
    output = header
    for i in range(0, len(mine_idlist)): # iterate over each mine
        mine_id = mine_idlist[i]
        this_s = s
        for j in range(1,10): # loop over the nine positions for this mine
            j = "%.0f" % j
            this_s = this_s.replace('$' + j + j, stub + '/' + mine_id + '_' + j + '.jpg')
            this_s = this_s.replace('$' + j, mine_id + '_' + j)
        this_s = this_s.replace('$0', mine_id)
        output = output + this_s
    output = output + footer
    return output

# get list of img filenames in this folder
minelist = output['scrape_id'].tolist()

# create our body html to include all mine image chips
html = replace_img(minelist, header)

# write out
with open(jpg_public + '/**mine_chips.html', 'w') as text_file:
    text_file.write(html)

# print out html location
print 'view html at http://caligari.dartmouth.edu/~' + myname + '/png/mine_chips/**mine_chips.html'


