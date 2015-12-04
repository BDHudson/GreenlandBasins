# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 09:36:22 2015

@author: Ben Hudson
"""
#IMPORT PYTHON PACKAGES/ ARCGIS TOOLBOXES
import pylab as plt
import numpy as np
import time
import datetime
import mpl_toolkits.basemap.pyproj as pyproj
import dask.array as da
import archook #The module which locates arcgis
archook.get_arcpy()
import arcpy # Hopefully will cut out ArcPy out of this at a later date

# projection information
#NSIDC Sea Ice Polar Stereographic North
PS_north = pyproj.Proj("+init=EPSG:3413")

# this is correct for greenland 
#utm22n = pyproj.Proj("+init=EPSG:32622")

# THIS IS WGS84
WGS84 = pyproj.Proj("+init=EPSG:4326")

# SET ENVIRONMENTS FOR ARCPY 
arcpy.env.overwriteOutput = True #this overwrites the old file when re-run
#TauDEM needs tiffs to be not compressed
arcpy.env.compression = "None"
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3413)

#Arcpy can now find TAUDEM routines
arcpy.ImportToolbox('C:\Program Files\TauDEM\TauDEM5Arc\TauDEM Tools.tbx')

# DATA SETS NEEDED 

#Give it bouds to work within (Subsbample bit all Greenland TIFF)
    # This is old, most of UO outlets
#Heimdal # Thule # NW # 
N = 65.5 #75.5#79.5 #72.5 #
S = 62.0 #63.0#76.0 #66.85 #
E = -40.0 #-28.0#-52.0 #-26. #
W = -48.0 #-56.2#-65.70 #-39. #

#outletLat = 72.827
#outletLon = -54.309

outletLat = 62.863 #68.545 #72.846
outletLon = -42.608 #-32.718 #-54.016

glacierName = "All"#"Heimdal"#"NW"#"ThuleArea" #"KangerD" 

# convery to 

E_PS, N_PS = pyproj.transform(WGS84,PS_north,E,N)
W_PS, S_PS = pyproj.transform(WGS84,PS_north,W,S)
outletLon_PS,outletLat_PS = pyproj.transform(WGS84,PS_north,outletLon,outletLat)


# THIS INCLUDES ICE FREE AREAS 

# TIFFS 
surface_tiff_ALL = 'Y:\\Documents\\DATA\\MORLIGHEM_NSIDC\\surface_Layer.tif'
bed_tiff_ALL = 'Y:\\Documents\\DATA\\MORLIGHEM_NSIDC\\bed_Layer.tif'
errbed_tiff_ALL = 'Y:\\Documents\\DATA\\MORLIGHEM_NSIDC\\errbed_Layer.tif'
mask_tiff_ALL = 'Y:\\Documents\\DATA\\MORLIGHEM_NSIDC\\mask_Layer.tif'

## TIFFS - the one was added because needed to write new tiff wiht no compression
#surface_tiff = 'Y:\\Documents\\DATA\\MORLIGHEM_NSIDC\\temp_surface_Layer.tif'
#bed_tiff = 'Y:\\Documents\\DATA\\MORLIGHEM_NSIDC\\temp_bed_Layer.tif'
#errbed_tiff = 'Y:\\Documents\\DATA\\MORLIGHEM_NSIDC\\temp_errbed_Layer.tif'
#mask_tiff = 'Y:\\Documents\\DATA\\MORLIGHEM_NSIDC\\temp_mask_Layer.tif'
#
## NOW CLIP OUT ONLY WHAT YOU NEED
##Clip Raster Dataset by known extent - Left Bottom Right Top
#arcpy.Clip_management(surface_tiff_ALL,str(W_PS)+" "+str(S_PS)+" "+str(E_PS)+" "+str(N_PS),surface_tiff, "#", "#", "NONE")
#arcpy.Clip_management(bed_tiff_ALL,str(W_PS)+" "+str(S_PS)+" "+str(E_PS)+" "+str(N_PS),bed_tiff, "#", "#", "NONE")
#arcpy.Clip_management(errbed_tiff_ALL,str(W_PS)+" "+str(S_PS)+" "+str(E_PS)+" "+str(N_PS),errbed_tiff, "#", "#", "NONE")
#arcpy.Clip_management(mask_tiff_ALL,str(W_PS)+" "+str(S_PS)+" "+str(E_PS)+" "+str(N_PS),mask_tiff, "#", "#", "NONE")

# clip to the region you care about, and work with temporary Tiffs from here

# TURN INTO NUMPY - Dask Arrays
surface_numpy = da.from_array(arcpy.RasterToNumPyArray(surface_tiff_ALL),chunks=(1000,1000))
bed_numpy = da.from_array(arcpy.RasterToNumPyArray(bed_tiff_ALL),chunks=(1000,1000))
errbed_numpy = da.from_array(arcpy.RasterToNumPyArray(errbed_tiff_ALL),chunks=(1000,1000))
mask_numpy = arcpy.RasterToNumPyArray(mask_tiff_ALL,nodata_to_value=0)

# for masked array 0 is ocean, 1 is land 2 is ice. tell it if want ice, land + ice, etc. 
# all values that are ice or land set to one
mask_numpy[mask_numpy >= 1] = 1

mask_numpy = da.from_array(mask_numpy,chunks=(1000,1000))

# Tell it some info about the file strucutre/namiming convetion you want.

#fileProcessingFolder = "W:\\MORLIGHEM_NSIDC\\ThuleArea\\"
#fileProcessingFolder = "W:\\MORLIGHEM_NSIDC\\NW_greenland\\"
fileProcessingFolder = "Y:\\Documents\\DATA\\MORLIGHEM_NSIDC\\All\\"

#set all masks to ice only 

#print "ICE SURFACE MASK IS OFF!"

iceSurface = surface_numpy * mask_numpy
iceBed = bed_numpy *  mask_numpy
iceErrBed = errbed_numpy * mask_numpy


def alterBedDEM(bed_numpy,errbed_numpy):
    
    #loading may not be the correct term 
    # uniform distrubution is maybe too agressive, see what Morgleghem says in paper    
    randomErrorLoading = np.random.uniform(-1.0,1.0,np.shape(errbed_numpy))
    
    #triangle distribution
    #randomErrorLoading = np.random.triangular(-1.0,0,1.0,np.shape(errbed_numpy))
    
    #multiply randomErrorLoading, by the error bed
    changeBedBy = randomErrorLoading * errbed_numpy
    
    # apply the error to alter the bed topography 
    newBed = bed_numpy + changeBedBy

    #change no data values back, they were randomly perturbed also    
    #newBed[bed_numpy == -9999] = -9999
    return newBed
    
def numpyGeoTiff(numpyarray,fileProcessingFolder,filePrefix,originalTiff):
    
    arcpy.env.overwriteOutput = True

    # get info about the raster
    r = arcpy.Raster(originalTiff)
    LL_corner_subset = r.extent.lowerLeft
    x_cell_size = r.meanCellWidth
    y_cell_size = r.meanCellHeight
    stringSpatialRef = r.spatialReference.exporttostring()
    
    # DUMP TO GEOTIFF
    
    outRaster1 = arcpy.NumPyArrayToRaster(numpyarray,LL_corner_subset,x_cell_size,y_cell_size,value_to_nodata=0)
    arcpy.DefineProjection_management(outRaster1,stringSpatialRef)    
    
    arcpy.RasterToOtherFormat_conversion(outRaster1,fileProcessingFolder,"TIFF")
    
    outRaster1.save(fileProcessingFolder+filePrefix+"DEM.tif")

def TauDEM_basins(fileProcessingFolder,filePrefix):
    
    arcpy.CheckOutExtension("Spatial")
    
    # StEP 1 - fill pits 
    arcpy.PitRemove(fileProcessingFolder+filePrefix+"DEM.tif",
                    "8",
                    fileProcessingFolder+filePrefix+"pitlessDEM.tif")   
    
    # STEP 2 - D8 FLOW DIRECTIONS 
    arcpy.D8FlowDir(fileProcessingFolder+filePrefix+"pitlessDEM.tif",
                    "8",
                    fileProcessingFolder+filePrefix+"p.tif",
                    fileProcessingFolder+filePrefix+"sd8.tif")
    
    #STEP 3 - Convert TauDEM directions to Arc directions 
    arcpy.gp.Reclassify_sa(fileProcessingFolder+filePrefix+"p.tif",
                           "Value",
                           "1 1;2 128;3 64;4 32;5 16;6 8;7 4;8 2",
                           fileProcessingFolder+filePrefix+"_Reclass.tif",
                           "DATA")
    
    #STEP 4 - Use arcpy to find basins. 
    arcpy.gp.Basin_sa(fileProcessingFolder+filePrefix+"_Reclass.tif",
                      fileProcessingFolder+filePrefix+"_basin.tif")
    
    # Step 5 - Delete the files you made except for basin. 
    arcpy.Delete_management(fileProcessingFolder+filePrefix+"DEM.tif")
    arcpy.Delete_management(fileProcessingFolder+filePrefix+"pitlessDEM.tif")
    arcpy.Delete_management(fileProcessingFolder+filePrefix+"p.tif")
    arcpy.Delete_management(fileProcessingFolder+filePrefix+"sd8.tif")
    arcpy.Delete_management(fileProcessingFolder+filePrefix+"_Reclass.tif")
    

for i in xrange(1):

    print i
    
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H%M%S_')
    
    filePrefix = glacierName+st
    # main code 
    #for i in xrange(0,10): 
    time1 = time.time()    
    newBed = alterBedDEM(iceBed,iceErrBed)
    
    # now calculate hydro potential 
    hydroPot = iceSurface + (.1 * newBed)
    
    hydroPot.compute()
    
    out = np.array(hydroPot)
    #plt.imshow(newBed, vmin=-2000,vmax=2500)
    #plt.show()
    
    # take numpy array to TauDEM suited tif
    numpyGeoTiff(out,fileProcessingFolder,filePrefix,surface_tiff_ALL)
    
    # calculate basins
    TauDEM_basins(fileProcessingFolder,filePrefix)
    
    # turn into Shapefile
    #arcpy.WaterShedGridToShapefile(fileProcessingFolder+filePrefix+"_basin.tif",fileProcessingFolder+glacierName+str(i+4)+".shp")
    
    time2 = time.time()
    
    print "time to run:", time2-time1, "seconds"
    print "or time to run:", (time2-time1)/60, "minutes"

# now make probability maps 

#ID all basin tifs
#import os
#
#basinRawSum = np.zeros(np.shape(mask_numpy))
#basinRawSum_loop = np.zeros(np.shape(mask_numpy))


## can only be one set of basin runs in 
#for files in os.listdir(fileProcessingFolder):
#    if files.endswith("_basin.tif"):
#        
#        
#        
#        #        
#            
#        # this is hardwired now, but should be fixed later
#        #basinID = basin_numpy[3729,1330] # nuuk
#        #basinID = basin_numpy[1422,1391] # isortoq
#        #basinID = basin_numpy[1540,1289] # watson
#        #basinID = basin_numpy[3605,1424] # watson
#        #basinID = basin_numpy[1676,1359] # umiiviit
#        #basinID = basin_numpy[2045,1143] # sarfartoq
#        
#        basin_numpy = arcpy.RasterToNumPyArray(fileProcessingFolder+files) 
#        
#        # everything else to 
#        #basin_numpy[basin_numpy != basinID] = 0
#        
#        # for basin ids that are allowed 
#        for ri in xrange(-150,150,150): 
#            for rj in xrange(-150,150,150): 
#                outletLon_PS_range = outletLon_PS + ri
#                outletLat_PS_range = outletLat_PS + rj
#                
#                #ptList.append([outletLon_PS_range,outletLat_PS_range])
#                
#                #basinID = arcpy.GetCellValue_management(fileProcessingFolder+files,str(outletLon_PS_range)+" "+str(outletLat_PS_range))           
#                basinID = 16762               
#                print basinID
#                # set basin to 1
#                basinRawSum_loop[basin_numpy == basinID] = 1
#        
#                basinRawSum = basinRawSum + basinRawSum_loop
#        plt.imshow(basinRawSum,vmin=0,vmax=1)
#        plt.show()        
#        #plt.imshow(basin_numpy)        
#
##for ri in xrange(-6000,6000,150): 
##        for rj in xrange(-6000,6000,150): 
##            outletLon_PS_range = outletLon_PS + ri
##            outletLat_PS_range = outletLat_PS + rj
##            
##            ptList.append([outletLon_PS_range,outletLat_PS_range])
##
##pt = arcpy.Point()
##ptGeoms = []
##
##for p in ptList:
##    pt.x = p[0]
##    pt.Y = p[1]
##    ptGeoms.append(arcpy.PointGeometry(pt))
#
##arcpy.CopyFeatures_management(ptGeoms, r"C:\Temp\test2.shp")
#
#basinPercent = basinRawSum / np.max(basinRawSum)
#
## sae
#numpyGeoTiff(basinPercent,fileProcessingFolder,glacierName,surface_tiff)
#    
