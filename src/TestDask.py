# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 12:47:57 2015

@author: Ben Hudson
"""
#import pylab as plt
import numpy as np
#import dask
import xray
import inspect
import os
import gdal
import archook #The module which locates arcgis
archook.get_arcpy()
import arcpy # Hopefully will cut out ArcPy out of this at a later date

#Arcpy can now find TAUDEM routines
arcpy.ImportToolbox('C:\Program Files\TauDEM\TauDEM5Arc\TauDEM Tools.tbx')

# -------------------------------------
# INPUTS 
# -------------------------------------

#goes to parentdirectory of project and gets net cdf file
#see http://stackoverflow.com/questions/3718657/how-to-properly-determine-current-script-directory-in-python


filename = inspect.getframeinfo(inspect.currentframe()).filename
dirNC =  os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(filename))))

ncFileName = dirNC +'\\GreenlandBasinsData\\MCdataset-2015-04-27.nc'
ncOutFileName = dirNC +'\\GreenlandBasinsData\\tempHydroPot.nc'


# -------------------------------------
# FUNCTIONS
# -------------------------------------

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

def netCDF_toGeoTiff(src_filename,dst_filename):
    
    import gdal
    
    #Open existing dataset
    src_ds = gdal.Open( src_filename )
    
    #Open output format driver, see gdal_translate --formats for list
    format = "GTiff"
    driver = gdal.GetDriverByName( format )
    
    #Output to new format
    dst_ds = driver.CreateCopy( dst_filename, src_ds, 0 )
    
    #Properly close the datasets to flush to disk
    dst_ds = None
    src_ds = None

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
    

# opening the dataset and setting up chunking scheme 
ds = xray.open_dataset(ncFileName,chunks={'x': 1000, 'y': 1000})

#trying out data array is this a necessary unpacking step ? 
mask = ds['mask']

# 2 is only land ice
simpleMask = mask.where((mask == 2))

# Split DataSet into more useful ? dataArrays
surface = ds['surface']
thickness = ds['thickness']
bed = ds['bed']
errbed = ds['errbed'] 

# START of monte carlo approach/loop
# alter the bed 

newBed = alterBedDEM(bed,errbed)
    
# now calculate hydro potential 
hydroPot = surface + (.1 * newBed)

hydroPotDS = hydroPot.to_dataset()
hydroPotDS.to_netcdf(ncOutFileName)


#
#
#for i in xrange(50):
#
#    print i
#    
#    ts = time.time()
#    st = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H%M%S_')
#    
#    filePrefix = glacierName+st
#    # main code 
#    #for i in xrange(0,10): 
#    time1 = time.time()    
#    newBed = alterBedDEM(iceBed,iceErrBed)
#    
#    # now calculate hydro potential 
#    hydroPot = iceSurface + (.1 * newBed)
#    
#    #plt.imshow(newBed, vmin=-2000,vmax=2500)
#    #plt.show()
#    
#    # take numpy array to TauDEM suited tif
#    numpyGeoTiff(hydroPot,fileProcessingFolder,filePrefix,surface_tiff)
#    
#    # calculate basins
#    TauDEM_basins(fileProcessingFolder,filePrefix)
#    
#    # turn into Shapefile
#    #arcpy.WaterShedGridToShapefile(fileProcessingFolder+filePrefix+"_basin.tif",fileProcessingFolder+glacierName+str(i+4)+".shp")
#    
#    time2 = time.time()
#    
#    print "time to run:", time2-time1, "seconds"
#    print "or time to run:", (time2-time1)/60, "minutes"
