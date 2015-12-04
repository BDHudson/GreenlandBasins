# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 15:19:03 2015

@author: Ben Hudson
"""
from __future__ import division

import os
#IMPORT PYTHON PACKAGES/ ARCGIS TOOLBOXES
import pylab as plt
import numpy as np
import arcpy, arcinfo
import mpl_toolkits.basemap.pyproj as pyproj

# projection information
#NSIDC Sea Ice Polar Stereographic North
PS_north = pyproj.Proj("+init=EPSG:3413")

# this is correct for greenland 
#utm22n = pyproj.Proj("+init=EPSG:32622")

# THIS IS WGS84
WGS84 = pyproj.Proj("+init=EPSG:4326")
#WGS84 = pyproj.Proj(proj='latlong')
# SET ENVIRONMENTS FOR ARCPY 
arcpy.env.overwriteOutput = True #this overwrites the old file when re-run
#TauDEM needs tiffs to be not compressed

arcpy.env.compression = "None"
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3413)

#Arcpy can now find TAUDEM routines
arcpy.ImportToolbox('C:\Program Files\TauDEM\TauDEM5Arc\TauDEM Tools.tbx')

# this is not my typical implmentation but was a quick and dirty work around.
# the main code block is a function so that all variables are cleared each loop 

def loopTheCode(loopNumber):
    # DATA SETS NEEDED 
    
    #Give it bouds to work within (Subsbample bit all Greenland TIFF)
    
    #N = 72.5 #75.5
    #S = 66.85 #63.0
    #E = -26. #-28.0
    #W = -39. #-56.2
    
    #outletLat = 72.827
    #outletLon = -54.309
    # RINK #Upernavik
      #[68.592]#[66.342] #[64.311] # [72.856] #[71.733] #[71.47] #[70.368]#[69.173] #,, # ,,68.545,
    
    outletLatList = [[77.498],[77.664],[74.399],[71.725],[68.592],[66.342],[64.234],[72.808],[71.781],[71.47],[70.368],[69.173]]
    # [77.509],[77.656]
      #[-32.943]#[-38.036] #[-49.655] # [-54.016] #[-51.675] #[-51.446] # [-50.616]#[-49.783] # ,,, # ,,-32.718
    
    outletLonList = [[-65.216],[-66.045],[-56.099],[-52.417],[-32.943],[-38.036],[-49.524],[-54.118],[-51.463],[-51.446],[-50.616],[-49.730]]
    #[-65.847],[-66.128]
      #["KangerD"]#["Sermilik"] #["KNS"] #["Upernavik"] #["Rink"] #["KS"]  # ["Store"]#["Jacobshavn"] #,,, #,,"Kangerdlussuaq",
    
    glacierNameList = [["Heilprin"],["Tracy"],["Alison"],["Umiamako"],["KangerD"],["Sermilik"],["KNS"],["Upernavik"],["Rink"],["KS"],["Store"],["Jacobshavn"]]
    
    # there is a better way to do this 
    outletLatList = outletLatList[loopNumber]
    outletLonList = outletLonList[loopNumber]
    glacierNameList = glacierNameList[loopNumber]
    
    fileProcessingFolderList = ["W:\\MORLIGHEM_NSIDC\\ThuleArea\\",
                                "W:\\MORLIGHEM_NSIDC\\ThuleArea\\",
                                "W:\\MORLIGHEM_NSIDC\\NW_greenland\\",
                                "W:\\MORLIGHEM_NSIDC\\NW_greenland\\",
                                "W:\\MORLIGHEM_NSIDC\\KangerD\\",
                                "W:\\MORLIGHEM_NSIDC\\NW_greenland\\",
                                "W:\\MORLIGHEM_NSIDC\\d8_tri_landIce\\",
                                "W:\\MORLIGHEM_NSIDC\\NW_greenland\\",
                                "W:\\MORLIGHEM_NSIDC\\NW_greenland\\",
                                "W:\\MORLIGHEM_NSIDC\\NW_greenland\\",
                                "W:\\MORLIGHEM_NSIDC\\NW_greenland\\",
                                "W:\\MORLIGHEM_NSIDC\\NW_greenland\\"]
    
    fileProcessingFolder = fileProcessingFolderList[loopNumber]
    
#    errbed_tiff = 'W:\\MORLIGHEM_NSIDC\\temp_errbed_Layer.tif'
    mask_tiff = 'W:\\MORLIGHEM_NSIDC\\temp_mask_Layer.tif'
    
    mask_numpy = arcpy.RasterToNumPyArray(mask_tiff,nodata_to_value=0)
    
    
    # for masked array 0 is ocean, 1 is land 2 is ice. tell it if want ice, land + ice, etc. 
    # all values that are ice or land set to one
    mask_numpy[mask_numpy >= 1] = 1
    
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
        
        outRaster1.save(fileProcessingFolder+filePrefix+".tif")
    

    # find first basin file and use it to set numpy arrays 
    
    ticker = 0
    ticker2 = 0 
    
    for files in os.listdir(fileProcessingFolder):
        if files.endswith("_basin.tif"):
            ticker = ticker + 1
            print files
            guideFile = arcpy.RasterToNumPyArray(fileProcessingFolder+files) 
        if ticker == 1:
            break
        
            
                         
    basinRawSum = np.zeros(np.shape(guideFile))
    basinRawSum_loop = np.zeros(np.shape(guideFile))
        
    #del guideFile 
    
    for locIndex in xrange(len(outletLatList)): 
        
        ptList = []
    
        outletLat = outletLatList[locIndex]
        outletLon = outletLonList[locIndex]
        
        glacierName = glacierNameList[locIndex]
       
        print glacierName
        print locIndex
        print outletLat
        print outletLon
        
        outletLon_PS,outletLat_PS = pyproj.transform(WGS84,PS_north,outletLon,outletLat)
        
        print "At start of script, location is:", outletLon_PS,outletLat_PS 
        
        
        
        # can only be one set of basin runs in 
        for files in os.listdir(fileProcessingFolder):
            if files.endswith("_basin.tif"):
            
                tiff_for_geoRef = files
                # this is hardwired now, but should be fixed later
                #basinID = basin_numpy[3729,1330] # nuuk
                #basinID = basin_numpy[1422,1391] # isortoq
                #basinID = basin_numpy[1540,1289] # watson
                #basinID = basin_numpy[3605,1424] # watson
                #basinID = basin_numpy[1676,1359] # umiiviit
                #basinID = basin_numpy[2045,1143] # sarfartoq
                
                basin_numpy = arcpy.RasterToNumPyArray(fileProcessingFolder+files) 
                
                # everything else to 
                #basin_numpy[basin_numpy != basinID] = 0
                
                # for basin ids that are allowed 
                # This defines the search radius for each catchment
                # It is then visually confirmed that this give the right catchments
                
                for ri in xrange(-6000,6001,3000): 
                    for rj in xrange(-6000,6001,3000): 
                        #ri = 0                
                        #rj = 0
                        outletLon_PS_range = outletLon_PS + ri
                        outletLat_PS_range = outletLat_PS + rj
                        #ptGeoms.append([outletLon_PS_range,outletLat_PS_range])
                        
                        ptList.append([outletLon_PS_range,outletLat_PS_range])
                        
                        # use the point(s) to get cell values
                        basinID_resultObject = arcpy.GetCellValue_management(fileProcessingFolder+files,str(outletLon_PS_range)+" "+str(outletLat_PS_range))           
                        
                        if str(basinID_resultObject.getOutput(0)) != 'NoData':                
                        
                            basinID = float(basinID_resultObject.getOutput(0))
                        
                            #basinID = 16762               
                            print basinID
                            print files
                            ticker2 = ticker2 +1 
                            print ticker2
                            # set basin to 1'
                        
                        
                            basinRawSum_loop[basin_numpy == basinID] = 1.0
                            #basinRawSum_loop[basin_numpy != basinID] = 0
                        
                            basinRawSum = basinRawSum + basinRawSum_loop
                            basinRawSum_loop = np.zeros(np.shape(guideFile))
        
                    
#                    # Dump Each step of basin raw sum
#                    originalTiff = fileProcessingFolder+tiff_for_geoRef 
#                    # get info about the raster
#                    r = arcpy.Raster(originalTiff)
#                    LL_corner_subset = r.extent.lowerLeft
#                    x_cell_size = r.meanCellWidth
#                    y_cell_size = r.meanCellHeight
#                    stringSpatialRef = r.spatialReference.exporttostring()
#                    
#                    shpName = files[0:-4]
#                    outRaster5 = arcpy.NumPyArrayToRaster(basinRawSum,LL_corner_subset,x_cell_size,y_cell_size) # ,value_to_nodata=0
#                    arcpy.DefineProjection_management(outRaster5,stringSpatialRef)    
#        
#                    #arcpy.RasterToOtherFormat_conversion(outRaster1,fileProcessingFolder,"TIFF")
#        
#                    outRaster5.save(fileProcessingFolder+glacierName+shpName+"STEP.tif")
                    
        pt = arcpy.Point()
        ptGeoms = []
        
        for p in ptList:
            pt.X = p[0]
            pt.Y = p[1]
            ptGeoms.append(arcpy.PointGeometry(pt))
        
        print "At the end of script, location is:", pt.X,pt.Y
        
        # This create a shapefile of the points that were used to 'grab' each catchment
        arcpy.CopyFeatures_management(ptGeoms,fileProcessingFolder+glacierName+"outlet.shp")
        
        basinPercent = basinRawSum / np.max(basinRawSum)
        
        
        
        #plt.imshow(basinRawSum,vmin=0,vmax=1)
        #plt.show()        
        #
        ## sae
        #numpyGeoTiff(basinPercent,fileProcessingFolder,glacierName,surface_tiff)
        #    
        
        numpyarray = basinPercent
        #fileProcessingFolder
        filePrefix = glacierName
        originalTiff = fileProcessingFolder+tiff_for_geoRef 
            
        arcpy.env.overwriteOutput = True
        arcpy.CheckOutExtension("Spatial")
        
        #Do min Max
        
        minRaster =np.zeros(np.shape(numpyarray))
        maxRaster =np.zeros(np.shape(numpyarray))
        
        maxRaster[numpyarray > 0.0] = 1
        minRaster[numpyarray > 0.95] = 1    
        
        #maxMinRaster = maxRaster + minRaster
        
        # get info about the raster
        r = arcpy.Raster(originalTiff)
        LL_corner_subset = r.extent.lowerLeft
        x_cell_size = r.meanCellWidth
        y_cell_size = r.meanCellHeight
        stringSpatialRef = r.spatialReference.exporttostring()
        
        # DUMP TO GEOTIFF
        # percent -------------------------
        arcpy.env.overwriteOutput = True #this overwrites the old file when re-run
        outRaster1 = arcpy.NumPyArrayToRaster(basinPercent,LL_corner_subset,x_cell_size,y_cell_size) # ,value_to_nodata=0
        arcpy.DefineProjection_management(outRaster1,stringSpatialRef)    
        
        #arcpy.RasterToOtherFormat_conversion(outRaster1,fileProcessingFolder,"TIFF")
        
        outRaster1.save(fileProcessingFolder+filePrefix+"FINAL_V02.tif")
        
        # MIN -------------------------
        
        outRaster2 = arcpy.NumPyArrayToRaster(minRaster,LL_corner_subset,x_cell_size,y_cell_size) # ,value_to_nodata=0
        arcpy.DefineProjection_management(outRaster2,stringSpatialRef)    
        
        #arcpy.RasterToOtherFormat_conversion(outRaster1,fileProcessingFolder,"TIFF")
        
        outRaster2.save(fileProcessingFolder+filePrefix+"min_V02.tif")    
        # The following inputs are layers or table views: "KangerDmin.tif"
        arcpy.gp.ContourList_sa(fileProcessingFolder+filePrefix+"min_V02.tif",fileProcessingFolder+filePrefix+"min_V02.shp","1")
        
        # MAX -------------------------
        
        outRaster3 = arcpy.NumPyArrayToRaster(maxRaster,LL_corner_subset,x_cell_size,y_cell_size) # ,value_to_nodata=0
        arcpy.DefineProjection_management(outRaster3,stringSpatialRef)    
        
        #arcpy.RasterToOtherFormat_conversion(outRaster1,fileProcessingFolder,"TIFF")
        
        outRaster3.save(fileProcessingFolder+filePrefix+"max_V02.tif") 
        arcpy.gp.ContourList_sa(fileProcessingFolder+filePrefix+"max_V02.tif",fileProcessingFolder+filePrefix+"max_V02.shp","1")
        
        #return basinRawSum_loop
       
# now run the loop. doing it in a function so that all variables are cleared each loop
for i in xrange(8,9):
    basinRawSum = loopTheCode(i)       
