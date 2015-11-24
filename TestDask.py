# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 12:47:57 2015

@author: Ben Hudson
"""
import pylab as plt
import numpy as np
import dask
import xray
import archook #The module which locates arcgis
archook.get_arcpy()
import arcpy # Hopefully will cut out ArcPy out of this at a later date

#Arcpy can now find TAUDEM routines
arcpy.ImportToolbox('C:\Program Files\TauDEM\TauDEM5Arc\TauDEM Tools.tbx')

# FUNCTIONS

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

# this was downloaded from NSIDC 
filePathNC = "Y:\Documents\DATA\MORLIGHEM_NSIDC\MCdataset-2014-11-19_COPY.nc"

# opening the dataset and setting up chunking scheme 
ds = xray.open_dataset(filePathNC,chunks={'x': 1000, 'y': 1000})

#trying out data array is this a necessary unpacking step ? 
mask = ds['mask']

# 2 is only land ice
simpleMask = mask.where((mask ==2))

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
