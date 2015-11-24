# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 12:47:57 2015

@author: Ben Hudson
"""

import dask
import xray

filePathNC = "Y:\Documents\DATA\MORLIGHEM_NSIDC\MCdataset-2014-11-19_COPY.nc"


ds = xray.open_dataset(filePathNC)