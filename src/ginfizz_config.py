#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# this file is intented to share global variables and configuration settings

# global variables coming from XML file
AcqNb = int(240)     # Acquisition number 
TR = 2.0          # TR in seconds   

# configuration settings
hpass =  0.01   # high pass for bandpassing
lpass = 0.1     # low pass for bandpassing

plotThreshold = 0.35    # the threshold to get 2D plots with nilearn