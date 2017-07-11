#!/usr/bin/env python
# coding: utf-8



    #---------------------------------------------------------------------------------------
 #
#        UTILS functions
#
#---------------------------------------------------------------------------------------

def vectorDerivative(v):
    import pandas as pd
    import numpy as np
    import os
    import ginfizz_config            
    
    dv = {}
    for i in range(ginfizz_config.AcqNb):
        # print mocodf['x'][i]
        if i== 0:
            dv[i]= 0
        elif i== ginfizz_config.AcqNb-1:
            dv[i]= v[i]-v[i-1]
        else:
            dv[i]= v[i]-v[i-1]
            #print 'derivative' + str(i)
            #print  v[i]
    return v

def plusDerivative(df):
    import pandas as pd
    import numpy as np
    import os
    import ginfizz_config 
    
    lg = len(df.columns.values)
    dg = df
    for j in list(df.columns.values):
        vprime = vectorDerivative(df[j])
        dg[lg+j]=vprime
    return dg

def plusSquare(df):
    import pandas as pd
    import numpy as np
    import os
    import ginfizz_config
    
    lg = len(df.columns.values)
    ds = df
    for j in range(6):
        print j
        vs = df[j]**2
        ds[lg+j]=vs
    return ds   

# end utils  