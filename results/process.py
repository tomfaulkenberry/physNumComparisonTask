#!/usr/bin/env python

print "\nProcessing..."

import os 
import glob

try:
    import numpy as np # Numeric calculation
    import pandas as pd # General purpose data analysis library
except:
    raise Exception("\
Whoops, you're missing some of the dependencies you need to run this script.\n\
You need to have numpy and pandas installed.")

this_dir = os.path.abspath('.')
print "Running in %s\n\
Checking for .csv files in %s" % (this_dir, os.path.join(this_dir, 'data'))

datafiles = glob.glob('data/*.csv')
print "%i files found:" % len(datafiles)
print '\n'.join(datafiles)



data = pd.concat(
    [pd.DataFrame(pd.read_csv(datafile)) 
     for datafile in datafiles])

print "Done!\n"

# Save data
data.to_csv('processed.csv', index=False)
print "Summary statistics saved to %s" % os.path.join(this_dir, 'processed.csv')
