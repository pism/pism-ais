#!/usr/bin/env python

"""
SCOPE: Fix time_bnds variable in NetCDF datasets for use with PISM (using NCO's ncap2).

USAGE: `python correct_time_bnds.py <ncfile> <year_ref> <year_start> <year_end>`

NOTE: - Depending on file size, this process can take a while.
      - The warning 'WARNING assign(): Var being read and written in ASSIGN time_bnds' might be ignored.
      - <year_ref> can for example be obtained using CDO ("RefTime" in `cdo sinfon <ncfile>` output).

FIXME: - Currently relying on time_bnds variable being present in the file. In future, we should create time_bnds, if not present.
       - Currently assumes a 365-day calendar with time units in days since <year_ref> as well as monthly timesteps of the data.

AUTHOR: julius.garbe@pik
"""

import sys
import subprocess

filename = sys.argv[1]

print "Correcting time_bnds for file:", filename

dayspermonth = [31,28,31,30,31,30,31,31,30,31,30,31] #365_day calendar

year_ref   = int(sys.argv[2])
year_start = int(sys.argv[3])
year_end   = int(sys.argv[4])

time_bnds = [[],[]]

for (y,yearafterref) in enumerate(range(year_start-year_ref,year_end-year_ref,1)): # years after reference year
    
    lower0 = yearafterref*sum(dayspermonth)
    
    for (m,days) in enumerate(dayspermonth):
        
        lower = lower0 if (y==0 and m==0) else time_bnds[1][-1]
        upper = lower+days
        
        time_bnds[0].append(lower)
        time_bnds[1].append(upper)
        
nmonths = len(time_bnds[0])
nyears = nmonths/12

print "Time units reference year: "+str(year_ref)
print "Years to be corrected: "+str(year_start)+"-"+str(year_end)+" ("+str(nyears)+" years / "+str(nmonths)+" months)"

lower_str = ",".join(map(str,time_bnds[0]))
upper_str = ",".join(map(str,time_bnds[1]))
ncap2_str = "ncap2 -O -s 'time_bnds(:,0)={"+lower_str+"};time_bnds(:,1)={"+upper_str+"};' "+filename+" "+filename
#print ncap2_str

print "Fixing time_bnds..."
subprocess.check_call(ncap2_str, shell=True)

