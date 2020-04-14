"""
A script that checks in step of 1 minute if a source is inside the defined window during the night for nu alert observations, starting from a given date.

The only thing that has to be changed is ra,dec of the source. The format of the coordinates can be also radians.

"""

from datetime import datetime, timedelta
import argparse
import time
import math
import ephem
import logging
import sys
import re
import numpy as np

from astropy.io import ascii
from astropy.time import Time

# CTA_N site parameters
CTA_N = ephem.Observer()
CTA_N.elevation = 2200
CTA_N.lat = ephem.degrees('28:45:43')
CTA_N.long = ephem.degrees('-17:53:24')
CTA_N.date = ephem.now()

# CTA_S site parameters
CTA_S = ephem.Observer()
CTA_S.elevation = 0
CTA_S.lat = ephem.degrees('00:00:00')
CTA_S.long = ephem.degrees('00:00:00')
CTA_S.date = ephem.now()

def two_coor(value):
    values = value.split()
    if len(values) != 2:
        raise argparse.ArgumentError
    return values

def parse_args_file():
    """
    Parse command line options and arguments.
    """

    parser = argparse.ArgumentParser(description="\nCheck in step of 1 minute if a list of sources is inside the sea window starting from a given date")
    parser.add_argument("file", metavar="File", type=str, help = "File containing a list of coordinates")

    return parser

def is_night(date=None):
    """
    Check if it is nighttime or daytime.
    """
    if date is None:
        date = datetime.utcnow()
    sun = ephem.Sun()
    CTA_N.date = date
    sun.compute(CTA_N)
    return float(sun.alt) < math.radians(-18)

# cut on Moon phase? which value?
def moonphase(date=None):
    """
    Return the phase of the moon at the given time.
    """
    if date is None:
        date = datetime.utcnow()
    moon = ephem.Moon()
    moon.compute(date)
    return moon.phase/100

def time_for_altitude(obj, alt, startvals, date):
    CTA_N.date = date
    def rel_pos(when):
        CTA_N.date = when
        obj.compute(CTA_N)
        return obj.alt - alt
    start = float(ephem.date(startvals))-0.08
    end = float(ephem.date(startvals))+0.08
    start, end = min(start,end), max(start,end)
    return ephem.date(ephem.newton(rel_pos,start,end))


def time_spent_inside_irf(ra=None, dec=None, date=None):
    if date is None:
        date = datetime.utcnow()
    obj = ephem.FixedBody()
    obj.name = "HET"
    obj._ra = ra
    obj._dec = dec
    obj._epoch = '2000'
    obj.compute(CTA_N)
    # here we need to adapt the alt values, I suggest we pass them as arguments depending on the IRF: time_spent_inside_irf(ra=None, dec=None, date=None, irf_alt0=None, irf_alt1=None)
    # the IRF will be set in the function checking_all_days_during_year_when_object_is_visible_inside_irf which calls this one
    alt0= math.radians(0.0)
    alt1= math.radians(-5.)
    startvals=date
    time_entering_irf=time_for_altitude(obj,alt0,startvals, date)
    time_leaving_irf=time_for_altitude(obj,alt1,startvals, date)
    difftime=float(time_leaving_irf-time_entering_irf)*24*60.
    print(' time entering irf {}'.format(time_entering_irf))
    print(' time leaving irf {}'.format(time_leaving_irf))
    return difftime

# we need to check the TS before entering this function, if TS is < 25 then we do not need to calculate anything
# maybe add an IRF identifier in the arguments? checking_all_days_during_year_when_object_is_visible_inside_irf(ra=None, dec=None, date=None, irf=None)
# then we can make: "case irf == XXX " then set limits to az and alt where it is valid and do the calculation of time
def checking_all_days_during_year_when_object_is_visible_inside_irf(ra=None, dec=None, date=None):
    if date is None:
        sys.exit(1)
    ephem.Date = date
    obj = ephem.FixedBody()
    obj.name = "nu_source"
    obj._ra = ra
    obj._dec = dec
    obj._epoch = '2000'
    obj.compute(CTA_N)
    # check for every minute - do we want it? maybe every 0.5 hour is enough? we need a rough estimate, not an actual schedule planner...
    for x in range(1, 5*24*60): #test check for 5 days
        CTA_N.date = date
        CTA_N.date=CTA_N.date + x/(24.0*60)
        obj.compute(CTA_N)
        alt = obj.alt
        az = obj.az
        altrad=repr(alt)
        altdeg=float(altrad)*180./math.pi
        azdeg=float(repr(az))*180./math.pi
        phase = moonphase(CTA_N.date)
        # here we need a case statement for different IRFs zd (alt) and az ranges for CTA_N and CTA_S;
        #BTW: what is the max zenith (and for which az range) on La Palma in the direction of the Roque?!
        # example values for MAGIC Sea window for tau neutrino observations ;)
        # we also need an if-continue for the Moon phase - obviously we won't observe with full moon ;) suggestions on the limits?
        altdegmin = -5
        altdegmax = 0
        azdegmin = 260
        azdegmax = 340
        if altdeg > altdegmin and altdeg < altdegmax and azdeg >azdegmin and azdeg <azdegmax and is_night(CTA_N.date):
            print(CTA_N.date, " JD", ephem.julian_date(CTA_N.date), " ALT [deg]: ", float(altrad)*180./math.pi," AZ [deg]:", float(repr(az))*180./math.pi," Time inside window [0,-5deg] in min:",time_spent_inside_irf(ra,dec,CTA_N.date)," MoonPhase:", phase) # obviously we do not need to print this out :) we just need the time_spent_inside_irf value summed up over the year to compare with other IRFs output
    
    return alt,az

def main(*args):
    """
    Main function.
    """
    parser = parse_args_file()
    flags = parser.parse_args(args)
    print(flags)
    
    file = flags.file

    # need to be adapted for our output alert files
    source_ids  = np.genfromtxt(file, dtype='int', comments="#", usecols=0, missing_values="NA", delimiter="|")
    ras         = np.genfromtxt(file, dtype='str', comments="#", usecols=1, missing_values="NA", delimiter="|")
    decs        = np.genfromtxt(file, dtype='str', comments="#", usecols=2, missing_values="NA", delimiter="|")

    # start time will be fixed
    start_time = 60676.00 # 2025-01-01 00:00:00  ??? someting else???  #= np.genfromtxt(file, dtype='float', comments="#", usecols=3, missing_values="NA", delimiter="|")
        #print(source_ids)
        #print(ras)
        #print(decs)
        #print(start_times)
    for (source_id, ra, dec) in zip(source_ids, ras, decs):
        # we should add a "if TS > 25" loop and select IRFs for which the source is detected only then call checking_all_days_during_year_when_object_is_visible_inside_irf with the selected IRF
        start_time = Time(start_time, format='mjd')
        #print(start_time.iso)
        print("***************************  Source ID {} ***************************".format(source_id))
        checking_all_days_during_year_when_object_is_visible_inside_irf(ephem.hours(ra), ephem.degrees(dec), start_time.iso)

if __name__ == '__main__':
    main(*sys.argv[1:])
