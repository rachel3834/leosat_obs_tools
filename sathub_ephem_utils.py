# sathub_ephem_utils

# Library of functions to query the SatHub Ephemeris service,
# which provides information on the expected positions of satellites
# in Earth orbit.

# Credit: This code is based on the demonstration notebook at
# https://github.com/iausathub/ephemeris_api_demo/blob/main/demo.ipynb
# by Anthony Rihani and Seigfried Eggl.
# Developer: R.A. Street

import requests
from os import path
from astropy.time import Time
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun
from astropy.coordinates import get_moon
from astropy import constants
import csv

def query(params):
    """Function to compose and send a query to SatHub's Ephemeris service

    Input:
    params  dict    Dictionary of required query parameters including:
                    name: string Name of the satellite in the Celestrak database
                    latitude: float  Latitude of observer in decimal degrees
                    longitude: float Longitude of observer in decimal degrees
                    elevation: float Height of observer above mean sealevel in m
                    JD: float        Time of observation
    Returns:
    results json    Dictionary of returned results
    """

    ROOT_URL = 'http://apexgroup.web.illinois.edu/ephemeris/'

    if 'jd' in params.keys():
        query = path.join(ROOT_URL, 'name', params['name'],
                            str(params['latitude']), str(params['longitude']),
                            str(params['elevation']), str(params['jd']))

    elif 'jdstart'in params.keys() and 'jdstop' in params.keys() \
            and 'jdstep' in params.keys():
        query = path.join(ROOT_URL, 'namejdstep', params['name'],
                            str(params['latitude']), str(params['longitude']),
                            str(params['elevation']),
                            str(params['jdstart']), str(params['jdstop']),
                            str(params['jdstep']))

    else:
        raiseIOError('Incorrect parameter set provided for query')

    results = requests.get(query)
    if results.status_code == 200:
        results = results.json()

    return results

def summarize_results(query_results):
    """Function to output a concise and user-friendly summary of observing
    opportunities"""

    if len(query_results) == 0:
        print('No observing windows within parameters given')
    else:
        for result in query_results:
            ts = Time(result['JULIAN_DATE'], format='jd')
            print(result['NAME']+' '+ts.isot+' RA='
                    +str(result['RIGHT_ASCENSION-DEG'])+' Dec='
                    +str(result['DECLINATION-DEG']))

def check_sat_illuminated(query_results, tel):
    """Function to test whether a satellite is likely to be illuminated, based
    on whether it falls in the Earth's shadow.  Based on code from Siegfried Eggl."""

    # Time of observation
    tobs = Time(query_results['JULIAN_DATE'], format='jd')

    # Calculate the satellite position wrt to the geocenter
    # sat = satellite.at(t).position.km
    #frame = GCRS(obstime=tobs)
    #field = SkyCoord(query_results['RIGHT_ASCENSION-DEG'],
    #                 query_results['DECLINATION-DEG'],
    #                 unit=(u.deg, u.deg))
    #sat = field.transform_to(frame)
    sat_range_from_geoc = 6371.146 + tel.elevation/1000.0 + query_results['RANGE-KM']

    #normalized
    satn = sat/np.linalg.norm(sat_range_from_geoc)

    # Calculate the Earth sun vector
    #frame = AltAz(obstime=tobs, location=tel.earthlocation)
    frame = GCRS(obstime=tobs)
    sun_altaz = get_sun(tobs).transform_to(frame)
    sun_altaz.distance*constants.au

    earthsun = earth.at(t).observe(sun).position.km
    earthsunn = earthsun/np.linalg.norm(earthsun)
    #satellite sun vector
    satsun =  sat - earthsun
    satsunn = satsun/np.linalg.norm(satsun)

    #Is the satellite in Earth's Shadow?
    r_parallel = np.dot(sat,earthsunn)*earthsunn
    r_tangential = sat-r_parallel

    illuminated = True
    if(np.linalg.norm(r_tangential)<rearthkm):

        #yes the satellite is in Earth's shadow
        illuminated = False
