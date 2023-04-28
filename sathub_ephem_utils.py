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

    return requests.get(query).json()
