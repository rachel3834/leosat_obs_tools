# Las Cumbres Observatory Module

# Functions to submit observation requests to the Las Cumbres Observatory APIs

from astropy.coordinates import Angle, SkyCoord
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import get_sun
from astropy.coordinates import get_moon
import json
from os import path
import requests
import numpy as np

# Developer: R.A. Street
class Telescope():
    def __init__(self):
        self.site = None
        self.enclosure = None
        self.latitude = None
        self.longitude = None
        self.elevation = None
        self.tel = None
        self.imager = None

    def build_location_dict(self):

        location = {
                    'telescope_class' : str(self.tel).replace('a',''),
                    'site':             str(self.site),
                    'enclosure':      str(self.enclosure)
                    }

        return location

    def get_instrument_type(self):
        if 'fa' in self.imager:
            self.instrument_type = '1M0-SCICAM-SINISTRO'
        elif 'mc' in self.imager:
            self.instrument_type = '2M0-SCICAM-MUSCAT'
        elif 'fs' in self.imager:
            self.instrument_type = '2M0-SCICAM-SPECTRAL'
        elif 'kb' in self.imager:
            self.instrument_type = '0M4-SCICAM-SBIG'
        else:
            raise IOError('Unrecognized instrument type, '+str(self.imager))

    def get_EarthLocation(self):
        self.earthlocation = EarthLocation(lat=self.latitude,
                                           lon=self.longitude,
                                           height=self.elevation)

    def check_visibility(self, target, time_observe,
                         min_alt=30.0*u.deg, time_format='JD'):
        """Method to check whether the given target is visible from this
        facility on the datetime given.

        Inputs:
            target dict {RA: float, Dec: float} Units of decimal degrees
            time_observe  float  if time_format=JD (default)
                          str    YYYY-MM-DDTHH:MM:SS

        Output:
            boolean
        """

        # Convert the pointing coordinates to a SkyCoord
        field = SkyCoord(target['RA'], target['Dec'], unit=(u.deg, u.deg))

        # Convert the time to observe to a Time object
        if time_format == 'JD':
            tobs = Time(time_observe, format='jd')
        else:
            tobs = Time(time_observe, format='isot', scale='utc')

        # Calculate the altitude of the target from this observatory, and
        # check whether the target rises above the minimum for observation
        frame = AltAz(obstime=tobs, location=self.earthlocation)
        altaz = field.transform_to(frame)
        altitude = altaz.alt
        if (altitude <= min_alt):
            return False, 'Target altitude too low'

        else:
            # Calculate the altitude of the Sun for the same time(s), and
            # check whether the Sun is lower than 12deg, meaning that the site
            # is in nighttime.
            sun_altaz = get_sun(tobs).transform_to(frame)
            sun_altitude = sun_altaz.alt
            if (sun_altitude >= 12.0*u.deg):
                return False, 'Target visible in daytime'

            # If both conditions are true, the target is visible
            else:
                return True, 'OK'


class LasCumbresNetwork():
    """Configure data on LCO Facilities.  Currently supports only imaging
    facilities on the 2m and 1m networks"""

    def __init__(self):
        self.get_imaging_facilities()

    def get_imaging_facilities(self):
        """Data on all imaging facilities of the Las Cumbres Observatory.
        Data courtesy of Tim Lister, LCO
        https://github.com/LCOGT/neoexchange/blob/
               71da4c6a23ba24e92398e7022fc83c1f9b83f9a4/neoexchange/
                astrometrics/ephem_subs.py#L1416
        Currently only imagers on the LCO 2m and 1m networks are supported.
        """

        #tel_code : site_code: latitude, longitude, elevation
        facilities = {
            'ogg-clma-2m0a': [Angle('20d42m25.5sN'), Angle('156d15m27.4sW'), 3055.0*u.m, 'mc03'],
            'ogg-clma-0m4b': [Angle('20d42m25.1sN'), Angle('156d15m27.11sW'), 3037.0*u.m, 'kb27'],
            #'ogg-clma-0m4c': [Angle('20d42m25.1sN'), Angle('156d15m27.12sW'), 3037.0*u.m, ''],
            'coj-clma-2m0a': [Angle('31d16m23.4sS'), Angle('149d4m13.0sE'), 1111.8*u.m, 'fs01'],
            'coj-doma-1m0a': [Angle('31d16m22.56sS'), Angle('149d4m14.33sE'), 1168.0*u.m, 'fa12'],
            'coj-domb-1m0a': [Angle('31d16m22.89sS'), Angle('149d4m14.75sE'), 1168.0*u.m, 'fa19'],
            'coj-clma-0m4a': [Angle('31d16m22.38sS'), Angle('149d4m15.05sE'), 1191.0*u.m, 'kb26'],
            'coj-clma-0m4b': [Angle('31d16m22.48sS'), Angle('149d4m14.91sE'), 1191.0*u.m, 'kb83'],
            'elp-doma-1m0a': [Angle('30d40m47.53sN'), Angle('104d0m54.63sW'), 2010.0*u.m, 'fa05'],
            'elp-domb-1m0a': [Angle('30d40m48.00sN'), Angle('104d0m55.74sW'), 2029.4*u.m, 'fa07'],
            #'elp-aqwa-0m4a': [Angle('30d40,48.15sN'), Angle('104d0m54.24sW'), 2027.0*u.m, ''],
            'lsc-doma-1m0a': [Angle('30d10m2.58sS'), Angle('70d48m17.24sW'), 2201.0*u.m, 'fa15'],
            'lsc-domb-1m0a': [Angle('30d10m2.39sS'), Angle('70d48m16.78sW'), 2201.0*u.m, 'fa04'],
            'lsc-domc-1m0a': [Angle('30d10m2.81sS'), Angle('70d48m16.85sW'), 2201.0*u.m, 'fa03'],
            'lsc-aqwa-0m4a': [Angle('30d10m3.79sS'), Angle('70d48m16.88sW'), 2202.5*u.m, 'kb98'],
            'lsc-aqwb-0m4a': [Angle('30d10m3.56sS'), Angle('70d48m16.74sW'), 2202.5*u.m, 'kb85'],
            'cpt-doma-1m0a': [Angle('32d22m50.0sS'), Angle('20d48m36.65sE'), 1807.0*u.m, 'fa14'],
            'cpt-domb-1m0a': [Angle('32d22m50.0sS'), Angle('20d48m36.13sE'), 1807.0*u.m, 'fa01'],
            'cpt-domc-1m0a': [Angle('32d22m50.38sS'), Angle('20d48m36.39sE'), 1807.0*u.m, 'fa06'],
            #'cpt-aqwa-0m4a': [Angle('32d22m50.25sS'), Angle('20d48m35.54sE'), 1804.0*u.m, ''],
            'tfn-doma-1m0a': [Angle('28d18m1.56sN'), Angle('16d30m41.82sE'), 2406.0*u.m, 'fa20'],
            'tfn-domb-1m0a': [Angle('28d18m1.8720sN'), Angle('16d30m41.4360sE'), 2406.0*u.m, 'fa11'],
            #'tfn-aqwa-0m4a': [Angle('28d18m1.11sN'), Angle('16d30m42.13sE'), 2390.0*u.m, 'kb95'],
            'tfn-aqwa-0m4a': [Angle('28d18m1.11sN'), Angle('16d30m42.21sE'), 2390.0*u.m, 'kb24'],
            }

        self.telescopes = {}
        for tel_code, tel_data in facilities.items():
            tel = Telescope()
            (tel.site, tel.enclosure, tel.tel) = tel_code.split('-')
            tel.latitude = tel_data[0]
            tel.longitude = tel_data[1]
            tel.elevation = tel_data[2]
            tel.imager = tel_data[3]
            tel.get_instrument_type()
            tel.get_EarthLocation()
            self.telescopes[tel_code] = tel

    def get_tel(self, tel_code):
        if tel_code in list(self.telescopes.keys()):
            return self.telescopes[tel_code]
        else:
            raise IOError('Error: '+tel_code
                    +' is not a recognized facility in the Las Cumbres Network')

class LasCumbresObservation():

    def __init__(self, params=None, facilities=None):
        self.group_id = None

        # Credentials
        self.submitter = None
        self.proposal_id = None

        # Facility
        self.facility = None

        # Target
        self.target_name = None
        self.target_type = 'ICRS'
        self.ra = None
        self.dec = None
        self.proper_motion_ra = None
        self.proper_motion_dec = None
        self.parallax = None
        self.epoch = 2000

        # Observation parameters
        self.exposure_times = []
        self.exposure_counts = []
        self.filters = []
        self.ipp = None
        self.obs_type = 'NORMAL'

        # Constraints
        self.max_airmass = None
        self.min_lunar_distance = None
        self.max_lunar_phase = None

        # Scheduling
        self.tstart = None
        self.tend = None

        # Assign parameters, if provided
        if params:
            for key, value in params.items():
                if hasattr(self, key):
                    if key in ['exposure_times', 'exposure_counts', 'filters']:
                        data = getattr(self, key)
                        for entry in value:
                            data.append(entry)
                        setattr(self, key, data)
                    else:
                        setattr(self, key, value)

        # Assign facility, if parameters given:
        if facilities and params:
            if params['tel_code'] in facilities.telescopes.keys():
                self.facility = facilities.telescopes[params['tel_code']]

    def build_target_dict(self):
        # Check target coordinates are in decimal degrees:
        if type(self.ra) == type(1.0):
            s = SkyCoord(self.ra, self.dec, frame='icrs', unit=(u.deg, u.deg))
        else:
            s = SkyCoord(self.ra, self.dec, frame='icrs', unit=(u.hourangle, u.deg))

        target =   {
                    'name': str(self.target_name),
                    'type': self.target_type,
                    'ra': s.ra.deg,
                    'dec': s.dec.deg,
                    'proper_motion_ra': 0,
                    'proper_motion_dec': 0,
                    'parallax': 0,
                    'epoch': 2000,
                    'extra_params': {}
                    }

        return target

    def build_constraints_dict(self):

        constraints = {}
        for key in ['max_airmass', 'min_lunar_distance', 'max_lunar_phase']:
            if getattr(self,key):
                constraints[key] = float(getattr(self, key))

        return constraints

    def build_instrument_configs(self, target, constraints):
        """Function to compose the instrument configuration dictionary for a
        set of exposures"""

        def parse_filter(f):
            filters = { 'SDSS-g': 'gp', 'SDSS-r': 'rp', 'SDSS-i': 'ip',
                       'Bessell-B': 'b', 'Bessell-V': 'v', 'Bessell-R': 'r',
                       'Cousins-Ic': 'i', 'Pan-STARRS-Z': 'zs'
                       }
            if f in filters.keys():
                return filters[f]
            else:
                raise ValueError('Unrecognized filter ('+f+') requested')

        config_list = []
        for i in range(0,len(self.exposure_times),1):
            config = {
                    'type': 'EXPOSE',
                    'instrument_type': self.facility.instrument_type,
                    'instrument_configs': [
                        {
                            'exposure_count': int(self.exposure_counts[i]),
                            'exposure_time': float(self.exposure_times[i]),
                            'mode': 'default',
                            'rotator_mode': '',
                            "extra_params": {
                                "defocus": 0,
                                "offset_ra": 0,
                                "offset_dec": 0
                            },
                            'optical_elements': {
                                'filter': parse_filter(self.filters[i])
                            }
                        }
                    ],
                    'acquisition_config': {
                        'mode': 'OFF',
                        'extra_params': {}
                    },
                    'guiding_config': {
                        'mode': 'ON',
                        'optional': 'true',
                        'extra_params': {}
                    },
                    'target': target,
                    'constraints': constraints,
                }
            config_list.append(config)

        return config_list

    def build_obs_request(self):

        request_group = {
                        "submitter": self.submitter,
                        "name": self.group_id,
                        "observation_type": self.obs_type,
                        "operator": "SINGLE",
                        "ipp_value": float(self.ipp),
                        "proposal": self.proposal_id,
                         }

        target = self.build_target_dict()
        location = self.facility.build_location_dict()
        constraints = self.build_constraints_dict()
        inst_config_list = self.build_instrument_configs(target, constraints)
        windows = {'start': self.tstart.strftime("%Y-%m-%d %H:%M:%S"),
                   'end': self.tend.strftime("%Y-%m-%d %H:%M:%S")}

        request_group['requests'] = [{
                                        "location": location,
                                        "configurations": inst_config_list,
                                        "windows": [windows],
                                        "observation_note": "",
                                        "acceptability_threshold": 90,
                                        "configuration_repeats": 1,
                                        "optimization_type": "TIME",
                                    }]

        self.request = request_group

def lco_api(request_group, credentials, end_point):
    """Function to communicate with various APIs of the LCO network.
    ur should be a user request in the form of a Python dictionary,
    while end_point is the URL string which
    should be concatenated to the observe portal path to complete the URL.
    Accepted end_points are:
        "requestgroups"
    Accepted methods are:
        POST
    """
    PORTAL_URL = 'https://observe.lco.global/api'

    jur = json.dumps(request_group)

    print('\n',jur,'\n')

    headers = {'Authorization': 'Token ' + credentials['lco_token']}

    if end_point[0:1] == '/':
        end_point = end_point[1:]
    if end_point[-1:] != '/':
        end_point = end_point+'/'
    url = path.join(PORTAL_URL,end_point)

    response = requests.post(url, headers=headers, json=request_group).json()

    return response

def load_lco_info(file_path):
    """Function to load a user's LCO API token and proposal ID code from a
    local file.  The contents of the file are expected to be JSON, of the form:
    {"submitter": <lco user ID>,
     "proposal_id": <proposal code>,
     "lco_token": <LCO token}
    """

    with open(file_path,'r') as f:
        credentials = json.load(f)

    return credentials
