# Las Cumbres Observatory Module

# Functions to submit observation requests to the Las Cumbres Observatory APIs

from astropy.coordinates import Angle, SkyCoord
from astropy import units as u

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
        else:
            raise IOError('Unrecognized instrument type, '+str(self.imager))

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
            'coj-clma-2m0a': [Angle('31d16m23.4sS'), Angle('149d4m13.0sE'), 1111.8*u.m, 'fs01'],
            'elp-doma-1m0a': [Angle('30d40m47.53sN'), Angle('104d0m54.63sW'), 2010.0*u.m, 'fa05'],
            'elp-domb-1m0a': [Angle('30d40m48.00sN'), Angle('104d0m55.74sW'), 2029.4*u.m, 'fa07'],
            'lsc-doma-1m0a': [Angle('30d10m2.58sS'), Angle('70d48m17.24sW'), 2201.0*u.m, 'fa15'],
            'lsc-domb-1m0a': [Angle('30d10m2.39sS'), Angle('70d48m16.78sW'), 2201.0*u.m, 'fa04'],
            'lsc-domc-1m0a': [Angle('30d10m2.81sS'), Angle('70d48m16.85sW'), 2201.0*u.m, 'fa03'],
            'cpt-doma-1m0a': [Angle('32d22m50.0sS'), Angle('20d48m36.65sE'), 1807.0*u.m, 'fa14'],
            'cpt-domb-1m0a': [Angle('32d22m50.0sS'), Angle('20d48m36.13sE'), 1807.0*u.m, 'fa01'],
            'cpt-domc-1m0a': [Angle('32d22m50.38sS'), Angle('20d48m36.39sE'), 1807.0*u.m, 'fa06'],
            'coj-doma-1m0a': [Angle('31d16m22.56sS'), Angle('149d4m14.33sE'), 1168.0*u.m, 'fa12'],
            'coj-domb-1m0a': [Angle('31d16m22.89sS'), Angle('149d4m14.75sE'), 1168.0*u.m, 'fa19'],
            'tfn-doma-1m0a': [Anglre('28d18m1.56sN'), Angle('16d30m41.82sE'), 2406.0*u.m, 'fa20'],
            'tfn-domb-1m0a': [Angle('28d18m1.8720sN'), Angle('16d30m41.4360sE'), 2406.0*u.m, 'fa11']
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
            self.telescopes[tel_code] = tel


class LasCumbresObservation():

    def __init__(self, params=None):
        # Credentials
        self.proposal_id = None

        # Facility
        self.facility = None

        # Target
        self.name = None
        self.target_type = 'ICRS'
        self.ra = None
        self.dec = None
        self.proper_motion_ra = None
        self.proper_motion_dec = None
        self.parallax = None
        self.epoch = 2000

        # Observation parameters
        self.exposure_sec = None
        self.nexposures = None
        self.filter = None
        self.ipp = None
        self.obs_type = 'NORMAL'

        # Constraints
        self.max_airmass = None
        self.min_lunar_distance = None
        self.max_lunar_phase = None

        # Scheduling
        self.tstart = None
        self.tend = None

        if params:
            for key, value in params.items():
                if hasattr(self, key):
                    setattr(self, key, value)

    def build_target_dict(self):
        # Check target coordinates are in decimal degrees:
        if type(self.ra) == type(1.0):
            s = SkyCoord(self.ra, self.dec, frame='icrs', unit=(u.deg, u.deg))
        else:
            s = SkyCoord(self.ra, self.dec, frame='icrs', unit=(u.hourangle, u.deg))

        target =   {
                    'name': str(self.name),
                    'type': self.target_type,
                    'ra': s.ra.deg,
                    'dec': s.dec.deg,
                    'proper_motion_ra': 0,
                    'proper_motion_dec': 0,
                    'parallax': 0,
                    'epoch': 2000,
                    }

        return target

    def build_constraints_dict(self):

        constraints = {}
        for key in ['max_airmass', 'min_lunar_distance', 'max_lunar_phase']:
            if getattr(self,key):
                constraints[key] = float(getattr(self, key))

        return constraints

    def build_instrument_configs(self):
        """Function to compose the instrument configuration dictionary for a
        set of exposures"""

        def parse_filter(f):
            filters = { 'SDSS-g': 'gp', 'SDSS-r': 'rp', 'SDSS-i': 'ip',
                       'Bessell-B': 'B', 'Bessell-V': 'V', 'Bessell-R': 'R',
                       'Cousins-Ic': 'I', 'Pan-STARRS-Z': 'zs'
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
                            'exposure_count': int(self.exposure_counts),
                            'exposure_time': float(self.exposure_times),
                            'mode': 'full_frame',
                            'rotator_mode': '',
                            'extra_params': {
                                'offset_ra': 0,
                                'offset_dec': 0,
                                'defocus': 0
                            },
                            'optical_elements': {
                                'filter': parse_filter(self.filter)
                            }
                        }
                    ],
                    'acquisition_config': {
                        'mode': 'OFF',
                        'extra_params': {}
                    },
                    'guiding_config': {
                        'mode': 'ON',
                        'optional': true,
                        'extra_params': {}
                    },
                    'target': target,
                    'constraints': constraints
                }
            config_list.append(config)

        return config_list

    def build_obs_request(self):

        if self.group_id == None:
            self.get_group_id()

        request_group = {
                         'name': self.group_id,
                         'proposal': self.proposal_id,
                         'ipp_value': float(self.ipp),
                         'operator': 'SINGLE',
                         'observation_type': self.obs_type
                         }

        target = self.build_target_dict()
        location = self.facility.build_location_dict()
        constraints = self.build_constraints_dict()
        inst_config_list = self.build_instrument_configs()
        windows = {'start': self.tstart.strftime("%Y-%m-%dT%H:%M:%S"),
                   'end': self.tend.strftime("%Y-%m-%dT%H:%M:%S")}

        request_group['requests'] = [{
                                    'acceptability_threshold': 90,
                                    'configuration_repeats': 1,
                                    'optimization_type': 'TIME',
                                    'configurations': inst_config_list,
                                    'windows': windows,
                                    'location': location
                                    }]

        self.request = request_group

def lco_api(self,request,credentials):
        """Function to communicate with various APIs of the LCO network.
        ur should be a user request in the form of a Python dictionary,
        while end_point is the URL string which
        should be concatenated to the observe portal path to complete the URL.
        Accepted end_points are:
            "userrequests"
        Accepted methods are:
            POST
        """
        PORTAL_URL = 'https://observe.lco.global/api'

        jur = json.dumps(request)

        headers = {'Authorization': 'Token ' + credentials['token']}

        if end_point[0:1] == '/':
            end_point = end_point[1:]
        if end_point[-1:] != '/':
            end_point = end_point+'/'
        url = path.join(PORTAL_URL,end_point)

        response = requests.post(url, headers=headers, json=ur).json()
        
        return response
