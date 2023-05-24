import pytest
from datetime import datetime
from astropy import units as u
from astropy.time import Time
import numpy as np

@pytest.mark.parametrize(
    "test, expected",
    [
        ({
          "group_id": "TEST1",
          "submitter": "rstreet",
          "proposal_id": "LCO2023A-001",
          "observation_type": "NORMAL",
          "target_name": "M31",
          "target_type": "ICRS",
          "ra": 121.174329,
          "dec": -21.573309,
          "max_airmass": 1.6,
          "min_lunar_distance": 30.0,
          "max_lunar_phase": 1.0,
          "tel_code": "ogg-clma-2m0a",
          "exposure_counts": [1],
          "exposure_times": [30.0],
          "filters": ["Bessell-R"],
          "ipp": 1.05,
          "tstart": datetime.strptime("2017-03-22 14:26:08","%Y-%m-%d %H:%M:%S"),
          "tend": datetime.strptime("2017-03-22 14:26:08","%Y-%m-%d %H:%M:%S")
          },
          {
            "submitter": "rstreet",
            "name": "TEST1",
            "observation_type": "NORMAL",
            "operator": "SINGLE",
            "ipp_value": 1.05,
            "proposal": "LCO2023A-001",
            "requests": [
                        {"location": {
                                     'telescope_class': '2m0',
                                     'site': 'ogg',
                                     'enclosure': 'clma'
                                     },
                        "configurations": [{
                                        "type": "EXPOSE",
                                        "instrument_type": "2M0-SCICAM-MUSCAT",
                                        "instrument_configs": [{
                                                                "exposure_time": 30.0,
                                                                "exposure_count": 1,
                                                                "mode": "default",
                                                                "rotator_mode": "",
                                                                "optical_elements": {
                                                                    "filter": "v"
                                                                },
                                                                "extra_params": {
                                                                    "defocus": 0,
                                                                    "offset_ra": 0,
                                                                    "offset_dec": 0
                                                                }
                                                            }],
                                        "acquisition_config": {
                                                                "mode": "OFF",
                                                                "extra_params": {}
                                                            },
                                        "guiding_config": {
                                                            "mode": "ON",
                                                            "optional": "true",
                                                            "extra_params": {}
                                                        },
                                        "constraints": {
                                                        "max_airmass": 1.6,
                                                        "min_lunar_distance": 30.0,
                                                        "max_lunar_phase": 1.0
                                                    },
                                        "target": {
                                                    "type": "ICRS",
                                                    "name": "M31",
                                                    "ra": 121.174329,
                                                    "dec": -21.573309,
                                                    "proper_motion_ra": 0,
                                                    "proper_motion_dec": 0,
                                                    "parallax": 0,
                                                    "epoch": 2000,
                                                    "extra_params": {}
                                                },
                                    }],
                        "windows": [{
                                    "start": "2017-03-22 14:26:08",
                                    "end": "2017-03-22 14:26:08"
                                    }],
                        "observation_note": "",
                        "state": "PENDING",
                        "acceptability_threshold": 90,
                        "configuration_repeats": 1,
                        "optimization_type": "TIME",
                        "extra_params": {}
                      }
            ]
        })
    ])
def test_build_obs_request(test, expected):

    from leosat_obs_tools.las_cumbres import LasCumbresObservation, LasCumbresNetwork

    facilities = LasCumbresNetwork()

    obs = LasCumbresObservation(test, facilities)
    obs.build_obs_request()

    for key, value in expected.items():
        assert(obs.request[key] == expected[key])

@pytest.mark.parametrize(
    "test, expected",
        [
            (
                {'tel_code': 'lsc-doma-1m0a',
                 'RA': 262.5*u.deg,
                 'Dec': -29.0*u.deg,
                 'time_observe': '2023-07-15T04:00:00'},
                 False
            ),
            (
                {'tel_code': 'ogg-clma-2m0a',
                 'RA': 262.5*u.deg,
                 'Dec': -29.0*u.deg,
                 'time_observe': '2023-07-15T18:00:00'},
                 False
            )
        ]
    )
def test_visibility(test, expected):

    from leosat_obs_tools.las_cumbres import LasCumbresNetwork

    facilities = LasCumbresNetwork()

    tel = facilities.get_tel(test['tel_code'])
    target = {'RA': test['RA'], 'Dec': test['Dec']}
    time_observe = Time(test['time_observe'], format='isot', scale='utc')
    (visible, status) = tel.check_visibility(target, time_observe.jd)
    print(visible,status)
    assert(visible == expected)

@pytest.mark.parametrize(
    "test, expected",
        [
            (
                {'tel_code': 'ogg-clma-2m0a',
                 'time_observe': '2023-07-15T18:00:00'},
                 2460142.125
            )
        ]
    )
def test_get_times_twilight(test, expected):

    from leosat_obs_tools.las_cumbres import LasCumbresNetwork

    facilities = LasCumbresNetwork()

    tel = facilities.get_tel(test['tel_code'])
    time_observe = Time(test['time_observe'], format='isot', scale='utc')
    twilight = tel.get_times_twilight(time_observe.jd, time_format='JD')

    assert(twilight[0].jd == expected)
