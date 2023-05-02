import pytest
from datetime import datetime

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
            "id": "TEST1",
            "submitter": "rstreet",
            "name": "TEST1",
            "observation_type": "NORMAL",
            "operator": "SINGLE",
            "ipp_value": 1.05,
            "state": "PENDING",
            "proposal": "LCO2023A-001",
            "requests": [
                        {"id": "TEST1",
                        "location": {
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
                                                                "mode": "full_frame",
                                                                "rotator_mode": "",
                                                                "optical_elements": {
                                                                    "filter": "R"
                                                                },
                                                                "extra_params": {
                                                                    "defocus": 0.0,
                                                                    "bin_x": 1,
                                                                    "bin_y": 1
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
                                                    "epoch": 2000
                                                },
                                        "extra_params" : {}
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
