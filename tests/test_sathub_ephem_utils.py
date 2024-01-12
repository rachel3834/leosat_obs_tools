import pytest

@pytest.mark.parametrize(
    "test, expected",
    [
        (
            {
                'api': 'name',
                'name': 'STARLINK-1600',
                'latitude': 40.1106,
                'longitude': -88.2073,
                'elevation': 222,
                'jd': 2460000.1
            },
            [{
                "ALTITUDE-DEG": -83.92701127488,
                "AZIMUTH-DEG": 74.23644169397,
                "DDEC-DEG_PER_SEC": -0.02013597761,
                "DECLINATION-DEG": -38.21698520948,
                "DRA_COSDEC-DEG_PER_SEC": 0.0273538424,
                "ILLUMINATED": 'true',
                "JULIAN_DATE": 2460000.1,
                "NAME": "STARLINK-1600",
                "PHASE_ANGLE-DEG": 75.47849642996,
                "RANGE-KM": 13236.231719560885,
                "RANGE_RATE-KM_PER_SEC": 0.656362193304,
                "RIGHT_ASCENSION-DEG": 94.32142232641,
                "TLE-DATE": "2023-09-05 16:20:37"
            }]
        ),
        (
                {
                    'api': 'catalog-number',
                    'catalog': '25544',
                    'latitude': 40.1106,
                    'longitude': -88.2073,
                    'elevation': 222,
                    'jd': 2460000.1
                },
            [{
                    "ALTITUDE-DEG": -59.42992120557,
                    "AZIMUTH-DEG": 288.04620638774,
                    "DDEC-DEG_PER_SEC": 0.02460147584,
                    "DECLINATION-DEG": -25.64785198072,
                    "DRA_COSDEC-DEG_PER_SEC": 0.02499960249,
                    "ILLUMINATED": 'true',
                    "JULIAN_DATE": 2460000.1,
                    "NAME": "ISS (ZARYA)",
                    "PHASE_ANGLE-DEG": 41.69217956408,
                    "RANGE-KM": 11477.324789805665,
                    "RANGE_RATE-KM_PER_SEC": -3.431545486776,
                    "RIGHT_ASCENSION-DEG": 134.21602941437,
                    "TLE-DATE": "2023-09-05 16:21:29"
            }]
            ),
        (
                {
                    'api': 'catalog-number-jdstep',
                    'catalog': '25544',
                    'latitude': 40.1106,
                    'longitude': -88.2073,
                    'elevation': 222,
                    'start_jd': 2460000.1,
                    'stop_jd': 2460000.3,
                    'step_jd': 0.1
                },
            [
            {
                "ALTITUDE-DEG": -59.42992120557,
                "AZIMUTH-DEG": 288.04620638774,
                "DDEC-DEG_PER_SEC": 0.02460147584,
                "DECLINATION-DEG": -25.64785198072,
                "DRA_COSDEC-DEG_PER_SEC": 0.02499960249,
                "ILLUMINATED": 'true',
                "JULIAN_DATE": 2460000.1,
                "NAME": "ISS (ZARYA)",
                "PHASE_ANGLE-DEG": 41.69217956408,
                "RANGE-KM": 11477.324789805665,
                "RANGE_RATE-KM_PER_SEC": -3.431545486776,
                "RIGHT_ASCENSION-DEG": 134.21602941437,
                "TLE-DATE": "2023-09-05 16:21:29"
            },
            {
                "ALTITUDE-DEG": -22.86735389391,
                "AZIMUTH-DEG": 142.33553116822,
                "DDEC-DEG_PER_SEC": -0.01420767889,
                "DECLINATION-DEG": -54.03105192755,
                "DRA_COSDEC-DEG_PER_SEC": 0.03650863588,
                "ILLUMINATED": 'true',
                "JULIAN_DATE": 2460000.2,
                "NAME": "ISS (ZARYA)",
                "PHASE_ANGLE-DEG": 118.54352293428,
                "RANGE-KM": 5908.636912798003,
                "RANGE_RATE-KM_PER_SEC": 6.290602878885,
                "RIGHT_ASCENSION-DEG": 30.83552022903,
                "TLE-DATE": "2023-09-05 16:21:29"
            }
            ]
        ),
        (
                {
                    'api': 'name-jdstep',
                    'name': 'STARLINK-1600',
                    'latitude': 40.1106,
                    'longitude': -88.2073,
                    'elevation': 222,
                    'start_jd': 2460000.1,
                    'stop_jd': 2460000.3,
                    'step_jd': 0.1
                },
                [
                    {
                        "ALTITUDE-DEG": -83.92701127488,
                        "AZIMUTH-DEG": 74.23644169397,
                        "DDEC-DEG_PER_SEC": -0.02013597761,
                        "DECLINATION-DEG": -38.21698520948,
                        "DRA_COSDEC-DEG_PER_SEC": 0.0273538424,
                        "ILLUMINATED": 'true',
                        "JULIAN_DATE": 2460000.1,
                        "NAME": "STARLINK-1600",
                        "PHASE_ANGLE-DEG": 75.47849642996,
                        "RANGE-KM": 13236.231719560885,
                        "RANGE_RATE-KM_PER_SEC": 0.656362193304,
                        "RIGHT_ASCENSION-DEG": 94.32142232641,
                        "TLE-DATE": "2023-09-05 16:20:37"
                    },
                    {
                        "ALTITUDE-DEG": -11.8036627367,
                        "AZIMUTH-DEG": 282.38507272541,
                        "DDEC-DEG_PER_SEC": 0.05433004435,
                        "DECLINATION-DEG": 1.75807790636,
                        "DRA_COSDEC-DEG_PER_SEC": 0.00760649602,
                        "ILLUMINATED": 'true',
                        "JULIAN_DATE": 2460000.2,
                        "NAME": "STARLINK-1600",
                        "PHASE_ANGLE-DEG": 53.73895247174,
                        "RANGE-KM": 4328.449597815868,
                        "RANGE_RATE-KM_PER_SEC": -6.016772535669,
                        "RIGHT_ASCENSION-DEG": 210.80053185868,
                        "TLE-DATE": "2023-09-05 16:20:37"
                    }
                ]
        )
    ])
def test_query(test, expected):
    """Test that a standard query to the SatHub Ephemeris service produces
    the expected results"""

    from leosat_obs_tools.sathub_ephem_utils import query

    assert(query(test), expected)
