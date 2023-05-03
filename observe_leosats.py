# Code to design and submit LEO satellite observations
# Developer: R. Street

# Test example commandlines:
# python observe_leosats.py ISS coj-doma-1m0a 2023-02-18 2023-02-20 0.2 <path_to_credentials_file> nogo


import argparse
import sathub_ephem_utils
import las_cumbres
from astropy.time import Time

def run():
    # Gather user input
    args = get_args()

    # Load information on LCO network facilities, and extract from it
    # the geographic location of the requested telescope
    lco_network = las_cumbres.LasCumbresNetwork()
    tel = lco_network.get_tel(args.tel_code)

    # Load the user's LCO credentials:
    lco_info = las_cumbres.load_lco_info(args.lco_info)

    # Query SatHub's ephemeris service
    ts_start = Time(args.date_start, format='isot', scale='utc')
    ts_stop = Time(args.date_stop, format='isot', scale='utc')
    params = {'name': args.sat_code,
              'latitude': tel.latitude.deg,
              'longitude': tel.longitude.deg,
              'elevation': tel.elevation.value,
              'jdstart': ts_start.jd,
              'jdstop': ts_stop.jd,
              'jdstep': args.jdstep}
    results = sathub_ephem_utils.query(params)

    # Compose observation requests for each entry, using a default strategy:
    nexp = 5
    overhead = 28.0
    expt = 1.0
    bandpass = 'Bessell-R'
    duration = (nexp*(expt + overhead))/(60.0*60.0*24.0)
    dt = duration / 2.0
    for entry in results:
        ts_obs = Time(entry['JULIAN_DATE'], format='jd', scale='utc')
        obs_start = ts_obs - dt
        obs_stop = ts_obs + dt
        obs_params = {
                      "group_id": args.sat_code+'_'+args.date_start,
                      "submitter": lco_info['submitter'],
                      "proposal_id": lco_info['proposal_id'],
                      "observation_type": "NORMAL",  # NEEDS TO CHANGE
                      "target_name": args.sat_code,
                      "target_type": "ICRS",
                      "ra": entry['RIGHT_ASCENSION-DEG'],
                      "dec": entry['DECLINATION-DEG'],
                      "max_airmass": 1.6,
                      "min_lunar_distance": 10.0,
                      "max_lunar_phase": 1.0,
                      "tel_code": args.tel_code,
                      "exposure_counts": [nexp],
                      "exposure_times": [expt],
                      "filters": [bandpass],
                      "ipp": 1.05,
                      "tstart": obs_start.datetime,
                      "tend": obs_stop.datetime
                      }

        obs = las_cumbres.LasCumbresObservation(obs_params, lco_network)
        obs.build_obs_request()
        print('Observing request details: ')
        print(obs.request)
        
        if args.submit == 'submit':
            response = las_cumbres.lco_api(obs.request, lco_info)
            print('Submitted observation to LCO Network with response: ')
            print(response)
        else:
            print('Status is '+args.submit+' so no observations submitted')

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('sat_code', type=str,
                    help='sat_code: Name of the satellite in the Celestrak database')
    parser.add_argument("tel_code", type=str,
                    help='tel_code: Code of the Las Cumbres Network facility to observe, e.g. lsc-doma-1m0a')
    parser.add_argument("date_start", type=str,
                    help='date_start: UTC Date to observe in YYYY-mm-ddTHH:MM:SS format')
    parser.add_argument("date_stop", type=str,
                    help='date_stop: UTC Date to observe in YYYY-mm-ddTHH:MM:SS format')
    parser.add_argument("jdstep", type=float,
                    help='jdstep: Float number of day intervals to calculate positions for')
    parser.add_argument("lco_info", type=str,
                    help='lco_info: Path to file containing the users LCO token and proposal ID')
    parser.add_argument("submit", type=str,
                    help='submit: Trigger to submit observations to LCO, either "nogo" or "submit"')

    args = parser.parse_args()

    if args.submit not in ["nogo", "submit"]:
        args.submit = "nogo"

    return args


if __name__ == '__main__':
    run()