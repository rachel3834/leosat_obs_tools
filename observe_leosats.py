# Code to design and submit LEO satellite observations
# Developer: R. Street

# Test example commandlines:
# python observe_leosats.py ISS coj-doma-1m0a 2023-02-18 2023-02-20 0.2 <path_to_credentials_file> nogo


import argparse
import sathub_ephem_utils
import las_cumbres
from astropy.time import Time, TimeDelta
import csv
import numpy as np

def run():
    # Gather user input
    args = get_args()

    # Resolve list of satellites to search for
    if args.sat_code == 'ALL':
        satcat = load_satcat()
    elif 'search_term=' in args.sat_code:
        search_term = args.sat_code.split('=')[-1]
        satcat = load_satcat(search_term)
    else:
        satcat = [args.sat_code]

    # Load information on LCO network facilities, and extract from it
    # the geographic location of the requested telescope
    lco_network = las_cumbres.LasCumbresNetwork()
    if 'aperture_class=' in args.tel_code:
        search_term = args.tel_code.split('=')[-1]
        tel_list = lco_network.get_aperture_class(search_term)
    else:
        tel_list = [lco_network.get_tel(args.tel_code)]
    print('Checking for observing opportunities for '
            +str(len(satcat))+' satellites from '
            +str(len(tel_list))+' telescope sites')

    # Resolve date range to search over
    ts_start = Time(args.date_start, format='isot', scale='utc')
    ts_stop = Time(args.date_stop, format='isot', scale='utc')

    # Load the user's LCO credentials:
    lco_info = las_cumbres.load_lco_info(args.lco_info)

    # Query SatHub's ephemeris service for all requested satellites, for
    # all selected locations:
    log = open('./visible_satellites.txt', 'w')

    visible_results = []

    # Query SatHub for all telescopes in the list
    for tel in tel_list:

        # Query for each requested satellite by name
        for sat_code in satcat:
            params = {
                'api': 'name-jdstep',
                'name': sat_code,
                'latitude': tel.latitude.deg,
                'longitude': tel.longitude.deg,
                'elevation': tel.elevation.value,
                'start_jd': ts_start.jd,
                'stop_jd': ts_stop.jd,
                'step_jd': args.jdstep
            }
            results = sathub_ephem_utils.query(params)

            # Check for visibility
            ntrails = 0
            for entry in results:
                if entry['ILLUMINATED']:
                    print(entry)
                    if not np.isnan(entry['RIGHT_ASCENSION-DEG']):
                        target = {'RA': entry['RIGHT_ASCENSION-DEG'],
                                  'Dec': entry['DECLINATION-DEG']}
                        (visible, status) = tel.check_visibility(target, entry['JULIAN_DATE'])

                        visible = True  ### REMOVE

                        if visible:
                            entry['VISIBLE'] = tel.tel_code+':True'
                            ntrails += 1

                            # Compose observation requests for each entry, using a default strategy:
                            if args.submit == 'submit':
                                nexp = 4
                                overhead_per_exp = 28.0
                                initial_overhead = 120.0
                                expt = 5.0
                                bandpass = 'Bessell-V'
                                duration = (initial_overhead + nexp*(expt + overhead_per_exp))/(60.0*60.0*24.0)
                                dt = duration / 2.0

                                # Take only the first entry, to avoid oversubmitting:
                                ts_obs = Time(entry['JULIAN_DATE'], format='jd', scale='utc')
                                obs_start = ts_obs - dt
                                obs_stop = ts_obs + dt
                                obs_params = {
                                              "group_id": sat_code.replace(' ','')+'_'+args.date_start,
                                              "submitter": lco_info['submitter'],
                                              "proposal_id": lco_info['proposal_id'],
                                              "obs_type": "TIME_CRITICAL",
                                              "target_name": sat_code,
                                              "target_type": "ICRS",
                                              "ra": entry['RIGHT_ASCENSION-DEG'],
                                              "dec": entry['DECLINATION-DEG'],
                                              "max_airmass": 2.0,
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
                                    response = las_cumbres.lco_api(obs.request, lco_info, 'requestgroups')
                                    print('Submitted observation to LCO Network with response: ')
                                    print(response)
                                else:
                                    print('Status is '+args.submit+' so no observations submitted')

                        else:
                            entry['VISIBLE'] = tel.tel_code+':False'
                    else:
                        entry['VISIBLE'] = tel.tel_code+':False'
                else:
                    entry['VISIBLE'] = tel.tel_code+':False'
                visible_results.append(entry)

                tout = Time(entry['JULIAN_DATE'], format='jd')
                log.write(entry['NAME']+' '+tout.isot+' RA='
                        +str(entry['RIGHT_ASCENSION-DEG'])+' Dec='
                        +str(entry['DECLINATION-DEG'])
                        +str(entry['VISIBLE'])+'\n')

            print(str(ntrails)+' opportunities to observe '+sat_code+' from '
                +tel.tel_code)

    if len(visible_results) == 0:
        log.write('No observing windows within parameters given\n')
    log.close()

def output_results_visibility(visible_results):
    """Function to output a concise and user-friendly summary of observing
    opportunities"""

    with open('./visible_satellites.txt', 'w') as file:
        if len(visible_results) == 0:
            file.write('No observing windows within parameters given\n')
        else:
            for result in visible_results:
                ts = Time(result['JULIAN_DATE'], format='jd')
                file.write(result['NAME']+' '+ts.isot+' RA='
                        +str(result['RIGHT_ASCENSION-DEG'])+' Dec='
                        +str(result['DECLINATION-DEG'])
                        +str(result['VISIBLE'])+'\n')

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('sat_code', type=str,
                    help='sat_code: Name of the satellite in the Celestrak database, ALL or search_term=')
    parser.add_argument("tel_code", type=str,
                    help='tel_code: Code of the Las Cumbres Network facility to observe, e.g. lsc-doma-1m0a or aperture_class=, e.g. aperture_class=0m4')
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

def load_satcat(search_term = None):
    satcat = []
    with open('./satcat.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for i,row in enumerate(reader):
            if i==0:
                keys = row[0].split(',')[1:]
            else:
                data = row[0].split(',')[0]
                if data != 'OBJECT' and not search_term:
                    satcat.append(data)
                elif search_term in data:
                    satcat.append(data)

    return satcat

if __name__ == '__main__':
    run()
