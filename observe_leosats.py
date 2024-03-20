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
import json
from os import path, remove
import logging

def run():
    # Gather user input
    args = get_args()
    log = start_log(args.data_dir)
    log.info('Started run of observe_leosats')

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
    log.info('Checking for observing opportunities for '
            +str(len(satcat))+' satellites from '
            +str(len(tel_list))+' telescope sites')

    # Resolve date range to search over
    ts_start = Time(args.date_start, format='isot', scale='utc')
    ts_stop = Time(args.date_stop, format='isot', scale='utc')

    # Load the user's LCO credentials:
    lco_info = las_cumbres.load_lco_info(args.lco_info)

    # Query SatHub's ephemeris service for all requested satellites, for
    # all selected locations:
    log2 = open(path.join(args.data_dir, 'visible_satellites.txt'), 'w')

    visible_results = []

    # Query SatHub for all telescopes in the list
    for tel in tel_list:

        # Query for each requested satellite by name
        for sat_code in satcat:
            params = {
                'api': args.endpoint,
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
                if 'Error' in str(entry):
                    print(entry)
                    exit()

                if entry['ILLUMINATED']:
                    log.info('Satellite illuminated on ' + str(entry['JULIAN_DATE']))
                    if not np.isnan(entry['RIGHT_ASCENSION-DEG']):
                        target = {'RA': entry['RIGHT_ASCENSION-DEG'],
                                  'Dec': entry['DECLINATION-DEG']}
                        (visible, status) = tel.check_visibility(target, entry['JULIAN_DATE'])

                        if visible:
                            entry['VISIBLE'] = tel.tel_code+':True'
                            ntrails += 1
                            log.info(entry)

                            # Compose observation requests for each entry, using a default strategy:
                            if args.submit == 'submit':
                                nexp = 11    # Must be an odd number
                                overhead_per_exp = 4.0
                                initial_overhead = 120.0
                                expt = 8.0
                                bandpass = 'Bessell-V'
                                duration = (nexp*(expt + overhead_per_exp))/(60.0*60.0*24.0)
                                dt = duration / 2.0

                                # Take only the first entry, to avoid oversubmitting:
                                ts_obs = Time(entry['JULIAN_DATE'], format='jd', scale='utc')
                                obs_start = ts_obs - TimeDelta(dt, format='sec') - TimeDelta(initial_overhead, format='sec')
                                obs_stop = obs_start + dt
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
                                              "min_lunar_distance": 20.0,
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
                                log.info('Observing request details: ')
                                log.info(obs.request)

                                if args.submit == 'submit':
                                    response = las_cumbres.lco_api(obs.request, lco_info, 'requestgroups')
                                    log.info('Submitted observation to LCO Network with response: ')
                                    log.info(response)
                            else:
                                log.info('Status is '+args.submit+' so no observations submitted')

                        else:
                            entry['VISIBLE'] = tel.tel_code+':False'
                            log.info('Satellite not visible from site ' + args.tel_code)

                    else:
                        entry['VISIBLE'] = tel.tel_code+':False'
                else:
                    entry['VISIBLE'] = tel.tel_code+':False'
                visible_results.append(entry)

                tout = Time(entry['JULIAN_DATE'], format='jd')
                log2.write(entry['NAME']+' '+tout.isot+' RA='
                        +str(entry['RIGHT_ASCENSION-DEG'])+' Dec='
                        +str(entry['DECLINATION-DEG'])
                        +str(entry['VISIBLE'])+'\n')

            log.info(str(ntrails)+' opportunities to observe '+sat_code+' from '
                +tel.tel_code)

    if len(visible_results) == 0:
        log.write('No observing windows within parameters given\n')
    log2.close()
    close_log(log)

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
    parser.add_argument('endpoint', type=str,
                        help='endpoint: SatChecker API to query, one of {name-jdstep, name, catalog-number-jdstep, catalog-number}')
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
    parser.add_argument("data_dir", type=str,
                    help='data_dir: Path to output logfile')
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


def start_log(log_dir, version=None):
    """Function to initialize a log file for a single stage of pyDANDIA.

    The naming convention for the file is [stage_name].log.

    The new file will automatically overwrite any previously-existing logfile
    for the given reduction.

    This function also configures the log file to provide timestamps for
    all entries.

    Parameters:
        log_dir   string        Directory path
                                log_root_name  Name of the log file
        stage_name  string      Name of the stage to be logged
                                Used as both the log file name and the name
                                of the logger Object
        version   string        [optional] Stage code version string
    Returns:
        log       open logger object
    """

    # Console output not captured, though code remains for testing purposes
    console = False
    stage_name = 'obsleosat'

    if path.isdir(log_dir) == False:
        makedirs(log_dir)

    log_file = path.join(log_dir, stage_name + '.log')
    if path.isfile(log_file) == True:
        remove(log_file)

    # To capture the logging stream from the whole script, create
    # a log instance together with a console handler.
    # Set formatting as appropriate.
    log = logging.getLogger(stage_name)

    if len(log.handlers) == 0:
        log.setLevel(logging.INFO)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)

        if console == True:
            console_handler = logging.StreamHandler()
            console_handler.setLevel(logging.INFO)

        formatter = logging.Formatter(fmt='%(asctime)s %(message)s', \
                                      datefmt='%Y-%m-%dT%H:%M:%S')
        file_handler.setFormatter(formatter)

        if console == True:
            console_handler.setFormatter(formatter)

        log.addHandler(file_handler)
        if console == True:
            log.addHandler(console_handler)

    log.info('Started run of ' + stage_name + '\n')
    if version != None:
        log.info('  Software version: ' + version + '\n')

    return log
def close_log(log):
    """Function to cleanly shutdown logging functions with a final timestamped
    entry.
    Parameters:
        log     logger Object
    Returns:
        None
    """

    log.info('Processing complete\n')
    logging.shutdown()


if __name__ == '__main__':
    run()
