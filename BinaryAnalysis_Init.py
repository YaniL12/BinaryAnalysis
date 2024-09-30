import os
working_directory = '/avatar/yanilach/PhD-Home/binaries_galah-main/spectrum_analysis/BinaryAnalysis'
os.chdir(working_directory)

# Basic packages
import numpy as np
import time
from pathlib import Path
import logging
import importlib
import sys
from datetime import datetime
import subprocess
import pandas as pd
import json
from astropy.io import fits

# Scipy
import scipy
from scipy.optimize import curve_fit
from scipy import signal

from astropy.table import Table
import AnalysisFunctions as af
import multiprocessing
from multiprocessing.pool import ThreadPool as Pool
# import mysql.connector

sys.path.append(os.path.join(working_directory, 'utils'))
import AstroPandas as ap

import stellarmodel
from stellarmodel import StellarModel

# Accepts mass, log(age), metallicity. Outputs Teff, logg, and log(L) bolometric (flux)
isochrone_table = Table.read(working_directory +  '/assets/parsec_isochrones_logt_8p00_0p01_10p17_mh_m2p75_0p25_m0p75_mh_m0p60_0p10_0p70_GaiaEDR3_2MASS.fits')
isochrone_interpolator = af.load_isochrones()

tracker_path = '/home/yanilach/public_html/avatar-tracker/'

file_lock = multiprocessing.Lock()


def edit_tracker(key, vals):
    # Step 1: Load existing data from JSON file (if it exists)
    try:
        with open("AnalysisTracker.json", "r") as f:
            data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        data = {}  # If the file doesn't exist or is empty, start with an empty dictionary

    # Step 2: 
    data[key] = vals

    # Step 3: Write the updated data back to the JSON file
    with file_lock:
        with open("AnalysisTracker.json", "w") as f:
            json.dump(data, f, indent=4)  # Write JSON with pretty formatting
        
        with open(tracker_path + "AnalysisTracker.json", "w") as f:
            json.dump(data, f, indent=4)  # Write JSON with pretty formatting


def update_tracker(object_ids, val=0, err=None):
    # Step 1: Load existing data from JSON file (if it exists)
    try:
        with open("AnalysisTracker.json", "r") as f:
            data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        data = {}  # If the file doesn't exist or is empty, start with an empty dictionary

    # Ensure the 'objects' key exists in the data
    if 'objects' not in data:
        data['objects'] = {}

    # Step 2: Update or add new entries with current date and time
    for s_id in object_ids:
        # Get the current date and time
        current_time = datetime.now().isoformat()  # e.g., '2024-09-11T14:23:45.123456'
        
        # Update the value and add the timestamp. Overides existing data.
        if val == 0:
            data['objects'][str(s_id)] = {
                'status': val
            }
        elif val == 1:
            # Assume the entry already exists
            data['objects'][str(s_id)]['status'] = val
            data['objects'][str(s_id)]['timestart'] = current_time
        else:
            # Update existing entry with new timestop
            # Assume the entry already exists
            data['objects'][str(s_id)]['status'] = val
            data['objects'][str(s_id)]['timestop'] = current_time

            if err:
                data['objects'][str(s_id)]['error'] = err

    # Step 3: Write the updated data back to the JSON file
    with file_lock:
        with open("AnalysisTracker.json", "w") as f:
            json.dump(data, f, indent=4)  # Write JSON with pretty formatting
        
        with open(tracker_path + "AnalysisTracker.json", "w") as f:
            json.dump(data, f, indent=4)  # Write JSON with pretty formatting



def run_script(args):
    object_id, tmass_id, ages, masses, m_hs = args

    # try:
    #     model, residual = fit_model(args)
    #     update_tracker([object_id], 1)

    #     print("Done")
    #     params = model.get_params(values_only=True)
    #     params = ', '.join(map(str, params))
    #     # print(params)

    #     with file_lock:
    #         with open("fit_results.txt", "a") as f:
    #             f.write(f"{object_id}, {residual}, {params}\n")
                

    # except Exception as e:
    #     print(f"Script failed for object_id {object_id}.")
    #     print("Error message:", e)
    #     update_tracker([object_id], -1)


    ### We HAVE to run this as an external script because of memory allocation issues!

    # Modify the command to run your script with the object ID argument
    print("Beginning script for object_id", object_id)
    command = ["python", "BinaryAnalysis.py", str(object_id), str(tmass_id), str(ages), str(masses), str(m_hs)]

    # Status codes:
    # 0 - Queued, 1 - Processing, 2 - Completed, -1 - Failed
    try:
        # Run the command and check for success
        update_tracker([object_id], 1)
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print(f"Script completed successfully for object_id {object_id}.")

        # Split the output by lines and get the last line
        output_lines = result.stdout.splitlines()
        # Just get the residual and parameters

        final_line = output_lines[-1] if output_lines else "No output received"
        print("Final Output:", final_line)

        #Check if final line contains the word' RV2.
        # TODO make this more robust
        if 'RV2' not in final_line:
            update_tracker([object_id], 2)

            with file_lock:
                with open("fit_results.txt", "a") as f:
                    f.write(f"{object_id}, {final_line}\n")
        else:
            update_tracker([object_id], -1, err=final_line)


    except subprocess.CalledProcessError as e:
        # Handle the error (non-zero return code)
        print(f"Script failed for object_id {object_id} with return code {e.returncode}.")
        print("Error message:", e.stderr)
        update_tracker([object_id], -1, err=e.stderr)

if __name__ == "__main__":

    # Remove pending items from the web interface - starting again

    # mydb = mysql.connector.connect(
    # host="139.99.208.176",
    # user="yanilach_analysis",
    # password="q9jjn=#RQX;Pq:A",
    # database="yanilach_binary_analysis"
    # )

    # mycursor = mydb.cursor()
    # sql = "INSERT INTO processing (SOBJECT_ID, PROCESSING) VALUES (%s, %s)"
    # val = (1704150015010971, 0)
    # mycursor.execute(sql, val)
    # mydb.commit()

   # Table data
    GALAH_DR4_dir = '/avatar/buder/GALAH_DR4/'
    GALAH_DR4 = ap.FitsToDF(GALAH_DR4_dir + "catalogs/galah_dr4_allspec_240207.fits")

    # Accepts mass, log(age), metallicity. Outputs Teff, logg, and log(L) bolometric (flux)
    isochrone_table = Table.read(working_directory +  '/assets/parsec_isochrones_logt_8p00_0p01_10p17_mh_m2p75_0p25_m0p75_mh_m0p60_0p10_0p70_GaiaEDR3_2MASS.fits')
    isochrone_interpolator = af.load_isochrones()

    # Results from CCF
    # results_text = pd.read_csv("/avatar/yanilach/PhD-Home/binaries_galah-main/spectrum_analysis" + "/CCF_results.txt", sep='\t', names=["sobject_id", "no_peaks", "RVs"])
    # results_text['index'] = results_text.index

    # Get all stars in GALAH DR4 where the sobject_id is in obvious_binaries array
    obvious_binaries = pd.read_csv(working_directory + "obvious_binaries.csv")

    # Check the star has a valid rv_2 value
    binary_stars = GALAH_DR4[GALAH_DR4['sobject_id'].isin(obvious_binaries['0'].values)]

    # Get the row where the object ID is 131216001101026
    # TESTING
    # binary_stars = GALAH_DR4[GALAH_DR4['sobject_id'] == 131216001101026]

    object_ids = binary_stars['sobject_id'].values
    tmass_ids = binary_stars['tmass_id'].values

     # Params
    age_min = (10**isochrone_table['logAge'].min()) / 1e9
    age_max = (10**isochrone_table['logAge'].max()) / 1e9

    ages = binary_stars['age'].values.clip(age_min, age_max)
    masses = binary_stars['mass'].values
    m_hs = binary_stars['fe_h'].values

    # print(binary_stars['age'].values)
    # print(ages)
    
    current_time = datetime.now().isoformat()  # e.g., '2024-09-11T14:23:45.123456'

    # Move the current tracker to a backup file in a sudirectory and delete the current tracker
    # Check if a tracker file already exists
    if os.path.exists(tracker_path + "AnalysisTracker.json"):
        backup_path = tracker_path + "runs/"
        Path(backup_path).mkdir(parents=True, exist_ok=True)
        os.rename(tracker_path + "AnalysisTracker.json", backup_path + "AnalysisTracker_" + current_time + ".json")

    if os.path.exists("fit_results.txt"):
        backup_path = "previous_fit_results/"
        Path(backup_path).mkdir(parents=True, exist_ok=True)
        os.rename("fit_results.txt", backup_path + "fit_results " + current_time + ".txt")



    val = {
        'timestart': current_time,
        'no_objects': len(object_ids),
    }
    edit_tracker('meta', val)

    # Append results to a file with safe access
    update_tracker(object_ids)


    # # Create a pool of worker processes. Max number here is 24 at home.
    num_cores_os = 30

    with Pool(processes=num_cores_os) as pool:
        # Run the scripts in parallel
        pool.map(run_script, zip(object_ids, tmass_ids, ages, masses, m_hs))

    current_time = datetime.now().isoformat()  # e.g., '2024-09-11T14:23:45.123456'
    val['timestop'] = current_time
    edit_tracker('meta', val)
