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
import math

# Scipy
import scipy
from scipy.optimize import curve_fit
from scipy import signal


# /pkg/linux/anaconda3/bin/python -m pip install pyswarm <-- Install here if requred. This is the python used by PBS script.
import pyswarm
from pyswarm import pso

from astropy.table import Table
import AnalysisFunctions as af
import multiprocessing
from multiprocessing.pool import ThreadPool as Pool
# import mysql.connector

from filelock import FileLock

sys.path.append(os.path.join(working_directory, 'utils'))
import AstroPandas as ap
import DataFunctions as df

import stellarmodel
from stellarmodel import StellarModel

# Accepts mass, log(age), metallicity. Outputs Teff, logg, and log(L) bolometric (flux)
isochrone_table = Table.read(working_directory +  '/assets/parsec_isochrones_logt_8p00_0p01_10p17_mh_m2p75_0p25_m0p75_mh_m0p60_0p10_0p70_GaiaEDR3_2MASS.fits')
isochrone_interpolator = af.load_isochrones(type='trilinear')

tracker_path = '/home/yanilach/public_html/avatar-tracker/tracking_files/'
tracking_file = ""
fit_results_file = ""

file_lock = multiprocessing.Lock()


# Are we running on the cluster?
cluster = True
if int(sys.argv[1]) == 1:
    cluster = True

def split_workload(data, num_chunks, chunk_index):
    """Split data into chunks for parallel processing."""
    chunk_size = math.ceil(len(data) / num_chunks)
    start = chunk_size * chunk_index
    end = start + chunk_size
    return data[start:end]


def edit_tracker(key, vals):
    # Step 1: Load existing data from JSON file (if it exists)
    try:
        with open(tracker_path + tracking_file, "r") as f:
            data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        data = {}  # If the file doesn't exist or is empty, start with an empty dictionary

    # Step 2: 
    data[key] = vals

    # Step 3: Ensure the directory for the web tracking file exists
    Path(tracker_path).mkdir(parents=True, exist_ok=True)

    # Step 4: Write the updated data back to the JSON file
    # lock = FileLock("AnalysisTracker.json")
    with file_lock:
    # with lock:
        # Local tracking file
        with open("AnalysisTracker.json", "w") as f:
            json.dump(data, f, indent=4)  # Write JSON with pretty formatting
        
        # Web tracking file
        with open(tracker_path + tracking_file, "w") as f:
            json.dump(data, f, indent=4)  # Write JSON with pretty formatting

            os.chmod(tracker_path + tracking_file, 0o644)


def update_tracker(object_ids, val=0, err=None):
    # Step 1: Load existing data from JSON file (if it exists)
    # -1 = Failed, 0 = Queued, 1 = Processing, 2 = Completed
    try:
        with open(tracker_path + tracking_file, "r") as f:
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

        # Check if this object has a timestart
        if 'timestart' not in data['objects'][str(s_id)]: 
            data['objects'][str(s_id)]['timestart'] = current_time

    # Step 3: Write the updated data back to the JSON file
    # lock = FileLock("AnalysisTracker.json")
    with file_lock:

        # Local tracking file
        with open("AnalysisTracker.json", "w") as f:
            json.dump(data, f, indent=4)  # Write JSON with pretty formatting

        # Web tracking file
        with open(tracker_path + tracking_file, "w") as f:
            json.dump(data, f, indent=4)  # Write JSON with pretty formatting
            
            os.chmod(tracker_path + tracking_file, 0o644)

def cleanup_PBS():
    # Define source and backup directories
    source_path = "/avatar/yanilach/"
    backup_path = "/avatar/yanilach/previous_PBS_runs/"
    
    # Create the backup directory if it doesn't exist
    Path(backup_path).mkdir(parents=True, exist_ok=True)
    
    # Iterate through files in the source directory
    for file in os.listdir(source_path):
        # Check if the filename contains '.o' or '.e'
        if '.o' in file or '.e' in file:
            source_file = os.path.join(source_path, file)
            destination_file = os.path.join(backup_path, file)
            
            # Move the file
            try:
                os.rename(source_file, destination_file)
                print(f"Moved {file} to {backup_path}")
            except Exception as e:
                print(f"Error moving {file}: {e}")


def run_script(args):
    object_id, tmass_id, ages, masses, m_hs = args
    ### We HAVE to run this as an external script because of memory allocation issues!

    # Modify the command to run your script with the object ID argument
    print("Beginning script for object_id", object_id)
    print("Arguments:", object_id, tmass_id, ages, masses, m_hs)
    command = ["python3", "BinaryAnalysis.py", str(object_id), str(tmass_id), str(ages), str(masses), str(m_hs)]
    
    

    # Status codes:
    # 0 - Queued, 1 - Processing, 2 - Completed, -1 - Failed
    try:
        # Run the command and check for success
        update_tracker([object_id], 1)
        # Timeout set to 2 hour
        result = subprocess.run(command, check=True, capture_output=True, text=True, timeout=7200)
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

            # lock = FileLock("fit_results.txt")
            with file_lock:
                with open(fit_results_file, "a") as f:
                    f.write(f"{object_id}, {final_line}\n")
        else:
            update_tracker([object_id], -1, err=final_line)


    except subprocess.TimeoutExpired:
        print(f"Script timed out after 60 minutes for object_id {object_id}.")
        update_tracker([object_id], -1, err="Timeout expired after 60 minutes")

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


    if cluster:
        # Get PBS_ARRAY_INDEX and total jobs from PBS script
        job_index = int(os.environ.get("PBS_ARRAY_INDEX", 0)) - 1
        num_jobs = int(os.environ.get("PBS_ARRAY_LENGTH", 1))
        num_jobs = 50

        print(f"Running job {job_index + 1} of {num_jobs}")

        if job_index == 1:
            cleanup_PBS()

   # Table data
    GALAH_DR4_dir = '/avatar/buder/GALAH_DR4/'
    GALAH_DR4 = df.FitsToDF(GALAH_DR4_dir + "catalogs/galah_dr4_allspec_240207.fits")

    # Accepts mass, log(age), metallicity. Outputs Teff, logg, and log(L) bolometric (flux)
    isochrone_table = Table.read(working_directory +  '/assets/parsec_isochrones_logt_8p00_0p01_10p17_mh_m2p75_0p25_m0p75_mh_m0p60_0p10_0p70_GaiaEDR3_2MASS.fits')
    isochrone_interpolator = af.load_isochrones(type='trilinear')

    # Results from CCF
    # results_text = pd.read_csv("/avatar/yanilach/PhD-Home/binaries_galah-main/spectrum_analysis" + "/CCF_results.txt", sep='\t', names=["sobject_id", "no_peaks", "RVs"])
    # results_text['index'] = results_text.index

    # Get all stars in GALAH DR4 where the sobject_id is in obvious_binaries array
    # obvious_binaries = pd.read_csv(working_directory + "obvious_binaries.csv")
    sample_selection = pd.read_csv(working_directory + "/sample_selection.csv")

    # Check the star has a valid rv_2 value
    binary_stars = GALAH_DR4[GALAH_DR4['sobject_id'].isin(sample_selection['sobject_id'].values)]

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

    # print(object_ids)

    # Split data into chunks for this PBS job
    if cluster:
        # Objects, number of total jobs, job index
        object_ids = split_workload(object_ids, num_jobs, job_index)
        tmass_ids = split_workload(tmass_ids, num_jobs, job_index)
        ages = split_workload(ages, num_jobs, job_index)
        masses = split_workload(masses, num_jobs, job_index)
        m_hs = split_workload(m_hs, num_jobs, job_index)

        if len(object_ids) < 1:
            print(f"No data to process for PBS job index {job_index}. Exiting.")
            exit(0)

    current_time = datetime.now().isoformat()  # e.g., '2024-09-11T14:23:45.123456'
    # Get date and hours minutes only
    current_time = current_time.split("T")[0] + "_" + current_time.split("T")[1].split(":")[0]

    # Move the current tracker to a backup file in a sudirectory and delete the current tracker
    # Check if a tracker file already exists
    if cluster:
        tracker_path_suffix = f"_{job_index}"
    else:
        tracker_path_suffix = ""

    tracking_file = "AnalysisTracker" + tracker_path_suffix + ".json"

    print("Checking for file: ", tracker_path + tracking_file)
    print(os.path.exists(tracker_path + tracking_file))
    if os.path.exists(tracker_path + tracking_file):
        backup_path = tracker_path + "previous_runs/" + current_time + "/"
        Path(backup_path).mkdir(parents=True, exist_ok=True)

        # We don't stricly need to move all files in a loop since this script is being called for each job.
        # However we may change the number of jobs in the future so it's better to be safe.
        if cluster:
            for i in range(100):
                if os.path.exists(tracker_path + tracking_file):
                    os.rename(tracker_path + "AnalysisTracker_" + str(job_index) + ".json", backup_path + current_time  + "_" + tracking_file)
            
            if os.path.exists(tracker_path + tracking_file):
                os.rename(tracker_path + tracking_file, backup_path + current_time + "_" + tracking_file)
        else:
            os.rename(tracker_path + tracking_file, backup_path + current_time + "_" + tracking_file)


    fit_results_file = "fit_results" + tracker_path_suffix + ".txt"
    if os.path.exists(fit_results_file):
        backup_path = "previous_fit_results/" + current_time + "/"
        Path(backup_path).mkdir(parents=True, exist_ok=True)

        if cluster:
            for i in range(100):
                if os.path.exists(fit_results_file):
                    os.rename(fit_results_file, backup_path + "fit_results_" + str(job_index) + ".txt")

            if os.path.exists("fit_results.txt"):
                os.rename("fit_results.txt", backup_path + "fit_results " + current_time + ".txt")

        else:
            os.rename("fit_results.txt", backup_path + "fit_results " + current_time + ".txt")

    current_time = datetime.now().isoformat()  # e.g., '2024-09-11T14:23:45.123456'
    val = {
        'timestart':  current_time,
        'no_objects': len(object_ids),
    }

    # Note we assume an existing tracker file here.
    # If it doesn't exist, we will create a new one.
    edit_tracker('meta', val)

    # Append results to a file with safe access
    update_tracker(object_ids)


    # # Create a pool of worker processes. Max number here is 24 at home.
    # Check this matches the PBS script. DO NOT MODIFY unless you know what you are doing.
    num_cores_os = 10

    with Pool(processes=num_cores_os) as pool:
        # Run the scripts in parallel
        pool.map(run_script, zip(object_ids, tmass_ids, ages, masses, m_hs))

    val['timestop'] = current_time
    edit_tracker('meta', val)
