import csv
import os
import sys
from pyteomics import mzml

# Define input variables from command-line arguments
input_mzml_file = sys.argv[1]  # Open mzML file 
output_csv_file = sys.argv[2]   # Save file name (should end with .csv)

# Ensure the output folder exists, or create it if it doesn't
output_folder = os.path.dirname(output_csv_file)  # Get the folder from the output path
#os.makedirs(output_folder, exist_ok=True)

# Open the mzML file for reading
with mzml.read(input_mzml_file) as mzml_file:
    # Create a CSV file for writing
    with open(output_csv_file, mode='w', newline='') as csv_file:  
        csv_writer = csv.writer(csv_file)
        
        # Write the header row
        csv_writer.writerow(['ScanNumber', 'RetentionTime', 'TotalIntensity'])

        # Initialize a dictionary to hold summed intensities
        intensity_sum = {}

        # Iterate through spectra in the mzML file
        for spectrum in mzml_file:
            scan_number = spectrum['id']
            retention_time = spectrum['scanList']['scan'][0]['scan start time']
            mzs = spectrum['m/z array']
            intensities = spectrum['intensity array']

            # Calculate the total intensity for the scan
            total_intensity = sum(intensities)

            # Store summed intensities for this scan ID
            intensity_sum[scan_number] = total_intensity

        # Write the summed intensity data to the CSV file
        for scan_number, total_intensity in intensity_sum.items():
            # Assuming the same retention time applies for each matched scan ID
            csv_writer.writerow([scan_number, retention_time, total_intensity])

print(f"Conversion completed for {input_mzml_file}. Data saved to {output_csv_file}")

