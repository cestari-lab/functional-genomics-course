import pandas as pd
import argparse

# Set up argument parser
def parse_arguments():
    parser = argparse.ArgumentParser(description='Process MS data and PSMs to calculate intensities.')
    parser.add_argument('ms_data_file', type=str, help='Path to the MS data file (ms-data.csv)')
    parser.add_argument('psm_file', type=str, help='Path to the validated PSM file (valid.mokapot.psms.csv)')
    parser.add_argument('output_file', type=str, help='Path to the output file (outputfile.csv)')
    return parser.parse_args()

# Function to process data
def process_data(ms_data_file, psm_file, output_file):
    # Load the original MS data
    input_df = pd.read_csv(ms_data_file)

    # Display the first few rows of the DataFrame for verification
    print("MS Data:")
    print(input_df.head())

    # Process ScanNumbers to extract numeric IDs from the ScanNumber column
    input_df['ScanNumber'] = input_df['ScanNumber'].str.extract(r'scan=(\d+)')[0]  # Extract the number after 'scan='
    input_df['ScanNumber'] = input_df['ScanNumber'].astype(int)  # Convert to integer

    # Display processed input data for verification
    print("Processed MS Data with Numeric ScanNumbers:")
    print(input_df[['ScanNumber']].head())

    # Load the validated PSM list from the CSV file
    psm_df = pd.read_csv(psm_file)

    # Display PSM Data for verification
    print("PSM Data:")
    print(psm_df.head())

    # Extract unique scan numbers from the PSM DataFrame as integers
    scan_numbers = psm_df['scan'].astype(int).unique()  # Ensure the 'scan' column matches your PSM dataset

    # Filter the MS data to only include relevant scans
    filtered_ms_data = input_df[input_df['ScanNumber'].isin(scan_numbers)]

    # Verify the filtered data
    print("Filtered MS Data based on scan numbers:")
    print(filtered_ms_data)

    # If filtered data is empty, handle accordingly
    if filtered_ms_data.empty:
        print("No matching scans found in the MS data.")
    else:
        # Sum the intensity values for each scan
        sum_intensities = filtered_ms_data.groupby('ScanNumber')['TotalIntensity'].sum().reset_index()

        # Rename the intensity column for clarity
        sum_intensities.rename(columns={'TotalIntensity': 'Intensity'}, inplace=True)

        # Merge the summed intensities back into the PSM DataFrame
        merged_data = pd.merge(psm_df, sum_intensities, how='left', left_on='scan', right_on='ScanNumber')

        # Fill NaN values with 0 for total intensity
        merged_data['Intensity'] = merged_data['Intensity'].fillna(0)

        # Save the final dataset to a new CSV file
        merged_data.to_csv(output_file, index=False)

        print(f"Final output dataset saved to {output_file}.")

# Main execution
if __name__ == "__main__":
    args = parse_arguments()  # Parse command line arguments
    process_data(args.ms_data_file, args.psm_file, args.output_file)  # Process the data based on inputs
