import numpy as np
import pandas as pd
import os
from google.cloud import storage

# Constants
number_of_PTP1B_molecules = 10**9  # 1 billion molecules
time_steps = 10  # 10 steps
chunk_size = 10**6  # Process molecules in chunks to handle large numbers

# Transition probabilities
P_oxidation_Cys215 = 0.1
P_oxidation_other = 0.028
P_reduction_Cys215 = 0.9
P_reduction_other = 0.1286

# Load proteoform data
proteoform_df = pd.read_csv('PTP1B_proteoforms_ordered.csv')
proteoform_df['k_value'] = proteoform_df.iloc[:, 1:].sum(axis=1)

# Group proteoforms by k value
proteoform_library = proteoform_df.groupby('k_value').apply(lambda x: x.iloc[:, 1:-1].values.tolist()).to_dict()

# Initialize distribution
proteoform_distribution = {k: [0] * len(proteoform_library[k]) for k in range(11)}
proteoform_distribution[0][0] = number_of_PTP1B_molecules  # All molecules start in k=0

# Monte Carlo simulation with chunking
for step in range(time_steps):
    new_distribution = {k: [0] * len(proteoform_library[k]) for k in range(11)}

    for k, proteoforms in proteoform_distribution.items():
        for i, count in enumerate(proteoforms):
            if count > 0:
                num_chunks = int(np.ceil(count / chunk_size))
                for _ in range(num_chunks):
                    chunk_count = min(chunk_size, count)
                    count -= chunk_count

                    # Oxidation: Move from k to k+1
                    if k < 10:
                        oxidation_prob = P_oxidation_Cys215 if k == 0 else P_oxidation_other
                        oxidized = np.random.binomial(chunk_count, oxidation_prob)

                        # Determine which proteoform to transition to
                        for j, next_proteoform in enumerate(proteoform_library[k + 1]):
                            if sum(np.array(next_proteoform) - np.array(proteoform_library[k][i])) == 1:
                                new_distribution[k + 1][j] += oxidized
                                break
                        new_distribution[k][i] += chunk_count - oxidized

                    # Reduction: Move from k to k-1
                    if k > 0:
                        reduction_prob = P_reduction_Cys215 if k == 1 else P_reduction_other
                        reduced = np.random.binomial(chunk_count, reduction_prob)

                        # Determine which proteoform to transition to
                        for j, prev_proteoform in enumerate(proteoform_library[k - 1]):
                            if sum(np.array(proteoform_library[k][i]) - np.array(prev_proteoform)) == 1:
                                new_distribution[k - 1][j] += reduced
                                break
                        new_distribution[k][i] += chunk_count - reduced

    proteoform_distribution = new_distribution

    # Add a progress indicator here
    print(f"Step {step + 1}/{time_steps} completed.")

# Final check
total_molecules_final = sum(sum(proteoforms) for proteoforms in proteoform_distribution.values())
print(f"Total molecules at the end of simulation: {total_molecules_final}")
if total_molecules_final == number_of_PTP1B_molecules:
    print("The total number of molecules matches the initial input.")
else:
    print("Warning: The total number of molecules does NOT match the initial input!")

# Generate Excel file
output_df = pd.DataFrame(columns=['Proteoform_ID', 'k_value'] + list(proteoform_df.columns[1:-1]) + ['Molecule_Count'])
for k, proteoforms in proteoform_distribution.items():
    for i, count in enumerate(proteoforms):
        if count > 0:
            proteoform_data = {
                'Proteoform_ID': proteoform_df[(proteoform_df.iloc[:, 1:-1].values.tolist() == proteoform_library[k][i]).all(axis=1)].index[0],
                'k_value': k,
                'Molecule_Count': count
            }
            proteoform_data.update({site: state for site, state in zip(proteoform_df.columns[1:-1], proteoform_library[k][i])})
            output_df = pd.concat([output_df, pd.DataFrame([proteoform_data])], ignore_index=True)

# Save to a local file first
local_file_path = 'PTP1B_proteoform_final_distribution_with_counts.xlsx'
output_df.to_excel(local_file_path, index=False)
print(f"File saved locally at {local_file_path}.")

# Function to upload to Google Cloud Storage
def upload_to_gcs(bucket_name, source_file_name, destination_blob_name):
    """Uploads a file to Google Cloud Storage."""
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(destination_blob_name)
    blob.upload_from_filename(source_file_name)

    print(f"File {source_file_name} uploaded to {destination_blob_name} in bucket {bucket_name}.")

# Use your actual bucket name
gcs_bucket_name = "gcs_bucket_jamesmontecarlo"
gcs_destination_blob_name = "PTP1B_proteoform_final_distribution_with_counts.xlsx"

# Upload the file to GCS
upload_to_gcs(gcs_bucket_name, local_file_path, gcs_destination_blob_name)
