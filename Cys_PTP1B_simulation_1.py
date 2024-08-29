import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
import seaborn as sns

# Constants
number_of_PTP1B_molecules = 6022141500000000  # Number of PTP1B molecules
time_steps = 173  # Number of steps
delta_t = 1  # Time per step in seconds

# Manually adjust this probability until the desired 10% is achieved
P_oxidation_Cys215 = 0.1  # Start with this value and adjust as needed
P_oxidation_other = 0.028  # This can remain constant or be adjusted if necessary
P_reduction_Cys215 = 0.9  # Adjust as needed
P_reduction_other = 0.1286  # Adjust as needed

# Load the proteoform data from the CSV file
proteoform_df = pd.read_csv('/content/PTP1B_proteoforms_ordered.csv')

# Calculate the k value (number of oxidized cysteines) for each row
proteoform_df['k_value'] = proteoform_df.iloc[:, 1:].sum(axis=1)  # Sum across all cysteine columns to get k value

# Group the proteoforms by k value to create a library
proteoform_library = proteoform_df.groupby('k_value').apply(lambda x: x.iloc[:, 1:-1].values.tolist()).to_dict()

# Initialize the distribution of PTP1B molecules
# 100% start at k = 0 (fully reduced)
proteoform_distribution = {k: 0 for k in range(11)}
proteoform_distribution[0] = number_of_PTP1B_molecules

# Monte Carlo simulation
for step in range(time_steps):
    new_distribution = {k: 0 for k in range(11)}
    
    for k, count in proteoform_distribution.items():
        if count > 0:
            # Oxidation: Move from k to k+1
            if k < 10:
                if k == 0:
                    oxidized = np.random.binomial(count, P_oxidation_Cys215)  # Cys215 oxidation
                else:
                    oxidized = np.random.binomial(count, P_oxidation_other)  # Other cysteines oxidation
                new_distribution[k + 1] += oxidized
                new_distribution[k] += count - oxidized
            else:
                new_distribution[k] += count
            
            # Reduction: Move from k to k-1
            if k > 0:
                if k == 1:
                    reduced = np.random.binomial(count, P_reduction_Cys215)  # Cys215 reduction
                else:
                    reduced = np.random.binomial(count, P_reduction_other)  # Other cysteines reduction
                new_distribution[k - 1] += reduced
                new_distribution[k] -= reduced

    proteoform_distribution = new_distribution
    
    # Print the current distribution of PTP1B molecules across k values
    print(f"Time Step {step}: Proteoform Distribution = {proteoform_distribution}")

# Final distribution of proteoforms
print(f"Final Proteoform Distribution: {proteoform_distribution}")

# Calculate the percentage of molecules in k = 1
percentage_k1 = (proteoform_distribution[1] / number_of_PTP1B_molecules) * 100
print(f"Percentage of molecules in k = 1: {percentage_k1:.4f}%")

# Calculate the overall weighted mean redox state
total_molecules = sum(proteoform_distribution.values())
weighted_mean_redox_state = sum(k * count for k, count in proteoform_distribution.items()) / total_molecules / 10
print(f"Overall weighted mean redox state: {weighted_mean_redox_state:.4f}")

# Generate the Excel file
site_columns = proteoform_df.columns[1:-1]  # Extract the cysteine site columns
output_df = pd.DataFrame(columns=['Proteoform_ID', 'k_value'] + list(site_columns))

# Populate the Excel file with the final distribution
rows = []
for k, count in proteoform_distribution.items():
    if count > 0:
        for proteoform in proteoform_library[k]:
            proteoform_data = {
                'Proteoform_ID': proteoform_df[(proteoform_df.iloc[:, 1:-1] == proteoform).all(axis=1)].index[0],
                'k_value': k
            }
            proteoform_data.update({site: state for site, state in zip(site_columns, proteoform)})
            rows.append(proteoform_data)

# Convert the list of rows to a DataFrame
output_df = pd.DataFrame(rows)

# Save the DataFrame to an Excel file
output_df.to_excel('PTP1B_proteoform_distribution.xlsx', index=False)

print("Proteoform distribution has been saved to 'PTP1B_proteoform_distribution.xlsx'")
