import numpy as np
import pandas as pd

# Constants
number_of_PTP1B_molecules = 6022141500000000  # Total number of PTP1B molecules
time_steps = 173  # Number of simulation steps
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

output_df.to_excel('PTP1B_proteoform_final_distribution_with_counts.xlsx', index=False)
print("Final proteoform distribution with molecule counts has been saved to 'PTP1B_proteoform_final_distribution_with_counts.xlsx'")
