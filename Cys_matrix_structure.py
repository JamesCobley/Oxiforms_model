import numpy as np
import pandas as pd
import itertools

# Define the cysteine residue positions for PTP1B
cysteine_positions = [32, 92, 121, 215, 226, 231, 324, 344, 414, 426]

# Number of cysteines
num_cysteines = len(cysteine_positions)

# Initialize an empty list to hold all the proteoforms in the correct order
proteoforms_ordered = []

# Iterate over k (the number of oxidized cysteines)
for k in range(num_cysteines + 1):
    # Generate all combinations of k oxidized cysteines
    for combo in itertools.combinations(range(num_cysteines), k):
        # Create a binary list representing the proteoform
        proteoform = [0] * num_cysteines
        for index in combo:
            proteoform[index] = 1
        proteoforms_ordered.append(proteoform)

# Convert the ordered proteoforms to a NumPy array
proteoforms_ordered = np.array(proteoforms_ordered, dtype=int)

# Generate proteoform names
proteoform_names = [f'PF{i+1:03d}' for i in range(len(proteoforms_ordered))]

# Convert binary arrays to strings
binary_structures = [''.join(map(str, row)) for row in proteoforms_ordered]

# Create a DataFrame with two columns: Proteoform and Structure
proteoform_df = pd.DataFrame({
    'Proteoform': proteoform_names,
    'Structure': binary_structures
})

# Display the first few rows of the DataFrame
print(proteoform_df.head(10))

# Export the DataFrame to a CSV file
proteoform_df.to_csv('PTP1B_proteoforms_ordered.csv', index=False)
