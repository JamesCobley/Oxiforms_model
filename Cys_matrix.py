import numpy as np
import pandas as pd

# Define the cysteine residue positions for PTP1B
cysteine_positions = [32, 92, 121, 215, 226, 231, 324, 344, 414, 426]

# Number of cysteines
num_cysteines = len(cysteine_positions)

# Generate all possible binary combinations for the given number of cysteines
proteoforms = np.array([list(f"{i:0{num_cysteines}b}") for i in range(2**num_cysteines)], dtype=int)

# Convert to a DataFrame for better visualization with cysteine positions as column names
proteoform_df = pd.DataFrame(proteoforms, columns=[f"Cys{pos}" for pos in cysteine_positions])

# Add a proteoform identifier column
proteoform_df.index = [f'PF{i+1:03d}' for i in range(len(proteoform_df))]

# Display the first few rows of the matrix for verification
print(proteoform_df.head(10))

# Export the DataFrame to a CSV file
proteoform_df.to_csv('PTP1B_proteoforms.csv', index=True)
