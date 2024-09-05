import numpy as np

# Define initial conditions
num_molecules = 70000  # Total number of PTP1B molecules
oxidation_probability = 0.1
reduction_probability = 0.9

# Starting with all molecules in the reduced state (k = 0)
molecules_in_k0 = num_molecules
molecules_in_k1 = 0

# Define the number of steps for the simulation
num_steps = 173

# Run the simulation
for step in range(num_steps):
    # Calculate the number of molecules transitioning from k = 0 to k = 1
    transitioning_to_k1 = np.random.binomial(molecules_in_k0, oxidation_probability)
    
    # Calculate the number of molecules transitioning from k = 1 to k = 0
    transitioning_to_k0 = np.random.binomial(molecules_in_k1, reduction_probability)
    
    # Update the counts
    molecules_in_k0 += transitioning_to_k0 - transitioning_to_k1
    molecules_in_k1 += transitioning_to_k1 - transitioning_to_k0

    # Print the results at each step
    print(f"Time Step {step+1}: k0 = {molecules_in_k0}, k1 = {molecules_in_k1}")

# Final distribution
total_molecules = molecules_in_k0 + molecules_in_k1
percentage_in_k1 = (molecules_in_k1 / total_molecules) * 100

# Print the final results
print(f"\nFinal Distribution after {num_steps} steps:")
print(f"Number of molecules in k = 0: {molecules_in_k0}")
print(f"Number of molecules in k = 1: {molecules_in_k1}")
print(f"Total number of molecules: {total_molecules}")
print(f"Percentage of molecules in k = 1: {percentage_in_k1:.4f}%")
