import numpy as np
import pandas as pd
import time

# Simulated smaller parameters for quick testing
number_of_molecules = 100  # Reduced for testing
time_steps = 10  # Reduced for testing

# Provided oxidation and reduction probabilities
P_oxidation_Cys215 = 0.1
P_oxidation_other = 0.028
P_reduction_Cys215 = 0.9
P_reduction_other = 0.1286

# Load the proteoform library CSV file (update with your file path in Colab)
proteoform_df = pd.read_csv('/content/PTP1B_proteoforms_ordered.csv')

# Initial proteoform state for all molecules
initial_proteoform = 'PF001'

# Initialize a route map DataFrame
route_map = pd.DataFrame(index=np.arange(number_of_molecules), columns=np.arange(time_steps + 1))

# Set all molecules to the initial proteoform at step 0
route_map.iloc[:, 0] = initial_proteoform

# Initialize variables to track unique proteoforms and k distribution
unique_proteoforms = set()
final_k_distribution = np.zeros(11)  # Assuming 10 cysteines, so k ranges from 0 to 10
final_proteoform_counts = {}

# Utility function to find the proteoform ID based on state vector
def find_proteoform_id(state_vector):
    return proteoform_df.loc[(proteoform_df.iloc[:, 1:-1] == state_vector).all(axis=1), 'Unnamed: 0'].values[0]

# Utility function to check valid transitions
def valid_transition(current_state, target_state):
    diff = np.array(target_state) - np.array(current_state)
    oxidations = np.where(diff == 1)[0]
    reductions = np.where(diff == -1)[0]

    # For oxidation, only allow if oxidizing a previously reduced site
    if len(oxidations) == 1 and len(reductions) == 0:
        return True

    # For reduction, only allow if reducing a previously oxidized site
    if len(oxidations) == 0 and len(reductions) == 1:
        return True

    return False

# Prepare a dictionary for quick lookup by k value
proteoform_library = proteoform_df.groupby(proteoform_df.iloc[:, 1:-1].sum(axis=1)).apply(lambda x: x.iloc[:, 1:-1].values.tolist()).to_dict()

# Function to print progress
def print_progress(step, total_steps):
    print(f"Step {step}/{total_steps} completed.")

# Start timing the simulation
start_time = time.time()

# Simulation loop
for step in range(1, time_steps + 1):
    for i in range(number_of_molecules):
        current_state = route_map.iloc[i, step - 1]  # Get the proteoform from the previous step
        
        # Convert current_state to the actual state (binary vector) based on the proteoform ID
        current_vector = proteoform_df.loc[proteoform_df['Unnamed: 0'] == current_state].iloc[0, 1:-1].tolist()
        
        k_value = sum(current_vector)
        oxidized = False
        reduced = False
        
        # Oxidation step
        if k_value < max(proteoform_library.keys()):  # Ensure we are not at the max k value
            oxidation_prob = P_oxidation_Cys215 if current_vector[3] == 0 else P_oxidation_other
            if np.random.rand() < oxidation_prob:
                possible_states = proteoform_library[k_value + 1]  # Get possible states for k + 1
                valid_states = [s for s in possible_states if valid_transition(current_vector, s)]
                if valid_states:
                    new_vector = valid_states[np.random.randint(len(valid_states))]
                    new_proteoform = find_proteoform_id(new_vector)
                    oxidized = True

        # Reduction step (if no oxidation occurred)
        if k_value > 0 and step > 1 and not oxidized:
            reduction_prob = P_reduction_Cys215 if current_vector[3] == 1 else P_reduction_other
            if np.random.rand() < reduction_prob:
                possible_states = proteoform_library[k_value - 1]  # Get possible states for k - 1
                valid_states = [s for s in possible_states if valid_transition(current_vector, s)]
                if valid_states:
                    new_vector = valid_states[np.random.randint(len(valid_states))]
                    new_proteoform = find_proteoform_id(new_vector)
                    reduced = True

        # If no change, retain the current proteoform
        if not oxidized and not reduced:
            new_proteoform = current_state

        # Update the route map for the current step
        route_map.iloc[i, step] = new_proteoform

        # Track unique proteoforms formed
        unique_proteoforms.add(new_proteoform)
        
        # Update final k value distribution
        if step == time_steps:
            final_k_distribution[sum(new_vector)] += 1
            if new_proteoform in final_proteoform_counts:
                final_proteoform_counts[new_proteoform] += 1
            else:
                final_proteoform_counts[new_proteoform] = 1

    # Print progress after each step
    print_progress(step, time_steps)

# End timing
end_time = time.time()
print(f"Simulation completed in {end_time - start_time:.2f} seconds.")

# Calculate the overall redox state of the population
overall_redox_state = sum(k * count for k, count in enumerate(final_k_distribution)) / number_of_molecules

# Print the results
print(f"Total unique proteoforms formed: {len(unique_proteoforms)} out of 1,024")
print(f"Final k distribution:\n{final_k_distribution}")
print(f"Overall redox state of the population: {overall_redox_state}")

# Save the route map to a CSV file
route_map.to_csv("/content/proteoform_route_map.csv", index=False)
print("Route map saved to /content/proteoform_route_map.csv")

# Create a DataFrame from the final proteoform counts and save to an Excel file
final_counts_df = pd.DataFrame(list(final_proteoform_counts.items()), columns=['Proteoform_ID', 'Count'])
final_counts_df['k_value'] = final_counts_df['Proteoform_ID'].apply(lambda pid: sum(proteoform_df.loc[proteoform_df['Unnamed: 0'] == pid].iloc[0, 1:-1]))

# Save the results to an Excel file
output_file = "/content/PTP1B_proteoform_final_distribution_with_counts.xlsx"
final_counts_df.to_excel(output_file, index=False)
print(f"Final distribution saved to {output_file}")
