import numpy as np
import pandas as pd

# Load the proteoform library from the provided CSV file
proteoform_df = pd.read_csv('/content/PTP1B_proteoforms_ordered.csv')

# Drop the first column which is likely an index or placeholder
proteoform_df = proteoform_df.drop(columns=['Unnamed: 0'])

# Constants
number_of_PTP1B_molecules = 10000  # Starting with 10,000 molecules
time_steps = 173  # 173 simulation steps

# Provided oxidation and reduction probabilities for 9 micromoles of H2O2
P_oxidation_Cys215 = 0.1  # Oxidation probability for Cys215
P_oxidation_other = 0.028  # Oxidation probability for other sites
P_reduction_Cys215 = 0.9  # Reduction probability for Cys215
P_reduction_other = 0.1286  # Reduction probability for other sites

# Initialize all molecules in PF001 (the first row of the CSV, assuming it is the fully reduced state)
initial_proteoform = proteoform_df.iloc[0, :].tolist()
proteoform_distribution = [initial_proteoform] * number_of_PTP1B_molecules

# Convert proteoform states to k values
proteoform_df['k_value'] = proteoform_df.iloc[:, :].sum(axis=1)
proteoform_library = proteoform_df.groupby('k_value').apply(lambda x: x.iloc[:, :-1].values.tolist()).to_dict()

# Function to find the proteoform ID given a state
def find_proteoform_id(state):
    for index, row in proteoform_df.iterrows():
        if np.array_equal(row[:-1].tolist(), state):  # Compare the state ignoring the k_value column
            return index
    return None

# Function to identify valid transitions
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

# Simulation loop
for step in range(time_steps):
    print(f"--- Step {step + 1} ---")
    new_distribution = []
    transitions = {'oxidized': 0, 'reduced': 0}
    for state in proteoform_distribution:
        original_state = state.copy()  # Track the original state for debugging
        k_value = sum(state)  # Determine k value
        current_library = proteoform_library[k_value]

        # Oxidation step
        oxidized = False
        if k_value < max(proteoform_library.keys()):
            oxidation_prob = P_oxidation_Cys215 if state[3] == 0 else P_oxidation_other  # Cys215 is at index 3
            if np.random.rand() < oxidation_prob:
                # Randomly choose a valid oxidation transition
                possible_states = proteoform_library[k_value + 1]
                valid_states = [s for s in possible_states if valid_transition(state, s)]
                if valid_states:
                    state = valid_states[np.random.randint(len(valid_states))]
                    oxidized = True
                    transitions['oxidized'] += 1

        # Reduction step (only if no oxidation occurred and step > 1)
        reduced = False
        if k_value > 0 and step > 1 and not oxidized:
            reduction_prob = P_reduction_Cys215 if state[3] == 1 else P_reduction_other  # Cys215 is at index 3
            if np.random.rand() < reduction_prob:
                # Randomly choose a valid reduction transition
                possible_states = proteoform_library[k_value - 1]
                valid_states = [s for s in possible_states if valid_transition(state, s)]
                if valid_states:
                    state = valid_states[np.random.randint(len(valid_states))]
                    reduced = True
                    transitions['reduced'] += 1

        if oxidized or reduced:
            print(f"Molecule transitioned at step {step + 1}: {original_state} -> {state}")

        new_distribution.append(state)
    
    proteoform_distribution = new_distribution  # Update distribution
    print(f"Transitions in step {step + 1}: {transitions}")

# Final distribution of proteoforms
final_distribution = pd.DataFrame(proteoform_distribution, columns=proteoform_df.columns[:-1])  # Exclude the k_value column
final_counts = final_distribution.groupby(final_distribution.columns.tolist(), as_index=False).size()

# Map back to proteoform IDs and k values
final_counts['Proteoform_ID'] = final_counts.iloc[:, :-1].apply(find_proteoform_id, axis=1)  # Use only the state columns for mapping
final_counts['k_value'] = final_counts.iloc[:, :-2].sum(axis=1)  # Calculate k_value based on state columns

# Output the results
output_file = '/content/PTP1B_proteoform_final_distribution_with_counts.xlsx'
final_counts.to_excel(output_file, index=False)
print(f"Final distribution saved to {output_file}")
