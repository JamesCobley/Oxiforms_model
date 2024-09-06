import pandas as pd

# Load both the Excel and CSV files
file_path_excel = '/content/PTP1B_proteoform_final_distribution_with_counts (12).xlsx'
file_path_csv = '/content/PTP1B_proteoforms_ordered.csv'

# Load the Excel file
df_excel = pd.read_excel(file_path_excel, sheet_name='Sheet1')  # Adjust 'Sheet1' if needed
print("Excel data loaded successfully.")

# Load the CSV file containing the proteoform matrix
df_matrix = pd.read_csv(file_path_csv)
print("Proteoform matrix loaded successfully.")

# Task 1: Define the number of proteoforms present and express as a percentage of 1,024
num_proteoforms = df_excel['Proteoform_ID'].nunique()
total_possible_proteoforms = 1024
proteoform_percentage = (num_proteoforms / total_possible_proteoforms) * 100
print(f"\nTask 1: Number of proteoforms: {num_proteoforms} ({proteoform_percentage:.2f}% of {total_possible_proteoforms})")

# Task 2: Print the number of molecules in each k_value
molecules_per_k_value = df_excel.groupby('k_value')['Count'].sum().reset_index()
print("\nTask 2: Number of molecules for each k_value:")
print(molecules_per_k_value)

# Task 3: Calculate the weighted mean redox state of the population based on k values (as a percentage)
total_molecules = df_excel['Count'].sum()
weighted_mean_redox = (df_excel['k_value'] * df_excel['Count']).sum() / total_molecules
weighted_mean_redox_percentage = weighted_mean_redox * 10  # Convert to percentage
print(f"\nTask 3: Weighted mean redox state of the population: {weighted_mean_redox_percentage:.2f}%")

# Task 4: Print the total number of molecules
print(f"\nTask 4: Total number of molecules: {total_molecules}")

# Task 5: Print the number of molecules in proteoform "PF005" and calculate its percentage of the total
proteoform_pf005 = df_excel[df_excel['Proteoform_ID'] == 'PF005']
if not proteoform_pf005.empty:
    molecules_pf005 = proteoform_pf005['Count'].sum()
    pf005_percentage = (molecules_pf005 / total_molecules) * 100
    print(f"\nTask 5: Number of molecules in proteoform 'PF005': {molecules_pf005} ({pf005_percentage:.2f}% of total molecules)")
else:
    print("\nProteoform 'PF005' not found in the data.")

# Task 6: Calculate cysteine oxidation state percentages based on the enumerated matrix
# Merge the matrix with the Excel data on "Proteoform_ID"
df_combined = pd.merge(df_excel, df_matrix, left_on='Proteoform_ID', right_on='Unnamed: 0')

# Calculate the oxidation percentage for each cysteine residue
cysteine_columns = df_matrix.columns[1:]  # All cysteine columns (e.g., 'Cys32', 'Cys92', ...)
oxidation_percentages = (df_combined[cysteine_columns].T * df_combined['Count']).sum(axis=1) / df_combined['Count'].sum() * 100

# Print the oxidation percentages for each cysteine
print("\nTask 6: Cysteine oxidation state percentages:")
for cysteine, percentage in oxidation_percentages.items():
    print(f"{cysteine}: {percentage:.2f}%")

# Task 7: Count proteoforms where Cys215 is oxidized (1)
proteoforms_with_cys215_oxidized = df_combined[df_combined['Cys215'] == 1]
num_proteoforms_cys215_oxidized = proteoforms_with_cys215_oxidized['Proteoform_ID'].nunique()
percentage_cys215_oxidized = (num_proteoforms_cys215_oxidized / 1024) * 100
print(f"\nTask 7: Number of proteoforms with Cys215 oxidized: {num_proteoforms_cys215_oxidized} ({percentage_cys215_oxidized:.2f}% of total possible proteoforms)")

# Task 8: Calculate percentage of PTP1B activation (proteoforms where Cys215 is oxidized as percentage of total population)
molecules_with_cys215_oxidized = proteoforms_with_cys215_oxidized['Count'].sum()
percentage_ptp1b_activation = (molecules_with_cys215_oxidized / total_molecules) * 100
print(f"\nTask 8: Percentage of PTP1B activation (Cys215 oxidized): {percentage_ptp1b_activation:.2f}%")

# Task 9: Print the top 10 most frequent proteoform IDs and their k_values
top_10_proteoforms = df_excel.nlargest(10, 'Count')[['Proteoform_ID', 'k_value']]
print("\nTask 9: Top 10 most frequent proteoform IDs and their k_values:")
print(top_10_proteoforms)
