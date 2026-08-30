import pandas as pd
import numpy as np
import sys # Import the sys module

def calculate_allelic_contribution(input_file, output_file, threshold=0.5):
    """
    Process a four-allele FPKM file, filter rows whose mean is below the threshold,
    and calculate the relative contribution of each allele (numeric range 0–1), adding the results as new columns.

    Args:
        input_file (str): Path to the input FPKM file.
        output_file (str): Path to the output file.
        threshold (float): Mean FPKM threshold used to filter rows.
    """

    # --- 1. Read and parse data ---
    with open(input_file, 'r') as f:
        header_line = f.readline().strip()
        data_lines = f.readlines()

    # Parse the header to identify each allele name and its value column name
    header_parts = header_line.split()
    allele_names = [] # For example, ['A', 'B', 'C', 'D']
    original_header_list = [] 
    for i in range(0, len(header_parts), 2):
        allele_letter = header_parts[i].replace('Allele', '')
        allele_names.append(allele_letter)
        original_header_list.append(header_parts[i]) # AlleleA
        original_header_list.append(header_parts[i+1]) # tetra_leaf_ave

    # Prepare storage for the parsed data
    parsed_data = []
    fpkm_value_cols_in_df = [] 
    
    # Build the initial list of DataFrame column names
    df_initial_columns = []
    for allele_letter in allele_names:
        df_initial_columns.append(f'Allele{allele_letter}_ID')
        df_initial_columns.append(f'Allele{allele_letter}_FPKM_Value')
        fpkm_value_cols_in_df.append(f'Allele{allele_letter}_FPKM_Value')

    for line in data_lines:
        parts = line.strip().split()
        if not parts: # Skip empty lines
            continue
        
        row_values = []
        for i, allele_letter in enumerate(allele_names):
            row_values.append(parts[2*i]) # Allele_ID
            fpkm_str = parts[2*i + 1]
            row_values.append(float(fpkm_str) if fpkm_str != 'NA' else 0.0) # FPKM_Value (handle NA)
        parsed_data.append(row_values)

    df = pd.DataFrame(parsed_data, columns=df_initial_columns)

    # --- 2. Clean data and perform calculations ---
    df['Total_FPKM_Sum'] = df[fpkm_value_cols_in_df].sum(axis=1)
    df['Average_FPKM'] = df[fpkm_value_cols_in_df].mean(axis=1)

    # --- 3. Filter rows ---
    df_filtered = df[df['Average_FPKM'] >= threshold].copy()

    if df_filtered.empty:
        print(f"No rows passed the mean FPKM threshold of {threshold}.")
        # Create and output an empty DataFrame containing only the original ID and FPKM columns plus the new relative-contribution columns
        final_output_cols = [col for col in df.columns if 'ID' in col or 'FPKM_Value' in col] + \
                            [f'Relative_{al}_FPKM' for al in allele_names] # Note that _Pct has been removed from the column names
        pd.DataFrame(columns=final_output_cols).to_csv(output_file, sep='\t', index=False)
        return

    # --- 4. Calculate relative contributions (numeric form) ---
    relative_contribution_output_cols = []
    for allele_letter in allele_names:
        fpkm_col_name = f'Allele{allele_letter}_FPKM_Value'
        # Modify the column name so that it no longer contains _Pct
        relative_col_name = f'Relative_{allele_letter}_FPKM' 
        
        # Avoid division by zero: if Total_FPKM_Sum is 0, set the relative contribution to 0
        df_filtered[relative_col_name] = np.where(
            df_filtered['Total_FPKM_Sum'] == 0,
            0,
            (df_filtered[fpkm_col_name] / df_filtered['Total_FPKM_Sum']) # Do not multiply by 100
        )
        relative_contribution_output_cols.append(relative_col_name)

    # --- 5. Prepare the output file ---
    # Output header: original header plus the new relative-contribution column names
    new_header_line = ' '.join(original_header_list) + ' ' + ' '.join(relative_contribution_output_cols)

    # Output DataFrame containing only the original ID/FPKM values and the newly calculated percentage columns
    output_df_final = df_filtered[df_initial_columns + relative_contribution_output_cols].copy()
    
    # --- 6. Write the new file ---
    with open(output_file, 'w') as f:
        f.write(new_header_line.replace(' ', '\t') + '\n') # Write the new header
        # Write the data without the DataFrame header and control the floating-point format
        output_df_final.to_csv(f, sep='\t', index=False, header=False, float_format='%.6f')

    print(f"Processing completed! The filtered data have been saved to '{output_file}'.")
    print(f"The original file contained {len(df)} rows; {len(df_filtered)} rows were retained after filtering.")

# --- Usage example ---
if __name__ == "__main__":
    # Use sys.argv to retrieve command-line arguments
    # sys.argv[0] is the script name
    # sys.argv[1] should be the input file
    # sys.argv[2] should be the output file
    
    if len(sys.argv) < 3:
        print("Usage: python your_script_name.py <input_fpkm_file> <output_file_name> [threshold]")
        print("Example: python process_fpkm.py input.txt output.txt 0.5")
        sys.exit(1) # Exit with an error status

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    
    # Check whether the optional threshold argument was provided
    threshold_value = 0.5 # Default value
    if len(sys.argv) > 3:
        try:
            threshold_value = float(sys.argv[3])
        except ValueError:
            print(f"Warning: Invalid threshold value '{sys.argv[3]}'. Using default threshold 0.5.")

    calculate_allelic_contribution(input_filename, output_filename, threshold=threshold_value)

    print(f"\nPlease check the generated file '{output_filename}'.")
