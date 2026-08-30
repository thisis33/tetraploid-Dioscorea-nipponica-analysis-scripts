import pandas as pd
import numpy as np
import sys

def classify_allelic_expression(input_file, defination_file, output_file):
    """
    Classify the observed expression pattern of each allele quartet using predefined ideal expression patterns.

    Args:
        input_file (str): Path to the input file containing relative-contribution values for allele quartets.
                           This should be the output file from the previous script.
        defination_file (str): Path to the file defining the nine ideal expression patterns.
        output_file (str): Path to the new output file containing classification results.
    """

    # --- 1. Read the ideal expression-pattern definition file ---
    # Define column names to match defination.txt
    def_cols = ['Pattern_Name', 'A', 'B', 'C', 'D']
    # Note: here we assume that defination.txt has no header, or that its header must be skipped
    # The first line in the sample file is "v A B C D", which appears to be a header, so header=0 is used
    def_df = pd.read_csv(defination_file, sep='\s+', header=0, names=def_cols)
    
    # Set Pattern_Name as the index for convenient lookup
    def_df = def_df.set_index('Pattern_Name')
    
    # Extract the numeric portion of the ideal patterns
    ideal_patterns = def_df[['A', 'B', 'C', 'D']].values
    ideal_pattern_names = def_df.index.tolist()

    print(f"Loaded {len(ideal_pattern_names)} ideal expression patterns.")
    # print("Ideal pattern names:", ideal_pattern_names)

    # --- 2. Read the input file containing observed relative contributions ---
    # Because the output header from the previous script was constructed manually, it can be read directly here
    # Assume that the input file is tab-delimited
    df_input = pd.read_csv(input_file, sep='\t')

    # Extract the relative-contribution columns for the observed quartets.
    # Assume that they are the last four columns and use names such as Relative_A_FPKM and Relative_B_FPKM.
    # These column names must exactly match the output from the previous script
    relative_cols = ['Relative_A_FPKM', 'Relative_B_FPKM', 'Relative_C_FPKM', 'Relative_D_FPKM']
    
    # Check whether these columns exist
    if not all(col in df_input.columns for col in relative_cols):
        print(f"Error: Not all expected relative-contribution columns were found in the input file '{input_file}': {relative_cols}")
        print("Please ensure that the input file is the correct output from the previous script and that the column names have not been modified.")
        sys.exit(1)

    actual_contributions = df_input[relative_cols].values

    # --- 3. Calculate Euclidean distances and classify ---
    classifications = []
    min_distances = []

    for i in range(len(actual_contributions)):
        actual_vec = actual_contributions[i]
        
        distances = []
        for j in range(len(ideal_patterns)):
            ideal_vec = ideal_patterns[j]
            
            # Calculate the Euclidean distance
            distance = np.sqrt(np.sum((actual_vec - ideal_vec)**2))
            distances.append(distance)
        
        # Find the prototype with the smallest distance
        min_dist_index = np.argmin(distances)
        min_distance = distances[min_dist_index]
        classified_pattern = ideal_pattern_names[min_dist_index]
        
        classifications.append(classified_pattern)
        min_distances.append(min_distance)

    # --- 4. Add classification results to the DataFrame and write the output ---
    df_input['Classified_Pattern'] = classifications
    # The minimum distance can also be included as an output column if needed
    df_input['Min_Euclidean_Distance'] = min_distances 

    df_input.to_csv(output_file, sep='\t', index=False, float_format='%.6f')

    print(f"Classification completed! Results have been saved to '{output_file}'.")
    print(f"Processed {len(df_input)} rows in total.")

# --- Usage example ---
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python classify_expression.py <input_processed_fpkm_file> <defination_file> <output_classified_file>")
        print("Example: python classify_expression.py output_relative_contribution_final.txt defination.txt classified_results.txt")
        sys.exit(1) # Exit with an error status

    input_filename = sys.argv[1]
    defination_filename = sys.argv[2]
    output_filename = sys.argv[3] # Assume that the third argument is the output filename

    classify_allelic_expression(input_filename, defination_filename, output_filename)

    print(f"\nPlease check the generated file '{output_filename}'.")
