from chembl_webresource_client.new_client import new_client
import pandas as pd

# Initialize the activity API client
activity = new_client.activity

# Define the target ChEMBL ID and filtering conditions
target_chembl_id = 'CHEMBL1862'  # ABL1 human
pchembl_threshold = 10  # Corresponds to pIC50 >= 10 (i.e., IC50 <= 0.1 nM)

print(f"Fetching IC50 data for target {target_chembl_id} with pChEMBL value >= {pchembl_threshold}...")

try:
    # Build the query
    # We want to retrieve activity data related to the target,
    # with standard_type 'IC50', and pChEMBL value matching the criteria.
    # Also ensure pChEMBL value is not null (pchembl_value__isnull=False).
    # __gte means "greater than or equal to".
    # __in means the value is in the provided list.
    res = activity.filter(
        target_chembl_id=target_chembl_id,
        standard_type__iexact="IC50",  # iexact for case-insensitive exact match
        pchembl_value__gte=pchembl_threshold,
        pchembl_value__isnull=False, # Ensure pchembl_value exists
        assay_type__in=['B', 'F'] # Filter for multiple assay_types using __in
    ).only(
        'molecule_chembl_id', 'standard_type', 'standard_relation', 'standard_value', 'standard_units',
        'pchembl_value', 'assay_chembl_id', 'assay_description', 'document_chembl_id', 'target_organism'
    ) # .only() specifies which fields to return, reducing data volume

    # Convert the results to a list
    activities_list = list(res)

    if activities_list:
        # Convert the list to a Pandas DataFrame for better viewing and handling
        df_activities = pd.DataFrame(activities_list)
        
        print(f"\nFound {len(df_activities)} activities matching the criteria.")
        
        # Print the first few records and some key columns
        print("\nSample of filtered activities:")
        print(df_activities[[
            'molecule_chembl_id', 'standard_type', 'standard_relation', 'standard_value', 'standard_units', 'pchembl_value',
        ]].head())

        # Optionally save to a CSV file
        output_filename = "./hw2/q4.csv"
        df_activities.to_csv(output_filename, index=False)
        print(f"\nFull results saved to {output_filename}")

    else:
        print("No activities found matching the specified criteria.")

except Exception as e:
    print(f"An error occurred: {e}")