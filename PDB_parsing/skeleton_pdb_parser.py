# Created by Ettore Bartalucci
# Test skeleton for parsing the PDB database and extract structures based on structural motif constraints
# Date: 29.09.23

# more infos at: https://search.rcsb.org/#search-api

# Modules
import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import requests # to make the HTTP requests to PDB RESTful API

import requests

def fetch_protein_structures(parameter):
    # Define the PDB RESTful API URL
    pdb_api_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # Define your criteria for searching (e.g., keyword, organism, resolution, etc.)
    # You can customize 'criteria' based on your specific requirements.
    parameter = f"your_search_criteria_here"
    
    # Create a query payload
    payload = {
        'queryType': 'org.pdb.query.simple.StructureIdQuery',
        'description': 'Query for Protein Structures',
        'query': parameter,
        'format': 'json',
    }
    
    # Make the GET request to the PDB API
    response = requests.get(pdb_api_url, params=payload)
    
    if response.status_code == 200:
        # Parse the JSON response
        data = response.json()
        
        # Extract the list of PDB IDs
        pdb_ids = data.get('result_set', [])
        
        # Print the PDB IDs
        for pdb_id in pdb_ids:
            print(f"Found PDB ID: {pdb_id}")
    else:
        print(f"Failed to retrieve data. Status code: {response.status_code}")

if __name__ == "__main__":
    # Specify your search criteria here (e.g., 'human', 'resolution < 2.0')
    search_criteria = "your_search_criteria_here"
    
    # Call the function to fetch protein structures based on the criteria
    fetch_protein_structures(search_criteria)
