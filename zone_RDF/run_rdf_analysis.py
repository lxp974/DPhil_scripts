# Import necessary modules
from pathlib import Path
import sys

# Path to your script
sys.path.append('.')

# Import the script's main function
from md_rdf_analysis import MDTrajAnalyser

# Define file paths
structure_file = Path('md.gro')
trajectory_file = Path('md_mol.xtc')

# Create analyzer and run
analyzer = MDTrajAnalyser(trajectory_file, structure_file)
analyzer.run_analysis([
    {
            'name': 'bulk',
            'ion_selection': 'name CL and prop 80 < z',
            'partner_selection': 'name OW',
            'frame_range': [20000,25000]
    },

        {
            'name': 'center',
            'ion_selection': 'name CL and cylayer 0 1.75 10 -10 (resname UNL)',
            'partner_selection': 'name OW',
            'frame_range': [20000,25000]
    },
    
        {
            'name': 'midinner',
            'ion_selection': 'name CL and cylayer 1.75 3.5 10 -10 (resname UNL)',
            'partner_selection': 'name OW',
            'frame_range': [20000,25000]
    },
    
        {
            'name': 'midouter',
            'ion_selection': 'name CL and cylayer 3.5 5.25 10 -10 (resname UNL)',
            'partner_selection': 'name OW',
            'frame_range': [20000,25000]
    },

        {
            'name': 'interface',
            'ion_selection': 'name CL and cylayer 5.25 7 10 -10 (resname UNL)',
            'partner_selection': 'name OW',
            'frame_range': [20000,25000],
    }
])
