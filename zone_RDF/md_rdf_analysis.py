#!/usr/bin/env python3
import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import MDAnalysis as mda
import MDAnalysis.analysis.rdf as mdaRDF
import numpy as np
import panda as pd

"""
Calculates RDFs of ions within radial sections of a carbon nanotube (CNT) pore.

This script parses a trajectory of a CNT system and identifies ions within user defined
radial sections of the pore. Ions are dynamically selected such that only ions 
residing in a specified region at a given time are select and the RDF is calculated.
I.e. once an ion has moved outside the boundary of the selected region, it is not 
included in the RDF calculation. 

Outputs the bins with RDF values and the chloride coordination numbers for each frame. 
"""
@dataclass
class MDTrajAnalyser:
    """
    Analysis of MD trajectories

    Calculates the RDF of (chloride) ions that reside in a given radial section of the internal CNT pore.
    """
    
    trajectory_file: Path
    structure_file: Path
    output_dir: Path = field(default_factory=lambda: Path('./analysis_results'))
    
    def __init__(self, trajectory_file: Path, structure_file: Path):
        self.trajectory_file = trajectory_file
        self.structure_file = structure_file
        self.output_dir = Path('./analysis_results')
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def analyze_rdf(
        self, 
        ion_selection: str, 
        partner_selection: str, 
        region: Optional[Tuple[float, float]] = None,
        frame_range: Optional[Tuple[int, int]] = None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        
        """
        Compute Radial Distribution Function (RDF) with dynamic ion selection for a given specified region
        
        Args:
            ion_selection: MDAnalysis atom selection for ions
            partner_selection: Atom selection for interaction partner i.e. Water-Oxygen
            region: Radial region boundaries (lower, upper)
            frame_range: Optional trajectory frame range
        
        Returns:
            Tuple of (RDF, cumulative RDF, coordination number, bins)
        """
        
        # Load universe
        u = mda.Universe(str(self.structure_file), str(self.trajectory_file))
        
        # Select frames
        if frame_range is None:
            frame_range = (20000, 25000)
        
        # Configurable region selection
        if region is not None:
            lower, upper = region
            ion_selection += f" and cylayer {lower} {upper} 10 -10"
    

        # RDF computation
        rdf_list, rdf_coord, rdf_cumu = [], [], []
        bins = None
        
        for frame in u.trajectory[slice(*frame_range)]:
            print(f"Analysing frame: {frame.frame}")

            # Selections
            ions = u.select_atoms(ion_selection)
            partners = u.select_atoms(partner_selection)

            # If there are no ions present
            if len(ions) == 0:
                continue  # Skip frames with no ions selected
            
            # Calculate RDF
            rdf = mdaRDF.InterRDF(ions, partners, nbins=150)
            rdf.run(frame.frame, frame.frame+1)

            if bins is None:
                bins = rdf.bins
                
            count = np.cumsum(rdf.count) / len(ions)
            rdf_list.append(rdf.rdf)
            rdf_coord.append(count[38])  # at 0.38 nm
            rdf_cumu.append(count)
        

        # Average results
        rdf_avg = np.mean(rdf_list, axis=0)

        coord_avg = np.mean(rdf_coord, axis=0)
        print (f'The average ion coordination number is {coord_avg}')
        
        cumu_avg = np.mean(rdf_cumu, axis=0)
        
        return bins, rdf_avg, cumu_avg, rdf_coord
    
    def plot_results(
        self, 
        bins: np.ndarray, 
        rdf: np.ndarray, 
        cumu_rdf: np.ndarray, 
        region_name: str
    ):
        """Generate and save RDF and cumulative RDF plots."""

        plt.figure(0)
        plt.plot(bins,rdf)
        plt.title(f'RDF - {region_name}')
        plt.xlabel('Radius (Å)')
        plt.ylabel('g(r)')
        plt.savefig(self.output_dir / f'{region_name}_rdf.png')
        plt.close()
        
        plt.figure(1)
        plt.plot(bins, cumu_rdf)
        plt.title(f'Cumulative RDF - {region_name}')
        plt.xlabel('Radius (Å)')
        plt.ylabel('Coordination Number')
        plt.savefig(self.output_dir / f'{region_name}_cumuav.png')
        plt.close()
    
    def run_analysis(self, configs: List[dict]):
        """
        Execute multiple RDF analyses based on configuration.
        
        Args:
            configs: List of configuration dictionaries for analysis
        """
        for config in configs:
            region_name = config.get('name', 'default')
            bins, rdf_avg, cumu_avg, rdf_coord = self.analyze_rdf(
                ion_selection=config['ion_selection'],
                partner_selection=config['partner_selection'],
                region=config.get('region'),
                frame_range=config.get('frame_range')
            )


            pd.DataFrame({
                'Radius (Å)': bins, 'RDF': rdf_avg
                }).to_csv(self.output_dir / f'{region_name}_rdf.csv', index=False)


            np.savetxt(
                self.output_dir / f'{region_name}_coordno.csv', 
                rdf_coord,
                header='Coordination number' 
            )
            
            # Plot graphs
            self.plot_results(bins, rdf_avg, cumu_avg, region_name)


def main():
    parser = argparse.ArgumentParser(description='Molecular Dynamics RDF Analysis')
    parser.add_argument('structure', type=Path, help='Input structure file (.gro)')
    parser.add_argument('trajectory', type=Path, help='Input trajectory file (.xtc)')
    
    args = parser.parse_args()
    
    analysis_configs = [
        {'name': 'bulk', 'ion_selection': 'name CL and prop 80 < z', 'partner_selection': 'name OW'},
        {'name': 'center', 'ion_selection': 'name CL and cylayer 0 1.75 10 -10 (resname UNL)', 'partner_selection': 'name OW'},
        {'name': 'midinner', 'ion_selection': 'name CL and cylayer 1.75 3.5 10 -10 (resname UNL)', 'partner_selection': 'name OW'},
        {'name': 'midouter', 'ion_selection': 'name CL and cylayer 3.5 5.25 10 -10 (resname UNL)', 'partner_selection': 'name OW'},
        {'name': 'interface', 'ion_selection': 'name CL and cylayer 5.25 7 10 -10 (resname UNL)', 'partner_selection': 'name OW'}
    ]
    
    analyzer = MDTrajAnalyser(args.trajectory, args.structure)
    analyzer.run_analysis(analysis_configs)

if __name__ == "__main__":
    main()
