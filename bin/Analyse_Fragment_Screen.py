import os, sys, glob

from PANDDAs.Main import multi_dataset_analyser

# ============================================================================>
#####
# Output Files
#####
# ============================================================================>

output_dir = '/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak'

output_points = os.path.join(output_dir, 'sample_points.csv')
output_map_values = os.path.join(output_dir, 'map_values.csv')
output_map_point_summary = os.path.join(output_dir, 'map_point_summary.csv')
output_z_values = os.path.join(output_dir, 'z_values.csv')
output_z_values_summary = os.path.join(output_dir, 'z_values_summary.csv')
output_modified_z_values_summary = os.path.join(output_dir, 'mod_z_values_summary.csv')
if os.path.exists(output_points):
    os.remove(output_points)
if os.path.exists(output_map_values):
    os.remove(output_map_values)
if os.path.exists(output_map_point_summary):
    os.remove(output_map_point_summary)
if os.path.exists(output_z_values):
    os.remove(output_z_values)
if os.path.exists(output_z_values_summary):
    os.remove(output_z_values_summary)

# ============================================================================>
#####
# Reference Dataset
#####
# ============================================================================>

input_dir = '/home/npearce/Analyses/FragmentSoaks/BAZ2BA-Feb13-Soak'
ref_pdb = os.path.join(input_dir,'reference.pdb')
ref_mtz = os.path.join(input_dir,'reference.mtz')

# ============================================================================>
#####
# Create list of data
#####
# ============================================================================>

raw_pdbs = sorted(glob.glob(os.path.join(input_dir,'ProcessedFragmentSoak/BAZ2BA-*/1-apo/apo-BAZ2BA-*-refmac.pdb')))
raw_file_pairs = []

for pdb in raw_pdbs:
    mtz = pdb.replace('.pdb','.mtz')
    if os.path.exists(mtz):
        raw_file_pairs.append((pdb,mtz))

#raw_file_pairs = raw_file_pairs[0:5]

if not raw_file_pairs: raise SystemExit('No Files Found!')

num_raw_datasets = len(raw_file_pairs)
file_pairs = []

print 'Number of Datasets: ', num_raw_datasets

# ============================================================================>
# ============================================================================>
# TODO Change the name of the main processor object

pandda = multi_dataset_analyser()
pandda.set_map_type(map_type='2mFo-DFc')

# ============================================================================>
#####
# Set Reference Dataset
#####
# Select the reference dataset
# TODO Choose the reference dataset better
# ============================================================================>

pandda.load_reference_dataset(ref_pdb=ref_pdb, ref_mtz=ref_mtz)

# ============================================================================>
#####
# Create Sample Grid
#####
# Choose the resolution limit for the maps
# Create reference grid based on the reference structure
#
# TODO Choose the resolution limit better
# TODO Add a buffer to the grid
# ============================================================================>

pandda.set_cut_resolution(d_min=pandda.get_reference_dataset().get_dataset_resolution())
pandda.create_reference_grid(res_factor=0.33, include_origin=True, buffer=None)
pandda.extract_reference_map_values()

# ============================================================================>
#####
# Create Local and Global Masks
#####
# Local mask used for forming groups of points around a grid point
# Global mask used for removing points in the bulk solvent regions
# ============================================================================>

pandda.get_reference_grid().create_local_mask(distance_cutoff=4)
pandda.get_reference_grid().create_global_mask(cart_sites=pandda.get_reference_dataset().get_calpha_sites(), distance_cutoff=15)
pandda.mask_resampled_reference_grid()

# ============================================================================>
#####
# Scale and Process All Data
#####
# Read in all of the files
# Scale and align maps to reference
# TODO Revisit Scaling
# ============================================================================>

pandda.add_input_files(raw_file_pairs)
pandda.load_input_datasets()
pandda.scale_datasets_and_load_maps()
pandda.align_and_filter_datasets()
pandda.extract_map_values()

# ============================================================================>
#####
# ANALYSIS
#####
# Calculate the moments of the distributions at the grid points
# Use the means and the stds to convert the maps to z-maps
# Use the local mask to look for groups of significant z-values
# ============================================================================>

pandda.calculate_map_statistics()
pandda.normalise_maps_to_z_maps()
pandda.post_process_z_maps()

# ============================================================================>

