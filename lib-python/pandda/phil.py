
pandda_phil_def = """
    pandda
    {
        input
            .help = "File names"
        {
            data_dirs = None
                .type = str
            pdb_style = final.pdb
                .type = str
            mtz_style = None
                .type = str
            lig_style = *.cif
                .type = str
            regex
                .help = "Advanced dataset labelling regexs"
            {
                dir_regex = None
                    .type = str
                pdb_regex = None
                    .type = str
                mtz_regex = None
                    .type = str
            }
            reference
                .help = "Manually define reference structure"
            {
                pdb = None
                    .type = path
                mtz = None
                    .type = path
                amp_label = None
                    .type = str
                pha_label = None
                    .type = str
            }
            flags
                .help = "Flags for individual datasets"
            {
                ignore_datasets = None
                    .help = 'Reject these datasets, don\'t even load them - comma separated list of dataset tags'
                    .type = str
                no_analyse = None
                    .help = 'Don\'t analyse these datasets, only use them to build the distributions - comma separated list of dataset tags'
                    .type = str
                no_build = None
                    .help = 'Don\'t use these datasets to build the distributions, only analyse them - comma separated list of dataset tags'
                    .type = str
            }

        }
        output
            .help = "Output directory"
        {
            out_dir = './pandda'
                .type = path
            dataset_prefix = ''
                .help = "Prefix to be attached to each dataset name"
                .type = str
            plot_graphs = True
                .help = "Output graphs using matplotlib"
                .type = bool
        }
        method
            .help = "High Level control of algorithm"
        {
            reload_existing_datasets = True
                .help = "Reload existing datasets? - if False, will only load new datasets (unprocessed datasets)"
                .type = bool
            reprocess_existing_datasets = False
                .help = "Reprocess existing datasets? (Time-consuming) - if False, will only calculate z-maps for new datasets (unprocessed datasets)"
                .type = bool
            reprocess_selected_datasets = None
                .help = "Reprocess selection of datasets (comma-separated list)"
                .type = str
            recalculate_statistical_maps = False
                .help = "Recalculate all of the statistical maps? (Time-consuming) - if False, looks for existing statistical maps and uses those (reverts to True if none are found)"
                .type = bool
        }
        params
            .help = "Algorithm Parameters"
        {
            alignment
                .help = "Settings to control the alignment of the structures"
            {
                method = global *local
                    .type = choice
            }
            filtering
                .help = "Settings to control when datasets are rejected from the analysis"
            {
                max_rmsd_to_reference = 1.5
                    .help = "Reject datasets that have a calpha rmsd of greater than this to the reference (after alignment)"
                    .type = float
                max_rfree = 0.4
                    .help = 'Maximum allowed rfree for a structure (datasets above this are rejected)'
                    .type = float
            }
            maps
                .help = "Settings to control the generation of the dataset maps"
            {
                ampl_label = FWT *2FOFCWT
                    .type = choice
                phas_label = PHWT *PH2FOFCWT
                    .type = choice
                scaling = none *sigma volume
                    .type = choice
                resolution_factor = 0.25
                    .help = 'Sampling factor for fft-ing the maps'
                    .type = float
                grid_spacing = 0.5
                    .help = 'Spacing of the grid points in the sampled maps (A) - fixed across resolutions'
                    .type = float
                padding = 3
                    .help = "Padding around the edge of the maps (A)"
                    .type = float
            }
            masks
                .help = "Parameters to control the masking of grid points around the protein"
            {
                inner_mask = 1.8
                    .help = "Points are masked within this distance of protein atoms"
                    .type = float
                outer_mask = 6
                    .help = "Points are masked outside this distance of protein atoms"
                    .type = float
            }
            analysis
                .help = "Settings to control the selection of datasets"
            {
                min_build_datasets = 40
                    .help = 'Minimum number of datasets needed to build distributions'
                    .type = int
                max_build_datasets = 75
                    .help = 'Maximum number of datasets used to build distributions'
                    .type = int
                dynamic_res_limits = True
                    .help = 'Allow the analysed resolution limits to change depending on the dataset resolution ranges'
                    .type = bool
                high_res_upper_limit = 0.0
                    .help = 'Highest resolution limit (maps are never calulcated above this limit)'
                    .type = float
                high_res_lower_limit = 4.0
                    .help = 'Lowest resolution limit (datasets below this are ignored)'
                    .type = float
                high_res_increment = 0.05
                    .help = 'Increment of resolution shell for map analysis'
                    .type = float
            }
            z_map
                .help = "Settings to control the calculation of z-maps"
            {
                map_type = naive adjusted uncertainty *adjusted+uncertainty
                    .help = 'Type of Z-map to calculate'
                    .type = choice
            }
            blob_search
                .help = "Settings to control the finding of blobs"
            {
                contour_level = 2.5
                    .help = 'Contour level when looking for blobs'
                    .type = float
                min_blob_volume = 10.0
                    .help = 'Blob volume filter for detecting blobs'
                    .type = float
                min_blob_z_peak = 3.0
                    .help = 'Blob Z-peak filter for detecting blobs'
                    .type = float
                blobs_to_image = 2
                    .help = 'Maximum number of blobs to image per dataset'
                    .type = int
            }
            background_correction
                .help = "Parameters to control the estimation of feature background corrections"
            {
                min_bdc = 0.0
                    .type = float
                    .help = 'Minimum background correction estimate'
                max_bdc = 1.0
                    .type = float
                    .help = 'Maximum background correction estimate'
                increment = 0.01
                    .type = float
                    .help = 'Resolution of background correction estimation'
            }
        }
        settings
            .help = "General Settings"
        {
            cpus = 1
                .type = int
            verbose = True
                .type = bool
            max_memory = 25.0
                .type = float
            exit_flags
                .help = "Flags for terminating the program early"
            {
                dry_run = False
                    .help = "Setup pandda, print working arguments, and exit"
                    .type = bool
                calculate_first_mean_map_only = False
                   .help = "Will calculate the highest resolution mean map and then exit - used for initial reference modelling"
                    .type = bool
            }
            developer
                .help = "Development Settings (Not needed by most users)"
            {
                write_maps_for_empty_datasets = False
                    .help = "Output maps for all datasets, not just datasets marked as interesting"
                    .type = bool
                write_grid_masks = False
                    .help = "Output the grid masks which control which areas are analysed"
                    .type = bool
                write_all_z_map_types = False
                    .help = "Output all possible types of Z-maps (Takes a lot of space on disk)"
                    .type = bool
                validate_output = False
                    .help = "Validate hits against a known list of ligands"
                    .type = bool
            }
            pickling
                .help = "Pickle Settings"
            {
                full_pickle = False
                    .help = "Pickle the entire PANDDA object for reloading - takes a lot of space"
                    .type = bool
                pickle_maps = False
                    .help = "Pickle the dataset maps as part of the dataset pickle object"
                    .type = bool
            }
        }
    }
"""
