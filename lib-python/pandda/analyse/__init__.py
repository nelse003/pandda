import os, sys, glob, time, gc, traceback

#################################
try:
    import matplotlib
    matplotlib.interactive(False)
    from matplotlib import pyplot
    pyplot.style.use('ggplot')
except Exception as e:
    print e
#################################

import numpy

from libtbx import easy_mp
from libtbx.utils import Sorry, Failure

from bamboo.common import Meta, Info
from bamboo.common.path import rel_symlink

from giant.jiffies import extract_params_default

from pandda import welcome, module_info
from pandda.phil import pandda_phil
from pandda.analyse.classes import PanddaMultiDatasetAnalyser, PanddaMapAnalyser, PanddaDataset
from pandda.analyse.functions import DatasetProcessor, NativeMapMaker, wrapper_run
from pandda.analyse.z_maps import PanddaZMapAnalyser
from pandda.analyse import graphs as analyse_graphs
from pandda.analyse import html as analyse_html

# ============================================================================>
#
###                 PanDDA Initialisation Functions
#
# ============================================================================>

def pandda_dataset_setup(pandda):
    """Initialise the pandda object and load the input datasets"""

    # ============================================================================>
    # Report
    # ============================================================================>
    pandda.log.heading('Running PanDDA Setup', blank=True)

    # ============================================================================>
    # Build list of files in data directories
    # ============================================================================>
    new_files = pandda.build_input_list()
    # Add the new files to the pandda object
    pandda.add_files(file_list=new_files)
    # Check that some datasets have been found or already loaded
    pandda.check_number_of_datasets()
    # ============================================================================>
    #####
    # Load and process input files
    #####
    # ============================================================================>
    if pandda.new_files():
        # ============================================================================>
        # Load and check datasets
        # ============================================================================>
        pandda.load_new_datasets()
        pandda.initialise_dataset_masks_and_tables()
        if pandda.args.flags.reprocess_existing_datasets or pandda.args.flags.reprocess_selected_datasets:
            pandda.check_loaded_datasets(datasets=pandda.datasets.all())
        else:
            pandda.check_loaded_datasets(datasets=pandda.datasets.mask(mask_name='old datasets', invert=True))
        # ============================================================================>
        # Select Reference Dataset
        # ============================================================================>
        if not pandda.datasets.reference():
            # Filter datasets against the provided filter pdb if given
            if pandda.args.input.filter.pdb is not None:
                pandda.log.bar()
                pandda.log('Filtering datasets against the provided pdb structure (defined by pandda.input.filter.pdb)', True)
                pandda.filter_datasets_1(filter_dataset=PanddaDataset.from_file(model_filename=pandda.args.input.filter.pdb))
            # Use given reference dataset, or select reference dataset
            if pandda.args.input.reference.pdb and pandda.args.input.reference.mtz:
                pandda.log.bar()
                pandda.log('Reference Provided by User', True)
                ref_pdb, ref_mtz = pandda.args.input.reference.pdb, pandda.args.input.reference.mtz
            else:
                pandda.log.bar()
                pandda.log('Selecting reference dataset from loaded datasets', True)
                ref_pdb, ref_mtz = pandda.select_reference_dataset(method='resolution')
            # Load the reference dataset
            pandda.load_reference_dataset(ref_pdb=ref_pdb, ref_mtz=ref_mtz)
            pandda.log.bar()
        # ============================================================================>
        # Align and Filter All Datasets
        # ============================================================================>
        # Filter out datasets with different protein structures
        pandda.filter_datasets_1()
        # Align structures to the reference
        pandda.align_datasets(method=pandda.params.alignment.method)
        # Pickle the new loaded datasets
        pandda.pickle_the_pandda(components=['datasets'], datasets=pandda.datasets.mask(mask_name='old datasets', invert=True))
    else:
        # ============================================================================>
        # Build masks for reloaded datasets based on current input parameters (quick)
        # ============================================================================>
        pandda.initialise_dataset_masks_and_tables()
        pandda.check_loaded_datasets(datasets=pandda.datasets.all())
        pandda.filter_datasets_1()

    # ============================================================================>
    # Check that enough VALID datasets have been found
    # ============================================================================>
    pandda.check_number_of_datasets()

    # ============================================================================>
    # If dry_run, exit after initial search
    # ============================================================================>
    if pandda.args.exit_flags.dry_run:
        raise SystemExit('Dry run only -exiting')

    # ============================================================================>
    #####
    # Filter and Analyse the Datasets
    #####
    # ============================================================================>
    # Collate dataset variables and parameters
    pandda.collate_dataset_variables()
    # Filter out the datasets that are not isomorphous and therefore incomparable
    pandda.filter_datasets_2()
    # Write output csvs
    pandda.write_output_csvs()

    return

def pandda_variation_analysis(pandda):
    """Analyse the variation in the loaded datasets"""

    # ============================================================================>
    # Report
    # ============================================================================>
    pandda.log.heading('Analysing dataset variation')

    # ============================================================================>
    #####
    # Perform variation analysis on the loaded datasets
    #####
    # ============================================================================>
    # Calculate the variation between the alignments of neighbouring C-alphas
    pandda.analyse_alignment_variation()
    pandda.analyse_overall_bfactors()

def pandda_grid_setup(pandda):
    """Generate the grid objects for the analysis"""

    # ============================================================================>
    # Report
    # ============================================================================>
    pandda.log.heading('Running Grid Setup')

    # ============================================================================>
    #####
    # Create Sampling Grid (for generated maps)
    #####
    # Create reference grid based on the reference structure
    # ============================================================================>
    if pandda.grid is None:
        pandda.create_reference_grid(
            dataset=pandda.datasets.reference(),
            grid_spacing=pandda.params.maps.grid_spacing)
        pandda.mask_reference_grid(
            dataset=pandda.datasets.reference())
        pandda.partition_reference_grid(
            dataset=pandda.datasets.reference())
    else:
        pandda.log('Grid loaded from previous analysis')
    # ============================================================================>
    # Store for reuse
    # ============================================================================>
    # Pickle all of the large arrays so they can be reloaded
    # ============================================================================>
    pandda.pickle_the_pandda(components=['grid'])

    # ============================================================================>
    # If setup_only, exit after initial search
    # ============================================================================>
    if pandda.args.exit_flags.setup_only:
        pandda.exit(error_msg=None)
        raise SystemExit('Setup Only: Exiting')

    return

# ============================================================================>
#
###                 PanDDA Processing Functions
#
# ============================================================================>

def pandda_main_loop(pandda):
    """Calculate the statistical maps, and then look for events in each dataset"""

    # ============================================================================>
    # Report
    # ============================================================================>
    pandda.log('Running Main PanDDA Analysis')

    # ============================================================================>
    # Validate/Reset the loaded datasets (those to be analysed)
    # ============================================================================>
    # Remove previous event information for datasets to be re-analysed
    pandda.reset_loaded_datasets()
    # Load reflection data for each dataset
    if pandda.args.flags.reprocess_existing_datasets or pandda.args.flags.reprocess_selected_datasets:
        pandda.load_and_scale_diffraction_data(datasets=pandda.datasets.all())
    else:
        pandda.load_and_scale_diffraction_data(datasets=pandda.datasets.mask(mask_name='old datasets', invert=True))

    # ============================================================================>
    #####
    # PRE-ANALYSIS ANALYSIS (DUMP OF DATASET PARAMETERS)
    #####
    # ============================================================================>
    pandda.write_output_csvs()
    analyse_graphs.write_dataset_summary_graphs(pandda)
    analyse_html.write_initial_html(pandda)

    # ============================================================================>
    #####
    # PREPARE VARIABLES TO LOOP OVER RESOLUTION SHELLS
    #####
    # ============================================================================>
    # Calculate cutoffs for resolution shells
    # ============================================================================>
    res_limits = pandda.select_resolution_limits()
    # ============================================================================>
    # Analyse all datasets from high_shell_limit -> cut_resolution (initialise high_shell_limit to 0)
    # ============================================================================>
    high_shell_limit = 0
    # ============================================================================>
    # Record how many datasets are processed at each resolution
    # ============================================================================>
    resolution_count = {}
    # ============================================================================>
    # Report
    # ============================================================================>
    pandda.log.bar()
    if len(res_limits)==1:
        pandda.log('Analysing All Maps at {!s}A'.format(res_limits[0]), True)
    else:
        pandda.log('Analysing Resolution Shells from {!s} -> {!s}A'.format(res_limits[0], res_limits[-1]), True)
        pandda.log('Limits: {!s}'.format(', '.join(map(str,res_limits))), True)

    # ============================================================================>
    #####
    # ANALYSE DATASETS - ITERATE THROUGH RESOLUTION SHELLS
    #####
    # ============================================================================>
    t_analysis_start = time.time()
    # ==================================================>
    pandda.log.bar()
    pandda.log('Dataset Analysis Started: {!s}'.format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(t_analysis_start))), True)

    for cut_resolution in res_limits:

        # ============================================================================>
        # Report
        # ============================================================================>
        pandda.log.heading('Processing Datasets from {!s}A -> {!s}A'.format(high_shell_limit, cut_resolution), spacer=True)

        # ============================================================================>
        # Select datasets to build the distributions
        # ============================================================================>
        if cut_resolution in pandda.stat_maps.get_resolutions():
            # Create a dummy mask as we won't be using any datasets for building
            pandda.log.bar()
            pandda.log('Using existing statistical maps @ {}A'.format(cut_resolution), True)
            building_mask_name = 'dummy mask @ {}A'.format(cut_resolution)
            building_mask = [False]*pandda.datasets.size()
            pandda.datasets.all_masks().add_mask(name=building_mask_name, values=building_mask)
        else:
            pandda.log.bar()
            pandda.log('Selecting Building Mask', True)
            pandda.log.bar()
            building_mask_name, building_mask = pandda.select_datasets_for_density_characterisation(high_res_cutoff=cut_resolution)
            # Check that we have enough datasets to build distributions
            if sum(building_mask) < pandda.params.analysis.min_build_datasets:
                pandda.log('Not enough datasets at this resolution to perform statistical density analysis ({!s}<{!s})'.format(sum(building_mask),pandda.params.analysis.min_build_datasets), True)
                pandda.log.bar()
                continue
            else:
                pandda.log('Sufficient datasets ar this resolution to perform statistical density analysis', True)
                pandda.log.bar()
                pandda.log('Using {} datasets for statistical electron density characterisation'.format(sum(building_mask)), True)
                pandda.log('Datasets to be used for density characterisation: {}'.format(','.join(['\n\t'*(not i%5)+d.tag for i,d in enumerate(pandda.datasets.mask(mask_name=building_mask_name))])))

        # ============================================================================>
        # Select the datasets to analyse
        # ============================================================================>
        pandda.log.bar()
        pandda.log('Selecting Analysis Mask', True)
        pandda.log.bar()
        analysis_mask_name, analysis_mask = pandda.select_datasets_for_analysis(high_res_large_cutoff=cut_resolution, high_res_small_cutoff=high_shell_limit)
        # Check that there're some datasets to analyse
        if sum(analysis_mask) == 0:
            pandda.log('No datasets to analyse at this resolution', True)
            pandda.log.bar()
            continue
        else:
            pandda.log('{!s} datasets selected for contrasting against the characterised electron density'.format(sum(analysis_mask)), True)
            pandda.log('Datasets to be analysed: {}'.format(','.join(['\n\t'*(not i%5)+d.tag for i,d in enumerate(pandda.datasets.mask(mask_name=analysis_mask_name))])))

        # ============================================================================>
        # Combine the masks as we will need to load maps for all datasets
        # ============================================================================>
        map_load_mask = pandda.datasets.all_masks().combine_masks(names=[analysis_mask_name, building_mask_name], invert=False, operation='or', invert_output=False)
        map_load_mask_name = 'Loading @ {!s}A'.format(cut_resolution)
        pandda.datasets.all_masks().add_mask(name=map_load_mask_name, values=map_load_mask.values)

        pandda.log.bar()
        pandda.log('Datasets for Building, Loading, and Analysis')
        pandda.log('Building: {!s} datasets'.format(sum(building_mask)))
        pandda.log('Analysis: {!s} datasets'.format(sum(analysis_mask)))
        pandda.log('Total No: {!s} datasets'.format(sum(map_load_mask)))
        pandda.log.bar()

        # ============================================================================>
        # Have now commited to processing this resoution - record start time
        # ============================================================================>
        t_loop_start = time.time()
        # ============================================================================>
        pandda.log.bar(True, False)
        pandda.log('{!s}A Analysis Started: {!s}'.format(cut_resolution, time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(t_loop_start))), True)
        pandda.log.bar(False, True)
        # ============================================================================>
        # Update limit (cut resolution becomes top limit in next shell)
        # ============================================================================>
        high_shell_limit = cut_resolution

        # ============================================================================>
        #####
        # Load, truncate and scale diffraction data
        #####
        # ============================================================================>
        # Report
        # ============================================================================>
        pandda.log.heading('Preparing selected datasets for analysis at {!s}A'.format(cut_resolution))

        # ============================================================================>
        #####
        # Truncate data, load and scale maps to reference dataset
        #####
        # ============================================================================>
        # Truncate the data to the same resolution and scale to reference data
        pandda.truncate_diffraction_data(datasets=pandda.datasets.mask(mask_name=map_load_mask_name), res_truncate=cut_resolution)
        # Get resolution range of the truncated data and plot histogram
        if pandda.settings.plot_graphs:
            analyse_graphs.write_truncated_data_plots(pandda     = pandda,
                                                      resolution = cut_resolution,
                                                      datasets   = pandda.datasets.mask(mask_name=map_load_mask_name))
        # ============================================================================>
        # Load the reference map so that we can re-scale the individual maps to this
        # ============================================================================>
        ref_map = pandda.load_reference_map(map_resolution=cut_resolution)

        # ============================================================================>
        #####
        # Perform the statistical density characterisation (if required)
        #####
        # ============================================================================>
        if cut_resolution not in pandda.stat_maps.get_resolutions():
            # ============================================================================>
            # Report
            # ============================================================================>
            pandda.log.heading('Performing statistical density characterisation at {!s}A'.format(cut_resolution))
            pandda.log.bar()
            pandda.log('Building density distributions from {!s} datasets at {!s}A'.format(pandda.datasets.size(mask_name=building_mask_name), cut_resolution), True)
            # ============================================================================>
            # Load maps for characterisation datasets
            # ============================================================================>
            map_holder_list = pandda.load_and_morph_maps(datasets       = pandda.datasets.mask(mask_name=building_mask_name),
                                                         ref_map        = ref_map,
                                                         map_resolution = cut_resolution)
            # ============================================================================>
            # Label the datasets
            # ============================================================================>
            map_holder_list.all_masks().add_mask(name=building_mask_name, values=True)
            map_holder_list.all_masks().add_mask(name=analysis_mask_name, values=False)
            for m in map_holder_list.all():
                if pandda.datasets.all_masks().get_value(name=analysis_mask_name, id=m.meta.tag):
                    map_holder_list.all_masks().set_value(name  = analysis_mask_name,
                                                          id    = m.meta.tag,
                                                          value = True)
            # ============================================================================>
            # Create an object to hold all of the maps, and can be used to calculate the mean maps, etc...
            # ============================================================================>
            map_analyser = PanddaMapAnalyser(   dataset_maps     = map_holder_list,
                                                meta             = Meta({'resolution'    : cut_resolution,
                                                                         'map_data_size' : pandda.grid.global_mask().outer_mask_indices().size()}),
                                                statistical_maps = None,
                                                parent           = pandda,
                                                log              = pandda.log   )
            # ============================================================================>
            # Add the meta data about which datasets are being used
            # ============================================================================>
            map_analyser.statistical_maps.meta.characterisation_datasets = [d.tag for d in pandda.datasets.mask(mask_name=building_mask_name)]

            # ============================================================================>
            # Calculate the mean map
            # ============================================================================>
            map_analyser.calculate_mean_map()
            # ============================================================================>
            # If only mean map requested, output and exit
            # ============================================================================>
            if pandda.args.exit_flags.calculate_first_mean_map_only:
                NativeMapMaker.process(dataset  = pandda.datasets.reference(),
                                       map      = map_analyser.statistical_maps['mean_map'],
                                       filename = pandda.file_manager.get_file('mean_map').format(cut_resolution),
                                       args     = pandda.args,
                                       verbose  = pandda.settings.verbose)
                raise SystemExit('Calculating first mean map only: Exiting')
            # ============================================================================>
            # Plot the reference dataset map against the mean map (sanity check)
            # ============================================================================>
            analyse_graphs.write_map_analyser_reference_dataset_graphs(pandda=pandda, map_analyser=map_analyser)
            # ============================================================================>
            # Calculate the uncertainty of all loaded maps (needs the mean map to have been calculated)
            # ============================================================================>
            inner_mask_idxs_reindx = pandda.grid.index_on_other(query=pandda.grid.global_mask().inner_mask_indices(), other=pandda.grid.global_mask().outer_mask_indices())
            map_analyser.calculate_map_uncertainties(masked_idxs=inner_mask_idxs_reindx, cpus=pandda.settings.cpus)
            # ============================================================================>
            # Plot uncertainties of maps
            # ============================================================================>
            try:
                from ascii_graph import Pyasciigraph
                g=Pyasciigraph(float_format='{0:.3f}')
                graph_data = [(m.meta.tag, round(m.meta.map_uncertainty,3)) for m in map_analyser.dataset_maps.all()]
                pandda.log.bar()
                for l in g.graph(label='Uncertainties of maps used for density characterisation', data=graph_data):
                    if l.startswith('#######'): continue
                    pandda.log(l.replace(u"\u2588", '=').replace('= ','> '), True)
                pandda.log.bar(False, True)
            except ImportError:
                print('IMPORT ERROR (ascii_graph) - cannot plot uncertainties')
            except:
                raise
            # ============================================================================>
            # Calculate the statistics of the maps
            # ============================================================================>
            map_analyser.calculate_statistical_maps(cpus=pandda.settings.cpus)
            # ============================================================================>
            # Write the statistical maps
            # ============================================================================>
            pandda.write_map_analyser_maps(map_analyser=map_analyser)
            # ============================================================================>
            # Store in main pandda object and pickle the statistical map objects
            # ============================================================================>
            pandda.log('Storing statistical maps @ {!s}A in main pandda object'.format(cut_resolution))
            pandda.stat_maps.add(stat_map_list=map_analyser.statistical_maps, resolution=cut_resolution, overwrite=False)
            pandda.pickle_the_pandda(components=['stat_maps'])
            # ============================================================================>
            # Pickle map analyser if requested
            # ============================================================================>
            if pandda.args.output.pickling.pickle_map_analysers:
                pandda.pickle(pickle_file=pandda.pickle_handler.get_file('map_analyser').format(cut_resolution), pickle_object=map_analyser, overwrite=True)

            # TODO TODO TODO
            # ============================================================================>
            # Perform analysis of the characterised maps
            # ============================================================================>
            # map_analyser.analyse_point_distributions()
            # TODO TODO TODO

            # ============================================================================>
            #####
            # Load the remaining maps and merge with those already loaded
            #####
            # ============================================================================>
            # Report
            # ============================================================================>
            pandda.log.heading('Loading remaining maps for analysis at {!s}A'.format(cut_resolution))
            # ============================================================================>
            # Select analysis datasets that are not building datasets and load them
            # ============================================================================>
            remain_mask = pandda.datasets.all_masks().combine_masks(names=[analysis_mask_name, building_mask_name],
                                                                    invert=[False,True],
                                                                    operation='and',
                                                                    invert_output=False)
            # Load these datasets
            map_holder_list = pandda.load_and_morph_maps(datasets       = pandda.datasets.mask(mask=remain_mask),
                                                         ref_map        = ref_map,
                                                         map_resolution = cut_resolution)
            # ============================================================================>
            # Extract the maps that are also to be analysed
            # ============================================================================>
            # Extract these maps from the loaded maps
            transfer_maps = map_analyser.dataset_maps.mask(mask_name=analysis_mask_name)
            # Append to the maps already loaded
            map_holder_list.add(transfer_maps)
            # ============================================================================>
            # Garbage collect to remove any remaining un-needed datasets
            # ============================================================================>
            map_analyser.parent = None; del map_analyser; map_analyser = None
            gc.collect(); gc.collect(); gc.collect();
        else:
            # ============================================================================>
            # Load maps for analysis
            # ============================================================================>
            map_holder_list = pandda.load_and_morph_maps(datasets       = pandda.datasets.mask(mask_name=analysis_mask_name),
                                                         ref_map        = ref_map,
                                                         map_resolution = cut_resolution)

        # ============================================================================>
        #####
        # Prepare datasets for comparison against the characterised electron density
        #####
        # ============================================================================>
        # Report
        # ============================================================================>
        pandda.log.heading('Preparing to analyse {!s} datasets at {!s}A'.format(pandda.datasets.size(mask_name=analysis_mask_name), cut_resolution))
        # ============================================================================>
        # Create new map analyser for analysis datasets
        # ============================================================================>
        map_analyser = PanddaMapAnalyser(   dataset_maps     = map_holder_list,
                                            meta             = Meta({'resolution'    : cut_resolution,
                                                                     'map_data_size' : pandda.grid.global_mask().outer_mask_indices().size()}),
                                            statistical_maps = pandda.stat_maps.get(cut_resolution),
                                            parent           = pandda,
                                            log              = pandda.log   )
        # ============================================================================>
        # Calculate the uncertainty of all loaded maps
        # ============================================================================>
        inner_mask_idxs_reindx = pandda.grid.index_on_other(query=pandda.grid.global_mask().inner_mask_indices(), other=pandda.grid.global_mask().outer_mask_indices())
        map_analyser.calculate_map_uncertainties(masked_idxs=inner_mask_idxs_reindx, cpus=pandda.settings.cpus)
        # ============================================================================>
        # Create a new "dummy" map analyser for the parallel steps
        # ============================================================================>
        dummy_map_analyser = PanddaMapAnalyser(dataset_maps = None,
                                               meta         = Meta({'resolution':cut_resolution,
                                                                    'map_data_size':pandda.grid.global_mask().outer_mask_indices().size()}),
                                               statistical_maps = pandda.stat_maps.get(cut_resolution),
                                               parent=None, log=None)
        # ============================================================================>
        # Print blob-search object settings
        # ============================================================================>
        dummy_blob_finder = PanddaZMapAnalyser(params=pandda.params.blob_search, grid=pandda.grid, log=pandda.log)
        dummy_blob_finder.print_settings()

        # ============================================================================>
        #####
        # Create Z-Maps processor objects
        #####
        # ============================================================================>
        t_anal_start = time.time()
        # List of objects for parallel calculation
        dataset_processor_list = []
        # Iterate through and prepare to calculate z-maps
        for i_d, dataset in enumerate(pandda.datasets.mask(mask_name=analysis_mask_name)):
            # ============================================================================>
            # Check/update dataset records/variables
            # ============================================================================>
            # Dataset should not have any events
            assert dataset.events == []
            # Record which resolution this dataset was analysed at
            resolution_count.setdefault(cut_resolution, []).append(dataset.tag)
            # Update datasets masks flag - for this analysis
            pandda.datasets.all_masks().set_value(name='analysed', id=dataset.tag, value=True)
            # Update the dataset meta object -- this is persistent
            dataset.meta.analysed = True
            # ============================================================================>
            # Extract the map for this dataset
            # ============================================================================>
            dataset_map = map_analyser.dataset_maps.get(tag=dataset.tag)
            # ============================================================================>
            # Compile arguments for this datasets
            # ============================================================================>
            dp = DatasetProcessor(dataset=dataset, dataset_map=dataset_map.make_sparse(),
                                  grid=pandda.grid, map_analyser=dummy_map_analyser,
                                  args=pandda.args, verbose=pandda.settings.verbose)
            dataset_processor_list.append(dp)

        # ===========================================================================>
        #####
        # Process the Z-maps for each dataset
        #####
        # ============================================================================>
        # Run processing in parallel
        # ============================================================================>
        pandda.log.heading('Calculating Z-maps and event maps for {!s} dataset(s) at {!s}A'.format(pandda.datasets.size(mask_name=analysis_mask_name), cut_resolution))
        proc_results = easy_mp.pool_map(func=wrapper_run, args=dataset_processor_list, processes=pandda.settings.cpus, chunksize=1)
        # ============================================================================>
        # Update the main thread objects with the results
        # ============================================================================>
        pandda.log.bar()
        pandda.log('Updating with results from analysing {!s} Dataset(s) at {!s}A'.format(pandda.datasets.size(mask_name=analysis_mask_name), cut_resolution), True)
        for results in proc_results:
            # ============================================================================>
            # Unpack results
            # ============================================================================>
            tmp_dataset, dataset_meta, log_strs = results
            # ============================================================================>
            # Report
            # ============================================================================>
            pandda.log('', True)
            pandda.log('======================================>>>', True)
            pandda.log('Z-Map Analysis Results for {}'.format(tmp_dataset.tag), True)
            pandda.log('======================================>>>', True)
            pandda.log('\n'.join(log_strs), True)
            # ============================================================================>
            # Store analysis data in dataset map table
            # ============================================================================>
            pandda.tables.dataset_map_info.set_value(tmp_dataset.tag, 'analysed_resolution', dataset_meta.resolution)
            pandda.tables.dataset_map_info.set_value(tmp_dataset.tag, 'map_uncertainty',     round(dataset_meta.map_uncertainty,3))
            pandda.tables.dataset_map_info.set_value(tmp_dataset.tag, 'obs_map_mean',        round(dataset_meta.obs_map_mean,3))
            pandda.tables.dataset_map_info.set_value(tmp_dataset.tag, 'obs_map_rms',         round(dataset_meta.obs_map_rms,3))
            pandda.tables.dataset_map_info.set_value(tmp_dataset.tag, 'scl_map_mean',        round(dataset_meta.scl_map_mean,3))
            pandda.tables.dataset_map_info.set_value(tmp_dataset.tag, 'scl_map_rms',         round(dataset_meta.scl_map_rms,3))
            pandda.tables.dataset_map_info.set_value(tmp_dataset.tag, 'z_map_mean',          round(dataset_meta.z_mean,3))
            pandda.tables.dataset_map_info.set_value(tmp_dataset.tag, 'z_map_std',           round(dataset_meta.z_stdv,3))
            pandda.tables.dataset_map_info.set_value(tmp_dataset.tag, 'z_map_skew',          round(dataset_meta.z_skew,3))
            pandda.tables.dataset_map_info.set_value(tmp_dataset.tag, 'z_map_kurt',          round(dataset_meta.z_kurt,3))
            # ============================================================================>
            # Write out dataset information to csv file
            # ============================================================================>
            out_list = pandda.tables.dataset_info.loc[tmp_dataset.tag].append(pandda.tables.dataset_map_info.loc[tmp_dataset.tag])
            out_list.to_csv(path=tmp_dataset.file_manager.get_file('dataset_info'), header=True, index_label='dtag')
            # ============================================================================>
            # Mark as interesting and add events to the event table
            # ============================================================================>
            if tmp_dataset.events:
                pandda.datasets.all_masks().set_value(name='interesting', id=tmp_dataset.tag, value=True)
                for e in tmp_dataset.events:
                    pandda.add_event_to_event_table(dataset=tmp_dataset, event=e)
            # ============================================================================>
            # Update the master copy of the dataset object
            # ============================================================================>
            master_dataset = pandda.datasets.get(tag=tmp_dataset.tag)
            master_dataset.events = tmp_dataset.events

        # ============================================================================>
        #####
        # Generate native-aligned maps (in the crystallographic unit cell)
        #####
        # ============================================================================>
        native_map_maker_list = []
        for i_d, dataset in enumerate(pandda.datasets.mask(mask_name=analysis_mask_name)):
            # ============================================================================>
            # Extract the map from the map analyser and store as child
            # ============================================================================>
            dataset.child = map_analyser.dataset_maps.get(tag=dataset.tag)
            # ============================================================================>
            # Make Z-map for each dataset (if events or write_z_maps_for_all_datasets)
            # ============================================================================>
            if (pandda.args.output.maps.write_z_maps=='interesting' and dataset.events) or (pandda.args.output.maps.write_z_maps=='all'):
                ref_z_map = map_analyser.calculate_z_map(map         = dataset.child,
                                                         uncertainty = dataset.child.meta.map_uncertainty,
                                                         method      = pandda.args.params.z_map.map_type)
                ref_z_map = ref_z_map.normalised_copy()
                map_maker = NativeMapMaker(dataset  = dataset,
                                           map      = ref_z_map,
                                           filename = dataset.file_manager.get_file('native_z_map'),
                                           args     = pandda.args,
                                           verbose  = pandda.settings.verbose)
                native_map_maker_list.append(map_maker)
            # ============================================================================>
            # Make dataset map for each dataset (if write_dataset_maps)
            # ============================================================================>
            if (pandda.args.output.maps.write_dataset_map=='interesting' and dataset.events) or (pandda.args.output.maps.write_dataset_map=='all'):
                map_maker = NativeMapMaker(dataset  = dataset,
                                           map      = dataset.child,
                                           filename = dataset.file_manager.get_file('native_map'),
                                           args     = pandda.args,
                                           verbose  = pandda.settings.verbose)
                native_map_maker_list.append(map_maker)
            # ============================================================================>
            # Make Event-map for each event
            # ============================================================================>
            for i,e in enumerate(dataset.events):
                ref_event_map = ( dataset.child - map_analyser.statistical_maps.mean_map * e.info.estimated_bdc )
                map_maker = NativeMapMaker(dataset  = dataset,
                                           map      = ref_event_map,
                                           filename = dataset.file_manager.get_file('native_event_map').format(e.id[1], 1-e.info.estimated_bdc),
                                           args     = pandda.args,
                                           verbose  = pandda.settings.verbose)
                native_map_maker_list.append(map_maker)
            # ============================================================================>
            # Mean (ground-state) map for this resolution
            # ============================================================================>
            if (pandda.args.output.maps.write_mean_map=='interesting' and dataset.events) or (pandda.args.output.maps.write_mean_map=='all'):
                map_maker = NativeMapMaker(dataset  = dataset,
                                           map      = map_analyser.statistical_maps.mean_map,
                                           filename = dataset.file_manager.get_file('native_mean_map'),
                                           args     = pandda.args,
                                           verbose  = pandda.settings.verbose)
                native_map_maker_list.append(map_maker)
        # ============================================================================>
        # Statistical maps in the native frame of the reference dataset
        # ============================================================================>
        if pandda.args.output.maps.write_statistical_maps:
            for m in ['mean_map', 'medn_map','stds_map','sadj_map','skew_map','kurt_map','bimo_map']:
                map_maker = NativeMapMaker(dataset  = pandda.datasets.reference(),
                                           map      = map_analyser.statistical_maps[m],
                                           filename = pandda.file_manager.get_file(m).format(cut_resolution),
                                           args     = pandda.args,
                                           verbose  = pandda.settings.verbose)
                native_map_maker_list.append(map_maker)
        # ============================================================================>
        # Write the compiled list of maps
        # ============================================================================>
        if native_map_maker_list:
            pandda.log.heading('Outputting {!s} native maps at {!s}A'.format(len(native_map_maker_list), cut_resolution))
            t_map_make_start = time.time()
            for n in native_map_maker_list:
                pandda.log('{:<30} -> {}'.format(n.data[0].tag, os.path.split(n.data[2])[-1]), True)
            proc_results = easy_mp.pool_map(func=wrapper_run, args=native_map_maker_list, processes=pandda.settings.cpus)
            t_map_make_end = time.time()
            pandda.log.bar(True, False)
            pandda.log('Map Generation Time: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(t_map_make_end-t_map_make_start))), True)

        # ============================================================================>
        # Ensure all data is synced between the datasets and the central tables
        # ============================================================================>
        pandda.sync_datasets(datasets=pandda.datasets.mask(mask_name=analysis_mask_name), overwrite_dataset_meta=True)
        # ============================================================================>
        #####
        # Write map_analyser summary
        #####
        # ============================================================================>
        if pandda.settings.plot_graphs:
            pandda.log.heading('Writing Statistical Maps Summary'.format(cut_resolution))
            analyse_graphs.write_map_analyser_graphs(pandda=pandda, resolution=cut_resolution,
                                                     analysis_mask_name=analysis_mask_name,
                                                     building_mask_name=building_mask_name)
            analyse_html.write_map_analyser_html(pandda=pandda, resolution=cut_resolution)

        # ============================================================================>
        #####
        # Print summaries
        #####
        # ============================================================================>
        t_loop_end = time.time()
        # ==================================================>
        pandda.log.heading('Finished Analysis at {!s}A'.format(cut_resolution))
        pandda.log.bar()
        pandda.log('{!s}A Z-Map Processing Time: {!s}'.format(cut_resolution, time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(t_loop_end - t_anal_start))), True)
        pandda.log('{!s}A Total Processing Time: {!s}'.format(cut_resolution, time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(t_loop_end - t_loop_start))), True)
        pandda.log.bar()
        pandda.log('@{!s}A:\t {!s}/{!s} New Datasets Analysed'.format(cut_resolution,
                                                                      pandda.datasets.all_masks().combine_masks(
                                                                                                    names=[analysis_mask_name, 'old datasets'],
                                                                                                    invert=[False,True],
                                                                                                    operation='and',
                                                                                                    invert_output=False).sum(),
                                                                      pandda.datasets.size(mask_name='valid - new')))
        pandda.log('@{!s}A:\t {!s}/{!s} Old Datasets Analysed'.format(cut_resolution,
                                                                      pandda.datasets.all_masks().combine_masks(
                                                                                                    names=[analysis_mask_name, 'old datasets'],
                                                                                                    invert=[False,False],
                                                                                                    operation='and',
                                                                                                    invert_output=False).sum(),
                                                                      pandda.datasets.size(mask_name='valid - old')))
        pandda.log('Total:\t {!s}/{!s} Datasets Analysed'.format(pandda.datasets.size(mask_name='analysed'), pandda.datasets.size(mask_name='rejected - total', invert=True)))
        pandda.log.bar()

        # ============================================================================>
        #####
        # PICKLE THE DATASETS THAT HAVE JUST BEEN PROCESSED
        #####
        # ============================================================================>
        # Clear the linked maps if requested
        if not pandda.args.output.pickling.pickle_dataset_maps:
            for dataset in pandda.datasets.mask(mask_name=analysis_mask_name):
                dataset.child = None
        # Pickle the datasets
        pandda.pickle_the_pandda(components=['datasets'], datasets=pandda.datasets.mask(mask_name=analysis_mask_name))

        # ============================================================================>
        # Delete the map analyser to free memory
        # ============================================================================>
        map_analyser.parent = None; del map_analyser; map_analyser = None
        # ============================================================================>
        # Remove the dataset maps from the dataset objects
        # ============================================================================>
        for dataset in pandda.datasets.mask(mask_name=map_load_mask_name):
            dataset.child = None
        # No maps should still be loaded for memory consumption reasons
        for dataset in pandda.datasets.all():
            assert dataset.child is None
        # ============================================================================>
        # Launch garbage-collection manually just to be sure
        # ============================================================================>
        gc.collect(); gc.collect(); gc.collect()

        # ============================================================================>
        #####
        # LIVE ANALYSIS - RUN AT THE END OF EACH LOOP
        #####
        # Analyse the processed data
        # ============================================================================>
        pandda.log.heading('Clustering identified events and updating central data tables')
        # Extract all events from datasets
        all_events=[]; [all_events.extend(d.events) for d in pandda.datasets.all()]
        # Process any identified events and update output files
        if all_events:
            # ==================================================>
            # Cluster events and add to pandda tables
            # ==================================================>
            site_list = pandda.cluster_events_and_update(events=all_events)
            # ==================================================>
            # Update the output analysis files
            # ==================================================>
            pandda.write_output_csvs()
            analyse_html.write_analyse_html(pandda)
            # ==================================================>

    # ============================================================================>
    # Ensure that the collation step happens even if no resolutions are processed
    # ============================================================================>
    if not resolution_count:
        # Extract all events from datasets
        all_events=[]; [all_events.extend(d.events) for d in pandda.datasets.all()]
        # Process any identified events and update output files
        if all_events:
            # ==================================================>
            # Cluster events and add to pandda tables
            # ==================================================>
            site_list = pandda.cluster_events_and_update(events=all_events)
            # ==================================================>
            # Update the output analysis files
            # ==================================================>
            pandda.write_output_csvs()
            analyse_html.write_analyse_html(pandda)
            # ==================================================>

    # Pickle statistical maps again just to be sure
    pandda.pickle_the_pandda(components=['stat_maps'])

    # ============================================================================>
    #####
    # END OF MAIN LOOP
    #####
    # ============================================================================>
    t_analysis_end = time.time()
    # ==================================================>
    pandda.log.heading('ANALYSIS COMPLETE')
    pandda.log('Total Analysis Time: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(t_analysis_end - t_analysis_start))), True)

    return

# ============================================================================>
#
###                 PanDDA Output + Wrap-Up Functions
#
# ============================================================================>

def pandda_end(pandda):
    """Write the output summaries and close"""

    # ============================================================================>
    # Report
    # ============================================================================>
    pandda.log.heading('Writing PanDDA End-of-Analysis Summary')

    # ============================================================================>
    #####
    # Final analysis - run once at the end of the program
    #####
    # ============================================================================>
    # Collate events
    # ============================================================================>
    event_total, event_num, all_dataset_events = pandda.collate_event_counts()
    # ============================================================================>
    # Print some slightly less important information
    # ============================================================================>
    pandda.log.bar()
    for d_tag, event_count in event_num:
        pandda.log('Dataset {!s}: {!s} Events'.format(d_tag, event_count), False)
    # ============================================================================>
    # Print a summary of the number of identified events
    # ============================================================================>
    pandda.log.bar()
    pandda.log('Total Datasets with Events: {!s}'.format(len(event_num)), True)
    pandda.log('Total Events: {!s}'.format(event_total), True)
    pandda.log.bar()
    if event_num:
        pandda.log.bar(True, False)
        pandda.log('Potentially Useful Shortcuts for Future Runs', True)
        pandda.log.bar()
        pandda.log('All datasets with events:')
        pandda.log('exclude_from_characterisation={!s}'.format(','.join(zip(*event_num)[0])), True)
        for site_num, event_list in pandda.tables.event_info.groupby('site_idx').groups.items():
            pandda.log.bar()
            pandda.log('Datasets with events at Site {}'.format(site_num))
            pandda.log('exclude_from_characterisation={!s}'.format(','.join(zip(*event_list)[0])))
        pandda.log.bar()
    # ============================================================================>
    # Record which datasets were used for each statistical map
    # ============================================================================>
    pandda.log.bar(True, False)
    pandda.log('Lists of datasets used to generate the statistical maps', True)
    for res in sorted(pandda.stat_maps.get_resolutions()):
        sm = pandda.stat_maps.get(res)
        pandda.log.bar()
        pandda.log('Statistical Electron Density Characterisation at {}A'.format(res))
        pandda.log('> Density characterised using {} datasets'.format(len(sm.meta.characterisation_datasets)))
        pandda.log('> Dataset IDs: {}'.format(','.join(['\n\t'*(not i%5)+d for i,d in enumerate(sm.meta.characterisation_datasets)])))
    # ============================================================================>
    # Summaries (csvs and html)
    # ============================================================================>
    pandda.log.bar()
    pandda.log('Writing final output files', True)
    pandda.write_output_csvs()
    analyse_html.write_analyse_html(pandda)
    pandda.log.bar()
    # ============================================================================>
    # Screen graphs of number of datasets analysed at each resolution
    # ============================================================================>
    resolution_counts = pandda.tables.dataset_map_info['analysed_resolution'].value_counts().sort_index()
    graph_data = [(str(r), c) for r,c in resolution_counts.iteritems()]
    if graph_data:
        try:
            from ascii_graph import Pyasciigraph
            g=Pyasciigraph()
            for l in g.graph(label='Datasets analysed at each resolution:', data=graph_data):
                if l.startswith('#######'): continue
                pandda.log(l.replace(u"\u2588", '=').replace('= ','> '), True)
        except ImportError: print('IMPORT ERROR (ascii_graph) - CANNOT GENERATE MAP ANALYSIS GRAPH')
        except:             raise
        pandda.log.bar(False, True)
    else:
        pandda.log('> No Resolutions Analysed')

    pandda.log.bar(True, False)
    pandda.log('Datasets Processed: {!s}'.format(sum(resolution_counts)))
    pandda.log('Datasets Loaded {!s}'.format(pandda.datasets.size(mask_name='rejected - total', invert=True)))
    pandda.log.bar(False, True)

    return

# ============================================================================>
#
###                 PanDDA Command-Line Function
#
# ============================================================================>

def pandda_analyse_main(args):
    """Run the PANDDA algorithm, using supplied args"""

    working_phil = extract_params_default(master_phil=pandda_phil, args=args, module_info=module_info)
    welcome()

    try:
        # ============================================================================>
        #####
        # Initialise
        #####
        # ============================================================================
        pandda = PanddaMultiDatasetAnalyser(params=working_phil.extract())
        pandda.run_analysis_init()
        # ============================================================================>
        #####
        # Load and pre-process datasets. Generate grid and grid masks.
        #####
        # ============================================================================>
        pandda_dataset_setup(pandda=pandda)
        pandda_grid_setup(pandda=pandda)
        # ============================================================================>
        #####
        # Run variation analysis on the loaded datasets
        #####
        # ============================================================================>
        #pandda_variation_analysis(pandda=pandda)
        # ============================================================================>
        #####
        # Run the main analysis loop
        #####
        # ============================================================================>
        pandda_main_loop(pandda=pandda)
        # ============================================================================>
        #####
        # Write summaries and post-process
        #####
        # ============================================================================>
        pandda_end(pandda=pandda)
        # ============================================================================>
        #####
        # End
        #####
        # ============================================================================>
    except KeyboardInterrupt:
        raise
    except SystemExit as s:
        # Non-error exit - print type and message
        try:
            pandda.exit(error_msg=traceback.format_exception(type(s), str(s), tb=False))
        except:
            print '<<< Pandda exited before being initialised >>>'
            raise
    except Sorry as s:
        # Recognised error - print type and message
        pandda.exit(error_msg=''.join(traceback.format_exception(type(s), str(s), tb=False)))
    except:
        # Unknown error - print full traceback
        if ("pandda" in locals()) and hasattr(pandda, 'log'):
            pandda.exit(error_msg=traceback.format_exc())
        else:
            raise
    else:
        pandda.exit(error_msg=None)

    return pandda

# ============================================================================>
#
#   COMMAND LINE RUN
#
# ============================================================================>

if __name__ == '__main__':

    pandda = pandda_analyse_main(args=sys.argv[1:])



