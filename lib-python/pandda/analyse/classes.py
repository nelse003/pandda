from __future__ import print_function

import os, sys, glob, time, re
import copy, warnings

import numpy, pandas, scipy.stats

from libtbx import easy_mp
from libtbx.utils import Sorry, Failure
from libtbx.math_utils import ifloor, iceil

import cctbx.maptbx
import scitbx.matrix

from scitbx.array_family import flex

from bamboo.common import Meta, Info
from bamboo.common.logs import Log
from bamboo.common.file import FileManager
from bamboo.common.path import easy_directory, rel_symlink, delete_with_glob
from bamboo.common.status import status_bar, status_bar_2
from bamboo.common.command import CommandManager
from bamboo.common.holders import HolderList

from bamboo.plot import bar
from bamboo.maths import round_no_fail

from giant.dataset import ModelAndData, ElectronDensityMap
from giant.manager import Program
from giant.grid import Grid
from giant.grid.masks import AtomicMask, GridMask
from giant.structure.align import GlobalAlignment
from giant.structure.select import calphas, protein, sel_altloc
from giant.structure.formatting import Labeller, ShortLabeller
from giant.xray.data import estimate_wilson_b_factor
from giant.xray.scaling import IsotropicBfactorScalingFactory

from pandda.phil import pandda_phil
from pandda.analyse.events import cluster_events
from pandda.analyse.functions import DatasetAligner, MapLoader, DensityStatistics, UncertaintyCalculator, wrapper_run
from pandda.constants import *
from pandda import HEADER_TEXT, VERSION


class PanddaDataset(ModelAndData):


    def __init__(self, model, data):
        """Subclass of ModelAndData used for PanDDA Analysis"""

        super(PanddaDataset, self).__init__(model=model, data=data)

        self.child = None
        self.events = []


class PanddaReferenceDataset(ModelAndData):


    _origin_shift = None

    def __init__(self, model, data):

        super(PanddaReferenceDataset, self).__init__(model=model, data=data)

        self.child = None

    def set_origin_shift(self, shift):
        """Creates an alignment corresponding to an origin shift"""

        self._origin_shift = shift
        r = scitbx.matrix.rec([1,0,0,0,1,0,0,0,1], (3,3))
        t = scitbx.matrix.rec(shift, (3,1))
        rt = scitbx.matrix.rt((r,t))
        self.model.alignment = GlobalAlignment(alignment_mx=rt, alignment_sites=calphas(self.model.hierarchy).atoms().extract_xyz(),  id='ref')
        return self.model.alignment

    def origin_shift(self):
        return self._origin_shift

    def nat2grid(self, *args, **kwargs):
        return self.model.alignment.nat2ref(*args, **kwargs)
    def grid2nat(self, *args, **kwargs):
        return self.model.alignment.ref2nat(*args, **kwargs)


class PanddaDatasetList(HolderList):


    _holder_class = PanddaDataset
    _reference_class = PanddaReferenceDataset
    _reference = None

    def _get_num(self, item):
        return item.num
    def _get_tag(self, item):
        return item.tag

    def set_reference(self, dataset):
        assert isinstance(dataset, self._reference_class), 'reference not of the right class: {}'.format(self._reference_class)
        self._reference = dataset

    def reference(self):
        return self._reference


class MapHolderList(HolderList):


    _holder_class = ElectronDensityMap

    def _get_num(self, item):
        return item.meta.num
    def _get_tag(self, item):
        return item.meta.tag


class MapList(Info):


    _map_names = []

    def __init__(self, map_names=None):
        if map_names is None: map_names=[]
        assert self._map_names+map_names, 'No Maps defined'
        for m in self._map_names+map_names:
            self.__dict__[m] = None
        self.meta = Meta()
        self._initialized = True

    def __get_item__(self, item):
        return self.__dict__[item]


class PanddaStatMapList(MapList):


    _map_names = ['mean_map','medn_map','stds_map','sadj_map','skew_map','kurt_map','bimo_map']


class PanddaMultipleStatMapList(object):


    def __init__(self):
        """Store lots of statistical maps"""
        self.map_lists = {}

    def get_resolutions(self):
        return sorted(self.map_lists.keys())

    def add(self, stat_map_list, resolution, overwrite=False):
        assert isinstance(resolution, float), 'Resolution of map must be of type float. Type given: {!s}'.format(type(resolution))
        if overwrite is not False: assert resolution not in self.map_lists.keys(), 'MAPS OF THIS RESOLUTION ALREADY ADDED'
        assert isinstance(stat_map_list, PanddaStatMapList), 'stat_map_list must be of type PanddaMultipleStatMapList. Type given: {!s}'.format(type(stat_map_list))
        self.map_lists[resolution] = stat_map_list

    def get(self, resolution):
        return self.map_lists[resolution]


class PanddaMultiDatasetAnalyser(Program):


    _NAME    = 'pandda.analyse'
    _DESC    = ''
    _VERSION = VERSION
    _TEXT    = HEADER_TEXT

    def __init__(self, params):
        """Main PanDDA Class"""

        # Log init time
        self._init_time = time.time()
        # ===============================================================================>
        # Process input arguments
        # ===============================================================================>
        # Store master phil for reference
        self.master_phil = pandda_phil
        # Store different levels of the params for easy access
        self._input_params = params
        self.args = self._input_params.pandda
        self.params = self._input_params.pandda.params
        self.settings = self._input_params.settings
        # ===============================================================================>
        # Output files stuff
        # ===============================================================================>
        assert self.args.output.out_dir, 'pandda.output.out_dir IS NOT DEFINED'
        self.out_dir = easy_directory(os.path.abspath(self.args.output.out_dir))
        # Create a couple of files and directories for the logs and the parameters
        self._log_dir = easy_directory(os.path.join(self.out_dir, 'logs'))
        self._log_filename = 'pandda-{}.log'.format(time.strftime("%Y-%m-%d-%H%M", time.gmtime()))
        self._def_filename = self._log_filename.replace('.log','.def')
        self._eff_filename = self._log_filename.replace('.log','.eff')
        # Create a log for the object
        self.log = Log(log_file=os.path.join(self._log_dir, self._log_filename), verbose=self.settings.verbose)
        # ===============================================================================>
        # Data and maps stuff
        # ===============================================================================>
        self._new_dataset_files = []
        self.grid = None
        self.datasets = PanddaDatasetList()
        self.pickled_dataset_meta = None
        self.stat_maps = PanddaMultipleStatMapList()
        # ===============================================================================>
        # Analysis objects
        # ===============================================================================>
        # Create tables object (to hold pandas dataframe objects)
        self.tables = Meta()
        # Record global information about the datasets
        self.tables.dataset_info        = pandas.DataFrame( data    = None,
                                                            index   = pandas.Index(data=[], name='dtag'),
                                                            columns = PanddaTableFields.all_dataset_fields      )
        # Record information about the created maps for each dataset
        self.tables.dataset_map_info    = pandas.DataFrame( data    = None,
                                                            index   = pandas.Index(data=[], name='dtag'),
                                                            columns = PanddaTableFields.all_dataset_map_fields  )
        # Record the events detected in each dataset
        self.tables.event_info          = pandas.DataFrame( data    = None,
                                                            index   = pandas.MultiIndex(levels=[[],[]], labels=[[],[]], names=['dtag','event_idx']),
                                                            columns = PanddaTableFields.all_event_fields        )
        # Record information about the clustered events (cluster of events = site)
        self.tables.site_info           = pandas.DataFrame( data    = None,
                                                            index   = pandas.Index(data=[], name='site_idx'),
                                                            columns = PanddaTableFields.all_site_fields         )
        # ===============================================================================>
        # Directory stuff
        # ===============================================================================>
        self._directory_setup()
        self._pickle_setup()

        if os.path.exists(self.pickle_handler.get_file('dataset_meta')):
            self._new_pandda = False
        else:
            self._new_pandda = True

    #########################################################################################################
    #                                                                                                       #
    #                                              Directory Setup                                          #
    #                                                                                                       #
    #########################################################################################################

    def _directory_setup(self):
        """Initialise the pandda directory system"""

        # ===============================================================================>
        # Create a file and directory organiser
        # ===============================================================================>
        fm = self.initialise_file_manager(rootdir=self.out_dir)
        # ===============================================================================>
        # Filename templates
        # ===============================================================================>
        f = PanddaAnalyserFilenames
        h = PanddaHtmlFilenames
        # ===============================================================================>
        # Global Directories that do not change from run to run
        # ===============================================================================>
        fm.add_dir(dir_name='logs',                 dir_tag='logs',                 top_dir_tag='root', create=False, exists=True)
        fm.add_dir(dir_name='processed_datasets',   dir_tag='processed_datasets',   top_dir_tag='root', create=False, exists=False)
        fm.add_dir(dir_name='analyses',             dir_tag='analyses',             top_dir_tag='root', create=False, exists=False)
        fm.add_dir(dir_name='aligned_structures',   dir_tag='aligned_structures',   top_dir_tag='root', create=False, exists=False)
        # ================================================>
        # Input + Status parameters
        # ================================================>
        fm.add_file(file_name='pandda.{}',         file_tag='status',           dir_tag='root')
        fm.add_file(file_name=self._def_filename,  file_tag='def_file',         dir_tag='logs')
        fm.add_file(file_name=self._eff_filename,  file_tag='eff_file',         dir_tag='logs')
        # ================================================>
        # Somewhere to store the analysis summaries - for the user
        # ================================================>
        fm.add_dir(dir_name='html_summaries', dir_tag='output_summaries', top_dir_tag='analyses', create=False, exists=False)
        fm.add_file(file_name=h.initial_html,                      file_tag='initial_html',            dir_tag='output_summaries')
        fm.add_file(file_name=h.map_html,                          file_tag='map_html',                dir_tag='output_summaries')
        fm.add_file(file_name=h.analyse_html,                      file_tag='analyse_html',            dir_tag='output_summaries')
        fm.add_file(file_name=h.analyse_site_graph,                file_tag='analyse_site_graph',      dir_tag='output_summaries')
        fm.add_file(file_name=h.analyse_site_graph_mult,           file_tag='analyse_site_graph_mult', dir_tag='output_summaries')
        fm.add_file(file_name=h.pymol_sites_py,                    file_tag='pymol_sites_py',          dir_tag='output_summaries')
        fm.add_file(file_name=h.pymol_sites_pml,                   file_tag='pymol_sites_pml',         dir_tag='output_summaries')
        fm.add_file(file_name=h.pymol_sites_png_1,                 file_tag='pymol_sites_png_1',       dir_tag='output_summaries')
        fm.add_file(file_name=h.pymol_sites_png_2,                 file_tag='pymol_sites_png_2',       dir_tag='output_summaries')
        # ================================================>
        # Dataset summary graphs
        # ================================================>
        fm.add_dir(dir_name='dataset_summary_graphs', dir_tag='d_graphs', top_dir_tag='analyses', create=False, exists=False)
        fm.add_file(file_name='dataset_resolutions.png',           file_tag='d_resolutions',           dir_tag='d_graphs')
        fm.add_file(file_name='dataset_rfactors.png',              file_tag='d_rfactors',              dir_tag='d_graphs')
        fm.add_file(file_name='dataset_global_rmsd_to_ref.png',    file_tag='d_global_rmsd_to_ref',    dir_tag='d_graphs')
        fm.add_file(file_name='dataset_cell_axes.png',             file_tag='d_cell_axes',             dir_tag='d_graphs')
        fm.add_file(file_name='dataset_cell_angles.png',           file_tag='d_cell_angles',           dir_tag='d_graphs')
        fm.add_file(file_name='dataset_cell_volumes.png',          file_tag='d_cell_volumes',          dir_tag='d_graphs')
        fm.add_file(file_name='dataset_unscaled_wilson_plots.png', file_tag='d_unscaled_wilson_plots', dir_tag='d_graphs')
        fm.add_file(file_name='dataset_unscaled_wilson_rmsds.png', file_tag='d_unscaled_wilson_rmsds', dir_tag='d_graphs')
        fm.add_file(file_name='dataset_scaled_wilson_plots.png',   file_tag='d_scaled_wilson_plots',   dir_tag='d_graphs')
        fm.add_file(file_name='dataset_scaled_wilson_rmsds.png',   file_tag='d_scaled_wilson_rmsds',   dir_tag='d_graphs')
        # ================================================>
        # Map analysis graphs
        # ================================================>
        fm.add_dir(dir_name='map_analysis_summary_graphs', dir_tag='m_graphs', top_dir_tag='analyses', create=False, exists=False)
        fm.add_file(file_name='{}A-truncated_data_resolutions.png', file_tag='dataset_res_hist',       dir_tag='m_graphs')
        fm.add_file(file_name='{}A-truncated_data_wilson_plot.png', file_tag='dataset_wilson_plot',    dir_tag='m_graphs')
        fm.add_file(file_name='{}A-reference_map_distribution.png', file_tag='ref_map_dist',           dir_tag='m_graphs')
        fm.add_file(file_name='{}A-reference_v_mean_unsorted.png',  file_tag='ref_v_mean_map_unsort',  dir_tag='m_graphs')
        fm.add_file(file_name='{}A-reference_v_mean_sorted.png',    file_tag='ref_v_mean_map_sort',    dir_tag='m_graphs')
        fm.add_file(file_name='{}A-maps_mean_medn_histograms.png',  file_tag='map_mean_median_hist',   dir_tag='m_graphs')
        fm.add_file(file_name='{}A-maps_mean_medn_scatter.png',     file_tag='map_mean_median_scat',   dir_tag='m_graphs')
        fm.add_file(file_name='{}A-maps_mean_medn_difference.png',  file_tag='map_mean_median_diff',   dir_tag='m_graphs')
        fm.add_file(file_name='{}A-maps_stds_sadj_histograms.png',  file_tag='map_stds_sadj_hist',     dir_tag='m_graphs')
        fm.add_file(file_name='{}A-maps_stds_sadj_scatter.png',     file_tag='map_stds_sadj_scat',     dir_tag='m_graphs')
        fm.add_file(file_name='{}A-dataset_unc_histograms.png',     file_tag='dataset_unc_hist',       dir_tag='m_graphs')
        fm.add_file(file_name='{}A-dataset_res_unc_scatter.png',    file_tag='dataset_res_unc_scat',   dir_tag='m_graphs')
        fm.add_file(file_name='{}A-dataset_res_rfree_scatter.png',  file_tag='dataset_res_rfree_scat', dir_tag='m_graphs')
        fm.add_file(file_name='{}A-dataset_unc_rfree_scatter.png',  file_tag='dataset_unc_rfree_scat', dir_tag='m_graphs')
        fm.add_file(file_name='{}A-zmap_mean_sadj_histograms.png',  file_tag='zmap_mean_sadj_hist',    dir_tag='m_graphs')
        fm.add_file(file_name='{}A-zmap_skew_kurt_scatter.png',     file_tag='zmap_skew_kurt_scat',    dir_tag='m_graphs')
        # ================================================>
        # Dataset information (general values)
        # ================================================>
        fm.add_file(file_name=f.dataset_info,                      file_tag='dataset_info',            dir_tag='analyses')
        fm.add_file(file_name=f.dataset_map_info,                  file_tag='dataset_map_info',        dir_tag='analyses')
        fm.add_file(file_name=f.dataset_combined_info,             file_tag='dataset_combined_info',   dir_tag='analyses')
        fm.add_file(file_name=f.dataset_masks,                     file_tag='dataset_masks',           dir_tag='analyses')
        # ================================================>
        # Dataset information (identified events + sites)
        # ================================================>
        fm.add_file(file_name=f.event_info,                        file_tag='event_info',              dir_tag='analyses')
        fm.add_file(file_name=f.site_info,                         file_tag='site_info',               dir_tag='analyses')
        fm.add_file(file_name='_point_distributions.csv',          file_tag='point_distributions',     dir_tag='analyses')
        # ================================================>
        # Pickled objects
        # ================================================>
        fm.add_dir(dir_name='pickled_data', dir_tag='pickles', top_dir_tag='root', create=False, exists=False)
        # ===============================================================================>
        # Reference Structure Files (should only be needed once for writing and then only for reloading)
        # ===============================================================================>
        fm.add_dir(dir_name='reference', dir_tag='reference', top_dir_tag='root', create=False, exists=False)
        fm.add_file(file_name=f.reference_structure,               file_tag='reference_structure', dir_tag='reference')
        fm.add_file(file_name=f.reference_dataset,                 file_tag='reference_dataset',   dir_tag='reference')
        fm.add_file(file_name=f.reference_on_origin,               file_tag='reference_on_origin', dir_tag='reference')
        fm.add_file(file_name=f.reference_symmetry,                file_tag='reference_symmetry',  dir_tag='reference')
        fm.add_file(file_name='grid-voronoi-{}.ccp4',              file_tag='grid_voronoi',        dir_tag='reference')
        # ===============================================================================>
        # Standard template files that will be populated when needed
        # ===============================================================================>
        fm.add_dir(dir_name='statistical_maps', dir_tag='statistical_maps', top_dir_tag='reference', create=False, exists=False)
        fm.add_file(file_name=f.mean_map,                          file_tag='mean_map',            dir_tag='statistical_maps')
        fm.add_file(file_name=f.medn_map,                          file_tag='medn_map',            dir_tag='statistical_maps')
        fm.add_file(file_name=f.stds_map,                          file_tag='stds_map',            dir_tag='statistical_maps')
        fm.add_file(file_name=f.sadj_map,                          file_tag='sadj_map',            dir_tag='statistical_maps')
        fm.add_file(file_name=f.skew_map,                          file_tag='skew_map',            dir_tag='statistical_maps')
        fm.add_file(file_name=f.kurt_map,                          file_tag='kurt_map',            dir_tag='statistical_maps')
        fm.add_file(file_name=f.bimo_map,                          file_tag='bimo_map',            dir_tag='statistical_maps')

    def _pickle_setup(self):
        """Initialise all of the pickle filenames"""

        # Pickle Handler
        self.pickle_handler = FileManager(rootdir=self.file_manager.get_dir('pickles'))
        self.pickle_handler.add_file(file_name='grid.pickle',              file_tag='grid')
        self.pickle_handler.add_file(file_name='reference_dataset.pickle', file_tag='reference_dataset')
        self.pickle_handler.add_file(file_name='dataset_masks.pickle',     file_tag='dataset_masks')
        self.pickle_handler.add_file(file_name='dataset_meta.pickle',      file_tag='dataset_meta')
        self.pickle_handler.add_file(file_name='statistical_maps.pickle',  file_tag='stat_maps')
        self.pickle_handler.add_file(file_name='map_analyser_{}A.pickle',  file_tag='map_analyser')
        self.pickle_handler.add_file(file_name='my_pandda.pickle',         file_tag='my_pandda')

    #########################################################################################################
    #                                                                                                       #
    #                                             Utility Functions                                         #
    #                                                                                                       #
    #########################################################################################################

    def is_new_pandda(self):
        """Is this the first time the program has been run?"""
        return self._new_pandda

    def pickle_the_pandda(self, components=None, all=False, datasets=None):
        """Pickles it's major components for quick loading..."""

        self.log.subheading('Pickling PanDDA objects')

        self.log.bar()
        if all == True:
            self.log('Pickling all components of PanDDA', True)
        elif not components:
            self.log('Pickling NOTHING', True)
            return
        else:
            self.log('Selective Pickling: {!s}'.format(', '.join(components).upper()), True)

        if all or ('grid' in components):
            self.log.bar()
            if self.grid is not None:
                self.log('Pickling Map Grid')
                self.pickle(pickle_file=self.pickle_handler.get_file('grid'),
                            pickle_object=self.grid, overwrite=False)
            else:
                self.log('No Reference Grid to Pickle')

        if all or ('datasets' in components):
            self.log.bar()
            if self.datasets.reference():
                self.log('Pickling Reference Dataset')
                self.pickle(pickle_file=self.pickle_handler.get_file('reference_dataset'),
                            pickle_object=self.datasets.reference().get_pickle_copy(), overwrite=True)
            # If no datasets given, pickle them all
            if not datasets:
                datasets = self.datasets.all()
            # Pickle the datasets (individual pickle files)
            if datasets:
                self.log('Pickling Datasets')
                for d in datasets:
                    self.pickle(pickle_file=d.file_manager.get_file('dataset_pickle'),
                                pickle_object=d.get_pickle_copy(), overwrite=True)
            else:
                self.log('No Datasets to Pickle')

        if all or ('stat_maps' in components):
            self.log.bar()
            if self.stat_maps is not None:
                self.log('Pickling Statistical Maps')
                self.pickle(pickle_file   = self.pickle_handler.get_file('stat_maps'),
                            pickle_object = self.stat_maps,
                            overwrite = True)
            else:
                self.log('No Statistical Maps to Pickle')

    def write_parameters_to_parameter_file(self, f_name):
        """Write the current parameter object to file f_name - overwrites current file"""

        # Write the used parameters to file
        with open(f_name, 'w') as out_file:
            out_file.write( '\n'.join([ '# Command Line Args',
                                        '# ', # The line below assures that strings are quoted for easy copy-pasting
                                        '# '+self._NAME+' '+' '.join(sys.argv[1:]).replace(' ','" ').replace('=','="')+'"',
                                        '',
                                        self.master_phil.format(python_object=self._input_params).as_str() ]))

    def exit(self, error_msg=None):
        """Exit the PANDDA, record runtime etc..."""

        self._finish_time = time.time()
        self.log('Runtime: {!s}'.format(time.strftime("%H hours:%M minutes:%S seconds", time.gmtime(self._finish_time - self._init_time))))

        # ===============================================================================>
        # If error, don't make meta or self pickle
        # ===============================================================================>
        if error_msg is not None:
            self.update_status('errored')
            self.log.heading('PanDDA exited with an error')
            self.log.bar(False, True)
            self.log(str(error_msg).strip())
            self.log.bar(True, False)
            self.log.heading('PanDDA exited with an error')
            sys.exit()
        else:
            self.update_status('done')
            self.log.heading('.. FINISHED! PANDDA EXITED NORMALLY ..')
        # ===============================================================================>
        # Write list of datasets loaded into the PanDDA
        # ===============================================================================>
        self.log.bar()
        self.log('Writing dataset information for any future runs', True)
        try:
            if self.pickled_dataset_meta and (not self.args.flags.reload_existing_datasets):
                # Combine with existing meta
                self.log('Combining old dataset meta with new meta for pickle')
                number_of_datasets  = self.pickled_dataset_meta.number_of_datasets  + self.datasets.size()
                dataset_labels      = self.pickled_dataset_meta.dataset_labels      + [d.tag for d in self.datasets.all()]
                dataset_pickle_list = self.pickled_dataset_meta.dataset_pickle_list + [os.path.relpath(d.file_manager.get_file('dataset_pickle'), start=self.out_dir) for d in self.datasets.all()]
            else:
                # Create new meta containing only information about current datasets
                self.log('Creating new meta for pickle')
                number_of_datasets  = self.datasets.size()
                dataset_labels      = [d.tag for d in self.datasets.all()]
                dataset_pickle_list = [os.path.relpath(d.file_manager.get_file('dataset_pickle'), start=self.out_dir) for d in self.datasets.all()]
            # Create a dictionary to be stored
            dataset_meta = Meta({'number_of_datasets'    : number_of_datasets,
                                 'dataset_labels'        : dataset_labels,
                                 'dataset_pickle_list'   : dataset_pickle_list})
            # Pickle the list of locations of the dataset pickles
            self.pickle(pickle_file=self.pickle_handler.get_file('dataset_meta'), pickle_object=dataset_meta, overwrite=True)
        except:
            self.log('FAILED TO PICKLE META')
        # ===============================================================================>
        # Lastly, pickle the main pandda object
        # ===============================================================================>
        try:
            if self.args.output.pickling.pickle_complete_pandda:
                self.log.bar()
                self.log('Pickling the whole PANDDA object (for developer access)')
                self.pickle(pickle_file=self.pickle_handler.get_file('my_pandda'), pickle_object=self, overwrite=True)
        except:
            self.log('FAILED TO PICKLE MYSELF')

        self.log('', True)

        sys.exit()

    #########################################################################################################
    #                                                                                                       #
    #                                       Parameter Validation/Updating                                   #
    #                                                                                                       #
    #########################################################################################################

    def apply_shortcuts(self):
        """Apply default settings where shortcuts are provided"""

        p = self.args

        if not (p.shortcuts.run_characterisation_for_all_resolutions or \
                p.shortcuts.run_in_single_dataset_mode):
            return

        self.log.heading('Applying shortcuts to parameters')

        # ================================================>
        # Seeding mode
        # ================================================>
        if p.shortcuts.run_characterisation_for_all_resolutions:
            self.log('Applying shortcut: run_characterisation_for_all_resolutions')

        # ================================================>
        # Seeding mode
        # ================================================>
        if p.shortcuts.run_in_single_dataset_mode:
            self.log('Applying shortcut: run_in_single_dataset_mode')
            p.params.analysis.min_build_datasets = 1
            p.params.analysis.max_build_datasets = 1

    def check_and_update_input_parameters(self):
        """Validate and preprocess the loaded parameters"""

        self.log.heading('Validating and updating input files')

        p = self.args

        # ================================================>
        # Validate input files and labelling
        # ================================================>
        if p.input.data_dirs is None: raise Sorry('No value has been given for pandda.input.data_dirs')
        if p.input.pdb_style is None: raise Sorry('No value has been given for pandda.input.pdb_style')
        if not (p.input.regex.pdb_regex or \
                p.input.regex.mtz_regex or \
                p.input.regex.dir_regex or \
                (p.input.pdb_style and ('*' in p.input.pdb_style)) or \
                (p.input.mtz_style and ('*' in p.input.mtz_style)) or \
                (p.input.data_dirs and ('*' in p.input.data_dirs))):
            raise Sorry('No method has been provided for labelling the datasets. Need to provide either pdb_regex or mtz_regex or dir_regex, or contain a * in one of pdb_style, mtz_style or data_dirs.')
        # ================================================>
        # Print information
        # ================================================>
        if p.input.mtz_style is None:
            p.input.mtz_style = p.input.pdb_style.replace('.pdb','.mtz')
            self.log('mtz_style not provided - using pdb style: {} -> {}'.format(p.input.pdb_style, p.input.mtz_style))
        if p.input.lig_style is None:
            self.log('pandda.input.lig_style has not been provided - ligand files will not be detected')
        # ================================================>
        # Make fullpaths so we can run on the eff file from anywhere and change directories without worrying about relative paths
        # ================================================>
        p.input.data_dirs = os.path.abspath(p.input.data_dirs)
        p.output.out_dir  = os.path.abspath(p.output.out_dir)
        if p.input.filter.pdb:
            p.input.filter.pdb = os.path.abspath(p.input.filter.pdb)
        if p.input.reference.pdb:
            p.input.reference.pdb = os.path.abspath(p.input.reference.pdb)
        if p.input.reference.mtz:
            p.input.reference.mtz = os.path.abspath(p.input.reference.mtz)
        # ================================================>
        # Trim any unnecessary '/'
        # ================================================>
        if p.input.pdb_style:
            p.input.pdb_style = p.input.pdb_style.strip('/')
        if p.input.mtz_style:
            p.input.mtz_style = p.input.mtz_style.strip('/')
        if p.input.lig_style:
            p.input.lig_style = p.input.lig_style.strip('/')
        # ================================================>
        # Whether or not to create a new ouput directory
        # ================================================>
        if self.is_new_pandda() or p.flags.reprocess_existing_datasets:
            self.log.bar()
            self.log('Setting output.new_analysis_dir = True')
            self.log('A new output analysis directory will be created')
            p.output.new_analysis_dir = True
        # ================================================>
        # Check input parameters are valid values
        # ================================================>
        if p.params.analysis.min_build_datasets < 1:
            raise Sorry('params.analysis.min_build_datasets must be greater than 0. Current value: {}'.format(p.params.analysis.min_build_datasets))
        if p.params.analysis.max_build_datasets < 1:
            raise Sorry('params.analysis.max_build_datasets must be greater than 0. Current value: {}'.format(p.params.analysis.max_build_datasets))
        if p.params.analysis.max_build_datasets < p.params.analysis.min_build_datasets:
            raise Sorry('params.analysis.max_build_datasets must be larger than (or the same size as) params.analysis.min_build_datasets. Currently params.analysis.min_build_datasets={} and params.analysis.max_build_datasets={}'.format(p.params.analysis.min_build_datasets, p.params.analysis.max_build_datasets))
        # ================================================>
        # Use default structure factors if none provided
        # ================================================>
        if not p.params.maps.structure_factors:
            self.log.bar()
            self.log('maps.structure_factors has not been provided.')
            p.params.maps.structure_factors.append('FWT,PHWT')
            p.params.maps.structure_factors.append('2FOFCWT_fill,PH2FOFCWT_fill')
            p.params.maps.structure_factors.append('2FOFCWT,PH2FOFCWT')
            self.log('Will look for the default PHENIX and REFMAC structure factor columns: {}'.format(' or '.join(p.params.maps.structure_factors)))
            self.log('(2FOFCWT_fill,PH2FOFCWT_fill are the "filled" phenix structure factors with missing values "filled" with DFc - if these are present it is prefereable to use these)')
        # ================================================>
        # Hardcoded limits
        # ================================================>
        # Require p.params.analysis.high_res_lower_limit to be set to less than 4A (0-4A)
        if p.params.analysis.high_res_lower_limit > 4.0:
            raise Sorry('params.analysis.high_res_lower_limit must be better/smaller than 4A. This is due to the need to scale the diffraction data '\
                        'correctly - only datasets with a better resolution than this can be analysed. Parameter params.analysis.high_res_lower_limit '\
                        'must be set to a value lower than 4.')
        # Check that the upper limit is smaller than the lower limit
        if p.params.analysis.high_res_upper_limit > p.params.analysis.high_res_lower_limit:
            raise Sorry('params.analysis.high_res_lower_limit must be a larger number than params.analysis.high_res_upper_limit. '+\
                        'Current values are \n\tparams.analysis.high_res_lower_limit={} '.format(p.params.analysis.high_res_lower_limit)+\
                        '\nand \n\tparams.analysis.high_res_upper_limit={}'.format(p.params.analysis.high_res_upper_limit))

        self.log.bar()

    def check_and_update_program_flags(self):
        """Validate and update program control flags"""

        self.log.heading('Validating and updating program flags')

        p = self.args

        # ================================================>
        # Check program flags for conflicts
        # ================================================>
        assert p.params.filtering.flags.same_space_group_only, 'PanDDA does not currently support the analysis of different spacegroups - sorry to get your hopes up...'
        # ================================================>
        # Check if matplotlib is load-able
        # ================================================>
        if self.settings.plot_graphs:
            self.log.bar()
            if self.check_for_matplotlib(backend=self.settings.plotting.backend, interactive=False):
                pass
            else:
                self.settings.plot_graphs = False
        # ================================================>
        # Check availability of existing statistical maps
        # ================================================>
        if not self.stat_maps.get_resolutions():
            self.log.bar()
            self.log('No statistical maps have been loaded from previous runs')
            self.log('Setting method.recalculate_statistical_maps to "Yes"')
            p.flags.recalculate_statistical_maps = "Yes"
        elif (p.flags.recalculate_statistical_maps == 'Yes') and self.stat_maps.get_resolutions():
            raise Sorry('method.recalculate_statistical_maps is set to Yes, but statistical maps have been reloaded when they should not have been.')
        # ================================================>
        # If any datasets are set to be reprocessed, reload all datasets (need to change this to allow for "reload_selected_datasets")
        # ================================================>
        if p.flags.reprocess_existing_datasets or p.flags.reprocess_selected_datasets:
            self.log.bar()
            self.log('method.reprocess_existing_datasets or p.flags.reprocess_selected_datasets are set to True')
            self.log('Old (previously processed) datasets will be reloaded')
            self.log('Setting method.reload_existing_datasets = True')
            p.flags.reload_existing_datasets = True
        # ================================================>
        # Developer flags
        # ================================================>
        if p.output.developer.write_all:
            self.log.bar()
            self.log('output.developer.write_all is set to True: updating developer flags')
            p.output.developer.write_reference_frame_maps = True
            p.output.developer.write_reference_frame_grid_masks = True
            p.output.developer.write_reference_frame_all_z_map_types = True
            self.log('output.developer.write_reference_frame_maps = {}'.format(p.output.developer.write_reference_frame_maps))
            self.log('output.developer.write_reference_frame_grid_masks = {}'.format(p.output.developer.write_reference_frame_grid_masks))
            self.log('output.developer.write_reference_frame_all_z_map_types = {}'.format(p.output.developer.write_reference_frame_all_z_map_types))

    def check_number_of_datasets(self):
        """Check enough datasets have been loaded for analysis (given the input parameters"""

        if (not self.datasets.all()) and (not self.new_files()):
            raise Sorry('No datasets have been found in data_dirs, and no datasets have been reloaded from previous runs')
        # ================================================>
        # Check to see if we're reusing statistical maps
        # ================================================>
        if (self.args.flags.recalculate_statistical_maps == "No") and self.stat_maps.get_resolutions():
            pass
        # ================================================>
        # Check that enough datasets have been found - before any filtering
        # ================================================>
        elif self.datasets.size()+len(self.new_files()) < self.params.analysis.min_build_datasets:
            self.log.bar()
            self.log('Not enough datasets are available for statistical map characterisation.', True)
            self.log('The minimum number required is controlled by changing analysis.min_build_datasets', True)
            self.log('Number loaded ({!s}) is less than the {!s} currently needed.'.format(self.datasets.size()+len(self.new_files()), self.params.analysis.min_build_datasets), True)
            raise Sorry('Not enough datasets are available for statistical map characterisation')
        # ================================================>
        # Check that enough VALID datasets have been loaded
        # ================================================>
        elif self.datasets.all_masks().has_mask('rejected - total') and (self.datasets.size(mask_name='rejected - total', invert=True) < self.params.analysis.min_build_datasets):
            self.log.bar()
            self.log('After filtering datasets, not enough datasets are available for statistical map characterisation.', True)
            self.log('The minimum number required is controlled by changing analysis.min_build_datasets', True)
            self.log('Number loaded ({!s}) is less than the {!s} currently needed.'.format(self.datasets.size(mask_name='rejected - total', invert=True), self.params.analysis.min_build_datasets), True)
            raise Sorry('Not enough datasets are available for statistical map characterisation')

    #########################################################################################################
    #                                                                                                       #
    #                                               Initialisation                                          #
    #                                                                                                       #
    #########################################################################################################

    def run_analysis_init(self):
        """Set up the pandda for a new analysis (doing this will override links to analyses)"""

        # ================================================>
        # Write the input parameters to file
        # ================================================>
        self.write_parameters_to_parameter_file(f_name=self.file_manager.get_file('def_file'))
        # ================================================>
        # Validate the input parameters
        # ================================================>
        self.check_and_update_input_parameters()
        # ================================================>
        # Create a new analysis directory for analyses/summaries (if required)
        # ================================================>
        if self.args.output.new_analysis_dir or (not os.path.exists(self.file_manager.get_dir('analyses'))):
            analysis_time_name = 'analyses-{}'.format(time.strftime("%Y-%m-%d-%H%M", time.gmtime(self._init_time)))
            analysis_time_path = easy_directory(os.path.join(self.file_manager.get_dir('root'), analysis_time_name))
            analysis_link_path = self.file_manager.get_dir('analyses')
            # Remove old analysis link if it exists and link in the new analysis directory
            if os.path.exists(analysis_link_path) and os.path.islink(analysis_link_path):
                os.unlink(analysis_link_path)
            rel_symlink(orig=analysis_time_path, link=analysis_link_path)
        assert os.path.exists(self.file_manager.get_dir('analyses')), 'Output analysis directory does not exist'
        # ===============================================================================>
        # Update the FileManager to make sure all directories are now created
        # ===============================================================================>
        self.file_manager.check_and_create_directories()
        # ===============================================================================>
        # Report
        # ===============================================================================>
        self.log.heading('Files/module information')
        self.log('Running from: {!s}'.format(sys.argv[0]), True)
        self.log.bar()
        self.log('Reading input from : {!s}'.format(self.args.input.data_dirs), True)
        self.log.bar()
        self.log('Writing output to: {!s}'.format(self.out_dir), True)
        # ===============================================================================>
        # Change into output directory
        # ===============================================================================>
        self.log.bar()
        self.log('Changing into output directory: {}'.format(self.out_dir), True)
        os.chdir(self.out_dir)
        # ===============================================================================>
        # Repopulate pandda from previous runs
        # ===============================================================================>
        # Load any objects from previous runs
        self.log.bar()
        self.log('Checking for existing analyses', True)
        self.load_pickled_objects()
        # Reload reference dataset
        if (not self.datasets.reference()) and os.path.exists(self.file_manager.get_file('reference_structure')) and os.path.exists(self.file_manager.get_file('reference_dataset')):
            self.log.bar()
            self.log('Loading Reference Dataset', True)
            self.load_reference_dataset(ref_pdb=self.file_manager.get_file('reference_structure'), ref_mtz=self.file_manager.get_file('reference_dataset'))
        # ================================================>
        # Validate the input parameters
        # ================================================>
        self.check_and_update_program_flags()
        # ===============================================================================>
        # Report
        # ===============================================================================>
        self.log.heading('Program setup complete - starting analysis.')
        # ===============================================================================>
        # Write the header to the log file and write the updated parameters
        # ===============================================================================>
        self.log(self._TEXT.format(program=self._NAME, description=self._DESC), True)
        self.write_running_parameters_to_log()
        self.write_parameters_to_parameter_file(f_name=self.file_manager.get_file('eff_file'))
        # ===============================================================================>
        # Log the start time and update status
        # ===============================================================================>
        self.log.bar()
        self.log('Analysis Started: {!s}'.format(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime(self._init_time))), True)
        self.log.bar()
        self.update_status('running')

    def load_pickled_objects(self):
        """Loads any pickled objects it finds"""

        self.log.bar()
        self.log('Looking for pickled files from previous runs in: {!s}'.format(os.path.relpath(self.pickle_handler.get_dir('root'))), True)

        # Record whether any pickled objects are loaded
        pickles_found = False

        # ==============================>
        # Load Grid
        # ==============================>
        if os.path.exists(self.pickle_handler.get_file('grid')):
            pickles_found = True
            self.log('-> Loading reference grid')
            self.grid = self.unpickle(self.pickle_handler.get_file('grid'))
        # ==============================>
        # Load Reference Dataset
        # ==============================>
        if os.path.exists(self.pickle_handler.get_file('reference_dataset')):
            pickles_found = True
            self.log('-> Loading reference dataset')
            self.datasets.set_reference(dataset=self.unpickle(self.pickle_handler.get_file('reference_dataset')))
        # ==============================>
        # Load the datasets
        # ==============================>
        if os.path.exists(self.pickle_handler.get_file('dataset_meta')):
            pickles_found = True
            self.log('-> Loading old dataset information (existing datasets)')
            self.pickled_dataset_meta = self.unpickle(self.pickle_handler.get_file('dataset_meta'))
            if self.args.flags.reload_existing_datasets:
                pickled_dataset_list = self.pickled_dataset_meta.dataset_pickle_list
                for filename in pickled_dataset_list:
                    assert os.path.isfile(os.path.join(self.out_dir, filename)), 'File does not exist: {!s}'.format(filename)
                self.log('-> Reloading old datasets')
                self.datasets.add([self.unpickle(os.path.join(self.out_dir,f)) for f in pickled_dataset_list])
            else:
                self.log('-> Not reloading old datasets')
        else:
            # No datasets to load - this must be False
            self.args.flags.reload_existing_datasets = False
            self.log('-> No old datasets found')
        # ==============================>
        # Load Statistical Maps
        # ==============================>
        if os.path.exists(self.pickle_handler.get_file('stat_maps')) and (self.args.flags.recalculate_statistical_maps != "Yes"):
            pickles_found = True
            self.log('-> Loading old statistical maps')
            self.stat_maps = self.unpickle(self.pickle_handler.get_file('stat_maps'))
        # ==============================>
        # Or Report
        # ==============================>
        if not pickles_found:
            self.log('-> No Pickles Found', True)
        self.log.bar()

    def initialise_dataset_masks_and_tables(self):
        """Add blank masks to the mask objects, based on how many datasets have been loaded"""

        # ==============================>
        # Initialise dataset masks
        # ==============================>
        self.log.bar()
        self.log('Initialising Dataset Masks.', True)
        # ==============================>
        # Initialise standard blank masks
        # ==============================>
        for mask_name in PanddaMaskNames.all_mask_names:
            self.datasets.all_masks().add_mask(name=mask_name, values=False)
        # ==============================>
        # Initialise masks for datasets that shouldn't be analysed or used for building
        # ==============================>
        if self.args.input.flags.exclude_from_zmap_analysis:
            no_analyse_tags = self.args.input.flags.exclude_from_zmap_analysis.split(',')
            self.log('Not analysing {!s} Datasets: \n\t{!s}'.format(len(no_analyse_tags), '\n\t'.join(no_analyse_tags)))
            no_analyse_mask = [True if d.tag in no_analyse_tags else False for d in self.datasets.all()]
            self.datasets.all_masks().add_mask(name='exclude_from_zmap_analysis', values=no_analyse_mask, overwrite=True)
        if self.args.input.flags.exclude_from_characterisation:
            no_build_tags = self.args.input.flags.exclude_from_characterisation.split(',')
            self.log('Not building distributions from {!s} Datasets: \n\t{!s}'.format(len(no_build_tags), '\n\t'.join(no_build_tags)))
            no_build_mask = [True if d.tag in no_build_tags else False for d in self.datasets.all()]
            self.datasets.all_masks().add_mask(name='exclude_from_characterisation', values=no_build_mask, overwrite=True)
        # ==============================>
        # Initialise mask for datasets that have been previously pickled ("old" datasets)
        # ==============================>
        self.datasets.all_masks().add_mask(name='old datasets', values=False)
        if self.pickled_dataset_meta and self.args.flags.reload_existing_datasets:
            for tag in self.pickled_dataset_meta.dataset_labels:
                self.datasets.all_masks().set_value(name='old datasets', id=tag, value=True)
            self.log('Considering {!s} datasets as "New Datasets"'.format(self.datasets.size(mask_name='old datasets', invert=True)))
            self.log('Considering {!s} datasets as "Old Datasets"'.format(self.datasets.size(mask_name='old datasets')))
        else:
            self.log('Considering all {!s} datasets as "New Datasets"'.format(self.datasets.size(mask_name='old datasets', invert=True)))
            assert self.datasets.size(mask_name='old datasets', invert=True) == self.datasets.size(), 'Total datasets should be same as total new datasets'
        # ==============================>
        # Initialise datasets log tables
        # ==============================>
        self.log.bar()
        self.log('Initialising dataset data tables.', True)
        self.tables.dataset_info     = self.tables.dataset_info.append(pandas.DataFrame(index=[d.tag for d in self.datasets.all()]), verify_integrity=True)
        self.tables.dataset_map_info = self.tables.dataset_map_info.append(pandas.DataFrame(index=[d.tag for d in self.datasets.all()]), verify_integrity=True)
        # ==============================>
        # Populate the event table with information from old datasets
        # ==============================>
        old_datasets = self.datasets.mask(mask_name='old datasets')
        if (not self.args.flags.reprocess_existing_datasets) and old_datasets:
            self.log('Syncing old dataset information to dataset tables.', True)
            self.sync_datasets(datasets=old_datasets)
            self.log('Syncing old dataset events to output tables.', True)
            for dataset in old_datasets:
                if dataset.events:
                    for e in dataset.events:
                        self.add_event_to_event_table(dataset=dataset, event=e)

    #########################################################################################################
    #                                                                                                       #
    #                                               Grid Functions                                          #
    #                                                                                                       #
    #########################################################################################################

    def create_reference_grid(self, dataset, grid_spacing):
        """Create a grid over the given dataset"""

        self.log.bar()
        self.log('Creating Reference Grid', True)

        # ==============================>
        # Create grid size based on reference dataset atoms and buffer zone
        # ==============================>
        sites_cart = protein(dataset.model.hierarchy).atoms().extract_xyz()
        buffer = self.params.masks.outer_mask + self.params.maps.padding
        grid_min = flex.double([s-buffer for s in sites_cart.min()])
        grid_max = flex.double([s+buffer for s in sites_cart.max()])

        # TODO origin -> grid_min, approx_max -> grid_max TODO

        # ==============================>
        # Create grid object
        # ==============================>
        self.grid = Grid(grid_spacing   = grid_spacing,
                         origin         = (0,0,0),
                         approx_max     = tuple(grid_max-grid_min),
                         verbose        = self.settings.verbose)
        self.log(self.grid.summary())

        return self.grid

    def mask_reference_grid(self, dataset):
        """Create masks for the reference grid based on distances from atoms in the reference structure"""

        self.log.bar()
        self.log('Masking Reference Grid', True)

        # ============================================================================>
        # Get main and neighbouring symmetry copies of the reference structures
        # ============================================================================>
        ref_sites_cart = dataset.model.alignment.nat2ref(protein(dataset.model.hierarchy).atoms().extract_xyz())
        sym_copies = dataset.model.crystal_contacts(distance_cutoff=self.args.params.masks.outer_mask+5, combine_copies=True)
        sym_sites_cart = dataset.model.alignment.nat2ref(protein(sym_copies).atoms().extract_xyz())
        # ============================================================================>
        # Global mask used for removing points in the bulk solvent regions
        # ============================================================================>
        if self.grid.global_mask() is None:
            self.log.bar()
            self.log('Generating Protein Mask')
            global_mask = AtomicMask(parent=self.grid, sites_cart=ref_sites_cart,
                                     max_dist=self.params.masks.outer_mask,
                                     min_dist=self.params.masks.inner_mask)
            self.grid.set_global_mask(global_mask)
        # ============================================================================>
        # Global mask used for removing points close to symmetry copies of the protein
        # ============================================================================>
        if self.grid.symmetry_mask() is None:
            self.log.bar()
            self.log('Generating Symmetry Mask')
            symmetry_mask = GridMask(parent=self.grid, sites_cart=sym_sites_cart,
                                     max_dist=self.params.masks.outer_mask,
                                     min_dist=self.params.masks.inner_mask_symmetry)
            self.grid.set_symmetry_mask(symmetry_mask)
        # ============================================================================>
        # Write masked maps
        # ============================================================================>
        # Write protein masked map
        self.grid.write_indices_as_map(indices=self.grid.global_mask().total_mask_indices(),
                                       f_name=self.file_manager.get_file('reference_dataset').replace('.mtz','.totalmask.ccp4'))
        # Write symmetry masked map
        self.grid.write_indices_as_map(indices=self.grid.symmetry_mask().total_mask_indices(),
                                       f_name=self.file_manager.get_file('reference_dataset').replace('.mtz','.symmask.ccp4'))

        return self.grid

    def partition_reference_grid(self, dataset, altlocs=['','A']):

        self.log.bar()
        self.log('Partitioning Reference Grid', True)

        # ============================================================================>
        # Select the sites for generating the voronoi alignments (calphas)
        # ============================================================================>
        partition_h = calphas(sel_altloc(dataset.model.hierarchy, altlocs=altlocs))
        site_cart_ca = dataset.nat2grid(partition_h.atoms().extract_xyz())
        # ============================================================================>
        # Create voronoi cells based on these atoms
        # ============================================================================>
        t1 = time.time()
        self.grid.create_grid_partition(sites_cart=site_cart_ca)
        self.grid.partition.partition(mask  = self.grid.global_mask(),
                                      cpus  = self.settings.cpus)
        t2 = time.time()
        self.log('> Grid partitioning complete > Time Taken: {!s} seconds'.format(int(t2-t1)))
        # ============================================================================>
        # Print cell-by-cell summary or the partitioning
        # ============================================================================>
        self.log.bar(True, False)
        self.log('Partition Summary:', True)
        self.log.bar(False, True)
        voronoi_counts = dict(zip(*numpy.unique(self.grid.partition.nn_groups, return_counts=True)))
        # Cell-by-Cell summary of the voronoi cells
        self.log.bar()
        self.log('CHN - RES -  RESID  - ATOM - ALT :    VORONOI VOLUME')
        self.log.bar()
        for i_atom, atom in enumerate(partition_h.atoms_with_labels()):
            self.log('{:<3} - {:<3} - {:<7} - {:<4} - {:<3} : {:>10} points'.format(atom.chain_id, atom.resname, atom.resid(), atom.name, atom.altloc, voronoi_counts.get(i_atom,0)))
        self.log.bar()
        self.log('Unpartitioned space: {} points'.format(voronoi_counts.get(-1,0)))
        self.log.bar()
        # Chain-by-chain summary of the voronoi cells
        for c in partition_h.chains():
            self.log('Chain {:1} - {:5} regions - ({:5} residues)'.format(c.id, len(c.atoms()), len(c.residue_groups())))
        self.log.bar()
        self.log('Total: {} regions ({} chains, {} residues)'.format(len(partition_h.atoms()), len(list(partition_h.chains())), len(list(partition_h.residue_groups()))), True)
        # ============================================================================>
        # Write grid summary for developer purposes
        # ============================================================================>
        if self.args.output.developer.write_reference_frame_grid_masks:
            self.log.bar()
            self.log('Writing Voronoi grid masks:', True)
            # Write out the un-partitioned section of the grid
            self.grid.write_array_as_map(array=(self.grid.partition.nn_groups==-1),
                                         f_name=self.file_manager.get_file('grid_voronoi').format('unpartitioned'))
            # Write out the voronoi masks for each atom
            for i_cell in range(0,10):
                self.grid.write_array_as_map(array=((self.grid.partition.nn_groups%10)==i_cell)*(self.grid.partition.nn_groups>=0).astype(int),
                                             f_name=self.file_manager.get_file('grid_voronoi').format('{:04}'.format(i_cell)))
            # Write out pymol script to allow results to be access easily
            from bamboo.pymol import PythonScript, Sphere
            pml = PythonScript()
            pml.set_normalise_maps(False)
            for i_atom, atom in enumerate(partition_h.atoms_with_labels()):
                ca_xyz = site_cart_ca[i_atom]
                ca_sph = Sphere(centre=tuple(ca_xyz), radius=min(5.0,max(0.1,0.0002*voronoi_counts.get(i_atom,0))))
                ca_nam = 'voronoi_centres'
                pml.add_shape(shape=ca_sph, obj=ca_nam)
                pml.colour(obj=ca_nam, colour='white')
            pdb = pml.load_pdb(f_name=self.file_manager.get_file('reference_on_origin'))
            pml.colour(obj=pdb, colour='grey')
            for grid_file in sorted(glob.glob(self.file_manager.get_file('grid_voronoi').format('*'))):
                if not grid_file.endswith('.ccp4'): continue
                map_name = pml.load_map(f_name=grid_file)
                mes_name = pml.make_mesh(map_name, contour_level=0)
                pml.colour(mes_name)
            pml.write_script(f_name=self.file_manager.get_file('grid_voronoi').format('cell-centres').replace('.ccp4','.py'))

        return self.grid

    #########################################################################################################
    #                                                                                                       #
    #                                 Dataset loading/alignment/filtering                                   #
    #                                                                                                       #
    #########################################################################################################

    def build_input_list(self):
        """Builds a list of input files from the command line arguments passed"""

        self.log.bar()
        self.log('Building List of Datasets')

        # ==============================>
        # Extract input styles from parameter object
        # ==============================>
        dir_style = self.args.input.data_dirs.strip('/')
        pdb_style = self.args.input.pdb_style.strip('/')
        mtz_style = self.args.input.mtz_style.strip('/')
        self.log('Looking for folders that match {}'.format(dir_style))
        self.log('...and for pdb files that match "{}" in each folder'.format(pdb_style))
        self.log('...and for mtz files that match "{}" in each folder'.format(mtz_style))
        # ==============================>
        # Find datasets in the input directories
        # ==============================>
        new_files = []
        empty_directories = []
        for dir in sorted(glob.glob(self.args.input.data_dirs)):
            pdb_files = [f for f in glob.glob(os.path.join(dir, pdb_style)) if os.path.exists(f)]
            mtz_files = [f for f in glob.glob(os.path.join(dir, mtz_style)) if os.path.exists(f)]
            if not (pdb_files and mtz_files):
                print('EMPTY DIRECTORY: {!s}'.format(dir))
                empty_directories.append(dir)
            elif not pdb_files:
                print('NO PDB IN DIRECTORY: {!s}'.format(dir))
                empty_directories.append(dir)
            elif not mtz_files:
                print('NO MTZ IN DIRECTORY: {!s}'.format(dir))
                empty_directories.append(dir)
            else:
                assert len(pdb_files) == 1, 'More than one matching PDB file found: {!s}'.format(os.path.join(dir, pdb_style))
                assert len(mtz_files) == 1, 'More than one matching MTZ file found: {!s}'.format(os.path.join(dir, mtz_style))
                # ==============================>
                # Found PDB anf MTZ file in directory
                # ==============================>
                new_pdb = pdb_files[0]
                new_mtz = mtz_files[0]
                dataset_tag = [None]
                # ==============================>
                # Regex Matching - PDB file
                # ==============================>
                if '*' in pdb_style:
                    pdb_base = os.path.basename(new_pdb)
                    if self.args.input.regex.pdb_regex:
                        pdb_regex = self.args.input.regex.pdb_regex
                    else:
                        pdb_regex = pdb_style.replace('*', '(.*)')
                    pdb_tag = re.findall(pdb_regex, pdb_base)
                    assert pdb_tag, 'NO PDB TAG FOUND: {!s} -> {!s}'.format(pdb_regex, pdb_base)
                    if isinstance(pdb_tag[0], tuple):
                        self.log('More than one PDB TAG found - choosing the first one of {!s}'.format(pdb_tag[0]))
                        pdb_tag = list(pdb_tag[0])[0:1]
                else: pdb_regex = pdb_tag = None
                # ==============================>
                # Regex Matching - MTZ file
                # ==============================>
                if '*' in mtz_style:
                    mtz_base = os.path.basename(new_mtz)
                    if self.args.input.regex.mtz_regex:
                        mtz_regex = self.args.input.regex.mtz_regex
                    else:
                        mtz_regex = mtz_style.replace('*', '(.*)')
                    mtz_tag = re.findall(mtz_regex, mtz_base)
                    assert mtz_tag, 'NO MTZ TAG FOUND: {!s} -> {!s}'.format(mtz_regex, mtz_base)
                    if isinstance(mtz_tag[0], tuple):
                        self.log('More than one MTZ TAG found - choosing the first one of {!s}'.format(mtz_tag[0]))
                        mtz_tag = list(mtz_tag[0])[0:1]
                else: mtz_regex = mtz_tag = None
                # ==============================>
                # Regex Matching - Directory
                # ==============================>
                if '*' in dir_style:
                    dir_base = os.path.dirname(pdb_files[0])
                    if self.args.input.regex.dir_regex:
                        dir_regex = self.args.input.regex.dir_regex
                    else:
                        dir_regex = dir_style.replace('*', '(.*)')
                    dir_tag = re.findall(dir_regex, dir_base)
                    assert dir_tag, 'NO DIR TAG FOUND: {!s} -> {!s}'.format(dir_regex, dir_base)
                    if isinstance(dir_tag[0], tuple):
                        self.log('More than one DIR TAG found - choosing the first one of {!s}'.format(dir_tag[0]))
                        dir_tag = list(dir_tag[0])[0:1]
                else: dir_regex = dir_tag = None
                # ==============================>
                # Check consistency
                # ==============================>
                if pdb_tag and mtz_tag: assert pdb_tag == mtz_tag, 'PDB-MTZ TAGS ARE NOT IDENTICAL: {} != {}'.format(pdb_tag, mtz_tag)
                if dir_tag and pdb_tag: assert dir_tag == pdb_tag, 'DIR-PDB TAGS ARE NOT IDENTICAL: {} != {}'.format(dir_tag, pdb_tag)
                if dir_tag and mtz_tag: assert dir_tag == mtz_tag, 'DIR-MTZ TAGS ARE NOT IDENTICAL: {} != {}'.format(dir_tag, mtz_tag)
                # ==============================>
                # Extract tag
                # ==============================>
                if   dir_tag: dataset_tag = dir_tag
                elif pdb_tag: dataset_tag = pdb_tag
                elif mtz_tag: dataset_tag = mtz_tag
                # ==============================>
                # Add prefix - slightly obsoleted
                # ==============================>
                if isinstance(dataset_tag[0], str): dataset_tag = [self.args.output.dataset_prefix + dataset_tag[0]]
                else:                               assert dataset_tag[0] is None

                new_files.append(pdb_files+mtz_files+dataset_tag)

        # ==============================>
        # Filter out the already added files
        # ==============================>
        if self.pickled_dataset_meta:
            filtered_new_files = []
            for i, (pdb, mtz, tag) in enumerate(new_files):
                if tag in self.pickled_dataset_meta.dataset_labels:
                    self.log('Dataset with this tag has already been loaded: {!s} - Not loading'.format(tag))
                else:
                    filtered_new_files.append(new_files[i])
        else:
            filtered_new_files = new_files
        # ==============================>
        # Filter out manually labelled datasets to ignore
        # ==============================>
        if self.args.input.flags.ignore_datasets:
            ignore_tags = self.args.input.flags.ignore_datasets.split(',')
            self.log('Ignoring {!s} Datasets: {!s}'.format(len(ignore_tags), ', '.join(ignore_tags)))
            re_filtered_new_files = []
            for i, (pdb, mtz, tag) in enumerate(filtered_new_files):
                if tag in ignore_tags:
                    self.log('Ignoring Dataset: {!s}'.format(tag))
                else:
                    re_filtered_new_files.append(filtered_new_files[i])
            filtered_new_files = re_filtered_new_files
        # ==============================>
        # Report number of empty datasets
        # ==============================>
        self.log.bar()
        self.log('{!s} EMPTY DIRECTORIES FOUND:'.format(len(empty_directories)), True)
        self.log.bar()
        for d in empty_directories:
            self.log('Empty Directory: {}'.format(d))
        # ==============================>
        # Report total number of datasets, and total number of new datasets
        # ==============================>
        if self.pickled_dataset_meta:
            num_old = self.pickled_dataset_meta.number_of_datasets
        else:
            num_old = 0
        self.log.bar()
        self.log('{!s} DATASETS FOUND (TOTAL)'.format(len(filtered_new_files)+num_old), True)
        self.log('{!s} DATASETS FOUND (NEW)'.format(len(filtered_new_files)), True)

        return filtered_new_files

    def add_files(self, file_list):
        """Add (pdb, mtz) file pairs to the datasets to be processed"""

        # ==============================>
        # Limit the number of files that can be added at once to a self
        # ==============================>
        if len(file_list) > self.args.input.max_new_datasets:
            self.log.bar()
            self.log('Limiting the number of new datasets that are added to the pandda analysis (controlled by input.max_new_datasets)')
            self.log('Limiting analysis to the first {} of {} datasets'.format(self.args.input.max_new_datasets, len(file_list)))
            file_list = file_list[:self.args.input.max_new_datasets]
            if len(file_list) != self.args.input.max_new_datasets:
                raise Sorry('Something has gone wrong, number of selected datasets ({}) is not equal to the maximum ({})'.format(len(file_list), self.args.input.max_new_datasets))
        # ==============================>
        # Append to input file list
        # ==============================>
        self._new_dataset_files += file_list

        self.log.bar()
        self.log('{!s} Datasets Added'.format(len(file_list)), True)

    def load_new_datasets(self):
        """Read in maps for the input datasets"""

        # ==============================>
        # Report
        # ==============================>
        if not self.datasets.all() and self.is_new_pandda():
            self.log('Adding First Datasets to PanDDA')
        else:
            self.log('Adding more datasets to PanDDA')
            self.log('{!s} previous datasets already loaded'.format(self.datasets.size()))
            self.log('{!s} previous datasets not loaded'.format(self.pickled_dataset_meta.number_of_datasets - self.datasets.size()))
        # ==============================>
        # Counting offset for dataset index
        # ==============================>
        if self.pickled_dataset_meta: n_offset = self.pickled_dataset_meta.number_of_datasets
        else:                         n_offset = 0
        # ==============================>
        # Load datasets in parallel
        # ==============================>
        start = time.time()
        self.log.bar()
        print('Loading Datasets...')
        loaded_datasets = [PanddaDataset.from_file(model_filename=pdb, data_filename=mtz).label(num=num+n_offset, tag=dtag) for num, (pdb, mtz, dtag) in enumerate(self.new_files())]
        finish = time.time()
        self.log('> Adding Datasets > Time Taken: {!s} seconds'.format(int(finish-start)), True)
        self.log.bar()
        # ==============================>
        # Initialise loaded datasets
        # ==============================>
        # Output Path Templates
        f = PanddaDatasetFilenames
        p = PanddaDatasetPNGFilenames
        for dataset in loaded_datasets:
            # ==============================>
            # Intialise the meta for the dataset
            # ==============================>
            dataset.meta.analysed = False
            dataset.meta.dataset_info = None
            dataset.meta.dataset_map_info = None
            # ==============================>
            # Create a file manager object
            # ==============================>
            dataset.initialise_output_directory(dir=os.path.join(self.file_manager.get_dir('processed_datasets'), dataset.tag))
            # ==============================>
            # Main input/output files
            # ==============================>
            dataset.file_manager.add_file(file_name=f.input_model.format(dataset.tag),                        file_tag='input_model'                  )
            dataset.file_manager.add_file(file_name=f.input_data.format(dataset.tag),                         file_tag='input_data'                   )
            dataset.file_manager.add_file(file_name=f.dataset_info.format(dataset.tag),                       file_tag='dataset_info'                 )
            dataset.file_manager.add_file(file_name=f.dataset_log.format(dataset.tag),                        file_tag='dataset_log'                  )
            dataset.file_manager.add_file(file_name=f.z_peaks_csv.format(dataset.tag),                        file_tag='z_peaks_csv')
            dataset.file_manager.add_file(file_name=f.aligned_model.format(dataset.tag),                      file_tag='aligned_model'                )
            dataset.file_manager.add_file(file_name=f.symmetry_copies.format(dataset.tag),                    file_tag='symmetry_copies'              )
            # ==============================>
            # Native (back-rotated/transformed) maps
            # ==============================>
            dataset.file_manager.add_file(file_name=f.native_obs_map.format(dataset.tag),                     file_tag='native_obs_map'       )
            dataset.file_manager.add_file(file_name=f.native_z_map.format(dataset.tag),                       file_tag='native_z_map'         )
            dataset.file_manager.add_file(file_name=f.native_event_map.format(dataset.tag,'{!s}','{!s}'),     file_tag='native_event_map'     )
            dataset.file_manager.add_file(file_name=f.native_mean_map.format(dataset.tag),                    file_tag='native_mean_map'      )
            # ==============================>
            # Map files (in reference frame)
            # ==============================>
            dataset.file_manager.add_file(file_name=f.sampled_map.format(dataset.tag),                        file_tag='sampled_map'                  )
            dataset.file_manager.add_file(file_name=f.mean_diff_map.format(dataset.tag),                      file_tag='mean_diff_map'                )
            dataset.file_manager.add_file(file_name=f.z_map.format(dataset.tag),                              file_tag='z_map'                        )
            dataset.file_manager.add_file(file_name=f.z_map_naive.format(dataset.tag),                        file_tag='z_map_naive'                  )
            dataset.file_manager.add_file(file_name=f.z_map_naive_norm.format(dataset.tag),                   file_tag='z_map_naive_normalised'       )
            dataset.file_manager.add_file(file_name=f.z_map_uncertainty.format(dataset.tag),                  file_tag='z_map_uncertainty'            )
            dataset.file_manager.add_file(file_name=f.z_map_uncertainty_norm.format(dataset.tag),             file_tag='z_map_uncertainty_normalised' )
            dataset.file_manager.add_file(file_name=f.z_map_corrected.format(dataset.tag),                    file_tag='z_map_corrected'              )
            dataset.file_manager.add_file(file_name=f.z_map_corrected_norm.format(dataset.tag),               file_tag='z_map_corrected_normalised'   )
            dataset.file_manager.add_file(file_name=f.event_map.format(dataset.tag, '{!s}', '{!s}'),          file_tag='event_map'                    )
            # ==============================>
            # Miscellaneous masks
            # ==============================>
            dataset.file_manager.add_file(file_name=f.high_z_mask.format(dataset.tag), file_tag='high_z_mask')
            dataset.file_manager.add_file(file_name=f.grid_mask.format(dataset.tag),   file_tag='grid_mask')
            # ==============================>
            # Links to ligand files (if they've been found)
            # ==============================>
            dataset.file_manager.add_dir(dir_name='ligand_files', dir_tag='ligand', top_dir_tag='root')
            # ==============================>
            # Fitted structures when modelled with pandda.inspect
            # ==============================>
            dataset.file_manager.add_dir(dir_name='modelled_structures', dir_tag='models', top_dir_tag='root')
            # ==============================>
            # Graphs
            # ==============================>
            dataset.file_manager.add_dir(dir_name='graphs', dir_tag='graphs', top_dir_tag='root')
            dataset.file_manager.add_file(file_name=p.s_map_png.format(dataset.tag),                          file_tag='s_map_png',                        dir_tag='graphs')
            dataset.file_manager.add_file(file_name=p.d_mean_map_png.format(dataset.tag),                     file_tag='d_mean_map_png',                   dir_tag='graphs')
            dataset.file_manager.add_file(file_name=p.z_map_naive_png.format(dataset.tag),                    file_tag='z_map_naive_png',                  dir_tag='graphs')
            dataset.file_manager.add_file(file_name=p.z_map_naive_norm_png.format(dataset.tag),               file_tag='z_map_naive_normalised_png',       dir_tag='graphs')
            dataset.file_manager.add_file(file_name=p.z_map_uncertainty_png.format(dataset.tag),              file_tag='z_map_uncertainty_png',            dir_tag='graphs')
            dataset.file_manager.add_file(file_name=p.z_map_uncertainty_norm_png.format(dataset.tag),         file_tag='z_map_uncertainty_normalised_png', dir_tag='graphs')
            dataset.file_manager.add_file(file_name=p.z_map_corrected_png.format(dataset.tag),                file_tag='z_map_corrected_png',              dir_tag='graphs')
            dataset.file_manager.add_file(file_name=p.z_map_corrected_norm_png.format(dataset.tag),           file_tag='z_map_corrected_normalised_png',   dir_tag='graphs')
            dataset.file_manager.add_file(file_name=p.z_map_qq_plot_png.format(dataset.tag),                  file_tag='z_map_qq_plot_png',                dir_tag='graphs')
            dataset.file_manager.add_file(file_name=p.bdc_est_png.format(dataset.tag, '{!s}'),                file_tag='bdc_est_png',                      dir_tag='graphs')
            dataset.file_manager.add_file(file_name=p.unc_qqplot_png.format(dataset.tag),                     file_tag='unc_qqplot_png',                   dir_tag='graphs')
            dataset.file_manager.add_file(file_name=p.obs_qqplot_sorted_png.format(dataset.tag),              file_tag='obs_qqplot_sorted_png',            dir_tag='graphs')
            dataset.file_manager.add_file(file_name=p.obs_qqplot_unsorted_png.format(dataset.tag),            file_tag='obs_qqplot_unsorted_png',          dir_tag='graphs')
            # ==============================>
            # Scripts
            # ==============================>
            dataset.file_manager.add_dir(dir_name='scripts', dir_tag='scripts', top_dir_tag='root')
            dataset.file_manager.add_file(file_name=f.pymol_script,     file_tag='pymol_script',    dir_tag='scripts')
            dataset.file_manager.add_file(file_name=f.ccp4mg_script,    file_tag='ccp4mg_script',   dir_tag='scripts')
            # ==============================>
            # Output images
            # ==============================>
            dataset.file_manager.add_dir(dir_name='images',  dir_tag='images', top_dir_tag='root')
#            dataset.file_manager.add_file(file_name=f.ccp4mg_png,       file_tag='ccp4mg_png',      dir_tag='images')
            # ==============================>
            # Pickled objects
            # ==============================>
            dataset.file_manager.add_dir(dir_name='pickles', dir_tag='pickles', top_dir_tag='root')
            dataset.file_manager.add_file(file_name=f.dataset_pickle,   file_tag='dataset_pickle',  dir_tag='pickles')

            # ==============================>
            # Create links to input files
            # ==============================>
            # Links for the dataset input files
            link_pdb = dataset.file_manager.get_file('input_model')
            link_mtz = dataset.file_manager.get_file('input_data')
            # Link the input files to the output folder
            if not os.path.exists(link_pdb): rel_symlink(orig=dataset.model.filename, link=link_pdb)
            if not os.path.exists(link_mtz): rel_symlink(orig=dataset.data.filename, link=link_mtz)
            # ==============================>
            # Search for ligand files and copy them to the output ligands folder
            # ==============================>
            lig_files = glob.glob(os.path.join(os.path.dirname(dataset.model.filename), self.args.input.lig_style))
            for lig_file in lig_files:
                # Find all files with the same basename but allowing for different extensions. Then link to output folder.
                lig_base = os.path.splitext(lig_file)[0] + '.*'
                lig_matches = glob.glob(lig_base)
                for lig in lig_matches:
                    out_path = os.path.join(dataset.file_manager.get_dir('ligand'), os.path.basename(lig))
                    if os.path.exists(lig) and (not os.path.exists(out_path)):
                        try: shutil.copy(lig, out_path)
                        except: pass
            # ==============================>
            # Lastly: Update the pointer to the new path (relative to the pandda directory)
            # ==============================>
            dataset.model.filename = os.path.relpath(link_pdb, start=self.out_dir)
            dataset.data.filename = os.path.relpath(link_mtz, start=self.out_dir)

        self.datasets.add(loaded_datasets)
        self.log('{!s} Datasets Loaded (New).          '.format(len(loaded_datasets), True))
        self.log('{!s} Datasets Loaded (Total).        '.format(self.datasets.size(), True))

        # Reset the list of input files
        self._new_dataset_files = []

    # TODO MOVE TO PANDDA DATASET LIST? TODO
    def select_reference_dataset(self, method='resolution', max_rfree=0.4, min_resolution=5):
        """Select dataset to act as the reference - scaling, aligning etc"""

        assert method in ['resolution','rfree'], 'METHOD FOR SELECTING THE REFERENCE DATASET NOT RECOGNISED: {!s}'.format(method)

        # ==============================>
        # Create a mask of the datasets that can be selected as the reference dataset
        # ==============================>
        potential_reference_mask = self.datasets.all_masks().combine_masks(names=['exclude_from_characterisation', 'rejected - total'], invert_output=True)
        self.datasets.all_masks().add_mask(name='potential reference datasets', values=potential_reference_mask)
        # ==============================>
        # Get the potential reference datasets
        # ==============================>
        filtered_datasets = self.datasets.mask(mask_name='potential reference datasets')
        if not filtered_datasets: raise Failure("Can't select a reference dataset - NO SUITABLE (NON-REJECTED) DATASETS REMAINING")
        # ==============================>
        # Select by either R-free or Resolution
        # ==============================>
        self.log.bar()
        self.log('Selecting Reference Dataset by: {!s}'.format(method), True)
        if method == 'rfree':
            # Get RFrees of datasets (set to dummy value of 999 if resolution is too high so that it is not selected)
            r_frees = [d.model.input.get_r_rfree_sigma().r_free if (d.data.mtz_object().max_min_resolution()[1] < min_resolution) else 999 for d in filtered_datasets]
            if len(r_frees) == 0: raise Exception('NO DATASETS BELOW RESOLUTION CUTOFF {!s}A - CANNOT SELECT REFERENCE DATASET'.format(min_resolution))
            ref_dataset_index = r_frees.index(min(r_frees))
        elif method == 'resolution':
            # Get Resolutions of datasets (set to dummy value of 999 if r-free is too high so that it is not selected)
            resolns = [d.data.mtz_object().max_min_resolution()[1] if (d.model.input.get_r_rfree_sigma().r_free < max_rfree) else 999 for d in filtered_datasets]
            if len(resolns) == 0: raise Exception('NO DATASETS BELOW RFREE CUTOFF {!s} - CANNOT SELECT REFERENCE DATASET'.format(max_rfree))
            ref_dataset_index = resolns.index(min(resolns))
        # ==============================>
        # Report and return
        # ==============================>
        reference = filtered_datasets[ref_dataset_index]
        self.log('Reference Selected: {!s}'.format(reference.tag), True)
        self.log('Resolution: {!s}, RFree: {!s}'.format(reference.data.mtz_object().max_min_resolution()[1], reference.model.input.get_r_rfree_sigma().r_free), True)

        return reference.model.filename, reference.data.filename

    def load_reference_dataset(self, ref_pdb, ref_mtz):
        """Set the reference dataset, to which all other datasets will be aligned and scaled"""

        self.log('---------->>>', True)
        self.log('Loading Reference Dataset: {!s}'.format(ref_mtz), True)

        # ==============================>
        # Output links to reference files
        # ==============================>
        link_ref_pdb = self.file_manager.get_file('reference_structure')
        link_ref_mtz = self.file_manager.get_file('reference_dataset')
        # ==============================>
        # Remove any old links to dataset
        # ==============================>
        if os.path.abspath(ref_pdb) != os.path.abspath(link_ref_pdb):
            if os.path.exists(link_ref_pdb): os.unlink(link_ref_pdb)
            if os.path.exists(link_ref_mtz): os.unlink(link_ref_mtz)
        # ==============================>
        # Create links to dataset
        # ==============================>
        if not os.path.exists(link_ref_pdb): rel_symlink(orig=ref_pdb, link=link_ref_pdb)
        if not os.path.exists(link_ref_mtz): rel_symlink(orig=ref_mtz, link=link_ref_mtz)
        # ==============================>
        # Create and set reference dataset
        # ==============================>
        ref_dataset = PanddaReferenceDataset.from_file(model_filename=os.path.relpath(link_ref_pdb, start=self.out_dir),
                                                       data_filename=os.path.relpath(link_ref_mtz, start=self.out_dir)).label(num=-1, tag='reference')

        # ==============================>
        # Calculate the shift required to move the reference structure into the positive quadrant
        # ==============================>
        buffer = self.params.masks.outer_mask + self.params.maps.padding
        sites_min = protein(ref_dataset.model.hierarchy).atoms().extract_xyz().min()
        ref_dataset.set_origin_shift(shift=-1*flex.double(sites_min)+buffer)
        self.log('Origin Shift for reference structure: {!s}'.format(tuple([round(s,3) for s in ref_dataset.origin_shift()])))
        # ==============================>
        # Set as the reference dataset for the analysis
        # ==============================>
        self.datasets.set_reference(dataset=ref_dataset)
        # ==============================>
        # Write out structure in reference coordinate frames
        # ==============================>
        tmp_r_hierarchy = ref_dataset.model.hierarchy.deep_copy()
        if not os.path.exists(self.file_manager.get_file('reference_on_origin')):
            tmp_r_hierarchy.atoms().set_xyz(ref_dataset.nat2grid(tmp_r_hierarchy.atoms().extract_xyz()))
            tmp_r_hierarchy.write_pdb_file(self.file_manager.get_file('reference_on_origin'))
        if not os.path.exists(self.file_manager.get_file('reference_symmetry')):
            ref_sym_copies = ref_dataset.model.crystal_contacts(distance_cutoff=self.args.params.masks.outer_mask+5, combine_copies=True)
            ref_sym_copies.atoms().set_xyz(ref_dataset.nat2grid(ref_sym_copies.atoms().extract_xyz()))
            ref_sym_copies.write_pdb_file(self.file_manager.get_file('reference_symmetry'))

        # ==============================>
        # Extract reference dataset SFs
        # ==============================>
        sf_cols = [self.args.input.reference.structure_factors.split(',')] if self.args.input.reference.structure_factors \
                    else [sf.split(',') for sf in self.args.params.maps.structure_factors]
        self.log('Checking for any of these structure factors in reference dataset:', True)
        for sf_pair in sf_cols:
            self.log('> {} and {}'.format(*sf_pair), True)
        self.log.bar()
        # Record when a pair is found
        dataset_sfs = None
        # Extract mtz object from the reference dataset
        mtz_obj = ref_dataset.data.mtz_object()
        # Iterate through possible structure factor pairs
        for sf_pair in sf_cols:
            # Check that the data contains the appropriate column
            if mtz_obj.has_column(sf_pair[0]) and mtz_obj.has_column(sf_pair[1]):
                dataset_sfs = ','.join(sf_pair)
                break
        # Raise error if no columns are identified
        if dataset_sfs is None:
            raise Sorry('No matching structure factors were found in the reflection data for reference dataset. \n'+\
                        'Looking for structure factors: \n\t{}\n'.format('\n\t'.join(map(' and '.join,sf_cols))) +\
                        'Structure factors in this dataset: \n\t{}\n'.format('\n\t'.join(mtz_obj.column_labels()))+\
                        'You may need to change the maps.structure_factors or the reference.structure_factors option.')
        # Store column labels for later
        ref_dataset.meta.column_labels = dataset_sfs
        # Load the diffraction data
        self.log('Loading structure factors for reference dataset: {}'.format(dataset_sfs))
        ref_dataset.data.miller_arrays[dataset_sfs] = self.datasets.reference().data.get_structure_factors(columns=dataset_sfs)

        return self.datasets.reference()

    def align_datasets(self, method):
        """Align each structure the reference structure"""

        self.log.heading('Aligning dataset structures')

        assert method in ['local','global'], 'METHOD NOT DEFINED: {!s}'.format(method)

        # ==============================>
        # Select the datasets for alignment
        # ==============================>
        datasets_for_alignment = [d for d in self.datasets.mask(mask_name='rejected - total', invert=True) if not d.model.alignment]
        # ==============================>
        # Report
        # ==============================>
        self.log.bar(True, False)
        if not datasets_for_alignment:
            self.log('All datasets are already aligned/No datasets to align')
            return
        self.log('Generating Alignments (using {!s} cores) for {} datasets'.format(self.settings.cpus, len(datasets_for_alignment)))
        # ==============================>
        # Create a shifted version of the reference model for alignment
        # ==============================>
        ref_model = copy.deepcopy(self.datasets.reference().model)
        ref_model.hierarchy.atoms().set_xyz(ref_model.alignment.nat2ref(ref_model.hierarchy.atoms().extract_xyz()))
        # ==============================>
        # Generate the alignments for each structure
        # ==============================>
        t1 = time.time()
        arg_list = [DatasetAligner(model=d.model, other=ref_model, method=method, id=d.tag) for d in datasets_for_alignment]
        dataset_alignments = easy_mp.pool_map(func=wrapper_run, args=arg_list, processes=self.settings.cpus)
        t2 = time.time()
        # ==============================>
        # Post-process the alignments and output the aligned structures
        # ==============================>
        self.log.bar()
        self.log('Alignment Summaries:', True)
        errors = []
        for dataset, alignment in zip(datasets_for_alignment, dataset_alignments):
            # If errored, print and record
            if isinstance(alignment, str):
                self.log.bar()
                self.log('Failed to align dataset {}'.format(dataset.tag))
                self.log(alignment)
                errors.append((dataset,alignment))
                continue
            # Attach alignment to dataset
            assert dataset.tag == alignment.id
            dataset.model.alignment = alignment
            # Output an aligned copy of the structure
            aligned_struc = dataset.model.hierarchy.deep_copy()
            aligned_struc.atoms().set_xyz(dataset.model.alignment.nat2ref(coordinates=dataset.model.hierarchy.atoms().extract_xyz()))
            aligned_struc.write_pdb_file(file_name=os.path.join(self.file_manager.get_dir('aligned_structures'), '{!s}-aligned.pdb'.format(dataset.tag)))
            aligned_struc.write_pdb_file(file_name=dataset.file_manager.get_file('aligned_model'))
            # Write alignment summary to log
            self.log.bar()
            self.log(dataset.model.alignment.summary())
        t3 = time.time()
        # ==============================>
        # Report Errors
        # ==============================>
        if errors:
            for dataset, error in errors:
                self.log.bar()
                self.log('Failed to align dataset {}'.format(dataset.tag))
                self.log.bar()
                self.log(error)
            raise Sorry('Failed to align {} datasets. Error messages printed above.'.format(len(errors)))
        # ==============================>
        # Report
        # ==============================>
        self.log.bar()
        self.log('> Generating Alignments > Time Taken: {!s} seconds'.format(int(t2-t1)), True)
        self.log('> Aligning Structures   > Time Taken: {!s} seconds'.format(int(t3-t2)), True)
        self.log.bar()

    def collate_dataset_variables(self):
        """Go through all of the datasets and collect lots of different characteristics of the datasets for identifying odd datasets"""

        self.log.bar()
        self.log('Collating Dataset Structure/Crystal Variables', True)

        # ==============================>
        # Extract information about each dataset
        # ==============================>
        for d in self.datasets.all():
            # Resolution info
            self.tables.dataset_info.set_value(d.tag, 'high_resolution', numpy.round(d.data.summary.high_res,3))
            self.tables.dataset_info.set_value(d.tag, 'low_resolution',  numpy.round(d.data.summary.low_res,3))
            # Unit cell info
            self.tables.dataset_info.set_value(d.tag, ['uc_a','uc_b','uc_c','uc_alpha','uc_beta','uc_gamma'],   numpy.round(d.data.summary.unit_cell.parameters(),3))
            self.tables.dataset_info.set_value(d.tag, 'uc_vol',                                                 numpy.round(d.data.summary.unit_cell.volume()),3)
            # Spacegroup info
            self.tables.dataset_info.set_value(d.tag, 'space_group', d.data.summary.space_group.info().type().lookup_symbol())
            # Quality info
            self.tables.dataset_info.set_value(d.tag, 'r_work', round_no_fail(d.model.input.get_r_rfree_sigma().r_work,3))
            self.tables.dataset_info.set_value(d.tag, 'r_free', round_no_fail(d.model.input.get_r_rfree_sigma().r_free,3))
            # Alignment info
            if d.model.alignment:
                self.tables.dataset_info.set_value(d.tag, 'rmsd_to_reference', numpy.round(d.model.alignment.alignment_rmsd(),3))

    def filter_datasets_1(self, filter_dataset=None):
        """Filter out the datasets which contain different protein models (i.e. protein length, sequence, etc)"""

        self.log.subheading('Filtering datasets (before alignment)')

        # ==============================>
        # Print list of possible error classes
        # ==============================>
        self.log.bar(True, False)
        self.log('Potential rejection classes:', True)
        for failure_class in PanddaMaskNames.reject_mask_names:
            self.log('\t{!s}'.format(failure_class), True)
        self.log.bar()
        # ==============================>
        # If no filtering dataset given, filter against the reference dataset
        # ==============================>
        if not filter_dataset: filter_dataset = self.datasets.reference()
        # ==============================>
        # Remove poor/incompatible datasets
        # ==============================>
        for dataset in self.datasets.mask(mask_name='rejected - total', invert=True):
            print('\rFiltering Dataset {!s}          '.format(dataset.tag), end=''); sys.stdout.flush()
            if self.params.filtering.flags.same_space_group_only and (dataset.model.space_group.info().symbol_and_number() != filter_dataset.model.space_group.info().symbol_and_number()):
                self.log('\rRejecting Dataset: {!s}          '.format(dataset.tag))
                self.log('Different Space Group')
                self.log('Reference: {!s}, {!s}: {!s}'.format(filter_dataset.model.space_group.info().symbol_and_number(),
                                                        dataset.tag, dataset.model.space_group.info().symbol_and_number()))
                self.log.bar()
                self.tables.dataset_info.set_value(dataset.tag, 'rejection_reason', 'Different Space Group to Reference')
                self.datasets.all_masks().set_value(name='rejected - different space group', id=dataset.tag, value=True)
            elif dataset.model.input.get_r_rfree_sigma().r_free > self.params.filtering.max_rfree:
                self.log('\rRejecting Dataset: {!s}          '.format(dataset.tag))
                self.log('RFree is higher than cutoff: {!s}'.format(self.params.filtering.max_rfree))
                self.log('High RFree: {!s}'.format(dataset.model.input.get_r_rfree_sigma().r_free))
                self.log.bar()
                self.tables.dataset_info.set_value(dataset.tag, 'rejection_reason', 'R-free is too high')
                self.datasets.all_masks().set_value(name='rejected - rfree', id=dataset.tag, value=True)
            elif self.params.filtering.flags.similar_models_only and (not dataset.model.hierarchy.is_similar_hierarchy(filter_dataset.model.hierarchy)):
                self.log('\rRejecting Dataset: {!s}          '.format(dataset.tag))
                self.log('Non-Identical Structure (Structures do not contain the same atoms)')
                self.log.bar()
                self.tables.dataset_info.set_value(dataset.tag, 'rejection_reason', 'Atoms present in the dataset are different to atoms present in the reference structure')
                self.datasets.all_masks().set_value(name='rejected - non-identical structures', id=dataset.tag, value=True)
            else:
                pass

        self.log('\rDatasets Filtered.                          ', True)
        # ==============================>
        # Update the masks
        # ==============================>
        self.update_masks()
        # ==============================>
        # Report
        # ==============================>
        self.log.bar()
        self.log('Rejected Datasets (Total):     {!s}'.format(sum(self.datasets.all_masks().get_mask(name='rejected - total'))), True)
        self.log.bar()
        # ==============================>
        # Print summary of rejected datasets
        # ==============================>
        reject_reasons = self.tables.dataset_info['rejection_reason'].value_counts().sort_index()
        if reject_reasons.any():
            self.log('Reasons for Rejection:')
            for reason, count in reject_reasons.iteritems():
                self.log('{} Dataset(s) - {}'.format(count, reason))

    def filter_datasets_2(self):
        """Filter out the non-isomorphous datasets"""

        self.log.subheading('Filtering datasets (after alignment).')

        # ==============================>
        # Print list of possible error classes
        # ==============================>
        self.log.bar(True, False)
        self.log('Potential rejection classes:', True)
        for failure_class in PanddaMaskNames.reject_mask_names:
            self.log('\t{!s}'.format(failure_class), True)
        self.log.bar()
        # ==============================>
        # Remove poor/incompatible datasets
        # ==============================>
        for dataset in self.datasets.mask(mask_name='rejected - total', invert=True):
            print('\rFiltering Dataset {!s}          '.format(dataset.tag), end=''); sys.stdout.flush()
            if dataset.model.alignment.alignment_rmsd() > self.params.filtering.max_rmsd_to_reference:
                self.log('\rRejecting Dataset: {!s}          '.format(dataset.tag))
                self.log('Alignment RMSD is too large')
                self.log('Aligned (Calpha) RMSD: {!s}'.format(dataset.model.alignment.alignment_rmsd()))
                self.log('---------->>>', True)
                self.tables.dataset_info.set_value(dataset.tag, 'rejection_reason', 'High RMSD to aligned reference structure')
                self.datasets.all_masks().set_value(name='rejected - rmsd to reference', id=dataset.tag, value=True)
            else:
                pass

        self.log('\rDatasets Filtered.                          ', True)
        # ==============================>
        # Update the masks
        # ==============================>
        self.update_masks()
        # ==============================>
        # Report
        # ==============================>
        self.log.bar()
        self.log('Rejected Datasets (Total):     {!s}'.format(sum(self.datasets.all_masks().get_mask(name='rejected - total'))), True)
        self.log.bar()
        # ==============================>
        # Print summary of rejected datasets
        # ==============================>
        reject_reasons = self.tables.dataset_info['rejection_reason'].value_counts().sort_index()
        if reject_reasons.any():
            self.log('Reasons for Rejection:')
            for reason, count in reject_reasons.iteritems():
                self.log('{} Dataset(s) - {}'.format(count, reason))

    def update_masks(self):
        """Update the combined masks"""

        self.log.bar()
        self.log('Updating dataset masks', True)

        # ==============================>
        # Combine all "rejected" masks into one combined mask
        # ==============================>
        combined_reject_mask = self.datasets.all_masks().combine_masks(names=PanddaMaskNames.reject_mask_names)
        self.datasets.all_masks().add_mask(name='rejected - total', values=combined_reject_mask, overwrite=True)

        # ==============================>
        # Combine the combined "rejected" mask with the "old datasets" mask
        # ==============================>
        valid_datasets_new = self.datasets.all_masks().combine_masks(names=['rejected - total', 'old datasets'],
                                                                     invert=[True, True],
                                                                     operation='and',
                                                                     invert_output=False)
        valid_datasets_old = self.datasets.all_masks().combine_masks(names=['rejected - total', 'old datasets'],
                                                                     invert=[True, False],
                                                                     operation='and',
                                                                     invert_output=False)
        self.datasets.all_masks().add_mask(name='valid - new', values=valid_datasets_new, overwrite=True)
        self.datasets.all_masks().add_mask(name='valid - old', values=valid_datasets_old, overwrite=True)

    #########################################################################################################
    #                                                                                                       #
    #                               Dataset diffraction data/map functions                                  #
    #                                                                                                       #
    #########################################################################################################

    def load_and_scale_diffraction_data(self, datasets=None):
        """Extract amplitudes and phases for creating map"""

        if datasets is None:
            datasets=self.datasets.mask(mask_name='rejected - total', invert=True)

        # ==============================>
        # Prepare objects for scaling
        # ==============================>
        # Extract reference miller array and prepare for scaling
        ref_dataset = self.datasets.reference()
        ref_miller = ref_dataset.data.miller_arrays[ref_dataset.meta.column_labels]
        factory = IsotropicBfactorScalingFactory(reference_miller_array=ref_miller.as_intensity_array())

        # ==============================>
        # Report
        # ==============================>
        t1 = time.time()
        self.log.heading('Loading and scaling structure factors for each dataset')

        for dataset in datasets:
            # Get the columns to be loaded for this dataset
            dataset_sfs = dataset.meta.column_labels
            # ==============================>
            # Load structure factors
            # ==============================>
            if dataset_sfs in dataset.data.miller_arrays.keys():
                self.log('Structure factors already loaded {:<20} - Dataset {!s}'.format(dataset_sfs, dataset.tag))
                ma_unscaled_com = dataset.data.miller_arrays[dataset_sfs]
            else:
                self.log('Loading structure factors {:<20} for dataset {!s}'.format(dataset_sfs, dataset.tag))
                ma_unscaled_com = dataset.data.get_structure_factors(columns=dataset_sfs)

            # ==============================>
            # Compare data to reference - allow spotting of outliers
            # ==============================>
            # Extract miller array and calculate amplitudes
            ma_unscaled_int = ma_unscaled_com.as_intensity_array()
            ma_unscaled_phs = ma_unscaled_com*(1.0/ma_unscaled_com.as_amplitude_array().data())
            # Scale to the reference dataset
            scaling = factory.calculate_scaling(miller_array=ma_unscaled_int)
            # Select high resolution and low resolution separately
            high_res_sel = scaling.x_values > 1/(4.0**2)
            low_res_sel = scaling.x_values <= 1/(4.0**2)
            # Log unscaled rmsds values
            self.tables.dataset_info.set_value(dataset.tag, 'unscaled_wilson_rmsd_all', numpy.round(scaling.unscaled_rmsd,3))
            self.tables.dataset_info.set_value(dataset.tag, 'unscaled_wilson_rmsd_<4A', numpy.round(scaling.rmsd_to_ref(values=scaling.scl_values, sel=high_res_sel),3))
            self.tables.dataset_info.set_value(dataset.tag, 'unscaled_wilson_rmsd_>4A', numpy.round(scaling.rmsd_to_ref(values=scaling.scl_values, sel=low_res_sel),3))
            self.log('RMSD to reference dataset: (before) {}'.format(numpy.round(scaling.unscaled_rmsd,3)))

            # ==============================>
            # Scale data to reference
            # ==============================>
            if not self.args.testing.perform_diffraction_data_scaling:
                # No scaling
                ma_scaled_com = ma_unscaled_com
            else:
                # Report values
                self.log('RMSD to reference dataset: (after)  {}'.format(numpy.round(scaling.scaled_rmsd,3)))
                self.log('Optimised B-factor Scaling Factor: {} ({})'.format(numpy.round(scaling.scaling_b_factor,3), 'sharpened' if scaling.scaling_b_factor<0 else 'blurred'))
                # Log scaled rmsds values
                self.tables.dataset_info.set_value(dataset.tag, 'applied_b_factor_scaling', numpy.round(scaling.scaling_b_factor,3))
                self.tables.dataset_info.set_value(dataset.tag, 'scaled_wilson_rmsd_all',   numpy.round(scaling.scaled_rmsd,3))
                self.tables.dataset_info.set_value(dataset.tag, 'scaled_wilson_rmsd_<4A',   numpy.round(scaling.rmsd_to_ref(values=scaling.out_values, sel=high_res_sel),3))
                self.tables.dataset_info.set_value(dataset.tag, 'scaled_wilson_rmsd_>4A',   numpy.round(scaling.rmsd_to_ref(values=scaling.out_values, sel=low_res_sel),3))
                # Apply scaling to diffraction data
                scaling.new_x_values(x_values=ma_unscaled_int.d_star_sq().data())
                ma_scaled_int = ma_unscaled_int.array(data=scaling.transform(ma_unscaled_int.data())).set_observation_type_xray_intensity()
                ma_scaled_com = ma_unscaled_phs * ma_scaled_int.as_amplitude_array().data()

            assert ma_unscaled_com.is_complex_array()
            assert ma_scaled_com.is_complex_array()
            # Store in the dataset object
            dataset.data.miller_arrays[dataset_sfs] = ma_unscaled_com
            dataset.data.miller_arrays['scaled']    = ma_scaled_com

            # ==============================>
            # Wilson B-factors
            # ==============================>
            dataset.meta.unscaled_wilson_b = estimate_wilson_b_factor(miller_array=ma_unscaled_com)
            dataset.meta.scaled_wilson_b   = estimate_wilson_b_factor(miller_array=ma_scaled_com)
            self.tables.dataset_info.set_value(index=dataset.tag, col='unscaled_wilson_B', value=dataset.meta.unscaled_wilson_b)
            self.tables.dataset_info.set_value(index=dataset.tag, col='scaled_wilson_B',   value=dataset.meta.scaled_wilson_b)

            self.log.bar()
        # ==============================>
        # Report
        # ==============================>
        t2 = time.time()
        self.log('\r> Structure factors extracted > Time taken: {!s} seconds'.format(int(t2-t1))+' '*30, True)
        # ==============================>
        # Update Z-score columns
        # ==============================>
        self.tables.dataset_info['unscaled_wilson_rmsd_all_z'] = scipy.stats.zscore(self.tables.dataset_info['unscaled_wilson_rmsd_all'])
        self.tables.dataset_info['unscaled_wilson_rmsd_<4A_z'] = scipy.stats.zscore(self.tables.dataset_info['unscaled_wilson_rmsd_<4A'])
        self.tables.dataset_info['unscaled_wilson_rmsd_>4A_z'] = scipy.stats.zscore(self.tables.dataset_info['unscaled_wilson_rmsd_>4A'])
        if self.args.testing.perform_diffraction_data_scaling:
            self.tables.dataset_info['scaled_wilson_rmsd_all_z'] = scipy.stats.zscore(self.tables.dataset_info['scaled_wilson_rmsd_all'])
            self.tables.dataset_info['scaled_wilson_rmsd_<4A_z'] = scipy.stats.zscore(self.tables.dataset_info['scaled_wilson_rmsd_<4A'])
            self.tables.dataset_info['scaled_wilson_rmsd_>4A_z'] = scipy.stats.zscore(self.tables.dataset_info['scaled_wilson_rmsd_>4A'])
        # ==============================>
        # Exclude from characterisation if poor quality diffraction
        # ==============================>
        self.log.subheading('Checking for datasets that exhibit large wilson plot RMSDs to the reference dataset')
        pref = '' if self.args.testing.perform_diffraction_data_scaling else 'un'
        for dataset in datasets:
            if (self.tables.dataset_info.get_value(index=dataset.tag, col=pref+'scaled_wilson_rmsd_all_z') > self.params.excluding.max_wilson_plot_rmsd_z_score) or \
               (self.tables.dataset_info.get_value(index=dataset.tag, col=pref+'scaled_wilson_rmsd_<4A_z') > self.params.excluding.max_wilson_plot_rmsd_z_score) or \
               (self.tables.dataset_info.get_value(index=dataset.tag, col=pref+'scaled_wilson_rmsd_>4A_z') > self.params.excluding.max_wilson_plot_rmsd_z_score):
                self.log('Dataset {} has a high wilson plot rmsd to the reference dataset (relative to other datasets) - excluding from characterisation'.format(dataset.tag))
                self.datasets.all_masks().set_value(name='exclude_from_characterisation', id=dataset.tag, value=True)

        return None

    def truncate_diffraction_data(self, datasets, res_truncate):
        """Truncate data at the same indices across all the datasets"""

        self.log.bar()
        self.log.heading('(Scaling and) truncating reflection data')

        # ==============================>
        # Find how many reflections are present in the reference dataset
        # ==============================>
        ref_cols = self.datasets.reference().meta.column_labels
        ref_size = self.datasets.reference().data.miller_arrays[ref_cols].set().size()
        self.log('Number of reflections in reference dataset: {!s}'.format(ref_size))

        # TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
        #
        # NEED TO REMOVE TRUNCATION STEP IF SPACEGROUPS ARE NOT IDENTICAL (OR UNIT CELLS ISOMORPHOUS)
        # 1) ONLY TRUNCATE IF IDENTICAL SPACEGROUP
        # 2) IF NOT IDENTICAL SPACEGROUP ASSERT THAT EACH DATASET HAS A FULL (COMPLETE) SET OF MILLER OF INDICES
        #
        # TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

        # ==============================>
        # Truncate miller indices to the common set (not including the reference dataset)
        # ==============================>
        common_set = datasets[0].data.miller_arrays['scaled'].set()
        for dataset in datasets[1:]:
            common_set = common_set.common_set(dataset.data.miller_arrays['scaled'], assert_is_similar_symmetry=False)
        # ==============================>
        # Report
        # ==============================>
        self.log.bar()
        self.log('Number of Common Reflections between Datasets: {!s} ({!s}% of reference)'.format(common_set.size(), int(100.0*common_set.size()/ref_size)))
        self.log('After Truncation - Reflections per dataset: {!s}'.format(common_set.size()))
        # ==============================>
        # Truncate diffraction data for all of the datasets (including the reference dataset)
        # ==============================>
        self.datasets.reference().data.miller_arrays['truncated'] = self.datasets.reference().data.miller_arrays[ref_cols].common_set(common_set, assert_is_similar_symmetry=False)
        for dataset in datasets:
            dataset.data.miller_arrays['truncated'] = dataset.data.miller_arrays['scaled'].common_set(common_set, assert_is_similar_symmetry=False)

    def load_reference_map(self, map_resolution=0):
        """Load the reference map, and calculate some map statistics"""

        # ==============================>
        # Get the reference dataset
        # ==============================>
        ref_dataset = self.datasets.reference()
        # ==============================>
        # Take the scaled diffraction data for the dataset and create fft
        # ==============================>
        fft_map = ref_dataset.data.miller_arrays['truncated'].fft_map(resolution_factor = self.params.maps.resolution_factor,
                                                                      d_min             = map_resolution,
                                                                      symmetry_flags    = cctbx.maptbx.use_space_group_symmetry)
        # ==============================>
        # Scale the map
        # ==============================>
        if   self.params.maps.scaling == 'none':   pass
        elif self.params.maps.scaling == 'sigma':  fft_map.apply_sigma_scaling()
        elif self.params.maps.scaling == 'volume': fft_map.apply_volume_scaling()
        # ==============================>
        # Transform to the reference frame
        # ==============================>
        # Extract the points for the map (in the grid frame)
        masked_cart = self.grid.grid2cart(self.grid.global_mask().outer_mask(), origin=True)
        # Transform to the frame of the reference dataset (this should become unnecessary in the future)
        masked_cart = self.datasets.reference().grid2nat(masked_cart)
        # Create map handler in the native frame and extract the map values
        ref_map_true = ElectronDensityMap.from_fft_map(fft_map).as_map()
        masked_vals = ref_map_true.get_cart_values(masked_cart)
        # ==============================>
        # Create a new electron density map object for the "grid map"
        # ==============================>
        ref_map = ElectronDensityMap(map_data=masked_vals, unit_cell=self.grid.unit_cell(),
                        map_indices=self.grid.global_mask().outer_mask_indices(),
                        map_size=self.grid.grid_size(), sparse=True)
        # Store the map as a child of the dataset
        ref_dataset.child = ref_map
        # ==============================>
        # Add some meta for debugging, etc
        # ==============================>
        ref_map.meta.type = 'reference-map'
        ref_map.meta.resolution = map_resolution
        ref_map.meta.map_mean = ref_map.get_map_data(sparse=True).min_max_mean().mean
        ref_map.meta.map_rms = ref_map.get_map_data(sparse=True).standard_deviation_of_the_sample()

        return ref_map

    def load_and_morph_maps(self, datasets, ref_map, map_resolution=0):
        """Create map from miller arrays. Transform map into the reference frame by sampling at the given points."""

        assert ref_map.is_sparse(), 'Reference map is not in sparse form'

        # ==============================>
        # Create holder for the output map objects
        # ==============================>
        map_list = MapHolderList()
        # ==============================>
        # Return empty list if no datasets
        # ==============================>
        if not datasets: return map_list

        # ==============================>
        # Report
        # ==============================>
        self.log.bar()
        self.log('Loading {} electron density maps @ {!s}A'.format(len(datasets), map_resolution), True)
        self.log.bar()
        # ==============================>
        # FFT diffraction data into maps
        # ==============================>
        self.log('Converting structure factors to electron density maps')
        start = time.time()
        for dataset in datasets:
            dataset.data.fft_maps['truncated'] = dataset.data.miller_arrays['truncated'].fft_map(resolution_factor = self.params.maps.resolution_factor,
                                                                                                 d_min             = map_resolution,
                                                                                                 symmetry_flags    = cctbx.maptbx.use_space_group_symmetry)
            # Remove diffraction data from dataset object to save memory
            dataset.data.miller_arrays['truncated'] = None
        # ==============================>
        # Report
        # ==============================>
        finish = time.time()
        self.log('> FFT-ing structure factors ({!s} Datasets) > Time Taken: {!s} seconds'.format(len(datasets), int(finish-start)), True)
        # ==============================>
        # Load maps in parallel
        # ==============================>
        self.log.bar()
        print('Loading maps (using {!s} cores)'.format(self.settings.cpus))
        start = time.time()
        arg_list = [MapLoader(dataset=d, grid=self.grid, reference_map=ref_map, args=self.args, verbose=self.settings.verbose) for d in datasets]
        # Print a sort of progress bar
        print('1'+''.join(['{:<5}'.format(i) for i in range(0,len(arg_list)+5,5)])[2:])
        print(' '*len(arg_list)+'|\r', end=''); sys.stdout.flush()
        dataset_maps = easy_mp.pool_map(func=wrapper_run, args=arg_list, processes=self.settings.cpus, chunksize=1)
        print('|')
        map_list.add(dataset_maps)
        # ==============================>
        # Clear fft map data to save memory
        # ==============================>
        for map in map_list.all():
            map_dataset = self.datasets.get(tag=map.meta.tag)
            map_dataset.data.fft_maps['truncated'] = None
        # ==============================>
        # Report
        # ==============================>
        finish = time.time()
        self.log('> Loading maps ({!s} Datasets) > Time Taken: {!s} seconds'.format(map_list.size(), int(finish-start)), True)
        self.log('----------------------------------->>>')

        return map_list

    #########################################################################################################
    #                                                                                                       #
    #                              Dataset utility functions (e.g. reset/sync)                              #
    #                                                                                                       #
    #########################################################################################################

    def new_files(self):
        """Get all of the files that were added on this run"""
        return self._new_dataset_files

    def reset_loaded_datasets(self):
        """Check that pickled datasets are ready for reprocessing, etc, if required"""

        if self.args.flags.reprocess_existing_datasets:   datasets_for_reprocessing = self.datasets.all()
        elif self.args.flags.reprocess_selected_datasets: datasets_for_reprocessing = [self.datasets.get(tag=t) for t in self.args.flags.reprocess_selected_datasets.split(',')]
        else:                                             datasets_for_reprocessing = []

        for dataset in datasets_for_reprocessing:
            # ==============================>
            # Reset the meta objects
            # ==============================>
            dataset.meta.analysed = False
            dataset.meta.dataset_map_info = None
            # ==============================>
            # Delete events from before
            # ==============================>
            dataset.events = []
            # ==============================>
            # Reset the map information for the dataset
            # ==============================>
            self.tables.dataset_map_info.loc[dataset.tag] = numpy.nan

    def check_loaded_datasets(self, datasets):
        """Check that the datasets are analysable (have the right mtz columns, etc)"""

        self.log.bar()
        self.log('Performing checks on the loaded datasets', True)
        self.log.bar()
        # ==============================>
        # Extract structure factor column names
        # ==============================>
        sf_cols = [c.split(',') for c in self.params.maps.structure_factors]
        # Print the columns we're looking for
        self.log('Checking for any of these structure factors in each dataset:', True)
        for sf_pair in sf_cols:
            self.log('> {} and {}'.format(*sf_pair), True)
        self.log.bar()
        # ==============================>
        # Check for these columns in each dataset
        # ==============================>
        for dataset in datasets:
            # Check that the input files exist
            if not os.path.exists(dataset.model.filename):
                raise Sorry('Model file does not exist for dataset {}: {}'.format(dataset.tag, dataset.model.filename))
            if not os.path.exists(dataset.data.filename):
                raise Sorry('Data file does not exist for dataset {}: {}'.format(dataset.tag, dataset.data.filename))
            # Extract reflection data
            mtz_obj = dataset.data.mtz_object()
            # The structure factors to be used in this dataset
            dataset_sfs = None
            # ==============================>
            # Find which structure factors are present in this dataset
            # ==============================>
            for sf_pair in sf_cols:
                # Check that the data contains the appropriate column
                if mtz_obj.has_column(sf_pair[0]) and mtz_obj.has_column(sf_pair[1]):
                    dataset_sfs = sf_pair
                    break
            # Raise error if no columns are identified
            if dataset_sfs is None:
                raise Sorry('No matching structure factors were found in the reflection data for dataset {}. \n'.format(dataset.tag)+\
                            'Looking for structure factors: \n\t{}\n'.format('\n\t'.join(map(' and '.join,sf_cols))) +\
                            'Structure factors in this dataset: \n\t{}\n'.format('\n\t'.join(mtz_obj.column_labels()))+\
                            'You may need to change the maps.structure_factors option.')
            # Store the column labels in the dataset object
            dataset.meta.column_labels = ','.join(dataset_sfs)
            # Report which are being used
            self.log('Validating structure factors for dataset {}: {}'.format(dataset.tag, dataset.meta.column_labels))
            # ==============================>
            # Check the data for the selected columns
            # ==============================>
            for c in dataset_sfs:
                # Extract column from the dataset, and associated miller_set
                col = mtz_obj.get_column(c)
                ms = col.mtz_crystal().miller_set(False)
                # Get the boolean selection for valid reflections from the column
                valid_selection = col.selection_valid()
                # ==============================>
                # Check that all data are valid values
                # ==============================>
                if self.params.checks.all_data_are_valid_values is True:
                    if sum(valid_selection) != len(valid_selection):
                        raise Sorry('Structure factor column "{}" in dataset {} has missing reflections (some values are set to N/A or zero). '.format(c, dataset.tag)+\
                                    '{} reflections have a value of zero. '.format(len(valid_selection)-sum(valid_selection))+\
                                    'You should populate the structure factors for these reflections with their estimated values. '\
                                    'Analysing maps with missing reflections (escepially low resolution reflections!) will degrade the quality of the analysis. '\
                                    'However, you can continue by setting checks.all_data_are_valid_values=None.')
                # ==============================>
                # Check that the data is complete up until a certain resolution
                # ==============================>
                if self.params.checks.low_resolution_completeness is not None:
                    # Find selection for the low resolution reflections
                    low_res_sel = ms.resolution_filter_selection(d_min=self.params.checks.low_resolution_completeness, d_max=999)
                    # Extract a complete set of miller indices up to cutoff resolution
                    ms_c = ms.complete_set(d_min_tolerance=0.0, d_min=self.params.checks.low_resolution_completeness, d_max=999)
                    # Check that there are the same number of reflections in this set as the other set
                    if ms_c.size() != sum(low_res_sel):
                        raise Sorry('Structure factor column "{}" in dataset {} has missing reflections below {}A (some reflections are not present in the miller set). '.format(c, dataset.tag, self.params.checks.low_resolution_completeness)+\
                                    '{} reflections are missing from the miller set in the reflection file. '.format(ms_c.size()-sum(low_res_sel))+\
                                    'You should add these missing reflections and populate the structure factors for these reflections with their estimated values. '\
                                    'Analysing maps with missing reflections (escepially low resolution reflections!) will degrade the quality of the analysis. '\
                                    'I really would not do this, but you can continue by setting checks.low_resolution_completeness to None.')
                    # Calculate overlap between low resolution set and valid set to ensure none missing from low resolution set
                    valid_low_res = valid_selection.select(low_res_sel)
                    # Check if any low resolution reflections are invlaid
                    if sum(valid_low_res) != len(valid_low_res):
                        raise Sorry('Structure factor column "{}" in dataset {} has missing reflections below {}A (some values are set to N/A or zero). '.format(c, dataset.tag, self.params.checks.low_resolution_completeness)+\
                                    '{} reflections have a value of zero. '.format(len(valid_low_res)-sum(valid_low_res))+\
                                    'You should populate the structure factors for these reflections with their estimated values. '\
                                    'Analysing maps with missing reflections (escepially low resolution reflections!) will degrade the quality of the analysis. '\
                                    'I really would not do this, but you can continue by setting checks.low_resolution_completeness=None.')

    def sync_datasets(self, datasets=None, overwrite_dataset_meta=False):
        """Sync the loaded datasets and the pandda dataset tables"""

        if not datasets: datasets = self.datasets.all()

        for dataset in datasets:
            # Copy data from pandda dataset tables to dataset
            if (dataset.meta.dataset_info is None) or overwrite_dataset_meta:
                dataset.meta.dataset_info = self.tables.dataset_info.loc[dataset.tag]
            # Copy data from dataset to pandda dataset tables
            else:
                for col,val in dataset.meta.dataset_info.iteritems():
                    self.tables.dataset_info.set_value(index=dataset.tag, col=col, value=val)

            # Copy data from pandda dataset tables to dataset
            if (dataset.meta.dataset_map_info is None) or overwrite_dataset_meta:
                dataset.meta.dataset_map_info = self.tables.dataset_map_info.loc[dataset.tag]
            # Copy data from dataset to pandda dataset tables
            else:
                for col,val in dataset.meta.dataset_map_info.iteritems():
                    self.tables.dataset_map_info.set_value(index=dataset.tag, col=col, value=val)

    #########################################################################################################
    #                                                                                                       #
    #                                  Dataset variation analysis functions                                 #
    #                                                                                                       #
    #########################################################################################################

    def analyse_unit_cell_variation(self):
        pass

    def analyse_alignment_variation(self):
        """Look at all of the rotation matrices for the local alignments and calculate the rms between neighbours"""

        # TODO TODO TODO
        return
        # TODO TODO TODO

        assert self.params.alignment.method in ['global', 'local']

        if self.params.alignment.method == 'global':
            self.log('GLOBAL ALIGNMENT SELECTED - NOT ANALYSING ROTATION MATRICES')
            return

        if self.settings.plot_graphs:
            import matplotlib
            matplotlib.interactive(False)
            from matplotlib import pyplot

        # Select datasets to analyse
        used_datasets = self.datasets.mask(mask_name='rejected - total', invert=True)

        ref_calpha_atoms = calphas(self.datasets.reference().model.hierarchy).atoms()
        ref_calpha_sites = ref_calpha_atoms.extract_xyz()
        ref_calpha_labels = [make_label(a) for a in ref_calpha_atoms]

        # Array to hold the output data
        num_datasets = len(used_datasets)
        num_pairs =  len(ref_c_alpha_labels)-1
        output_diffs = numpy.zeros((num_datasets, num_pairs, 2))

        # Iterate through the datasets and pull out the alignment matrices
        for d_num, dataset in enumerate(used_datasets):
            # Extract and sort dataset alignments
            alignments = dataset.local_alignment_transforms()
            alignment_keys = sorted(alignments.keys())
            assert alignment_keys == ref_c_alpha_labels

            # Iterate through adjacent pairs of matrices
            for i in range(0, num_pairs):
                # Label and lsq fit for the current calpha
                calpha_1 = alignment_keys[i]
                rt_1 = alignments[calpha_1]
                # And for the next calpha
                calpha_2 = alignment_keys[i+1]
                rt_2 = alignments[calpha_2]

                assert calpha_1 == ref_c_alpha_labels[i]
                assert calpha_2 == ref_c_alpha_labels[i+1]

                # Calculate the difference in the angles of the alignment matrices
                theta_1 = scitbx.math.math.acos((rt_1.r.trace()-1)/2.0)
                theta_2 = scitbx.math.math.acos((rt_2.r.trace()-1)/2.0)
                # XXX Should we calculate the absolute of the difference?
                theta_rad = theta_2-theta_1
                theta_deg = theta_rad * 180.0/scitbx.math.math.pi
                # Calculate the difference in the translation
                t_shift = (rt_2.t-rt_1.t).norm_sq()**0.5

#                # Calculate the angles from the multiplication of one by the inverse of the other
#                rt_1_2 = rt_1 * rt_2.inverse()
#                # Calculate the angle of the rotation matrix
#                theta_rad = scitbx.math.math.acos((rt_1_2.r.trace()-1)/2.0)
#                theta_deg = theta_rad * 180.0/scitbx.math.math.pi
#                # Calculate the length of the shift
#                t_shift =  rt_1_2.t.norm_sq()**0.5

                # Append to the array
                output_diffs[d_num, i, :] = theta_deg, t_shift

        # Directory to write the output to
        var_out_dir = self.file_manager.get_dir('analyses')
        # Write out to file
        numpy.savetxt(  fname = os.path.join(var_out_dir, 'calpha_rt_r_variation.csv'), X=output_diffs[:,:,0], delimiter=',', newline='\n' )
        numpy.savetxt(  fname = os.path.join(var_out_dir, 'calpha_rt_t_variation.csv'), X=output_diffs[:,:,1], delimiter=',', newline='\n' )

        # Write out graphs
        if self.settings.plot_graphs:

            # Create labels
            labels = ['']*num_pairs
            for i in range(0, num_pairs, 5)+[num_pairs-1]:
                labels[i] = i+1
            # Clear the last n before the last one
            n = 4
            labels[-1-n:-1] = ['']*n

            # BOX PLOT OF ROTATION AND TRANSLATION SHIFTS
            fig = pyplot.figure()
            pyplot.title('Rotation-translation alignment matrix variation between adjacent C-alpha')
            # ADJACENT ANGLE VARIATION
            pyplot.subplot(2, 1, 1)
            pyplot.boxplot(x=output_diffs[:,:,0], notch=True, sym='.', widths=0.5, whis=[5,95], whiskerprops={'ls':'-'}, flierprops={'ms':1}, labels=labels) # whis='range'
            pyplot.xlabel('C-alpha index')
            pyplot.ylabel('Angle Difference\n(degrees)')
            # ADJACENT SHIFT VARIATION
            pyplot.subplot(2, 1, 2)
            pyplot.boxplot(x=output_diffs[:,:,1], notch=True, sym='.', widths=0.5, whis=[5,95], whiskerprops={'ls':'-'}, flierprops={'ms':1}, labels=labels) # whis='range'
            pyplot.xlabel('C-alpha Index')
            pyplot.ylabel('Translation Difference\n(angstroms)')
            # Apply tight layout to prevent overlaps
            fig.set_tight_layout(True)
            # Save both
            pyplot.savefig(os.path.join(var_out_dir, 'calpha_rt_variation.png'), format='png')
            pyplot.close(fig)

    #########################################################################################################
    #                                                                                                       #
    #                                       Dataset selection functions                                     #
    #                                                                                                       #
    #########################################################################################################

    def select_resolution_limits(self):
        """Generate a set of resolution limits for the data analysis"""

        self.log.heading('Selecting resolution limits for analysis')

        # ================================================>
        # Use resolution of already-loaded maps
        # ================================================>
        if (self.args.flags.recalculate_statistical_maps == "No"):
            res_limits = self.stat_maps.get_resolutions()
            assert len(self.stat_maps.get_resolutions()) > 0
            self.log('----------------------------------->>>', True)
            self.log('Using previously loaded statistical maps')
            self.log('Resolutions of reloaded maps: {!s}'.format(', '.join(map(str,res_limits))), True)
            return res_limits

        # ================================================>
        # Select resolution range based on input parameters
        # ================================================>
        # Update the resolution limits using the resolution limits from the datasets supplied
        if self.params.analysis.dynamic_res_limits:
            self.log('Updating resolution limits based on the resolution of the loaded datasets')
            small_limit = max(self.params.analysis.high_res_upper_limit,
                              self.datasets.reference().data.summary.high_res,
                              min(self.tables.dataset_info['high_resolution']))
            large_limit = min(self.params.analysis.high_res_lower_limit,
                              max(self.tables.dataset_info['high_resolution']))
            # Round the limits up and down to create sensible limits
            small_limit = round(small_limit - 0.005, 2)    # i.e. 1.344 -> 1.34
            large_limit = round(large_limit + 0.005, 2)    # i.e. 3.423 -> 3.43
            self.log('Input limits: ({}-{})'.format(self.params.analysis.high_res_upper_limit,self.params.analysis.high_res_lower_limit))
            self.log('Updated limits:  ({}-{})'.format(small_limit, large_limit))
        else:
            self.log('Using the resolution limits defined in the input parameters')
            small_limit = self.params.analysis.high_res_upper_limit
            large_limit = self.params.analysis.high_res_lower_limit
            self.log('Resolution Limits: ({}-{})'.format(small_limit, large_limit))

        # ============================================================================>
        # Determine the range of resolutions already covered
        # ============================================================================>
        if (self.args.flags.recalculate_statistical_maps == "Extend"):
            # Extract the resolution limits from previous runs
            curr_res_limits = self.stat_maps.get_resolutions()
            curr_small_limit = min(curr_res_limits)
            curr_large_limit = max(curr_res_limits)
        else:
            curr_res_limits = []
            curr_small_limit = small_limit
            curr_large_limit = small_limit
        # ============================================================================>
        # Combine to create new set of resolution limits
        # ============================================================================>
        if not self.params.analysis.high_res_increment:
            # No variable cutoff - select all
            res_limits = curr_res_limits + [large_limit]  # i.e. [2]
        else:
            small_res_limits = [round(x, 4) for x in numpy.arange(small_limit, curr_small_limit, self.params.analysis.high_res_increment).tolist()]
            large_res_limits = [round(x, 4) for x in numpy.arange(curr_large_limit, large_limit, self.params.analysis.high_res_increment).tolist()]
            res_limits = small_res_limits + curr_res_limits + large_res_limits + [large_limit]

        return sorted(set(res_limits))

    def select_datasets_for_density_characterisation(self, high_res_cutoff):
        """Select all datasets with resolution better than high_res_cutoff"""

        building_mask_name = 'characterisation @ {!s}A'.format(high_res_cutoff)

        # Create empty mask
        self.datasets.all_masks().add_mask(name=building_mask_name, values=False)
        # Counter for the number of datasets to select
        total_build = 0
        # Select from the datasets that haven't been rejected
        for dataset in self.datasets.mask(mask_name='rejected - total', invert=True):
            # Check the resolution of the dataset
            if self.tables.dataset_info.get_value(index=dataset.tag, col='high_resolution') > high_res_cutoff:
                continue
            # Check to see if this has been excluded from building
            elif self.datasets.all_masks().get_value(name='exclude_from_characterisation', id=dataset.tag) == True:
                self.log('Rejecting Dataset {!s}: Excluded from building'.format(dataset.tag))
                continue
            else:
                self.datasets.all_masks().set_value(name=building_mask_name, id=dataset.tag, value=True)
                # Check to see if the number of datasets to use in building has been reached
                total_build += 1
                if total_build >= self.params.analysis.max_build_datasets:
                    self.log('Maximum number of datasets for building reached: {!s}={!s}'.format(total_build, self.params.analysis.max_build_datasets))
                    break
        return building_mask_name, self.datasets.all_masks().get_mask(name=building_mask_name)

    def select_datasets_for_analysis(self, high_res_large_cutoff, high_res_small_cutoff):
        """Select all datasets with resolution between high and low limits"""

        assert high_res_large_cutoff > high_res_small_cutoff, '{!s} must be larger than {!s}'.format(high_res_large_cutoff, high_res_small_cutoff)

        analysis_mask_name = 'analysis @ {!s}A'.format(high_res_large_cutoff)

        if self.args.flags.reprocess_selected_datasets: datasets_for_reprocessing = self.args.flags.reprocess_selected_datasets.split(',')
        else:                                           datasets_for_reprocessing = []

        # Create empty mask
        self.datasets.all_masks().add_mask(name=analysis_mask_name, values=False)
        # Select from the datasets that haven't been rejected
        for dataset in self.datasets.mask(mask_name='rejected - total', invert=True):
            # Check the resolution of the dataset (is not too low)
            if self.tables.dataset_info.get_value(index=dataset.tag, col='high_resolution') > high_res_large_cutoff:
                continue
            # Check the resolution of the dataset (is not too high)
            elif self.tables.dataset_info.get_value(index=dataset.tag, col='high_resolution') <= high_res_small_cutoff:
                continue
            # Check to see if this has been excluded from building
            elif self.datasets.all_masks().get_value(name='exclude_from_zmap_analysis', id=dataset.tag) == True:
                self.log('Not selecting dataset {!s}: excluded from analysis'.format(dataset.tag))
                continue
            elif self.datasets.all_masks().get_value(name='old datasets', id=dataset.tag) and (not self.args.flags.reprocess_existing_datasets) and (dataset.tag not in datasets_for_reprocessing):
                self.log('Not selecting dataset {!s}: already processed (old dataset)'.format(dataset.tag))
                continue
            else:
                self.datasets.all_masks().set_value(name=analysis_mask_name, id=dataset.tag, value=True)
        return analysis_mask_name, self.datasets.all_masks().get_mask(name=analysis_mask_name)

    #########################################################################################################
    #                                                                                                       #
    #                                               Event functions                                         #
    #                                                                                                       #
    #########################################################################################################

    def collate_event_counts(self):
        """Collate events from all of the datasets"""

        self.log('----------------------------------->>>', False)
        self.log('Collating Clusters', False)

        # List of points to be returned
        all_dataset_events = dict([(d.tag, d.events) for d in self.datasets.all()])

        # Print Cluster Summaries
        event_num = [(k, len(all_dataset_events[k])) for k in sorted(all_dataset_events.keys()) if all_dataset_events[k]]
        event_total = sum([a[1] for a in event_num])

        return event_total, event_num, all_dataset_events

    def cluster_events_and_update(self, events=[], update_tables=True, update_output=True):
        """Cluster events to sites and add information to the pandda tables"""

        if not events:
            print('No Events Found')
            return None

        self.log('----------------------------------->>>', True)
        self.log('Clustering identified events: {} Event(s)'.format(len(events)))
        site_list = cluster_events(events=events, cutoff=15.0/self.grid.grid_spacing(), linkage='average')
        site_list.sort(key=lambda s: (s.info.num_events, max([e.cluster.max for e in s.children])), reverse=True).renumber()
        # Add meta to the site list TODO implement this function -- blank at the moment TODO
        [s.find_protein_context(hierarchy=self.datasets.reference().model.hierarchy) for s in site_list.children]
        # Update the pandda tables?
        if update_tables:
            self.update_site_table(site_list=site_list, clear_table=True)
            self.update_event_table_site_info(events=events)
        # Generate output images and graphs?
        if update_output:
            # Plot output graph of site list
            self.log('Deleting old images: ')
            delete_with_glob(glob_str=self.file_manager.get_file('analyse_site_graph_mult').format('*'))
            bar.multiple_bar_plot_over_several_images(
                                    f_template = self.file_manager.get_file('analyse_site_graph_mult'),
                                    plot_vals  = [sorted([e.cluster.max for e in s.children],reverse=True) for s in site_list.children]   )
            # Create pictures of the sites on the protein
            self.make_pymol_site_image_and_scripts(site_list=site_list, make_images=True)

        return site_list

#    def image_blob(self, script, image, dataset, point, point_no, towards=[10,10,10]):
#        """Take pictures of the maps with ccp4mg"""
#
#        from giant.graphics import calculate_view_quaternion, multiply_quaternions
#
#        # Get the template to be filled in
#        template = PANDDA_HTML_ENV.get_template('ccp4mg-pic.py')
#
#        orientation = calculate_view_quaternion(towards, point)
#        rotate_1 = multiply_quaternions(orientation, (0.0, 0.5**0.5, 0.0, 0.5**0.5))
#        rotate_2 = multiply_quaternions(orientation, (0.5**0.5, 0.0, 0.0, 0.5**0.5))
#
#        for view_no, view in enumerate([orientation, rotate_1, rotate_2]):
#
#            view_script = script.format(point_no, view_no)
#            view_image  = image.format(point_no, view_no)
#
#            ccp4mg_script = template.render({
#                                                'view'  :{
#                                                                'camera_centre' : [-1*c for c in point],
#                                                                'orientation'   : list(view)
#                                                            },
#                                                'mol'   :{
#                                                                'path'  : dataset.file_manager.get_file('aligned_structure'),
#                                                                'name'  : 'aligned_structure'
#                                                            },
#                                         #       'map'   :{
#                                         #                       'path'    : dataset.file_manager.get_file('sampled_map'),
#                                         #                       'name'    : 'sampled_map',
#                                         #                       'contour' : [1]
#                                         #                   },
#                                                'diff_map' :{
#                                                                'path'    : dataset.file_manager.get_file('z_map_corrected_normalised'),
#                                                                'name'    : 'diff_map',
#                                                            #    'neg-contour' : -3,
#                                                                'pos-contour' : [2,3,4,5]
#                                                            }
#                                            })
#
#            # Write out the ccp4mg script to the dataset's scripts folder
#            with open(view_script, 'w') as fh:
#                fh.write(ccp4mg_script)
#
#            # Make the images
#            c = CommandManager('ccp4mg')
#            c.SetArguments(['-norestore','-picture', view_script, '-R', view_image, '-RO', """'{"size":"1500x1500"}'""", '-quit'])
#            c.Run()
#
#            if not os.path.exists(view_image):
#                print('FAILED TO MAKE IMAGES')
#                print(c.err)

    #########################################################################################################
    #                                                                                                       #
    #                                         Summary/output functions                                      #
    #                                                                                                       #
    #########################################################################################################

    def write_map_analyser_maps(self, map_analyser):
        """Write statistical maps for a map_analyser object"""

        # Resolution for file naming
        map_res = map_analyser.meta.resolution
        # Extract stat maps for this resolution
        stat_maps = map_analyser.statistical_maps
        # Set symmetry mask values to zero
        mask = flex.size_t(self.grid.symmetry_mask().inner_mask_indices())

        self.log('----------------------------------->>>')
        self.log('=> Writing characterised statistical maps @ {!s}A'.format(map_res))

        self.grid.write_array_as_map(array  = stat_maps.mean_map.get_map_data(sparse=False).as_1d().set_selected(mask, 0.0),
                                     f_name = self.file_manager.get_file('mean_map').format(map_res))
        self.grid.write_array_as_map(array  = stat_maps.stds_map.get_map_data(sparse=False).as_1d().set_selected(mask, 0.0),
                                     f_name = self.file_manager.get_file('stds_map').format(map_res))
        self.grid.write_array_as_map(array  = stat_maps.sadj_map.get_map_data(sparse=False).as_1d().set_selected(mask, 0.0),
                                     f_name = self.file_manager.get_file('sadj_map').format(map_res))
        self.grid.write_array_as_map(array  = stat_maps.skew_map.get_map_data(sparse=False).as_1d().set_selected(mask, 0.0),
                                     f_name = self.file_manager.get_file('skew_map').format(map_res))
        self.grid.write_array_as_map(array  = stat_maps.kurt_map.get_map_data(sparse=False).as_1d().set_selected(mask, 0.0),
                                     f_name = self.file_manager.get_file('kurt_map').format(map_res))
        self.grid.write_array_as_map(array  = stat_maps.bimo_map.get_map_data(sparse=False).as_1d().set_selected(mask, 0.0),
                                     f_name = self.file_manager.get_file('bimo_map').format(map_res))

    def write_output_csvs(self):
        """Write CSV file of dataset variables"""

        self.log.subheading('Writing output CSVs', False, True)

        # Write the dataset information to csv file
        self.log.bar()
        self.log('Writing Dataset + Dataset Map Summary CSVs')
        self.log('\t'+self.file_manager.get_file('dataset_info'))
        self.tables.dataset_info.dropna(axis=1, how='all').to_csv(path_or_buf=self.file_manager.get_file('dataset_info'), index_label='dtag')
        self.log('\t'+self.file_manager.get_file('dataset_map_info'))
        self.tables.dataset_map_info.dropna(axis=1, how='all').to_csv(path_or_buf=self.file_manager.get_file('dataset_map_info'), index_label='dtag')
        self.log('\t'+self.file_manager.get_file('dataset_masks'))
        self.datasets.all_masks().table.to_csv(path_or_buf=self.file_manager.get_file('dataset_masks'), index_label='id')

        # Join the tables on the index of the main table
        self.log.bar()
        self.log('Writing COMBINED Dataset Summary CSV')
        comb_tab = (
            self.tables.dataset_info
                .join(self.tables.dataset_map_info, how='outer')
                .join(self.datasets.all_masks().table[PanddaMaskNames.write_mask_names], how='outer')
        )
        self.log('\t'+self.file_manager.get_file('dataset_combined_info'))
        comb_tab.dropna(axis=1, how='all').to_csv(path_or_buf=self.file_manager.get_file('dataset_combined_info'), index_label='dtag')

        # Write the event data only once events have been recorded
        if len(self.tables.event_info.index):
            self.log.bar()
            self.log('Writing Event+Site Summary CSVs')
            # Sort the event data by z-peak and write out
            sort_eve = self.tables.event_info.sort_values(by=['site_idx',self.args.results.events.order_by], ascending=[1,0])
            sort_eve = sort_eve.join(comb_tab, how='right')
            self.log('\t'+self.file_manager.get_file('event_info'))
            sort_eve.to_csv(path_or_buf=self.file_manager.get_file('event_info'))
            # Sort the sites by number of events and write out
            sort_sit = self.tables.site_info.sort_values(by=[self.args.results.sites.order_by],ascending=[0])
            self.log('\t'+self.file_manager.get_file('site_info'))
            sort_sit.to_csv( path_or_buf=self.file_manager.get_file('site_info'))

    def update_site_table(self, site_list, clear_table=True):
        """Add site entries to the site table"""

        # Clear an existing table
        if clear_table:
            self.tables.site_info = pandas.DataFrame(data    = None,
                                                     index   = self.tables.site_info.index.reindex([])[0],
                                                     columns = self.tables.site_info.columns)
        # Go through and update the site information
        for site in site_list.children:
            self.tables.site_info.loc[site.id,:] = None
            centroid_cart = tuple(flex.double(site.info.centroid)*self.grid.grid_spacing())
            self.tables.site_info.set_value(site.id, 'centroid', centroid_cart)
            self.tables.site_info.set_value(site.id, 'native_centroid', tuple(self.datasets.reference().model.alignment.ref2nat(coordinates=[centroid_cart])[0]))

    def make_pymol_site_image_and_scripts(self, site_list, make_images=True):
        """Generate pymol script to mark the location of identified sites"""

        pymol_str =  '# Mark the identified sites on the protein\n'
        pymol_str += 'from pymol import cmd\n'
        pymol_str += 'from pymol.cgo import *\n'
        pymol_str += 'cmd.load("{}", "reference")\n'.format(os.path.relpath(self.file_manager.get_file('reference_on_origin'), start=self.file_manager.get_dir('output_summaries')))
        pymol_str += 'cmd.show_as("cartoon", "reference")\n'
        pymol_str += 'cmd.color("cyan", "reference")\n'
        # Add sphere at each of the sites
        for site in site_list.children:
            # Only print the site if it has more than one event
            if len(site.children) > 1:
                lab = 'site_{}'.format(site.id)
                com = tuple(flex.double(site.info.centroid)*self.grid.grid_spacing())
                pymol_str += 'cmd.pseudoatom("{}", pos={}, vdw=2.5)\n'.format(lab, com)
                pymol_str += 'cmd.show("sphere", "{}")\n'.format(lab)
                pymol_str += 'cmd.label("{}", "{}")\n'.format(lab, site.id)
                pymol_str += 'cmd.color("deepteal", "{}")\n'.format(lab)
                pymol_str += 'cmd.set("label_color", "white", "{}")\n'.format(lab)
            # Label events as smaller spheres
            for event in site.children:
                lab = 'event'
                com = tuple(flex.double(event.cluster.centroid)*self.grid.grid_spacing())
                pymol_str += 'cmd.pseudoatom("{}", pos={}, vdw=0.5)\n'.format(lab, com)
                pymol_str += 'cmd.show("sphere", "{}")\n'.format(lab)
                pymol_str += 'cmd.color("blue", "{}")\n'.format(lab)
        # Set label things...
        pymol_str += 'cmd.set("label_size", 25)\n'
        pymol_str += 'cmd.set("label_position", (0,0,4))\n'
        pymol_str += 'cmd.bg_color(color="white")\n'
        # Write as python script
        with open(self.file_manager.get_file(file_tag='pymol_sites_py'), 'w') as fh:
            fh.write(pymol_str)

        # Run Pymol to generate images and output to pngs
        if make_images:
            pymol_str =  '# Load the protein representation and output images of sites\n'
            pymol_str += 'run {}\n'.format(os.path.relpath(self.file_manager.get_file(file_tag='pymol_sites_py'), start=self.file_manager.get_dir('output_summaries')))
            pymol_str += 'set ray_opaque_background, off\n'
            pymol_str += 'set specular, off\n'
            pymol_str += 'orient\n'
            pymol_str += 'png {}, width=1200, dpi=300, ray=1\n'.format(os.path.relpath(self.file_manager.get_file(file_tag='pymol_sites_png_1'), start=self.file_manager.get_dir('output_summaries')))
            pymol_str += 'rotate y, 180\n'
            pymol_str += 'png {}, width=1200, dpi=300, ray=1\n'.format(os.path.relpath(self.file_manager.get_file(file_tag='pymol_sites_png_2'), start=self.file_manager.get_dir('output_summaries')))
            pymol_str += 'quit'

            with open(self.file_manager.get_file(file_tag='pymol_sites_pml'), 'w') as fh:
                fh.write(pymol_str)

            # Change into directory as script runs off of relative paths
            os.chdir(self.file_manager.get_dir('output_summaries'))
            c = CommandManager('pymol')
            c.add_command_line_arguments(['-k', '-q', '-c', self.file_manager.get_file(file_tag='pymol_sites_pml')])
            try:    c.run()
            except: print("Failed to start pymol - maybe it's not available?")
            # Change back to top directory
            os.chdir(self.out_dir)

        #os.remove(self.file_manager.get_file(file_tag='pymol_sites_pml'))
        #os.remove(self.file_manager.get_file(file_tag='pymol_sites_py'))

    def add_event_to_event_table(self, dataset, event):
        """Add event entries to the event table"""

        # Check event has not been added previously
        assert event.id not in self.tables.event_info.index.values.tolist(), 'Event Already Added!: {!s}'.format(event.id)
        # Add values to a new row in the table
        self.tables.event_info.loc[event.id,:] = None
        # Default to site_idx of 0 if no site given
        if event.parent:    site_idx = event.parent.id
        else:               site_idx = 0
        self.tables.event_info.set_value(event.id, 'site_idx', site_idx)
        # Event and cluster information
        self.tables.event_info.set_value(event.id, '1-BDC',  round(1.0-event.info.estimated_bdc,2))
        self.tables.event_info.set_value(event.id, 'z_peak', round(event.cluster.max,2))
        self.tables.event_info.set_value(event.id, 'z_mean', round(event.cluster.mean,2))
        self.tables.event_info.set_value(event.id, 'cluster_size', event.cluster.size)
        self.tables.event_info.set_value(event.id, ['refx','refy','refz'], list(self.grid.grid2cart([event.cluster.peak],origin=False)[0]))
        self.tables.event_info.set_value(event.id, ['x','y','z'], list(dataset.model.alignment.ref2nat(coordinates=self.grid.grid2cart([event.cluster.peak],origin=False))[0]))

    def update_event_table_site_info(self, events):
        """Update the event table for pre-existing events"""
        for e in events:
            assert e.id, 'NO ID GIVEN: {!s}'.format(e.id)
            assert e.parent, 'EVENT HAS NO PARENT: {!s}'.format(e.parent)
            self.tables.event_info.set_value(e.id, 'site_idx', e.parent.id)


class PanddaMapAnalyser(object):


    def __init__(self, dataset_maps, meta=None, statistical_maps=None, parent=None, log=None):
        """Class to hold dataset maps, statistical maps and meta data for a set of related maps. Also holds functions for analysing the maps."""

        # Validate the dataset maps
        assert (dataset_maps is None) or isinstance(dataset_maps, MapHolderList), 'dataset_maps must be stored in a MapHolderList. Type given: {!s}'.format(type(dataset_maps))
        # Validate the meta
        if meta:
            assert isinstance(meta, Meta), 'meta must be of type Meta. Type given: {!s}'.format(type(meta))
        else:
            meta = Meta()
        # Validate the statistical maps
        if statistical_maps:
            assert isinstance(statistical_maps, PanddaStatMapList), 'statistical_maps must be of type MapList. Type given: {!s}'.format(type(statistical_maps))
        else:
            statistical_maps = PanddaStatMapList()
        # Validate the parent object (main pandda object)
        if parent:
            assert isinstance(parent, PanddaMultiDatasetAnalyser), 'parent must be of type PanddaMultiDatasetAnalyser. Type given: {!s}'.format(type(parent))
        # Create log as appropriate
        if log:
            pass
        elif parent:
            log = parent.log
        else:
            log = Log(verbose=False)

        self.dataset_maps = dataset_maps
        self.meta = meta
        self.statistical_maps = statistical_maps
        self.parent = parent
        self.log = log

        if self.dataset_maps is not None: self._validate()

    def _validate(self):
        """Check that all of the added maps are the same size etc..."""

        m_ref = None
        for m in self.dataset_maps.all():
            if m_ref is None: m_ref = m

            assert m.data.all()  == m_ref.data.all(), 'Not all maps are the same shape'
            assert m.data.size() == m_ref.data.size(), 'Not all maps are the same size'

            if not hasattr(m.meta, 'map_uncertainty'):
                m.meta.map_uncertainty = None

        self.meta.map_data_size  = m_ref.data.size()
        self.meta.map_data_shape = m_ref.data.all()

    def calculate_mean_map(self, mask_name=None):
        """Calculate the mean map from all of the different observations"""

        # Extract the maps to be used for averaging
        dataset_maps = self.dataset_maps.mask(mask_name=mask_name)

        if len(dataset_maps) == 1:
            self.log.bar(True, False)
            self.log('One dataset map provided for mean map calculation -- using this map.', True)

            # Extract the map from the list
            m = dataset_maps[0]
            # Mean and median are simply the map value -- copy directly to the statistical maps
            mean_map_vals = medn_map_vals = numpy.array(m.data)

        else:
            self.log.bar(True, False)
            self.log('Calculating mean map from {} dataset maps'.format(len(dataset_maps)), True)

            # Chunk the points into groups - Compromise between cpu time and memory usage - ~200 dataset -> chunksize of 5000
            chunk_size = 500*iceil(1000.0/len(dataset_maps))
            chunk_idxs = [i for i in range(0, self.meta.map_data_size, chunk_size)]
            num_chunks = len(chunk_idxs)

            self.log('Iterating through {!s} points in {!s} chunks'.format(self.meta.map_data_size, num_chunks), True)
            t1 = time.time()

            mean_map_vals = numpy.zeros(self.meta.map_data_size)
            medn_map_vals = numpy.zeros(self.meta.map_data_size)

            for i_chunk, chunk_start in enumerate(chunk_idxs):
                status_bar_2(n=i_chunk, n_max=num_chunks)

                tmp_map_vals = numpy.array([m.data[chunk_start:chunk_start+chunk_size] for m in dataset_maps])

                # Check that the output values are the expected dimensions
                if i_chunk+1 < num_chunks:
                    assert len(tmp_map_vals) == len(dataset_maps)
                    assert len(tmp_map_vals.T) == chunk_size

                tmp_map_means = numpy.mean(tmp_map_vals, axis=0)
                mean_map_vals[chunk_start:chunk_start+chunk_size] = tmp_map_means
                tmp_map_medns = numpy.median(tmp_map_vals, axis=0)
                medn_map_vals[chunk_start:chunk_start+chunk_size] = tmp_map_medns

            status_bar_2(n=num_chunks, n_max=num_chunks)

            t2 = time.time()
            self.log('> Calculation of mean map > Time Taken: {!s} seconds'.format(int(t2-t1)))

        # Store the mean and median map values in the statistical maps object
        self.statistical_maps.mean_map = m.new_from_template(map_data=flex.double(mean_map_vals.flatten()), sparse=m.is_sparse())
        self.statistical_maps.medn_map = m.new_from_template(map_data=flex.double(medn_map_vals.flatten()), sparse=m.is_sparse())

        return self.statistical_maps.mean_map

    def calculate_map_uncertainties(self, masked_idxs=None, mask_name=None, q_cut=1.5, cpus=1):
        """Calculate the uncertainty in each of the different maps"""

        if masked_idxs is None:
            masked_idxs = flex.size_t(range(0, self.meta.map_data_size))
        else:
            assert max(masked_idxs) < self.meta.map_data_size, 'masked_idxs out of range of map'
            masked_idxs = flex.size_t(masked_idxs)

        # Extract masked map values from the mean map... and sort them
        mean_vals = self.statistical_maps.mean_map.data.select(masked_idxs)

        self.log.bar()
        self.log('Selecting datasets for uncertainty calculation')
        self.log.bar()

        arg_list = []

        for i_m, m in enumerate(self.dataset_maps.mask(mask_name=mask_name)):

            if m.meta.map_uncertainty is not None:
                print('Uncertainty already calculated for dataset {!s} ({!s}/{!s})'.format(m.meta.tag, i_m+1, self.dataset_maps.size(mask_name=mask_name)))
                arg_list.append(None)
                continue

            # Get the file manager for the dataset to allow outputting of files
            if self.parent:
                file_manager = self.parent.datasets.get(tag=m.meta.tag).file_manager
            else:
                fine_manager = None

            u = UncertaintyCalculator(query_values=m.data.select(masked_idxs), ref_values=mean_vals, file_manager=file_manager)
            arg_list.append(u)

        t1 = time.time()
        num_to_process = len(arg_list)-arg_list.count(None)
        self.log('Calculating uncertainties of {!s} maps (using {!s} cores)'.format(num_to_process, cpus))
        self.log('----------------------------------->>>')
        print('1'+''.join(['{:<5}'.format(i) for i in range(0,num_to_process+5,5)])[2:])
        print(' '*num_to_process+'|\r', end=''); sys.stdout.flush()
        map_uncertainties = easy_mp.pool_map(func=wrapper_run, args=arg_list, processes=cpus, chunksize=1)
        print('|')
        t2 = time.time()
        self.log('> Calculation of map uncertainties > Time Taken: {!s} seconds'.format(int(t2-t1)))
        self.log('----------------------------------->>>')

        for i_m, m in enumerate(self.dataset_maps.mask(mask_name=mask_name)):
            map_unc = map_uncertainties[i_m]
            if m.meta.map_uncertainty is not None:
                assert map_unc is None
            else:
                assert map_unc is not None
                m.meta.map_uncertainty = map_unc
                self.log('Dataset {}: Uncertainty {:.4f}'.format(m.meta.tag, m.meta.map_uncertainty))

        return [m.meta.map_uncertainty for m in self.dataset_maps.mask(mask_name=mask_name)]

    def calculate_statistical_maps(self, mask_name=None, ignore_warnings=True, cpus=1):
        """Take the sampled maps and calculate statistics for each grid point across the datasets"""

        # Extract the maps to be used for averaging
        dataset_maps = self.dataset_maps.mask(mask_name=mask_name)

        if len(dataset_maps) == 1:
            self.log.bar(True, False)
            self.log('One dataset map provided for statistical map calculation -- setting all statistical values to zero.', True)

            # Output array of the 5 statistics for each map point
            point_statistics = numpy.zeros((self.meta.map_data_size, 5))

        else:
            self.log.bar()
            self.log('Calculating statistics of {} grid points (using {!s} cores)'.format(self.meta.map_data_size, cpus))

            # Create statistics objects for each grid point
            if ignore_warnings:
                self.log('Suppressing Numpy Iteration Warnings... (Nothing to worry about)')
                warnings.simplefilter('ignore', category=RuntimeWarning)

            # Extract the map uncertainties
            uncertainties = [m.meta.map_uncertainty for m in dataset_maps]
            assert uncertainties.count(None) == 0, 'some maps have not got associated uncertainties'

            # Chunk the points into groups - Compromise between cpu time and memory usage - 1000 per cpu at 50 datasets
            chunk_size = iceil(1000.0*cpus*50.0/len(dataset_maps))
            chunk_idxs = [i for i in range(0, self.meta.map_data_size, chunk_size)]
            num_chunks = len(chunk_idxs)

            # Second level of iteration - split the first chunk level between the cpus
            chunk_size_2 = iceil(1.0*chunk_size/cpus)
            chunk_idxs_2 = [i for i in range(0, chunk_size, chunk_size_2)]
            num_chunks_2 = len(chunk_idxs_2)

            self.log('Iterating through {!s} points in blocks of {!s} ({!s} chunks)'.format(self.meta.map_data_size, chunk_size, num_chunks))

            t1 = time.time()

            # Output array of the 5 statistics for each map point
            point_statistics = numpy.zeros((self.meta.map_data_size, 5))

            tot = 0
            for i_chunk, chunk_start in enumerate(chunk_idxs):
                status_bar_2(n=i_chunk, n_max=num_chunks)

                # Argument list for multiprocessing
                arg_list = []

                # Loop through the secondary chunks and send for multi-core processing
                for i_chunk_2, chunk_start_2 in enumerate(chunk_idxs_2):

    #                print('\n--------------------->')
    #                print('Chunk {}-{}'.format(i_chunk, i_chunk_2))
    #                print('Getting {} - {}'.format(chunk_start+chunk_start_2, chunk_start+chunk_start_2+chunk_size_2))

                    # Lower limit - always the beginning of the chunk
                    l1 = chunk_start+chunk_start_2
                    # Upper limit - full chunk size, limited by the larger chunk size, or by map size
                    l2 = min(chunk_start+chunk_start_2+chunk_size_2, chunk_start+chunk_size, self.meta.map_data_size)

                    if l1 >= l2:
                        continue

                    # Extract map values from the maps
                    map_vals = [m.data[l1:l2] for m in dataset_maps]
                    # Want to iterate over grid points not datasets
                    map_vals = numpy.transpose(map_vals)
                    assert map_vals.shape[1] == len(dataset_maps)

                    # Create DensityStatistics object for analysis of the density variation
                    arg_list.append(DensityStatistics(observations_array=map_vals, uncertainties=uncertainties))

                if not arg_list: continue

                # Calculate the statistis of the grid points
                tmp_point_statistics = easy_mp.pool_map(func=wrapper_run, args=arg_list, processes=cpus)

                # Put values into the output array
                offset = 0
                for point_vals in tmp_point_statistics:
                    assert point_vals.shape[1] == 5
                    l1 = chunk_start+offset
                    l2 = l1+point_vals.shape[0]
                    if not (point_statistics[l1:l2,:] == 0.0).all():
                        print('Overwriting data?!')
                        print(point_statistics[l1-10:l2+10,:])
                        assert point_statistics[l1:l2,:].shape == point_vals.shape, '{} != {}'.format(point_statistics[l1:l2,:].shape, point_vals.shape)
                    point_statistics[l1:l2,:] = point_vals
                    offset += point_vals.shape[0]
                tot += offset

            status_bar_2(n=num_chunks, n_max=num_chunks)

            # Check that we've calculated the right number of things
            assert tot == self.meta.map_data_size, 'tot {}, map size {}'.format(tot, self.meta.map_data_size)

            t2 = time.time()
            self.log('> Calculation of map variation statistics > Time Taken: {!s} seconds'.format(int(t2-t1)))

        # Use the mean map as the template for the other maps
        mean_map = self.statistical_maps.mean_map
        # Create the other statistical maps
        self.statistical_maps.stds_map = mean_map.new_from_template(map_data=flex.double(point_statistics[:,0].tolist()), sparse=mean_map.is_sparse())
        self.statistical_maps.sadj_map = mean_map.new_from_template(map_data=flex.double(point_statistics[:,1].tolist()), sparse=mean_map.is_sparse())
        self.statistical_maps.skew_map = mean_map.new_from_template(map_data=flex.double(point_statistics[:,2].tolist()), sparse=mean_map.is_sparse())
        self.statistical_maps.kurt_map = mean_map.new_from_template(map_data=flex.double(point_statistics[:,3].tolist()), sparse=mean_map.is_sparse())
        self.statistical_maps.bimo_map = mean_map.new_from_template(map_data=flex.double(point_statistics[:,4].tolist()), sparse=mean_map.is_sparse())

        return self.statistical_maps

    def calculate_z_map(self, method, map=None, tag=None, uncertainty=None):
        """Calculate the z-map relative to the mean and std map"""

        assert method in ['naive','adjusted','uncertainty','adjusted+uncertainty']
        assert [tag, map].count(None) == 1, 'Must provide tag OR map'

        if tag:
            map = self.dataset_maps.get(tag=tag)

        if uncertainty is not None:
            uncertainty = map.meta.map_uncertainty

        if 'uncertainty' in method:
            assert uncertainty is not None

        # Extract maps in the right sparseness
        is_sparse = map.is_sparse()
        # Extract mean values (for subtraction)
        mean_vals = self.statistical_maps.mean_map.get_map_data(sparse=is_sparse)

        # Extract the normalisation values (for division)
        if method == 'naive':
            norm_vals = self.statistical_maps.stds_map.get_map_data(sparse=is_sparse)
        elif method == 'adjusted':
            norm_vals = self.statistical_maps.sadj_map.get_map_data(sparse=is_sparse)
        elif method == 'uncertainty':
            norm_vals = uncertainty
        elif method == 'adjusted+uncertainty':
            norm_vals = flex.sqrt(self.statistical_maps.sadj_map.get_map_data(sparse=is_sparse)**2 + uncertainty**2)
        else:
            raise Exception('method not found: {!s}'.format(method))

        return (map - mean_vals)*(1.0/norm_vals)
