# Adapted from Code by Sebastian Kelm

import os, sys, tempfile
import pandas

from libtbx.utils import Sorry, Failure

from bamboo.common.command import CommandManager
from bamboo.ccp4_utils import MtzSummary

class Edstats(object):


    def __init__(self, mtz_file, pdb_file, f_label=None):
        scores, command = score_with_edstats_to_dict(mtz_file=mtz_file, pdb_file=pdb_file, f_label=f_label)
        self.scores = pandas.DataFrame.from_dict(scores)
        self._command = command

        if self.scores.empty and self._command.error:
            raise Failure('EDSTATS has failed to run (error message below)\n=========>\n{!s}\n=========>'.format(self._command.error))

    def extract_residue_group_scores(self, residue_group, data_table=None, rg_label=None, column_suffix=''):
        """Extract density quality metrics for a residue group from precalculated edstats scores"""

        rg = residue_group
        # Set defaults
        if rg_label is None:    rg_label = (rg.unique_resnames()[0]+'-'+rg.parent().id+'-'+rg.resseq+rg.icode).replace(' ','')
        if data_table is None:  data_table = pandas.DataFrame(index=[rg_label], column=[])
        # Check validity
        if len(rg.unique_resnames()) != 1: raise Failure(rg_label+': More than one residue name associated with residue group -- cannot process')

        # Extract residue scores
        ed_scores = self.scores[(rg.unique_resnames()[0], rg.parent().id, rg.resseq_as_int(), rg.icode)]
        # Append scores to data_table
        data_table.set_value(index=rg_label, col='RSCC'+column_suffix, value=ed_scores['CCSa'])
        data_table.set_value(index=rg_label, col='RSR' +column_suffix, value=ed_scores['Ra']  )
        data_table.set_value(index=rg_label, col='B_AV'+column_suffix, value=ed_scores['BAa'] )
        data_table.set_value(index=rg_label, col='RSZO'+column_suffix, value=ed_scores['ZOa'] )
        data_table.set_value(index=rg_label, col='RSZD'+column_suffix, value=ed_scores['ZDa'] )

        return data_table


class EdstatsLogSummary:
    """Class to process, store and present the output from edstats.pl"""


    def __init__(self, logtext):
        self.logtext = logtext
        self.mean_biso = None
        self.rms_scores = {'atm_all_ZD+':None,'atm_all_ZD-':None,'atm_all_ZO':None,'atm_main_ZD+':None,'atm_main_ZD-':None,'atm_main_ZO':None,'atm_side_ZD+':None,'atm_side_ZD-':None,'atm_side_ZO':None,'het_all_ZD+':None,'het_all_ZD-':None,'het_all_ZO':None,'het_main_ZD+':None,'het_main_ZD-':None,'het_main_ZO':None,'het_side_ZD+':None,'het_side_ZD-':None,'het_side_ZO':None}

        # Populate fields
        self._parselogtext()

    def _parselogtext(self):
        """Read and sort through log input"""

        output = [line for line in self.logtext.split('\n') if line]
        for i, line in enumerate(output):
            if line.startswith('Overall mean Biso:'):
                self.mean_biso = float(line.replace('Overall mean Biso:',''))
            elif line.startswith('RMS Z-scores for protein residues:'):
                assert ['Main','Side','All'] == output[i+1].strip().split(), 'COLUMN HEADINGS ARE NOT WHAT IS EXPECTED! {!s}'.format(line)
                for incr in [2,3,4]:
                    score, main_val, side_val, all_val = output[i+incr].strip().split()
                    self.rms_scores['_'.join(['atm','main',score])] = float(main_val)
                    self.rms_scores['_'.join(['atm','side',score])] = float(side_val)
                    self.rms_scores['_'.join(['atm','all',score])] = float(all_val)
            elif line.startswith('RMS Z-scores for hetero residues:'):
                assert ['Main','Side','All'] == output[i+1].strip().split(), 'COLUMN HEADINGS ARE NOT WHAT IS EXPECTED! {!s}'.format(line)
                for incr in [2,3,4]:
                    score, main_val, side_val, all_val = output[i+incr].strip().split()
                    self.rms_scores['_'.join(['het','main',score])] = float(main_val)
                    self.rms_scores['_'.join(['het','side',score])] = float(side_val)
                    self.rms_scores['_'.join(['het','all',score])] = float(all_val)

########################################################################################################

def score_file_with_edstats(mtz_file, pdb_file):
    """Score pdb file against electron density"""

    # Score the complex with Edstats
    edstats = Edstats(mtz_file=mtz_file, pdb_file=pdb_file)
    # Process the std out of the program
    summary = EdstatsLogSummary(edstats._command.output)

    return edstats, summary

def score_with_edstats_to_dict(mtz_file, pdb_file, f_label=None):
    """Scores residues against density, then converts to dict"""

    # Generate the list of the EDSTATS scores for each residue
    scores, header, command = score_with_edstats_to_list(mtz_file, pdb_file, f_label=f_label)
    # Create dict for the mapping between residue_id and scores
    mapping = {}
    for label, values in scores:
        hash = {}
        for k,v in zip(header,values):
            hash[k] = v
        mapping[label] = hash

    return mapping, command

def score_with_edstats_to_list(mtz_file, pdb_file, f_label=None):
    """Scores residues against density, then returns list"""

    assert os.path.exists(mtz_file), 'MTZ file for edstats does not exist! {!s}'.format(mtz_file)
    assert os.path.exists(pdb_file), 'PDB file for edstats does not exist! {!s}'.format(mtz_file)

    # Create a file handle and path for the output
    temp_handle, temp_path = tempfile.mkstemp(suffix='.table', prefix='edstats_')

    # Collate summary of MTZ file
    m_summ = MtzSummary(mtz_file)

    # Use column labels if given
    if (f_label is not None) and (f_label not in m_summ.summary['colheadings']):
        raise Sorry('Selected f_label ({}) not found in mtz file ({}) -- mtz contains columns {}'.format(f_label, mtz_file, m_summ.summary['colheadings']))
    # else guess the labels in the mtzfile
    else:
        f_label = m_summ.label.f

    # Check for f_label
    if not f_label:
        raise Sorry('No F label selected/found in mtz file: {!s} -- mtz contains columns {}'.format(mtz_file, m_summ.summary['colheadings']))

    # Run EDSTATS on the files
    try:
        # Initialise Command Manager to run edstats
        command = CommandManager('edstats.pl')
        command.add_command_line_arguments(['-hklin',mtz_file,'-xyzin',pdb_file,'-output',temp_path,'-noerror','-flabel',f_label])
        command.set_timeout(timeout=600)
        command.run()
        # Read the output
        with os.fdopen(temp_handle) as f:
            output = f.read().strip().replace('\r\n','\n').replace('\r','\n').splitlines()
        command.file_output = output
    finally:
        os.remove(temp_path)

    # Process the output header
    if output:
        # Check header and then remove the first three columns
        header = output.pop(0).split()
        assert header[:3] == ['RT', 'CI', 'RN'], 'edstats output headers are not as expected! {!s}'.format(output)
        num_fields = len(header)
        header = header[3:]
    else:
        header = []

    # List to be returned
    outputdata = []

    # Process the rest of the data
    for line in output:
        line = line.strip()
        if not line:
            continue

        fields = line.split()
        if len(fields) != num_fields:
            raise ValueError("Error Parsing EDSTATS output: Header & Data rows have different numbers of fields")

        # Get and process the residue information - TODO CI column can include alternate conformer?! TODO
        residue, chain, resnum = fields[:3]
        try:
            resnum = int(resnum)
            inscode = ' '
        except ValueError:
            inscode = resnum[-1:]
            resnum = int(resnum[:-1])

        # Remove the processed columns
        fields = fields[3:]

        # Process the other columns (changing n/a to None and value to int)
        for i, x in enumerate(fields):
            if x == 'n/a':
                fields[i] = None
            else:
                try:
                    fields[i] = int(x)
                except ValueError:
                    try:
                        fields[i] = float(x)
                    except ValueError:
                        pass

        outputdata.append([(residue, chain, resnum, inscode),fields])

    return outputdata, header, command


