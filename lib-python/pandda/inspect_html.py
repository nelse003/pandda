import os, glob

from pandda.settings import PANDDA_TOP, PANDDA_TEXT
from pandda.html import PANDDA_HTML_ENV
from Bamboo.html import path2url
from pandda.constants import PanddaAnalyserFilenames

def write_inspect_html(top_dir, inspector):

    # Get template to be filled in
    template = PANDDA_HTML_ENV.get_template('pandda_summary.html')
    # Output file
    out_file = os.path.join(top_dir, 'results_summaries', 'pandda_inspect.html')
    # Output directory (for relative symlinks)
    out_dir  = os.path.abspath(os.path.dirname(out_file))

    all_data = inspector.log_table
    len_data = len(all_data.index)

    num_fitted = sum(all_data['Ligand Placed'])
    num_viewed = sum(all_data['Viewed'])
    num_empty  = num_viewed - num_fitted
    num_unviewed = len_data - num_viewed

    num_blobs = len_data
    num_sites = len(set(all_data['site_idx']))
    # Number of datasets with hits
    try:    num_d_hit = len(set(zip(*all_data.index[all_data['Ligand Placed'] == True])[0]))
    except: num_d_hit = 0

    # ===========================================================>
    # Construct the data object to populate the template
    output_data = {}
    output_data['header'] = 'PANDDA Inspect Summary'
    output_data['title'] = 'PANDDA Inspect Summary'
    output_data['introduction'] = 'Summary of Inspection of Datasets'
    # ===========================================================>
    # Header Images
    output_data['top_images'] = []
    output_data['top_images'].append({ 'path': './'+os.path.relpath(path=os.path.join(top_dir, 'results_summaries', PanddaAnalyserFilenames.pymol_sites_png_1), start=out_dir),
                                       'title': 'Identified Sites (Front)' })
    output_data['top_images'].append({ 'path': './'+os.path.relpath(path=os.path.join(top_dir, 'results_summaries', PanddaAnalyserFilenames.pymol_sites_png_2), start=out_dir),
                                       'title': 'Identified Sites (Back)' })
    for i_png, png in enumerate(sorted(glob.glob(os.path.join(top_dir, 'results_summaries', PanddaAnalyserFilenames.inspect_site_graph_mult).format('*')))):
        output_data['top_images'].append({ 'path': './'+os.path.relpath(path=png, start=out_dir),
                                           'title': 'Identified Site Summary ({})'.format(i_png+1) })
    # ===========================================================>
    # Summary Bar
    output_data['summary_bar'] = []
    output_data['summary_bar'].append({'colour':'info',    'text':'PANDDA-Identified Blobs: {}'.format(num_blobs)})
    output_data['summary_bar'].append({'colour':'info',    'text':'PANDDA-Identified Sites: {}'.format(num_sites)})
    output_data['summary_bar'].append({'colour':'success', 'text':'Number of Ligands Fitted: {}'.format(num_fitted)})
    output_data['summary_bar'].append({'colour':'success', 'text':'Datasets w. Fitted Ligands: {}'.format(num_d_hit)})
    # ===========================================================>
    # Progress Bars
    output_data['progress_bar'] = []
    output_data['progress_bar'].append({'title':'Fitting Progress', 'data':[]})
    output_data['progress_bar'][0]['data'].append({'text':'Fitted',           'colour':'success', 'size':100.0*num_fitted/len_data})
    output_data['progress_bar'][0]['data'].append({'text':'Unviewed',         'colour':'warning', 'size':100.0*num_unviewed/len_data})
    output_data['progress_bar'][0]['data'].append({'text':'No Ligand Fitted', 'colour':'danger',  'size':100.0*num_empty/len_data})
    # ===========================================================>
    # Tables
    output_data['table'] = {}
    output_data['table']['column_headings'] = ['Dataset','Viewed','Interesting','Lig. Placed','Event','Site','1 - BDC','Z-Peak','Map Res.','Map Unc.','Confidence','Comment','']
    output_data['table']['rows'] = []
    # Add the datasets as rows
    for i_d in range(len(all_data.index)):

        # colour choices - 'success', 'muted', 'danger'
        # icon choices   - 'ok', 'flag', 'remove'

        d_data = all_data.iloc[i_d].to_dict()
        d_tag, d_event = all_data.index[i_d]

        columns = []

        if d_data['Viewed']:        columns.append({'colour':'success', 'icon':'ok',     'message':d_data['Viewed']})
        else:                       columns.append({'colour':'danger',  'icon':'remove', 'message':d_data['Viewed']})

        if d_data['Interesting']:   columns.append({'colour':'success', 'icon':'ok',     'message':d_data['Interesting']})
        else:                       columns.append({'colour':'danger',  'icon':'remove', 'message':d_data['Interesting']})

        if d_data['Ligand Placed']: columns.append({'colour':'success', 'icon':'ok',     'message':d_data['Ligand Placed']})
        else:                       columns.append({'colour':'danger',  'icon':'remove', 'message':d_data['Ligand Placed']})

        columns.append({'message':d_event})
        columns.append({'message':d_data['site_idx']})
        columns.append({'message':round(d_data['1-BDC'],3)})
        columns.append({'message':round(d_data['z_peak'],3)})
        columns.append({'message':d_data['analysed_resolution']})
        columns.append({'message':round(d_data['map_uncertainty'],3)})

        columns.append({'message':d_data['Ligand Confidence']})
        columns.append({'message':d_data['Comment']})

        row_message = 'Hit' if d_data['Ligand Placed'] else \
                      ''
        row_colour  = 'success' if d_data['Ligand Placed'] else \
                      'danger' if not d_data['Viewed'] else \
                      'info'

        output_data['table']['rows'].append({'heading' : d_tag,
                                             'colour'  : row_colour,
                                             'message' : row_message,
                                             'columns' : columns})

    with open(out_file, 'w') as out_html:
        out_html.write(template.render(output_data))