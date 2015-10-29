import os, sys, copy, re

import libtbx.phil

import numpy
import pandas

from Giant.Xray.Data import crystalSummary
from Giant.Xray.Data.cluster import crystalGroup

blank_arg_prepend = 'mtz='

master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .type = path
        .multiple = True
    mtz_regex = None
        .type = str
}
check_for{
    labels = None
        .type = str
        .multiple = True
}
output {
    out_dir = None
        .type = path
}
""")

def run(params):

    if params.output.out_dir and (not os.path.exists(params.output.out_dir)):
        os.mkdir(params.output.out_dir)
    # Dump images into the current directory if no directory is given
    if not params.output.out_dir: img_dir = './'
    else: img_dir = params.output.out_dir

    #########################################################################################################

    if params.input.mtz_regex: labels = [re.findall(params.input.mtz_regex, m)[0] for m in params.input.mtz]
    else:                      labels = ['MTZ-{:06d}'.format(i) for i in range(len(params.input.mtz))]

    #########################################################################################################

    print '============================================>'
    print 'GROUPING BY SPACE GROUP'

    cg = crystalGroup([crystalSummary.from_mtz(mtz_file=f, id=lab) for f,lab in zip(params.input.mtz, labels)])

    #########################################################################################################

    print '==============================================================================>'
    print 'SPACE GROUPS: {} - {} DATASETS'.format(','.join(cg.space_groups), len(cg.crystals))

    #########################################################################################################

    res_sorted = sorted(cg.crystals, key=lambda c: c.high_res)
    high_res_crystal = res_sorted[0]
    low_res_crystal  = res_sorted[-1]

    print '===========================>'
    print 'HIGHEST RESOLUTION: {:.2f}'.format(high_res_crystal.high_res)
    print 'DATASET + PATH: {} - {}'.format(high_res_crystal.id, high_res_crystal.mtz_file)
    print 'LOWEST RESOLUTION: {:.2f}'.format(low_res_crystal.high_res)
    print 'DATASET + PATH: {} - {}'.format(low_res_crystal.id, low_res_crystal.mtz_file)

    #########################################################################################################

    for c in cg.crystals:

        if params.check_for.labels:
            for lab in params.check_for.labels:
                if lab not in c.column_labels:
                    print 'COLUMN NOT IN FILE ({}): {}'.format(lab, c.mtz_file)




