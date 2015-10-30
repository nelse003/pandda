#!/usr/bin/env pandda.python

import os, sys, copy, re

import libtbx.phil
from Bamboo.Common.command import commandManager
from Giant.jiffies import parse_phil_args
from Giant.Xray.Data import crystalSummary

blank_arg_prepend = 'mtz='

master_phil = libtbx.phil.parse("""
input {
    mtz = None
        .type = path
        .multiple = True
    columns = F,SIGF
        .type = str
}
method {
    package = *phenix ccp4
        .type = choice
}
output {
    column = IMEAN
        .type = str
    file_prefix = intensities_
        .type = str
}
""")

def run(params):

    input_cols = params.input.columns.split(',')

    for mtz_file in params.input.mtz:

        output_mtz = os.path.join(os.path.dirname(mtz_file), params.output.file_prefix+os.path.basename(mtz_file))

        cs = crystalSummary.from_mtz(mtz_file=mtz_file)

        if params.output.column in cs.column_labels:
            print 'OUTPUT COLUMN ALREADY PRESENT IN MTZ_FILE: {}'.format(mtz_file)
            if not os.path.exists(output_mtz):
                os.symlink(os.path.basename(mtz_file), output_mtz)
            continue

        cols_correct  = True
        for col in input_cols:
            if col not in cs.column_labels:
                cols_correct = False
                break
        if cols_correct == False:
            print 'Column not present in mtz_file:'
            print cs.column_labels
            print 'SKIPPING'
            continue
       
        if params.method.package == 'phenix':
            cm = commandManager('phenix.reflection_file_converter')
            cm.add_command_line_arguments([ mtz_file,
                                            '--mtz-root-label='+params.output.column,
                                            '--label='+params.input.columns,
                                            '--write-mtz-intensities',
                                            '--mtz='+output_mtz     ])
        elif params.method.package == 'ccp4':
            raise Exception('NOT IMPLEMENTED')
            
        cm.print_settings()
        cm.run()
        print '============================>'
        print cm.output
        print '============================>'
        print cm.error
        print '============================>'

if __name__ == '__main__':

    # Show Defaults (just values)
    if '--show-defaults' in sys.argv:
        master_phil.show(attributes_level=0)
    # Show Defaults (including information)
    elif '--help' in sys.argv:
        master_phil.show(attributes_level=2)
    # ... or just run ...
    elif '--expert' in sys.argv:
        master_phil.show(attributes_level=4)
    # ... or just run ...
    else:
        working_phil = parse_phil_args(master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)
        out = run(params=working_phil.extract())
    # Exit (unnecessary, but eh)
    sys.exit()
