#!/usr/bin/env pandda.python

import os, sys, copy, re, shutil

import libtbx.phil
import libtbx.easy_mp
from Bamboo.Common.command import commandManager
from Giant.jiffies import parse_phil_args
from Giant.Xray.Data import crystalSummary

blank_arg_prepend = {'.pdb':'reference_pdb=', '.mtz':'mtz='}

master_phil = libtbx.phil.parse("""
input {
    reference_pdb = None
        .type = path
        .multiple = True
    reference_mtz = None
        .type = path
    mtz = None
        .type = path
        .multiple = True
}
options {
    use_reference_mtz = False
        .type = bool
    cpus = 1
        .type = int
}
output {
    overwrite = True
        .type = bool
    dir_suffix = '_dimple'
        .type = str
}
""")

def run(params):

    def run_dimple(c):
        print '============================>'
        print 'Running New Dimple'
        cm.print_settings()
        ret_code = c.run()
        if ret_code != 0:
            print '============================>'
            print 'Dimple returned with an error'
            print '============================>'
            cm.print_settings()
            print '============================>'
            print cm.output
            print '============================>'
            print cm.error
            print '============================>'
        return ret_code

    ref_mtz_file = params.input.reference_mtz
    if ref_mtz_file: ref_mtz_file = os.path.abspath(ref_mtz_file)

    for ref_pdb_file in params.input.reference_pdb:
            
        commands = []

        ref_pdb_file = os.path.abspath(ref_pdb_file)
        
        if params.options.use_reference_mtz:
            if not ref_mtz_file: ref_mtz_file = os.path.abspath(ref_pdb_file.replace('.pdb','.mtz'))
            assert os.path.exists(ref_mtz_file), 'Reference MTZ does not exist: {}'.format(ref_mtz_file)

        output_dirname = os.path.basename(ref_pdb_file).replace('.pdb', params.output.dir_suffix)
        print 'Placing all dimple runs for {} into directories called {}'.format(ref_pdb_file, output_dirname)

        for mtz_file in params.input.mtz:

            mtz_file = os.path.abspath(mtz_file)

            output_dir = os.path.join(os.path.dirname(mtz_file), output_dirname)

            if os.path.exists(output_dir):
                if params.output.overwrite:
                    shutil.rmtree(output_dir)
                else:
                    raise Exception('Output Directory already exists: {}'.format(output_dir))

            # Need to precreate the directory
            os.mkdir(output_dir)

            cm = commandManager('dimple')
            if params.options.use_reference_mtz:
                cm.add_command_line_arguments([ '--free-r-flags', ref_mtz_file ])
            cm.add_command_line_arguments([ mtz_file,
                                            ref_pdb_file,
                                            output_dir,
                                        ])
            cm.print_settings()
            commands.append(cm)

        # Run Dimples for this reference pdb file
        print '============================>'
        print 'Running Dimple for {} Datasets'.format(len(commands))
        returned = libtbx.easy_mp.pool_map(fixed_func=run_dimple, args=commands, processes=params.options.cpus)

        print '============================>'
        print 'Return Code Total:'
        print sum(returned)

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

