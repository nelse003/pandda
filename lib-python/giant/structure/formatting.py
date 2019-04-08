import iotbx.pdb.hierarchy as iotbx_pdbh


class _Selection:


    _labels = ('model','chain','resseq','icode','resname','altloc','name')


    @classmethod
    def format(cls, obj):
        return cls.join_and([s for s in cls._format_any_to_list(obj) if s not in cls.remove])

    @classmethod
    def join_or(cls, strings, extra_join=''):
        return (cls._join_or+extra_join).join([cls._sep_or.format(s) for s in strings])


    @classmethod
    def join_and(cls, strings, extra_join=''):
        return (cls._join_and+extra_join).join([cls._sep_and.format(s) for s in strings])


    @classmethod
    def join_custom(cls, strings, join):
        if '\n' in join:
            return join.join([s.replace('\n', join) for s in strings])
        else:
            return join.join(strings)


    @classmethod
    def _format_any_to_list(cls, obj):
        """Return the relevant label for a supplied hierarchy/atom object"""
        if   isinstance(obj, iotbx_pdbh.model):
            s = cls._format_mo(obj)
        elif isinstance(obj, iotbx_pdbh.chain):
            s = cls._format_ch(obj)
        elif isinstance(obj, iotbx_pdbh.residue_group):
            s = cls._format_rg(obj)
        elif isinstance(obj, iotbx_pdbh.atom_group):
            s = cls._format_ag(obj)
        elif isinstance(obj, iotbx_pdbh.conformer):
            s = cls._format_co(obj)
        elif isinstance(obj, iotbx_pdbh.residue):
            s = cls._format_re(obj)
        elif isinstance(obj, iotbx_pdbh.atom):
            if hasattr(obj, 'chain_id'): s = cls._format_al(obj)
            else:                        s = cls._format_at(obj)
        elif isinstance(obj, iotbx_pdbh.atom_with_labels):
            s = cls._format_al(obj)
        elif isinstance(obj, dict):
            s = cls._format_dict(obj)
        else:
            raise Exception('Invalid object type provided: {}'.format(type(obj)))
        return s


    @classmethod
    def _format_mo(cls, obj):
        return [ cls.model.format(obj.id) ]


    @classmethod
    def _format_ch(cls, obj):
        return cls._format_any_to_list(obj.parent()) + \
                [ cls.chain.format(obj.id),
                  ]


    @classmethod
    def _format_rg(cls, obj):
        return cls._format_any_to_list(obj.parent()) + \
                [ cls.resseq.format(obj.resseq),
                  cls.icode.format(obj.icode),
                ]


    @classmethod
    def _format_ag(cls, obj):
        return cls._format_any_to_list(obj.parent()) + \
                [ cls.resname.format(obj.resname),
                  cls.altloc.format(obj.altloc),
                ]


    @classmethod
    def _format_co(cls, obj):
        return cls._format_any_to_list(obj.parent()) + \
                [ cls.altloc.format(obj.altloc),
                ]


    @classmethod
    def _format_re(cls, obj):
        return cls._format_any_to_list(obj.parent()) + \
                [ cls.resname.format(obj.resname),
                  cls.resseq.format(obj.resseq),
                  cls.icode.format(obj.icode),
                ]


    @classmethod
    def _format_at(cls, obj):
        return cls._format_any_to_list(obj.parent()) + \
                [ cls.name.format(obj.name),
                ]


    @classmethod
    def _format_al(cls, obj):
        return [cls.model.format(obj.model_id),
                cls.chain.format(obj.chain_id),
                cls.resseq.format(obj.resseq),
                cls.icode.format(obj.icode),
                cls.resname.format(obj.resname),
                cls.altloc.format(obj.altloc),
                cls.name.format(obj.name),
               ]


    @classmethod
    def _format_dict(cls, obj_dict):
        return [ cls.__dict__.get(l).format(obj_dict.get(l))
                 for l
                 in cls._labels if (l in obj_dict.keys()) ]


class Labeller(_Selection):


    _join_and = '-'
    _join_or  = None    # Refmac has no equivalent atom selection syntax
    _sep_and  = '{}'
    _sep_or   = None    # Refmac has no equivalent atom selection syntax

    remove = ['',' ']

    model       = 'mod({})'
    chain       = 'chn({})'
    resseq      = 'res({})'
    icode       = 'ins({})'
    resname     = '{}'
    altloc      = 'alt({})'
    name        = '[{}]'


class ShortLabeller(Labeller):


    model       = ''
    chain       = '{}'
    resseq      = '{}'
    icode       = '{}'
    resname     = '{}'
    altloc      = '({})'
    name        = '[{}]'


class GenericSelection(_Selection):


    _join_and = '-'
    _join_or  = None
    _sep_and  = '{}'
    _sep_or   = None

    remove = ['']

    model       = '{}'
    chain       = '{}'
    resseq      = '{}'
    icode       = '{}'
    resname     = '{}'
    altloc      = '{}'
    name        = '{}'

    @classmethod
    def to_str(cls, obj):
        obj_dict = cls.to_dict(obj)
        obj_list = ['{}({})'.format(l, obj_dict.get(l,'*')) for l in cls._labels]
        return ','.join(obj_list)

    @classmethod
    def to_dict(cls, obj):
        """Format any object to a dictionary (except atoms)"""

        if   isinstance(obj, iotbx_pdbh.model):
            labs = ('model')
            info = cls._format_mo(obj)
        elif isinstance(obj, iotbx_pdbh.chain):
            labs = ('model','chain')
            info = cls._format_ch(obj)
        elif isinstance(obj, iotbx_pdbh.residue_group):
            labs = ('model','chain','resseq','icode')
            info = cls._format_rg(obj)
        elif isinstance(obj, iotbx_pdbh.atom_group):
            labs = ('model','chain','resseq','icode','resname','altloc')
            info = cls._format_ag(obj)
        elif isinstance(obj, iotbx_pdbh.conformer):
            labs = ('model','chain','altloc')
            info = cls._format_co(obj)
        elif isinstance(obj, iotbx_pdbh.residue):
            labs = ('model','chain','altloc','resname','resseq','icode')
            info = cls._format_re(obj)
        elif isinstance(obj, iotbx_pdbh.atom):
            raise Exception('Not implemented')
            labs = ('model','chain','resseq','icode','resname','altloc','name')
            if hasattr(obj, 'chain_id'): info = cls._format_al(obj)
            else:                        info = cls._format_at(obj)
        elif isinstance(obj, iotbx_pdbh.atom_with_labels):
            labs = ('model','chain','resseq','icode','resname','altloc','name')
            info = cls._format_al(obj)
        else:
            raise Exception('Invalid object type provided: {}'.format(type(obj)))

        assert len(labs) == len(info)
        return dict(zip(labs, info))


class RefmacSelection(_Selection):


    _join_and = ' '
    _join_or  = None    # Refmac has no equivalent atom selection syntax
    _sep_and  = '{}'
    _sep_or   = None    # Refmac has no equivalent atom selection syntax

    remove = ['','model ','alte ', 'insc  ']

    model       = 'model {}'
    chain       = 'chain {}'
    resseq      = 'resi {}'
    icode       = 'insc {}'
    resname     = ''
    altloc      = 'alte {}'
    name        = 'atom {}'


class BusterSelection(_Selection):

    """
    Selection syntax for buster

    Attributes
    ----------
    _join_and: ''
        provides string to join other attributes
    _sep_and: {}
        provides formatting seperator,
        needs to include a curly brace format

    Notes
    ------
    Atoms to be formatted as <chain>|<residue_id>:<atom_type>.<altloc>
    """

    _join_and = ''
    _join_or  = None
    _sep_and  = '{}'
    _sep_or   = None

    remove = ['resname','model ','resseq','icode','name']

    model       = '{}'
    chain       = '{}|'
    resseq      = '{}:'
    icode       = '{}'
    resname     = '{}'
    altloc      = '.{}'
    name        = '{}'


    @classmethod
    def format(cls, obj):
        """
        Override to allow formatting as <chain>|<residue_id>:<atom_type>.<altloc>
        """

        atm_description = [s.strip(' ') for s in cls._format_any_to_list(obj) if s not in cls.remove]
        print(atm_description)

        # when obj is an atoms with labels
        if len(atm_description) == 7:
            out = cls.join_and([atm_description[1],
                          atm_description[2],
                          atm_description[6],
                          atm_description[5]])
        else:
            out = cls.join_and([s for s in cls._format_any_to_list(obj) if s not in cls.remove])
        return out


class PhenixSelection(_Selection):


    _join_and = ' and '
    _join_or  = ' or '
    _sep_and  = '{}'
    _sep_or   = '({})'

    remove = ['','model ']

    model       = "model {}"
    chain       = "chain '{:1}'"
    resseq      = "resseq {}"
    icode       = "icode '{:1}'"
    resname     = "resname '{}'"
    altloc      = "altid '{:1}'"
    name        = "name {}"


class PymolSelection(_Selection):


    _join_and = ' and '
    _join_or  = ' or '
    _sep_and  = '{}'
    _sep_or   = '({})'

    remove = ['','model ']

    model       = "model {}"
    chain       = "chain {:1}"
    resseq      = "resi {}"
    icode       = ""
    resname     = "resn '{}'"
    altloc      = "alt '{}'"
    name        = "name {}"


class _Formatter:


    selection = None

    @classmethod
    def _unimplemented(cls, *args, **kw_args):
        raise Exception('Not implemented for this class')

    distance_restraint  = _unimplemented
    peptide_bond        = _unimplemented

    @classmethod
    def format(cls, *args, **kwargs):
        return cls.selection.format(*args, **kwargs)

    @classmethod
    def format_distance_restraints(cls, restraint_list):
        return cls._distance_restraint_format.format(
            cls.selection.join_custom(strings=restraint_list,
                                      join=cls._distance_restraint_format_join))

    @classmethod
    def format_occupancy_restraints(cls, restraint_list):
        return cls._occupancy_restraint_format.format(
            cls.selection.join_custom(strings=restraint_list,
                                      join=cls._occupancy_restraint_format_join))

    @classmethod
    def make_occupancy_restraints(cls, list_of_lists_of_groups, group_completeness=None, idx=1):
        if group_completeness is None:
            group_completeness = [False]*len(list_of_lists_of_groups)
        assert len(group_completeness) == len(list_of_lists_of_groups)
        r_list = []
        for list_of_groups, complete in zip(list_of_lists_of_groups, group_completeness):
            r_list.append(cls.occupancy_restraint(list_of_groups=list_of_groups, complete=complete, idx=idx))
            idx += len(list_of_groups)
        return r_list


class RefmacFormatter(_Formatter):


    selection = RefmacSelection

    _distance_restraint_format = """{}"""
    _distance_restraint_format_join = "\n"
    _distance_restraint = 'exte dist first {} second {} value {} sigma {} type {}'

    @classmethod
    def make_distance_restraint(cls, atm_1, atm_2, value, sigma, add=True):
        s1 = cls.selection.format(atm_1)
        s2 = cls.selection.format(atm_2)
        if add:  add = 1
        else:    add = 0
        return cls._distance_restraint.format(s1, s2, value, sigma, add)

    _occupancy_restraint_format = """{}\noccupancy refine"""
    _occupancy_restraint_format_join = "\n"
    _occupancy_restraint = '{}\noccupancy group alts {} {}'
    _occupancy_restraint_join = '\n'
    _occupancy_group = "occupancy group id {} {}"

    @classmethod
    def occupancy_restraint(cls, list_of_groups, complete=False, idx=1):
        r_list = [cls.occupancy_group(objects=g, idx=idx+i) for i,g in enumerate(list_of_groups)]
        return cls._occupancy_restraint.format(cls.selection.join_custom(strings=r_list, join=cls._occupancy_restraint_join),
                                               'complete' if complete is True else 'incomplete',
                                               ' '.join(map(str,range(idx,idx+len(list_of_groups)))))
    @classmethod
    def occupancy_group(cls, objects, idx):
        return '\n'.join(sorted([cls._occupancy_group.format(idx, cls.selection.format(o)) for o in objects]))


class BusterFormatter(_Formatter):
    """
    Produce formatting for Buster:

    Notes
    ------

    Refmac distance restraint:

    exte dist first chain A resi  110 alte A atom  N   second chain A resi  110 alte C atom  N   value 0.0 sigma 0.02 type 1

    In buster is:

    NOTE BUSTER_DISTANCE =0.0        0.02       A|110:N.A     A|110:N.C

    {}|{}:{}.{}

    chain|residue_id:Atom.altloc

    somewhere for formatting s1 and s1

    Occupancy refinement:

    NOTE BUSTER_RESET_CONSTANT_COMBINE

    This command will Reset the refinement state;
    by default occupancies are all fixed to their initial values,
    NOTE BUSTER_RESET_CONSTANT_COMBINE turns this off.
    It also turns off any grouping of positions or B-factors,
    not sure currently what this implies.

    NOTE BUSTER_SET AltOccAll = Empty

    This defines an empty set of atoms to add all altlocs
    that are being refined to.

    NOTE BUSTER_SET OccZeroH = OccZero & Hydrogen

    This sets the occupancy of hydrogens to zero

    NOTE BUSTER_SET FixOcc = All

    This is a set of all atoms

    NOTE BUSTER_SET FixOcc = FixOcc \ AltOccAll

    This removes the atoms that have been added to the "AltOccAll" set
    from the set of all atoms

    NOTE BUSTER_CONSTANT OCC FixOcc

    This then fixes the occupancy of those atoms,
    which are not in the "AltOccAll" set

    """

    # This instantiates a selection object
    # TODO implement correct selection object
    selection = BusterSelection

    # This is the header line for the whole section of distance restraints.
    # currently empty
    _distance_restraint_format = "{}"
    # join between restraints (new line currently)
    _distance_restraint_format_join = "\n"
    # individual restraint lines
    _distance_restraint = "NOTE BUSTER_DISTANCE = {0}    {1}    {2}    {3}"

    @classmethod
    def make_distance_restraint(cls, atm_1, atm_2, value, sigma, add=True):
        s1 = cls.selection.format(atm_1)
        s2 = cls.selection.format(atm_2)
        if not add:  raise Exception('Not implemented')
        return cls._distance_restraint.format(value, sigma, s1, s2)

    _occupancy_restraint_format = \
"""
NOTE BUSTER_RESET_CONSTANT_COMBINE
NOTE BUSTER_SET AltOccAll = Empty
NOTE BUSTER_SET OccZeroH = OccZero & Hydroen
{}
NOTE BUSTER_SET FixOcc = All
NOTE BUSTER_SET FixOcc = FixOcc \\ AltOccAll
NOTE BUSTER_CONSTANT OCC FixOcc
"""
    _occupancy_restraint_join = "\n"
    _occupancy_restraint_format_join = ""

    # TODO Test with multiple altlocs
    # TODO This needs to be format as D|1:S2.A E|1:S1.B
    _occupancy_restraint = "NOTE BUSTER_OCCSUM 1.0 0.005 {} {}"

    #TODO formatting selection for format A|282:*.D or chain|resid:*.Altloc
    _occupancy_group = \
"""
NOTE BUSTER_SET AltOcc{0} =  {{{1}}}
NOTE BUSTER_SET AltOcc{0} = AltOcc{0} \\ OccZeroH
NOTE BUSTER_SET AltOccAll  = AltOccAll  + AltOcc{0} 
"""


    @classmethod
    def occupancy_restraint(cls, list_of_groups, complete=True, idx=1):
        r_list = [cls.occupancy_group(objects=g, idx=idx + i) for i, g in enumerate(list_of_groups)]
        # TODO How does this work?
        return cls._occupancy_restraint.format(
            cls.selection.join_custom(strings=r_list,
                              join=cls._occupancy_restraint_join),
                              ''.join(map(str, range(idx, idx + len(list_of_groups)))))

    # TODO How does this work?
    @classmethod
    def occupancy_group(cls, objects, idx):
        return '\n'.join(sorted([cls._occupancy_group.format(idx, cls.selection.format(o)) for o in objects]))


class PhenixFormatter(_Formatter):


    selection = PhenixSelection

    _distance_restraint_format = """refinement.geometry_restraints.edits {{\n    {}\n}}"""
    _distance_restraint_format_join = "\n    "
    _distance_restraint = """bond {{
    action = *add
    atom_selection_1 = {}
    atom_selection_2 = {}
    distance_ideal = {}\n    sigma = {}\n    slack = None\n}}"""

    @classmethod
    def make_distance_restraint(cls, atm_1, atm_2, value, sigma, add=True):
        s1 = cls.selection.format(atm_1)
        s2 = cls.selection.format(atm_2)
        if not add:  raise Exception('Not implemented')
        return cls._distance_restraint.format(s1, s2, value, sigma)

    _occupancy_restraint_format = """refinement.refine.occupancies {{\n    {}\n}}"""
    _occupancy_restraint_format_join = "\n    "
    _occupancy_restraint = """constrained_group {{\n    {!s}\n}}"""
    _occupancy_restraint_join = """\n    """
    _occupancy_group = """selection = {!s}"""

    @classmethod
    def occupancy_restraint(cls, list_of_groups, complete=False, idx=1):
        r_list = [cls.occupancy_group(objects=g) for g in list_of_groups]
        return cls._occupancy_restraint.format(cls.selection.join_custom(strings=r_list, join=cls._occupancy_restraint_join))
    @classmethod
    def occupancy_group(cls, objects, idx=0):
        return cls._occupancy_group.format(cls.selection.join_or(sorted([cls.selection.format(obj) for obj in objects]), extra_join='\\\n'+' '*12))

    _peptide_bond = """refinement.pdb_interpretation.apply_cif_link {{\n    data_link = {}
    residue_selection_1 = {}
    residue_selection_2 = {}\n}}"""

    @classmethod
    def peptide_bond(cls, atm_1, atm_2, chain_1=None, chain_2=None, config='TRANS'):
        atm_1 = atm_1.fetch_labels()
        atm_2 = atm_2.fetch_labels()
        if chain_1 is not None: atm_1.chain_id=chain_1
        if chain_2 is not None: atm_2.chain_id=chain_2
        s1 = cls.selection.format(atm_1)
        s2 = cls.selection.format(atm_2)
        return cls._peptide_bond.format(config, s1, s2)

