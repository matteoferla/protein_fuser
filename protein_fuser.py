# prevent pymol from launching in normal mode.
if __name__ == '__main__':
    pymol_argv = ['pymol', '-qc']
else:
    import __main__

    __main__.pymol_argv = ['pymol', '-qc']

import pymol

pymol.finish_launching()

from pprint import PrettyPrinter

pprint = PrettyPrinter().pprint

import os, requests, time, json, uuid
import pandas as pd
import numpy as np
from IPython.display import display, HTML, Image
from Bio.SeqUtils import seq3
from collections import namedtuple
from warnings import warn

##resn = seq3(self.seq[self.end]).upper()

_Coords = namedtuple('Coords', ['name', 'x', 'y', 'z'])


class Coords(_Coords):
    def to_np(self, coords):
        return np.array([coords.x, coords.y, coords.z])


def coords_to_np(self, coords):
    np.array([coords.x, coords.y, coords.z])


####################################################################################
class Workshop:
    """
    >>> ws = Workshop(pdbs, debug).order().save()
    debug: verbosity boolean
    pdbs: [{'x': stated_start_in_uniprot,
            'y': stated_end_in_uniprot,
            'id': 'code_chain', 
            'tier': opt. int. low => max importance (e.g. crystal), high => lower importance (e.g. models)
            'description': opt. 'Not used.'}]
    In the case of PDBs the metadata is taken in order to find where the true start is without cloning scars.
    The code can be a PDB code. or 'PDB:code'  or 'SWISSMODEL:id' or 'LOCAL:filename' (without underscores!)
    For PDB codes, the metadata is fetch to verify what the true values are.
    See model for more about the number madness
    order() readies and aligns or projects the non-redudant models.
    the init step actually calls `.categorise()` which categorises the models based on redundancy.
    ready_models() loads models
    align_models() aligns overlapping ones
    project_models() places models with gaps along the x-axis.
    `ws.models` list of models.
    Note that there is not a create step to fuse the parts as simply saving as pdb will do it.
    But saving as pse to check things is nice.
    `pymol.exporting.multisave('tmp.pdb')` is a handy command if there is a version mismathc with pymol.
    
    Model has many bound methods that use pymol cmds, that affect only one model, say `protein.roll(40)`,
    while workshop uses more thatn one say, `ws.roll_free(protein_A, protein_B)`.
    Some methods seems repeated but aren't. the `_models` operate on all.
    ws.ready_models() --> model.fetch_n_clean()
    ws.align_models() --> ws.align(A,B)
    ws.project_models() --> ws.project(A,B)
    """

    def __init__(self, pdbs, debug=False):
        self.debug = bool(debug)
        self.pdbs = [{'tier': 1, 'description': 'NA', **p} for p in pdbs]
        self.metadata = {}
        self.fetch_metadata()
        if self.debug:
            print(self.metadata)
        self.models = []
        self.catergorise()
        ### pymol operations.
        self.last_het_index = 9000
        self.joints = []  ##{'joint': 'resi_fore': int, 'fore_model': 'aft'..., 'type': 'end|alignment|projection', 'overlap' distance 'note'}

    def log_joint(self, form='ERROR', fore=None, aft=None):
        if fore is None and aft is None:
            raise ValueError
        elif fore is None:  # nothing before
            self.joints.append({'type': 'end',
                                'aft_resi_start': aft.start,
                                'aft_source': f'{aft.mode}:{aft.code}'})
        elif aft is None:  # nothing after
            self.joints.append({'type': 'end',
                                'fore_resi_start': fore.end,
                                'fore_source': f'{fore.mode}:{fore.code}'})
        else:
            f = fore.get_resi_coords(fore.end)['C']
            a = aft.get_resi_coords(aft.start)['N']
            d = fore.distance(f, a)
            o = pymol.cmd.overlap(f"model {fore.name}", f"model {aft.name}")
            self.joints.append({'type': form,
                                'fore_resi_end': fore.end,
                                'fore_source': f'{fore.mode}:{fore.code}',
                                'aft_resi_start': aft.start,
                                'aft_source': f'{aft.mode}:{aft.code}',
                                'distance': d,
                                'overlap': o
                                })
        return self

    def order(self):
        """
        clear pymol, ready models, align overlapping, project gaps.
        """
        pymol.cmd.delete('all')
        self.ready_models()
        pymol.cmd.delete('cylinder*')
        if self.debug:
            self.show_load_checks()
            self.show_before()
        if self.debug:
            pymol.cmd.zoom(self.models[0].name)
            self.show_pose('mid.png')
            print('### Aligning!')
        self.align_models()
        if self.debug:
            print('### Projecting!')
        self.project_models()
        self.log_joint('end', aft=self.models[0])
        self.log_joint('end', fore=self.models[-1])
        if self.debug:
            self.show_after()
        return self

    def save(self, filename='fusion.pdb'):
        pymol.cmd.save(filename)
        return self

    def align_models(self):
        """
        Fuses all overlapping models. The name will be a uuid of 6 letters
        """
        for i in range(1, len(self.models)):
            fore = self.models[i - 1]
            aft = self.models[i]
            if fore.end > aft.start:
                self.align(fore, aft)
                self.log_joint('alignment', fore, aft)
                # changing aft to a fusion.
                self.models[i].start = fore.start
                self.models[i].length = fore.length + aft.length
                new_name = uuid.uuid1().hex.replace('-', '')[:6]
                pymol.cmd.create(new_name, f'{fore.name} or {aft.name}')
                pymol.cmd.delete(fore.name)
                pymol.cmd.delete(aft.name)
                self.models[i].name = new_name
                self.models[i - 1] = None  ## kill fore from list.
        self.models = [m for m in self.models if m is not None]
        return self

    def project_models(self):
        """
        places models with gaps along the x-axis.
        """
        for m in self.models:  ##no idea why this is needed.
            m.angle_fix()
        for i in range(1, len(self.models)):
            fore = self.models[i - 1]
            aft = self.models[i]
            self.project(fore, aft)
            self.roll_free(fore, aft)
            self.log_joint('projection', fore, aft)
        return self

    def show_pose(self, filename='temp.png'):
        """
        Jupyter notebook specific.
        """
        if '.png' not in filename:
            filename += '.png'
        if os.path.exists(filename):
            os.remove(filename)
        pymol.cmd.png(filename)
        while not os.path.exists(filename):
            time.sleep(1)
        display(Image(filename=filename))

    def show_before(self):
        self.add_xyz()
        modsele = f'{self.models[0].sele(self.models[0].end)} or {self.models[1].sele(self.models[1].start)} or {self.models[1].sele(self.models[1].end)}'
        pymol.cmd.color('yellow', modsele)
        pymol.cmd.show(representation="sticks", selection=modsele)
        pymol.cmd.zoom(modsele)
        print('Before')
        self.show_pose('before.png')

    def show_after(self):
        print('after')
        pymol.cmd.zoom('all')
        self.show_pose('after.png')
        pymol.cmd.save('after.pdb', 'all')

    def show_load_checks(self):
        pymol.cmd.zoom('all')
        for model in self.models:
            pymol.cmd.color('cyan', model.name)
            pymol.cmd.color('red', f'{model.name} and chain {model.chain}')
            modsele = model.sele(f'{model.start}-{model.end}')
            pymol.cmd.color('white', modsele)
        self.show_pose('all.png')

    def roll_free(self, fore, aft, cycle=1):
        if cycle == 0:
            return self
        if self.has_overlap(fore, aft):
            for i in range(6):
                print(f'Overlap between {fore.name} and {aft.name}')
                aft.roll(60)
                if not self.has_overlap(fore, aft):
                    return self
            else:
                print(f'Failed cycle {cycle}')
                self.project(fore, aft, extra_d=5 * cycle)
                self.roll_free(fore, aft, cycle=cycle + 1)
        return self

    def align(self, fore, aft, cycle=0):
        mobile = aft.sele(f"{aft.start}-{fore.end}")
        m_match = self.safe_select(mobile)
        target = fore.sele(f"{aft.start}-{fore.end}")
        t_match = self.safe_select(target)
        if m_match == 0 or t_match == 0:
            warn(f'Possible major issue...\n{mobile} ... {m_match}\n{target} ... {t_match}')
            return self
        pymol.cmd.align(mobile, target)
        scores = {}
        if fore.tier <= aft.tier:
            ranger = range(fore.end, aft.start - 1, -1)
        else:
            ranger = range(aft.start, fore.end + 1, +1)
        for i in ranger:
            pymol.cmd.pair_fit(aft.sele(i, 'N'), fore.sele(i, 'N'),
                               aft.sele(i, 'CA'), fore.sele(i, 'CA'),
                               aft.sele(i, 'C'), fore.sele(i, 'C'))
            if self.debug:
                print(f'... {mobile} --> {target}')
                pymol.cmd.zoom(mobile)
                pymol.cmd.color('red', mobile)
                pymol.cmd.color('green', target)
                self.show_pose('align.png')
            if cycle > 0:
                aft.roll(60 * cycle)
            o = pymol.cmd.overlap(fore.sele(f"{fore.start}-{i}"), aft.sele(f"{i + 1}-{aft.end}"))
            if o == 0:
                break
            else:
                if self.debug:
                    print(f'!!! There is an overlap of {o}!')
                scores[i] = o
        else:
            i = list(scores.keys())[list(scores.values()).index(min(scores.values()))]
            if scores[i] < 100:
                if self.debug:
                    print(f'Aligning via residue {i}. score={scores[i]}')
                pymol.cmd.pair_fit(fore.sele(i, 'N'), aft.sele(i, 'N'),
                                   fore.sele(i, 'CA'), aft.sele(i, 'CA'),
                                   fore.sele(i, 'C'), aft.sele(i, 'C'))
                if cycle > 0:
                    aft.roll(60 * cycle)
            else:
                print('Good alignment impossible.')
                if cycle < 6:
                    self.align(fore, aft, cycle=cycle + 1)
                else:
                    raise Exception('Really?? Overlap after alignment even with rotations?? Sigh.')
        if fore.tier <= aft.tier:
            pymol.cmd.remove(fore.sele(f'{i + 1}-{fore.end}'))
            fore.end = i
            pymol.cmd.remove(aft.sele(f'{aft.start}-{i}'))
            aft.start = i + 1
        else:
            pymol.cmd.remove(fore.sele(f'{i}-{fore.end}'))
            fore.end = i - 1
            pymol.cmd.remove(aft.sele(f'{aft.start}-{i - 1}'))
            aft.start = i
        return self

    def fetch_metadata(self):
        if self.debug: print('## Fetching metadata')
        ### fetch the metadata
        for entry in self.pdbs:
            mode, code, chain = Model.split_id(entry)
            if mode == 'PDB':
                reply = requests.get(f'http://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{code}').json()
                self.metadata[code] = reply[code.lower()]
            else:
                self.metadata[code] = {}
        return self

    def catergorise(self):
        """
        Clusters the proteins
        """
        if self.debug: print('## Running categorisation')
        ## find longest!
        tdata = []
        from collections import defaultdict
        starters = defaultdict(list)
        enders = defaultdict(list)
        for entry in self.pdbs:
            print(entry)
            mode, code, chain = Model.split_id(entry)
            if mode == 'PDB':
                c_entity = [entity for entity in self.metadata[code] if chain in entity['in_chains']][0]
                mappings = c_entity['source'][0]['mappings']
                if len(c_entity['source'][0]['mappings']) > 1:
                    raise ValueError('MULTIPLE MAPPINGS!')
                s = mappings[0]['start']['residue_number']
                e = mappings[0]['end']['residue_number']
                seq = c_entity['sequence']
            else:
                # Model.infer_metadata(entry)
                # burden on user.
                s = 1
                e = entry['y'] - entry['x'] + 1
                seq = None
            tier = float(entry['tier'])
            tdata.append({'code': code,
                          'mode': mode,
                          'chain': chain,
                          'tier': tier,
                          'x': entry['x'],
                          'y': entry['y'],
                          'seq': seq,
                          'start': s,
                          'end': e,
                          'length': e - s})
        df = pd.DataFrame.from_records(tdata).sort_values(['x', 'length', 'tier'], ascending=[True, False, True])
        if self.debug:
            print('There are all the chains:')
            display(df)
        ## thin!
        new = df.iloc[0]
        get_similars = lambda: list(set(df.loc[(df.y <= new.y) & (df.x >= new.x)].code.values) - {new.code})
        new.at['similars'] = get_similars()
        o = df.loc[df.y > new.y]
        parts = [new]
        old = new
        while len(o):
            new = o.iloc[0]
            new.at['similars'] = get_similars()
            if new.x < old.y:  # mark overlap for debug.
                new.at['subtender'] = old.code
            ## curious case of better models within makes this weird.
            # so something is better is pulled out of the similars (i.e. smaller so excluded.) group.
            # but to be better it has to be of a lower tier and not present in the selection.
            betters = df[(df.code.isin(new.similars)) & (df.tier < new.tier) & (~df.code.isin([p.code for p in parts]))]
            bi = len(betters.index)
            if bi == 0:
                parts.append(new)  ## has to happen!!
            else:
                better = betters.iloc[0]
                print(f'Nested case {new.code} and {betters.code.values}')
                ## dangerous hack. It assumes that a good match can be found and won't overrun.
                ## as the priority differs it will chew off one side then the other, hence the double.
                if better.x > new.x:
                    parts.append(new)
                parts.append(better)
                if better.y < new.y:
                    parts.append(new)
            ## cleanup.
            o = o.loc[o.y > new.y]
            old = new
        if self.debug:
            print('These are the non-redudant chains:')
            display(pd.DataFrame(parts))
        self.models = [Model(workshop=self, metadata=self.metadata[row.code], **row) for row in parts]
        return self

    def ready_models(self):
        """
        Runs model.fetch_n_clean() on all models.
        """
        if self.debug: print('## Reading models')
        for model in self.models:
            model.fetch_n_clean()
        self.models[0].camera_fix()
        return self

    def project(self, fore, aft, extra_d=0):
        """extra_d is to move it away when there is on clean rotational solution."""
        aft.angle_fix()
        (ref, pointer, θ, φ, internal_d) = fore.angles
        (ref_aft, pointer_aft, θ_aft, φ_aft, internal_d_aft) = aft.angles
        ## translate
        current_end = aft.get_resi_vector(aft.start)['N']  # re-get as it moved.
        missing = aft.end - fore.start
        missing_d = abs(missing * 3.5) ** 0.5  # the James Ross fudge distance
        if self.debug: print(f'Projecting {missing_d} + {extra_d} Å')
        d = missing_d + internal_d + extra_d
        target_end = np.array([d * np.cos(-θ) * np.cos(φ) + ref.x,
                               d * np.sin(-θ) * np.cos(φ) + ref.y,
                               d * np.sin(-θ) * np.sin(φ) + ref.z])
        pymol.cmd.translate((target_end - current_end).tolist(), aft.name, camera=0)
        if self.debug:
            fore.add_NC_arrow('cylinderF', r=0, g=0, b=1)
            aft.add_NC_arrow('cylinderA', r=0, g=1, b=0)

    def has_overlap(self, fore, aft):
        o = pymol.cmd.overlap(f"model {fore.name} and chain {fore.chain}", f"model {aft.name} and chain{aft.chain}")
        if o > 0:
            if self.debug:
                print(f"Overlap between {fore.name} chain {fore.chain} and {aft.name} chain{aft.chain} of {o}")
            return True
        o = pymol.cmd.overlap(f"model {fore.name}", f"model {aft.name}")
        if o > 0:
            if self.debug:
                print(f"Overlap between {fore.name} (whole) and {aft.name} (whole)")
            return True
        return False

    def add_xyz(self):
        pymol.cmd.load_cgo([9.0, 0, 0, 0, 100, 0, 0, 3, 1, 1, 1, 0, 0, 0], "cylinderx")
        pymol.cmd.load_cgo([9.0, 0, 0, 0, 0, 100, 0, 3, 1, 1, 1, 0, 0, 0], "cylindery")
        pymol.cmd.load_cgo([9.0, 0, 0, 0, 0, 0, 100, 3, 1, 1, 1, 0, 0, 0], "cylinderz")

    def safe_select(self, name, sele=None, silent=False):
        """
        There is something not quite right with the selection. It throws an error..."""
        try:
            if not sele:
                sele = name
                name = 'test'
            n = pymol.cmd.select(name, sele)
            return n
        except:
            if '{' in sele:
                raise ValueError(f'There seems to be a curly bracket in the selection ({sele}): bad formatting?')
            else:
                if not silent:
                    warn(f'Failed selection "{sele}"')
                return 0


###############################################################################################

class Model:
    def __init__(self, mode, chain, code, x, y, start=1, end=None, seq=None, length=None, workshop=None, tier=0,
                 metadata=None, similars=None, subtender=None):
        """
        Init just stores the variables. but .fetch_and_clean() does all the work by calling:
        * model.kill_states()     ---> use only state 1.
        * model.move_away_hetatm()  ---> protect the hetatms ligands!
        * model.clean_Nterminus()   ---> trim the linkers!
        * model.clean_Cterminus()
        * model.fix_numbers()   --> the numbering off quite frequently...
        * model.shift_to_true_Nterminus()  ---> gets rid of unsolved terminus
        * model.shift_to_true_Cterminus()  ---> gets rid of unsolved terminus
        While other methods are called by the project method of workshop. Say all the angle operations.
        e.g.model.roll(n) or model.get_resi_coords(resi)
        
        The code can be a PDB code. or... PDB:code  or SWISSMODEL:id or LOCAL:filename (without underscores!)
        
        There are three sets of numbers.
        uniprot start (`x`) which is the position in the whole protein of the first real residue based on the stated sequence regardless of density.
        start (`start`) which is the first residue stated in PDBe query based on the stated sequence regardless of density.
        Namely... the mapping start is at what residue index within the model does the match start,
        while uniprot start is what the position that residues has in the whole sequence.
        A method to keep an eye for is `.angle_fix()`. This ensures that the the N and C termini are on the x axis.
        """
        self.chain = chain
        self.mode = mode
        self.code = code
        self.name = code  ## a backup value in case of duplicates...
        self.tier = tier
        self.metadata = metadata
        self.x = x
        self.y = y
        self.start = start
        if end:
            self.end = end
        else:
            self.end = y - x
        if length:
            self.length = length
        else:
            self.length = y - x
        self.seq = seq
        self.subtender = subtender
        if workshop:
            self.workshop = workshop
        else:
            warn('No workshop?')
            self.workshop = namedtuple('fakeshop', ['debug'])(False)
        self.length = length

    @classmethod
    def from_ds(cls, ds):
        cls(**ds)

    def sele(self, resi, atom=None, no_hetatm=True):
        if isinstance(resi, float):
            resi = int(resi)
        if no_hetatm:
            extra = ' and not hetatm'
        else:
            extra = ''
        if atom:
            return f'/{self.name}//{self.chain}/{resi}/{atom}' + extra
        else:
            return f'/{self.name}//{self.chain}/{resi}/' + extra

    def __str__(self):
        return f'<Model name={self.name} chain={self.chain} start={self.start} end={self.end}>'

    def describe(self):
        (ref, pointer, θ, φ, internal_d) = self.angles
        print(f'Angles: θ={np.rad2deg(θ)} and φ={np.rad2deg(φ)})')
        print('Nterm:', self.get_resi_coords(self.start)['N'])
        print('Cterm:', self.get_resi_coords(self.end)['C'])

    @staticmethod
    def split_id(ID):
        """
        splits an entry id into mode, code, chain
        """
        if isinstance(ID, dict):
            ID = ID['id']
        if ':' not in ID:
            mode = 'PDB'
        else:
            mode, ID = ID.split(':')
        if '_' in ID:
            code, chain = ID.split('_')
        else:
            code, chain = ID, 'A'
        return mode, code, chain

    @classmethod
    def infer_metadata(cls, ID):
        mode, code, chain = cls.split(ID)
        m = cls(mode=mode, chain=chain, code=code, start=1, end=1, x=1, y=1, seq='', length=1).load()
        residex = m.resi_list()
        pymol.cmd.delete(m.name)
        return {'start': residex[0], 'end': residex[-1], 'length': residex[-1] - residex[0]}

    def resi_list(self):
        return sorted(set([int(at.resi) for at in pymol.cmd.get_model(self.name).atom]))

    def residues(self):
        rset = set([(at.resn, int(at.resi), True) for at in pymol.cmd.get_model(self.name).atom])
        r = sorted(rset, key=lambda x: x[1])
        s = ''.join(pymol.cmd.get_fastastr(f'{self.name} and chain {self.chain} and not hetatm').split('\n')[1:])
        first_stated_resn = seq3(s[0]).upper()
        first_structured_resi = r[0][1]
        first_structured_resn = r[0][0]

        if first_structured_resi == 1:
            if first_stated_resn == first_structured_resn:  ## case 1. all is good.
                pass
            else:  ## case 2. negative numbers??!
                raise NotImplementedError(f'First residues are {r[0:3]} but the sequence is {s}')
        elif first_structured_resi > 1:
            if first_stated_resn == first_structured_resn:  ## case 3. The first residue exist is not 1.
                pass
            elif seq3(s[first_structured_resi - 1]).upper() == first_structured_resn:
                ## case 4. the first stated residue is 1, but does not exist.
                r = [(seq3(s[i]).upper(), i + 1, False) for i in range(1, first_structured_resi)] + r
            else:
                raise NotImplementedError(f'First residues are {r[0:3]} but the sequence is {s}')
        else:
            raise NotImplementedError(f'First residues are {r[0:3]} but the sequence is {s}')

        return r

        ###################### METHODS FOR LOADING ####################################

    def fetch_n_clean(self):
        if self.workshop.debug: print(f'### Fetching for {self.name}')
        c = self.code
        self.load()
        self.force_chain_A()
        self.start_from_one()
        self.move_away_hetatm()
        self.clean_Nterminus()
        self.clean_Cterminus()
        self.fix_numbers()
        self.shift_to_true_Nterminus()
        self.shift_to_true_Cterminus()

    def load(self):
        print(f'loading {self.name}')
        if self.name in pymol.cmd.get_object_list():
            preexisting = self.name
            self.name += 'X'
            warn(f'{preexisting} is already loaded. Copying as {self.name}.')
            pymol.cmd.create(self.name, preexisting)
            return self
        if self.mode == 'PDB':
            print('Fetching biological assembly.')
            self.name = pymol.cmd.fetch(self.code, type='pdb1')
            if self.name == -1:
                self.name = pymol.cmd.fetch(self.code)
                print('\t... solved by fetching assimetric unit.')
            if self.name == -1:
                raise ValueError(f'Could not load {self.code}')
            self.kill_states()
        elif self.mode == 'SWISSMODEL':
            url = f'https://swissmodel.expasy.org/repository/{self.code}.pdb'
            if self.workshop.debug:
                print(url)
            x = requests.get(url).text
            if 'END' in x:
                pymol.cmd.read_pdbstr(x, self.code)
            else:
                raise ValueError(f'{url} failed.')
            chains = sorted({at.chain for at in pymol.cmd.get_model(self.name).atom})
            if len(chains) == 1:
                if chains[0] == '':
                    pymol.alter.cmd(self.name, 'chain="A"')
                    self.chain = 'A'
                else:
                    self.chain = chains[0]
            else:
                raise NotImplementedError
        elif self.mode == 'LOCAL':
            pymol.cmd.load(self.code)
            self.name = os.path.splitext(self.code)[0]
        else:
            raise ValueError(f'UNKNOWN mode {self.mode}. PDB: SWISSMODEL: LOCAL:')
        if self.seq is None:
            self.seq = ''.join(
                pymol.cmd.get_fastastr(f'{self.name} and chain {self.chain} and not hetatm').split('\n')[1:])
        return self

    def force_chain_A(self):
        get_untaken = lambda c: list(set('ABCDEFGHIJKLMNOPQRSTUVWXYZ').difference(c))[0]
        if self.chain != 'A':
            chains = sorted({at.chain for at in pymol.cmd.get_model(self.name).atom})
            if 'A' in chains:
                # move chain A away.
                new = get_untaken(chains)
                pymol.cmd.alter(f'model {self.name} and chain A', f'chain="{new}"')
                pymol.cmd.sort()
            pymol.cmd.alter(f'model {self.name} and chain {self.chain}', f'chain="A"')
            pymol.cmd.sort()
        segis = sorted({at.chain for at in pymol.cmd.get_model(f'model {self.name} and chain A').atom})
        if len(segis) > 1:
            if '' not in segis:
                subs = segis[1:]
            else:
                i = segis.index('')
                subs = segis[:i] + segis[i + 1:]
            for s in subs:
                new = get_untaken(chains)
                pymol.cmd.alter(f'model {self.name} and chain A and segi {s}', f'chain="{new}"')
                pymol.cmd.sort()
        elif segis[0] != '':
            pymol.cmd.alter(f'model {self.name} and chain A and segi {segis[0]}', f'segi=""')
            pymol.cmd.sort()
        else:
            pass
        self.chain = 'A'
        return self

    def kill_states(self):
        # fetch(.., state=1) does not work. !?
        pymol.cmd.split_states(self.name, 1, 1)
        pymol.cmd.set_name(f'{self.name}_0001', self.code)
        self.name = self.code
        return self

    def move_away_hetatm(self):
        if self.workshop.debug: print(f'### Protecting hetatms for {self.name}')
        """
        Renumbers the hetatm residues to 9000+ values, unique in the whole workshop.
        """
        hetatms = {self.sele(at.resi, no_hetatm=False) for at in pymol.cmd.get_model(self.name).atom if at.hetatm}
        for h in hetatms:
            pymol.cmd.alter(h, f"resi={self.workshop.last_het_index}")
            self.workshop.last_het_index += 1
        pymol.cmd.sort()
        return self

    def start_from_one(self):
        """
        The first residue may exist or not.
        """
        residues = self.residues()
        if residues[0][1] != 1:
            Δ = residues[0][1] - 1
            pymol.cmd.alter(f'{self.name} and chain {self.chain}', f'resi=str(int(resi)-{Δ})')
            pymol.cmd.sort()
        return self

    def clean_Nterminus(self):
        if self.workshop.debug: print(f'### NTerm cleaning for {self.name}')
        pymol.cmd.remove(self.sele(f'0-{self.start - 1}'))
        pymol.cmd.sort()
        return self

    def clean_Cterminus(self):
        if self.workshop.debug: print(f'### CTerm cleaning for {self.name}')
        """
        convert the N of the junk resi to OXT.
        """
        has_trailing = pymol.cmd.alter(self.sele(self.end + 1, 'N'), 'name="OXT"')
        if has_trailing:
            pymol.cmd.alter(self.sele(self.end + 1, 'N'), f'resi="{self.end}"')
            pymol.cmd.sort()
            resn = pymol.cmd.get_model(self.sele(self.end)).atom[0].resn
            pymol.cmd.alter(self.sele(self.end + 1, 'OXT'), f'resn="{resn}"')
            pymol.cmd.sort()
        pymol.cmd.remove(self.sele(f'{self.end + 1}-9999'))
        pymol.cmd.sort()
        return self

    def fix_numbers(self):
        if self.workshop.debug:
            print(f'### Fixing numbering for {self.name}')
        if self.start != self.x:
            Δ = self.x - self.start
            pymol.cmd.alter(f'{self.name} and chain {self.chain}', f'resi=str(int(resi) + {Δ})')
            pymol.cmd.sort()
            self.start = self.start + Δ
            self.end = self.end + Δ
        return self

    def shift_to_true_Nterminus(self):
        """
        The model may have missing residues at the NTerm, but in reality it is as if they were not there.
        """
        if self.workshop.debug: print(f'### Shifting Nterm for {self.name}')
        ## verify that the start exists.
        coords = None
        l = self.start
        while coords is None:
            coords = pymol.cmd.get_coords(self.sele(l))
            l += 1
        self.start = l - 1

    def shift_to_true_Cterminus(self):
        """
        The model may have missing residues at the CTerm, but in reality it is as if they were not there.
        """
        if self.workshop.debug: print(f'### Shifting Cterm for {self.name}')
        ## verify that the end exists.
        coords = None
        l = self.end
        while coords is None:
            coords = pymol.cmd.get_coords(self.sele(l))
            l -= 1
            if l < 0:
                raise ValueError(
                    f'Trimming missing residues went haywire! start={self.x}/{self.start} end={self.y}/{self.end}. Died at {self.sele(l)}')
        self.end = l + 1

    ###################### METHODS FOR ADJUSTING ####################################
    def get_resi_coords(self, resi):
        match_coords = lambda c: {'x': c[0], 'y': c[1], 'z': c[2]}
        coordsdex = self.get_resi_vector(resi)
        return {atom: Coords(name=atom, **match_coords(coordsdex[atom]))
                for atom in coordsdex}

    def distance(self, A, B):
        return ((A.x - B.x) ** 2 + (A.y - B.y) ** 2 + (A.z - B.z) ** 2) ** 0.5

    def get_resi_vector(self, resi):
        # if an error occurs here, it is because there is no self.sele(resi), which should not be possible
        # as it is dealt with upstream in the trimming. Export the scene with pymol.exporting.multisave('tmp.pdb')
        # as see what is wrong. It is somethign weird. Maybe a hetatm not marked as such?
        return {atom: pymol.cmd.get_coords(self.sele(resi, atom))[0] for atom in ('C', 'N', 'CA')}

    @property
    def angles(self):
        ref = self.get_resi_coords(self.start)['N']
        pointer = self.get_resi_coords(self.end)['C']
        internal_d = self.distance(pointer, ref)
        θ = np.arctan((pointer.y - ref.y) / (pointer.x - ref.x))
        φ = np.arctan((pointer.z - ref.z) / (pointer.x - ref.x))
        return (ref, pointer, θ, φ, internal_d)

    def add_NC_arrow(self, name, r=0, g=0, b=0):
        N = self.get_resi_coords(self.start)['N']
        C = self.get_resi_coords(self.end)['C']
        pymol.cmd.load_cgo([9.0, N.x, N.y, N.z, C.x, C.y, C.z, 3, r, g, b, r / 3, g / 3, b / 3], name)

    def camera_fix(self):
        # pymol.cmd.alter('all','segi=""')
        self.angle_fix()

    def angle_fix(self, cycles=3):
        if self.workshop.debug:
            print(f'Initial position for model {self.name}.....')
            self.describe()
        ## rotate the xy plane
        for i in range(cycles):
            (ref, pointer, θ, φ, internal_d) = self.angles
            pymol.cmd.translate([-ref.x, -ref.y, -ref.z], self.name, camera=0)
            pymol.cmd.rotate([0, 0, 1], -np.rad2deg(θ), selection=self.name, camera=0)
            ## rotate the xz plane
            (ref, pointer, θ, φ, internal_d) = self.angles
            pymol.cmd.translate([-ref.x, -ref.y, -ref.z], self.name, camera=0)
            pymol.cmd.rotate([0, 1, 0], np.rad2deg(φ), selection=self.name, camera=0)
            (ref, pointer, θ, φ, internal_d) = self.angles
            pymol.cmd.translate([-ref.x, -ref.y, -ref.z], self.name, camera=0)
            tail = self.get_resi_coords(self.end)['C']
            if tail.x < 0:
                if self.workshop.debug:
                    print('It is off by 180... Odd.')
                pymol.cmd.rotate([0, 0, 1], 180, selection=self.name, camera=0)  ## rotate the xy plane
                (ref, pointer, θ, φ, internal_d) = self.angles
                pymol.cmd.translate([-ref.x, -ref.y, -ref.z], self.name, camera=0)
        if self.workshop.debug:
            print(f'Final position for model {self.name}.....')
            self.describe()
        return self

    def roll(self, φ):
        """
        Whereas they are aligned on the x-axis with the y,z coordinates of C and N termini near 0,0.
        They can rotate around the X. roll. This is good to remove clashes.
        """
        pymol.cmd.rotate([1, 0, 0], φ, selection=self.name, camera=0)
        self.angle_fix()
        return self
