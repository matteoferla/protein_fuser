from __future__ import annotations
import os
import requests
import time
import uuid
from collections import namedtuple
from pprint import PrettyPrinter
from typing import List, Optional, Tuple
from warnings import warn

import numpy as np
import pandas as pd
from Bio.SeqUtils import seq3
from IPython.display import display, Image

pprint = PrettyPrinter().pprint

# ======================================================================================================================

_Coords = namedtuple('Coords', ['name', 'x', 'y', 'z'])


class Coords(_Coords):
    def to_np(self):
        return np.array([self.x, self.y, self.z])


# ======================================================================================================================

class Fuser:
    """
    >>> fuser = Fuser(pdbs=[{},{}], debug=False).order().save()
    debug: verbosity boolean
    pdbs: [{'uniprot_start': stated_start_in_uniprot,
            'uniprot_end': stated_end_in_uniprot,
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
    But saving as pse to check things is nice. `pymol.cmd.save('test.pse')`.
    `pymol.exporting.multisave('tmp.pdb')` is a handy command if there is a version mismatch with pymol.
    
    Model has many bound methods that use pymol cmds, that affect only one model, say `protein.roll(40)`,
    while fuser uses more thatn one say, `ws.roll_free(protein_A, protein_B)`.
    Some methods seems repeated but aren't. the `_models` operate on all.
    ws.ready_models() --> model.fetch_n_clean()
    ws.align_models() --> ws.align(A,B)
    ws.project_models() --> ws.project(A,B)
    """
    pymol = None

    def __init__(self, pdbs: List[dict], debug: bool = False):
        self.debug = bool(debug)
        # adding tier and description defaults
        self.pdbs = [{'tier': 1, 'description': 'NA', **p} for p in pdbs]
        self.metadata = {}
        self.fetch_metadata()
        if self.debug:
            print(self.metadata)
        self.models = []
        self.categorise()
        # pymol operations.
        self.last_het_index = 9000
        self.joints = []
        self.old2new = {}
        # {'joint': 'resi_fore': int, 'fore_model': 'aft'...,
        # 'type': 'terminus|alignment|projection', 'overlap', distance 'note'}

    def log_joint(self, form='ERROR', fore=None, aft=None):
        if fore is None and aft is None:
            raise ValueError
        elif fore is None:  # nothing before
            self.joints.append({'type': 'terminus',
                                'aft_resi_start': aft.pdb_start,
                                'aft_source': f'{aft.mode}:{aft.code}'})
        elif aft is None:  # nothing after
            self.joints.append({'type': 'terminus',
                                'fore_resi_start': fore.pdb_end,
                                'fore_source': f'{fore.mode}:{fore.code}'})
        else:
            try:
                f = fore.get_resi_coords(fore.pdb_end)['C']
                a = aft.get_resi_coords(aft.pdb_start)['N']
                d = fore.distance(f, a)
            except:
                f = float('nan')
                a = float('nan')
                d = float('nan')
            o = self.pymol.cmd.overlap(f"model {fore.name}", f"model {aft.name}")
            self.joints.append({'type': form,
                                'fore_resi_end': fore.pdb_end,
                                'fore_source': f'{fore.mode}:{fore.code}',
                                'aft_resi_start': aft.pdb_start,
                                'aft_source': f'{aft.mode}:{aft.code}',
                                'distance': d,
                                'overlap': o
                                })
        return self

    def order(self):
        """
        clear pymol, ready models, align overlapping, project gaps.
        """
        self.pymol.cmd.delete('all')
        self.ready_models()
        self.pymol.cmd.delete('cylinder*')
        if self.debug:
            self.show_load_checks()
            self.show_before()
        if self.debug:
            self.pymol.cmd.zoom(self.models[0].name)
            self.show_pose('mid.png')
            print('# Aligning!')
        self.align_models()
        if self.debug:
            print('# Projecting!')
        self.project_models()
        self.log_joint('terminus', aft=self.models[0])
        self.log_joint('terminus', fore=self.models[-1])
        if self.debug:
            self.show_after()
        return self

    def save(self, filename='fusion.pdb'):
        self.pymol.cmd.create('fusion', ' or '.join(self.pymol.cmd.get_names(selection='polymer')))
        self.pymol.cmd.save(filename, 'fusion')
        return self

    def get_model_pair(self, i) -> Tuple[Model, Model]:
        fore = self.models[i]
        aft = self.models[i + 1]
        if aft is None and aft.name in self.old2new:
            warn(f'deleted object requested')
            aft = self.old2new[aft.name]
        if fore is None and fore.name in self.old2new:
            warn(f'deleted object requested')
            fore = self.old2new[fore.name]
            if fore.name == aft.name:
                warn('Fore and aft are the same: ignoring.')
                raise ValueError
            else:
                self.models[i] = fore
        if fore.name == aft.name:
            warn('Fore and aft are the same: ignoring.')
            raise ValueError
        else:
            self.models[i] = fore
            self.models[i + 1] = aft
            return fore, aft

    def align_models(self):
        """
        Fuses all overlapping models. The name will be a uuid of 6 letters
        """
        for i in range(0, len(self.models) - 1):
            try:
                fore, aft = self.get_model_pair(i)
            except ValueError:
                continue
            if fore.pdb_end > aft.pdb_start:
                self.align(fore, aft)
                self.log_joint('alignment', fore, aft)
                # changing aft to a fusion.
                new_name = uuid.uuid1().hex.replace('-', '')[:6]
                if self.debug:
                    print(f'Combining {fore.name} and {aft.name} => {new_name}')
                # make aft the combo.
                aft.pdb_start = fore.pdb_start
                aft.length = fore.length + aft.length # wrong...
                self.pymol.cmd.create(new_name, f'{fore.name} or {aft.name}')
                self.pymol.cmd.delete(fore.name)
                self.pymol.cmd.delete(aft.name)
                self.old2new[fore.name] = self.models[i]
                self.old2new[aft.name] = self.models[i + 1]
                self.models[i + 1].name = new_name # aft becomes the combo.
                #self.models[i + 1].tier #the Cterminus will be the tier of the previous.
                self.models[i] = None  # kill fore from list. # fore may
        self.models = [m for m in self.models if m is not None]
        return self

    def project_models(self):
        """
        places models with gaps along the x-axis.
        """
        for m in self.models:  # no idea why this is needed.
            m.angle_fix()
        for i in range(0, len(self.models) - 1):
            try:
                fore, aft = self.get_model_pair(i)
            except ValueError:
                continue
            self.project(fore, aft)
            self.roll_free(fore, aft)
            self.log_joint('projection', fore, aft)
        return self

    def show_pose(self, filename='temp.png'):
        """
        Jupyter notebook specific. Change to NGLview!
        BAD NAME. Not pyrosetta object. Will change.
        """
        if '.png' not in filename:
            filename += '.png'
        if os.path.exists(filename):
            os.remove(filename)
        self.pymol.cmd.png(filename)
        while not os.path.exists(filename):
            time.sleep(1)
        display(Image(filename=filename))

    def show_before(self):
        self.add_xyz()
        modsele = f'{self.models[0].make_selection(self.models[0].pdb_end)} or {self.models[1].make_selection(self.models[1].pdb_start)} or {self.models[1].make_selection(self.models[1].pdb_end)}'
        self.pymol.cmd.color('yellow', modsele)
        self.pymol.cmd.show(representation="sticks", selection=modsele)
        self.pymol.cmd.zoom(modsele)
        print('Before')
        self.show_pose('before.png')

    def show_after(self):
        print('after')
        self.pymol.cmd.zoom('all')
        self.show_pose('after.png')
        self.pymol.cmd.save('after.pdb', 'all')

    def show_load_checks(self):
        self.pymol.cmd.zoom('all')
        for model in self.models:
            self.pymol.cmd.color('cyan', model.name)
            self.pymol.cmd.color('red', f'{model.name} and chain {model.chain}')
            modsele = model.make_selection(f'{model.pdb_start}-{model.pdb_end}')
            self.pymol.cmd.color('white', modsele)
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
        """
        * Aligns them
        * Iterate across the overlap pair_fitting the backbone and getting a score of the overlap
        * get best overlap.
        * on bad overlaps rotate the aft by 60deg. Not sure why, this was done.
        """
        mobile = aft.make_selection(f"{aft.pdb_start}-{fore.pdb_end}")
        m_match = self.safe_select(mobile)
        target = fore.make_selection(f"{aft.pdb_start}-{fore.pdb_end}")
        t_match = self.safe_select(target)
        if m_match == 0 or t_match == 0:
            warn(f'Possible major issue...\n{mobile} ... {m_match}\n{target} ... {t_match}')
            return self
        self.pymol.cmd.align(mobile, target)
        scores = {}
        # Align all combinations by mapping the backbone and find the one with least clashes.
        if fore.tier <= aft.tier:
            ranger = range(fore.pdb_end, aft.pdb_start - 1, -1)
        else:
            ranger = range(aft.pdb_start, fore.pdb_end + 1, +1)
        for over_i in ranger:
            for model in (aft, fore):
                n_ca = self.pymol.cmd.select(aft.make_selection(over_i, 'CA'))
                assert n_ca == 1, f'There are {n_ca} atoms for {model}'
            self.pymol.cmd.pair_fit(aft.make_selection(over_i, 'N'), fore.make_selection(over_i, 'N'),
                                    aft.make_selection(over_i, 'CA'), fore.make_selection(over_i, 'CA'),
                                    aft.make_selection(over_i, 'C'), fore.make_selection(over_i, 'C'))
            if self.debug:
                print(f'... {mobile} --> {target}')
                self.pymol.cmd.zoom(mobile)
                self.pymol.cmd.color('red', mobile)
                self.pymol.cmd.color('green', target)
                self.show_pose('align.png')
            if cycle > 0:
                aft.roll(60 * cycle)
            scores[over_i] = self.pymol.cmd.overlap(fore.make_selection(f"{fore.pdb_start}-{over_i}"), aft.make_selection(f"{over_i + 1}-{aft.pdb_end}"))
            if scores[over_i] == 0: # no overlap.
                break
            else:
                if self.debug:
                    print(f'!!! There is an overlap of {scores[over_i]}!')
        else:
            over_i = list(scores.keys())[list(scores.values()).index(min(scores.values()))]
            if scores[over_i] < 100:
                if self.debug:
                    print(f'Aligning via residue {over_i}. score={scores[over_i]}')
                self.pymol.cmd.pair_fit(fore.make_selection(over_i, 'N'), aft.make_selection(over_i, 'N'),
                                        fore.make_selection(over_i, 'CA'), aft.make_selection(over_i, 'CA'),
                                        fore.make_selection(over_i, 'C'), aft.make_selection(over_i, 'C'))
                if cycle > 0:
                    aft.roll(60 * cycle)
            else:
                print('Good alignment impossible.')
                if cycle < 6:
                    self.align(fore, aft, cycle=cycle + 1)
                else:
                    raise Exception('Really?? Overlap after alignment even with rotations?? Sigh.')
        if fore.tier <= aft.tier:
            if self.debug:
                print(f'overlap removal')
                print('fore', fore)
                print('aft', aft)
            self.pymol.cmd.remove(fore.make_selection(f'{over_i + 1}-{fore.pdb_end}'))
            fore.pdb_end = over_i
            self.pymol.cmd.remove(aft.make_selection(f'{aft.pdb_start}-{over_i}'))
            aft.pdb_start = over_i + 1
            assert aft.make_selection(aft.pdb_start) != 0, f'Model start is in gap????  {aft.pdb_start}'
        else:
            if self.debug:
                print(f'overlap: {over_i}')
                print('fore', fore)
                print('aft', aft)
            self.pymol.cmd.remove(fore.make_selection(f'{over_i}-{fore.pdb_end}'))
            fore.pdb_end = over_i - 1
            self.pymol.cmd.remove(aft.make_selection(f'{aft.pdb_start}-{over_i - 1}'))
            aft.pdb_start = over_i
            while aft.make_selection(aft.pdb_start) == 0:
                aft.pdb_start += 1
                warn(f'Model start is in gap????  {aft.pdb_start}')
        return self

    def fetch_metadata(self):
        if self.debug:
            print('# Fetching metadata')
        # fetch the metadata
        for entry in self.pdbs:
            mode, code, chain = Model.split_id(entry)
            if mode.upper() == 'PDB':
                reply = requests.get(f'http://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{code}').json()
                self.metadata[code] = reply[code.lower()]
            else:
                self.metadata[code] = {}
        return self

    def categorise(self):
        """
        Clusters the proteins
        """
        if self.debug: print('# Running categorisation')
        # find longest!
        tdata = []
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
                e = entry['uniprot_end'] - entry['uniprot_start'] + 1
                seq = None
            tier = float(entry['tier'])
            tdata.append({'code': code,
                          'mode': mode,
                          'chain': chain,
                          'tier': tier,
                          'uniprot_start': entry['uniprot_start'],
                          'uniprot_end': entry['uniprot_end'],
                          'seq': seq,
                          'url': entry['url'],
                          'pdb_start': s,
                          'pdb_end': e,
                          'length': e - s})
        df = pd.DataFrame.from_records(tdata).sort_values(['uniprot_start', 'length', 'tier'],
                                                          ascending=[True, False, True])
        if self.debug:
            print('There are all the chains:')
            display(df)
        # thin!
        new = df.iloc[0]
        get_similars = lambda: list(
            set(df.loc[(df.uniprot_end <= new.uniprot_end) & (df.uniprot_start >= new.uniprot_start)].code.values) - {
                new.code})
        new.at['similars'] = get_similars()
        o = df.loc[df.uniprot_end > new.uniprot_end]
        parts = [new]
        old = new
        while len(o):
            new = o.iloc[0]
            new.at['similars'] = get_similars()
            if new.uniprot_start < old.uniprot_end:  # mark overlap for debug.
                new.at['subtender'] = old.code
            # curious case of better models within makes this weird.
            # so something is better is pulled out of the similars (i.e. smaller so excluded.) group.
            # but to be better it has to be of a lower tier and not present in the selection.
            betters = df[(df.code.isin(new.similars)) & (df.tier < new.tier) & (~df.code.isin([p.code for p in parts]))]
            bi = len(betters.index)
            if bi == 0:
                parts.append(new)  # has to happen!!
            else:
                better = betters.iloc[0]
                print(f'Nested case {new.code} and {betters.code.values}')
                # dangerous hack. It assumes that a good match can be found and won't overrun.
                # as the priority differs it will chew off one side then the other, hence the double.
                if better.uniprot_start > new.uniprot_start:
                    parts.append(new)
                parts.append(better)
                if better.uniprot_end < new.uniprot_end:
                    parts.append(new)
            # cleanup.
            o = o.loc[o.uniprot_end > new.uniprot_end]
            old = new
        if self.debug:
            print('These are the non-redudant chains:')
            display(pd.DataFrame(parts))
        self.models = [Model(fuser=self, metadata=self.metadata[row.code], **row) for row in parts]
        return self

    def ready_models(self):
        """
        Runs model.fetch_n_clean() on all models.
        """
        if self.debug: print('# Reading models')
        for model in self.models:
            model.fetch_n_clean()
        self.models[0].camera_fix()
        return self

    def project(self, fore, aft, extra_d=0):
        """extra_d is to move it away when there is on clean rotational solution."""
        aft.angle_fix()
        (ref, pointer, θ, φ, internal_d) = fore.angles
        (ref_aft, pointer_aft, θ_aft, φ_aft, internal_d_aft) = aft.angles
        # translate
        current_end = aft.get_resi_vector(aft.pdb_start)['N']  # re-get as it moved.
        missing = aft.pdb_end - fore.pdb_start
        missing_d = abs(missing * 3.5) ** 0.5  # the James Ross fudge distance
        if self.debug: print(f'Projecting {missing_d} + {extra_d} Å')
        d = missing_d + internal_d + extra_d
        target_end = np.array([d * np.cos(-θ) * np.cos(φ) + ref.x,
                               d * np.sin(-θ) * np.cos(φ) + ref.y,
                               d * np.sin(-θ) * np.sin(φ) + ref.z])
        self.pymol.cmd.translate((target_end - current_end).tolist(), aft.name, camera=0)
        if self.debug:
            fore.add_NC_arrow('cylinderF', r=0, g=0, b=1)
            aft.add_NC_arrow('cylinderA', r=0, g=1, b=0)

    def has_overlap(self, fore, aft):
        o = self.pymol.cmd.overlap(f"model {fore.name} and chain {fore.chain}",
                                   f"model {aft.name} and chain {aft.chain}")
        if o > 0:
            if self.debug:
                print(f"Overlap between {fore.name} chain {fore.chain} and {aft.name} chain {aft.chain} of {o}")
            return True
        o = self.pymol.cmd.overlap(f"model {fore.name}", f"model {aft.name}")
        if o > 0:
            if self.debug:
                print(f"Overlap between {fore.name} (whole) and {aft.name} (whole)")
            return True
        return False

    def add_xyz(self):
        self.pymol.cmd.load_cgo([9.0, 0, 0, 0, 100, 0, 0, 3, 1, 1, 1, 0, 0, 0], "cylinderx")
        self.pymol.cmd.load_cgo([9.0, 0, 0, 0, 0, 100, 0, 3, 1, 1, 1, 0, 0, 0], "cylindery")
        self.pymol.cmd.load_cgo([9.0, 0, 0, 0, 0, 0, 100, 3, 1, 1, 1, 0, 0, 0], "cylinderz")

    def safe_select(self, name:str, sele:str=None, silent=False):
        """
        There is something not quite right with the selection. It throws an error..."""
        try:
            if not sele:
                sele = name
                name = 'test'
            n = self.pymol.cmd.select(name, sele) # name is the selection name (sele)
            return n
        except:
            if '{' in sele:
                raise ValueError(f'There seems to be a curly bracket in the selection ({sele}): bad formatting?')
            else:
                if not silent:
                    warn(f'Failed selection "{sele}"')
                return 0

    # ================= Uniprot ========================================================================================
    @classmethod
    def from_uniprot(cls, uniprot: str, debug: bool = False):
        data = cls.get_converted_data(uniprot)
        self = cls(data, debug)
        return self

    @classmethod
    def get_converted_data(cls, uniprot: str):
        return [cls.convert_SM_entry(p) for p in cls.get_data(uniprot, provider='pdb')] + \
               [cls.convert_SM_entry(p) for p in cls.get_data(uniprot, provider='swissmodel')]

    @staticmethod
    def convert_SM_entry(data):
        """
        From swissmodel api to this jumble of a format...
        """
        assert len(data['chains']) == 1, f'There are multiple chains! {data["chains"]}'
        try:
            convert = dict(
                chain=data['chains'][0]['id'],
                url=data['coordinates'],
                uniprot_start=data['chains'][0]['segments'][0]['uniprot']['from'],
                uniprot_end=data['chains'][0]['segments'][0]['uniprot']['to']
            )
            if data['provider'] == 'PDB':
                convert['tier'] = 1
                convert['id'] = data['template'] + '_' + data['chains'][0]['id']
                convert['description'] = data['chains'][0]['segments'][0]['pdb']['description']
                convert['pdb_start'] = data['chains'][0]['segments'][0]['pdb']['from']
                convert['pdb_end'] = data['chains'][0]['segments'][0]['pdb']['to']
            else:
                convert['tier'] = 2
                span = str(data['chains'][0]['segments'][0]['uniprot']['from']) + \
                        str(data['chains'][0]['segments'][0]['uniprot']['to'])
                convert['id'] = data['provider'] + ':' + data['template'] + span
                descr = data['chains'][0]['segments'][0]['smtl']['description']
                convert['description'] = f"{data['template']} ({descr}, {data['similarity']:.1}%)"
                convert['pdb_start'] = data['chains'][0]['segments'][0]['uniprot']['from']
                convert['pdb_end'] = data['chains'][0]['segments'][0]['uniprot']['to']
        except KeyError as error:
            raise KeyError(f'The key {error} could not be found in {data}')
        return convert

    @staticmethod
    def get_data(uniprot: str, provider='swissmodel') -> List[dict]:
        """
        Gets data from Swissmodel. Convert each entry with convert_SM_entry
        """
        assert provider in ('swissmodel', 'pdb'), 'Provider has to be swissmodel or pdb'
        url = f'https://swissmodel.expasy.org/repository/uniprot/{uniprot.strip()}.json?provider={provider.strip()}'
        data = requests.get(url).json()
        return data['result']['structures']  # [0] # from to coordinates description


# ======================================================================================================================

class Model:
    pymol = None

    def __init__(self, mode, chain, code, uniprot_start, uniprot_end, pdb_start=1, pdb_end=None, seq=None, length=None,
                 fuser=None, tier=0,
                 metadata=None, similars=None, subtender=None, url=None):
        """
        Init just stores the variables. but .fetch_and_clean() does all the work by calling:
        * model.kill_states()     ---> use only state 1.
        * model.move_away_hetatm()  ---> protect the hetatms ligands!
        * model.clean_Nterminus()   ---> trim the linkers!
        * model.clean_Cterminus()
        * model.fix_numbers()   --> the numbering off quite frequently...
        * model.shift_to_true_Nterminus()  ---> gets rid of unsolved terminus
        * model.shift_to_true_Cterminus()  ---> gets rid of unsolved terminus
        While other methods are called by the project method of fuser. Say all the angle operations.
        e.g.model.roll(n) or model.get_resi_coords(resi)
        
        The code can be a PDB code. or... PDB:code  or SWISSMODEL:id or LOCAL:filename (without underscores!)
        
        There are three sets of numbers.
        uniprot start (`uniprot_start`) which is the position in the whole protein of the first real residue based on the stated sequence regardless of density.
        pdb_start (`pdb_start`) which is the first residue stated in PDBe query based on the stated sequence regardless of density.
        Namely... the mapping start is at what residue index within the model does the match start,
        while uniprot start is what the position that residues has in the whole sequence.
        A method to keep an eye for is `.angle_fix()`. This ensures that the the N and C termini are on the x axis.
        """
        self.chain = chain
        self.mode = mode
        self.code = code
        self.name = code.replace('.','')[:8]  # a name for pymol. and backup value in case of duplicates...
        self.tier = tier
        self.metadata = metadata
        self.uniprot_start = uniprot_start  # uniprot
        self.uniprot_end = uniprot_end
        self.pdb_start = pdb_start  # pdb numbering
        if pdb_end:
            self.pdb_end = pdb_end
        else:
            self.pdb_end = uniprot_end - uniprot_start
        if length:
            self.length = length
        else:
            self.length = uniprot_end - uniprot_start
        self.url = url
        self.seq = seq
        self.subtender = subtender
        if fuser:
            self.fuser = fuser
            self.pymol = fuser.pymol
        elif self.pymol is None:
            raise EnvironmentError('Please assign a pymol session to `cls.pymol`.')
        else:
            warn('No fuser?')
            self.fuser = namedtuple('faux', ['debug'])(False)
        self.length = length

    @classmethod
    def from_ds(cls, ds):
        cls(**ds)

    def make_selection(self, resi:Optional[int]=None, atom:Optional[str]=None, no_hetatm:bool=True):
        chain = f'chain {self.chain}'
        if no_hetatm:
            extra = 'not hetatm'
        else:
            extra = ''
        if resi is None:
            resi = ''
        elif isinstance(resi, float):
            resi = f'resi {int(resi)}'
        else:
            resi = f'resi {resi}'
        if atom is None:
            atom = ''
        else:
            atom = f'name {atom}'
        return f'({" and ".join([x for x in (self.name, chain, resi, atom, extra) if x])})'

    def __str__(self):
        detail = f'name={self.name} chain={self.chain} start={self.pdb_start} end={self.pdb_end} tier={self.tier}'
        return f'<Model({detail}) object at 0x{id(self):02x}>'

    def describe(self):
        (ref, pointer, θ, φ, internal_d) = self.angles
        print(f'Angles: θ={np.rad2deg(θ)} and φ={np.rad2deg(φ)})')
        print('Nterm:', self.get_resi_coords(self.pdb_start)['N'])
        print('Cterm:', self.get_resi_coords(self.pdb_end)['C'])

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
        mode, code, chain = cls.split_id(ID)
        m = cls(mode=mode, chain=chain, code=code, pdb_start=1, pdb_end=1, uniprot_start=1, uniprot_end=1, seq='',
                length=1).load()
        residex = m.resi_list()
        cls.pymol.cmd.delete(m.name)
        return {'pdb_start': residex[0], 'pdb_end': residex[-1], 'length': residex[-1] - residex[0]}

    def resi_list(self):
        return sorted(set([int(at.resi) for at in self.pymol.cmd.get_model(self.name).atom]))

    def residues(self):
        rset = set([(at.resn, int(at.resi), True) for at in self.pymol.cmd.get_model(self.name).atom])
        r = sorted(rset, key=lambda x: x[1])
        s = ''.join(self.pymol.cmd.get_fastastr(f'{self.name} and chain {self.chain} and not hetatm').split('\n')[1:])
        first_stated_resn = seq3(s[0]).upper()
        first_structured_resi = r[0][1]
        first_structured_resn = r[0][0]

        if first_structured_resi == 1:
            if first_stated_resn == first_structured_resn:  # case 1. all is good.
                pass
            else:  # case 2. negative numbers??!
                raise NotImplementedError(f'First residues are {r[0:3]} but the sequence is {s}')
        elif first_structured_resi > 1:
            if first_stated_resn == first_structured_resn:  # case 3. The first residue exist is not 1.
                pass
            elif seq3(s[first_structured_resi - 1]).upper() == first_structured_resn:
                # case 4. the first stated residue is 1, but does not exist.
                r = [(seq3(s[i]).upper(), i + 1, False) for i in range(1, first_structured_resi)] + r
            else:
                raise NotImplementedError(f'First residues are {r[0:3]} but the sequence is {s}')
        else:
            raise NotImplementedError(f'First residues are {r[0:3]} but the sequence is {s}')

        return r

        # ================ METHODS FOR LOADING =========================================================================

    def fetch_n_clean(self):
        if self.fuser.debug: print(f'# Fetching for {self.name}')
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
        self.name = self.name[:8]
        if self.name in self.pymol.cmd.get_object_list():
            warn(f'{self.name} is already loaded. Copying as {self.name}X.')
            self.pymol.cmd.create(self.name + 'X', self.name)
            self.name += 'X'
            assert self.pymol.cmd.select(self.name), f'{self.name} has no atoms.'
            return self
        if self.mode == 'PDB':
            print('Fetching biological assembly.')
            self.name = self.pymol.cmd.fetch(self.code, type='pdb1')
            if self.name == -1:
                self.name = self.pymol.cmd.fetch(self.code)
                print('\t... solved by fetching assimetric unit.')
            if self.name == -1:
                raise ValueError(f'Could not load {self.code}')
            self.kill_states()
        elif self.mode == 'SWISSMODEL':
            if self.url is None:
                self.url = f'https://swissmodel.expasy.org/repository/{self.code}.pdb'
            if self.fuser.debug:
                print(self.url)
            pdbblock = requests.get(self.url).text
            if 'END' in pdbblock:
                self.pymol.cmd.read_pdbstr(pdbblock, self.name)
            else:
                print(pdbblock)
                raise ValueError(f'{self.url} failed: not a PDB')
            chains = sorted({at.chain for at in self.pymol.cmd.get_model(self.name).atom})
            if len(chains) == 1:
                if chains[0] == '':
                    self.pymol.alter.cmd(self.name, 'chain="A"')
                    self.chain = 'A'
                else:
                    self.chain = chains[0]
            else:
                raise NotImplementedError
        elif self.mode == 'LOCAL':
            self.name = os.path.splitext(self.code)[0][:8]
            self.pymol.cmd.load(self.code, self.name)
        else:
            raise ValueError(f'UNKNOWN mode {self.mode}. PDB: SWISSMODEL: LOCAL:')
        if self.seq is None:
            self.seq = ''.join(
                self.pymol.cmd.get_fastastr(f'{self.name} and chain {self.chain} and not hetatm').split('\n')[1:])
        return self

    def force_chain_A(self):
        get_untaken = lambda c: list(set('ABCDEFGHIJKLMNOPQRSTUVWXYZ').difference(c))[0]
        chains = sorted({at.chain for at in self.pymol.cmd.get_model(self.name).atom})
        if self.chain != 'A':
            if 'A' in chains:
                # move chain A away.
                new = get_untaken(chains)
                self.pymol.cmd.alter(f'model {self.name} and chain A', f'chain="{new}"')
                self.pymol.cmd.sort()
            self.pymol.cmd.alter(f'model {self.name} and chain {self.chain}', f'chain="A"')
            self.pymol.cmd.sort()
        else:
            pass
        segis = sorted({at.chain for at in self.pymol.cmd.get_model(f'model {self.name} and chain A').atom})
        if len(segis) > 1:
            if '' not in segis:
                subs = segis[1:]
            else:
                i = segis.index('')
                subs = segis[:i] + segis[i + 1:]
            for s in subs:
                new = get_untaken(chains)
                self.pymol.cmd.alter(f'model {self.name} and chain A and segi {s}', f'chain="{new}"')
                self.pymol.cmd.sort()
        elif segis[0] != '':
            self.pymol.cmd.alter(f'model {self.name} and chain A and segi {segis[0]}', f'segi=""')
            self.pymol.cmd.sort()
        else:
            pass
        self.chain = 'A'
        return self

    def kill_states(self):
        # fetch(.., state=1) does not work. !?
        self.pymol.cmd.split_states(self.name, 1, 1)
        self.pymol.cmd.set_name(f'{self.name}_0001', self.name)
        return self

    def move_away_hetatm(self):
        if self.fuser.debug:
            print(f'# Protecting hetatms for {self.name}')
        """
        Renumbers the hetatm residues to 9000+ values, unique in the whole fuser.
        """
        hetatms = {self.make_selection(at.resi, no_hetatm=False) for at in self.pymol.cmd.get_model(self.name).atom if
                   at.hetatm}
        for h in hetatms:
            self.pymol.cmd.alter(h, f"resi={self.fuser.last_het_index}")
            self.fuser.last_het_index += 1
        self.pymol.cmd.sort()
        return self

    def start_from_one(self):
        """
        The first residue may exist or not.
        """
        residues = self.residues()
        if residues[0][1] != 1:
            Δ = residues[0][1] - 1
            self.pymol.cmd.alter(f'{self.name} and chain {self.chain}', f'resv-={Δ}')
            self.pymol.cmd.sort()
        return self

    def clean_Nterminus(self):
        sele = self.make_selection(f'0-{self.pdb_start - 1}')
        if self.pymol.cmd.select(sele) == 0:
            return self
        if self.fuser.debug:
            print(f'# NTerm cleaning for {self.name}')
        self.pymol.cmd.remove(sele)
        self.pymol.cmd.sort()
        return self

    def clean_Cterminus(self):
        if self.fuser.debug:
            print(f'# CTerm cleaning for {self.name}')
        """
        convert the N of the junk resi to OXT.
        """
        has_trailing = self.pymol.cmd.alter(self.make_selection(self.pdb_end + 1, 'N'), 'name="OXT"')
        if has_trailing:
            self.pymol.cmd.alter(self.make_selection(self.pdb_end + 1, 'N'), f'resi="{self.pdb_end}"')
            self.pymol.cmd.sort()
            resn = self.pymol.cmd.get_model(self.make_selection(self.pdb_end)).atom[0].resn
            self.pymol.cmd.alter(self.make_selection(self.pdb_end + 1, 'OXT'), f'resn="{resn}"')
            self.pymol.cmd.sort()
        self.pymol.cmd.remove(self.make_selection(f'{self.pdb_end + 1}-9999'))
        self.pymol.cmd.sort()
        return self

    def fix_numbers(self):

        if self.pdb_start != self.uniprot_start:
            Δ = self.uniprot_start - self.pdb_start
            if self.fuser.debug:
                print(f'# Fixing numbering for {self.name} (offset={Δ}')
            self.pymol.cmd.alter(f'{self.name} and chain {self.chain}', f'resv+={Δ}')
            self.pymol.cmd.sort()
            self.pdb_start = self.pdb_start + Δ
            self.pdb_end = self.pdb_end + Δ
        return self

    def shift_to_true_Nterminus(self):
        """
        The model may have missing residues at the NTerm, but in reality it is as if they were not there.
        """
        if self.fuser.debug: print(f'# Shifting Nterm for {self.name}')
        # verify that the start exists.
        coords = None
        current_start = self.pdb_start
        while coords is None:
            coords = self.pymol.cmd.get_coords(self.make_selection(current_start))
            current_start += 1
        self.pdb_start = current_start - 1

    def shift_to_true_Cterminus(self):
        """
        The model may have missing residues at the CTerm, but in reality it is as if they were not there.
        """
        if self.fuser.debug: print(f'# Shifting Cterm for {self.name}')
        # verify that the end exists.
        coords = None
        l = self.pdb_end
        while coords is None:
            coords = self.pymol.cmd.get_coords(self.make_selection(l))
            l -= 1
            if l < 0:
                raise ValueError(
                    f'Trimming missing residues went haywire! start={self.uniprot_start}/{self.pdb_start} end={self.uniprot_end}/{self.pdb_end}. Died at {self.make_selection(l)}')
        self.pdb_end = l + 1

    # ============ METHODS FOR ADJUSTING ===============================================================================
    def get_resi_coords(self, resi):
        match_coords = lambda c: {'x': c[0], 'y': c[1], 'z': c[2]}
        coordsdex = self.get_resi_vector(resi)
        return {atom: Coords(name=atom, **match_coords(coordsdex[atom]))
                for atom in coordsdex}

    def distance(self, A, B):
        return ((A.x - B.x) ** 2 + (A.y - B.y) ** 2 + (A.z - B.z) ** 2) ** 0.5

    def get_resi_vector(self, resi):
        # if an error occurs here, it is because there is no self.make_selection(resi), which should not be possible
        # as it is dealt with upstream in the trimming. Export the scene with pymol.exporting.multisave('tmp.pdb')
        # as see what is wrong. It is somethign weird. Maybe a hetatm not marked as such?
        return {atom: self.pymol.cmd.get_coords(self.make_selection(resi, atom))[0] for atom in ('C', 'N', 'CA')}

    @property
    def angles(self):
        ref = self.get_resi_coords(self.pdb_start)['N']
        pointer = self.get_resi_coords(self.pdb_end)['C']
        internal_d = self.distance(pointer, ref)
        θ = np.arctan((pointer.y - ref.y) / (pointer.x - ref.x))
        φ = np.arctan((pointer.z - ref.z) / (pointer.x - ref.x))
        return ref, pointer, θ, φ, internal_d

    def add_NC_arrow(self, name, r=0, g=0, b=0):
        N = self.get_resi_coords(self.pdb_start)['N']
        C = self.get_resi_coords(self.pdb_end)['C']
        self.pymol.cmd.load_cgo([9.0, N.x, N.y, N.z, C.x, C.y, C.z, 3, r, g, b, r / 3, g / 3, b / 3], name)

    def camera_fix(self):
        # self.pymol.cmd.alter('all','segi=""')
        self.angle_fix()

    def angle_fix(self, cycles=3):
        if self.fuser.debug:
            print(f'Initial position for model {self.name}.....')
            self.describe()
        # rotate the xy plane
        for i in range(cycles):
            (ref, pointer, θ, φ, internal_d) = self.angles
            self.pymol.cmd.translate([-ref.x, -ref.y, -ref.z], self.name, camera=0)
            self.pymol.cmd.rotate([0, 0, 1], -np.rad2deg(θ), selection=self.name, camera=0)
            # rotate the xz plane
            (ref, pointer, θ, φ, internal_d) = self.angles
            self.pymol.cmd.translate([-ref.x, -ref.y, -ref.z], self.name, camera=0)
            self.pymol.cmd.rotate([0, 1, 0], np.rad2deg(φ), selection=self.name, camera=0)
            (ref, pointer, θ, φ, internal_d) = self.angles
            self.pymol.cmd.translate([-ref.x, -ref.y, -ref.z], self.name, camera=0)
            tail = self.get_resi_coords(self.pdb_end)['C']
            if tail.x < 0:
                if self.fuser.debug:
                    print('It is off by 180... Odd.')
                self.pymol.cmd.rotate([0, 0, 1], 180, selection=self.name, camera=0)  # rotate the xy plane
                (ref, pointer, θ, φ, internal_d) = self.angles
                self.pymol.cmd.translate([-ref.x, -ref.y, -ref.z], self.name, camera=0)
        if self.fuser.debug:
            print(f'Final position for model {self.name}.....')
            self.describe()
        return self

    def roll(self, φ):
        """
        Whereas they are aligned on the x-axis with the y,z coordinates of C and N termini near 0,0.
        They can rotate around the X. roll. This is good to remove clashes.
        """
        self.pymol.cmd.rotate([1, 0, 0], φ, selection=self.name, camera=0)
        self.angle_fix()
        return self

if __name__ == '__main__':
    import pymol2, argparse

    parser = argparse.ArgumentParser(description='Combine structures for a Uniprot')
    parser.add_argument('uniprot', type=str, help='Uniprot to generate a combined model')
    args = parser.parse_args()

    with pymol2.SingletonPyMOL() as pymol:
        Fuser.pymol = pymol
        w = Fuser.from_uniprot(args.uniprot, debug=False)
        w.order()
        w.save()
        #w.add_xyz()