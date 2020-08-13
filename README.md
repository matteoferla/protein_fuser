# protein_fuser
A py3 script to fuse structures together.

![fig](protein_fixer-01.png)

This script, best run in a Jupyter notebook, requires pymol within your py environment, which can be installed with conda.
This is not a PyMOL script, but a Python3 script that uses `module pymol`.

## Raison d'etre

Human protein, in particular, are often huge multidomain protein, like beads on a string. Each domain crystallised separately.
Although it is important to use the primary source, it is sometimes beneficial to stitch domains together for a **simple** illustration without the confusion of separate images of say model A with domain X solved at 2.1 &Aring; while model B is an NMR of domain Y bound to protein 2, _etc._

If Rosetta remodel or the ITasser DEMO tool is used, the protein ends up looking like a tangle and not like beads...

A nice analogy are the planets. When pictured they are in syzygy, 
which is actually extremely rare or impossible depending on how precisely one expects them to align.

This script does the following:

* Finds which bits to use based on length, start, stop and tier &mdash; resolution and missing loops are not taken into account.
* Aligns models that can be aligned due to overlap, with shifts and, if need be, rotation, to find the best two overlapping amino acids to stitch the protein together.
* Places the N and C termini of each model-agglomerate on the x axis `(x, 0, 0)`.
* If there is a gap it projects the N terminus of the next model by missingAA &times; 3.5^0.5 &Aring; away from the C terminus of the previous. If there is a clash, the model is rotate/pushed out until it is okay.

The resulting protein will be chain A. With any other chains being from the biologically assembly if present.
At present it does not play with the partners, although in future it would be nice if it gathered binding partners too. The metadata comes from the PDBe API and albeit the gene symbols are a bit of a mess they are definitely matchable.

The attribute `ws.joints` contains information of what came from where.

**CAVEAT** It is essential that you let any viewer know that your protein is a Frankenstein protein and it is not known how the bits fit together.

## Next

* I have code to make Rosetta add the missing loops that may be nice to add.
* I have code to make Rosetta add the post translation modifications that may be nice to add.
* I would like to automated choose the binding partners
* I would like to run the script on the whole of the human PDB, with SWISSMODEL, Phyre2 (they have a nice amount thanks to Miscast!) and I-TASSER*

&lowast; I have written a [scraper for I-Tasser](https://github.com/matteoferla/ITasser_miner), which I think is okay, but I need to find out.

## Input

The input is based on my Uniprot parser (see [protein-module-for-VENUS](https://github.com/matteoferla/protein-module-for-VENUS)).
This is slightly odd in that it is for use with the NextProt featureViewer ([ref](https://github.com/calipho-sib/feature-viewer)), so the pdb is a list of dictionaries with keys `id`, `x`, `y`, `description` &mdash;description is not used, while id value is `(?PDB|SWISSMODEL|LOCAL:)<code>_<chain>` format, _e.g._ `1WG7_A` or `PDB:1WG7_A` or `1WG7` or `SWISSMODEL:5bf8103b02efd001eeb29f09` or `LOCAL:myfile.pdb`. Additionally, the field `tier` can be added to specify priority. _E.g._ a Phyre2 model is less reliable than a PDB model, so the PDB model should be trusted more, so the user would specify a lower tier (e.g. `1` vs. `2`).


## Workshop class

	ws = Workshop(pdbs, debug).order().save()

* `debug`: verbosity boolean
* `pdbs`: `[{'x': stated_start_in_uniprot,
    'y': stated_end_in_uniprot,
    'id': 'code_chain', 
    'tier': opt. int. low => max importance (e.g. crystal), high => lower importance (e.g. models)
    'description': opt. 'Not used.'}]`

In the case of PDBs the metadata is taken in order to find where the true start is without cloning scars.
The code can be a PDB code. or `PDB:code`  or `SWISSMODEL:id` or `LOCAL:filename` (without underscores!)

For PDB codes, the metadata is fetched to verify what the true values are. While for swissmodel and local it is assumed that the first residue that exists in the model is the `x` in the input.
Note that Uniprot and PDBe metadata `x` refers to the first residue in the cif file sequence, regardless of existence in the solved structure as is often not the case.



`ws.order()` readies and aligns and projects the non-redudant models. The all in one workhorse.
the init step actually calls `wa.categorise()` which categorises the models based on redundancy.
`ws.ready_models()` loads models
`ws.align_models()` aligns overlapping ones
`ws.project_models()` places models with gaps along the x-axis.
`ws.models` list of models.
For various reasons, `model.code` is what is used to fetch/load/pdbstr, while `model.name` is the name of the model within PyMol.
See model for more about the madness.

Note that there is not a `pymol.cmd.create` step to fuse the parts as simply saving as pdb will do it and instead saving as pse to check things is nice or if there is a version mismathc with pymol (_i.e._ system = Pymol1.8, python pymol2) `pymol.exporting.multisave('tmp.pdb')` is a handy command.

Model has many bound methods that use pymol cmds, that affect only one model, say `protein.roll(40)`. while workshop uses more thatn one say, `ws.roll_free(protein_A, protein_B)`.

Some methods act on all models and will have `_models` in the name.

* `ws.ready_models()` --> `model.fetch_n_clean()`
* `ws.align_models()` --> `ws.align(A,B)`
* `ws.project_models()` --> `ws.project(A,B)`

`ws.show_pose()` is a Jupyter notebook specific thing of mine as I am running in headless mode. Namely, it calls `pymol.cmd.save`, waits for the async function to do its thing and load the image as a `display(Image(filename))`.

## Model

 Init just stores the variables. but `model.fetch_and_clean()` does all the work loading the structure by calling:

* `model.kill_states()`    ---> use only state 1.
* `model.move_away_hetatm()`  ---> protect the hetatms ligands!
* `model.clean_Nterminus()`   ---> trim the linkers!
* `model.clean_Cterminus()`
* `model.fix_numbers()`   --> the numbering off quite frequently...
* `model.shift_to_true_Nterminus()`  ---> gets rid of unsolved terminus
* `model.shift_to_true_Cterminus()`  ---> gets rid of unsolved terminus

While other methods are called by the project method of workshop. Say all the angle operations.
e.g. `model.roll(n)` or `model.get_resi_coords(resi)`. 
A method to keep an eye for is `model.angle_fix()`. This ensures that the the N and C termini are on the x axis.
        
The code can be a PDB code. or... PDB:code  or SWISSMODEL:id or LOCAL:filename (without underscores!)

There are three sets of numbers.
uniprot start (`x`) which is the position in the whole protein of the first real residue based on the stated sequence regardless of density.
start (`start`) which is the first residue stated in PDBe query based on the stated sequence regardless of density.
Namely... the mapping start is at what residue index within the model does the match start,
while uniprot start is what the position that residues has in the whole sequence.

## Example
	
	uniprot = 'Q00341' ## vigilin
	pdbs = [{'x': 432, 'y': 502, 'id': '1VIG_A', 'description': '1VIG'}, 
		{'x': 432, 'y': 502, 'id': '1VIH_A', 'description': '1VIH'}, 
		{'x': 142, 'y': 222, 'id': '2CTE_A', 'description': '2CTE'},
		{'x': 346, 'y': 434, 'id': '2CTF_A', 'description': '2CTF'}, 
		{'x': 645, 'y': 726, 'id': '2CTJ_A', 'description': '2CTJ'},
		{'x': 964, 'y': 1054, 'id': '2CTK_A', 'description': '2CTK'},
		{'x': 1044, 'y': 1127, 'id': '2CTL_A', 'description': '2CTL'},
		{'x': 1119, 'y': 1200, 'id': '2CTM_A', 'description': '2CTM'}]
	swissmodel = [{'x': 581, 'y': 724, 'id': '5d5620399ffd12cf2f734256', 'description': '2n8m.1.A (id:2e+01%)'},
		      {'x': 1053, 'y': 1198, 'id': '5c5f33ed8fd6f9f9ae62f475', 'description': '2jvz.1.A (id:2e+01%)'},
		      {'x': 583, 'y': 723, 'id': '5c5f33ed8fd6f9f9ae62f465', 'description': '3aev.1.B (id:3e+01%)'},
		      {'x': 726, 'y': 1068, 'id': '5d5620399ffd12cf2f73425e', 'description': '3n89.1.A (id:1e+01%)'},
		      {'x': 346, 'y': 434, 'id': '5c5f33ed8fd6f9f9ae62f469', 'description': '2ctf.1.A (id:1e+02%)'},
		      {'x': 1053, 'y': 1199, 'id': '5d56203a9ffd12cf2f734266', 'description': '3krm.1.A (id:2e+01%)'},
		      {'x': 508, 'y': 648, 'id': '5d5620399ffd12cf2f734246', 'description': '5wyj.31.A (id:2e+01%)'},
		      {'x': 584, 'y': 720, 'id': '5c5f33ed8fd6f9f9ae62f485', 'description': '1j4w.1.B (id:2e+01%)'},
		      {'x': 860, 'y': 1194, 'id': '5d5620399ffd12cf2f73424e', 'description': '3n89.1.A (id:2e+01%)'},
		      {'x': 655, 'y': 798, 'id': '5c5f33ed8fd6f9f9ae62f491', 'description': '2jvz.1.A (id:2e+01%)'},
		      {'x': 653, 'y': 796, 'id': '5c5f33ed8fd6f9f9ae62f499', 'description': '6fai.1.F (id:2e+01%)'},
		      {'x': 651, 'y': 799, 'id': '5d5620399ffd12cf2f73425a', 'description': '2n8m.1.A (id:2e+01%)'},
		      {'x': 73, 'y': 219, 'id': '5d5620399ffd12cf2f73423a', 'description': '3i82.1.A (id:1e+01%)'},
		      {'x': 1055, 'y': 1193, 'id': '5c5f33ed8fd6f9f9ae62f471', 'description': '1j4w.1.B (id:2e+01%)'},
		      {'x': 432, 'y': 502, 'id': '5be4a0a502efd0cd3370539f', 'description': '1vih (id:%)'},
		      {'x': 432, 'y': 502, 'id': '5be49e6e02efd0c70a281b5c', 'description': '1vig (id:%)'},
		      {'x': 346, 'y': 434, 'id': '5be49cbe02efd0c17498a975', 'description': '2ctf (id:%)'},
		      {'x': 1119, 'y': 1200, 'id': '5be4adbd02efd0f036d2a8c4', 'description': '2ctm (id:%)'},
		       {'x': 964, 'y': 1054, 'id': '5be4a89802efd0e17f75807a', 'description': '2ctk (id:%)'},
		      {'x': 1044, 'y': 1127, 'id': '5be4ab5f02efd0e942407c0c', 'description': '2ctl (id:%)'},
		      {'x': 142, 'y': 222, 'id': '5be4999502efd0b97cdb5dd6', 'description': '2cte (id:%)'},
		      {'x': 645, 'y': 726, 'id': '5be4a62b02efd0db3ee53fd0', 'description': '2ctj (id:%)'}]
	

	p = [{'tier': 1, **p} for p in pdbs] + [{'tier': 2, **p, 'id': 'SWISSMODEL:'+p['id']} for p in swissmodel]
	
	w = Workshop(p, debug=False).order().save()
	
	##what operations were done??
	display(pd.DataFrame.from_records(w.joints).sort_values(by='aft_resi_start'))

	print('done')




