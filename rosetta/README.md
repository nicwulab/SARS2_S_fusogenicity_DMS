# Rosetta modelling of SARS-CoV-2 spike D994Q mutation

## Dependencies

* [Rosetta Software Suite](https://www.rosettacommons.org/software/license-and-download)
    * Either an academic or commercial license is required. One can request a license in the link above.
* [PyMOL](https://pymol.org/2/)
* [pdb-tools](http://www.bonvinlab.org/pdb-tools/)

## Input Files

* D994Q.resfile in the **mutagenesis** folder for mutating the specified residues
* **spike.pdb** containing the structure of prefusion-stabilized, uncleavable SARS-CoV-2 spike for mutagenesis. The original PDB is [6ZGE](https://www.rcsb.org/structure/6ZGE); water and N-acetyl-D-glucosamine (NAG) molecules were removed in PyMOL using 'remove solvent', and 'select nag, ////NAG' followed by 'remove (nag)', respectively. The new molecule was exported as a pdb file with CONECT records written for all bonds, multiple bonds as duplicate CONECT records, and segment identifiers written.

## Steps

1. Renumber the PDB using `python3 pdb_reres.py spike.pdb > spike_renum.pdb` via pdb-tools. Use the renumbered file for subsequent steps.

2. Point mutagenesis was run using `<rosetta_location>/main/source/bin/fixbb.static.macosclangrelease -s spike_renum.pdb -resfile <Mutation>.resfile -nstruct 25` while in the folder where `spike_renum.pdb` is located. One-hundred poses in total were generated.
    * Change `macosclangrelease` depending on your operating system.

3. Examine the .sc score file and select the pose that has the lowest overall score. Copy this file to another folder for relax.

4. Generate a constraint file using `<rosetta_location>/main/source/bin/minimize_with_cst.static.linuxgccrelease -s <pdb file with lowest score after point mutagenesis>.pdb -in:file:fullatom -ignore_unrecognized_res -fa_max_dis 9.0 -database <rosetta_location>/main/database/ -ddg::harmonic_ca_tether 0.5 -score:weights <rosetta_location>/main/database/scoring/weights/pre_talaris_2013_standard.wts -restore_pre_talaris_2013_behavior -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix min_cst_0.5 -ddg::sc_min_only false -score:patch <rosetta_location>/main/database/scoring/weights/score12.wts_patch > mincst.log`.

5. Convert the .log file to a .cst file using `bash <rosetta_location>/main/source/src/apps/public/ddg/convert_to_cst_file.sh ./mincst.log > ./constraint.cst`.

6. With the **constraint.cst** file and the pdb file of the lowest scoring pose from point mutagenesis in the same folder, perform fast relax using `<rosetta_location>/main/source/bin/relax.static.linuxgccrelease -database <rosetta_location>/main/database/ -in:file:s <pdb file with lowest score after point mutagenesis>.pdb -in:file:fullatom -relax:fast -constraints:cst_file constraint.cst -relax:ramp_constraints false -nstruct 5`. Thirty poses in total were generated.
    * Change `macosclangrelease` depending on your operating system.

## Notes

* The parameter _nstruct_ can be varied depending on the number of output poses. In our case, we chose 25 output poses for fixed backbone point mutagenesis (per replicate, for 4 replicates) and 5 output poses (per replicate, for 6 replicates) for fast relax. The lowest scoring pose after fixed backbone mutagenesis is in the **mutagenesis** folder. The lowest scoring pose after fast relax is in the **relax** folder.
