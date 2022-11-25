# Sidekick
Sidekick software package for automating transmembrane helix simulation 

http://sbcb.bioch.ox.ac.uk/Sidekick/

## Licence

Free for non-commercial use. Check Sidekick_4.5/licence.doc for licence information specifics.

## Requirements

Install a version of gromacs into Sidekick_4.5/gromacs. This is available from gromacs.org under its own licence.

`https://ftp.gromacs.org/pub/gromacs/`

## Usage

Example of basic usage for a single helix simulation
CG_Helix.py -s GWWLALALALALALALWWA -l 100 -r 5 -t MARTINI_1.1.2.b --special ntermini

## Install tips

On ubuntu bash links to dash- this can be worked around by switching manually with

`dkpg-reconfigure dash`

Set variable values in HAConf.py to point to your install and data locations

You get a password error if you don't set a stored password- this is legacy (used in Xgrid) and can be generated as below

python2 set_password.py sometext
