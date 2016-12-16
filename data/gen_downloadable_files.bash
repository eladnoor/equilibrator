#!/usr/bin/bash

python data/export_formation_energies.py --ph 5.0 -i 0.0 -o compound_formation_energies
python data/export_formation_energies.py --ph 5.5 -i 0.0 -o compound_formation_energies
python data/export_formation_energies.py --ph 6.0 -i 0.0 -o compound_formation_energies
python data/export_formation_energies.py --ph 6.5 -i 0.0 -o compound_formation_energies
python data/export_formation_energies.py --ph 7.0 -i 0.0 -o compound_formation_energies
python data/export_formation_energies.py --ph 7.5 -i 0.0 -o compound_formation_energies
python data/export_formation_energies.py --ph 8.0 -i 0.0 -o compound_formation_energies
python data/export_formation_energies.py --ph 8.5 -i 0.0 -o compound_formation_energies
python data/export_formation_energies.py --ph 9.0 -i 0.0 -o compound_formation_energies

mv compound_formation_energies*.json static/downloads/
mv compound_formation_energies*.tsv static/downloads/


