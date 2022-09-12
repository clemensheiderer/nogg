
EggNOG:  4.5
author: clemens




positional arguments:
  speciesfile - eggnog4.species_list.tsv
  membersfile - meNOG.members.tsv
  annotations - meNOG.members.tsv

optional arguments:
  -h, --help            show this help message and exit
  -i, --intersection
                            load it with: python3 finalclemens.py eggnog4.species_list.tsv meNOG.members.tsv meNOG.annotations.tsv -i

                            Enter 2 mammal species your choice.
                            All genes/ProteinIds of first species are calculated, which have at least one homolog in the second species.

  -t, --three_species
                            load it with: python3 finalclemens.py eggnog4.species_list.tsv meNOG.members.tsv meNOG.annotations.tsv -t

                            The relative compliment of Homo sapiens in Mu musculus is made, intersection with Pan troglodytes generates
                            all needed indices.

  -f, --family_specific_genes

                            load it with: python3 finalclemens.py eggnog4.species_list.tsv meNOG.members.tsv meNOG.annotations.tsv -i -f

                            Mus Musculus and Rattus norvegicus share the same Muridae family.