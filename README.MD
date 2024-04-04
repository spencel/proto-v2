
# Todo
* Get deflines
* Figure out a best way to add ways to execute CLI commands for a module.
* Fasta analyzer: create taxon analyzer, check for allowed and prohibited clades, species, etc.


# Commands
Use [Python Fire](https://google.github.io/python-fire/guide/) to call functions and class methods.
```
# Example A
python -m fire ./modules/fasta.py Fasta get_deflines --fpath=/media/sl/T7-Shield/PaRTISeq/blast_db_genomes/Step1_20211004.fa

# Example B
python -m fire ./commands/parti_seq/validate_db_fastas.py validate_db_fastas --arg-fasta-fpath=/media/sl/T7-Shield/PaRTISeq/blast_db_genomes/Step1_20211004.fa
```