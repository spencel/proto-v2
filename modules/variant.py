
import json
import os

from devtools import debug


# Config
bio_config = json.load(open(os.path.join("data", "bio.json"), 'r'))


class Variant():


    @classmethod
    def reconstruct_sequence_alignment(cls,
        reference_sequence = None, # The reference nucleic acid sequence
        reference_offset_origin = None, # Where the reference sequence starts if it's a subsection of the whole sequence
        reference_segments = None, # A list of these
        variant_segments   = None, # A list of these
        variant_segment_start_loci_on_reference_sequence = None, # A list of these
        mutation_types = None # A list of these
    ):
        scope_pref = f"{cls.__name__}.reconstruct_sequence_alignment"

        # Case 1
        if reference_sequence and variant_segments and variant_segment_start_loci_on_reference_sequence and mutation_types:

            # Handle arguments
            na_sequence = reference_sequence
            tech_seq_start = reference_offset_origin
            ref_sequences = reference_segments
            var_sequences = variant_segments
            ref_starts = variant_segment_start_loci_on_reference_sequence
            mut_types = mutation_types

            # Loop through variant segments and create reference and variant technical sequence alignment
            # Set the offset to account for insertions
            ref_offset = 0
            var_offset = 0
            # Initialize nucleotide sequences
            ref_tech_sequence = na_sequence
            var_tech_sequence = na_sequence
            for i in range(len(ref_starts)):
                # Alias variant start relative to technical sequence, not start on variant genome
                ref_start = ref_starts[i] - tech_seq_start + ref_offset
                var_start = ref_starts[i] - tech_seq_start + var_offset
                ref_sequence = ref_sequences[i]
                var_sequence = var_sequences[i]
                mut_type = mut_types[i]

                # Handle insertion
                if mut_type == bio_config["mutation_types"]["insertion"]:
                    ref_tech_sequence = ref_tech_sequence[:ref_start] + ref_sequence + ref_tech_sequence[ref_start:]
                    var_tech_sequence = var_tech_sequence[:var_start] + var_sequence + var_tech_sequence[var_start:]
                    ref_offset += len(ref_sequence)

                # Handle deletion
                if mut_type == bio_config["mutation_types"]["deletion"]:
                    var_tech_sequence = \
                        var_tech_sequence[:ref_start] \
                        + var_sequence \
                        + var_tech_sequence[ref_start + len(var_sequence):]

                # Handle substitution
                elif mut_type == bio_config["mutation_types"]["substitution"]:
                    var_tech_sequence = \
                          var_tech_sequence[:ref_start] \
                        + var_sequence \
                        + var_tech_sequence[ref_start + len(var_sequence):]

            return ref_tech_sequence, var_tech_sequence

