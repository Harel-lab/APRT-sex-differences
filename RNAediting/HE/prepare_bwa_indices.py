# TODO fix - didn't use it that way in the Killifish project


# from pathlib import  Path
# import argparse
#
#
#
#
#
# class GeneralHyperEditingTransformer():
#     """HyperEditingTransformer transforms genome for Hyper Editing according to given transformation.
#
#         For example, if transformation == 'A2G', then a seq 'ATCCT' will become 'GTCCT'."""
#
#     def __init__(self, genes, transformation):
#         self.genes = genes
#         self.prev_base, self.new_base = transformation.upper().split("2")  # transformation should already be uppercase
#         self.transformed = []
#
#     def process(self):
#         for isoform in self.genes:
#             new_gene = FastaGene(isoform.get_header(), str(isoform).replace(self.prev_base, self.new_base))
#             self.transformed.append(new_gene)
#
#     def get_results(self):
#         for trans in self.transformed:
#             yield trans
#
#
#
#
# def general_he_transformations(args):
#     genome = args.input_files.pop()
#     output_dir = args.output_dir
#     genome_prefix = args.transcriptome_prefix
#     _, genome_name = os.path.split(genome)
#     genome_name = genome_name.split(".")[0]
#     transformed_paths = [os.path.join(output_dir, f"{genome_name}.{trans.lower()}") for trans in
#                          HE_TRANSFORMATIONS]
#     writers = [FileWriter(path) for path in transformed_paths]
#     for writer, transformation in zip(writers, HE_TRANSFORMATIONS):
#         fasta_reader = FastaReader(input_file_path=genome)
#         genes = fasta_reader.objects_from_file()  # genes == generator
#         transformer = GeneralHyperEditingTransformer(genes, transformation)
#         transformer.process()
#         transformed_genes = transformer.get_results()
#         writer.write_objects(transformed_genes)
#     func_params = {"transcriptome_prefix": genome_prefix}
#     create_bwa_index(output_dir, func_params=func_params, max_sub_processes=6)
#     print(FINAL_WORDS)
#
#
# class FileWriter:
#     """docstring for FileWriter"""
#     def __init__(self, output_file, fields=None):
#         self.output_file = output_file
#         self.fields = fields
#
#     def write_lines(self, lines):
#         with open(self.output_file, mode="wt") as new_file:
#             for line in lines:
#                 new_file.write(line + "\n")   # presumably line doesn't end with '\n'
#
# class GenomeReader:
#     """GenomeReader is a simpler and more general version of FastaReader."""
#
#     @staticmethod
#     def parse_genome(genome_file):
#         """
#         Parse genome to dict of {chromosome_name: sequence}.
#
#         @param genome_file: path to a genome fasta file
#         @type genome_file: str
#         @return: genome, dict of {chromosome_name: sequence}
#         @rtype: dict
#         """
#         genome = {}
#         region_name = ""
#         seqs = []
#         with Path(genome_file).open("r") as genome_file:
#             for line_num, line in enumerate(genome_file, start=1):
#                 line = line.rstrip("\n")
#                 if not line:  # skip empty lines
#                     continue
#                 if line.startswith(">"):    # identifier line
#                     if region_name and seqs:
#                         merged_seq = "".join(seq for seq in seqs)
#                         genome[region_name] = merged_seq
#                         seqs = []
#                     region_name = line.split(" ")[0][1:]
#                 else:  # seq line
#                     if not region_name:
#                         raise ValueError(
#                             f"Error at line number {line_num} of genome_file, a seq is present before " +
#                             f"identifier")
#                     seqs.append(line)
#         if region_name and seqs:
#             merged_seq = "".join(seq for seq in seqs)
#             genome[region_name] = merged_seq
#         return genome
#
#
#     @staticmethod
#     def verbose_parse_genome(genome_file):
#         """
#         Parse genome to dict of header_line, seq.
#
#         @param genome_file: path to genome file
#         @type genome_file: str
#         @return: genome, dict of {header_line: seq}
#         @rtype: dict
#         """
#         genome = {}
#         header_line = ""
#         partial_seqs = []
#         with Path(genome_file).open("r") as genome_file:
#             for line_num, line in enumerate(genome_file, start=1):
#                 line = line.rstrip("\n")
#                 if not line:  # skip empty lines
#                     continue
#                 if line.startswith(">"):    # identifier line
#                     if header_line and partial_seqs:
#                         seq = "".join(partial_seq for partial_seq in partial_seqs)
#                         genome[header_line] = seq
#                         partial_seqs = []
#                     header_line = line
#                 else:  # seq line
#                     if not header_line:
#                         raise ValueError(
#                             f"Error at line number {line_num} of genome_file, a seq is present before " +
#                             f"identifier")
#                     partial_seqs.append(line)
#         if header_line and partial_seqs:
#             seq = "".join(partial_seq for partial_seq in partial_seqs)
#             genome[header_line] = seq
#         return genome
#
#
# class Seq(ABC):
#     """Seq is an abstract class for an object representing a sequence."""
#
#     def __init__(self, header, seq):
#         self.header = header
#         self.seq = seq  # pay attention - sometimes lower-case letters might be used by purpose
#
#     def __repr__(self):
#         if len(self.seq) <= 12:
#             seq_repr = self.seq
#         else:
#             seq_repr = f"{self.seq[0:6]}...{self.seq[len(self.seq) - 6:len(self.seq)]}"
#         return f"<{self.get_header()}   {seq_repr}>"
#
#     def __str__(self):
#         """Return the sequence of the object."""
#         return self.seq
#
#     def __len__(self):
#         """Return the the length of the sequence."""
#         return len(self.seq)
#
#     def get_header(self):
#         """Return the header of the object."""
#         return self.header
#
#     def get_lines_for_file(self):
#         """Return header and seq, in separate lines, so they can be conveniently written to file."""
#         return self.header + "\n" + self.seq + "\n"
#
#     def complement(self, seq=None):
#         """Return a string of the complementary sequence, e.g., for 'ATCG' return 'TAGC'."""
#         if seq is None:
#             seq = self.seq
#         replacement_0 = seq.upper()
#         replacement_1 = replacement_0.replace('A', 't')
#         replacement_2 = replacement_1.replace('T', 'a')
#         replacement_3 = replacement_2.replace('C', 'g')
#         replacement_4 = replacement_3.replace('G', 'c')
#         return replacement_4.upper()
#
#     def reverse_complement(self, seq=None):
#         """Return a string of the reversed complementary sequence, e.g., for 'ATCG' return 'CGAT'."""
#         if seq is None:
#             return "".join(reversed(self.complement()))
#         return "".join(reversed(self.complement(seq)))
#
#     def transcribe(self, seq=None):
#         """Return a string of the equivalent RNA seq, e.g., for 'ATCG' return 'AUCG'."""
#         if seq is None:
#             seq = self.seq
#         seq = seq.upper()
#         return seq.replace("T", "U")
#
#     def config_iteration_params(self, first_base=None, last_base=None, inclusive=True, first_base_is_1=True,
#                                 seq=None, positive_strand=True, transcribe=False):
#         """
#         Configure parameters for iterating over the seq (usually as triplets/codons for translation).
#
#         Note to supply both first_base and last_base together or not at all. If they are not provided then the whole seq
#         is used (as is frame +1). If they are provided, it's assumed (inclusive=True) they are both inclusive. Else,
#         if inclusive==False, last_base is exclusive.
#         In addition, if they are provided, it's assumed (first_base_is_1=True) that the first base whose index is 1
#         is found at index 0 in the string. This default behavior is suited for dealing with blast results.
#
#         If seq is not provided then the object's seq (self.seq) is used.
#         In that case, if first_base and last_base are not provided too, it's assumed (positive_strand=True) that the
#         positive strand is required.
#
#         The length of the returned seq isn't changed, only the order in case transcribe==True.
#         Regardless of the input flags, the returned first_base is always inclusive and returned last_base is
#         always exclusive, and first_base <= last_base. Also, if first_base points to the first base of the seq, its
#         index will be 0 (first_base_is_1==False).
#         """
#         strand = Strand.FORWARD
#         if seq is None:
#             seq = self.seq
#         length = len(seq)
#         if first_base is not None and last_base is not None:
#             # convert to 0-index notation
#             if first_base_is_1:
#                 first_base -= 1
#                 last_base -= 1
#             # swap first and last indices
#             if last_base < first_base:
#                 strand = Strand.REVERSE
#                 # temporarily convert to inclusive notation
#                 if not inclusive:
#                     first_base -= 1
#                     inclusive = True
#                 first_base = length - first_base - 1
#                 last_base = length - last_base - 1
#             # finally convert to exclusive notation
#             if inclusive:
#                 last_base += 1
#         elif (first_base is None and last_base is not None) or (first_base is not None and last_base is None):
#             raise Exception("first_base and last_base must be used together, or not at all")
#         else:  # using self.seq
#             if not positive_strand:
#                 strand = strand.REVERSE
#             first_base = 0
#             last_base = length
#         if strand == strand.REVERSE:
#             seq = self.reverse_complement(seq)
#         if transcribe:
#             seq = self.transcribe(seq)
#         return seq, first_base, last_base
#
#     def iterate_codons(self, first_base=None, last_base=None, inclusive=True, first_base_is_1=True,
#                        seq=None, positive_strand=True):
#         """Generate RNA codons, e.g., for 'ATCGGGTTC' yield 'AUC','GGG','UUC'."""
#         params = self.config_iteration_params(first_base=first_base, last_base=last_base, inclusive=inclusive,
#                                               first_base_is_1=first_base_is_1, seq=seq, positive_strand=positive_strand,
#                                               transcribe=True)
#         seq, first_base, last_base = params
#         current_base = first_base
#         while current_base + 3 <= last_base:
#             yield seq[current_base:current_base + 3]
#             current_base += 3
#
#     def translate(self, first_base=None, last_base=None, inclusive=True, first_base_is_1=True,
#                   seq=None, positive_strand=True):
#         """Return a string representing transcription and translation of DNA to AAs.
#
#         For example, if the sequence is 'ATCGGGTTCTAG' then the method returns 'IGF*', where '*' is a stop codon"""
#         codon_iterator = self.iterate_codons(first_base=first_base, last_base=last_base, inclusive=inclusive,
#                                              first_base_is_1=first_base_is_1, seq=seq, positive_strand=positive_strand)
#         return "".join([RNA_TO_AA_DICT[codon] for codon in codon_iterator])
#
#     def set_seq(self, seq):
#         self.seq = seq
#
#
# class FastaGene(Seq):
#     def __init__(self, header, seq):
#         super().__init__(header, seq)