from pathlib import Path


class GenomeReader:
    @staticmethod
    def parse_genome(genome_file):
        """
        Parse genome to dict of {chromosome_name: sequence}.

        @param genome_file: path to a genome fasta file
        @type genome_file: str or Path
        @return: genome, dict of {chromosome_name: sequence}
        @rtype: dict
        """
        genome = {}
        region_name = ""
        seqs = []
        with Path(genome_file).open("r") as genome_file:
            for line_num, line in enumerate(genome_file, start=1):
                line = line.rstrip("\n")
                if not line:  # skip empty lines
                    continue
                if line.startswith(">"):    # identifier line
                    if region_name and seqs:
                        merged_seq = "".join(seq for seq in seqs)
                        genome[region_name] = merged_seq
                        seqs = []
                    region_name = line.split(" ")[0][1:]
                else:  # seq line
                    if not region_name:
                        raise ValueError(
                            f"Error at line number {line_num} of genome_file, a seq is present before " +
                            f"identifier")
                    seqs.append(line)
        if region_name and seqs:
            merged_seq = "".join(seq for seq in seqs)
            genome[region_name] = merged_seq
        return genome
