import subprocess


def index(root_dir, bam_files_suffix, follow_links, stranded, paired_end, gene_expression_file, refseq_file,
          groups_file, genome_file, regions_file, sub_out_dir, mismatches, sample_threads,
          sample_strands_threads, index_path):

    index_cmd = f"{index_path} -d {root_dir} -f {bam_files_suffix}"

    bool_flags_values = [follow_links, stranded, paired_end]
    if any(bool_flags_values):
        bool_flags = ["--follow_links", "--paired_end", "--stranded "]
        given_bool_flags = " ".join(f"{bool_flag}" for bool_flag, bool_flag_value in zip(bool_flags, bool_flags_values)
                                    if bool_flag_value)
        index_cmd += f" {given_bool_flags}"

    index_cmd += f" --genes_expression {gene_expression_file} --genome 'UserProvided' --refseq {refseq_file} " \
                 f"-g {groups_file} -gf {genome_file} -rb {regions_file} --per_region_output " \
                 f"--per_sample_output -mm {mismatches} -o {sub_out_dir} -os {sub_out_dir} -l {sub_out_dir} " \
                 f"--ts {sample_threads} --tsd {sample_strands_threads}"

    subprocess.run(index_cmd, shell=True)
