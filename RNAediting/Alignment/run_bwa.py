from multiprocessing import Pool
from pathlib import Path
import subprocess
import argparse
import inspect
import sys
import re

sys.path.append(str(Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent))
from General.os_utils import  find_files, group_pe_fastq_files, delete_folder_with_files, decompress_file
from General.consts import final_words
from Alignment.alignment_utils import sort_alignments, index_alignments, base_alignment_stats


def bwa(in_dir, postfix, parity, mate_prefix, recursive, decompress_cmd, out_dir, index_dir, index_name, bwa_path,
        bwa_algorithm, processes, threads, additional_bwa_flags):

    index_dir = index_dir.absolute()

    # define inputs and alignment function
    fastq_files = find_files(in_dir, postfix, recursive) # individual fastq files
    if parity == "pe":
        fastq_files = group_pe_fastq_files(fastq_files, postfix, mate_prefix)  # tuples of paired fastq file

    if bwa_algorithm == "mem":
        additional_bwa_mem_flags = additional_bwa_flags.get("mem", "")
        if parity == "se":
            starmap_inputs = [(fastq_file, postfix, decompress_cmd, mate_prefix, out_dir, threads, index_name,
                               index_dir, bwa_path, additional_bwa_mem_flags)
                              for fastq_file in fastq_files]
            bwa_func = se_bwa_mem_one_sample
        else:
            starmap_inputs = [(fastq_pair[0],fastq_pair[1], postfix, decompress_cmd, mate_prefix, out_dir, threads,
                               index_name, index_dir, bwa_path, additional_bwa_mem_flags)
                              for fastq_pair in fastq_files]
            bwa_func = pe_bwa_mem_one_sample

    else:   # bwa_algorithm == "backtrack"
        additional_aln_flags = additional_bwa_flags.get("aln", "")
        if parity == "se":
            additional_samse_flags = additional_bwa_flags.get("samse", "")
            starmap_inputs = [(fastq_file, postfix, decompress_cmd, mate_prefix, out_dir, threads, index_name,
                               index_dir, bwa_path, additional_aln_flags, additional_samse_flags)
                              for fastq_file in fastq_files]
            bwa_func = se_bwa_backtrack_one_sample
        else:  # parity is Parity.PE
            additional_sampe_flags = additional_bwa_flags.get("sampe", "")
            starmap_inputs = [(fastq_pair[0], fastq_pair[1], postfix, decompress_cmd, mate_prefix, out_dir, threads,
                               index_name, index_dir, bwa_path, additional_aln_flags, additional_sampe_flags)
                              for fastq_pair in fastq_files]
            bwa_func = pe_bwa_backtrack_one_sample

    # align
    with Pool(processes=processes) as pool:
        pool.starmap(func=bwa_func, iterable=starmap_inputs)


def se_bwa_mem_one_sample(fastq_file, postfix, decompress_cmd, mate_prefix, out_dir, threads, index_name, index_dir,
                          bwa_path, additional_bwa_mem_flags):
    sample_name = fastq_file.name.rstrip(postfix).split(mate_prefix)[0]
    if decompress_cmd:
        decompressed_tmp_dir = Path(out_dir, f"{sample_name}_decompressed_tmp")
        decompressed_tmp_dir.mkdir(exist_ok=True)
        decompressed_fastq_file = Path(decompressed_tmp_dir, ".".join(part for part in fastq_file.name.split(".")[:-1]))
        decompress_file(fastq_file, decompressed_fastq_file, decompress_cmd=decompress_cmd)
        fastq_file = decompressed_fastq_file
    alignment_file = Path(out_dir, f"{sample_name}.Aligned.sam").absolute()
    fastq_file = fastq_file.absolute()
    if additional_bwa_mem_flags:
        alignment_cmd = f"{bwa_path} mem {additional_bwa_mem_flags} -t {threads} {index_name} {fastq_file} " \
                        f"-o {alignment_file}"
    else:
        alignment_cmd = f"{bwa_path} mem -t {threads} {index_name} {fastq_file} -o {alignment_file}"
    subprocess.run(alignment_cmd, shell=True, cwd=index_dir)
    if decompressed_tmp_dir.exists():  # if decompress_cmd
        delete_folder_with_files(decompressed_tmp_dir)


def pe_bwa_mem_one_sample(fastq_file_1, fastq_file_2, postfix, decompress_cmd, mate_prefix, out_dir, threads, index_name,
                          index_dir, bwa_path, additional_bwa_mem_flags):
    sample_name = fastq_file_1.name.rstrip(postfix).split(mate_prefix)[0]  # fastq_file_2 might be used just as well
    if decompress_cmd:
        decompressed_tmp_dir = Path(out_dir, f"{sample_name}_decompressed_tmp")
        decompressed_tmp_dir.mkdir(exist_ok=True)
        decompressed_fastq_file_1 = Path(decompressed_tmp_dir,
                                         ".".join(part for part in fastq_file_1.name.split(".")[:-1]))
        decompressed_fastq_file_2 = Path(decompressed_tmp_dir,
                                         ".".join(part for part in fastq_file_2.name.split(".")[:-1]))
        decompress_file(fastq_file_1, decompressed_fastq_file_1, decompress_cmd=decompress_cmd)
        decompress_file(fastq_file_2, decompressed_fastq_file_2, decompress_cmd=decompress_cmd)
        fastq_file_1 = decompressed_fastq_file_1
        fastq_file_2 = decompressed_fastq_file_2
    alignment_file = Path(out_dir, f"{sample_name}.Aligned.sam").absolute()  # fastq_file_2 might be used just as well
    fastq_file_1 = fastq_file_1.absolute()
    fastq_file_2 = fastq_file_2.absolute()
    if additional_bwa_mem_flags:
        alignment_command = f"{bwa_path} mem {additional_bwa_mem_flags } -t {threads} {index_name} {fastq_file_1} " \
                            f"{fastq_file_2} -o {alignment_file}"
    else:
        alignment_command = f"{bwa_path} mem -t {threads} {index_name} {fastq_file_1} {fastq_file_2} " \
                            f"-o {alignment_file}"
    subprocess.run(alignment_command, shell=True, cwd=index_dir)
    if decompressed_tmp_dir.exists():  # if decompress_cmd
        delete_folder_with_files(decompressed_tmp_dir)


def se_bwa_backtrack_one_sample(fastq_file, postfix, decompress_cmd, mate_prefix, out_dir, threads, index_name, index_dir,
                                bwa_path, additional_aln_flags, additional_samse_flags):
    # (1) bwa aln ref.fa reads.fq > reads.sai
    # (2) bwa samse ref.fa reads.sai reads.fq > aln-se.sam
    # define paths
    sample_name = fastq_file.name.rstrip(postfix).split(mate_prefix)[0]
    if decompress_cmd:
        decompressed_tmp_dir = Path(out_dir, f"{sample_name}_decompressed_tmp")
        decompressed_tmp_dir.mkdir(exist_ok=True)
        decompressed_fastq_file = Path(decompressed_tmp_dir, ".".join(part for part in fastq_file.name.split(".")[:-1]))
        decompress_file(fastq_file, decompressed_fastq_file, decompress_cmd=decompress_cmd)
        fastq_file = decompressed_fastq_file
    temp_sai_alignment_file = Path(out_dir, f"{fastq_file.name}.sai").absolute()
    final_sam_alignment_file = Path(out_dir, f"{fastq_file.name}.Aligned.sam").absolute()
    fastq_file = fastq_file.absolute()
    # define alignment commands
    if additional_aln_flags:
        aln_cmd = f"{bwa_path} aln -t {threads} {additional_aln_flags} {index_name} {fastq_file} > " \
                  f"{temp_sai_alignment_file}"
    else:
        aln_cmd = f"{bwa_path} aln -t {threads} {index_name} {fastq_file} > {temp_sai_alignment_file}"
    if additional_samse_flags:
        samse_cmd = f"{bwa_path} samse {additional_samse_flags} {index_name} {temp_sai_alignment_file} {fastq_file} " \
                    f"> {final_sam_alignment_file}"
    else:
        samse_cmd = f"{bwa_path} samse {index_name} {temp_sai_alignment_file} {fastq_file} > {final_sam_alignment_file}"
    # run
    subprocess.run(aln_cmd, shell=True, cwd=index_dir)
    subprocess.run(samse_cmd, shell=True, cwd=index_dir)
    # delete temp sai file
    temp_sai_alignment_file.unlink()
    if decompressed_tmp_dir.exists():  # if decompress_cmd
        delete_folder_with_files(decompressed_tmp_dir)


def pe_bwa_backtrack_one_sample(fastq_file_1, fastq_file_2, postfix, decompress_cmd, mate_prefix, out_dir, threads, index_name,
                                bwa_path, index_dir, additional_aln_flags, additional_sampe_flags):
    # (1) bwa aln ref.fa read1.fq > read1.sai; bwa aln ref.fa read2.fq > read2.sai
    # (2) bwa sampe ref.fa read1.sai read2.sai read1.fq read2.fq > aln-pe.sam
    # define paths
    sample_name = fastq_file_1.name.rstrip(postfix).split(mate_prefix)[0]  # fastq_file_2 might be used just as well
    if decompress_cmd:
        decompressed_tmp_dir = Path(out_dir, f"{sample_name}_decompressed_tmp")
        decompressed_tmp_dir.mkdir(exist_ok=True)
        decompressed_fastq_file_1 = Path(decompressed_tmp_dir,
                                         ".".join(part for part in fastq_file_1.name.split(".")[:-1]))
        decompressed_fastq_file_2 = Path(decompressed_tmp_dir,
                                         ".".join(part for part in fastq_file_2.name.split(".")[:-1]))
        decompress_file(fastq_file_1, decompressed_fastq_file_1, decompress_cmd=decompress_cmd)
        decompress_file(fastq_file_2, decompressed_fastq_file_2, decompress_cmd=decompress_cmd)
        fastq_file_1 = decompressed_fastq_file_1
        fastq_file_2 = decompressed_fastq_file_2
    temp_sai_1_alignment_file = Path(out_dir, f"{fastq_file_1.name}.sai").absolute()
    temp_sai_2_alignment_file = Path(out_dir, f"{fastq_file_2.name}.sai").absolute()
    final_sam_alignment_file = Path(out_dir, f"{sample_name}.Aligned.sam").absolute()
    fastq_file_1 = fastq_file_1.absolute()
    fastq_file_2 = fastq_file_2.absolute()
    # define alignment commands
    if additional_aln_flags:
        aln_1_cmd = f"{bwa_path} aln -t {threads} {additional_aln_flags} {index_name} {fastq_file_1} > " \
                    f"{temp_sai_1_alignment_file}"
        aln_2_cmd = f"{bwa_path} aln -t {threads} {additional_aln_flags} {index_name} {fastq_file_2} > " \
                    f"{temp_sai_2_alignment_file}"
    else:
        aln_1_cmd = f"{bwa_path} aln -t {threads} {index_name} {fastq_file_1} > {temp_sai_1_alignment_file}"
        aln_2_cmd = f"{bwa_path} aln -t {threads} {index_name} {fastq_file_2} > {temp_sai_2_alignment_file}"
    if additional_sampe_flags:
        sampe_cmd = f"{bwa_path} sampe {additional_sampe_flags} {index_name} {temp_sai_1_alignment_file} " \
                    f"{temp_sai_2_alignment_file} {fastq_file_1} {fastq_file_2} > {final_sam_alignment_file}"
    else:
        sampe_cmd = f"{bwa_path} sampe {index_name} {temp_sai_1_alignment_file} {temp_sai_2_alignment_file} " \
                    f"{fastq_file_1} {fastq_file_2} > {final_sam_alignment_file}"
    # run
    subprocess.run(aln_1_cmd, shell=True, cwd=index_dir)
    subprocess.run(aln_2_cmd, shell=True, cwd=index_dir)
    subprocess.run(sampe_cmd, shell=True, cwd=index_dir)
    # delete temp sai files
    temp_sai_1_alignment_file.unlink()
    temp_sai_2_alignment_file.unlink()
    if decompressed_tmp_dir.exists():  # if decompress_cmd
        delete_folder_with_files(decompressed_tmp_dir)


def main(*, in_dir, postfix, parity, mate_prefix, recursive, decompress_cmd, out_dir, index_name, index_dir,
         bwa_path, bwa_algorithm, processes, threads, additional_bwa_flags, samtools_path):
    out_dir.mkdir(exist_ok=True)
    bwa_path = bwa_path.expanduser()
    samtools_path = samtools_path.expanduser()
    # parse additional flags for the different algorithms
    additional_bwa_flags = parse_additional_bwa_flags(additional_bwa_flags)
    # align
    bwa(in_dir, postfix, parity, mate_prefix, recursive, decompress_cmd, out_dir, index_dir, index_name, bwa_path,
        bwa_algorithm, processes, threads, additional_bwa_flags)
    # sort alignments
    sort_alignments(in_dir=out_dir, postfix=".Aligned.sam", recursive=False, out_dir=out_dir, processes=processes,
                    threads=threads, samtools_path=samtools_path)
    # base alignment stats with samtools
    base_alignment_stats(in_dir=out_dir, postfix=".Aligned.SortedByCoords.bam", recursive=False, out_dir=out_dir,
                         processes=processes, threads=threads, samtools_path=samtools_path)
    # index alignment (create .bai files)
    index_alignments(in_dir=out_dir, postfix=".Aligned.SortedByCoords.bam", recursive=False, out_dir=out_dir,
                     processes=processes, threads=threads, samtools_path=samtools_path)


def parse_additional_bwa_flags(additional_bwa_flags):
    """
    Parse additional_bwa_flags, which come as a list with dashes prefixed by "\\", to a dictionary of the bwa programs.

    The order of "mem", "aln", "samse" and "sampe" is not important, but of course all additional args of "mem" (for
    example) should proceed those of "aln".

    @param additional_bwa_flags: a list of strings, which might look something like this:
    ['mem', '\\-a', '\\-w', '1', 'aln', '\\-n', '0.7']
    @type additional_bwa_flags: list[str]
    @return: additional_bwa_flags
    @rtype: dict
    """
    # check if additional_bwa_flags was supplied, and return empty dictionary if not
    if not additional_bwa_flags:
        return {}
    # current programs wrapped in this module
    programs = ("mem", "aln", "samse", "sampe")
    # remove "\\" prefixing "-", which was needed for argparse
    additional_flags_revived_dashes = []
    for arg in additional_bwa_flags:
        arg = re.sub(r"\\", "", arg)
        additional_flags_revived_dashes.append(arg)
    # find start locations of programs in the list
    programs_start_location = {program: None for program in programs}
    for program in programs:
        try:
            program_start_location = additional_flags_revived_dashes.index(program)
            programs_start_location[program] = program_start_location
        except ValueError:
            pass
    # define present programs
    present_programs = [program for program in programs_start_location if programs_start_location[program]]
    # sort present_programs_start_location according to their start location
    present_programs_start_location = [(program, programs_start_location[program]) for program in present_programs]
    present_programs_start_location.sort(key=lambda x: x[1])
    # define range of args in additional_flags_revived_dashes, according to each program in
    # present_programs_start_location
    present_programs_start_end_locations = {}
    for program_num, (program, start_location) in enumerate(present_programs_start_location, start=1):
        start_location += 1    # as the first arg is the name of program itself, and we already know it
        if program_num == len(present_programs_start_location):
            end_location = len(additional_flags_revived_dashes)
        else:
            _, next_start_location = present_programs_start_location[program_num + 1]  # _ == next_program_name
            end_location = next_start_location
        present_programs_start_end_locations[program] = (start_location, end_location)
    # final assembly of additional_bwa_flags according to parsed present_programs_start_end_locations
    final_additional_bwa_flags = {}
    for program in present_programs_start_end_locations:
        start_location, end_location = present_programs_start_end_locations[program]
        program_args = additional_flags_revived_dashes[start_location:end_location]
        program_args = " ".join(program_args)
        final_additional_bwa_flags[program] = program_args
    # return
    additional_bwa_flags = final_additional_bwa_flags
    return additional_bwa_flags


if __name__ == "__main__":

    # create parser
    parser = argparse.ArgumentParser()

    # define args
    parser.add_argument("--in_dir",
                        required=True,
                        type=Path,
                        help="A folder with fastq files, where each sample's name format is something like "
                             "`$sample_1.fastq.lzma` and `$sample_2.fastq.lzma` (or just the first if the files are "
                             "single-ended).")
    parser.add_argument("--postfix",
                        default=".fastq",
                        help="Postfix of wanted files in `in_dir`.")
    parser.add_argument("--parity",
                        default="se",
                        choices=["se", "pe"],
                        help="Reads' parity, whether single-ended (`se`) or paired-ended (`pe`).")
    parser.add_argument("--mate_prefix",
                        default="_",
                        help="See `in_dir` documentation.")
    parser.add_argument("--recursive",
                        action="store_true",
                        help="Whether to search recursively in subdirectories of `in_dir` for fastq files.")
    parser.add_argument("--decompress_cmd",
                        help="Any command that decompresses the file to stdout, such as `lzcat`, `xzcat`, or "
                             "`gunzip --stdout`. If given, it is assumed files' names end with `.lzma` (for "
                             "example). If absent, the file is uncompressed.")

    parser.add_argument("--out_dir",
                        type=Path,
                        required=True)

    parser.add_argument("--index_dir",
                        type=Path,
                        required=True,
                        help="Where the bwa index resides.")
    parser.add_argument("--index_name",
                        required=True,
                        help="The name of the index, as created by the `bwa index` command.")
    parser.add_argument("--bwa_path",
                        type=Path,
                        default="~/anaconda3/envs/Killifish-APRT-RNAEditing/bin/bwa")
    parser.add_argument("--bwa_algorithm",
                        default="mem",
                        choices=["mem", "backtrack"],
                        help="Use `mem` for reads longer than ~70bp (default), or `backtrack` for reads shorter "
                             "than ~70bp")
    parser.add_argument("--additional_bwa_flags", nargs="+",
                        help="Carefully read the documentation of 'parse_additional_bwa_flags' in order to understand"
                             " how additional_bwa_flags is parsed, and therefore how to construct it")
    parser.add_argument("--samtools_path",
                        type=Path,
                        default="~/anaconda3/envs/Killifish-APRT-RNAEditing/bin/samtools")

    parser.add_argument("--processes",
                        type=int,
                        default=6,
                        help="Number of processes to run in parallel.")
    parser.add_argument("--threads",
                        type=int,
                        default=10,
                        help="Threads used in each subprocess.")

    # parse args
    args = parser.parse_args()

    # run
    main(in_dir=args.in_dir, postfix=args.postfix, parity=args.parity, mate_prefix=args.mate_prefix,
         recursive=args.recursive, decompress_cmd=args.decompress_cmd, out_dir=args.out_dir, index_name=args.index_name,
         index_dir=args.index_dir, bwa_path=args.bwa_path, bwa_algorithm=args.bwa_algorithm,
         samtools_path=args.samtools_path, processes=args.processes, threads=args.threads,
         additional_bwa_flags=args.additional_bwa_flags)

    # end
    final_words()
