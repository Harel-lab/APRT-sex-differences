"""Run salmon in a mapping-based mode."""

from multiprocessing import Pool
from pathlib import Path
import subprocess
import argparse
import inspect
import sys

sys.path.append(str(Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent))
from General.os_utils import decompress_file, delete_folder_with_files, find_files, group_pe_fastq_files
from General.consts import final_words


def run_salmon(*, gene_map, salmon_index, in_dir, postfix, parity, mate_prefix, recursive, decompress_cmd, out_dir,
               processes, threads, lib_type):

    # find fastq files
    fastq_files = find_files(in_dir, postfix, recursive) # individual fastq files

    if parity == "se":
        samples_names = [fastq_file.name.split(mate_prefix)[0] for fastq_file in fastq_files]
    else:
        fastq_files = group_pe_fastq_files(fastq_files, postfix, mate_prefix)  # tuples of paired fastq file
        samples_names = [fastq_pair[0].name.split(mate_prefix)[0] for fastq_pair in fastq_files]

    samples_out_dirs = [Path(out_dir, sample_name) for sample_name in samples_names]

    run_one_sample_inputs = [(lib_type, decompress_cmd, fastq, sample_out_dir, salmon_index, gene_map, threads)
                             for fastq, sample_out_dir in zip(fastq_files, samples_out_dirs)]

    if parity == "se":
        one_sample_salmon_func = run_se_salmon_sample
    else:
        one_sample_salmon_func = run_pe_salmon_sample

    # quantify each sample in a different process
    with Pool(processes=processes) as pool:
        pool.starmap(func=one_sample_salmon_func, iterable=run_one_sample_inputs)


def run_se_salmon_sample(lib_type, decompress_cmd, fastq_file, sample_out_dir, salmon_index, gene_map, threads):
    if decompress_cmd:
        sample_out_dir.mkdir(exist_ok=True)
        decompressed_tmp_dir = Path(sample_out_dir, "decompressed_tmp")
        decompressed_tmp_dir.mkdir(exist_ok=True)
        decompressed_fastq_file = Path(decompressed_tmp_dir, ".".join(part for part in fastq_file.name.split(".")[:-1]))
        decompress_file(fastq_file, decompressed_fastq_file, decompress_cmd=decompress_cmd)
        fastq_file = decompressed_fastq_file
    salmon_cmd = f"salmon quant -l {lib_type} -r {fastq_file} -o {sample_out_dir} -i {salmon_index} -g {gene_map} " \
                 f"-p {threads}"
    subprocess.run(salmon_cmd, shell=True)
    if decompressed_tmp_dir.exists():  # if decompress_cmd
        delete_folder_with_files(decompressed_tmp_dir)


def run_pe_salmon_sample(lib_type, decompress_cmd, fastq_files, sample_out_dir, salmon_index, gene_map, threads):
    fastq_file_1 = fastq_files[0]
    fastq_file_2 = fastq_files[1]
    if decompress_cmd:
        sample_out_dir.mkdir(exist_ok=True)
        decompressed_tmp_dir = Path(sample_out_dir, "decompressed_tmp")
        decompressed_tmp_dir.mkdir(exist_ok=True)
        decompressed_fastq_file_1 = Path(decompressed_tmp_dir,
                                         ".".join(part for part in fastq_file_1.name.split(".")[:-1]))
        decompressed_fastq_file_2 = Path(decompressed_tmp_dir,
                                         ".".join(part for part in fastq_file_2.name.split(".")[:-1]))
        decompress_file(fastq_file_1, decompressed_fastq_file_1, decompress_cmd=decompress_cmd)
        decompress_file(fastq_file_2, decompressed_fastq_file_2, decompress_cmd=decompress_cmd)
        fastq_file_1 = decompressed_fastq_file_1
        fastq_file_2 = decompressed_fastq_file_2
    salmon_cmd = f"salmon quant -l {lib_type} -1 {fastq_file_1} -2 {fastq_file_2} -o {sample_out_dir} " \
                 f"-i {salmon_index} -g {gene_map} -p {threads}"
    subprocess.run(salmon_cmd, shell=True)
    if decompressed_tmp_dir.exists():  # if decompress_cmd
        delete_folder_with_files(decompressed_tmp_dir)


if __name__ == "__main__":
    # define args
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene_map",
                        required=True,
                        help="Salmon's gene map")
    parser.add_argument("--salmon_index",
                        required=True,
                        help="Salmon's index dir")
    parser.add_argument("--in_dir",
                        required=True,
                        type=Path,
                        help="A folder with fastq files, where each sample's name format is `$sample_1.fastq.lzma` "
                             "and `$sample_2.fastq.lzma` (or just the first if the files are single-ended")
    parser.add_argument("--postfix",
                        default=".fastq.lzma",
                        help="Postfix of wanted files in `in_dir`")
    parser.add_argument("--parity",
                        default="se",
                        choices=["se", "pe"],
                        help="Reads' parity, whether single-ended (`se`) or paired-ended (`pe`)")
    parser.add_argument("--mate_prefix",
                        default="_",
                        help="See `in_dir` documentation")
    parser.add_argument("--recursive",
                        action="store_true",
                        help="Whether to search recursively in subdirectories of `in_dir` for fastq files")
    parser.add_argument("--decompress_cmd",
                        help="Any command that decompresses the file to stdout, such as `lzcat`, `xzcat`, or "
                             "`gunzip --stdout`. If given, it is assumed files' names end with `.lzma` (for "
                             "example). If absent, the file is uncompressed.")
    parser.add_argument("--out_dir",
                        required=True,
                        type=Path)
    parser.add_argument("--processes",
                        default=6,
                        type=int,
                        help="Number of processes to run in parallel")
    parser.add_argument("--threads",
                        default=8,
                        type=int,
                        help="Number of threads for each process")
    parser.add_argument("--lib_type",
                        default="SR",
                        help="https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype")
    # parse args
    args = parser.parse_args()
    # run
    run_salmon(gene_map=args.gene_map, salmon_index=args.salmon_index, in_dir=args.in_dir, postfix=args.postfix,
               parity=args.parity, mate_prefix=args.mate_prefix, recursive=args.recursive,
               decompress_cmd=args.decompress_cmd, out_dir=args.out_dir, processes=args.processes,
               threads=args.threads, lib_type=args.lib_type)
    # end
    final_words()
