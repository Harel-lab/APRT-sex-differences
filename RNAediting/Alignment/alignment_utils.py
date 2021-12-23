from multiprocessing import Pool
from pathlib import Path
import subprocess
import inspect
import sys

sys.path.append(str(Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent))
from General.os_utils import delete_folder_with_files, find_files


def sort_alignments(in_dir, postfix, recursive, out_dir, processes, threads, samtools_path):
    """
    Sort alignment sam files and save them as bam files.

    @param in_dir: where the alignment files reside
    @type in_dir: Path
    @param postfix: find alignment files that end with postfix
    @type postfix: str
    @param recursive: whether to search recursively for alignments files in sub dirs of in_dir
    @type recursive: bool
    @param out_dir: write sorted bam files to out_dir
    @type out_dir: str or Path
    @param processes: number of processes to run in parallel
    @type processes: int
    @param threads: number of threads to use in each process
    @type threads: int
    @param samtools_path: path to a samtools' executable
    @type samtools_path: str or Path
    """
    alignment_files = find_files(in_dir, postfix, recursive)
    starmap_inputs = [(alignment_file, out_dir, threads, samtools_path) for alignment_file in alignment_files]
    with Pool(processes=processes) as pool:
        pool.starmap(func=sort_one_alignment_file, iterable=starmap_inputs)


def sort_one_alignment_file(alignment_file, out_dir, threads, samtools_path):
    """
    Sort alignment_file and save it as a bam file in out_dir.

    @param alignment_file: an alignment sam file
    @type alignment_file: Path
    @param out_dir: write sorted bam file to out_dir
    @type out_dir: Path
    @param threads: number of threads to use
    @type threads: int
    @param samtools_path: path to a samtools' executable
    @type samtools_path: str or Path
    """
    sample_name = ".".join(part for part in alignment_file.name.split(".")[:-1])  # remove ".sam" or ".bam"
    temp_sorting_dir = Path(out_dir, f"{sample_name}.TempSortingFiles")
    temp_sorting_dir.mkdir(exist_ok=True)
    sorted_alignment_file = Path(out_dir, f"{sample_name}.SortedByCoords.bam")
    sorting_command = f"{samtools_path} sort -@ {threads} -T {temp_sorting_dir} -o {sorted_alignment_file} " \
                      f"{alignment_file}"
    subprocess.run(sorting_command, shell=True)
    delete_folder_with_files(temp_sorting_dir)
    alignment_file.unlink()


def index_alignments(in_dir, postfix, recursive, out_dir, processes, threads, samtools_path):
    """
    Run samtools' flagstat on all files in in_dir that end with postfix.

    @param in_dir: where the aligned and sorted files reside
    @type in_dir: Path
    @param postfix: find alignment files that end with postfix
    @type postfix: str
    @param recursive: whether to search recursively for alignments files in sub dirs of in_dir
    @type recursive: bool
    @param out_dir: write output to out_dir
    @type out_dir: Path
    @param processes: number of processes to run in parallel
    @type processes: int
    @param threads: number of threads to use in each process
    @type threads: int
    @param samtools_path: path to a samtools' executable
    @type samtools_path: str or Path
    """
    alignment_files = find_files(in_dir, postfix, recursive)
    starmap_inputs = [(alignment_file, out_dir, threads, samtools_path) for alignment_file in alignment_files]
    with Pool(processes=processes) as pool:
        pool.starmap(func=index_one_alignment, iterable=starmap_inputs)


def index_one_alignment(alignment_file, out_dir, threads, samtools_path):
    """
    Run samtools' flagstat on alignment_file.

    @param alignment_file: sorted alignment file
    @type alignment_file: Path
    @param out_dir: write output to out_dir
    @type out_dir: Path
    @param threads:
    @type threads:
    @param samtools_path:
    @type samtools_path:
    """
    index_file = Path(out_dir, f"{alignment_file.name}.bai")
    index_cmd = f"{samtools_path} index -@ {threads} {alignment_file} {index_file}"
    subprocess.run(index_cmd, shell=True)


def base_alignment_stats(in_dir, postfix, recursive, out_dir, processes, threads, samtools_path):
    """
    Run samtools' flagstat on all files in in_dir that end with postfix, and write the outputs to out_dir.

    @param in_dir: folder with sam/bam alignment files
    @type in_dir: str or Path
    @param out_dir: folder to which the stats will be written
    @type out_dir: str or Path
    @param postfix: postfix of required files in in_dir
    @type postfix: str
    @param recursive: whether to search recursively for alignments files in sub dirs of in_dir
    @type recursive: bool
    @param processes: number of processes to run in parallel
    @type processes: int
    @param threads: number of threads to use in each process
    @type threads: int
    @param samtools_path: path to samtools' executable
    @type samtools_path: str or Path
    """
    # set paths
    in_dir = Path(in_dir)
    out_dir = Path(out_dir)
    alignment_files = find_files(in_dir, postfix, recursive)
    # run samtools' flagstat on each alignment_file
    starmap_inputs = [(alignment_file, out_dir, threads, samtools_path) for alignment_file in alignment_files]
    with Pool(processes=processes) as pool:
        pool.starmap(func=one_file_base_alignment_stats, iterable=starmap_inputs)


def one_file_base_alignment_stats(alignment_file, out_dir, threads, samtools_path):
    """
    Run samtools' flagstat on in_file and write its output to out_dir.

    @param alignment_file: sam/bam alignment file
    @type alignment_file: Path
    @param out_dir: the dir in which the alignment stats file will be written
    @type out_dir: Path
    @param threads: number of threads to use
    @type threads: int
    @param samtools_path: path to a samtools' executable
    @type samtools_path: str or Path
    """
    out_file = Path(out_dir, f"{alignment_file.name}.stats")
    flagstat_cmd = f"{samtools_path} flagstat -@ {threads} {alignment_file} > {out_file}"
    subprocess.run(flagstat_cmd, shell=True)
