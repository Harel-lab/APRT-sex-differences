import subprocess


def find_files(in_dir, postfix, recursive):
    if recursive:
        fastq_files = list(in_dir.glob(f"**/*{postfix}"))
    else:
        fastq_files = list(in_dir.glob(f"*{postfix}"))
    return fastq_files


def group_pe_fastq_files(fastq_files, postfix, mate_prefix):
    """
    Sort files in fastq_files into tuples of paired-ended files.

    :param fastq_files: all
    :type fastq_files: list[Path]
    :param postfix:
    :type postfix: str
    :param mate_prefix:
    :type mate_prefix: str
    :return: paired_fastq_files
    :rtype: list[tuple[Path]]
    """
    sorted_fastq_files = sorted(fastq_files, key=lambda fastq_file: fastq_file.name.rstrip(postfix).split(mate_prefix))
    paired_fastq_files = [(sorted_fastq_files[x], sorted_fastq_files[x + 1])
                          for x in range(0, len(sorted_fastq_files), 2)]
    return paired_fastq_files


def decompress_file(in_file, out_file, decompress_cmd="lzcat"):
    """
    Decompress a single file using any command that decompresses the file to stdout, such as `lzcat` or
    `gunzip --stdout`.
    """
    final_decompress_cmd = f"{decompress_cmd} {in_file} > {out_file}"
    subprocess.run(final_decompress_cmd, shell=True)


def copy_text_file(src_file, dest_file):
    """
    Copy text from src_file to dest_file.

    @param src_file: source text file
    @type src_file: Path
    @param dest_file: target text file
    @type dest_file: Path
    """
    dest_file.write_text(src_file.read_text())


def delete_folders_with_files(folders):
    """
    Delete each folder and its files.

    @param folders: folders with files
    @type folders: iterable
    """
    for folder in folders:
        delete_folder_with_files(folder)


def delete_folder_with_files(folder):
    """
    Delete folder and its files.

    @param folder: folder with files to delete
    @type folder: str or Path
    """
    subprocess.run(f"rm -rf {folder}", shell=True)
