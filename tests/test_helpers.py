from typing import List
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from make_prg import MSA
from make_prg.prg_builder import PrgBuilderZipDatabase
import os
import shutil
import filecmp
from pathlib import Path
from zipfile import ZipFile


def make_alignment(seqs: List[str], ids: List[str] = None) -> MSA:
    seq_lengths = set(map(len, seqs))
    assert (
        len(seq_lengths) == 1
    ), "Sequences are not the same length, does not represent an alignment"
    if ids is None:
        seqrecords = [SeqRecord(Seq(seq), id=f"s{i}", description="") for i, seq in enumerate(seqs)]
    else:
        seqrecords = [SeqRecord(Seq(seq), id=ID, description="") for seq, ID in zip(seqs, ids)]
    return MSA(seqrecords)


sample_prg = " 5 AATAGGCCG 7  9 GATGCAGTTCAA 10 GATGCGGCGTA 9 AACGCCTTATCCGGCATACGA 11 ATTTATT 12 TTTTATT 11  8  13 G " \
             "14 A 13 ATGCGGCGTACGAATTTAT 15 T 16 C 15  7 CGGCCTGGCTCCCCGTAGGCCG 17 A 18 G 17 " \
             "ATAAGATGCGCCAGCATCGCATCCGGCTATAATGC 19 G 20 A 19  6  21 TTCATTGG 22 TTCAATG 23 G 24 A 23  21 " \
             "TTTATAATGCCTGATAAACGCACGGTCGATCCCCTCGCCCCTTCGGGGAGAGGATTAGGGTGAGGGGGTACAAGCCAGCCAGAGACCAGGCAA 25 " \
             "TGACATG 26 CGACATG 25  5 CACATAACC 27 TCT 28 ACC 27 TGAAACT 29  31 CTTT 32 CGTC 32 CTTC 32 CATC 31 " \
             "CCCAGAGCCTCTT 33 CAGC 34 TAGC 34 CAGT 33 CATCTATT 35 CA 36 TG 35  30  37 ACATCTCTTCA 38 ACATTTCTTCA 37  " \
             "29 GGAGCAAACAATTTCAT 39 GCCAACTC 40 TCCAACTC 40 TCCAACTT 40 ACCAACTC 39 ATAACCCCAGCATATAAATCCAG 41 T 42 " \
             "A 41 TGGTAACTTTT 43 A 44 C 43 TTTAACCT 45 G 46 A 45 AAACCAGTTT 47 TATCCAC 49 T 50 C 49  48 AATCCACC 47 " \
             "ATTTATAAAATTATGTGAAGCATTTCATAGAAGAAAAATCACTGGC 51 C 52 T 51 TAAACATTAT 53 C 54 T 53 CCCTTTTTGC 55 CTGG " \
             "56 CTGA 56 ATGA 56 CTAG 56 CTTA 56 CTGT 55 TTTTTGACCATTTCCG 57 C 58 T 57 " \
             "GATTTGTTACACATTGAAATATCACTTTTGCTGTGCGTAATATGGCTATTCGTTAGC 59 C 60 A 59 AAAAAATAAGAAAAGAT 61 T 62 A 61 "

def first_dict_contained_in_second(dict_1, dict_2) -> bool:
    return dict_1.items() <= dict_2.items()

def remove_dir_if_exists(directory):
    if os.path.exists(directory) and os.path.isdir(directory):
        shutil.rmtree(directory)


# Adapted from https://stackoverflow.com/a/6681395/5264075
def are_dir_trees_equal(dir1, dir2):
    """
    Compare two directories recursively. Files in each directory are
    assumed to be equal if their names and contents are equal.
    Zip files are "unzipped" and the contents are compared byte by byte,
    to avoid mismatching metadata issues.

    @param dir1: First directory path
    @param dir2: Second directory path

    @return: True if the directory trees are the same and
        there were no errors while accessing the directories or files,
        False otherwise.
   """
    # compare directory listings
    dirs_cmp = filecmp.dircmp(dir1, dir2)
    if len(dirs_cmp.left_only)>0 or len(dirs_cmp.right_only)>0 or \
        len(dirs_cmp.funny_files)>0:
        return False

    # compare non-zip files
    common_files_without_zip = list(filter(lambda filename: not filename.endswith(".zip"), dirs_cmp.common_files))
    (_, mismatch, errors) = filecmp.cmpfiles(dir1, dir2, common_files_without_zip, shallow=False)
    if len(mismatch) > 0:
        print(f"[Dir comparison] File mismatches: {mismatch}")
        return False
    if len(errors)>0:
        print(f"[Dir comparison] Errors: {errors}")
        return False

    # compare zip files
    common_files_with_zip = list(filter(lambda filename: filename.endswith(".zip"), dirs_cmp.common_files))
    for file in common_files_with_zip:
        zip_file_1 = os.path.join(dir1, file)
        zip_file_2 = os.path.join(dir2, file)
        if not are_zip_files_equal(zip_file_1, zip_file_2):
            print(f"[Dir comparison] Zip file mismatch: {zip_file_1} and {zip_file_2}")
            return False

    # recursively compare subdirs
    for common_dir in dirs_cmp.common_dirs:
        new_dir1 = os.path.join(dir1, common_dir)
        new_dir2 = os.path.join(dir2, common_dir)
        if not are_dir_trees_equal(new_dir1, new_dir2):
            return False

    return True


def are_zip_files_equal(file_1: [Path, str], file_2: [Path, str]) -> bool:
    with ZipFile(file_1) as zip_file_1, ZipFile(file_2) as zip_file_2:
        if zip_file_1.namelist() != zip_file_2.namelist():
            print(f"Error: {file_1} and {file_2} differ due to namelist: ")
            print(f"File 1 namelist: {' '.join(zip_file_1.namelist())}")
            print(f"File 2 namelist: {' '.join(zip_file_2.namelist())}")
            return False

        is_update_DS_zip = str(file_1).endswith(".update_DS.zip") and str(file_2).endswith(".update_DS.zip")
        if is_update_DS_zip:
            try:
                zip_db_1 = PrgBuilderZipDatabase(Path(file_1))
                zip_db_1.load()
                zip_db_2 = PrgBuilderZipDatabase(Path(file_2))
                zip_db_2.load()
                return zip_db_1 == zip_db_2
            finally:
                zip_db_1.close()
                zip_db_2.close()
        else:
            for file in zip_file_1.namelist():
                bytes_from_zip_file_1 = zip_file_1.read(file)
                bytes_from_zip_file_2 = zip_file_2.read(file)
                if bytes_from_zip_file_1 != bytes_from_zip_file_2:
                    print(f"Error: file {file} differs in bytes when comparing {file_1} and {file_2}")
                    return False

        return True
