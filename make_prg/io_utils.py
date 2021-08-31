import logging
import re
import gzip
from pathlib import Path
from argparse import ArgumentParser

from Bio import AlignIO

from make_prg.from_msa import MSA
from make_prg.prg_encoder import PrgEncoder, PRG_Ints


def load_alignment_file(msa_file: str, alignment_format: str) -> MSA:
    msa_file = str(msa_file)
    logging.info("Read from MSA file %s", msa_file)
    if msa_file.endswith(".gz"):
        logging.debug("MSA is gzipped")
        handle = gzip.open(msa_file, "rt")
        alignment = AlignIO.read(handle, alignment_format)
        handle.close()
    else:
        alignment = AlignIO.read(msa_file, alignment_format)
    for record in alignment:
        record.seq = record.seq.upper()
    return alignment


# ************/
# GFA code */
# ***********/
class GFA_Output:
    """
    A simple class for converting a PRG string into a GFA file
    TODO: Update to GFA2 format
    """

    def __init__(self, gfa_string="", gfa_id=0, gfa_site=5):
        self.gfa_string = gfa_string
        self.gfa_id = gfa_id
        self.gfa_site = gfa_site
        self.delim_char = " "  # This mirrors the AlignedSeq class.

    def split_on_site(self, prg_string, site_num):
        site_coords = [
            (a.start(), a.end())
            for a in list(
                re.finditer(
                    "%s%d%s" % (self.delim_char, site_num, self.delim_char), prg_string
                )
            )
        ]
        last_pos = None
        split_strings = []
        for (start, end) in site_coords:
            split_strings.append(prg_string[last_pos:start])
            last_pos = end
        split_strings.append(prg_string[last_pos:])
        delim = "%s%d%s" % (self.delim_char, site_num, self.delim_char)
        check_string = delim.join(split_strings)
        assert check_string == prg_string, (
            "Something has gone wrong with the string split for site %d\nsplit_"
            "strings: %s" % (site_num, split_strings)
        )
        return split_strings

    def build_gfa_string(self, prg_string, pre_var_id=None):
        """Takes prg_string and builds a gfa_string with fragments
        from the prg_string."""
        end_ids = []
        # iterate through sites present, updating gfa_string with each in turn
        while str(self.gfa_site) in prg_string:
            logging.debug("gfa_site: %d", self.gfa_site)
            prgs = self.split_on_site(prg_string, self.gfa_site)
            logging.debug("prgs: %s", prgs)
            assert len(prgs) == 3, "Invalid prg sequence %s for site %d and id %d" % (
                prg_string,
                self.gfa_site,
                self.gfa_id,
            )

            # add pre-var site string and links from previous seq fragments
            if prgs[0] != "":
                self.gfa_string += "S\t%d\t%s\tRC:i:0\n" % (self.gfa_id, prgs[0])
            else:
                # adds an empty node for empty pre var site seqs
                self.gfa_string += "S\t%d\t%s\tRC:i:0\n" % (self.gfa_id, "*")
            pre_var_id = self.gfa_id
            self.gfa_id += 1
            for id in end_ids:
                self.gfa_string += "L\t%d\t+\t%d\t+\t0M\n" % (id, pre_var_id)
                end_ids = []

            # recursively add segments for each of the variant haplotypes at
            # this site, saving the end id for each haplotype
            vars = self.split_on_site(prgs[1], self.gfa_site + 1)
            assert len(vars) > 1, "Invalid prg sequence %s for site %d and id %d" % (
                prg_string,
                self.gfa_site + 1,
                self.gfa_id,
            )
            logging.debug("vars: %s", vars)
            self.gfa_site += 2
            logging.debug("gfa_site: %d", self.gfa_site)
            for var_string in vars:
                if pre_var_id != None:
                    self.gfa_string += "L\t%d\t+\t%d\t+\t0M\n" % (
                        pre_var_id,
                        self.gfa_id,
                    )
                var_end_ids = self.build_gfa_string(
                    prg_string=var_string, pre_var_id=pre_var_id
                )
                end_ids.extend(var_end_ids)

            prg_string = prgs[2]
            pre_var_id = None

        # finally add the final bit of sequence after variant site
        if prg_string != "":
            self.gfa_string += "S\t%d\t%s\tRC:i:0\n" % (self.gfa_id, prg_string)
        else:
            self.gfa_string += "S\t%d\t%s\tRC:i:0\n" % (self.gfa_id, "*")
        for id in end_ids:
            self.gfa_string += "L\t%d\t+\t%d\t+\t0M\n" % (id, self.gfa_id)
        end_ids = []
        return_id = [self.gfa_id]
        self.gfa_id += 1
        return return_id


def write_gfa(outfile, prg_string):
    """
    Writes a gfa file from prg string.
    """
    with open(outfile, "w") as f:
        # initialize gfa_string, id and site, then update string with the prg
        gfa_string = "H\tVN:Z:1.0\tbn:Z:--linear --singlearr\n"
        gfa_id = 0
        gfa_site = 5
        gfa_obj = GFA_Output(gfa_string)
        gfa_obj.build_gfa_string(prg_string=prg_string)
        f.write(gfa_obj.gfa_string)


# ******************/
# Write PRG code */
# *****************/


def write_prg(prg_fname: Path, prg_string: str, options: ArgumentParser):
    """
    Writes th prg to `output_file`.
    Writes it as a human readable string, and also as an integer vector
    """
    seqid = options.seqid or options.prg_name

    if options.output_type.prg:
        with prg_fname.open("w") as prg:
            header = f">{seqid} max_nest={options.max_nesting} min_match={options.min_match_length}"
            print(f"{header}\n{prg_string}", file=prg)

    if options.output_type.binary:
        prg_ints_fpath = prg_fname.with_suffix(".bin")
        prg_encoder = PrgEncoder()
        prg_ints: PRG_Ints = prg_encoder.encode(prg_string)

        with prg_ints_fpath.open("wb") as ostream:
            prg_encoder.write(prg_ints, ostream)
