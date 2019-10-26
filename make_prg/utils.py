import logging
import re
from typing import Generator, Sequence


def remove_duplicates(seqs: Sequence) -> Generator:
    seen = set()
    for x in seqs:
        if x in seen:
            continue
        seen.add(x)
        yield x


def remove_gaps(sequence: str) -> str:
    return sequence.replace("-", "")


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


class Integer_Encoder:
    """
    A class that converts a prg string as produced by this program into an integer vector.
    Public interface:
        Production: run `run()`
        Serialisation: run `write()`
    """

    DNA = {"A": 1, "C": 2, "G": 3, "T": 4}
    NUM_BYTES = 4  # How many bytes to use per serialised integer?

    def __init__(self, prg_string):
        self.marker_units = prg_string.split()
        self.vector = []

    def run(self):
        for unit in self.marker_units:
            self.vector.extend(self._encode_unit(unit))

    def write(self, fpath):
        with open(fpath, "wb") as f:
            for integer in self.vector:
                f.write(integer.to_bytes(self.NUM_BYTES, "little"))  # Little endian

    def _DNA_to_int(self, input_char):
        assert input_char in self.DNA, logging.error(
            f"Conversion error: char {input_char} is not in {self.DNA}"
        )
        return self.DNA[input_char]

    def _encode_unit(self, input_string):
        if input_string[0] in self.DNA:
            output = list(map(self._DNA_to_int, input_string))
        else:
            output = [int(input_string)]
        return output


def write_prg(outf_prefix, prg_string):
    """
    Writes the prg to outfile.
    Writes it as a human readable string, and also as an integer vector
    """
    out_name = f"{outf_prefix}.prg"
    with open(out_name, "w") as f:
        f.write(prg_string)

    out_name = f"{outf_prefix}.bin"
    int_enc = Integer_Encoder(prg_string)
    int_enc.run()
    # print(int_enc.vector.elements)
    int_enc.write(out_name)
