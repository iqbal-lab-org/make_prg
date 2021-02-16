from typing import Dict, List, BinaryIO


class ConversionError(Exception):
    pass


class EncodeError(Exception):
    pass


PRG_Ints = List[int]

BYTES_PER_INT = 4
ENDIANNESS = "little"


def to_bytes(integer: int):
    return integer.to_bytes(BYTES_PER_INT, ENDIANNESS)


class PrgEncoder:
    """
    A class that converts a prg string as produced by this program into an integer vector.
    Note that gramtools uses the following encoding for an A/T SNP: 5A6T6
    Ie, requires an odd marker at the beginning of a site and an even marker at the end of a site.
    """

    encoding = {"A": 1, "C": 2, "G": 3, "T": 4}

    def __init__(self, encoding: Dict[str, int] = None):
        """An encoder to convert a PRG string produced by this program into a vector of
        integers as specified by `encoding`.

        :param encoding: Specifies which integer characters in the PRG are to be
        encoded as.
        """
        if encoding is not None:
            self.encoding = encoding
        self._seen_site_markers = set()

    def encode(self, prg: str) -> PRG_Ints:
        marker_units = prg.split()
        encoding = []
        for unit in marker_units:
            encoding.extend(self._encode_unit(unit))
        return encoding

    @staticmethod
    def write(
        encoding: List[int],
        ostream: BinaryIO,
    ):
        ostream.write(b"".join(map(to_bytes, encoding)))

    def _dna_to_int(self, input_char: str) -> int:
        input_char = input_char.upper()
        if input_char not in self.encoding:
            raise ConversionError(f"Char '{input_char}' is not in {self.encoding}")
        return self.encoding[input_char]

    def _encode_unit(self, unit: str) -> List[int]:
        is_empty_string = not unit
        if is_empty_string:
            raise EncodeError("Cannot encode an empty string")

        chars_are_valid = all(c.upper() in self.encoding for c in unit)
        if chars_are_valid:
            output = [self._dna_to_int(char) for char in unit]
        elif unit.isdigit():
            site_marker = int(unit)
            if site_marker in self._seen_site_markers:
                output = [site_marker + 1]
                self._seen_site_markers.remove(site_marker)
            else:
                self._seen_site_markers.add(site_marker)
                output = [site_marker]
        else:
            raise EncodeError("Unit {} contains invalid characters".format(unit))

        return output
