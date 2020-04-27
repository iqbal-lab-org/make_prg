from typing import Dict, List, BinaryIO


class ConversionError(Exception):
    pass


class EncodeError(Exception):
    pass


PRG_Ints = List[int]

BYTES_PER_INT = 4
ENDIANNESS = "little"


class PrgEncoder:
    """
    A class that converts a prg string as produced by this program into an integer vector.
    Public interface:
        Production: run `run()`
        Serialisation: run `write()`
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
        byteorder: str = ENDIANNESS,
        num_bytes: int = BYTES_PER_INT,
    ):
        for integer in encoding:
            ostream.write(integer.to_bytes(num_bytes, byteorder=byteorder))

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
            output = [int(unit)]
        else:
            raise EncodeError("Unit {} contains invalid characters".format(unit))

        return output
