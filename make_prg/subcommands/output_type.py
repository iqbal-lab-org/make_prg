class UnknownOutputTypeError(Exception):
    pass


class OutputType:
    BINARY = "b"
    ALL = "a"
    PRG = "p"
    GFA = "g"

    def __init__(self, value: str):
        self.type = set(value.lower())

        if not self._any():
            raise UnknownOutputTypeError(f"{value} is an unknown output type")

    def _any(self) -> bool:
        return self._all() or self.prg or self.binary or self.gfa

    def _all(self) -> bool:
        return self.ALL in self.type

    @property
    def prg(self) -> bool:
        return self._all() or self.PRG in self.type

    @property
    def binary(self) -> bool:
        return self._all() or self.BINARY in self.type

    @property
    def gfa(self) -> bool:
        return self._all() or self.GFA in self.type
