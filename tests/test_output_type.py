from make_prg.subcommands.output_type import OutputType, UnknownOutputTypeError

import pytest


class TestOutputType:
    def test_case_insensitive(self):
        value = "P"
        otype = OutputType(value)

        assert otype.prg

    def test_unknown_raises_error(self):
        value = "x"

        with pytest.raises(UnknownOutputTypeError):
            OutputType(value)

    def test_multiple_types(self):
        value = "pb"
        otype = OutputType(value)

        assert otype.prg
        assert otype.binary
        assert not otype.gfa

    def test_multiple_with_unknown_and_known_ignores_unknown(self):
        value = "gqqqq"
        otype = OutputType(value)

        assert otype.gfa

    def test_all(self):
        value = "a"
        otype = OutputType(value)

        assert otype.prg
        assert otype.binary
        assert otype.gfa

    def test_all_and_single_ignores_single(self):
        value = "ap"
        otype = OutputType(value)

        assert otype.prg
        assert otype.binary
        assert otype.gfa
