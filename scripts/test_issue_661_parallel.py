#!/usr/bin/env python3
"""Regression tests for paired-end filename flags in parallel.py (issue #661)."""

import contextlib
import io
from pathlib import Path
import shlex
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import parallel


class ParallelPairingTest(unittest.TestCase):
    def test_complete_suffix_flags(self):
        for extension in (".fq", ".fastq", ".fq.gz", ".fastq.gz"):
            for read in ("_1", "_2"):
                flag = read + extension
                with self.subTest(flag=flag):
                    self.assertTrue(parallel.matchFlag("sample" + flag, flag))

    def test_partial_flags_and_delimiters(self):
        for flag in ("R1", "_1", "_1.fastq"):
            for delimiter in (".", "_", "-"):
                with self.subTest(flag=flag, delimiter=delimiter):
                    self.assertTrue(parallel.matchFlag("sample" + flag + delimiter + "gz", flag))
                    self.assertTrue(parallel.matchFlag("sample" + flag + delimiter + "fq", flag + delimiter))

    def test_does_not_match_longer_read_numbers(self):
        for flag in ("R1", "_1"):
            for extension in (".fq", ".fastq", ".fq.gz", ".fastq.gz"):
                with self.subTest(flag=flag, extension=extension):
                    self.assertFalse(parallel.matchFlag("sample" + flag + "0" + extension, flag))

    def assert_command(self, filenames, flags, read1, read2=None):
        # Read flags in the parent directory must not affect mate lookup.
        with tempfile.TemporaryDirectory(prefix="fastp_R1__1_") as directory:
            folder = Path(directory)
            for filename in filenames:
                (folder / filename).touch()
            options = SimpleNamespace(
                read1_flag=flags[0], read2_flag=flags[1], command=None,
                args=None, parallel=1, out_dir=str(folder / "out"),
                report_dir=str(folder / "reports"),
            )
            with patch.object(parallel.os, "listdir", return_value=filenames), \
                    patch.object(parallel, "run_command", return_value="") as run_command, \
                    contextlib.redirect_stdout(io.StringIO()):
                parallel.processDir(str(folder), options)

            run_command.assert_called_once()
            command = shlex.split(run_command.call_args.args[0])
            expected = ["fastp", "-i", str(folder / read1)]
            if read2 is not None:
                expected += ["-I", str(folder / read2)]
            expected += ["-o", str(folder / "out" / (read1.split(".")[0] + ".clean.fastq.gz"))]
            if read2 is not None:
                expected += ["-O", str(folder / "out" / (read2.split(".")[0] + ".clean.fastq.gz"))]
            else:
                self.assertNotIn("-I", command)
                self.assertNotIn("-O", command)
            self.assertEqual(command[:len(expected)], expected)

    def test_complete_flags_pair_reads_in_either_order(self):
        for extension in (".fq", ".fastq", ".fq.gz", ".fastq.gz"):
            read1, read2 = "sample_1" + extension, "sample_2" + extension
            flags = ("_1" + extension, "_2" + extension)
            for filenames in ([read1, read2], [read2, read1]):
                with self.subTest(extension=extension, filenames=filenames):
                    self.assert_command(filenames, flags, read1, read2)

    def test_partial_flags_still_pair_reads(self):
        for flags in (("R1", "R2"), ("_1", "_2"), ("_1.fastq", "_2.fastq")):
            read1 = "sample_" + flags[0] + (".gz" if flags[0].endswith(".fastq") else ".fastq.gz")
            read2 = read1.replace(flags[0], flags[1])
            with self.subTest(flags=flags):
                self.assert_command([read2, read1], flags, read1, read2)

    def test_missing_mate_remains_single_end(self):
        for flags in (("_1.fastq.gz", "_2.fastq.gz"), ("_1", "_2")):
            with self.subTest(flags=flags):
                self.assert_command(["sample_1.fastq.gz"], flags, "sample_1.fastq.gz")


if __name__ == "__main__":
    unittest.main()
