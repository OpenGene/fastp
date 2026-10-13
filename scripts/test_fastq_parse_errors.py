#!/usr/bin/env python3
"""Run: python3 scripts/test_fastq_parse_errors.py ./fastp."""
import gzip
import pathlib
import subprocess
import sys
import tempfile


def main():
    exe = str(pathlib.Path(sys.argv[1]).resolve())
    record = b'@valid\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n'
    # Place errors beyond initial sampling as well as at the first record.
    cases = [
        ('valid', record*4, None),
        ('no_final_newline', (record*4).rstrip(b'\n'), None),
        ('invalid_separator', b'@bad\nACGT\nwrong\nIIII\n', 'Invalid FASTQ separator'),
        ('missing_separator', b'@bad\nACGT\n', 'Invalid FASTQ separator'),
        ('short_quality', b'@bad\nACGT\n+\nIII\n', 'FASTQ sequence/quality length mismatch'),
        ('long_quality', b'@bad\nACGT\n+\nIIIII\n', 'FASTQ sequence/quality length mismatch'),
        ('missing_quality', b'@bad\nACGT\n+\n', 'FASTQ sequence/quality length mismatch'),
    ]
    failures = 0
    with tempfile.TemporaryDirectory(prefix='fastp-parse-test-') as tmp:
        root = pathlib.Path(tmp)
        for zipped in (False, True):
            for name, data, expected_error in cases:
                for late in ((False, True) if expected_error else (False,)):
                    label = '%s-%s-%s' % (name, zipped, late)
                    content = record*40000 + data if late else data
                    source = root/(label + ('.fq.gz' if zipped else '.fq'))
                    output = root/(label+'.out.fq')
                    source.write_bytes(gzip.compress(content) if zipped else content)
                    cmd = [exe, '-i', str(source), '-o', str(output), '-w', '2',
                           '-A', '-Q', '-L', '-G', '--dont_eval_duplication',
                           '-j', '/dev/null', '-h', '/dev/null']
                    try:
                        run = subprocess.run(cmd, capture_output=True, timeout=30)
                        if expected_error:
                            ok = (run.returncode > 0 and
                                  expected_error.encode() in run.stderr and
                                  str(source).encode() in run.stderr)
                        else:
                            ok = (run.returncode == 0 and output.exists() and
                                  output.read_bytes() == content.rstrip(b'\n')+b'\n')
                        detail = 'exit=%d' % run.returncode
                    except subprocess.TimeoutExpired:
                        ok, detail = False, 'timeout'
                    print('%s %s %s' % ('PASS' if ok else 'FAIL', label, detail))
                    failures += not ok
    return bool(failures)


if __name__ == '__main__':
    sys.exit(main())
