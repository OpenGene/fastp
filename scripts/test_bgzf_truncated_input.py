#!/usr/bin/env python3
"""Run: python3 scripts/test_bgzf_truncated_input.py ./fastp

Synthetic BGZF fixtures; no external data or compression tools required.
"""
import gzip
import pathlib
import struct
import subprocess
import sys
import tempfile
import zlib


def block(data):
    compressor = zlib.compressobj(wbits=-15)
    payload = compressor.compress(data) + compressor.flush()
    size = 18 + len(payload) + 8
    assert size <= 65536
    return (b'\x1f\x8b\x08\x04' + b'\x00'*4 + b'\x00\xff' +
            struct.pack('<H', 6) + b'BC' + struct.pack('<HH', 2, size-1) +
            payload + struct.pack('<II', zlib.crc32(data), len(data)))


def main():
    exe = str(pathlib.Path(sys.argv[1]).resolve())
    first = b'@first\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n'
    second = first.replace(b'first', b'second')
    a, b, eof = block(first), block(second), block(b'')
    malformed_size = bytearray(b)
    malformed_size[16:18] = struct.pack('<H', 0)
    cases = [
        ('valid', a+b+eof, first+second, None),
        ('no_eof_marker', a+b, first+second, None),
        ('concatenated', a+eof+b+eof, first+second, None),
        ('partial_header', a+b[:9], None, 'Truncated BGZF block header'),
        ('partial_payload', a+b[:20], None, 'Truncated BGZF block body'),
        ('partial_trailer', a+b[:-3], None, 'Truncated BGZF block body'),
        ('undersized_block', a+bytes(malformed_size), None, 'Invalid BGZF block size'),
    ]
    failed = 0
    with tempfile.TemporaryDirectory(prefix='fastp-bgzf-test-') as tmp:
        root = pathlib.Path(tmp)
        for threads in (1, 4):
            for name, content, expected, diagnostic in cases:
                source = root/('%s-%d.fq.gz' % (name, threads))
                output = root/('%s-%d.out.fq' % (name, threads))
                source.write_bytes(content)
                if expected is not None:
                    assert gzip.decompress(content) == expected
                cmd = [exe, '-i', str(source), '-o', str(output),
                       '-w', str(threads), '-A', '-Q', '-L', '-G',
                       '--dont_eval_duplication', '-j', '/dev/null', '-h', '/dev/null']
                try:
                    run = subprocess.run(cmd, capture_output=True, timeout=30)
                    if diagnostic:
                        ok = run.returncode > 0 and diagnostic.encode() in run.stderr
                    else:
                        ok = (run.returncode == 0 and output.exists() and
                              output.read_bytes() == expected)
                    detail = 'exit=%d' % run.returncode
                except subprocess.TimeoutExpired:
                    ok, detail = False, 'timeout'
                print('%s %s threads=%d %s' %
                      ('PASS' if ok else 'FAIL', name, threads, detail))
                failed += not ok
    return bool(failed)


if __name__ == '__main__':
    sys.exit(main())
