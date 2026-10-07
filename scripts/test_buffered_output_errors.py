#!/usr/bin/env python3
"""POSIX regression test: python3 scripts/test_buffered_output_errors.py ./fastp.

Only temporary synthetic data is used. RLIMIT_FSIZE simulates write failures
without filling a disk; SIGXFSZ is ignored so fastp must detect the I/O error.
"""
import gzip
import pathlib
import random
import resource
import signal
import subprocess
import sys
import tempfile


def restrict_output():
    signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
    resource.setrlimit(resource.RLIMIT_FSIZE, (4096, 4096))


def main():
    exe = str(pathlib.Path(sys.argv[1]).resolve())
    rng = random.Random(917)
    data = ''.join('@read%d\n%s\n+\n%s\n' %
                   (i, ''.join(rng.choices('ACGT', k=100)), 'I'*100)
                   for i in range(2000)).encode()
    failures = []
    with tempfile.TemporaryDirectory(prefix='fastp-writer-test-') as tmp:
        root = pathlib.Path(tmp)
        source = root/'input.fq'
        source.write_bytes(data)
        for mode in ('plain', 'gzip', 'stdout', 'split'):
            for limited in (False, True):
                case = root/('%s-%s' % (mode, limited))
                case.mkdir()
                output = case/('out.fq.gz' if mode == 'gzip' else 'out.fq')
                cmd = [exe, '-i', str(source), '-w', '1', '-A', '-Q', '-L',
                       '-G', '--dont_eval_duplication',
                       '-j', '/dev/null', '-h', '/dev/null']
                if mode == 'stdout':
                    cmd.append('--stdout')
                else:
                    cmd += ['-o', str(output)]
                if mode == 'split':
                    cmd += ['--split', '2']
                with (case/'stdout.fq').open('wb') as stream:
                    run = subprocess.run(cmd, stdout=stream, stderr=subprocess.PIPE,
                                         timeout=30, restore_signals=False,
                                         preexec_fn=restrict_output if limited else None)
                error = run.stderr.decode(errors='replace')
                if limited:
                    ok = (run.returncode > 0 and
                          ('Failed to write complete output to:' in error or
                           'Failed to finalise output to:' in error))
                else:
                    if mode == 'stdout':
                        actual = (case/'stdout.fq').read_bytes()
                    elif mode == 'gzip':
                        actual = gzip.decompress(output.read_bytes())
                    elif mode == 'split':
                        actual = b''.join(p.read_bytes() for p in sorted(case.glob('*.out.fq')))
                    else:
                        actual = output.read_bytes()
                    ok = run.returncode == 0 and actual == data
                print('%s %s limited=%s exit=%d' %
                      ('PASS' if ok else 'FAIL', mode, limited, run.returncode))
                if not ok:
                    failures.append((mode, limited, error))
    for mode, limited, error in failures:
        print('%s limited=%s:\n%s' % (mode, limited, error), file=sys.stderr)
    return bool(failures)


if __name__ == '__main__':
    sys.exit(main())
