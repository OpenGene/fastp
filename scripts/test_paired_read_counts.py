"""Run: python3 scripts/test_paired_read_counts.py ./fastp (POSIX)."""
import pathlib
import subprocess
import tempfile
import sys

exe = str(pathlib.Path(sys.argv[1]).resolve())
failures = 0
with tempfile.TemporaryDirectory(prefix='fastp-pairs-') as tmp:
    root = pathlib.Path(tmp)
    for threads in (1, 4):
        for n, m, limit in ((999,999,0),(1000,1000,0),(1001,1001,0),
                            (999,1000,0),(1000,999,0),(1000,1001,0),
                            (1001,1000,0),(41000,41001,0),(41001,41000,0),
                            (1000,1001,500)):
            paths = [root/'r1.fq', root/'r2.fq']
            for mate, count in enumerate((n,m)):
                paths[mate].write_text(''.join('@r%d/%d\n%s\n+\n%s\n' %
                    (i,mate+1,'ACGT'*10,'I'*40) for i in range(count)))
            cmd = [exe,'-i',str(paths[0]),'-I',str(paths[1]),'-o',str(root/'o1.fq'),
                   '-O',str(root/'o2.fq'),'-w',str(threads),'-A','-Q','-L','-G',
                   '--dont_eval_duplication','-j','/dev/null','-h','/dev/null']
            if limit:
                cmd += ['--reads_to_process',str(limit)]
            try:
                run = subprocess.run(cmd,capture_output=True,timeout=5)
                if n != m and not limit:
                    ok = (run.returncode > 0 and
                          b'Paired-end input files contain different numbers of reads' in run.stderr)
                else:
                    expected = limit or n
                    ok = run.returncode == 0
                    for mate in (1, 2):
                        output = root/('o%d.fq' % mate)
                        ok = ok and output.exists() and output.read_bytes().count(b'\n') == expected*4
                failures += not ok
                print('PASS' if ok else 'FAIL',n,m,threads,limit,'exit',run.returncode,flush=True)
            except subprocess.TimeoutExpired:
                failures += 1
                print(n,m,threads,limit,'TIMEOUT',flush=True)
sys.exit(bool(failures))
