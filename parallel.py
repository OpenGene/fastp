#!/usr/bin/env python

# This script is used to process FASTQ files in a folder in parallel.
# It uses the fastp/fastplong command to preprocess the FASTQ files.
# It can also generate a summary HTML report of the QC metrics.

import os,sys
from optparse import OptionParser
import time
from multiprocessing import Process, Queue
import copy
import subprocess
from concurrent.futures import ThreadPoolExecutor
import json

FASTP_PY_VERSION = "0.0.1"

def parseCommand():
    usage = "A python script to use fastp/fastplong to preprocess all FASTQ files within a folder"
    parser = OptionParser(usage = usage, version = FASTP_PY_VERSION)
    parser.add_option("-i", "--input_dir", dest = "input_dir", default = ".",
        help = "the folder contains the FASTQ files to be preprocessed, by default is current dir (.)")
    parser.add_option("-o", "--out_dir", dest = "out_dir", default = None,
        help = "the folder to store the clean FASTQ. If not specified, then there will be no output files.")
    parser.add_option("-r", "--report_dir", dest = "report_dir", default = None,
        help = "the folder to store QC reports. If not specified, use out_dir if out_dir is specified, otherwise use input_dir.")
    parser.add_option("-c", "--command", dest = "command", default = None,
        help = "the path to fastp/fastplong command, if not specified, then it will use 'fastp' in PATH")
    parser.add_option("-a", "--args", dest = "args", default = None,
        help = "the arguments that will be passed to fastp. Enclose in quotation marks. Like --args='-f 3 -t 3' ")
    parser.add_option("-p", "--parallel", dest = "parallel", default = None, type = "int",
        help = "the number of fastp processes can be run in parallel, if not specified, then it will be CPU_Core/4")
    parser.add_option("-1", "--read1_flag", dest = "read1_flag", default = "R1",
        help = "specify the name flag of read1, default is R1, which means a file with name *R1* is read1 file")
    parser.add_option("-2", "--read2_flag", dest = "read2_flag", default = "R2",
        help = "specify the name flag of read2, default is R2, which means a file with name *R2* is read2 file")
    return parser.parse_args()

def matchFlag(filename, flag):
    if flag.endswith('.') or flag.endswith('_') or flag.endswith('-'):
        return flag in filename
    else:
        return (flag+"." in filename) or (flag+"_" in filename) or (flag+"-" in filename)
    
def getBaseName(filename):
    fqext = (".fq.gz", ".fastq.gz", ".fq", ".fastq")
    for ext in fqext:
        if filename.endswith(ext):
            return filename[:-len(ext)]

def processDir(folder, options):
    fqext = (".fq", ".fastq", ".fq.gz", ".fastq.gz")
    read1name = options.read1_flag
    read2name = options.read2_flag
    
    #is not a dir
    if not os.path.isdir(folder):
        return
        
    options_list = []
    processed =  set()  # to avoid processing the same file multiple times
    
    files = os.listdir(folder)
    for f in files:
        path = os.path.join(folder, f)
        if os.path.isdir(path):
            continue
        
        isfq = False
        for ext in fqext:
            if f.endswith(ext):
                isfq = True
        if isfq == False:
            continue

        if processed.__contains__(path):
            continue

        # skip read2 files here, we will find read1 files first
        if matchFlag(f, read2name):
            continue

        processed.add(path)

        # here we skip those files with name starting with Undetermined
        # because these files are usually with unknown barcode and have no need to be processed
        if f.startswith("Undetermined"):
            continue
        
        #find read1 file, try to find read2 file
        if matchFlag(f, read1name):
            opt = copy.copy(options)
            read1 = path
            opt.read1_file = read1
            read2 = read1.replace(read1name, read2name)
            if os.path.exists(read2):
                opt.read2_file = read2
                processed.add(read2)
            options_list.append(opt)
        else:
            # if not read1 file, then run fastp as single-end
            opt = copy.copy(options)
            opt.read1_file = path
            options_list.append(opt)

    commands = []
    for opt in options_list:
        cmd = ""
        if opt.command:
            cmd = opt.command
            if not os.path.exists(cmd):
                print(f"Error: {cmd} not found, please specify the correct path to fastp/fastplong with -c option")
                sys.exit(1)
        else:
            cmd = "fastp"
        
        cmd = cmd + " -i " + opt.read1_file
        if hasattr(opt, 'read2_file'):
            cmd += " -I " + opt.read2_file
        if opt.out_dir:
            if not os.path.exists(opt.out_dir):
                os.makedirs(opt.out_dir)
            out_prefix1 = os.path.join(opt.out_dir, os.path.basename(getBaseName(opt.read1_file)))
            cmd += " -o " + out_prefix1 + ".clean.fastq.gz"
            if hasattr(opt, 'read2_file'):
                out_prefix2 = os.path.join(opt.out_dir, os.path.basename(getBaseName(opt.read2_file)))
                cmd += " -O " + out_prefix2 + ".clean.fastq.gz"
        
        if opt.args:
            cmd += " " + opt.args

        if opt.report_dir:
            if not os.path.exists(opt.report_dir):
                os.makedirs(opt.report_dir)
        
        report_file = os.path.join(opt.report_dir, os.path.basename(opt.read1_file).replace(opt.read1_flag, "report"))
        cmd += " --html=" + report_file + ".html --json=" + report_file + ".json"
        
        commands.append(cmd)
    
    if len(options_list) == 0:
        print("No FASTQ file found, do you call the program correctly?")
        print("See -h for help")
        return

    if options.parallel is None:
        options.parallel = max(1, os.cpu_count() // 4)

    with ThreadPoolExecutor(max_workers=opt.parallel) as executor:
        futures = [executor.submit(run_command, cmd) for cmd in commands]
        
        for future in futures:
            print(future.result())

def run_command(command):
    print("Running command: " + command)
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return result.stdout
    
def generate_summary_html(report_dir):
    # Collect all JSON report files
    json_files = [f for f in os.listdir(report_dir) if f.endswith('.json')]
    stats = []
    for jf in json_files:
        path = os.path.join(report_dir, jf)
        with open(path) as f:
            data = json.load(f)
            summary = data.get('summary', {})
            # Extract before and after stats, fallback to 0 if missing
            before = summary.get('before_filtering', {})
            after = summary.get('after_filtering', {})
            stat = {
                'file': jf.replace('.json', ''),
                'total_reads_before': before.get('total_reads', 0),
                'total_reads_after': after.get('total_reads', 0),
                'total_bases_before': before.get('total_bases', 0),
                'total_bases_after': after.get('total_bases', 0),
                'q20_rate_before': before.get('q20_rate', 0) * 100,
                'q20_rate_after': after.get('q20_rate', 0) * 100,
                'q30_rate_before': before.get('q30_rate', 0) * 100,
                'q30_rate_after': after.get('q30_rate', 0) * 100,
                'gc_content_before': before.get('gc_content', 0) * 100,
                'gc_content_after': after.get('gc_content', 0) * 100,
                'html_report': jf.replace('.json', '.html')
            }
            stats.append(stat)
    # Generate HTML
    html = '''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>FASTQ Summary Report</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/modern-normalize/2.0.0/modern-normalize.min.css">
    <style>
        body { font-family: 'Segoe UI', Arial, sans-serif; background: #f8f9fa; margin: 0; padding: 2em; }
        h1 { color: #2c3e50; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 2em; background: #fff; }
        th, td { border: 1px solid #e1e4e8; padding: 0.75em 1em; text-align: center; }
        th { background: #f3f6fa; color: #34495e; }
        tr:nth-child(even) { background: #f9fafb; }
        a { color: #2980b9; text-decoration: none; }
        a:hover { text-decoration: underline; }
        .chart-container { width: 100%; max-width: 900px; margin: 2em auto; background: #fff; padding: 2em; border-radius: 8px; box-shadow: 0 2px 8px #0001; }
    </style>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
</head>
<body>
    <h1>FASTQ Aggregate Summary</h1>
    <div class="chart-container">
        <canvas id="readsChart"></canvas>
    </div>
    <div class="chart-container">
        <canvas id="basesChart"></canvas>
    </div>
    <div class="chart-container">
        <canvas id="q20Chart"></canvas>
    </div>
    <div class="chart-container">
        <canvas id="q30Chart"></canvas>
    </div>
    <div class="chart-container">
        <canvas id="gcChart"></canvas>
    </div>
    <table>
        <thead>
            <tr>
                <th>File</th>
                <th>Total Reads (Before)</th>
                <th>Total Reads (After)</th>
                <th>Total Bases (Before)</th>
                <th>Total Bases (After)</th>
                <th>Q20 Rate (Before)</th>
                <th>Q20 Rate (After)</th>
                <th>Q30 Rate (Before)</th>
                <th>Q30 Rate (After)</th>
                <th>GC Content (Before)</th>
                <th>GC Content (After)</th>
                <th>HTML Report</th>
            </tr>
        </thead>
        <tbody>
'''
    for s in stats:
        html += f'<tr>'
        html += f'<td>{s["file"]}</td>'
        html += f'<td>{s["total_reads_before"]:,}</td>'
        html += f'<td>{s["total_reads_after"]:,}</td>'
        html += f'<td>{s["total_bases_before"]:,}</td>'
        html += f'<td>{s["total_bases_after"]:,}</td>'
        html += f'<td>{s["q20_rate_before"]:.2f}</td>'
        html += f'<td>{s["q20_rate_after"]:.2f}</td>'
        html += f'<td>{s["q30_rate_before"]:.2f}</td>'
        html += f'<td>{s["q30_rate_after"]:.2f}</td>'
        html += f'<td>{s["gc_content_before"]:.2f}</td>'
        html += f'<td>{s["gc_content_after"]:.2f}</td>'
        html += f'<td><a href="{s["html_report"]}">View</a></td>'
        html += '</tr>'
    html += '''
        </tbody>
    </table>
    <script>
        const files = ''' + json.dumps([s['file'] for s in stats]) + ''';
        const totalReadsBefore = ''' + json.dumps([s['total_reads_before'] for s in stats]) + ''';
        const totalReadsAfter = ''' + json.dumps([s['total_reads_after'] for s in stats]) + ''';
        const totalBasesBefore = ''' + json.dumps([s['total_bases_before'] for s in stats]) + ''';
        const totalBasesAfter = ''' + json.dumps([s['total_bases_after'] for s in stats]) + ''';
        const q20Before = ''' + json.dumps([s['q20_rate_before'] for s in stats]) + ''';
        const q20After = ''' + json.dumps([s['q20_rate_after'] for s in stats]) + ''';
        const q30Before = ''' + json.dumps([s['q30_rate_before'] for s in stats]) + ''';
        const q30After = ''' + json.dumps([s['q30_rate_after'] for s in stats]) + ''';
        const gcBefore = ''' + json.dumps([s['gc_content_before'] for s in stats]) + ''';
        const gcAfter = ''' + json.dumps([s['gc_content_after'] for s in stats]) + ''';
        // Reads chart (grouped bar)
        new Chart(document.getElementById('readsChart'), {
            type: 'bar',
            data: {
                labels: files,
                datasets: [
                    { label: 'Total Reads (Before)', data: totalReadsBefore, backgroundColor: '#95a5a6' },
                    { label: 'Total Reads (After)', data: totalReadsAfter, backgroundColor: '#3498db' }
                ]
            },
            options: {
                responsive: true,
                plugins: { legend: { position: 'top' } },
                scales: { y: { beginAtZero: true } }
            }
        });
        // Bases chart (grouped bar)
        new Chart(document.getElementById('basesChart'), {
            type: 'bar',
            data: {
                labels: files,
                datasets: [
                    { label: 'Total Bases (Before)', data: totalBasesBefore, backgroundColor: '#dfe6e9' },
                    { label: 'Total Bases (After)', data: totalBasesAfter, backgroundColor: '#636e72' }
                ]
            },
            options: {
                responsive: true,
                plugins: { legend: { position: 'top' } },
                scales: { y: { beginAtZero: true } }
            }
        });
        // Q20 chart (grouped bar)
        new Chart(document.getElementById('q20Chart'), {
            type: 'bar',
            data: {
                labels: files,
                datasets: [
                    { label: 'Q20 Rate (Before)', data: q20Before, backgroundColor: '#fab1a0' },
                    { label: 'Q20 Rate (After)', data: q20After, backgroundColor: '#e17055' }
                ]
            },
            options: {
                responsive: true,
                plugins: { legend: { position: 'top' } },
                scales: { y: { beginAtZero: true, max: 100 } }
            }
        });
        // Q30 chart (grouped bar)
        new Chart(document.getElementById('q30Chart'), {
            type: 'bar',
            data: {
                labels: files,
                datasets: [
                    { label: 'Q30 Rate (Before)', data: q30Before, backgroundColor: '#b2bec3' },
                    { label: 'Q30 Rate (After)', data: q30After, backgroundColor: '#2ecc71' }
                ]
            },
            options: {
                responsive: true,
                plugins: { legend: { position: 'top' } },
                scales: { y: { beginAtZero: true, max: 100 } }
            }
        });
        // GC chart (grouped bar)
        new Chart(document.getElementById('gcChart'), {
            type: 'bar',
            data: {
                labels: files,
                datasets: [
                    { label: 'GC Content (Before)', data: gcBefore, backgroundColor: '#ffeaa7' },
                    { label: 'GC Content (After)', data: gcAfter, backgroundColor: '#e67e22' }
                ]
            },
            options: {
                responsive: true,
                plugins: { legend: { position: 'top' } },
                scales: { y: { beginAtZero: true, max: 100 } }
            }
        });
    </script>
</body>
</html>
'''
    with open(os.path.join(report_dir, 'overall.html'), 'w') as f:
        f.write(html)

def main():
    time1 = time.time()
    
    (options, args) = parseCommand()
    options.version = FASTP_PY_VERSION

    if options.input_dir == None:
        options.input_dir="."

    if options.report_dir == None:
        if options.out_dir:
            options.report_dir = options.out_dir
        else:
            # if out_dir is not specified, use input_dir as report_dir
            options.report_dir = options.input_dir
    
    processDir(options.input_dir, options)
    # After processing, generate summary
    if options.report_dir:
        generate_summary_html(options.report_dir)
    time2 = time.time()
    print('Time used: ' + str(time2-time1))
    
if __name__  == "__main__":
    main()
