#!/usr/bin/env python

# This script is used to process FASTQ files in a folder in parallel.
# It uses the fastp command to preprocess the FASTQ files.
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
    usage = "A python script to use fastp to preprocess all FASTQ files within a folder"
    parser = OptionParser(usage = usage, version = FASTP_PY_VERSION)
    parser.add_option("-i", "--input_dir", dest = "input_dir", default = ".",
        help = "the folder contains the FASTQ files to be preprocessed, by default is current dir (.)")
    parser.add_option("-o", "--out_dir", dest = "out_dir", default = None,
        help = "the folder to store the clean FASTQ. If not specified, then there will be no output files.")
    parser.add_option("-r", "--report_dir", dest = "report_dir", default = None,
        help = "the folder to store QC reports. If not specified, use out_dir if out_dir is specified, otherwise use input_dir.")
    parser.add_option("-c", "--command", dest = "command", default = None,
        help = "the path to fastp command, if not specified, then it will use 'fastp' in PATH")
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
                print(f"Error: {cmd} not found, please specify the correct path to fastp with -c option")
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
        
        report_file = os.path.join(opt.report_dir, os.path.basename(opt.read1_file).replace(opt.read1_flag, "pe"))
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
    
def generate_summary_html(report_dir, fastp_cmd=None):
    # Get fastp version from JSON files
    fastp_version = "fastp"
    if fastp_cmd:
        fastp_version = fastp_cmd

    json_files = [f for f in os.listdir(report_dir) if f.endswith('.json')]
    if json_files:
        try:
            # Use the first JSON file to get the fastp version
            first_json = os.path.join(report_dir, json_files[0])
            with open(first_json) as f:
                data = json.load(f)
                fastp_version = fastp_version + " "+ data.get('summary', {}).get('fastp_version', 'unknown')
        except:
            pass
    
    # Collect all JSON report files
    json_files = [f for f in os.listdir(report_dir) if f.endswith('.json')]
    stats = []
    mean_qual_curves = []
    max_len_before = 0
    max_len_after = 0
    gc_curves = []
    max_len_gc_before = 0
    max_len_gc_after = 0
    mean_qual_curves_r2 = []
    max_len_before_r2 = 0
    max_len_after_r2 = 0
    gc_curves_r2 = []
    max_len_gc_before_r2 = 0
    max_len_gc_after_r2 = 0
    for jf in json_files:
        path = os.path.join(report_dir, jf)
        with open(path) as f:
            data = json.load(f)
            summary = data.get('summary', {})
            before = summary.get('before_filtering', {})
            after = summary.get('after_filtering', {})
            # Extract quality and GC curves for read1
            qual_curve_before = data.get('read1_before_filtering', {}).get('quality_curves', {}).get('mean', [])
            qual_curve_after = data.get('read1_after_filtering', {}).get('quality_curves', {}).get('mean', [])
            gc_curve_before = data.get('read1_before_filtering', {}).get('content_curves', {}).get('GC', [])
            gc_curve_after = data.get('read1_after_filtering', {}).get('content_curves', {}).get('GC', [])
            mean_qual_curves.append({
                'file': jf.replace('.json', ''),
                'curve_before': qual_curve_before,
                'curve_after': qual_curve_after
            })
            gc_curves.append({
                'file': jf.replace('.json', ''),
                'curve_before': gc_curve_before,
                'curve_after': gc_curve_after
            })
            if len(qual_curve_before) > max_len_before:
                max_len_before = len(qual_curve_before)
            if len(qual_curve_after) > max_len_after:
                max_len_after = len(qual_curve_after)
            if len(gc_curve_before) > max_len_gc_before:
                max_len_gc_before = len(gc_curve_before)
            if len(gc_curve_after) > max_len_gc_after:
                max_len_gc_after = len(gc_curve_after)
            # Extract quality and GC curves for read2 if present
            qual_curve_before_r2 = data.get('read2_before_filtering', {}).get('quality_curves', {}).get('mean', [])
            qual_curve_after_r2 = data.get('read2_after_filtering', {}).get('quality_curves', {}).get('mean', [])
            gc_curve_before_r2 = data.get('read2_before_filtering', {}).get('content_curves', {}).get('GC', [])
            gc_curve_after_r2 = data.get('read2_after_filtering', {}).get('content_curves', {}).get('GC', [])
            if qual_curve_before_r2 or qual_curve_after_r2 or gc_curve_before_r2 or gc_curve_after_r2:
                mean_qual_curves_r2.append({
                    'file': jf.replace('.json', ''),
                    'curve_before': qual_curve_before_r2,
                    'curve_after': qual_curve_after_r2
                })
                gc_curves_r2.append({
                    'file': jf.replace('.json', ''),
                    'curve_before': gc_curve_before_r2,
                    'curve_after': gc_curve_after_r2
                })
                if len(qual_curve_before_r2) > max_len_before_r2:
                    max_len_before_r2 = len(qual_curve_before_r2)
                if len(qual_curve_after_r2) > max_len_after_r2:
                    max_len_after_r2 = len(qual_curve_after_r2)
                if len(gc_curve_before_r2) > max_len_gc_before_r2:
                    max_len_gc_before_r2 = len(gc_curve_before_r2)
                if len(gc_curve_after_r2) > max_len_gc_after_r2:
                    max_len_gc_after_r2 = len(gc_curve_after_r2)
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
        .chart-container { width: 100%; max-width: none; aspect-ratio: 4/1; }
        .row-charts-table { width: 100%; margin-bottom: 2em; background: none; border: none; }
        .row-charts-table td { border: none; vertical-align: top; width: 50%; }
    </style>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
    <h2>FASTQ Aggregate Summary (''' + fastp_version + ''')</h2>
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
    def human_format(num):
        if num >= 1e9:
            return f"{num/1e9:.2f}G"
        elif num >= 1e6:
            return f"{num/1e6:.2f}M"
        elif num >= 1e3:
            return f"{num/1e3:.2f}K"
        else:
            return str(num)
    for s in stats:
        html += f'<tr>'
        html += f'<td>{s["file"]}</td>'
        html += f'<td>{human_format(s["total_reads_before"])}</td>'
        html += f'<td>{human_format(s["total_reads_after"])}</td>'
        html += f'<td>{human_format(s["total_bases_before"])}</td>'
        html += f'<td>{human_format(s["total_bases_after"])}</td>'
        html += f'<td>{s["q20_rate_before"]:.2f}%</td>'
        html += f'<td>{s["q20_rate_after"]:.2f}%</td>'
        html += f'<td>{s["q30_rate_before"]:.2f}%</td>'
        html += f'<td>{s["q30_rate_after"]:.2f}%</td>'
        html += f'<td>{s["gc_content_before"]:.2f}%</td>'
        html += f'<td>{s["gc_content_after"]:.2f}%</td>'
        html += f'<td><a href="{s["html_report"]}">View</a></td>'
        html += '</tr>'
    html += '''
        </tbody>
    </table>
    <table class="row-charts-table">
        <tr>
            <td><div id="meanQualPlotBefore" style="width:100%;height:400px;"></div></td>
            <td><div id="meanQualPlotAfter" style="width:100%;height:400px;"></div></td>
        </tr>
'''
    if mean_qual_curves_r2:
        html += '        <tr>\n'
        html += '            <td><div id="meanQualPlotBeforeR2" style="width:100%;height:400px;"></div></td>\n'
        html += '            <td><div id="meanQualPlotAfterR2" style="width:100%;height:400px;"></div></td>\n'
        html += '        </tr>\n'
    html += '        <tr>\n'
    html += '            <td><div id="gcCurvePlotBefore" style="width:100%;height:400px;"></div></td>\n'
    html += '            <td><div id="gcCurvePlotAfter" style="width:100%;height:400px;"></div></td>\n'
    html += '        </tr>\n'
    if gc_curves_r2:
        html += '        <tr>\n'
        html += '            <td><div id="gcCurvePlotBeforeR2" style="width:100%;height:400px;"></div></td>\n'
        html += '            <td><div id="gcCurvePlotAfterR2" style="width:100%;height:400px;"></div></td>\n'
        html += '        </tr>\n'
    html += '    </table>\n'
    html += '''
    <div class="chart-container" style="width:100%; max-width:none; aspect-ratio: 4/1;">
        <canvas id="qRateChart" style="height:200px;"></canvas>
    </div>
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
        // Plotly mean quality curves (before)
        const meanQualCurves = ''' + json.dumps(mean_qual_curves) + ''';
        const maxLenBefore = ''' + str(max_len_before) + ''';
        const maxLenAfter = ''' + str(max_len_after) + ''';
        const qualLabelsBefore = Array.from({length: maxLenBefore}, (_, i) => i + 1);
        const qualLabelsAfter = Array.from({length: maxLenAfter}, (_, i) => i + 1);
        const plotlyTracesBefore = meanQualCurves.map((item, idx) => {
            const before = item.curve_before || [];
            const beforePad = before.concat(Array(maxLenBefore - before.length).fill(null));
            return {
                x: qualLabelsBefore,
                y: beforePad,
                mode: 'lines',
                name: item.file,
                line: { width: 1 }
            };
        });
        Plotly.newPlot('meanQualPlotBefore', plotlyTracesBefore, {
            title: 'Mean Quality Curve (Read1, Before Filtering)',
            xaxis: { title: '' },
            yaxis: { title: 'Mean Quality', rangemode: 'tozero' },
            legend: { orientation: 'h' },
            margin: { t: 50, l: 60, r: 30, b: 60 }
        }, {responsive: true});
        // Plotly mean quality curves (after)
        const plotlyTracesAfter = meanQualCurves.map((item, idx) => {
            const after = item.curve_after || [];
            const afterPad = after.concat(Array(maxLenAfter - after.length).fill(null));
            return {
                x: qualLabelsAfter,
                y: afterPad,
                mode: 'lines',
                name: item.file,
                line: { width: 1 }
            };
        });
        Plotly.newPlot('meanQualPlotAfter', plotlyTracesAfter, {
            title: 'Mean Quality Curve (Read1, After Filtering)',
            xaxis: { title: '' },
            yaxis: { title: 'Mean Quality', rangemode: 'tozero' },
            legend: { orientation: 'h' },
            margin: { t: 50, l: 60, r: 30, b: 60 }
        }, {responsive: true});
        // Plotly GC content curves (before)
        const gcCurves = ''' + json.dumps(gc_curves) + ''';
        const maxLenGCBefore = ''' + str(max_len_gc_before) + ''';
        const maxLenGCAfter = ''' + str(max_len_gc_after) + ''';
        const gcLabelsBefore = Array.from({length: maxLenGCBefore}, (_, i) => i + 1);
        const gcLabelsAfter = Array.from({length: maxLenGCAfter}, (_, i) => i + 1);
        const plotlyTracesGCBefore = gcCurves.map((item, idx) => {
            const before = item.curve_before || [];
            const beforePad = before.concat(Array(maxLenGCBefore - before.length).fill(null));
            return {
                x: gcLabelsBefore,
                y: beforePad,
                mode: 'lines',
                name: item.file,
                line: { width: 1 }
            };
        });
        Plotly.newPlot('gcCurvePlotBefore', plotlyTracesGCBefore, {
            title: 'GC Content Curve (Read1, Before Filtering)',
            xaxis: { title: '' },
            yaxis: { title: 'GC %', rangemode: 'tozero' },
            legend: { orientation: 'h' },
            margin: { t: 50, l: 60, r: 30, b: 60 }
        }, {responsive: true});
        // Plotly GC content curves (after)
        const plotlyTracesGCAfter = gcCurves.map((item, idx) => {
            const after = item.curve_after || [];
            const afterPad = after.concat(Array(maxLenGCAfter - after.length).fill(null));
            return {
                x: gcLabelsAfter,
                y: afterPad,
                mode: 'lines',
                name: item.file,
                line: { width: 1 }
            };
        });
        Plotly.newPlot('gcCurvePlotAfter', plotlyTracesGCAfter, {
            title: 'GC Content Curve (Read1, After Filtering)',
            xaxis: { title: '' },
            yaxis: { title: 'GC %', rangemode: 'tozero' },
            legend: { orientation: 'h' },
            margin: { t: 50, l: 60, r: 30, b: 60 }
        }, {responsive: true});
        // Plotly mean quality and GC curves for read2 (if any)
        const meanQualCurvesR2 = ''' + json.dumps(mean_qual_curves_r2) + ''';
        const gcCurvesR2 = ''' + json.dumps(gc_curves_r2) + ''';
        const maxLenBeforeR2 = ''' + str(max_len_before_r2) + ''';
        const maxLenAfterR2 = ''' + str(max_len_after_r2) + ''';
        const maxLenGCBeforeR2 = ''' + str(max_len_gc_before_r2) + ''';
        const maxLenGCAfterR2 = ''' + str(max_len_gc_after_r2) + ''';
        if (meanQualCurvesR2.length > 0 || gcCurvesR2.length > 0) {
            const qualLabelsBeforeR2 = Array.from({length: maxLenBeforeR2}, (_, i) => i + 1);
            const qualLabelsAfterR2 = Array.from({length: maxLenAfterR2}, (_, i) => i + 1);
            const gcLabelsBeforeR2 = Array.from({length: maxLenGCBeforeR2}, (_, i) => i + 1);
            const gcLabelsAfterR2 = Array.from({length: maxLenGCAfterR2}, (_, i) => i + 1);
            const plotlyTracesBeforeR2 = meanQualCurvesR2.map((item, idx) => {
                const before = item.curve_before || [];
                const beforePad = before.concat(Array(maxLenBeforeR2 - before.length).fill(null));
                return {
                    x: qualLabelsBeforeR2,
                    y: beforePad,
                    mode: 'lines',
                    name: item.file,
                    line: { width: 1 }
                };
            });
            Plotly.newPlot('meanQualPlotBeforeR2', plotlyTracesBeforeR2, {
                title: 'Mean Quality Curve (Read2, Before Filtering)',
                xaxis: { title: '' },
                yaxis: { title: 'Mean Quality', rangemode: 'tozero' },
                legend: { orientation: 'h' },
                margin: { t: 50, l: 60, r: 30, b: 60 }
            }, {responsive: true});
            const plotlyTracesAfterR2 = meanQualCurvesR2.map((item, idx) => {
                const after = item.curve_after || [];
                const afterPad = after.concat(Array(maxLenAfterR2 - after.length).fill(null));
                return {
                    x: qualLabelsAfterR2,
                    y: afterPad,
                    mode: 'lines',
                    name: item.file,
                    line: { width: 1 }
                };
            });
            Plotly.newPlot('meanQualPlotAfterR2', plotlyTracesAfterR2, {
                title: 'Mean Quality Curve (Read2, After Filtering)',
                xaxis: { title: '' },
                yaxis: { title: 'Mean Quality', rangemode: 'tozero' },
                legend: { orientation: 'h' },
                margin: { t: 50, l: 60, r: 30, b: 60 }
            }, {responsive: true});
            // GC content curves for read2
            const plotlyTracesGCBeforeR2 = gcCurvesR2.map((item, idx) => {
                const before = item.curve_before || [];
                const beforePad = before.concat(Array(maxLenGCBeforeR2 - before.length).fill(null));
                return {
                    x: gcLabelsBeforeR2,
                    y: beforePad,
                    mode: 'lines',
                    name: item.file,
                    line: { width: 1 }
                };
            });
            Plotly.newPlot('gcCurvePlotBeforeR2', plotlyTracesGCBeforeR2, {
                title: 'GC Content Curve (Read2, Before Filtering)',
                xaxis: { title: '' },
                yaxis: { title: 'GC %', rangemode: 'tozero' },
                legend: { orientation: 'h' },
                margin: { t: 50, l: 60, r: 30, b: 60 }
            }, {responsive: true});
            const plotlyTracesGCAfterR2 = gcCurvesR2.map((item, idx) => {
                const after = item.curve_after || [];
                const afterPad = after.concat(Array(maxLenGCAfterR2 - after.length).fill(null));
                return {
                    x: gcLabelsAfterR2,
                    y: afterPad,
                    mode: 'lines',
                    name: item.file,
                    line: { width: 1 }
                };
            });
            Plotly.newPlot('gcCurvePlotAfterR2', plotlyTracesGCAfterR2, {
                title: 'GC Content Curve (Read2, After Filtering)',
                xaxis: { title: '' },
                yaxis: { title: 'GC %', rangemode: 'tozero' },
                legend: { orientation: 'h' },
                margin: { t: 50, l: 60, r: 30, b: 60 }
            }, {responsive: true});
        }
        // Q20/Q30/Q40 chart (grouped bar)
        new Chart(document.getElementById('qRateChart'), {
            type: 'bar',
            data: {
                labels: files,
                datasets: [
                    { label: 'Q20 Rate (Before)', data: q20Before, backgroundColor: '#fab1a0' },
                    { label: 'Q20 Rate (After)', data: q20After, backgroundColor: '#e17055' },
                    { label: 'Q30 Rate (Before)', data: q30Before, backgroundColor: '#b2bec3' },
                    { label: 'Q30 Rate (After)', data: q30After, backgroundColor: '#2ecc71' }
                ]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: { legend: { position: 'top' } },
                scales: { y: { beginAtZero: true, max: 100 } },
                animation: false
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
        generate_summary_html(options.report_dir, options.command)
    time2 = time.time()
    print('Time used: ' + str(time2-time1))
    
if __name__  == "__main__":
    main()
