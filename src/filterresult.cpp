#include <stdlib.h>
#include "filterresult.h"
#include "stats.h"
#include "htmlreporter.h"

FilterResult::FilterResult(Options* opt, bool paired){
    mOptions = opt;
    mPaired = paired;
    mTrimmedAdapterRead = 0;
    mTrimmedAdapterBases = 0;
    for(int i=0; i<FILTER_RESULT_TYPES; i++) {
        mFilterReadStats[i] = 0;
    }
}

FilterResult::~FilterResult() {
}

void FilterResult::addFilterResult(int result) {
    if(result < PASS_FILTER || result >= FILTER_RESULT_TYPES)
        return ;
    // for paired end data, both reads are filtered together
    if(mPaired)
        mFilterReadStats[result] += 2;
    else
        mFilterReadStats[result]++;
}

FilterResult* FilterResult::merge(vector<FilterResult*>& list) {
    if(list.size() == 0)
        return NULL;
    FilterResult* result = new FilterResult(list[0]->mOptions, list[0]->mPaired);

    long* target = result->getFilterReadStats();
    
    for(int i=0; i<list.size(); i++) {
        long* current = list[i]->getFilterReadStats();
        for(int j=0; j<FILTER_RESULT_TYPES; j++) {
            target[j] += current[j];
        }
        result->mTrimmedAdapterRead += list[i]->mTrimmedAdapterRead;
        result->mTrimmedAdapterBases += list[i]->mTrimmedAdapterBases;

        // merge adapter stats
        map<string, long>::iterator iter;
        for(iter = list[i]->mAdapter1.begin(); iter != list[i]->mAdapter1.end(); iter++) {
            if(result->mAdapter1.count(iter->first) > 0)
                result->mAdapter1[iter->first] += iter->second;
            else
                result->mAdapter1[iter->first] = iter->second;
        }
        for(iter = list[i]->mAdapter2.begin(); iter != list[i]->mAdapter2.end(); iter++) {
            if(result->mAdapter2.count(iter->first) > 0)
                result->mAdapter2[iter->first] += iter->second;
            else
                result->mAdapter2[iter->first] = iter->second;
        }
    }

    // sort adapters list by adapter length from short to long

    return result;
}

void FilterResult::addAdapterTrimmed(string adapter) {
    if(adapter.empty())
        return;
    mTrimmedAdapterRead++;
    mTrimmedAdapterBases += adapter.length();
    if(mAdapter1.count(adapter) >0 )
        mAdapter1[adapter]++;
    else
        mAdapter1[adapter] = 1;
}

void FilterResult::addAdapterTrimmed(string adapter1, string adapter2) {
    // paired
    mTrimmedAdapterRead += 2;
    mTrimmedAdapterBases += adapter1.length() + adapter2.length();
    if(!adapter1.empty()){
        if(mAdapter1.count(adapter1) >0 )
            mAdapter1[adapter1]++;
        else
            mAdapter1[adapter1] = 1;
    }
    if(!adapter2.empty()) {
        if(mAdapter2.count(adapter2) >0 )
            mAdapter2[adapter2]++;
        else
            mAdapter2[adapter2] = 1;
    }
}

void FilterResult::print() {
    cout <<  "reads passed filter: " << mFilterReadStats[PASS_FILTER] << endl;
    cout <<  "reads failed due to low quality: " << mFilterReadStats[FAIL_QUALITY] << endl;
    cout <<  "reads failed due to too many N: " << mFilterReadStats[FAIL_N_BASE] << endl;
    cout <<  "reads failed due to too short: " << mFilterReadStats[FAIL_LENGTH] << endl;
    cout <<  "reads with adapter trimmed: " << mTrimmedAdapterRead << endl;
    cout <<  "bases trimmed due to adapters: " << mTrimmedAdapterBases << endl;
}

void FilterResult::reportJson(ofstream& ofs, string padding) {
    ofs << "{" << endl;

    ofs << padding << "\t" << "\"passed_filter_reads\": " << mFilterReadStats[PASS_FILTER] << "," << endl;
    ofs << padding << "\t" << "\"low_quality_reads\": " << mFilterReadStats[FAIL_QUALITY] << "," << endl;
    ofs << padding << "\t" << "\"too_many_N_reads\": " << mFilterReadStats[FAIL_N_BASE] << "," << endl;
    ofs << padding << "\t" << "\"too_short_reads\": " << mFilterReadStats[FAIL_LENGTH] << endl;

    ofs << padding << "}," << endl;
}

void FilterResult::outputAdaptersJson(ofstream& ofs, map<string, long, classcomp>& adapterCounts) {
    map<string, long>::iterator iter;

    long total = 0;
    for(iter = adapterCounts.begin(); iter!=adapterCounts.end(); iter++) {
        total += iter->second;
    }

    if(total == 0)
        return ;

    const double reportThreshold = 0.01;
    const double dTotal = (double)total;
    bool firstItem = true;
    long reported = 0;
    for(iter = adapterCounts.begin(); iter!=adapterCounts.end(); iter++) {
        if(iter->second /dTotal < reportThreshold )
            continue;

        if(!firstItem)
            ofs << ", ";
        else
            firstItem = false;
        ofs << "\"" << iter->first << "\":" << iter->second;

        reported += iter->second;
    }

    long unreported = total - reported;

    if(unreported > 0) {
        if(!firstItem)
            ofs << ", ";
        ofs << "\"" << "others" << "\":" << unreported;
    }
}

void FilterResult::reportAdapterJson(ofstream& ofs, string padding) {
    ofs << "{" << endl;

    ofs << padding << "\t" << "\"adapter_trimmed_reads\": " << mTrimmedAdapterRead << "," << endl;
    ofs << padding << "\t" << "\"adapter_trimmed_bases\": " << mTrimmedAdapterBases << "," << endl;

    ofs << padding << "\t" << "\"read1_adapter_counts\": " << "{";
        outputAdaptersJson(ofs, mAdapter1);
    ofs << "}";
    if(mOptions->isPaired())
        ofs << ",";
    ofs << endl;

    if(mOptions->isPaired()) {
        ofs << padding << "\t" << "\"read2_adapter_counts\": " << "{";
            outputAdaptersJson(ofs, mAdapter2);
        ofs << "}" << endl;
    }

    ofs << padding << "}," << endl;
}

/*void FilterResult::reportHtml(ofstream& ofs, long totalReads) {
    const int types = 4;
    const string divName = "filtering_result";
    string labels[4] = {"good_reads", "low_quality_reads", "too_many_N_reads", "too_short_reads"};
    long counts[4] = {mFilterReadStats[PASS_FILTER], mFilterReadStats[FAIL_QUALITY], mFilterReadStats[FAIL_N_BASE], mFilterReadStats[FAIL_LENGTH]};
    
    string json_str = "var data=[";
    json_str += "{values:[" + Stats::list2string(counts, types) + "],";
    json_str += "labels:['good_reads', 'low_quality_reads', 'too_many_N_reads', 'too_short_reads'],";
    json_str += "textinfo: 'none',";
    json_str += "type:'pie'}];\n";
    string title = "Filtering statistics of sampled " + to_string(totalReads) + " reads";
    json_str += "var layout={title:'" + title + "', width:800, height:400};\n";
    json_str += "Plotly.newPlot('" + divName + "', data, layout);\n";

    ofs << "<div class='figure' id='" + divName + "'></div>\n";
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << json_str;
    ofs << "</script>" << endl;
} */

void FilterResult::reportHtml(ofstream& ofs, long totalReads) {
    double total = (double)totalReads;
    ofs << "<table class='summary_table'>\n";
    HtmlReporter::outputRow(ofs, "reads passed filters:", HtmlReporter::formatNumber(mFilterReadStats[PASS_FILTER]) + " (" + to_string(mFilterReadStats[PASS_FILTER] * 100.0 / total) + "%)");
    HtmlReporter::outputRow(ofs, "reads with low quality:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_QUALITY]) + " (" + to_string(mFilterReadStats[FAIL_QUALITY] * 100.0 / total) + "%)");
    HtmlReporter::outputRow(ofs, "reads with too many N:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_N_BASE]) + " (" + to_string(mFilterReadStats[FAIL_N_BASE] * 100.0 / total) + "%)");
    HtmlReporter::outputRow(ofs, "reads too short:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_LENGTH]) + " (" + to_string(mFilterReadStats[FAIL_LENGTH] * 100.0 / total) + "%)");
    ofs << "</table>\n";
}

void FilterResult::reportAdapterHtml(ofstream& ofs) {
    ofs << "<div class='subsection_title'>Adapter or bad ligation of read1</div>\n";
    ofs << "<table class='summary_table'>\n";
    outputAdaptersHtml(ofs, mAdapter1);
    ofs << "</table>\n";
    if(mOptions->isPaired()) {
        ofs << "<div class='subsection_title'>Adapter or bad ligation of read2</div>\n";
        ofs << "<table class='summary_table'>\n";
        outputAdaptersHtml(ofs, mAdapter2);
        ofs << "</table>\n";
    }
}



void FilterResult::outputAdaptersHtml(ofstream& ofs, map<string, long, classcomp>& adapterCounts) {
    map<string, long>::iterator iter;

    long total = 0;
    for(iter = adapterCounts.begin(); iter!=adapterCounts.end(); iter++) {
        total += iter->second;
    }

    if(total == 0)
        return ;

    ofs << "<tr><td class='adapter_col' style='font-size:14px;color:#ffffff;background:#556699'>" << "Sequence" << "</td><td class='col2' style='font-size:14px;color:#ffffff;background:#556699'>" << "Count" << "</td></tr>\n";

    const double reportThreshold = 0.01;
    const double dTotal = (double)total;
    long reported = 0;
    for(iter = adapterCounts.begin(); iter!=adapterCounts.end(); iter++) {
        if(iter->second /dTotal < reportThreshold )
            continue;

        ofs << "<tr><td class='adapter_col'>" << iter->first << "</td><td class='col2'>" << iter->second << "</td></tr>\n";

        reported += iter->second;
    }

    long unreported = total - reported;

    if(unreported > 0) {
        ofs << "<tr><td class='adapter_col'>" << "others" << "</td><td class='col2'>" << unreported << "</td></tr>\n";
    }
}