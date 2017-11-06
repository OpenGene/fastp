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
    }
    return result;
}

void FilterResult::addAdapterTrimmed(string adapter) {
    mTrimmedAdapterRead++;
    mTrimmedAdapterBases += adapter.length();
}

void FilterResult::addAdapterTrimmed(string adapter1, string adapter2) {
    // paired
    mTrimmedAdapterRead += 2;
    mTrimmedAdapterBases += adapter1.length() + adapter2.length();
}

void FilterResult::print() {
    cout << (mPaired?"Read pairs":"Reads") << " passed filter: " << mFilterReadStats[PASS_FILTER] << endl;
    cout << (mPaired?"Read pairs":"Reads") << " failed due to low quality: " << mFilterReadStats[FAIL_QUALITY] << endl;
    cout << (mPaired?"Read pairs":"Reads") << " failed due to too many N: " << mFilterReadStats[FAIL_N_BASE] << endl;
    cout << (mPaired?"Read pairs":"Reads") << " failed due to too short: " << mFilterReadStats[FAIL_LENGTH] << endl;
    cout << (mPaired?"Read pairs":"Reads") << " with adapter trimmed: " << mTrimmedAdapterRead << endl;
    cout << "Bases" << " trimmed for adapters: " << mTrimmedAdapterBases << endl;
}

void FilterResult::reportJson(ofstream& ofs, string padding) {
    ofs << "{" << endl;

    ofs << padding << "\t" << "\"passed_filter_reads\": " << mFilterReadStats[PASS_FILTER] << "," << endl;
    ofs << padding << "\t" << "\"low_quality_reads\": " << mFilterReadStats[FAIL_QUALITY] << "," << endl;
    ofs << padding << "\t" << "\"too_many_N_reads\": " << mFilterReadStats[FAIL_N_BASE] << "," << endl;
    ofs << padding << "\t" << "\"too_short_reads\": " << mFilterReadStats[FAIL_LENGTH] << "," << endl;
    ofs << padding << "\t" << "\"adapter_trimmed_reads\": " << mTrimmedAdapterRead << "," << endl;
    ofs << padding << "\t" << "\"adapter_trimmed_bases\": " << mTrimmedAdapterBases << endl;

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