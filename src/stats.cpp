#include "stats.h"
#include <memory.h>
#include <sstream>
#include "util.h"

#define KMER_LEN 5

Stats::Stats(Options* opt, bool isRead2, int guessedCycles, int bufferMargin){
    mOptions = opt;
    mIsRead2 = isRead2;
    mReads = 0;
    mLengthSum = 0;

    mEvaluatedSeqLen = mOptions->seqLen1;
    if(mIsRead2)
        mEvaluatedSeqLen = mOptions->seqLen2;

    if(guessedCycles == 0) {
        guessedCycles = mEvaluatedSeqLen;
    }

    mCycles = guessedCycles;
    mBases = 0;
    mQ20Total = 0;
    mQ30Total = 0;
    summarized = false;
    mKmerMin = 0;
    mKmerMax = 0;

    // extend the buffer to make sure it's long enough
    mBufLen = guessedCycles + bufferMargin;

    for(int i=0; i<8; i++){
        mQ20Bases[i] = 0;
        mQ30Bases[i] = 0;
        mBaseContents[i] = 0;

        mCycleQ30Bases[i] = new long[mBufLen];
        memset(mCycleQ30Bases[i], 0, sizeof(long) * mBufLen);

        mCycleQ20Bases[i] = new long[mBufLen];
        memset(mCycleQ20Bases[i], 0, sizeof(long) * mBufLen);

        mCycleBaseContents[i] = new long[mBufLen];
        memset(mCycleBaseContents[i], 0, sizeof(long) * mBufLen);

        mCycleBaseQual[i] = new long[mBufLen];
        memset(mCycleBaseQual[i], 0, sizeof(long) * mBufLen);
    }
    mCycleTotalBase = new long[mBufLen];
    memset(mCycleTotalBase, 0, sizeof(long)*mBufLen);

    mCycleTotalQual = new long[mBufLen];
    memset(mCycleTotalQual, 0, sizeof(long)*mBufLen);

    mKmerBufLen = 2<<(KMER_LEN * 2);
    mKmer = new long[mKmerBufLen];
    memset(mKmer, 0, sizeof(long)*mKmerBufLen);

    initOverRepSeq();
}

void Stats::extendBuffer(int newBufLen){
    if(newBufLen <= mBufLen)
        return ;

    long* newBuf = NULL;

    for(int i=0; i<8; i++){
        newBuf = new long[newBufLen];
        memset(newBuf, 0, sizeof(long)*newBufLen);
        memcpy(newBuf, mCycleQ30Bases[i], sizeof(long) * mBufLen);
        delete mCycleQ30Bases[i];
        mCycleQ30Bases[i] = newBuf;

        newBuf = new long[newBufLen];
        memset(newBuf, 0, sizeof(long)*newBufLen);
        memcpy(newBuf, mCycleQ20Bases[i], sizeof(long) * mBufLen);
        delete mCycleQ20Bases[i];
        mCycleQ20Bases[i] = newBuf;

        newBuf = new long[newBufLen];
        memset(newBuf, 0, sizeof(long)*newBufLen);
        memcpy(newBuf, mCycleBaseContents[i], sizeof(long) * mBufLen);
        delete mCycleBaseContents[i];
        mCycleBaseContents[i] = newBuf;

        newBuf = new long[newBufLen];
        memset(newBuf, 0, sizeof(long)*newBufLen);
        memcpy(newBuf, mCycleBaseQual[i], sizeof(long) * mBufLen);
        delete mCycleBaseQual[i];
        mCycleBaseQual[i] = newBuf;
    }
    newBuf = new long[newBufLen];
    memset(newBuf, 0, sizeof(long)*newBufLen);
    memcpy(newBuf, mCycleTotalBase, sizeof(long)*mBufLen);
    delete mCycleTotalBase;
    mCycleTotalBase = newBuf;

    newBuf = new long[newBufLen];
    memset(newBuf, 0, sizeof(long)*newBufLen);
    memcpy(newBuf, mCycleTotalQual, sizeof(long)*mBufLen);
    delete mCycleTotalQual;
    mCycleTotalQual = newBuf;

    mBufLen = newBufLen;
}

Stats::~Stats() {
    for(int i=0; i<8; i++){
        delete mCycleQ30Bases[i];
        mCycleQ30Bases[i] = NULL;

        delete mCycleQ20Bases[i];
        mCycleQ20Bases[i] = NULL;

        delete mCycleBaseContents[i];
        mCycleBaseContents[i] = NULL;

        delete mCycleBaseQual[i];
        mCycleBaseQual[i] = NULL;
    }

    delete mCycleTotalBase;
    delete mCycleTotalQual;

    // delete memory of curves
    map<string, double*>::iterator iter;
    for(iter = mQualityCurves.begin(); iter != mQualityCurves.end(); iter++) {
        delete iter->second;
    }
    for(iter = mContentCurves.begin(); iter != mContentCurves.end(); iter++) {
        delete iter->second;
    }
    delete mKmer;

    deleteOverRepSeqDist();
}

void Stats::summarize(bool forced) {
    if(summarized && !forced)
        return;

    // first get the cycle and count total bases
    for(int c=0; c<mBufLen; c++) {
        mBases += mCycleTotalBase[c];
        if (mCycleTotalBase[c] == 0){
            mCycles = c;
            break;
        }
    }
    if(mCycleTotalBase[mBufLen-1]>0)
        mCycles = mBufLen;

    // Q20, Q30, base content
    for(int i=0; i<8; i++) {
        for(int c=0; c<mCycles; c++) {
            mQ20Bases[i] += mCycleQ20Bases[i][c];
            mQ30Bases[i] += mCycleQ30Bases[i][c];
            mBaseContents[i] += mCycleBaseContents[i][c];
        }
        mQ20Total += mQ20Bases[i];
        mQ30Total += mQ30Bases[i];
    }


    // quality curve for mean qual
    double* meanQualCurve = new double[mCycles];
    memset(meanQualCurve, 0, sizeof(double)*mCycles);
    for(int c=0; c<mCycles; c++) {
        meanQualCurve[c] = (double)mCycleTotalQual[c] / (double)mCycleTotalBase[c];
    }
    mQualityCurves["mean"] = meanQualCurve;

    // quality curves and base content curves for different nucleotides
    char alphabets[5] = {'A', 'T', 'C', 'G', 'N'};
    for(int i=0; i<5; i++) {
        char base = alphabets[i];
        // get last 3 bits
        char b = base & 0x07;
        double* qualCurve = new double[mCycles];
        memset(qualCurve, 0, sizeof(double)*mCycles);
        double* contentCurve = new double[mCycles];
        memset(contentCurve, 0, sizeof(double)*mCycles);
        for(int c=0; c<mCycles; c++) {
            if(mCycleBaseContents[b][c] == 0)
                qualCurve[c] = meanQualCurve[c];
            else
                qualCurve[c] = (double)mCycleBaseQual[b][c] / (double)mCycleBaseContents[b][c];
            contentCurve[c] = (double)mCycleBaseContents[b][c] / (double)mCycleTotalBase[c];
        }
        mQualityCurves[string(1, base)] = qualCurve;
        mContentCurves[string(1, base)] = contentCurve;
    }

    // GC content curve
    double* gcContentCurve = new double[mCycles];
    memset(gcContentCurve, 0, sizeof(double)*mCycles);
    char gBase = 'G' & 0x07;
    char cBase = 'C' & 0x07;
    for(int c=0; c<mCycles; c++) {
        gcContentCurve[c] = (double)(mCycleBaseContents[gBase][c] + mCycleBaseContents[cBase][c]) / (double)mCycleTotalBase[c];
    }
    mContentCurves["GC"] = gcContentCurve;

    mKmerMin = mKmer[0];
    mKmerMax = mKmer[0];
    for(int i=0; i<mKmerBufLen; i++) {
        if(mKmer[i] > mKmerMax)
            mKmerMax = mKmer[i];
        if(mKmer[i] < mKmerMin)
            mKmerMin = mKmer[i];
    }

    summarized = true;
}

int Stats::getMeanLength() {
    if(mReads == 0)
        return 0.0;
    else
        return mLengthSum/mReads;
}

void Stats::statRead(Read* r) {
    int len = r->length();

    mLengthSum += len;

    if(mBufLen < len) {
        extendBuffer(max(len + 100, (int)(len * 1.5)));
    }
    const char* seqstr = r->mSeq.mStr.c_str();
    const char* qualstr = r->mQuality.c_str();

    int kmer = 0;
    bool needFullCompute = true;
    for(int i=0; i<len; i++) {
        char base = seqstr[i];
        char qual = qualstr[i];
        // get last 3 bits
        char b = base & 0x07;

        const char q20 = '5';
        const char q30 = '?';

        if(qual >= q30) {
            mCycleQ30Bases[b][i]++;
            mCycleQ20Bases[b][i]++;
        } else if(qual >= q20) {
            mCycleQ20Bases[b][i]++;
        }

        mCycleBaseContents[b][i]++;
        mCycleBaseQual[b][i] += (qual-33);

        mCycleTotalBase[i]++;
        mCycleTotalQual[i] += (qual-33);

        if(base == 'N'){
            needFullCompute = true;
            continue;
        }

        // 5 bases required for kmer computing
        if(i<4)
            continue;

        // calc 5 KMER
        // 0x3FC == 0011 1111 1100
        if(!needFullCompute){
            int val = base2val(base);
            if(val < 0){
                needFullCompute = true;
                continue;
            } else {
                kmer = ((kmer<<2) & 0x3FC ) | val;
                mKmer[kmer]++;
            }
        } else {
            bool valid = true;
            kmer = 0;
            for(int k=0; k<5; k++) {
                int val = base2val(seqstr[i - 4 + k]);
                if(val < 0) {
                    valid = false;
                    break;
                }
                kmer = ((kmer<<2) & 0x3FC ) | val;
            }
            if(!valid) {
                needFullCompute = true;
                continue;
            } else {
                mKmer[kmer]++;
                needFullCompute = false;
            }
        }

    }

    // do overrepresentation analysis for 1 of every 100 reads
    if(mOptions->overRepAnalysis.enabled) {
        if(mReads % mOptions->overRepAnalysis.sampling == 0) {
            const int steps[5] = {10, 20, 40, 100, min(150, mEvaluatedSeqLen-2)};
            for(int s=0; s<5; s++) {
                int step = steps[s];
                for(int i=0; i<len-step; i++) {
                    string seq = r->mSeq.mStr.substr(i, step);
                    if(mOverRepSeq.count(seq)>0) {
                        mOverRepSeq[seq]++;
                        for(int p = i; p < seq.length() + i && p < mEvaluatedSeqLen; p++) {
                            mOverRepSeqDist[seq][p]++;
                        }
                        i+=step;
                    }
                }
            }
        }
    }

    mReads++;
}

int Stats::base2val(char base) {
    switch(base){
        case 'A':
            return 0;
        case 'T':
            return 1;
        case 'C':
            return 2;
        case 'G':
            return 3;
        default:
            return -1;
    }
}

int Stats::getCycles() {
    if(!summarized)
        summarize();
    return mCycles;
}

long Stats::getReads() {
    if(!summarized)
        summarize();
    return mReads;
}

long Stats::getBases() {
    if(!summarized)
        summarize();
    return mBases;
}

long Stats::getQ20() {
    if(!summarized)
        summarize();
    return mQ20Total;
}

long Stats::getQ30() {
    if(!summarized)
        summarize();
    return mQ30Total;
}

long Stats::getGCNumber() {
    if(!summarized)
        summarize();
    return mBaseContents['G' & 0x07] + mBaseContents['C' & 0x07];
}

void Stats::print() {
    if(!summarized) {
        summarize();
    }
    cerr << "total reads: " << mReads << endl;
    cerr << "total bases: " << mBases << endl;
    cerr << "Q20 bases: " << mQ20Total << "(" << (mQ20Total*100.0)/mBases << "%)" << endl;
    cerr << "Q30 bases: " << mQ30Total << "(" << (mQ30Total*100.0)/mBases << "%)" << endl;
}

void Stats::reportJson(ofstream& ofs, string padding) {
    ofs << "{" << endl;

    ofs << padding << "\t" << "\"total_reads\": " << mReads << "," << endl;
    ofs << padding << "\t" << "\"total_bases\": " << mBases << "," << endl;
    ofs << padding << "\t" << "\"q20_bases\": " << mQ20Total << "," << endl;
    ofs << padding << "\t" << "\"q30_bases\": " << mQ30Total << "," << endl;
    ofs << padding << "\t" << "\"total_cycles\": " << mCycles << "," << endl;

    // quality curves
    string qualNames[5] = {"A", "T", "C", "G", "mean"};
    ofs << padding << "\t" << "\"quality_curves\": {" << endl;
    for(int i=0 ;i<5; i++) {
        string name=qualNames[i];
        double* curve = mQualityCurves[name];
        ofs << padding << "\t\t" << "\"" << name << "\":[";
        for(int c = 0; c<mCycles; c++) {
            ofs << curve[c];
            // not the end
            if(c != mCycles - 1)
                ofs << ",";
        }
        ofs << "]";
        // not the end;
        if(i != 5-1)
            ofs << ",";
        ofs << endl; 
    }
    ofs << padding << "\t" << "}," << endl;

    // content curves
    string contentNames[6] = {"A", "T", "C", "G", "N", "GC"};
    ofs << padding << "\t" << "\"content_curves\": {" << endl;
    for(int i=0 ;i<6; i++) {
        string name=contentNames[i];
        double* curve = mContentCurves[name];
        ofs << padding << "\t\t" << "\"" << name << "\":[";
        for(int c = 0; c<mCycles; c++) {
            ofs << curve[c];
            // not the end
            if(c != mCycles - 1)
                ofs << ",";
        }
        ofs << "]";
        // not the end;
        if(i != 6-1)
            ofs << ",";
        ofs << endl; 
    }
    ofs << padding << "\t" << "}," << endl;

    // KMER counting
    ofs << padding << "\t" << "\"kmer_count\": {" << endl;
    for(int i=0; i<64; i++) {
        string first = kmer3(i);
        for(int j=0; j<16; j++) {
            int target = (i<<4) + j;
            long count = mKmer[target];
            string last = kmer2(j);
            ofs << padding << "\t\t\"" << first << last << "\":" << count;
            if(j != 16-1)
                ofs << ",";
        }
        if(i != 64-1)
            ofs << "," << endl;
        else
            ofs << endl;
    }
    ofs << padding << "\t" << "}," << endl;

    // over represented seqs
    map<string, long>::iterator iter;
    bool first = true;
    ofs << padding << "\t" << "\"overrepresented_sequences\": {" << endl;
    for(iter=mOverRepSeq.begin(); iter!=mOverRepSeq.end(); iter++) {
        string seq = iter->first;
        long count = iter->second;
        if(!overRepPassed(seq, count))
            continue;
        if(!first) {
            ofs << "," << endl;
        } else
            first = false;
        ofs << padding << "\t\t\"" << seq <<  "\":" << count;
    }
    ofs << padding << "\t" << "}" << endl;

    ofs << padding << "}," << endl;
}

string Stats::list2string(double* list, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << list[i];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

string Stats::list2string(double* list, int size, long* coords) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        // coords is 1,2,3,...
        long start = 0;
        if(i>0)
            start = coords[i-1];
        long end = coords[i];

        double total = 0.0;
        for(int k=start; k<end; k++)
            total += list[k];

        // get average
        if(end == start)
            ss << "0.0";
        else
            ss << total / (end - start);
        //ss << list[coords[i]-1];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

string Stats::list2string(long* list, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << list[i];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

void Stats::reportHtml(ofstream& ofs, string filteringType, string readName) {
    reportHtmlQuality(ofs, filteringType, readName);
    reportHtmlContents(ofs, filteringType, readName);
    reportHtmlKMER(ofs, filteringType, readName);
    if(mOptions->overRepAnalysis.enabled) {
        reportHtmlORA(ofs, filteringType, readName);
    }
}

bool Stats::overRepPassed(string& seq, long count) {
    int s = mOptions->overRepAnalysis.sampling;
    switch(seq.length()) {
        case 10:
            return s * count > 500;
        case 20:
            return s * count > 200;
        case 40:
            return s * count > 100;
        case 100:
            return s * count > 50;
        default:
            return s * count > 20;
    }
}

void Stats::reportHtmlORA(ofstream& ofs, string filteringType, string readName) {
    // over represented seqs
    double dBases = mBases;
    map<string, long>::iterator iter;
    int displayed = 0;

    // KMER
    string subsection = filteringType + ": " + readName + ": overrepresented sequences";
    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    string title = "";

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div  id='" << divName << "'>\n";
    ofs << "<div class='sub_section_tips'>Sampling rate: 1 / " << mOptions->overRepAnalysis.sampling << "</div>\n";
    ofs << "<table class='summary_table'>\n";
    ofs << "<tr style='font-weight:bold;'><td>overrepresented sequence</td><td>count (% of bases)</td><td>distribution: cycle 1 ~ cycle " << mEvaluatedSeqLen << "</td></tr>"<<endl;
    int found = 0;
    for(iter=mOverRepSeq.begin(); iter!=mOverRepSeq.end(); iter++) {
        string seq = iter->first;
        long count = iter->second;
        if(!overRepPassed(seq, count))
            continue;
        found++;
        double percent = (100.0 * count * seq.length() * mOptions->overRepAnalysis.sampling)/dBases;
        ofs << "<tr>";
        ofs << "<td width='400' style='word-break:break-all;font-size:8px;'>" << seq << "</td>";
        ofs << "<td width='200'>" << count << " (" << to_string(percent) <<"%)</td>";
        ofs << "<td width='250'><canvas id='" << divName << "_" << seq << "' width='240' height='20'></td>";
        ofs << "</tr>" << endl;
    }
    if(found == 0)
        ofs << "<tr><td style='text-align:center' colspan='3'>not found</td></tr>" << endl;
    ofs << "</table>\n";
    ofs << "</div>\n";

    // output the JS
    ofs << "<script language='javascript'>" << endl;
    ofs << "var seqlen = " << mEvaluatedSeqLen << ";" << endl;
    ofs << "var orp_dist = {" << endl;
    bool first = true;
    for(iter=mOverRepSeq.begin(); iter!=mOverRepSeq.end(); iter++) {
        string seq = iter->first;
        long count = iter->second;
        if(!overRepPassed(seq, count))
            continue;

        if(!first) {
            ofs << "," << endl;
        } else
            first = false;
        ofs << "\t\"" << divName << "_" << seq << "\":[";
        for(int i=0; i<mEvaluatedSeqLen; i++){
            if(i !=0 )
                ofs << ",";
            ofs << mOverRepSeqDist[seq][i];
        }
        ofs << "]";
    }
    ofs << "\n};" << endl;

    ofs << "for (seq in orp_dist) {"<< endl;
    ofs << "    var cvs = document.getElementById(seq);"<< endl;
    ofs << "    var ctx = cvs.getContext('2d'); "<< endl;
    ofs << "    var data = orp_dist[seq];"<< endl;
    ofs << "    var w = 240;"<< endl;
    ofs << "    var h = 20;"<< endl;
    ofs << "    ctx.fillStyle='#cccccc';"<< endl;
    ofs << "    ctx.fillRect(0, 0, w, h);"<< endl;
    ofs << "    ctx.fillStyle='#0000FF';"<< endl;
    ofs << "    var maxVal = 0;"<< endl;
    ofs << "    for(d=0; d<seqlen; d++) {"<< endl;
    ofs << "        if(data[d]>maxVal) maxVal = data[d];"<< endl;
    ofs << "    }"<< endl;
    ofs << "    var step = (seqlen-1) /  (w-1);"<< endl;
    ofs << "    for(x=0; x<w; x++){"<< endl;
    ofs << "        var target = step * x;"<< endl;
    ofs << "        var val = data[Math.floor(target)];"<< endl;
    ofs << "        var y = Math.floor((val / maxVal) * h);"<< endl;
    ofs << "        ctx.fillRect(x,h-1, 1, -y);"<< endl;
    ofs << "    }"<< endl;
    ofs << "}"<< endl;
    ofs << "</script>"<< endl;
}

bool Stats::isLongRead() {
    return mCycles > 300;
}

void Stats::reportHtmlKMER(ofstream& ofs, string filteringType, string readName) {

    // KMER
    string subsection = filteringType + ": " + readName + ": KMER counting";
    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    string title = "";

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div  id='" << divName << "'>\n";
    ofs << "<div class='sub_section_tips'>Darker background means larger counts. The count will be shown on mouse over.</div>\n";
    ofs << "<table class='kmer_table' style='width:680px;'>\n";
    ofs << "<tr>";
    ofs << "<td></td>";
    // the heading row
    for(int h=0; h<16; h++) 
        ofs << "<td style='color:#333333'>" << kmer2(h) << "</td>";
    ofs << "</tr>\n";
    // content
    for(int i=0; i<64; i++) {
        ofs << "<tr>";

        ofs << "<td style='color:#333333'>" << kmer3(i) << "</td>";
        for(int j=0; j<16; j++) {
            ofs << makeKmerTD(i,j) ;
        }
        ofs << "</tr>\n";
    }
    ofs << "</table>\n";
    ofs << "</div>\n";
}

string Stats::makeKmerTD(int i, int j) {
    int target = (i<<4) + j;
    long val = mKmer[target];
    // 3bp + 2bp = 5bp
    string first = kmer3(i);
    string last = kmer2(j);
    string kmer = first+last;
    double meanBases = (double)(mBases+1) / mKmerBufLen;
    double prop = val / meanBases;
    double frac = 0.5;
    if(prop > 2.0) 
        frac = (prop-2.0)/20.0 + 0.5;
    else if(prop< 0.5)
        frac = prop;

    frac = max(0.01, min(1.0, frac));
    int r = (1.0-frac) * 255;
    int g = r;
    int b = r;
    stringstream ss;
    ss << "<td style='background:#"; 
    if(r<16)
        ss << "0";
    ss<<hex<<r;
    if(g<16)
        ss << "0";
    ss<<hex<<g;
    if(b<16)
        ss << "0";
    ss<<hex<<b;
    ss << dec << "' title='"<< kmer << ": " << val << "\n" << prop << " times as mean value'>";
    ss << kmer << "</td>";
    return ss.str();
}

string Stats::kmer3(int val) {
    const char bases[4] = {'A', 'T', 'C', 'G'};
    string ret(3, ' ');
    ret[0] = bases[(val & 0x30) >> 4];
    ret[1] = bases[(val & 0x0C) >> 2];
    ret[2] = bases[(val & 0x03)];
    return ret;
}

string Stats::kmer2(int val) {
    const char bases[4] = {'A', 'T', 'C', 'G'};
    string ret(2, ' ');
    ret[0] = bases[(val & 0x0C) >> 2];
    ret[1] = bases[(val & 0x03)];
    return ret;
}

void Stats::reportHtmlQuality(ofstream& ofs, string filteringType, string readName) {

    // quality
    string subsection = filteringType + ": " + readName + ": quality";
    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    string title = "";

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";
    
    string alphabets[5] = {"A", "T", "C", "G", "mean"};
    string colors[5] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    long *x = new long[mCycles];
    int total = 0;
    if(!isLongRead()) {
        for(int i=0; i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
    } else {
        const int fullSampling = 40;
        for(int i=0; i<fullSampling && i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
        // down sampling if it's too long
        if(mCycles>fullSampling) {
            double pos = fullSampling;
            while(true){
                pos *= 1.05;
                if(pos >= mCycles)
                    break;
                x[total] = (int)pos;
                total++;
            }
            // make sure lsat one is contained
            if(x[total-1] != mCycles){
                x[total] = mCycles;
                total++;
            }
        }
    }
    // four bases
    for (int b = 0; b<5; b++) {
        string base = alphabets[b];
        json_str += "{";
        json_str += "x:[" + list2string(x, total) + "],";
        json_str += "y:[" + list2string(mQualityCurves[base], total, x) + "],";
        json_str += "name: '" + base + "',";
        json_str += "mode:'lines',";
        json_str += "line:{color:'" + colors[b] + "', width:1}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'position'";
    // use log plot if it's too long
    if(isLongRead()) {
        json_str += ",type:'log'";
    }
    json_str += "}, yaxis:{title:'quality'}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
}

void Stats::reportHtmlContents(ofstream& ofs, string filteringType, string readName) {

    // content
    string subsection = filteringType + ": " + readName + ": base contents";
    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    string title = "";

    ofs << "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('" << divName << "')>" + subsection + "</a></div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";
    
    string alphabets[6] = {"A", "T", "C", "G", "N", "GC"};
    string colors[6] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(255, 0, 0, 1.0)", "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    long *x = new long[mCycles];
    int total = 0;
    if(!isLongRead()) {
        for(int i=0; i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
    } else {
        const int fullSampling = 40;
        for(int i=0; i<fullSampling && i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
        // down sampling if it's too long
        if(mCycles>fullSampling) {
            double pos = fullSampling;
            while(true){
                pos *= 1.05;
                if(pos >= mCycles)
                    break;
                x[total] = (int)pos;
                total++;
            }
            // make sure lsat one is contained
            if(x[total-1] != mCycles){
                x[total] = mCycles;
                total++;
            }
        }
    }
    // four bases
    for (int b = 0; b<6; b++) {
        string base = alphabets[b];
        long count = 0;
        if(base.size()==1) {
            char b = base[0] & 0x07;
            count = mBaseContents[b];
        } else {
            count = mBaseContents['G' & 0x07] + mBaseContents['C' & 0x07] ;
        }
        string percentage = to_string((double)count * 100.0 / mBases);
        if(percentage.length()>5)
            percentage = percentage.substr(0,5);
        string name = base + "(" + percentage + "%)"; 

        json_str += "{";
        json_str += "x:[" + list2string(x, total) + "],";
        json_str += "y:[" + list2string(mContentCurves[base], total, x) + "],";
        json_str += "name: '" + name + "',";
        json_str += "mode:'lines',";
        json_str += "line:{color:'" + colors[b] + "', width:1}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'position'";
    // use log plot if it's too long
    if(isLongRead()) {
        json_str += ",type:'log'";
    }
    json_str += "}, yaxis:{title:'base content ratios'}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
}

Stats* Stats::merge(vector<Stats*>& list) {
    if(list.size() == 0)
        return NULL;

    //get the most long cycles
    int cycles = 0;
    for(int t=0; t<list.size(); t++) {
        list[t]->summarize();
        cycles = max(cycles, list[t]->getCycles());
    }

    Stats* s = new Stats(list[0]->mOptions, list[0]->mIsRead2, cycles, 0);

    // init overrepresented seq maps
    map<string, long>::iterator iter;

    for(int t=0; t<list.size(); t++) {
        int curCycles =  list[t]->getCycles();
        // merge read number
        s->mReads += list[t]->mReads;
        s->mLengthSum += list[t]->mLengthSum;

        // merge per cycle counting for different bases
        for(int i=0; i<8; i++){
            for(int j=0; j<cycles && j<curCycles; j++) {
                s->mCycleQ30Bases[i][j] += list[t]->mCycleQ30Bases[i][j];
                s->mCycleQ20Bases[i][j] += list[t]->mCycleQ20Bases[i][j];
                s->mCycleBaseContents[i][j] += list[t]->mCycleBaseContents[i][j];
                s->mCycleBaseQual[i][j] += list[t]->mCycleBaseQual[i][j];
            }
        }

        // merge per cycle counting for all bases
        for(int j=0; j<cycles && j<curCycles; j++) {
            s->mCycleTotalBase[j] += list[t]->mCycleTotalBase[j];
            s->mCycleTotalQual[j] += list[t]->mCycleTotalQual[j];
        }

        // merge kMer
        for(int i=0; i<s->mKmerBufLen; i++) {
            s->mKmer[i] += list[t]->mKmer[i];
        }

        // merge over rep seq
        for(iter = s->mOverRepSeq.begin(); iter != s->mOverRepSeq.end(); iter++) {
            string seq = iter->first;
            s->mOverRepSeq[seq] += list[t]->mOverRepSeq[seq];
            if(s->mIsRead2 != list[t]->mIsRead2 || list[t]->mOverRepSeqDist[seq] == NULL)
                cerr << t <<seq<< ":" << (s->mIsRead2?2:1 ) << "," << (list[t]->mIsRead2?2:1 ) <<endl;
            for(int i=0; i<s->mEvaluatedSeqLen; i++) {
                s->mOverRepSeqDist[seq][i] += list[t]->mOverRepSeqDist[seq][i];
            }
        }
    }

    s->summarize();

    return s;
}

void Stats::initOverRepSeq() {
    map<string, long> overRepSeq;
    if(mIsRead2)
        overRepSeq = mOptions->overRepSeqs2;
    else
        overRepSeq = mOptions->overRepSeqs1;

    map<string, long>::iterator iter;
    for(iter = overRepSeq.begin(); iter!=overRepSeq.end(); iter++) {
        string seq = iter->first;
        mOverRepSeq[seq] = 0;
        long* distBuf = new long[mEvaluatedSeqLen];
        memset(distBuf, 0, sizeof(long)*mEvaluatedSeqLen);
        mOverRepSeqDist[seq] = distBuf;
    }
}

void Stats::deleteOverRepSeqDist() {
    map<string, long>::iterator iter;
    for(iter = mOverRepSeq.begin(); iter!=mOverRepSeq.end(); iter++) {
        string seq = iter->first;
        delete mOverRepSeqDist[seq];
        mOverRepSeqDist[seq] = NULL;
    }
}

