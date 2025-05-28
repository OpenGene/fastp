/*
MIT License

Copyright (c) 2021 Shifu Chen <chen@haplox.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "writer.h"
#include "util.h"
#include <string.h>

Writer::Writer(Options* opt, string filename, int compression, bool isSTDOUT){
	mCompression = compression;
	mFilename = filename;
	mCompressor = NULL;
	mZipped = false;
	haveToClose = true;
	mBuffer = NULL;
	mBufDataLen = 0;
	mOptions = opt;
	mBufSize = mOptions->writerBufferSize;
	mSTDOUT = isSTDOUT;
	init();
}

Writer::~Writer(){
	flush();
	if(haveToClose) {
		close();
	}
}

void Writer::flush() {
	if(mBufDataLen > 0) {
		writeInternal(mBuffer, mBufDataLen);
		mBufDataLen = 0;
	}
}

string Writer::filename(){
	return mFilename;
}

void Writer::init(){
	mBuffer = (char*) malloc(mBufSize);
	if(mBuffer == NULL) {
		error_exit("Failed to allocate write buffer with size: " + to_string(mBufSize));
	}
	if(mSTDOUT) {
		mFP = stdout;
		return ;
	}
	if (ends_with(mFilename, ".gz")){
		mCompressor = libdeflate_alloc_compressor(mCompression);
		if(mCompressor == NULL) {
			error_exit("Failed to alloc libdeflate_alloc_compressor, please check the libdeflate library.");
		}
		mZipped = true;
		mFP = fopen(mFilename.c_str(), "wb");
		if(mFP == NULL) {
			error_exit("Failed to write: " + mFilename);
		}
	} else {
		mFP = fopen(mFilename.c_str(), "wb");
		if(mFP == NULL) {
			error_exit("Failed to write: " + mFilename);
		}
		//mOutStream = new ofstream();
		//mOutStream->open(mFilename.c_str(), ifstream::out);
	}
}

bool Writer::writeString(const string& str) {
	return write(str.data(), str.length());
}

bool Writer::writeString(string* str) {
	return write(str->data(), str->length());
}

bool Writer::write(const char* strdata, size_t size) {
	if(size + mBufDataLen > mBufSize)
		flush();
	if(size > mBufSize)
		return writeInternal(strdata, size);
	else {
		memcpy(mBuffer + mBufDataLen, strdata, size);
		mBufDataLen += size;
	}
	return true;
}

bool Writer::writeInternal(const char* strdata, size_t size) {
	size_t written;
	bool status;
	
	if(mZipped){
		size_t bound = libdeflate_gzip_compress_bound(mCompressor, size);
		void* out = malloc(bound);
		size_t outsize = libdeflate_gzip_compress(mCompressor, strdata, size, out, bound);
		if(outsize == 0)
			status = false;
		else {
			size_t ret = fwrite(out, 1, outsize, mFP );
			status = ret>0;
			//mOutStream->write((char*)out, outsize);
			//status = !mOutStream->fail();
		}
		free(out);
	}
	else{
		size_t ret = fwrite(strdata, 1, size, mFP );
		status = ret>0;
	}
	return status;
}

void Writer::close(){
	if (mZipped){
		if (mCompressor){
			libdeflate_free_compressor(mCompressor);
			mCompressor = NULL;
		}
	}
	if(mBuffer) {
		free(mBuffer);
		mBuffer = NULL;
	}
	if(mFP && !mSTDOUT) {
		fclose(mFP);
		mFP = NULL;
	}
}

bool Writer::isZipped(){
	return mZipped;
}