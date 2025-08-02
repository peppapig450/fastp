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

#include "fastqreader.h"
#include "read.h"
#include "util.h"
#include <string.h>
#include <sys/stat.h>
#include <cassert>
#include <cstddef>
#include <memory>
#include <utility>

#define FQ_BUF_SIZE (1<<23)
#define IGZIP_IN_BUF_SIZE (1<<22)
#define GZIP_HEADER_BYTES_REQ (1<<16)

namespace {

// Return file size or 0 on failure/stdin (and other special files)
inline auto stat_file_size(const std::string& filename) -> std::size_t {
	// never try to stat the std-in/out/err pseudo-files
	if (starts_with(filename, "/dev/std")) {
		return 0;
	}

	struct stat fileStatus {};
	if (stat(filename.c_str(), &fileStatus) != 0) {
		return 0;
	}

	// only report a size for regular files
	if (!S_ISREG(fileStatus.st_mode)) {
		return 0;
	}

	return static_cast<std::size_t>(fileStatus.st_size);
}
}

FastqReader::FastqReader(string filename, bool hasQuality, bool phred64){
	mFilename = filename;
	mZipped = false;
	mFile = NULL;
	mStdinMode = false;
	mFastqBuf = new char[FQ_BUF_SIZE];
	mBufDataLen = 0;
	mBufUsedLen = 0;
	mHasNoLineBreakAtEnd = false;
	mGzipInputBufferSize = IGZIP_IN_BUF_SIZE;
	mGzipInputBuffer = new unsigned char[mGzipInputBufferSize];
	mGzipOutputBufferSize = FQ_BUF_SIZE;
	mGzipOutputBuffer = (unsigned char*)mFastqBuf;
	mCounter = 0;
	mPhred64 = phred64;
	mHasQuality = hasQuality;
	mHasNoLineBreakAtEnd = false;
	mGzipInputUsedBytes = 0;
	mReadPool = NULL;
	mFileSize = stat_file_size(filename);
	init();
}

FastqReader::~FastqReader(){
	close();
	delete[] mFastqBuf;
	delete[] mGzipInputBuffer;
}

bool FastqReader::hasNoLineBreakAtEnd() {
	return mHasNoLineBreakAtEnd;
}

void FastqReader::setReadPool(ReadPool* rp) {
	mReadPool = rp;
}


bool FastqReader::bufferFinished() {
	if(mZipped) {
		return eof() && mGzipState.avail_in == 0;
	} else {
		return eof();
	}
}

void FastqReader::readToBufIgzip(){
	mBufDataLen = 0;
	while(mBufDataLen == 0) {
		if(eof() && mGzipState.avail_in==0)
			return;
		if (mGzipState.avail_in == 0) {
			mGzipState.next_in = mGzipInputBuffer;
			mGzipState.avail_in = fread(mGzipState.next_in, 1, mGzipInputBufferSize, mFile);
			mGzipInputUsedBytes += mGzipState.avail_in;
		}
		mGzipState.next_out = mGzipOutputBuffer;
		mGzipState.avail_out = mGzipOutputBufferSize;

		int ret = isal_inflate(&mGzipState);
		if (ret != ISAL_DECOMP_OK) {
			error_exit("igzip: encountered while decompressing file: " + mFilename);
		}
		mBufDataLen = mGzipState.next_out - mGzipOutputBuffer;
		if(eof() || mGzipState.avail_in>0)
			break;
	}
	// this block is finished
	if(mGzipState.block_state == ISAL_BLOCK_FINISH) {
		// a new block begins
		if(!eof() || mGzipState.avail_in > 0) {
			if (mGzipState.avail_in == 0) {
				isal_inflate_reset(&mGzipState);
				mGzipState.next_in = mGzipInputBuffer;
				mGzipState.avail_in = fread(mGzipState.next_in, 1, mGzipInputBufferSize, mFile);
				mGzipInputUsedBytes += mGzipState.avail_in;
			} else if (mGzipState.avail_in >= GZIP_HEADER_BYTES_REQ){
				unsigned char* old_next_in = mGzipState.next_in;
				size_t old_avail_in = mGzipState.avail_in;
				isal_inflate_reset(&mGzipState);
				mGzipState.avail_in = old_avail_in;
				mGzipState.next_in = old_next_in;
			} else {
				size_t old_avail_in = mGzipState.avail_in;
				memmove(mGzipInputBuffer, mGzipState.next_in, mGzipState.avail_in);
				size_t added = 0;
				if(!eof()) {
					added = fread(mGzipInputBuffer + mGzipState.avail_in, 1, mGzipInputBufferSize - mGzipState.avail_in, mFile);
					mGzipInputUsedBytes += added;
				}
				isal_inflate_reset(&mGzipState);
				mGzipState.next_in = mGzipInputBuffer;
				mGzipState.avail_in = old_avail_in + added;
			}
			int ret = isal_read_gzip_header(&mGzipState, &mGzipHeader);
			if (ret != ISAL_DECOMP_OK) {
				error_exit("igzip: invalid gzip header found");
			}
		}
	}

	if(eof() && mGzipState.avail_in == 0) {
		// all data was processed - fail if not at logical end of zip file (truncated?)
		if (mGzipState.block_state != ISAL_BLOCK_FINISH || !mGzipState.bfinal) {
			error_exit("igzip: unexpected eof");
		}
	}
}

void FastqReader::readToBuf() {
	mBufDataLen = 0;
	if(mZipped) {
		readToBufIgzip();
	} else {
		if(!eof())
			mBufDataLen = fread(mFastqBuf, 1, FQ_BUF_SIZE, mFile);
	}
	mBufUsedLen = 0;

	if(bufferFinished() && mBufDataLen>0) {
		if(mFastqBuf[mBufDataLen-1] != '\n')
			mHasNoLineBreakAtEnd = true;
	}
}

void FastqReader::init(){
	if (ends_with(mFilename, ".gz")){
		mFile = fopen(mFilename.c_str(), "rb");
		if(mFile == NULL) {
			error_exit("Failed to open file: " + mFilename);
		}
		isal_gzip_header_init(&mGzipHeader);
		isal_inflate_init(&mGzipState);
		mGzipState.crc_flag = ISAL_GZIP_NO_HDR_VER;
		mGzipState.next_in = mGzipInputBuffer;
		mGzipState.avail_in = fread(mGzipState.next_in, 1, mGzipInputBufferSize, mFile);
		mGzipInputUsedBytes += mGzipState.avail_in;
		int ret = isal_read_gzip_header(&mGzipState, &mGzipHeader);
		if (ret != ISAL_DECOMP_OK) {
			error_exit("igzip: Error invalid gzip header found: "  + mFilename);
		}
		mZipped = true;
	}
	else {
		if(mFilename == "/dev/stdin") {
			mFile = stdin;
		}
		else
			mFile = fopen(mFilename.c_str(), "rb");
		if(mFile == NULL) {
			error_exit("Failed to open file: " + mFilename);
		}
		mZipped = false;
	}
	readToBuf();
}

void FastqReader::getBytes(size_t& bytesRead, size_t& bytesTotal) {
	if(mZipped) {
		bytesRead = mGzipInputUsedBytes - mGzipState.avail_in;
	} else {
		bytesRead = ftell(mFile);//mFile.tellg();
	}

	bytesTotal = mFileSize;
}

void FastqReader::clearLineBreaks(char* line) {

	// trim \n, \r or \r\n in the tail
	int readed = strlen(line);
	if(readed >=2 ){
		if(line[readed-1] == '\n' || line[readed-1] == '\r'){
			line[readed-1] = '\0';
			if(line[readed-2] == '\r')
				line[readed-2] = '\0';
		}
	}
}

bool FastqReader::eof() {
	return feof(mFile);//mFile.eof();
}

void FastqReader::getLine(string* line){
	int copied = 0;

	int start = mBufUsedLen;
	int end = start;

	while(end < mBufDataLen) {
		if(mFastqBuf[end] != '\r' && mFastqBuf[end] != '\n')
			end++;
		else
			break;
	}

	// this line well contained in this buf, or this is the last buf
	if(end < mBufDataLen || bufferFinished()) {
		int len = end - start;
		line->assign(mFastqBuf+start, len);

		// skip \n or \r
		end++;
		// handle \r\n
		if(end < mBufDataLen-1 && mFastqBuf[end-1]=='\r' && mFastqBuf[end] == '\n')
			end++;

		mBufUsedLen = end;

		return ;
	}

	// this line is not contained in this buf, we need to read new buf
	line->assign(mFastqBuf+start, mBufDataLen - start);

	while(true) {
		readToBuf();
		start = 0;
		end = 0;
		// handle the case that \r or \n in the start of buf
		if(line->empty()) {
			while(start < mBufDataLen && (mFastqBuf[start] == '\r' || mFastqBuf[start] == '\n'))
				start++;
			end = start;
		}
		while(end < mBufDataLen) {
			if(mFastqBuf[end] != '\r' && mFastqBuf[end] != '\n')
				end++;
			else
				break;
		}
		// this line well contained in this buf
		if(end < mBufDataLen || bufferFinished()) {
			int len = end - start;
			line->append(mFastqBuf+start, len);

			// skip \n or \r
			end++;
			// handle \r\n
			if(end < mBufDataLen-1 && mFastqBuf[end] == '\n')
				end++;

			mBufUsedLen = end;
			return;
		}
		// even this new buf is not enough, although impossible
		line->append(mFastqBuf+start, mBufDataLen);
	}

	return;
}

Read* FastqReader::read(){
	if(mBufUsedLen >= mBufDataLen && bufferFinished()) {
		return nullptr;
	}

	std::string name;
	std::string sequence;
	std::string strand;
	std::string quality;

	getLine(&name);
	// name should start with @
	while((name.empty() && (mBufUsedLen < mBufDataLen || !bufferFinished())) || (!name.empty() && name[0]!='@')){
		getLine(&name);
	}

	if (name.empty()) {
		return nullptr;
	}

	getLine(&sequence);
	getLine(&strand);
	getLine(&quality);

	if (strand.empty() || strand[0]!='+') {
		cerr << name << endl;
		cerr << "Expected '+', got " << strand << endl;
		cerr << "Your FASTQ may be invalid, please check the tail of your FASTQ file" << endl;

		return nullptr;
	}

	if(quality.length() != sequence.length()) {
		cerr << "ERROR: sequence and quality have different length:" << endl;
		cerr << name << endl;
		cerr << sequence << endl;
		cerr << strand << endl;
		cerr << quality << endl;
		cerr << "Your FASTQ may be invalid, please check the tail of your FASTQ file" << endl;

		return nullptr;
	}

	Read* readInPool = nullptr;
	if (mReadPool != nullptr) {
		readInPool = mReadPool->getOne();
	}

	if (readInPool != nullptr) {
		// Update the pooled read's contents
		*readInPool = Read(std::move(name),
						   std::move(sequence),
						   std::move(strand),
						   std::move(quality),
						   mPhred64);
		return readInPool;
	}

	return new Read(std::move(name),
					std::move(sequence),
					std::move(strand),
					std::move(quality),
					mPhred64);
}

void FastqReader::close(){
	if (mFile){
		fclose(mFile);
		mFile = NULL;
	}
}

bool FastqReader::isZipFastq(string filename) {
	if (ends_with(filename, ".fastq.gz"))
		return true;
	else if (ends_with(filename, ".fq.gz"))
		return true;
	else if (ends_with(filename, ".fasta.gz"))
		return true;
	else if (ends_with(filename, ".fa.gz"))
		return true;
	else
		return false;
}

bool FastqReader::isFastq(string filename) {
	if (ends_with(filename, ".fastq"))
		return true;
	else if (ends_with(filename, ".fq"))
		return true;
	else if (ends_with(filename, ".fasta"))
		return true;
	else if (ends_with(filename, ".fa"))
		return true;
	else
		return false;
}

bool FastqReader::isZipped(){
	return mZipped;
}

bool FastqReader::test(){
	FastqReader reader1("testdata/R1.fq");
	FastqReader reader2("testdata/R1.fq");
	int i=0;
	while(true){
		i++;
		std::unique_ptr<Read> r1(reader1.read());
		std::unique_ptr<Read> r2(reader2.read());
		if(r1 == nullptr || r2 == nullptr) {
			break;
		}
		r1->print();
		r2->print();
	}
	return true;
}

FastqReaderPair::FastqReaderPair(FastqReader* left, FastqReader* right){
	mLeft = left;
	mRight = right;
}

FastqReaderPair::FastqReaderPair(string leftName, string rightName, bool hasQuality, bool phred64, bool interleaved){
	mInterleaved = interleaved;
	mLeft = new FastqReader(leftName, hasQuality, phred64);
	if(mInterleaved)
		mRight = NULL;
	else
		mRight = new FastqReader(rightName, hasQuality, phred64);
}

FastqReaderPair::~FastqReaderPair(){
	if(mLeft){
		delete mLeft;
		mLeft = NULL;
	}
	if(mRight){
		delete mRight;
		mRight = NULL;
	}
}

ReadPair* FastqReaderPair::read(){
	std::unique_ptr<Read> left(mLeft->read());
	std::unique_ptr<Read> right;
	if (mInterleaved) {
		right.reset(mLeft->read());
	} else {
		right.reset(mRight->read());
	}

	if (left == nullptr || right == nullptr) {
		return nullptr;
	}

	return new ReadPair(std::move(left), std::move(right));
}
