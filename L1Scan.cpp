#include <iostream>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <seqan/seeds.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>


using namespace std;
using namespace seqan;
typedef Seed<Simple> TSeed;
typedef SeedSet<TSeed> TSeedSet;
typedef Iterator<TSeedSet>::Type TIterator;


// Function to convert a DNA string (up to 16 bases) to an integer
unsigned int dnaStringToInt(char* dna, int k) {
    unsigned int result = 0;
    
    // Iterate over the DNA string and convert each base to its 2-bit representation
    for (size_t i = 0; i < k; ++i) {
        result <<= 2;  // Shift the result by 2 bits
        switch (dna[i]) {
            case 'A': case 'a': result |= 0b00; break;  // A = 00
            case 'C': case 'c': result |= 0b01; break;  // C = 01
            case 'G': case 'g': result |= 0b10; break;  // G = 10
            case 'T': case 't': result |= 0b11; break;  // T = 11
        }
    }

    return result;
}

typedef unsigned int kmer_t;
typedef map<kmer_t, vector<int> > kmerMap_t;

void StoreKmerToPos(string &seq, int k, kmerMap_t &map) {
  for (int i=0; i + k  < seq.size(); i++) {
    kmer_t kmer = dnaStringToInt(&seq[i], k);
    if (map.find(kmer) == map.end()) {
      map[kmer] = vector<int>();
    }
    map[kmer].push_back(i);
  }
}

// Function to fetch reference sequence near insertion point
std::string fetch_reference_sequence(faidx_t* fai, const std::string& chrom, int start, int length) {
  int len = 0;
  char* seq = fai_fetch(fai, (chrom + ":" + std::to_string(start + 1) + "-" + std::to_string(start + length)).c_str(), &len);
  std::string ref_seq(seq, len);
  free(seq);
  return ref_seq;
}

// Function to extract query sequence from BAM alignment
std::string get_query_sequence(const bam1_t* aln) {
  const uint8_t* seq = bam_get_seq(aln);
  int len = aln->core.l_qseq;
  std::string query(len, 'N');
  for (int i = 0; i < len; ++i) {
    query[i] = seq_nt16_str[bam_seqi(seq, i)];
  }
  return query;
}


// Function to perform global alignment with traceback
int globalAlignmentWithTraceback(
    const std::string &seq1,
    const std::string &seq2,
    int matchScore,
    int mismatchPenalty,
    int gapPenalty) 
{
    size_t n = seq1.size() + 1; // Rows for seq1 (length + 1 for gap)
    size_t m = seq2.size() + 1; // Columns for seq2 (length + 1 for gap)

    // DP matrix to store alignment scores
    std::vector<std::vector<int>> dp(n, std::vector<int>(m, 0));
    // Matrix to store traceback pointers (0: diagonal, 1: up, 2: left)
    std::vector<std::vector<int>> traceback(n, std::vector<int>(m, -1));

    // Initialize the first row and column for global alignment
    for (size_t i = 0; i < n; ++i) {
        dp[i][0] = i * gapPenalty; // Gap penalties for seq1
        traceback[i][0] = 1;      // Pointer: up
    }
    for (size_t j = 0; j < m; ++j) {
        dp[0][j] = j * gapPenalty; // Gap penalties for seq2
        traceback[0][j] = 2;      // Pointer: left
    }

    // Variables to track the maximum score and its position
    int maxScore = dp[0][0];
    size_t maxRow = 0, maxCol = 0;

    // Fill the DP matrix
    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 1; j < m; ++j) {
            // Match/mismatch score
            int matchMismatch = dp[i - 1][j - 1] +
                                (seq1[i - 1] == seq2[j - 1] ? matchScore : mismatchPenalty);
            // Gap penalties
            int gapSeq1 = dp[i - 1][j] + gapPenalty; // Gap in seq2
            int gapSeq2 = dp[i][j - 1] + gapPenalty; // Gap in seq1

            // Calculate the max score for this cell and set the traceback pointer
            if (matchMismatch >= gapSeq1 && matchMismatch >= gapSeq2) {
                dp[i][j] = matchMismatch;
                traceback[i][j] = 0; // Diagonal
            } else if (gapSeq1 >= gapSeq2) {
                dp[i][j] = gapSeq1;
                traceback[i][j] = 1; // Up
            } else {
                dp[i][j] = gapSeq2;
                traceback[i][j] = 2; // Left
            }

            // Update max score and position if needed
            if (dp[i][j] > maxScore) {
                maxScore = dp[i][j];
                maxRow = i;
                maxCol = j;
            }
        }
    }

    // Traceback to build the alignment strings
    std::string alignedSeq1, alignedSeq2;
    size_t i = maxRow, j = maxCol;
    while (i > 0 || j > 0) {
        if (traceback[i][j] == 0) { // Diagonal
            alignedSeq1 += seq1[i - 1];
            alignedSeq2 += seq2[j - 1];
            --i;
            --j;
        } else if (traceback[i][j] == 1) { // Up
            alignedSeq1 += seq1[i - 1];
            alignedSeq2 += '-';
            --i;
        } else { // Left
            alignedSeq1 += '-';
            alignedSeq2 += seq2[j - 1];
            --j;
        }
    }

    // Reverse the aligned sequences (traceback builds them backwards)
    std::reverse(alignedSeq1.begin(), alignedSeq1.end());
    std::reverse(alignedSeq2.begin(), alignedSeq2.end());


    // Print the alignment
    /*
    std::cout << "\nAlignment:\n";
    std::cout << alignedSeq1 << "\n";
    std::cout << alignedSeq2 << "\n";
    */
    return maxRow;
}

// Function to find TSD length by comparing reference and inserted sequence
int find_before_tsd_length(const std::string& ref_flank, const std::string& ins_seq, int max_tsd_len, string &tsd) {
  int endPolyA=ins_seq.size();
  max_tsd_len=min(max_tsd_len, endPolyA-1);
  while (endPolyA > 0 and toupper(ins_seq[endPolyA-1]) == 'A') { endPolyA--;}
  
  string insSub = ins_seq.substr(endPolyA - max_tsd_len, max_tsd_len);
  /*
  cout << "before tsd: " << endl;
  cout << ref_flank << endl;
  cout << insSub << endl;
  */
  int tsdLen=0;
  string revRef(ref_flank);
  reverse(insSub.begin(), insSub.end());
  reverse(revRef.begin(), revRef.end());
  std::transform(revRef.begin(), revRef.end(), revRef.begin(), ::toupper);  
  tsdLen = globalAlignmentWithTraceback(insSub, revRef, 2, -6, -4);
  tsd = insSub.substr(0,tsdLen);
  return tsdLen;
}

// Function to find TSD length by comparing reference and inserted sequence
int find_after_tsd_length(const std::string& ref_flank, const std::string& ins_seq, int max_tsd_len, string &tsd) {
  int tsdLen;
  string insSub = ins_seq.substr(0,max_tsd_len);
  string refSub = ref_flank;
  std::transform(refSub.begin(), refSub.end(), refSub.begin(), ::toupper);
  /*  
  cout << "After:" << endl;
  cout << insSub << endl;
  cout << refSub << endl;
  */
  tsdLen = globalAlignmentWithTraceback(insSub, refSub, 2, -6, -4);
  tsd = ins_seq.substr(0, tsdLen);
  return tsdLen;
}


// Function to read sequences from a FASTA file
std::vector<std::pair<std::string, std::string>> readFasta(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error opening FASTA file: " + filename);
    }

    std::vector<std::pair<std::string, std::string>> sequences;
    std::string line, header, seq;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!header.empty() && !seq.empty()) {
                sequences.emplace_back(header, seq);
                seq.clear();
            }
            header = line.substr(1); // Remove '>'
        } else {
            seq += line;
        }
    }

    if (!header.empty() && !seq.empty()) {
        sequences.emplace_back(header, seq);
    }

    file.close();
    return sequences;
}



typedef unsigned int kmer_t;
typedef map<kmer_t, vector<int> > kmerMap_t;


void StoreKmerMatches( kmerMap_t &queryMap, kmerMap_t &targetMap, int k, TSeedSet &seedSet) {
  for (auto it : queryMap) {
    if (targetMap.find(it.first) != targetMap.end()) {
      for (auto qp: queryMap[it.first]) {
	for (auto tp: targetMap[it.first]) {
	  TSeed seed(qp,tp,k);
	  addSeed(seedSet, seed, Single());
	}
      }
    }
  }
}

void  CalculateAlignmentIdentity(Align<DnaString> const & align, float &alnIdentity, int &matchCount) {
    // Get the number of columns (length of the alignment)
    size_t numColumns = length(row(align, 0));

    // Iterate through the alignment and compare bases at each column
    int firstMatch =-1, lastMatch=0;
    for (size_t i = 0; i < numColumns; ++i) {
        // Check if both sequences have a base at the same position (not gaps)
        if (row(align, 0)[i] != '-' && row(align, 1)[i] != '-') {
	  
            // Increment match count if both bases are the same
            if (row(align, 0)[i] == row(align, 1)[i]) {
	      if (firstMatch == -1) { firstMatch = i;}
	      lastMatch =i;
                ++matchCount;
            }
        }
    }

    // Calculate alignment identity as the fraction of matching positions
    numColumns = lastMatch - firstMatch + 1;
    
    alnIdentity = static_cast<float>(matchCount) / numColumns;
    
}

// Function to reverse complement a DNA sequence
std::string MakeReverseComplement(const std::string& seq) {
    std::string revComp = seq;
    std::reverse(revComp.begin(), revComp.end());
    for (char& base : revComp) {
        switch (base) {
            case 'A': base = 'T'; break;
            case 'T': base = 'A'; break;
            case 'C': base = 'G'; break;
            case 'G': base = 'C'; break;
        }
    }
    return revComp;
}

// Function to compute alignments and output differences
void AlignBidirectional( std::string& query, std::string& ref,
			 int k,
			 kmerMap_t &refMap,
			 int &alnDir,
			 float &alnIdentity,
			 int &nMatches,
			 vector<int> &opPos, vector<int> &opLens, vector<string> &opSeq) {
    
    // Initialize the WFA aligner
  int match = 0;
  int mismatch = 1;
  int gap_opening = 3;
  int gap_extension = 1;

  map<kmer_t, vector<int> > queryMap, revQueryMap;;
  std::string revQuery = MakeReverseComplement(query);  
  StoreKmerToPos(query, k, queryMap);
  StoreKmerToPos(revQuery, k, revQueryMap); 

  TSeedSet queryMatches, revQueryMatches;
  StoreKmerMatches(queryMap, refMap, k, queryMatches);
  StoreKmerMatches(revQueryMap, refMap, k, revQueryMatches);  

  String<Seed<Simple> > chain, revChain;
  Align<DnaString, ArrayGaps> alignment, revAlignment;
  int result=0, revResult=0;
  Score<int, Simple> scoringScheme(2, -1, -2);
  //  AlignConfig<false, false, false, false> alignConfig;
  AlignConfig<true, true, true, true> alignConfig;  
    int chainIndex =0;  

    /*
  for (i=0; i + 1 < length(chain); i++) {
    cout << beginPositionH(chain[i]) <<  "\t" << beginPositionH(chain[i+1]) - beginPositionH(chain[i]) - k <<   "\t" << beginPositionV(chain[i]) << "\t" << beginPositionV(chain[i+1]) - beginPositionV(chain[i]) - k << endl;
  }
    */
    int ph, pv;
  if (length(queryMatches) > 0) {
    chainSeedsGlobally(chain, queryMatches, SparseChaining());
    int pri =0;
    /*
    cout << "Before pruning:\t" << length(queryMatches) << "\t" << length(chain) << endl;
    ph=0; pv=0;
    for (auto c: chain) {
      cout << pri << "\t" << beginPositionH(c) << "\t" << beginPositionV(c) <<  "\t" << beginPositionH(c) - ph << "\t" << beginPositionV(c) - pv << endl;
      ph=beginPositionH(c);
      pv=beginPositionV(c);
      pri++;
    }
    */
    int lc=length(chain);
    chainIndex=length(chain);
    chainIndex=lc;
    for (; chainIndex > 1; chainIndex--) {
      if (beginPositionH(chain[chainIndex-1]) - beginPositionH(chain[chainIndex-2]) - k < 30 and
	  beginPositionV(chain[chainIndex-1]) - beginPositionV(chain[chainIndex-2]) - k < 30) {
	break;
      }
      else {
	eraseBack(chain);
	lc=length(chain);
      }
    }
    for (chainIndex=0; chainIndex + 1 < length(chain); chainIndex++) {
      if (beginPositionH(chain[chainIndex+1]) - beginPositionH(chain[chainIndex]) - k < 30 and
	  beginPositionV(chain[chainIndex+1]) - beginPositionV(chain[chainIndex]) - k < 30) {
	break;
      }
    }
    int chainLength =length(chain);
    if (chainIndex > 0) {
      erase(chain, 0, chainIndex+1);
      chainLength=length(chain);
    }
    /*
    cout << "After pruning " << endl;
    ph=0; pv=0;
    pri=0;
    for (auto c: chain) {
      cout << pri << "\t" << beginPositionH(c) << "\t" << beginPositionV(c) <<  "\t" << beginPositionH(c) - ph << "\t" << beginPositionV(c) - pv << endl;
      ph=beginPositionH(c);
      pv=beginPositionV(c);
      pri++;
    }
    */
    //    cout << "After pruning: " << "\t" << length(chain) << endl;
    /*
  for (i=0; i + 1 < length(chain); i++) {
    cout << beginPositionH(chain[i]) <<  "\t" << beginPositionH(chain[i+1]) - beginPositionH(chain[i]) - k <<   "\t" << beginPositionV(chain[i]) << "\t" << beginPositionV(chain[i+1]) - beginPositionV(chain[i]) - k << endl;
    }*/

    /*    cout << "Forward chain " << endl;
    i=0;
    for (auto c: chain) {
      cout << i << "\t" << beginPositionH(c) << "\t" << beginPositionV(c) << endl;
      i++;
    }
    */
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), query);
    assignSource(row(alignment, 1), ref);
    chainLength = length(chain);
    ph=0,pv=0;    
    if (  chainLength > 2) {
      int lastChain = chainLength - 1;
      int querySpan = beginPositionV(chain[lastChain])  - beginPositionV(chain[0]);
      int refSpan   = beginPositionH(chain[lastChain])  - beginPositionH(chain[0]);
      int span = min(querySpan, refSpan);
      int minAnchors = (int)(0.5 * ( span / k ));

      
      if (chainLength > minAnchors) {
	int lc=length(chain);
	int pci=0;
	/*
	for (auto c: chain) {
	  cout << pci << "\t" << beginPositionH(c) << "\t" << beginPositionV(c) <<  "\t" << beginPositionH(c) - ph << "\t" << beginPositionV(c) - pv << endl;
	  ph=beginPositionH(c);
	  pv=beginPositionV(c);
	  pci++;
	}
	*/
	try {
	  result = bandedChainAlignment(alignment, chain, scoringScheme, scoringScheme,  alignConfig, 10);
	  //	  cout << "Alignment " << alignment << endl;
	  //	  cout << query << endl;
	}
	catch(seqan::ClassTest::AssertionFailedException e) {
	  result=0;
	  ofstream queryOut("query_out.seq");
	  queryOut << query << endl;
	  queryOut.close();
	  ofstream refOut("ref_out.seq");
	  refOut << ref << endl;
	  refOut.close();	  
	}

      }
      else {
	//	cout << "Skipping for banded alignment " << length(chain) << "\t" << minAnchors << endl;
      }
      /*
      cout << "Forward banded chain: " << endl;
      cout << alignment << endl;
      */
    }
  }
  if (length(revQueryMatches) > 0) {
    chainSeedsGlobally(revChain, revQueryMatches, SparseChaining());
    resize(rows(revAlignment), 2);
    assignSource(row(revAlignment, 0), revQuery);
    assignSource(row(revAlignment, 1), ref);
    int i =0;
    
    int nh, ch, nv, cv, ph, pv;
    ph=0; pv=0;
    /*
    cout << "Rev chain before pruning: " << length(revChain) << endl;    
    for (auto c: revChain) {
      cout << i << "\t" << beginPositionH(c) << "\t" << beginPositionV(c) << "\t" << beginPositionH(c) - ph << "\t" << beginPositionV(c) - pv << endl;
      ph = beginPositionH(c);
      pv = beginPositionV(c);
      i++;
    }
    */
    int revChainLength = length(revChain);    
    for (i=length(revChain); i > 1; i--) {
      ph = beginPositionH(revChain[i-1]);
      ch = beginPositionH(revChain[i-2]);
      pv = beginPositionV(revChain[i-1]);
      cv = beginPositionV(revChain[i-2]);
      
      if (beginPositionH(revChain[i-1]) - beginPositionH(revChain[i-2]) - k < 30 and
	  beginPositionV(revChain[i-1]) - beginPositionV(revChain[i-2]) - k < 30) {
	break;
      }
      else {
	eraseBack(revChain);
      }
    }
    for (i=0; i + 1 < length(revChain); i++) {
      nh = beginPositionH(revChain[i+1]);
      ch = beginPositionH(revChain[i]);
      nv = beginPositionV(revChain[i+1]);
      cv = beginPositionV(revChain[i]);
      if (beginPositionH(revChain[i+1]) - beginPositionH(revChain[i]) - k < 30 and
	  beginPositionV(revChain[i+1]) - beginPositionV(revChain[i]) - k < 30) {
	break;
      }
    }
    if (i > 0) {
      erase(revChain, 0, i+1);
    }
    //    cout << "Rev chain after pruning: " << length(revChain) << endl;    
    if (length(revChain) > 2) {
      int lastRevChain = length(revChain)-1;
      int querySpan = beginPositionV(revChain[lastRevChain])  - beginPositionV(revChain[0]);
      int refSpan   = beginPositionH(revChain[lastRevChain])  - beginPositionH(revChain[0]);
      int span = min(querySpan, refSpan);
      int minAnchors = (int)(0.5 * ( span / k ));
      revChainLength = length(revChain);
      if (revChainLength > minAnchors) {
	try {
	  revResult = bandedChainAlignment(revAlignment, revChain, scoringScheme, scoringScheme,  alignConfig, 10);
	  //	  cout << "rev alignment " << revAlignment << endl;
	}
	catch(seqan::ClassTest::AssertionFailedException e) {	  
	  revResult =0;
	}
      }
      else {
	//	cout << "Skipping rev banded alignment " << length(revChain) << "\t" << minAnchors << endl;
      }
      /*
      cout << "Reverse banded chain: " << endl;
      cout << revAlignment << endl;
      */
    }      
  }


  if (result > revResult and result > 0 and length(alignment) > 0) {
    CalculateAlignmentIdentity(alignment, alnIdentity, nMatches);
  }
  else if (revResult > result and revResult > 0 and length(revAlignment) > 0) {
    CalculateAlignmentIdentity(revAlignment, alnIdentity, nMatches);
  }
}


int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <input.bam> <reference.fa> <L1HS.fasta" << std::endl;
    return 1;
  }

  // Open BAM file
  samFile* bamFile = sam_open(argv[1], "r");
  if (!bamFile) {
    std::cerr << "Error: Unable to open BAM file." << std::endl;
    return 1;
  }

  // Load BAM header
  bam_hdr_t* header = sam_hdr_read(bamFile);

  // Open reference fasta and create an index
  faidx_t* fai = fai_load(argv[2]);
  if (!fai) {
    std::cerr << "Error: Unable to load reference fasta file." << std::endl;
    sam_close(bamFile);
    return 1;
  }
  auto l1hsAll = readFasta(argv[3]);
  auto l1hs = l1hsAll[0].second;
  // Initialize alignment
  bam1_t* aln = bam_init1();
  kmerMap_t refMap;
  int k=7;
  StoreKmerToPos(l1hs, k, refMap);
  // Maximum length to search for TSDs
  const int max_tsd_len = 50;

  // Iterate through alignments
  int nReads=0;

  cout << "#chrom" << "\t" << "refPos" << "\t" << "readName" << "\t" << "strand" << "\t" << "insLen" << "\t" << "alnIdentity" << "\t" << "before_tsd_len" << "\t" << "beforeTsd" << "\t" << "after_tsd_len" << "\t" << "afterTsd" << "\t" << "regionAlnIdentity" << "\t" << "regionNumMatches" << "\t" << "inserted_seq" << endl;  
  while (sam_read1(bamFile, header, aln) >= 0) {
    ++nReads;
    if (nReads % 10000 == 0) {
      cerr << "Processed " << nReads << endl;
    }
    uint32_t* cigar = bam_get_cigar(aln);
    int32_t refPos = aln->core.pos; // 0-based alignment position
    std::string chrom = header->target_name[aln->core.tid];
    std::string query_seq = get_query_sequence(aln);
    if (query_seq.size() < 2) {
      continue;
    }
    std::string read_name(bam_get_qname(aln));
    bool isReverseStrand = aln->core.flag & BAM_FREVERSE;
    char strand = '+';
    if (isReverseStrand) { strand = '-';}
    
    int32_t queryPos = 0;            // Query position within the read
	
    for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
      uint32_t op = bam_cigar_op(cigar[i]);
      uint32_t opLen = bam_cigar_oplen(cigar[i]);
      if (op == BAM_CINS and opLen > 500 and opLen < 9000) { // Check for insertion
	// Extract inserted sequence from the query
	int query_start = queryPos;
	std::string inserted_seq = query_seq.substr(query_start, opLen);
	int alnDir;
	float alnIdentity = 0;
	vector<int> opPos, opLens;
	vector<string> opSeq;

	if (inserted_seq.size() > 500) {
	  int numL1Matches=0;
	  AlignBidirectional( inserted_seq, l1hs, k, refMap,alnDir,
			      alnIdentity, numL1Matches,
			      opPos, opLens, opSeq);

	  /*	  fastaOut << ">" << nReads<< "/" << queryPos << "/" << inserted_seq.size() << endl;
	  fastaOut << inserted_seq << endl;
	  */

	  
	  if (alnIdentity > 0.9) {
	    // Fetch reference flank
	    //	    cout << "Read " << nReads << " chrom: " << chrom << " pos: " << refPos << " Query of len " << inserted_seq.size() << " Identity " << alnIdentity << endl;	  	    
	    int flank_length = std::min(std::min(max_tsd_len, refPos), (int) inserted_seq.size());
	    std::string ref_before_flank = fetch_reference_sequence(fai, chrom, refPos - flank_length, flank_length);
	    std::string ref_after_flank = fetch_reference_sequence(fai, chrom, refPos, flank_length);


	  
	    // Compare reference flank with start of inserted sequence
	    string beforeTsd, afterTsd;
	    int before_tsd_len = find_before_tsd_length(ref_before_flank, inserted_seq, flank_length, beforeTsd);
	    int after_tsd_len = find_after_tsd_length(ref_after_flank, inserted_seq, flank_length, afterTsd);
	    int tsdLen = max(before_tsd_len, after_tsd_len);
	    /*
	    if (tsdLen > 0) {
	     
	      std::cout << "TSD detected at " << chrom << ":" << refPos << ", Insertion length: " << opLen
			<< ", TSD length: "   << tsdLen << "\t" << before_tsd_len << " " << beforeTsd << "\t" << after_tsd_len << " " << afterTsd  << std::endl;
	      std::string insertedSeqSuffix(inserted_seq, inserted_seq.size() - flank_length, inserted_seq.size());
	      std::cout << insertedSeqSuffix << std::endl;
	    }
	    */
	    int refRegionStart = max(0, (int) refPos - (int) inserted_seq.size());


	    int refChromLen = faidx_seq_len(fai, chrom.c_str()); // Replace "chr1" with the desired chromosome
	    int refRegionEnd   = min((int) refChromLen, (int) (refPos + inserted_seq.size()));
	    std::string refRegion = fetch_reference_sequence(fai, chrom, refRegionStart, refRegionEnd - refRegionStart);
	    std::transform(refRegion.begin(), refRegion.end(), refRegion.begin(), ::toupper);	    
	    
	    if (beforeTsd == "") { beforeTsd = "NA";}
	    if (afterTsd == "") { afterTsd = "NA";}
	    //	    cout << chrom << "\t" << refPos << "\t" << read_name << "\t" << strand << "\t" << inserted_seq.size() << "\t" << alnIdentity << "\t" << before_tsd_len << "\t" << beforeTsd << "\t" << after_tsd_len << "\t" << afterTsd << endl;
	    int numMatches=0;
	    kmerMap_t refRegionMap;
	    int k=7;
	    float regionAlnIdentity=0;
	    int regionNumMatches=0;
	    StoreKmerToPos(refRegion, k, refRegionMap);
	    AlignBidirectional( inserted_seq, refRegion, k, refRegionMap, alnDir, regionAlnIdentity, regionNumMatches,
				opPos, opLens, opSeq);
	    cout << chrom << "\t" << refPos << "\t" << read_name << "\t" << strand << "\t" << inserted_seq.size() << "\t" << alnIdentity << "\t" << before_tsd_len << "\t" << beforeTsd << "\t" << after_tsd_len << "\t" << afterTsd << "\t" << regionAlnIdentity << "\t" << regionNumMatches << "\t" << inserted_seq << endl;
	    
	  }
	}
      }
      switch (op) {
      case BAM_CSOFT_CLIP:
	queryPos += opLen;
	break;
      case BAM_CMATCH:
      case BAM_CEQUAL: // Match
	queryPos += opLen;
	refPos   += opLen;
	break;
      case BAM_CDIFF: // Mismatch
	queryPos += opLen; 
	refPos += opLen;
	break;
      case BAM_CINS: // Insertion
	queryPos+= opLen;
	break;
      case BAM_CDEL: // Deletion
	refPos+=opLen;
	break;
      }
    }
  }

  // Cleanup
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(bamFile);
  fai_destroy(fai);

  return 0;
}
