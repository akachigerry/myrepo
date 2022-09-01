#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cstddef>
#include <ctime>
#include <map>
#include <unordered_map>
#include <sstream>		//allows to use istringstream
#include <memory>		//required for the smart pointers
#include <future>
#include <bitset>
#include <functional>
#include <thread>


struct SourceAttrib
{
	int depth = 0;
	std::unordered_map<std::vector<bool>, int> src_of_src; //map of 1st base of source of source kmer and no of times they overlap. 
	std::map<std::string, int> read_id;
};

typedef std::unordered_map<std::vector<bool>, std::unordered_map<std::vector<bool>, SourceAttrib>> hashMap;
typedef std::unordered_map<char, std::pair<bool, bool>> CharBoolMap;
typedef std::unordered_map<std::string, char> BoolCharMap;

std::vector<std::string> split1(char const *s, char delimiter)
{
	std::vector<std::string> tokens;
	std::string token;
	std::istringstream tokenStream(s);
	while (std::getline(tokenStream, token, delimiter))
	{
		tokens.push_back(token);
	}
	return tokens;
}


std::string uniq_str_gen(std::string &str)
{
	//sort string
	std::string str1 = str;
	std::sort(str1.begin(), str1.end());
	//then use std::unique to generate iterator of set of unique elements. 
	//Iterator points to position beyond last item and doesn't resize string afer duplicates removed.
	std::string::iterator str_iter;
	str_iter = std::unique(str1.begin(), str1.end());
	str1.resize(std::distance(str1.begin(), str_iter)); //std::distance determines number of items within range specified by both iterator arguments
	return str1;
}

std::string cannonseq(std::string &read)
{
	std::string rev_read;
	std::string cannonical_read;
	std::vector<std::string> readVec;
	readVec.push_back(read);
	for (std::string::reverse_iterator rit = read.rbegin(); rit != read.rend(); ++rit)
	{
		char base = *rit;
		if (base == 'A')
		{
			base = 'T';
			rev_read.push_back(base);
			continue;
		}
		if (base == 'C')
		{
			base = 'G';
			rev_read.push_back(base);
			continue;
		}
		if (base == 'G')
		{
			base = 'C';
			rev_read.push_back(base);
			continue;
		}
		if (base == 'T')
		{
			base = 'A';
			rev_read.push_back(base);
			continue;
		}		
	}
	readVec.push_back(rev_read);
	std::sort(readVec.begin(), readVec.end());
	cannonical_read = readVec.front();
	return cannonical_read;
}

int fileLineCounter(std::string fastqname)
{
	std::ifstream fastqfile(fastqname);
	std::string fastqline;
	int lineCount = 0;
	while (std::getline(fastqfile, fastqline))
	{
		++lineCount;
	}
	return lineCount;
}

void decompress(BoolCharMap &bool2char, std::vector<bool> &boolvec, std::string &str)
{
	int boolvecsz = boolvec.size();
	for (int i = 0; i < boolvecsz; i += 2)
	{
		std::string boolstr;
		if (boolvec[i] == true)
		{
			boolstr = "true";
		}
		else
		{
			boolstr = "false";
		}
		if (boolvec[i+1] == true)
		{
			boolstr = boolstr + "_" + "true";
		}
		else
		{
			boolstr = boolstr + "_" + "false";
		}
		str.push_back(bool2char[boolstr]);
	}
}



void fastqfilebreak(int num_of_procs, int fileLineCnt, std::string filename)
{
	std::cout << "Breaking " << filename << " into " << num_of_procs << " chunks\n";
	int num_of_reads;
	num_of_reads = (int) (fileLineCnt / 4);
	int read_chunkSize = (int) (num_of_reads / num_of_procs);
	int chunkSize_tracker = 0;

	std::ifstream fastqfile(filename);

	for (int i = 1; i <= num_of_procs; ++i)
	{
		std::string chunk_num = std::to_string(i);
		std::string chunkfileName = filename + "_" + chunk_num;
		std::ofstream chunkfileOut(chunkfileName);

		int linecount = 0;
		int seqlineFound = 0;
		std::string fastqline;
		while (std::getline(fastqfile, fastqline)) //this file reading loop is repeated 'num_of_procs' times. But each time,
			//the file reading continues from the line where it broke off.
		{			
			if (seqlineFound == 0 && fastqline[0] != '@')
			{
				linecount += 1;
				continue;
			}

			if (seqlineFound == 0 && fastqline[0] == '@' && std::count(fastqline.begin(), fastqline.end(), '@') == 1)
			{
				seqlineFound += 1;
			}			
			linecount += 1;
			if (linecount <= read_chunkSize)
			{
				chunkfileOut << fastqline << "\n";
			}
			else
			{
				break;
			}
		}
	}
	std::cout << filename << " broken into " << num_of_procs << " chunks\n";
}


hashMap graphBuild(std::string &&filename, int &&k, CharBoolMap &&Char2bool, BoolCharMap &&bool2char)
{
	int loopcnt = 0, loopcnt1 = 1;
	std::string readId;
	std::string fastqline;
	hashMap dbgraph;
	std::ifstream fastqfile1(filename);
	
	while (std::getline(fastqfile1, fastqline))
	{
		if (fastqfile1.fail() && !fastqfile1.eof())
		{
			std::cout << "First Fastq file could not be read " << "\n";
			EXIT_FAILURE;
		}
		//if only end of file is reached, and there's no file read failure,
		//simply break out from while loop.
		if (!fastqfile1.fail() && fastqfile1.eof())
		{
			break;
		}
		std::string uniq_fastqline;
		
		if (loopcnt == 0 || (loopcnt % 4 == 0 && fastqline.front() == '@' &&\
			std::count(fastqline.begin(), fastqline.end(), '@') == 1)) //for read id line
			//which is 4th 0-base index line from begining
		{											
			std::cout << "seq title line: " << fastqline << "\n";
			std::vector<std::string> fastqlineVec = split1(fastqline.c_str(), ':');
			if (fastqlineVec.size() == 3) //means reads from ion torrent
			{
				readId = fastqlineVec[1] + fastqlineVec[2];
			}
			else if (fastqlineVec.size() == 5 && fastqlineVec.back().find("#") != std::string::npos) //means reads from old illumina
			{
				std::string ycoord = fastqlineVec.back();
				ycoord = split1(ycoord.c_str(), '#')[0];
				readId = fastqlineVec[2] + fastqlineVec[3] + ycoord;
			}
			else if (fastqlineVec.size() == 10 && fastqlineVec[6].find(" ") != std::string::npos)
			{
				std::string ycoord = fastqlineVec[6];
				ycoord = split1(ycoord.c_str(), ' ')[0];
				readId = fastqlineVec[4] + fastqlineVec[5] + ycoord;
			}
			else  //every other read id condition
			{
				std::cout << "Fastq reads are not of the right format.\n";
				EXIT_FAILURE;
			}
			loopcnt += 1;
			continue;
		}

		if (fastqline.size() < k + 2)
		{
			loopcnt += 1;
			continue;
		}	
		uniq_fastqline = uniq_str_gen(fastqline);		
		if (uniq_fastqline != "ACGT" && uniq_fastqline != "ACGNT")
		{				
			loopcnt += 1;
			continue;
		}

		loopcnt += 1;
		std::cout << "unique read string: " << uniq_fastqline << "\n";
		fastqline = cannonseq(fastqline); //better to cannonise entire read seq than each kmer from it.
										  //Reasons to be explained later.
		int fastqlineSize = fastqline.size();	
		std::cout << "seq length: " << fastqlineSize << "\n";
		std::vector<bool> source1;
		std::vector<bool> sink;
	
		for (int j = 0; j < fastqlineSize - (k + 1); ++j)
		{
			char base, sink_base;
			bool firstbit, secondbit;

			if (j == 0)
			{
				//loop to simultaneosly take ksize substrings of fastqline as source and string and also
				//compress them into vector of bools	
				for (int idx = 0; idx < k; ++idx)
				{
					base = fastqline[idx];
					firstbit = Char2bool[base].first;
					secondbit = Char2bool[base].second;
					if (idx == k - 1) //last loop cycle ensures sink is last base of kmer and not found in source
					{
						sink.push_back(firstbit);
						sink.push_back(secondbit);
						break;
					}
					source1.push_back(firstbit);
					source1.push_back(secondbit);					
				}
			}
						
			int sourceSz = k - 1;
			int sinkSz = sourceSz;
			SourceAttrib sourceData;
			std::vector<bool> src_of_source;
			std::string src_of_source_1stBase;
			std::vector<bool> src_lastbase;
			std::vector<bool> src_of_src_1stbase_bool;						
			if (j > 0)
			{
				src_lastbase = sink; //prev. source is current 'src_of_source', therefore
				//prev. sink would last base of current source or sink of src_of_source'
				sink_base = fastqline[(j + k) - 1];
				firstbit = Char2bool[sink_base].first;
				secondbit = Char2bool[sink_base].second;
				//clear contents of prev 'sink' vector before inserting new 'sink'
				sink.clear();
				sink.push_back(firstbit);
				sink.push_back(secondbit);

				src_of_source = source1; //source1 is the previous one
				//get 1st base of src_of_source
				firstbit = src_of_source[0];
				secondbit = src_of_source[1];
				src_of_src_1stbase_bool.emplace_back(firstbit);
				src_of_src_1stbase_bool.emplace_back(secondbit);
				src_of_source_1stBase = fastqline[j - 1];			
				//Now to get current source
				source1.erase(source1.begin(), source1.begin() + 2); //..remove 1st 2 bools (1st base)
				//..then add last 2 bools (last base) of previous sink.
				source1.insert(source1.end(), src_lastbase.begin(), src_lastbase.end());
			}

			if (dbgraph.count(source1) > 0)
			{
				if (dbgraph[source1].count(sink) > 0)
				{
					sourceData = dbgraph[source1][sink];
					sourceData.depth += 1;

					if (sourceData.read_id.count(readId) > 0) //if read is already in sourceData
					{
						sourceData.read_id[readId] += 1;
					}
					else
					{
						sourceData.read_id[readId] = 1;
					}

					//steps below to take care of connection btw current source and its own source
					if (dbgraph.count(src_of_source) > 0) //src_of_source is expected to already be in graph b4 this stage
					{
						//last base of current source was sink in preceding source ('src_of_source')
						if (dbgraph[src_of_source].count(src_lastbase) > 0)
						{
							int srcOfsrc2srcOvelapCnt = dbgraph[src_of_source][src_lastbase].depth;

							//adds srcOfsrc to src overlap cnt to src_of_src data if non-existent, 
							//updates it to new cnt if already there. It is really an update cos 'srcOfsrc2srcOvelapCnt'
							//is derived depth of srcOfsrc to src connecxn which is incremented by 1 
							//each time a kmer of both are found
							sourceData.src_of_src[src_of_src_1stbase_bool] = srcOfsrc2srcOvelapCnt;

							std::cout << "Src of src 1st char: " << src_of_source_1stBase << " overlap cnt with src: " << srcOfsrc2srcOvelapCnt << "\n";
						}
					}
					dbgraph[source1][sink] = sourceData; //replace old 'sourceData' in dbgraph[source][sink] with updated one
					continue;
				}
				else //if sink is new and extra connecxn from source
				{
					//build new sourceAttrib data ('sourceData')
					sourceData.depth += 1; //Nb. increment is ok here cos 'depth' member already
					//initialised to 0 in structure.
					sourceData.read_id[readId] = 1;

					//steps below to take care of connection btw current source and its own source
					if (dbgraph.count(src_of_source) > 0) //src_of_source is expected to already be in graph b4 this stage
					{
						if (dbgraph[src_of_source].count(source1) > 0)
						{
							int srcOfsrc2srcOvelapCnt = dbgraph[src_of_source][src_lastbase].depth;
							sourceData.src_of_src[src_of_src_1stbase_bool] = srcOfsrc2srcOvelapCnt;
							std::cout << "Src of src 1st char: " << src_of_source_1stBase << " overlap cnt with src: " << srcOfsrc2srcOvelapCnt << "\n";
						}
					}
					dbgraph[source1][sink] = sourceData;
					continue;
				}
			}
			else //new source encountered
			{
				//first generate full kmer bool from source and sink bools, reverse complement
				//and check for availability of source of reversed kmer in 'dbgraph'. If
				//available, update depth of link to sink if already exists or create new sink 
				//if not
				std::vector<bool> kmer = source1; 
				kmer.insert(kmer.end(), sink.begin(), sink.end());
				int kmerSz = kmer.size();
				std::string kmerStr;
				decompress(bool2char, kmer, kmerStr);
				std::string revkmerStr = cannonseq(kmerStr);
				std::vector<bool> revkmer;
				for (auto ch : revkmerStr)
				{
					std::pair<bool, bool> bitpair = Char2bool[ch];
					revkmer.push_back(bitpair.first);
					revkmer.push_back(bit2pair.second);
				}
				
				std::vector<bool> revkmer1 = revkmer;
				//removes the last two bools representing last element of revkmer.
				//this in order to get the source of revkmer.
				revkmer1.erase(revkmer1.begin() + (k - 2), revkmer1.begin() + (k - 1)); 

				std::vector<bool> revsource1 = revkmer1;
				//revsink is last two bools of revkmer
				std::vector<bool> revsink(revkmer.begin() + (k - 2), revkmer.begin() + (k - 1));
				if (dbgraph.count(revsource1) > 0)
				{
					SourceAttrib sourceData;
					if (dbgraph[revsource1].count(revsink) > 0)
					{					
						dbgraph[revsource1][revsink].depth += 1;

						if (dbgraph[revsource1][revsink].read_id.count(readId) > 0) //if read is already in sourceData
						{
							dbgraph[revsource1][revsink].read_id[readId] += 1;
						}
						else
						{
							dbgraph[revsource1][revsink].read_id[readId] = 1;
						}						
						continue;
					}
					else //if sink is new and extra connecxn from source
					{
						sourceData.depth = 1;
						sourceData.read_id[readId] = 1;
						dbgraph[revsource1][revsink];
						dbgraph[revsource1][revsink] = sourceData;
						continue;
					}
				}
				else
				{
					sourceData.depth += 1;
					sourceData.read_id[readId] = 1;
					if (dbgraph.count(src_of_source) > 0) //though current source is not yet added, its own source would have been added 
														  //(if this is not 1st pass thru loop) with curent source added as its sink
					{
						if (dbgraph[src_of_source].count(src_lastbase) > 0)
						{
							int srcOfsrc2srcOvelapCnt = dbgraph[src_of_source][src_lastbase].depth;
							sourceData.src_of_src[src_of_src_1stbase_bool] = srcOfsrc2srcOvelapCnt;
							std::cout << "Src of src 1st char: " << src_of_source_1stBase << " overlap cnt with src: " << srcOfsrc2srcOvelapCnt << "\n";
						}
					}
					std::unordered_map<std::vector<bool>, SourceAttrib> innerMap;
					innerMap[sink] = sourceData;
					dbgraph[source1] = innerMap;
				}				
			}

		}
	}
	std::cout << "part graph size: " << dbgraph.size() << "\n";
	return dbgraph;
}


hashMap graphMerger(std::vector<std::future<hashMap>> &dbgraphVec)
{
	hashMap mergedgraph;
	int cnt = 0;
	std::cout << "waiting to build 1st graph...\n";
	auto start = std::chrono::steady_clock::now(); //current time.
	while (dbgraphVec.size() < 1) //waits until at least one process is done building graph.
								//Nb that this function is run immediately after graphBuild func
							 //w/o waiting for the latter to finish.
	{
		continue;
	}
	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	
	std::cout << "First graph built...\n";
	
	for (auto &graph : dbgraphVec)
	{		
		if (cnt == 0)
		{					
			mergedgraph = graph.get(); //calling process (main) blocks here until
			//content of future object 'graph' is obtained			
			std::cout << "Graph building completed for fastq chunk " << cnt + 1 << " ..now merging" << "\n";
			cnt += 1;
			continue;
		}
		std::cout << "Awaiting graph building processes\n";
		
		hashMap graph1 = graph.get();		
		cnt += 1;
		std::cout << "Graph building completed for fastq chunk " << cnt+1 << "\n";

		for (hashMap::iterator itr = graph1.begin(); itr != graph1.end(); ++itr)
		{
			std::vector<bool> sorce = itr->first;
			std::unordered_map<std::vector<bool>, SourceAttrib> sorceData;
			sorceData = itr->second;
			
			if (mergedgraph.count(sorce) > 0)
			{
				for (std::unordered_map<std::vector<bool>, SourceAttrib>::iterator itr1 = sorceData.begin(); itr1 != sorceData.end(); ++itr1)
				{
					std::vector<bool> sinke = itr1->first;
					SourceAttrib src_attrib = itr1->second;
					if (mergedgraph[sorce].count(sinke) > 0)
					{
						mergedgraph[sorce][sinke].depth += src_attrib.depth;
					
						std::map<std::string, int> readIdMap = src_attrib.read_id;
						for (std::map<std::string, int>::iterator iter1 = readIdMap.begin(); iter1 != readIdMap.end(); ++iter1)
						{
							std::string readid = iter1->first;
							int readcnt = iter1->second;
							if (mergedgraph[sorce][sinke].read_id.count(readid) > 0) //don't expect this to ever be true
							{
								mergedgraph[sorce][sinke].read_id[readid] += 1;
							}
							else
							{
								mergedgraph[sorce][sinke].read_id[readid] = 1;
							}
						}
						
						//accessing data for source's link to its source
						std::unordered_map<std::vector<bool>, int> srcOfsrc_data = src_attrib.src_of_src;						
						if (srcOfsrc_data.size() > 0)
						{
							for (std::unordered_map<std::vector<bool>, int>::iterator it2 = srcOfsrc_data.begin(); it2 != srcOfsrc_data.end(); ++it2)
							{
								std::vector<bool> srcofsrc = it2->first;
								int src2srcOfsrcCov = it2->second;
								if ((mergedgraph[sorce][sinke].src_of_src).count(srcofsrc) > 0)
								{
									mergedgraph[sorce][sinke].src_of_src[srcofsrc] += src2srcOfsrcCov;
								}
								else //if current graph section contains source of source not found in existing merged graph
								{
									mergedgraph[sorce][sinke].src_of_src[srcofsrc] = src2srcOfsrcCov;
								}
							}
						}
					}
					else
					{
						mergedgraph[sorce][sinke] = src_attrib;
					}
				}
			}
			else 
			{
				mergedgraph[sorce] = std::move(sorceData);
				continue;
			}
		}
		std::cout << "Graph segment: " << cnt << " integrated into merged graph" << "\n";
	}
	return mergedgraph;
}


int modalCovCalc(hashMap &dbgraph1)
{
	std::map<int, int> covfreq;
	for (hashMap::iterator it = dbgraph1.begin(); it != dbgraph1.end(); ++it)
	{
		std::vector<bool> source = it->first;
		std::unordered_map<std::vector<bool>, SourceAttrib> sourceData = it->second;

		if (sourceData.size() == 1)
		{
			std::unordered_map<std::vector<bool>, SourceAttrib>::iterator it1 = sourceData.begin();
			int cov = (it1->second).depth;
			std::cout << "coverage: " << cov << "\n";
			//coverage approximated to nearest 5x interval value
			int remainder = cov % 5;
			if (remainder < 3)
			{				
				if (cov < 3)
				{
					cov = 1;
				}
				else
				{
					cov -= remainder;
				}
			}
			else
			{
				cov += (5 - remainder);
			}
			if (covfreq.count(cov) > 1)
			{
				covfreq[cov] += 1;
			}
			else
			{
				covfreq[cov] = 1;
			}
		}
		if (sourceData.size() > 1)
		{
			int locusCov = 0;
			for (std::unordered_map<std::vector<bool>, SourceAttrib>::iterator it1 = sourceData.begin(); it1 != sourceData.end(); ++it1)
			{
				int cov = (it1->second).depth;
				locusCov += cov;
			}
			int remainder = locusCov % 5;
			if (remainder % 5 < 3)
			{
				if (locusCov < 3)
				{
					locusCov = 1;
				}
				else
				{
					locusCov -= remainder;
				}
			}
			else
			{
				locusCov += remainder;
			}
			if (covfreq.count(locusCov) > 1)
			{
				covfreq[locusCov] += 1;
			}
			else
			{
				covfreq[locusCov] = 1;
			}
		}
	}
	int modalcov = 0;
	int maxfreq = 0;
	for (std::map<int, int>::iterator it2 = covfreq.begin(); it2 != covfreq.end(); ++it2)
	{
		int cov = it2->first;
		if (cov < 2)
		{
			continue;
		}		
		int freq = it2->second;
		if (freq > maxfreq)
		{
			maxfreq = freq;
			modalcov = cov;
		}
	}
	std::cout << "modal coverage: " << modalcov << "\n";
	return modalcov;
}


//below function implements procedure which finds the longest tandem repeat substring in input string
std::string tandemRepeatFinder(std::string &str)
{
	std::map <std::string, int> repeatMap;
	int strSz = str.size();
	for (int i = 0; i < strSz; ++i)
	{
		int minrepeatSz = 2; //initialised with 2 cos i dont want repeat unit size < 2
		int maxrepeatSz = (int)strSz / 2;
		for (int s = minrepeatSz; s < maxrepeatSz; ++s)
		{
			int repeatSz = s;
			if (i + repeatSz >= strSz) //if below substr step is allowed under this condition segfault will result
			{
				break;
			}
			std::string repeat = str.substr(i, repeatSz);
			if (repeatMap.count(repeat) > 0)
			{
				repeatMap[repeat] += 1;
			}
			else
			{
				repeatMap[repeat] = 1;
			}
		}
	}
	std::string longestRepeat("");
	std::vector<std::string> tandemRepeats;
	for (std::map <std::string, int>::iterator it = repeatMap.begin(); it != repeatMap.end(); ++it)
	{
		std::string repeat = it->first;
		int repeatSz = repeat.size();
		int repeatCnt = it->second;
		if (repeatCnt < 2) //selected repeat must be found at least twice within input string
		{
			continue;
		}
		else
		{
			//collate all tandemly repeated substrings (found at least twice)
			std::size_t matchPos = str.find(repeat); //matchpos points to pos of 1st char of the 1st 'repeat' in 'str'
			std::size_t prematchPos = matchPos;
			int tandemRepeatFound = 0;
			while (strSz - matchPos > repeatSz) //continue search if unsearched section of 'str' is at least the size of 'repeat'
			{
				matchPos = str.find(repeat, matchPos+repeatSz); //search 4 another match of 'repeat' from pos followng pos of last char of prev match
				if (matchPos - prematchPos == repeatSz) //must be true for tandemness to be confirmed
				{
					tandemRepeatFound += 1;
					prematchPos = matchPos;
				}
				else
				{
					tandemRepeatFound = 0;
					break;
				}
			}
			if (tandemRepeatFound > 1)
			{
				if (repeat.size() > longestRepeat.size())
				{
					longestRepeat = repeat;
				}
			}
			
		}
	}
	return longestRepeat;
}


std::vector<std::pair<std::string, std::string> > contigGen(hashMap &mergedgraph, int modalcov, BoolCharMap &bool2Char)
{
	std::vector<std::pair<std::string, std::string> > contigs;
	std::string hap1_contigs;
	std::string hap2_contigs;	
	hashMap::iterator src_it = mergedgraph.begin();
	
	std::unordered_map<std::vector<bool>, SourceAttrib> sourcedata = src_it->second;
	std::unordered_map<std::vector<bool>, SourceAttrib>::iterator sink_it = sourcedata.begin();
	std::unordered_map<std::vector<bool>, int> sources_of_source = sink_it->second.src_of_src;

	int repeatFound = 0;
	//ensure there's just a single path from first source and that it's not within repeat
	while (sourcedata.size() > 1 || sources_of_source.size() > 1)
	{
		++src_it;
		sourcedata = std::move(src_it->second);
		sink_it = sourcedata.begin();
		sources_of_source = sink_it->second.src_of_src;
	}
	std::vector<bool> source_bool = std::move(src_it->first);
	int srcBoolSz = source_bool.size();

	
	//activebubble, tuggles btw 0 and 1 to indicate if if there's active bubble
	int newcontig = 0, bubblecnt = 0, bubblelen = 0;
	std::map<std::string, int> hap1_readsData;
	std::map<std::string, int> hap2_readsData;
	std::map<std::string, int> hap1BubbleReads;
	std::map<std::string, int> hap2BubbleReads;
	//contig generator loop. Doesn't end until all paths in graph traversed.
	int breakit = 0;
	int potMiddleRepeat = 0;
	while (1)
	{
		std::vector<bool> hap1Src_bool;
		std::vector<bool> hap2Src_bool;

		//get string format of source
		std::string source_str("");
		for (int y = 0; y<srcBoolSz; y += 2)
		{
			bool firstBit = source_bool[y];
			bool secondBit = source_bool[y + 1];
			std::string firstBit_str;
			std::string secondBit_str;
			std::string tru("true");
			std::string fals("false");
			if (firstBit == true)
			{
				firstBit_str = tru;
			}
			else
			{
				firstBit_str = fals;
			}
			if (secondBit == true)
			{
				secondBit_str = tru;
			}
			else
			{
				secondBit_str = fals;
			}
			std::string Key = firstBit_str + "_" + secondBit_str;
			char base = bool2Char[Key];
			std::string base1;
			base1.push_back(base);
			source_str += base1;
		}

		int terminateContig = 0;
		if (mergedgraph.count(source_bool) > 0)
		{
			int souceSz = source_str.size();
			std::unordered_map<std::vector<bool>, SourceAttrib> sourcedata = mergedgraph[source_bool];
			std::unordered_map<std::vector<bool>, SourceAttrib>::iterator itr1 = sourcedata.begin();
			int locuscov = 0;
			if (newcontig == 0) //new contig is 1st primed with entire source seq. last base of sink seq are added subsequently
			{
				SourceAttrib Src_attrib = itr1->second;
				int cov = Src_attrib.depth;
				double locuscov2modeRatio = (double)cov / modalcov;
				if (locuscov2modeRatio >= 2.0)
				{
					continue; //can't start contig within potential repeat.
				}
				hap1_contigs += source_str;
				hap2_contigs += source_str;
				newcontig += 1;
			}			
			std::vector<bool> sink = itr1->first;
			SourceAttrib Src_attrib = itr1->second;
			
			int cov = Src_attrib.depth;
			int sourceDataSz = sourcedata.size();
			std::unordered_map<std::vector<bool>, int> sources_of_src;
			if (sourceDataSz > 1) //where there's a bubble, cov will take into acct all sinks
			{
				sources_of_src = Src_attrib.src_of_src; //NB. 'sources_of_src' has been taken from 1 sink though there's more
				//than a single sink from current source. This is ok cos all the sinks have same 'sources_of_src' data
				for (std::unordered_map<std::vector<bool>, SourceAttrib>::iterator it2 = sourcedata.begin(); it2 != sourcedata.end(); ++it2)
				{
					std::vector<bool> sink1 = it2->first;
					int hapcov = (it2->second).depth;
					locuscov += hapcov;
				}
			}
			else
			{
				locuscov = cov;
			}
			double locuscov2modeRatio = (double)locuscov / modalcov;

			int falseRepeat = 0;
			if (locuscov2modeRatio > 2.0) //potential repeat
			{
				std::unordered_map<std::vector<bool>, SourceAttrib> srcData = mergedgraph[source_bool];
				std::vector<bool> firstSource;
				int srcDataSz = srcData.size();
				int srcNotFound = 0;
				int srcOfsrcSz = sources_of_src.size();
				std::cout << "Locus cov > 2x modal cov. potential repeat\n";

				//if >2x modal cov is not accompanied by two inward connections (edges) to source,
				//and 2 outward connexns, then its not considered true repeat.
				if (potMiddleRepeat == 0 && srcDataSz < 2 && srcOfsrcSz <= 1)
				{
					std::cout << "false repeat found\n";
					falseRepeat += 1;
				}

				//any potential middle of tandem repeat found is 1st marked by increment of 
				//potMiddleRepeat to 1, then graph traversal is continued w/o contig extention.
				//However, potential end of the repeat picked from the middle is always
				//checked and when found (on basis of > 2x modal cov and 2 sinks), the code
				//enters the repeat resolution section below
				int maxcov = 0;
				int mincov = 0;
				if (potMiddleRepeat > 0 && srcDataSz >= 2)
				{
					std::cout << "potential repeat likely found middle of repeat\n";
					std::vector<bool> loopsink;
					potMiddleRepeat = 0;
					int i = 0;
					for (std::unordered_map<std::vector<bool>, SourceAttrib>::iterator sit = srcData.begin(); sit != srcData.end();++sit)
					{
						int cov = (sit->second).depth;
						if (i == 0)
						{
							maxcov = cov;
							mincov = cov;
							loopsink = sit->first;
							i += 1;
							continue;
						}
						if (cov > maxcov) //sink with maxcov at end of tadem repeat unit
							//theoretically loops back to connect to 1st source of repeat unit
						{
							maxcov = cov;
							loopsink = sit->first;
						}
						else
						{
							mincov = cov;
						}
					}
					//maxcov represent cov of sink potentially looped back to start of tandem repeat
					//mincov represent cov of sink exiting repeat end.
					if ((double)maxcov / modalcov >= 2 && ((double)mincov / modalcov > 0.5 && (double)mincov / modalcov < 2))
					{
						srcOfsrcSz += 1; //cos we assume a real repeat is found and a real repeat 
						//has > 1 inward edges to repeat start. the increment also signals to
						//section of code below dealing with repeats
						firstSource = source_bool;
						firstSource.erase(firstSource.begin(), firstSource.begin() + 2);
						//potential 1st source of potential repeat obtained
						firstSource.insert(firstSource.end(), loopsink.begin(), loopsink.end());
					}
					else //sinks of fairly similar cov could indicate end of non-tandem 
						//duplicated segments. Otherwise if mincov is too small, might be error
					{
						if ((double)mincov / modalcov <= 0.5) //if mincov is too small, it could 
							//be error. In expected repeat end is false, assume no repeat
						{
							std::cout << "repeat end is judged false. Likely the result of error\n";
							falseRepeat += 1;
						}
						else //likely non-tandem duplicate but has to be proven below
						{
							srcOfsrcSz += 1;
						}
					}
				}
				

				if (srcOfsrcSz > 1) //if repeat is suspected (tandem repeats and duplcate segments)
									//Nb. for duplicated segments sources_of_src might be > 2 while it must be 2 for tandem repeats
				{
					std::cout << "potential repeat found from repeat start\n";
					std::cout << "modal cov: " << modalcov << "; Repeat cov: " << locuscov << "\n";
					int modeRepeatCov = 0;
					std::map<int, int> repeatCovDist;
					int repeatcov_av = 0;
					std::string repeatUnit;

					std::string tandemRepeatUnit("");
					int tandemRepeatUnitSz = 0;
					std::vector<bool> repeat_1stbase;

					if (firstSource.size() == 0)
					{
						firstSource = source_bool; //1st source of repeat
					}
					decompress(bool2Char, firstSource, repeatUnit);

					repeat_1stbase.push_back(source_bool[0]); //1st bool(bit) from source_bool used to initialise repeat_1stbase..
					repeat_1stbase.push_back(source_bool[1]);//..2nd bool  completes the bool representatn of 1st base 

					int tandemRepeat = 0, prematureRepeatTerminatn = 0, longTandemUnit = 0;
					int cnt1 = 0, repeatUnitComplete = 0;
					int kmerPerRead = 0, repeatcovFail = 0;
					for (std::map<std::string, int>::iterator readIter = Src_attrib.read_id.begin(); readIter != Src_attrib.read_id.end(); ++readIter)
					{
						if (readIter->second > 1) //kmer found more than once in same read or fragment
						{
							kmerPerRead = readIter->second;
							break;
						}
					}

					if (kmerPerRead > 1) //more than a single src-sink unit (kmer) in a read/fragment suggests tandem repeat
					{
						std::cout << "Repeat is probably short tandem\n";
						tandemRepeat += 1;
						//determine if src-sink (kmer unit) size is > repeat unit and find repeat unit if it is
						tandemRepeatUnit = tandemRepeatFinder(source_str);
						tandemRepeatUnitSz = tandemRepeatUnit.size();
						if (tandemRepeatUnitSz > 0)
						{
							std::cout << "Repeat is definitely short tandem\n";
						}
					}
					std::vector<bool> sinke;
					std::map<std::string, int> read_Ids;
					std::vector<bool> lastsink;
					std::vector<bool> lastsource;
					std::vector <std::string> duplicSegReads;
					int tandemRepeatEnd = 0, duplicateRepeatEnd = 0;
					int duplicateSegment = 0, bubblelength = 0, bubblecnt1 = 0;
					for (int i = 0; i < 10000; ++i) //its assumed no repeat/paralog/duplicated segment unit could be > 10kb
					{
						if (mergedgraph.count(source_bool) > 0)
						{
							srcData = mergedgraph[source_bool];
							srcDataSz = srcData.size();
							std::unordered_map<std::vector<bool>, SourceAttrib>::iterator iter = srcData.begin();

							SourceAttrib src_attrib = iter->second;

							int src2sink_cov = src_attrib.depth;
							int locuscov = 0;
							std::vector<std::vector<bool> > sinkPaths; //collates the diff sinks leading from source 

							if (srcDataSz > 1) //If not error, suggests SNV bubble within repeat or end of repeat.   
											   //Snv bubble here could only suggest duplicated segments /paralogs/gene families
											   //if its not a repeat end indicator
							{
								bubblecnt1 += 1;
								bubblelength = 0; //bubblelength>0 indicate bubble due to dupl segm
												  //hence is reinitialised for every new variant in readiness for 
												  //finding a duplicate segment
								int count = 1;
								int nonRepeatSink = 0;

								//obtain total cov for all variants at locus
								for (std::unordered_map<std::vector<bool>, SourceAttrib>::iterator sit = srcData.begin(); sit != srcData.end();++sit)
								{
									locuscov += (sit->second).depth;
								}
								if (cnt1 == 0) //cnt1=0 for only for 1st repeat nuceotide 
								{
									repeatcov_av = locuscov; //repeatcov_av is moving cov ave over repeat unit
									if (tandemRepeat == 0 && tandemRepeatUnitSz == 0) //if repeat is not short tandem, then it could be
																					  //duplicate segment or long tandem unit. Impt to determine this with 1st repeat base as traversed sources
																					  //won't be deleted for duplicate segments. This would allow traversal from alternate paths 
																					  //(often nearby seqs of trans-located duplicates) to the duplicate segment.
									{
										duplicateSegment += 1;
									}
									cnt1 += 1;
								}
								else
								{
									repeatcov_av = (int)(repeatcov_av + locuscov) / 2;
								}
								cnt1 += 1;
								int truehap = 0;
								for (std::unordered_map<std::vector<bool>, SourceAttrib>::iterator sitr = srcData.begin(); sitr != srcData.end();++sitr)
								{
									std::vector<bool> sink1 = sitr->first;

									int hapcov = (sitr->second).depth;
									std::vector<bool> fullsink;
									fullsink.push_back(source_bool[2]); //we need all bases of source except 1st one to form
																		//fullsink. hence initiate fullsink with 3rd bool w starts bit code for 2nd base
									fullsink.insert(fullsink.end(), source_bool.begin() + 3, source_bool.end()); //insert rest of source bools
									fullsink.insert(fullsink.end(), sink1.begin(), sink1.end()); //finally add bools of last sink base	
									if (fullsink == firstSource)//end of tandem repeat unit confirmed by sink 
																//equating 1st repeat source
									{
										tandemRepeatEnd += 1;
										longTandemUnit += 1;
										duplicateSegment = 0; //if repeat is tandem, then its not duplicateSegment.
									}

									if ((double)(hapcov / locuscov) < 1 / (2 * srcDataSz))//hapcov to locuscov ratio must be at least 
																						  //1/2 of 1/srcDataSz (if srcDataSz=2, equiv to 25% of tot cov (locuscov) across duplicate loci). In this case
																						  //it's not so the hap allele classed as error and removed.
									{
										if (fullsink != firstSource) //sink can be proven error and then deleted
																	 //only if it doesn't indicate end of tandem repeat 
										{
											mergedgraph[source_bool].erase(sink1);
										}
										else //update 'read_Ids' with read data of sinks of current bubble if 1st bubble
										{
											if (bubblecnt1 == 1 && read_Ids.size() == 0)
											{
												std::map<std::string, int> currentReads = mergedgraph[source_bool][sink1].read_id;
												read_Ids.insert(currentReads.begin(), currentReads.end());
												sinke = sink1;
											}
										}
									}
									else
									{
										sinkPaths.push_back(sink1); //sink1 collated in list of paths 4rm currnt source cos it
																	//passes above test
									}

									//at least 1 of the sinks connectng from source with cov closer to modalcov (representative of 
									//non-repeat coverages) is indicative of end of duplicate segm if tandemRepeat == 0
									if (duplicateSegment > 0 && mergedgraph[source_bool].size()>1) //duplicateSegment > 0 means short tandemRepeat=0
									{
										if ((double)modalcov / hapcov < 1.7 && (double)hapcov / modalcov < 1.7)
										{
											truehap += 1;
										}
										if (truehap > 1)
										{
											duplicateRepeatEnd += 1;
										}
									}
								}
								if (mergedgraph[source_bool].size() == 1) //if error trimming of bubble leaves only 1 remaining sink,
																		  //no chance of duplicate segment ending here. Hence its reverted to 0
																		  //Nb tandemRepeatEnd and longTandemUnit would remain 0 so no equating to 0.
								{
									duplicateRepeatEnd = 0;
									bubblecnt = 0; //since there's no longer a bubble
								}
								
								if (mergedgraph[source_bool].size() > 2 && duplicateSegment > 0) //source connexn to > 2 sinks only
																								 //possible for duplicate segments. updating traversed bubble length only relevant for 
																								 //duplicate segments
								{
									bubblelength += 1;
								}
							}
							else //if there no bubble at current pos of repeat unit
							{
								std::unordered_map<std::vector<bool>, SourceAttrib>::iterator srcData_it = srcData.begin();
								sinke = srcData_it->first;
								locuscov = (srcData_it->second).depth;
								if (cnt1 == 0) //cnt1=0 for only for 1st repeat nuceotide 
								{
									repeatcov_av = locuscov;
									cnt1 += 1;
								}
								else
								{
									repeatcov_av = (int)(repeatcov_av + locuscov) / 2;
								}
								cnt1 += 1;
								if (tandemRepeat == 0 && tandemRepeatUnitSz == 0)
								{
									duplicateSegment += 1;
								}
							}


							//for the very 1st bubble under repeat conditions, whether repeat is
							//duplicate segment or not, collate reads of 
							//remember, bubblelength > 0 is indicative of duplicate seg bubble						
							if (bubblelength > 0)//while variant-induced divergence 
												 //of source-sink persists for entire (kmer-1 or source) length, collate read ids for path of bubble  
												 //followed by current segment. Will be useful in determing traversal path should another bubble be  
												 //found within duplicated segment
							{
								//bubblecnt1is indicative of any bubble
								if (bubblecnt1 == 1 && sinkPaths.size() > 0) //the very 1st variant in duplicate segment. Sink path may be randomly selected
																			 //only for this variant.
								{
									sinke = sinkPaths[0];
									SourceAttrib srcAttrib_data = mergedgraph[source_bool][sinke];
									read_Ids = srcAttrib_data.read_id;
									//now erase the sink path taken here so as not to take it again
									mergedgraph[source_bool].erase(sinke);
								}
								if (bubblecnt1 > 1) //From 2nd varaint pos, sink path selection must consider read ids of preceeding variants
								{
									if (srcDataSz > 1)
									{
										int loopcnt = 0;
										int sinkPathsSz = sinkPaths.size();
										for (auto sink : sinkPaths) //for every variant at current pos within duplicate segment
										{
											SourceAttrib src_attribData = mergedgraph[source_bool][sink];
											std::map<std::string, int> sinkReads = src_attribData.read_id;
											int matchfound = 0;
											for (std::map<std::string, int>::iterator it = sinkReads.begin(); it != sinkReads.end(); ++it)
											{
												std::string read = it->first;
												if (read_Ids.count(read) > 0) //if id of any sinkPath read matches with read ids of
																			  //previous segment variant loci/locus
												{
													sinke = sink;
													//now erase the sink path taken here so as not to take it again. Impt for
													//traversing putative duplicated segment sections of graph
													mergedgraph[source_bool].erase(sinke);
													matchfound += 1;
													break;
												}
											}
											if (matchfound > 0)
											{
												read_Ids.insert(sinkReads.begin(), sinkReads.end());
											}
											if (matchfound == 0 && loopcnt == sinkPathsSz - 1) //not finding linkage with prev variant, exit repeat loop and 
																							   //terminate contig extention
											{
												terminateContig += 1;
												//as the current source is an unresolved multi-path node of mergedgraph, delete it
												mergedgraph.erase(source_bool);
											}
											loopcnt += 1;
										}
									}
									else //for non-variant pos within variant-generated bubble, just update 'read_Ids'
										 //with their read ids
									{
										SourceAttrib src_Attrib = (srcData.begin())->second;
										read_Ids.insert((src_Attrib.read_id).begin(), (src_Attrib.read_id).end());
									}
								}
							}
							//update repeat unit with above sink ('sinke') (chosen from 1 of multiple 
							//graph paths if there's multiple sinks (bubble) from current source
							std::string Sink_str("");
							decompress(bool2Char, sinke, Sink_str);
							repeatUnit += Sink_str;
							source_bool.erase(source_bool.begin(), source_bool.begin() + 1);
							source_bool.push_back(sinke[0]);
							source_bool.push_back(sinke[1]);

							//next 3 steps to approximate locuscov to nearest multiple of five
							int rem = locuscov % 5;
							if (rem > 2)
							{
								locuscov = locuscov + (5 - rem);
							}
							else
							{
								locuscov = locuscov - rem;
							}

							if (repeatCovDist.count(locuscov) > 0)
							{
								repeatCovDist[locuscov] += 1;
							}
							else
							{
								repeatCovDist[locuscov] = 1;
							}

							//number of src-sink connexns (depth) at each repeat pos must be no more be within 2x > or < than 
							//cov ave up to the pos
							if ((double)locuscov / repeatcov_av > 2 || (double)repeatcov_av / locuscov > 2)
							{
								repeatcovFail += 1;
								int repeatUnitSz = repeatUnit.size();
								if (repeatUnitSz > 0 && (double)repeatcovFail > (0.3*repeatUnitSz))
								{
									//if the 30% or more of loci are rather closer 'modalcov' (representative non-repeat cov)
									//then we might be dealing with a non-repeat. In this case, abort repeat unit and continue
									//normal contig extension
									if (((double)locuscov / modalcov > 1.0 && (double)locuscov / modalcov <= 1.5) || \
										(double)modalcov / locuscov > 1.0 && (double)modalcov / locuscov <= 1.5)
									{
										prematureRepeatTerminatn += 1;
										lastsink = sinke;
										break;
									}
									else
									{
										terminateContig += 1;
										lastsink = sinke;
										break; //get out of repeat loop.
									}
								}
							}

							//for confirmed tandem repeats, delete mergedgraph[source] as no other path is expected to lead to it
							//(i.e. it is not expected to be duplicated in other location in genome) 
							//then exit loop once num of loops = repeat sz is completed
							if (tandemRepeatUnitSz > 0 && i == tandemRepeatUnitSz - 1)
							{
								std::cout << "short tandem repeat unit sequence obtained\n";
								repeatUnit = tandemRepeatUnit;
								mergedgraph.erase(source_bool);
							}
							if (duplicateSegment == 0 && i >= 500) //i don't expect > 500bp sz for tandem repeat unit
							{
								terminateContig += 1; //repeat was not resolved cos last repeat base not found
								lastsink = sinke;
								std::cout << "Repeat not resolved as last repeat base not found\n";
								break;
							}
							if (tandemRepeatEnd > 0 || longTandemUnit > 0)
							{
								std::cout << "Reached end of long tandem repeat unit\n";
								lastsink = sinke;
								break;
							}
							//only duplicate segments could have more than 2 sink connections from a source. The diff sinks 
							//represent the diff variants of the segments at current pos
							if (mergedgraph[source_bool].size() >= 2 && duplicateRepeatEnd > 0)
							{
								std::cout << "Reached end of duplicate repeat unit\n";
								lastsink = sinke;
								break;
							}
							if (duplicateSegment > 0 && duplicateRepeatEnd == 0 && i >= 9900) //if after 29kbp, and duplicate segm/ paralog is 
																							   //hasn't ended, stop repeat loop
							{
								terminateContig += 1; //duplicate repeat building inconclusive
								lastsink = sinke;
								std::cout << "Repeat not resolved as duplicate segment expected but not found\n";
								break;
							}
						}
						else //if source is not found in graph increment 'terminatContig'. This 
							 //would wrap up the contig and transfer control to outermost loop.
							 //Also, a 'srcNotFound' signal is incremented. This prevents attempt to
							 //update sink into next source ensuring the non-exitent source checked
							 //at the outermost loop, from where control is handed to section dealing 
							 //with cases of non-existent sources.
						{
							terminateContig += 1;
							srcNotFound += 1;
							lastsink = sinke;
						}
					} //repeat loop ends here
					if (prematureRepeatTerminatn > 0)
					{
						prematureRepeatTerminatn = 0;
						hap1_contigs += repeatUnit;
						hap2_contigs += repeatUnit;
						source_bool.erase(source_bool.begin(), source_bool.begin() + 2);
						source_bool.push_back(lastsink[0]);
						source_bool.push_back(lastsink[1]);
						continue;
					}

					if (terminateContig > 0)
					{
						terminateContig = 0;
						hap1_contigs += repeatUnit;
						hap2_contigs += repeatUnit;
						std::pair<std::string, std::string> hap1hap2;
						hap1hap2.first = hap1_contigs;
						hap1hap2.second = hap2_contigs;
						contigs.push_back(hap1hap2); //extending contig terminated and contig collated
						hap1_contigs = std::string();
						hap2_contigs = std::string();
						if (srcNotFound > 0)
						{
							srcNotFound = 0;
							continue;
						}
						source_bool.erase(source_bool.begin(), source_bool.begin() + 2);
						source_bool.push_back(lastsink[0]);
						source_bool.push_back(lastsink[1]);
						continue;
					}

					//find modal repeat cov
					int maxfreq = 0;
					int modalRepeatCov = 0;
					for (std::map<int, int>::iterator itr = repeatCovDist.begin(); itr != repeatCovDist.end(); ++itr)
					{
						int freq = itr->second;
						int cov = itr->first;
						if (freq > maxfreq && cov != 0)
						{
							maxfreq = freq;
							modalRepeatCov = cov;
						}
					}

					if (tandemRepeat > 0 || longTandemUnit > 0 || tandemRepeatEnd > 0)
					{
						std::string tandemRepeatSeg("");
						int numOfRepeats = (int)modalRepeatCov / modalcov;
						if (tandemRepeatUnit.size() > 0)
						{
							std::string shorttandemRepeatSeg("");
							if (numOfRepeats > 1)
							{
								for (int i = 0; i < numOfRepeats; ++i)
								{
									shorttandemRepeatSeg += tandemRepeatUnit;
								}
								hap1_contigs += shorttandemRepeatSeg;
								hap2_contigs += shorttandemRepeatSeg;

								source_bool.erase(source_bool.begin(), source_bool.begin() + 2);
								source_bool.insert(source_bool.end(), lastsink.begin(), lastsink.end());
								continue;
							}
							else if (numOfRepeats == 1) //approx to numOfRepeats=2
							{
								numOfRepeats = 2;
								for (int i = 0; i < numOfRepeats; ++i)
								{
									shorttandemRepeatSeg += tandemRepeatUnit;
								}
								hap1_contigs += shorttandemRepeatSeg;
								hap2_contigs += shorttandemRepeatSeg;

								source_bool.erase(source_bool.begin(), source_bool.begin() + 2);
								source_bool.insert(source_bool.end(), lastsink.begin(), lastsink.end());
								continue;
							}
							else //not a real repeat, extend contigs with repeat units and
								 //continue normal grpah traversal.
							{
								hap1_contigs += tandemRepeatUnit;
								hap2_contigs += tandemRepeatUnit;
								source_bool.erase(source_bool.begin(), source_bool.begin() + 2);
								source_bool.insert(source_bool.end(), lastsink.begin(), lastsink.end());
								continue;
							}
						}
						if (tandemRepeatEnd > 0 || longTandemUnit > 0)
						{
							std::string tandemRepeatSeg("");
							if (numOfRepeats > 1)
							{
								for (int i = 0; i < numOfRepeats; ++i)
								{
									tandemRepeatSeg += repeatUnit;
								}
								hap1_contigs += tandemRepeatSeg;
								hap2_contigs += tandemRepeatSeg;
								std::vector<bool> fullsink;
								fullsink.push_back(source_bool[2]);
								fullsink.insert(fullsink.end(), source_bool.begin() + 3, source_bool.end()); //insert rest of source bools
								fullsink.insert(fullsink.end(), lastsink.begin(), lastsink.end());
								source_bool = fullsink;
							}
							else if (numOfRepeats == 1) //approx to numOfRepeats=2
							{
								numOfRepeats = 2;
								for (int i = 0; i < numOfRepeats; ++i)
								{
									tandemRepeatSeg += repeatUnit;
								}
								hap1_contigs += tandemRepeatSeg;
								hap2_contigs += tandemRepeatSeg;

								source_bool.erase(source_bool.begin(), source_bool.begin() + 2);
								source_bool.insert(source_bool.end(), lastsink.begin(), lastsink.end());
								continue;
							}
							else //not a real repeat, extend contigs with repeat units and
								 //continue normal grpah traversal.
							{
								hap1_contigs += repeatUnit;
								hap2_contigs += repeatUnit;
								source_bool.erase(source_bool.begin(), source_bool.begin() + 2);
								source_bool.insert(source_bool.end(), lastsink.begin(), lastsink.end());
								continue;
							}
						}
					}
					if (duplicateRepeatEnd > 0)
					{
						hap1_contigs += repeatUnit;
						hap2_contigs += repeatUnit;
						std::vector<bool> fullsink;
						fullsink.push_back(source_bool[2]);
						fullsink.insert(fullsink.end(), source_bool.begin() + 3, source_bool.end()); //insert rest of source bools
						fullsink.insert(fullsink.end(), lastsink.begin(), lastsink.end());
						source_bool = fullsink;
					}
				}
			}


			if ((sources_of_src.size() < 2 && locuscov2modeRatio <= 2.0) || falseRepeat > 0) //if there's no suspected repeat
			{	
				if (falseRepeat > 0)
				{
					std::cout << "..reverting to non-repeat assembly\n";
				}
				falseRepeat = 0;
				//check if the sink leading from source is start of bubble (SNP or error) prior to haplotisation 				
				if (sourceDataSz > 1) //if there was a SNP or error
				{
					std::vector<bool> Sink1;
					std::vector<bool> Sink2;
					if (sourceDataSz > 2) //if there's > 2 branches in bubble, select only the 2 with largest coverage
					{
						std::map<int, std::vector<bool>> covs;
						for (std::unordered_map<std::vector<bool>, SourceAttrib>::iterator iter = sourcedata.begin(); iter!= sourcedata.end(); ++iter)
						{
							int cov = iter->second.depth;
							covs[cov] = iter->first;
						}
						std::map<int, std::vector<bool>>::reverse_iterator riT = covs.rbegin(); //the first item of reverse map has the 
																 //largest coverage since st::map is sorted in ascending order
						Sink1 = riT->second;
						std::cout << "largest allele cov at snp locus: " << riT->first << "\n";
						++riT; //moves pointer to index with second-largest coverage
						Sink2 = riT->second;
						std::cout << "2nd largest allele cov at snp locus: " << riT->first << "\n";
					}
					
					if (sourceDataSz == 2)
					{
						std::unordered_map<std::vector<bool>, SourceAttrib>::iterator iter = sourcedata.begin();
						Sink1 = iter->first;
						++iter;
						Sink2 = iter->first;
					}

					//check likelihood the sink is an error based on its coverage relative to total cov
					//sink wi cov < 25% of locus depth is error, tru snp picked 4rm other sink w cov > 25%
					int hap1cov = sourcedata[Sink1].depth;
					int hap2cov = sourcedata[Sink2].depth;
					locuscov = hap1cov + hap2cov;

					double hap1covRatio = (double)hap1cov / locuscov;
					double hap2covRatio = (double)hap2cov / locuscov;					

					int hap1_read1 = 0, hap1_read2 = 0;
					int hap2_read1 = 0, hap2_read2 = 0;
					if (bubblecnt == 0) //for the first ever bubble found while extending a contig
					{
						hap1BubbleReads = sourcedata[Sink1].read_id;
						hap2BubbleReads = sourcedata[Sink2].read_id;
						hap1_read1 += 1; //cos reads of sink1 partition into hap1
						hap2_read2 += 1; //cos reads of sink2 partition into hap2						
					}
					else //if there's been previous bubbles in this contig, partion alleles into
						//existing haplots based on matching reads
					{
						std::map<std::string, int> reads1 = sourcedata[Sink1].read_id;
						std::map<std::string, int> reads2 = sourcedata[Sink2].read_id;
						for (std::map<std::string, int>::iterator Itr = reads1.begin(); Itr != reads1.end(); ++Itr)
						{
							std::string Read = Itr->first;
							if (hap1BubbleReads.count(Read) > 0 && hap2BubbleReads.count(Read) == 0) //if sink1 reads overlap with hap1
							{
								std::cout << "Read 1 partions to hap1\n";
								hap1_read1 += 1;
								break;
							}
							if (hap2BubbleReads.count(Read) > 0 && hap1BubbleReads.count(Read) == 0)
							{
								std::cout << "Read 1 partions to hap2\n";
								hap2_read1 += 1;
								break;
							}
						}
						for (std::map<std::string, int>::iterator itr = reads2.begin(); itr != reads2.end(); ++itr)
						{
							std::string Read = itr->first;
							if (hap1BubbleReads.count(Read) > 0 && hap2BubbleReads.count(Read) == 0)
							{
								std::cout << "Read 2 partions to hap1\n";
								hap1_read2 += 1;
								break;
							}
							if (hap2BubbleReads.count(Read) > 0 && hap1BubbleReads.count(Read) == 0)
							{
								std::cout << "Read 2 partions to hap2\n";
								hap2_read2 += 1;
								break;
							}
						}
						//Now check if the bubble found here is the result of a true snp or error.
						//Use proportion of coverage of each 'allele'. Even if this bubble is the
						//result of error, still partition its reads into 'hap1BubbleReads' and 
						//'hap2BubbleReads' as you would for a real snp. This is okay cos the 
						//errorneous 'allele' is still found in only a fraction of the reads 
						//of the locus (albeit a very small fraction) as in a real snp and the 
						//total coverage at the locus is within modal coverage range (meaning
						//bubble does not arise from duplicate segments) as it is a snp.
						//if Sink1 is judged error use only Sink2, hence no variant here
						if (hap1covRatio < 0.25)
						{
							std::string Sinke("");
							std::string boolstr;
							if (Sink2[0] == true)
							{
								boolstr = "true";
							}
							else
							{
								boolstr = "false";
							}
							if (Sink2[1] == true)
							{
								boolstr = boolstr + "_" + "true";
							}
							else
							{
								boolstr = boolstr + "_" + "false";
							}
							Sinke.push_back(bool2Char[boolstr]);
							hap1_contigs += Sinke;
							std::cout << "hap1 contig: " << hap1_contigs << "\n";
							hap2_contigs += Sinke;
							std::cout << "hap2 contig: " << hap2_contigs << "\n";
							mergedgraph.erase(source_bool);

							std::vector<bool>::iterator src_it_beg = source_bool.begin();
							source_bool.erase(src_it_beg, src_it_beg + 2);
							source_bool.insert(source_bool.end(), Sink2.begin(), Sink2.end());
							hap1Src_bool = source_bool;
							hap2Src_bool = source_bool;
							continue;
						}
						if (hap2covRatio < 0.25)
						{
							std::string Sinke("");
							std::string boolstr;
							if (Sink1[0] == true)
							{
								boolstr = "true";
							}
							else
							{
								boolstr = "false";
							}
							if (Sink1[1] == true)
							{
								boolstr = boolstr + "_" + "true";
							}
							else
							{
								boolstr = boolstr + "_" + "false";
							}
							Sinke.push_back(bool2Char[boolstr]);
						
							hap1_contigs += Sinke;
							std::cout << "hap1 contig: " << hap1_contigs << "\n";
							hap2_contigs += Sinke;
							std::cout << "hap2 contig: " << hap2_contigs << "\n";
							mergedgraph.erase(source_bool);

							std::vector<bool>::iterator src_it_beg = source_bool.begin();
							source_bool.erase(src_it_beg, src_it_beg + 2);
							source_bool.insert(source_bool.end(), Sink1.begin(), Sink1.end());
							hap1Src_bool = source_bool;
							hap2Src_bool = source_bool;
							continue;
						}
						//partion reads at 'error' locus
						if (hap1covRatio < 0.25 || hap2covRatio < 0.25)
						{
							if (hap1_read1 > 0 || hap2_read2 > 0)
							{
								for (std::map<std::string, int>::iterator iter = reads1.begin(); iter != reads1.end(); ++iter)
								{
									std::string readid = iter->first;
									if (hap1BubbleReads.count(readid) > 0)
									{
										hap1BubbleReads[readid] += 1;
									}
									else
									{
										hap1BubbleReads[readid] = iter->second;
									}
								}
								for (std::map<std::string, int>::iterator iter = reads2.begin(); iter != reads2.end(); ++iter)
								{
									std::string readid = iter->first;
									if (hap2BubbleReads.count(readid) > 0)
									{
										hap2BubbleReads[readid] += 1;
									}
									else
									{
										hap2BubbleReads[readid] = iter->second;
									}
								}
							}
							if ((hap1_read2 > 0 || hap2_read1 > 0) && (hap1_read1 == 0 || hap2_read2 == 0))
							{
								for (std::map<std::string, int>::iterator iter = reads2.begin(); iter != reads2.end(); ++iter)
								{
									std::string readid = iter->first;
									if (hap1BubbleReads.count(readid) > 0)
									{
										hap1BubbleReads[readid] += 1;
									}
									else
									{
										hap1BubbleReads[readid] = iter->second;
									}
								}
								for (std::map<std::string, int>::iterator iter = reads1.begin(); iter != reads1.end(); ++iter)
								{
									std::string readid = iter->first;
									if (hap2BubbleReads.count(readid) > 0)
									{
										hap2BubbleReads[readid] += 1;
									}
									else
									{
										hap2BubbleReads[readid] = iter->second;
									}
								}
							}			
						}

						else //if snp is not error
						{
							hap1BubbleReads.clear(); //empty 'hap1BubbleReads' and repalce with matching current read data 
							hap2BubbleReads.clear();
							if (hap1_read1 > 0 || hap2_read2 > 0)
							{
								hap1BubbleReads = reads1;
								hap2BubbleReads = reads2;
							}
							if (hap1_read2 > 0 || hap2_read1 > 0)
							{
								hap1BubbleReads = reads2;
								hap2BubbleReads = reads1;
							}														
						}						
					}
										
					//loop kmer-1 (source size) times and the variant moves up (5' to 3') the entire length of 'sources' same as
					//length of bubble					
					std::string hap1_source;
					std::string hap2_source;

					int newSnpInBubble = 0;
					std::string Snk1("");
					std::string Snk2("");
					int maxbubbleSz = 2 * souceSz; //would allow max 2 SNPs per snp bubble
					int loopcnt = 0;
					std::vector<bool> source_bool1 = source_bool;
					std::vector<bool> source_bool2 = source_bool;
					std::vector<bool>::iterator src_it_beg = source_bool1.begin();
					source_bool1.erase(src_it_beg, src_it_beg + 2);
					source_bool1.insert(source_bool1.end(), Sink1.begin(), Sink1.end());
					hap1Src_bool = source_bool1;

					src_it_beg = source_bool2.begin();
					source_bool2.erase(src_it_beg, src_it_beg + 2);
					source_bool2.insert(source_bool2.end(), Sink2.begin(), Sink2.end());
					hap2Src_bool = source_bool2;

					std::string boolstr1;
					std::string boolstr2;
					std::string tru("true");
					std::string fals("false");
					if (Sink1[0] == true)
					{
						boolstr1 = "true";
					}
					else
					{
						boolstr1 = "false";
					}
					if (Sink1[1] == true)
					{
						boolstr1 = boolstr1 + "_" + "true";
					}
					else
					{
						boolstr1 = boolstr1 + "_" + "false";
					}

					if (Sink2[0] == true)
					{
						boolstr2 = tru;
					}
					else
					{
						boolstr2 = fals;
					}
					if (Sink2[1] == true)
					{
						boolstr2 = boolstr2 + "_" + "true";
					}
					else
					{
						boolstr2 = boolstr2 + "_" + "false";
					}
					Snk1.push_back(bool2Char[boolstr1]);
					Snk2.push_back(bool2Char[boolstr2]);
					SourceAttrib src_attrib_1 = mergedgraph[source_bool][Sink1];
					SourceAttrib src_attrib_2 = mergedgraph[source_bool][Sink2];
					std::map<std::string, int> read_data1 = src_attrib_1.read_id;
					std::map<std::string, int> read_data2 = src_attrib_2.read_id;
					if (hap1_read1 > 0 || hap2_read2 > 0)
					{						
						hap1_contigs += Snk1;
						hap2_contigs += Snk2;							
						std::cout << "hap1 contig: " << hap1_contigs << "\n";
						std::cout << "hap2 contig: " << hap2_contigs << "\n";
						hap1BubbleReads.insert(read_data1.begin(), read_data1.end());
						hap2BubbleReads.insert(read_data2.begin(), read_data2.end());
						mergedgraph.erase(source_bool);
						source_bool = hap1Src_bool;
					}
					if (hap1_read2 > 0 || hap2_read1 > 0)
					{
						hap1_contigs += Snk2;
						hap2_contigs += Snk1;
						hap1BubbleReads.insert(read_data2.begin(), read_data2.end());
						hap2BubbleReads.insert(read_data1.begin(), read_data1.end());
						mergedgraph.erase(source_bool);
						source_bool = hap1Src_bool;
					}
					//bubble loop. Nb. alleles for the variant pos (start of bubble) have already
					//been partioned into existing haplots above. Overlaps btw current allele 
					//reads and haplot reads were used for this. This loop will just extend each
					//hapsource-hapsink pair (the 1st pair already extended above) 
					//simultaneously until (hopefully) kmers for the haplots become 
					//same (converge) again at end of bubble.
					while (hap1Src_bool != hap2Src_bool)//for std snp, this is expected to be false 
						//after number of loops equal to size of source
					{
						if (loopcnt < maxbubbleSz)
						{
							loopcnt += 1;
						}
						else
						{
							break;
						}
						std::unordered_map<std::vector<bool>, SourceAttrib> hap1srcData;
						std::unordered_map<std::vector<bool>, SourceAttrib> hap2srcData;
						if (mergedgraph.count(hap1Src_bool) > 0)
						{
							hap1srcData = mergedgraph[hap1Src_bool];
						}
						else //if hap1 source not found, equate to hap2 source (if found)
							//and exit bubble loop. Nb. snp would be been added at this stage
							//and is not removed cos it met set criteria. In exiting the bubble
							//loop we do not assume the snp isn't real. Rather, the absense of
							//hap1 source here is interpreted as lack coverage 
							//over 3' neigbourhood of the snp. Same interpretation is hap2
							//source is missing.
						{
							if (mergedgraph.count(hap2Src_bool) > 0)
							{
								hap1Src_bool = hap2Src_bool;
								source_bool = hap1Src_bool;
								breakit += 1;
								break;
							}
							else
							{
								//equate source_bool to non-existent hap1Src_bool, increment
								//'breakit' and exit bubble loop. A -ve 'breakit' takes us to
								//the outermost loop where presence of hap1Src_bool is
								//tested again and when not found, contigGen is prepared for
								//termination.
								source_bool = hap1Src_bool;
								breakit += 1;
								break;
							}
						}
						if (mergedgraph.count(hap2Src_bool) > 0)
						{
							hap2srcData = mergedgraph[hap2Src_bool];
						}
						else
						{
							{
								if (mergedgraph.count(hap1Src_bool) > 0)
								{
									hap2Src_bool = hap1Src_bool;
									source_bool = hap1Src_bool;
									breakit += 1;
									break;
								}
								else
								{
									source_bool = hap1Src_bool;
									breakit += 1;
									break;
								}
							}
						}

						//ensure there's only a single sink from each hap source. More than sink
						//suggests one must be an error or perhaps the segmental duplication.
						//I'd go with the former and simply clip the 'allele' with lower cov.
						if (hap1srcData.size() == 1)
						{
							std::unordered_map<std::vector<bool>, SourceAttrib>::iterator it = hap1srcData.begin();
							Sink1 = it->first;
						}
						if (hap2srcData.size() == 1)
						{
							std::unordered_map<std::vector<bool>, SourceAttrib>::iterator it = hap2srcData.begin();
							Sink2 = it->first;
						}

						if (hap1srcData.size() > 1)
						{
							int maxcov = 0;
							for (std::unordered_map<std::vector<bool>, SourceAttrib>::iterator it = hap1srcData.begin(); it != hap1srcData.end(); ++it)
							{
								SourceAttrib sinkData = it->second;
								std::vector<bool> sink = it->first;
								if (sinkData.depth > maxcov)
								{
									maxcov = sinkData.depth;
									Sink1 = sink;
								}
								else //clip sink with lower than max allele cov 
								{
									mergedgraph[hap1Src_bool].erase(sink);

								}
							}
						}
						if (hap2srcData.size() > 1)
						{
							int maxcov = 0;
							for (std::unordered_map<std::vector<bool>, SourceAttrib>::iterator it = hap2srcData.begin(); it != hap2srcData.end(); ++it)
							{
								SourceAttrib sinkData = it->second;
								std::vector<bool> sink = it->first;
								if (sinkData.depth > maxcov)
								{
									maxcov = sinkData.depth;
									Sink2 = sink;
								}
								else //clip sink with lower than max allele cov 
								{
									mergedgraph[hap2Src_bool].erase(sink);

								}
							}
						}

						std::string tru("true");
						std::string fals("false");
						if (Sink1[0] == true)
						{
							boolstr1 = "true";
						}
						else
						{
							boolstr1 = "false";
						}
						if (Sink1[1] == true)
						{
							boolstr1 = boolstr1 + "_" + "true";
						}
						else
						{
							boolstr1 = boolstr1 + "_" + "false";
						}

						if (Sink2[0] == true)
						{
							boolstr2 = tru;
						}
						else
						{
							boolstr2 = fals;
						}
						if (Sink2[1] == true)
						{
							boolstr2 = boolstr2 + "_" + "true";
						}
						else
						{
							boolstr2 = boolstr2 + "_" + "false";
						}
						Snk1.push_back(bool2Char[boolstr1]);
						Snk2.push_back(bool2Char[boolstr2]);

						if (hap1_read1 > 0 || hap2_read2 > 0)
						{
							hap1_contigs += Snk1; //'hap1_read1' > 0 means 'reads1', the read set overlapng sink1,
												  //also contains reads found read set of hap1.
							hap2_contigs += Snk2;

							std::vector<bool>::iterator src_itr_beg = hap1Src_bool.begin();
							hap1Src_bool.erase(src_itr_beg, src_itr_beg + 2);
							hap1Src_bool.insert(hap1Src_bool.end(), Sink1.begin(), Sink1.end());
							src_itr_beg = hap2Src_bool.begin();
							hap2Src_bool.erase(src_itr_beg, src_itr_beg + 2);
							hap2Src_bool.insert(hap2Src_bool.end(), Sink2.begin(), Sink2.end());

							SourceAttrib src_attrib_1 = mergedgraph[source_bool][Sink1];
							SourceAttrib src_attrib_2 = mergedgraph[source_bool][Sink2];
							std::map<std::string, int> read_data1 = src_attrib_1.read_id;
							std::map<std::string, int> read_data2 = src_attrib_2.read_id;
							hap1BubbleReads.insert(read_data1.begin(), read_data1.end());
							hap2BubbleReads.insert(read_data2.begin(), read_data2.end());
							mergedgraph.erase(source_bool);
							source_bool = hap1Src_bool;
							continue;
						}
						if (hap1_read2 > 0 || hap2_read1 > 0)
						{
							hap1_contigs += Snk2;
							hap2_contigs += Snk1;

							std::vector<bool>::iterator src_itr_beg = hap1Src_bool.begin();
							hap1Src_bool.erase(src_itr_beg, src_itr_beg + 2);
							hap1Src_bool.insert(hap1Src_bool.end(), Sink2.begin(), Sink2.end());
							src_itr_beg = hap2Src_bool.begin();
							hap2Src_bool.erase(src_itr_beg, src_itr_beg + 2);
							hap2Src_bool.insert(hap2Src_bool.end(), Sink1.begin(), Sink1.end());

							SourceAttrib src_attrib_1 = mergedgraph[source_bool][Sink1];
							SourceAttrib src_attrib_2 = mergedgraph[source_bool][Sink2];
							std::map<std::string, int> read_data1 = src_attrib_1.read_id;
							std::map<std::string, int> read_data2 = src_attrib_2.read_id;
							hap1BubbleReads.insert(read_data2.begin(), read_data2.end());
							hap2BubbleReads.insert(read_data1.begin(), read_data1.end());
							mergedgraph.erase(source_bool);
							source_bool = hap1Src_bool;
							continue;
						}
					}
												
					//after exiting bubble loop, impt to check if the loop was completed naturally or if there was 
					//an abrupt exit. Also, /full sinks of 'hap1_source' (hap1_source+'sink_1') and 'hap2_source' 
					//(hap2_source+'sink_2') are expected to be same after completion of bubble loop. 
					//If this isn't the case, then there must have been an extra wihtin-bubble SNP unaccounted for. 
					//Best to terminate contig chain here.
					if (breakit > 0)
					{
						breakit = 0;
						std::pair<std::string, std::string> hap1hap2;
						hap1hap2.first = hap1_contigs;
						hap1hap2.second = hap2_contigs;
						contigs.push_back(hap1hap2);
						bubblecnt = 0; //re-intialise in readiness for building new contig
						hap1BubbleReads.clear(); //cleared for same reason as above
						hap2BubbleReads.clear();
						continue;
					}					
				}
				else
				{
					std::string Sinke("");
					std::string boolstr;
					
					if (sink[0] == true)
					{
						boolstr = "true";
					}
					else
					{
						boolstr = "false";
					}
					if (sink[1] == true)
					{
						boolstr = boolstr + "_" + "true";
					}
					else
					{
						boolstr = boolstr + "_" + "false";
					}
					Sinke.push_back(bool2Char[boolstr]);
					hap1_contigs += Sinke;
					hap2_contigs += Sinke;
					mergedgraph.erase(source_bool);

					std::vector<bool>::iterator src_it_beg = source_bool.begin();
					source_bool.erase(src_it_beg, src_it_beg + 2);
					source_bool.insert(source_bool.end(), sink.begin(), sink.end());
					continue;
				}
			}

						
			if (sources_of_src.size() < 2 && locuscov2modeRatio >= 2.0) //for suspect middle of repeat unit
			{
				potMiddleRepeat += 1;
				breakit += 1;
			}
	    }
		else //if source is not found, extending contig terminated and collated
		{
			std::pair<std::string, std::string> hap1hap2;
			hap1hap2.first = hap1_contigs;
			hap1hap2.second = hap2_contigs;
			contigs.push_back(hap1hap2); 	
			//hap contigs re-initialised
			hap1_contigs = std::string(""); 
			hap2_contigs = std::string("");

			if (mergedgraph.size() > 0)
			{
				//note that mergedgraph here is not of same size as one b4 contig building.
				hashMap::iterator src_it1 = mergedgraph.begin();

				std::unordered_map<std::vector<bool>, SourceAttrib> sourceData1 = src_it1->second;
				std::unordered_map<std::vector<bool>, SourceAttrib>::iterator sink_it1 = sourceData1.begin();
				std::unordered_map<std::vector<bool>, int> sources_of_source1 = sink_it1->second.src_of_src;

				//ensure there's just a single path from next source and that it's not within repeat
				int srcDataSz1 = sourceData1.size();
				int sourcesOfsrcSize = sources_of_source1.size();
				while (srcDataSz1 < 1 || sourcesOfsrcSize > 1) //loop also to continue when no sink is found
				{
					//srcDataSz1 > 1 indicates end of repeat in general or a bubble
					//srcDataSz1 < 1 for non-tandem duplicates ofr which sinks are removed
					//sourcesOfsrcSize > 1 indicates start of repeat
					if (srcDataSz1 < 1) //if source exists w/o sink connexn (cos removd e.g @ duplicate segments), delete source
					{
						std::vector<bool> src = src_it1->first;
						mergedgraph.erase(src);
					}
					++src_it1;
					sourceData1 = src_it1->second;
					sink_it1 = sourceData1.begin();
					srcDataSz1 = sourceData1.size();
					sources_of_source1 = sink_it1->second.src_of_src;
					sourcesOfsrcSize = sources_of_source1.size();
				}
				if (src_it1 == mergedgraph.end()) //if a non-repeat and non-bubble source can't be found
					//exit the contig-building while loop. this ensures contig-building loop can still be exited
					//even though kmers source data within duplicate segments are not deleted (not deleted so 
					//as to allow multi-pass by trans-located versions of the segment)
				{
					break;
				}
				source_bool = src_it1->first;
			}
			else
			{
				break;
			}
		}
	}
	return contigs;
}


int main(int argc, char *argv[]) {
	if (argc < 2)
	{
		std::cout << "usage: " << argv[0] << "-k WordSize" << "-n number of processes" << "<FastqFile1>" << "<FastqFile2>" << std::endl;
		return 1;
	}
	//below, data stxs are devised to allow mapping seq char to a bool representation and vice-versa
	CharBoolMap char2bool;
	std::pair<bool, bool> falsefalse;
	std::pair<bool, bool> falsetru;
	std::pair<bool, bool> trufalse;
	std::pair<bool, bool> trutru;

	falsefalse.first = false;
	falsefalse.second = false;
	falsetru.first = false;
	falsetru.second = true;
	trufalse.first = true;
	trufalse.second = false;
	trutru.first = true;
	trutru.second = true;
	char2bool['A'] = falsefalse;
	char2bool['C'] = falsetru;
	char2bool['G'] = trufalse;
	char2bool['T'] = trutru;

	BoolCharMap bool2char;
	bool2char["true_true"] = 'T';
	bool2char["true_false"] = 'G';
	bool2char["false_true"] = 'C';
	bool2char["false_false"] = 'A';

	std::string fastqline;
	std::vector<std::string> fastqdata;

	hashMap dbGraph;
	int defaultProcs = 1; //default num of processes
	int defaultk = 31;
	int singlefastq = 0;

	std::string fastqname1_1;
	std::string fastqname2_1;

	if (argc == 5)
	{
		if (argv[1][1] == 'k' && argv[2][1] == 'n')
		{
			std::string ksize(argv[1]);
			ksize = ksize.substr(2);
			defaultk = std::stoi(ksize);

			std::string nprocs(argv[2]);
			nprocs = nprocs.substr(2);
			defaultProcs = std::stoi(nprocs);
		}
		else if (argv[1][1] == 'n' && argv[2][1] == 'k')
		{
			std::string nprocs(argv[1]);
			nprocs = nprocs.substr(2);
			defaultProcs = std::stoi(nprocs);

			std::string ksize(argv[2]);
			ksize = ksize.substr(2);
			defaultk = std::stoi(ksize);
		}
		else
		{
			std::cout << "You entered 5 arguments. The second and third must be k size and number of processes in no order\n";
			return 1;
		}
		std::ifstream fastqfile1(argv[3]);
		std::ifstream fastqfile2(argv[4]);
		std::string fastqname1(argv[3]);
		std::string fastqname2(argv[4]);
		fastqname1_1 = fastqname1;
		fastqname2_1 = fastqname2;

		int numOflines1 = fileLineCounter(fastqname1);
		int numOflines2 = fileLineCounter(fastqname2);
		
		fastqfilebreak(defaultProcs, numOflines1, fastqname1);
		fastqfilebreak(defaultProcs, numOflines2, fastqname2);
	}

	if (argc == 4)
	{
		if (argv[1][0] == '-' && argv[1][1] == 'k')
		{
			std::string ksize(argv[1]);
			ksize = ksize.substr(2);
			std::cout << "ksize: " << ksize << "\n";
			defaultk = std::stoi(ksize);
			if (argv[2][0] != '-' && argv[2][1] != 'n')
			{
				std::ifstream fastqfile1(argv[2]);
				std::ifstream fastqfile2(argv[3]);
				std::string fastqname1(argv[2]);
				std::string fastqname2(argv[3]);
				fastqname1_1 = fastqname1;
				fastqname2_1 = fastqname2;

				int numOflines1 = fileLineCounter(fastqname1);
				int numOflines2 = fileLineCounter(fastqname2);				
				fastqfilebreak(defaultProcs, numOflines1, fastqname1);
				fastqfilebreak(defaultProcs, numOflines2, fastqname2);
			}
			if (argv[2][0] == '-' && argv[2][1] == 'n')
			{
				std::ifstream fastqfile1(argv[3]);
				std::string nprocs(argv[2]);
				nprocs = nprocs.substr(2);
				defaultProcs = std::stoi(nprocs);
				std::string fastqname1(argv[3]);
				fastqname1_1 = fastqname1;
				
				int numOflines1 = fileLineCounter(fastqname1);
				
				fastqfilebreak(defaultProcs, numOflines1, fastqname1);
				singlefastq += 1;
			}
			else
			{
				std::cout << "four arguments supplied. At least one must be kmer or number of processor parameter\n";
				return 1;
			}
		}
		else
		{
			if (argv[1][0] == '-' && argv[1][1] == 'n')
			{
				std::string nprocs(argv[2]);
				nprocs = nprocs.substr(2);
				defaultProcs = std::stoi(nprocs);

				if (argv[2][0] == '-' && argv[2][1] == 'k')
				{
					std::ifstream fastqfile1(argv[3]);
					std::string ksize(argv[2]);
					ksize = ksize.substr(2);
					defaultk = std::stoi(ksize);
					std::string fastqname1(argv[3]);
					fastqname1_1 = fastqname1;
					
					int numOflines1 = fileLineCounter(fastqname1);
					
					fastqfilebreak(defaultProcs, numOflines1, fastqname1);
					singlefastq += 1;
				}
				else
				{
					std::ifstream fastqfile1(argv[2]);
					std::ifstream fastqfile2(argv[3]);
					std::string fastqname1(argv[2]);
					std::string fastqname2(argv[3]);
					fastqname1_1 = fastqname1;
					fastqname2_1 = fastqname2;

					int numOflines1 = fileLineCounter(fastqname1);
					int numOflines2 = fileLineCounter(fastqname2);
					
					fastqfilebreak(defaultProcs, numOflines1, fastqname1);
					fastqfilebreak(defaultProcs, numOflines2, fastqname2);
				}
			}
			else
			{
				std::cout << "four arguments supplied. At least one must be kmer or number of processor parameter\n";
				return 1;
			}
		}
	}

	if (argc == 3)
	{
		if (argv[1][0] != '-')
		{
			std::ifstream fastqfile1(argv[1]);
			std::ifstream fastqfile2(argv[2]);
			std::string fastqname1(argv[1]);
			std::string fastqname2(argv[2]);
			fastqname1_1 = fastqname1;
			fastqname2_1 = fastqname2;

			int numOflines1 = fileLineCounter(fastqname1);
			int numOflines2 = fileLineCounter(fastqname2);
			
			fastqfilebreak(defaultProcs, numOflines1, fastqname1);
			fastqfilebreak(defaultProcs, numOflines2, fastqname2);
		}
		if (argv[1][0] == '-' && argv[1][1] == 'k')
		{
			std::string ksize(argv[1]);
			ksize = ksize.substr(2);
			defaultk = std::stoi(ksize);
			std::ifstream fastqfile1(argv[2]);
			std::string fastqname1(argv[2]);
			fastqname1_1 = fastqname1;
			
			int numOflines1 = fileLineCounter(fastqname1);			
			fastqfilebreak(defaultProcs, numOflines1, fastqname1);
			singlefastq += 1;
		}
		if (argv[1][0] == '-' && argv[1][1] == 'n')
		{
			std::string nprocs(argv[1]);
			nprocs = nprocs.substr(2);
			defaultProcs = std::stoi(nprocs);
			std::ifstream fastqfile1(argv[2]);
			std::string fastqname1(argv[2]);
			fastqname1_1 = fastqname1;
			
			int numOflines1 = fileLineCounter(fastqname1);
			fastqfilebreak(defaultProcs, numOflines1, fastqname1);
			singlefastq += 1;
		}
		else
		{
			std::cout << "arguments have been entered incorrectly\n";
			return 1;
		}
	}
	if (argc == 2)
	{
		if ((argv[1][0] == '-' && argv[1][1] == 'k') || (argv[1][0] == '-' && argv[1][1] == 'k'))
		{
			std::cout << "please provide at least one fastqfile\n";
			return 1;
		}
		if (argv[1][0] != '-')
		{
			std::ifstream fastqfile1(argv[1]);
			std::string fastqname1(argv[1]);
			fastqname1_1 = fastqname1;
			
			int numOflines1 = fileLineCounter(fastqname1);
			fastqfilebreak(defaultProcs, numOflines1, fastqname1);
			singlefastq += 1;
		}
	}

	std::vector<std::future<hashMap>> graphFutures;
	for (int i = 1; i <= defaultProcs; ++i)
	{
		std::string idx = std::to_string(i);
		std::cout << "Commencing building of graph " << i << "\n";
		if (singlefastq > 0)
		{
			std::string fileName = fastqname1_1 + "_" + idx;
			std::ifstream fastqFile(fileName.c_str());
			graphFutures.emplace_back(std::async(std::launch::async, std::move(graphBuild), std::move(fileName), std::move(defaultk), std::move(char2bool), std::move(bool2char))); //will start a new process each time its run
		}
		else
		{
			std::string fileName1 = fastqname1_1 + "_" + idx;
			std::string fileName2 = fastqname2_1 + "_" + idx;
			std::ifstream fastqFile1(fileName1.c_str());
			std::ifstream fastqFile2(fileName2.c_str());
			std::future<hashMap> graph1 = std::async(std::launch::async, std::move(graphBuild), std::move(fileName1), std::move(defaultk), std::move(char2bool), std::move(bool2char));
			std::future<hashMap> graph2 = std::async(std::launch::async, std::move(graphBuild), std::move(fileName2), std::move(defaultk), std::move(char2bool), std::move(bool2char));
			graphFutures.push_back(std::move(graph1));
			graphFutures.push_back(std::move(graph2));
		}
	}
	std::cout << "Merging assembly graphs\n";
	hashMap mergedGraph = graphMerger(graphFutures);
	std::cout << "Estimating sequencing depth\n";
	int modalcoverage = modalCovCalc(mergedGraph);
	std::cout << "merged graph size: " << mergedGraph.size() << "\n";

	std::cout << "Generating contigs\n";
	std::vector<std::pair<std::string, std::string>> contigData = contigGen(mergedGraph, modalcoverage, bool2char);
	std::ofstream contig1_out("Contig1_seqs.fasta");
	std::ofstream contig2_out("Contig2_seqs.fasta");
	int seqcnt = 1;
	for (auto hapPair : contigData)
	{
		std::string contig1 = hapPair.first;
		std::string contig2 = hapPair.second;
		contig1_out << "contig1_" << seqcnt << "\n";
		contig1_out << contig1 << "\n";
		contig2_out << "contig2_" << seqcnt << "\n";
		contig2_out << contig2 << "\n";
		seqcnt += 1;
	}

	return 0;
}



	


