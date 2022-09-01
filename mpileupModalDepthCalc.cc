#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>		//allows to use istringstream
#include <stdlib.h> 	//for atoi
#include <stdio.h>
#include <string.h>
#include <utility>		//for std::pair
#include <array>
#include <csignal>
#include <signal.h>
#include <map>
#include <locale>		//for std::tolower


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

std::string lowerCase(std::string str)
{
	int str_length = str.length();
	std::locale loc;
	for (std::string::size_type i = 0; i < str_length; ++i)
	{
		char base;
		base = std::tolower(str[i], loc);
		str[i] = base;
	}
	return str;
}

int main(int argc, char *argv[])
{
	std::vector<std::string> args(argv, argv + argc); /*converts array of cstrings into vector of strings*/
	std::vector<std::string>::reverse_iterator outputPtr = args.rbegin();
	std::string outputName = *outputPtr;
	std::vector<std::string>::reverse_iterator pileupPtr = args.rbegin() + 1;
	std::string pileupName = *pileupPtr;
	
	//if arg count < minimum or the parameter ('-') qualifier is found as part of pileup file name or output file name, 
	//quit with output of tool usage info
	if (argc < 4 || outputName.find('-') != std::string::npos || pileupName.find('-') != std::string::npos)
	{
		std::cout << "usage: " << argv[0] << "-D <min locus depth> " << "-p min alt allele proportion " << "-f <min alt allele freq " << \
			"<mpileup file>" << "<sample list in mpileup order>" << "<SNP matrix otput>" << std::endl;
		std::cout << "quiting now..\n";
		return 1;
	}
	if (argc == 7 && (args[1][1] != 'D' || args[2][1] != 'p' || args[3][1] != 'f'))
	{
		std::cout << "usage: " << argv[0] << "-D <min locus depth> " << "-p min alt allele proportion " << "-f <min alt allele freq " << \
			"<mpileup file>" << "<sample list in mpileup order>" << "<SNP matrix otput>" << std::endl;
		std::cout << "quiting now1..\n";
		return 1;
	}

	if (argc == 6 && ((args[1][1] != 'D'  && args[1][1] != 'p' && args[1][1] != 'f') || \
		(args[2][1] != 'D'  && args[2][1] != 'p' && args[2][1] != 'f')))
	{
		std::cout << "usage: " << argv[0] << "-D <min locus depth> " << "-p min alt allele proportion " << "-f <min alt allele freq " << \
			"<mpileup file>" << "<sample list in mpileup order>" << "<SNP matrix otput>" << std::endl;
		std::cout << "quiting now2..\n";
		return 1;
	}

	if (argc == 5 && (args[1][1] != 'D'  && args[1][1] != 'p' && args[1][1] != 'f'))
	{
		std::cout << "usage: " << argv[0] << "-D <min locus depth> " << "-p min alt allele proportion " << "-f <min alt allele freq " << \
			"<mpileup file>" << "<sample list in mpileup order>" << "<SNP matrix otput>" << std::endl;
		std::cout << "quiting now3..\n";
		return 1;
	}

	int defaultMinDepth = 12;
	std::string suppliedMinDepth;
	if (args[1][1] == 'D' && args[1].size() > 2)
	{
		suppliedMinDepth = args[1].substr(2);
		int suppliedMinDepth1 = std::stoi(suppliedMinDepth);
		defaultMinDepth = suppliedMinDepth1;
	}

	double defaultMinAltProp = 0.3;
	std::string suppliedMinAltProp;
	if ((args[2][1] == 'p' && args[2].size() > 2) || (args[1][1] == 'p' && args[1].size() > 2))
	{
		if (args[2][1] == 'p')
		{
			suppliedMinAltProp = args[2].substr(2);
		}
		if (args[1][1] == 'p')
		{
			suppliedMinAltProp = args[1].substr(2);
		}
		double suppliedMinAltProp1 = std::stod(suppliedMinAltProp);
		defaultMinAltProp = suppliedMinAltProp1;
	}

	double defaultMinAltfreq = 0.1;
	std::string suppliedMinAltfreq;
	std::string defaultMinAltProp1 = std::to_string(defaultMinAltProp);
	
	if ((args[3][1] == 'f' && args[3].size() > 2) || (args[2][1] == 'f' && args[2].size() > 2) \
			|| (args[1][1] == 'f' && args[1].size() > 2))
	{
		if (args[3][1] == 'f')
		{
			suppliedMinAltfreq = args[3].substr(2);
		}
		if (args[2][1] == 'f')
		{
			suppliedMinAltfreq = args[2].substr(2);
		}
		if (args[1][1] == 'f')
		{
			suppliedMinAltfreq = args[1].substr(2);
		}
		double suppliedMinAltfreq1 = std::stod(suppliedMinAltfreq);
		defaultMinAltfreq = suppliedMinAltfreq1;
	}

	std::ifstream mpileupFile(argv[argc-2]);
	std::ifstream sampleListFile(argv[argc-1]);

	std::cout << "Reading mpileup file...\n";
	std::string line;
	std::string sampleline;
	int linecnt = 0;

	std::vector<std::string> sampleListVec;
	std::map<std::string, std::map<std::string, std::string> > SNPmatrix;

	while (std::getline(sampleListFile, sampleline))
	{
		if (sampleListFile.fail() && !sampleListFile.eof())
		{
			std::cout << "sample list file could not be read " << "\n";
			return 0;
		}
		if (!sampleListFile.fail() && sampleListFile.eof())
		{
			break;
		}
		sampleListVec.push_back(sampleline);
		std::cout << "sample " << sampleline << "\n";
	}

	std::ofstream totalSnpCntOut("totalSnpCountData");
	std::map<std::string, int> totalSnpCnt;
	std::map<std::string, std::map<int, int>> locusDepthDist;
	while (std::getline(mpileupFile, line))
	{
		if (mpileupFile.fail() && !mpileupFile.eof())
		{
			std::cout << "mpileup file could not be read " << "\n";
			return 0;
		}
		if (!mpileupFile.fail() && mpileupFile.eof())
		{
			break;
		}
		std::cout << "mpileup line: " << linecnt << " read\n";
		linecnt += 1;

		std::vector<std::string> lineVec = split1(line.c_str(), '\t');
		

		int depthIndexStart = 3;
		int locusBasesIndexStart = 4;

		int lineVecSize = lineVec.size();
		std::vector<std::string> locusReadBases;
		int locusDepthPassed = 0;
		int SNPfound = 0;
		int sampleindex = 0;
		std::map <std::string, std::string > innerMap1;
		int depth = 0;
		std::string chrid = lineVec[0];
		std::string pos = lineVec[1];
		std::string coord = chrid + "$" + pos;
		std::string base = lineVec[2];

		std::string depth1;
		for (int i = 0; i < lineVecSize; ++i)
		{
			std::string parent = sampleListVec[sampleindex];

			if (i == depthIndexStart)
			{
				depth = std::stoi(lineVec[i]);
				depth1 = std::to_string(depth);
				std::string baseNdepth = base + "$" + depth1;

				if (depth < defaultMinDepth)
				{
					depthIndexStart += 3; //jumps to depthIndex of next sample
										  //if we've  gone beyond the 1st sample
					i += 1;
					continue;
				}
				
				if (SNPmatrix.count(coord) > 0)
				{
					SNPmatrix[coord][parent] = baseNdepth;
				}
				else
				{
					std::map<std::string, std::string> innerMap;
					innerMap[parent] = baseNdepth;
					SNPmatrix[coord] = innerMap;
				}

				depthIndexStart += 3; //jumps to depthIndex of next sample
				//if we've  gone beyond the 1st sample
				if (i > 3)
				{
					sampleindex += 1;
					parent = sampleListVec[sampleindex];
				}

				if (locusDepthDist.count(parent) > 0)
				{
					if (locusDepthDist[parent].count(depth) > 0 && depth > 0)
					{
						locusDepthDist[parent][depth] += 1;
						int freq = locusDepthDist[parent][depth];
					}
					else
					{
						locusDepthDist[parent][depth] = 1;
					}
				}
				else
				{
					std::map<int, int> innerMap;
					innerMap[depth] = 1;
					locusDepthDist[parent] = innerMap;
				}
			}

			if (i == locusBasesIndexStart)
			{
				i += 1;
				continue;
			}
		}
	}


	int targetSize = 0;
	
	std::ofstream depthOutput("perSnpPerParentDepth2");
	std::ofstream modaldepth("perParentModalDepth2");

	depthOutput << "Chrom\tPosition";
	
	for (auto parent : sampleListVec)
	{
		depthOutput << "\t" << parent;
		
	}
	depthOutput << "\n";
	


	int cnt = 0;
	for (auto parent : sampleListVec)
	{
		if (locusDepthDist.count(parent) > 0)
		{
			modaldepth << parent << "\n";
			modaldepth << "depth" << "\t" << "freq\n";
			std::map<int, int> parentDepthData = locusDepthDist[parent];
			int maxfreq = 0;
			int mode = 0;
			int depthRange = 5;
			std::vector<int> rangeStartVec;
			for (std::map<int, int>::iterator it1 = parentDepthData.begin(); it1 != parentDepthData.end(); ++it1)
			{
				
				int depth = it1->first;
				
				int freq = it1->second;
				modaldepth << depth << freq << "\n";
			}
			modaldepth << "\n";
		}
	}
	
	

	std::cout << "SNPmatrix size" << SNPmatrix.size() << "\n";
	std::string prevChrid("");
	std::cout << "sampleList size " << sampleListVec.size() << "\n";
	for (std::map<std::string, std::map<std::string, std::string> >::iterator iter = SNPmatrix.begin(); iter != SNPmatrix.end(); ++iter)
	{
		std::map<std::string, std::string > InnerMap = iter->second;
		std::string snp_coord = iter->first;
		
		std::vector<std::string> snpPosVec = split1(snp_coord.c_str(), '$');
		std::string chrid = snpPosVec[0];
		std::string snpPos = snpPosVec[1];
		int snpPos1 = std::stoi(snpPos);
		std::cout << "Inner Map size: " << InnerMap.size() << "\n";
		depthOutput << chrid << "\t" << snpPos;
		for (auto parent : sampleListVec)
		{
			if (InnerMap.count(parent) > 0)
			{
				std::vector<std::string> baseNdepth = split1((InnerMap[parent]).c_str(), '$');
				std::string base = baseNdepth.front();
				char base1 = base.front();
				std::cout << "base is: " << base << "\n";
				std::cout << "base1 is: " << base1 << "\n";

				std::string depth = baseNdepth[1];
				std::cout << "depth is: " << depth << "\n";
				int depth1 = std::stoi(depth);
				depthOutput << "\t" << depth1;
			}
			else
			{
				depthOutput << "\t" << ".";
			}
		}
		depthOutput << "\n";
	}
	return 0;
}

