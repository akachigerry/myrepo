#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <map>
#include <sstream>		//allows to use istringstream

int main(int argc, char *argv[]) {
	if (argc != 3)
	{
		std::cout << "usage: " << argv[0] << "<parent ids>" << "<vcf file>" << std::endl;
	}
	//read in parent ids and vcf file
	std::ifstream parentsFile(argv[1]);
	std::ifstream vcfFile(argv[2]);
	std::ofstream vcfOutFile("filtered_variants.vcf");
	std::string varline;
	std::string parent;
	std::vector<std::string> parentList;
	//this initialises a nested map of string key and map value which
	//in turn has a string key and a vector as value
	//declares a vector of vectors to hold portion of vcf we want
	std::map<std::string, std::map<std::string, std::vector<std::string> > > vcfMap;
	std::vector<std::vector<std::string> > vcf4out;

	if (!parentsFile.is_open() || !vcfFile.is_open())
	{
		std::cout << "Could not open file " << argv[1] << "or " << argv[1] << "or both " << std::endl;
		return 1;
	}
	else
	{
		while (parentsFile >> parent)
		{
			if (!parentsFile.eof() && parentsFile.fail())
			{
				std::cout << "Reading " << argv[0] << "failed\n";
				return 1;
			}
			if (parentsFile.eof() && !parentsFile.fail())
			{
				break;
			}
			parentList.push_back(parent);
		}
		//vector to record index positions of our parents of interest in vcf column header list
		std::vector<int> index;
		while (std::getline(vcfFile, varline))
		{
			if (!vcfFile.eof() && vcfFile.fail())
			{
				std::cout << "Reading " << argv[1] << "failed\n";
				return 1;
			}
			if (vcfFile.eof() && !vcfFile.fail())
			{
				break;
			}

			std::map<std::string, std::vector<std::string> > innerVcfMap;
			std::string wordStr;
			//if "#CHROM" in variant line
			if (varline.find("#CHROM") != std::string::npos)
			{
				std::vector<std::string> innerVcf4out;
				std::istringstream varWordStream(varline);

				//create another stream, this time not of the opened file
				//but of each line of vcf file initially read by getline
				int count = 0;					//initialise a count variable for use in tracking index number of parent in vcf column header list.	
				while (varWordStream >> wordStr)
				{
					if (count < 9)				//condition to ensure addition of first 9 columns of column title line to vector and nothing more.
					{
						innerVcf4out.push_back(wordStr);
						//if wordStr in parentList, then add it to innerVcf4out array
					}
					std::vector<std::string>::const_iterator it = std::find(parentList.begin(), parentList.end(), wordStr);
					if (it != parentList.end())
					{
						innerVcf4out.push_back(wordStr);
						std::cout << count << "\n";
						index.push_back(count);
					}
					count += 1;
				}
				vcf4out.push_back(innerVcf4out);	
			}
			if ((varline.find("Salix_viminalis_TGAC_ROT_v2_scaffold") != std::string::npos) && (varline.find("#") == std::string::npos))
			{
				std::vector<std::string> innerVcf4out;
				std::istringstream varWordStream(varline);
				int count = 0;
				while (varWordStream >> wordStr)
				{
					if (count < 9)
					{
						innerVcf4out.push_back(wordStr);
					}
					else
					{
						if (std::find(index.begin(), index.end(), count) != index.end())
						{
							innerVcf4out.push_back(wordStr);
						}
					}
					count += 1;
				}
				vcf4out.push_back(innerVcf4out);
			}
		}
		


		if (!vcfOutFile.is_open())
		{
			std::cout << "Could not open output vcf file " << std::endl;
			return 1;
		}

		for (int i; i < vcf4out.size(); i++)
		{
			int count = 0;
			for (std::vector<std::string>::const_iterator it = vcf4out[i].begin(); it != vcf4out[i].end(); ++it)
			{
				//before getting to last item of the inner vector, vcf4output[i], separate 
				//words (strings) with tab. After the last word, add a new line.
				if (count == vcf4out[i].size() - 1)
				{
					vcfOutFile << *it << "\n";
				}
				else 
				{
					vcfOutFile << *it << "\t";
				}
				count+=1;
			}
		}
	}
	return 0;
}
