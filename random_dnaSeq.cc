#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdlib>
#include <ctime>


int randRange(int low, int high)
{
	//get a random number, get it to be between 0 and difference between
	//high and low, then add the lowest possible value.
	return rand() % (high - low) + low;
}

int get_filesize(std::string filename)				//function used to get size of file
	{
		std::streampos begin, end;
		std::ifstream partialFile(filename, std::ios::binary);
		begin = partialFile.tellg();
		partialFile.seekg(0, std::ios::end);
		end = partialFile.tellg();
		int filesize = end - begin;
		partialFile.close();
		return filesize;
	}
		
int main(int argc, char *argv[]) {
	if (argc != 2)
	{
		std::cout << "usage: " << argv[0] << "<filename>" << std::endl;
	}

	//Assume argv[1] is a filename to open
	std::ifstream file(argv[1]);
	//Always check to see if file opening succeeded
	if (!file.is_open())
	{
		std::cout << "Could not open file " << argv[1] << std::endl;
		return 1;
	}
	std::string line;
	std::vector<std::string> seqdata;
	std::vector<std::string> seqtitle;
	//Use an infinite loop to iteratively read in each line of file.
	//Break once end of file is reached.
	//Use while(file<<line) which stops automatically once
	//end of file is reached. 
	
	while (file >> line)
	{
		//if read failure occurs before end of file, exit as this indicates
		//something faulty with file reading process
		if (file.fail() && !file.eof())
		{
			std::cout << "Faulty file input " << "\n";
			return 0;
		}
		//if only end of file is reached, and there's no file read failure,
		//simply break out from while loop.
		if (!file.fail() && file.eof())
		{
			break;
		}
		
		if (line[0] == '>')
		{
			seqtitle.push_back(line.substr(1, line.size() - 1));		//add substring of sequence tile excluding the 
																		//'<' character to back of vector
		}
		else
		{
			seqdata.push_back(line);
		}
	}
	
	srand(1);		//sets random seed to current time
					//keep generating random sequences until max seq threshold is reached

	std::ofstream outstream1("randomseqFile1.fasta");
	std::ofstream outstream2("randomseqFile2.fasta");
	std::ofstream outstream3("randomseqFile3.fasta");
	
	if (!outstream1.is_open())
		{
			std::cout << "Could not open " << "randomseqFile1.fasta" << std::endl;
			return 1;
		}
	if (!outstream2.is_open())
		{
			std::cout << "Could not open " << "randomseqFile2.fasta" << std::endl;
			return 1;
		}
	if (!outstream3.is_open())
		{
			std::cout << "Could not open " << "randomseqFile3.fasta" << std::endl;
			return 1;
		}
	
	int randomIndex;
	for (int i = 0; i<seqdata.size(); i++)
	{
		randomIndex = randRange(0, seqdata.size() - 1);
		std::cout << randomIndex << "\n";
		outstream1 << ">" << seqtitle[randomIndex] << "\n";
		outstream1 << seqdata[randomIndex] << "\n";
		
		if (get_filesize("randomseqFile1.fasta") >= 80000000)
		{
			outstream1.close();
		}
		
		randomIndex = randRange(0, seqdata.size() - 1);
		std::cout << randomIndex << "\n";
		outstream2 << ">" << seqtitle[randomIndex] << "\n";
		outstream2 << seqdata[randomIndex] << "\n";
		
		if (get_filesize("randomseqFile2.fasta") >= 80000000)
		{
			outstream2.close();
		}

		randomIndex = randRange(0, seqdata.size() - 1);
		std::cout << randomIndex << "\n";
		outstream3 << ">" << seqtitle[randomIndex] << "\n";
		outstream3 << seqdata[randomIndex] << "\n";
		
		if (get_filesize("randomseqFile3.fasta") >= 80000000)
		{
			outstream3.close();
			break;
		}
	}

	return 0;
}
