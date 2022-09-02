#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdlib>
#include <map>
#include <unordered_map>
#include <sstream>		//allows to use istringstream
#include <stdlib.h> 	//for atoi
#include <stdio.h>
#include <cstdio> //for open
#include <string.h>
#include <utility>		//for std::pair
#include <numeric>
#include <csignal>
#include <signal.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/unistd.h>
#include <memory>
#include <sys/mman.h>  /*For mmap*/
#include <json.hpp>


using json = nlohmann::json;

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


void attachShm (key_t ShmKey, const char* Shmptr, std::string shmName)
{
	int shmid = shmget(ShmKey, 0, 0); //opens existing shm for read access only
	shmError(shmid);
	std::cout << "Existing shared memory opened\n";
	//attach shm to processes data space

	Shmptr = (const char*)shmat(shmid, 0, 0);
	if (Shmptr == (const char *)-1)
	{
		std::cout << "shmat(): " << shmName << "\"Shmptr\" failed";
		exit(EXIT_FAILURE);
	}
	std::cout << shmName << " shared memory attached\n";
}

void shmError(int shm)
{
	if (shm < 0) //if shmget generates error
	{
		perror("shm failed to open");
	}
	if (errno == EACCES)
	{
		perror("access to shared memory denied");
		exit(EXIT_FAILURE);
	}
	if (errno == EEXIST)
	{
		perror("Shared memory already exists");
		exit(EXIT_FAILURE);
	}
	if (errno == EINVAL)
	{
		perror("Shared memory outside allowed range");
		exit(EXIT_FAILURE);
	}
	if (errno == ENFILE)
	{
		perror("Too many shared memory opened");
		exit(EXIT_FAILURE);
	}
	if (errno == ENOENT)
	{
		perror("IPC_CREAT not set and named memory segement doesn't exist");
		exit(EXIT_FAILURE);
	}
	if (errno == ENOMEM)
	{
		perror("No memory could be allocated for segment overhead");
		exit(EXIT_FAILURE);
	}
}


int main(int argc, char *argv[])
{
	
	if (argc < 4) //4 is min number of args
	{
		std::cout << "insufficient number of arguments supplied\n";
		std::cout << "usage: " << argv[0] << "-v" << "<vcf db file>" << "-a" << "<annot db file>" << "-s" << "<seq db file>" << "<region bed file>" << "\n";
		return 1;
	}
	
	if (!(argv[1][0] == "-" && argv[1][1] == "v"))
	{
			std::cout << "The second argument must be vcf db flag\n;
			return 1;
	}
	else
	{
		/*client side read of shared memory (shm) contents.*/
		//First, open shm
		key_t vcfShmKey = 0;
		key_t annotShmKey = 0;
		key_t seqShmKey = 0;
		int vcfshmid = 0;
		int annotshmid = 0;
		int seqshmid = 0;
		const char * vcfShmptr;
		const char * annotShmptr;
		const char * seqShmptr;
		std::string vcfShmName("Vcf");
		std::string annotShmName("Annotion");
		std::string seqShmName("Ref sequence");
		
		if ((fd = open(argv[2], O_CREAT | O_RDWR, 0666)) < 0) {
			perror("Cannot open vcf file");
			exit(1);
		}
		//determine file size
		std::ifstream(argv[2], std::ios::binary);
		jsonInput.seekg(0, jsonInput.end); //moves file pointer (needle) from beginning to last byte
		int jsonLength = jsonInput.tellg() + 1; //+1 to allow space for trailing '\0' char
		jsonInput.seekg(0, jsonInput.beg); //move needle back to begin of file
		//ftruncate. Would add null byte ('\0') to end of file cos of extra byte added in file length above
		if (ftruncate(fd, jsonLength) == -1)
		{
			perror("ftruncate failure with error code -1");
			exit(1);
		}
		if ((vcfShmptr = (const char*)mmap(0, jsonLength, PROT_READ, MAP_SHARED, fd, 0)) == (void *)-1)
		{
			perror("mmap failure");
			exit(5);

		vcfShmKey = 1234;
		attachShm(vcfShmKey, vcfShmptr, vcfShmName);
		
		if (argc == 5 || argc == 7)
		{
			if ((argv[3][0] == "-" && argv[3][1] == "a") && !(argv[3][0] == "-" && argv[3][1] == "s"))
			{
				annotShmKey = 1236; //annot db shm key
				attachShm(annotShmKey, annotShmptr, annotShmName);
			}
			if ((argv[3][0] == "-" && argv[3][1] == "s") && !(argv[3][0] == "-" && argv[3][1] == "a"))
			{
				seqShmKey = 1235; //ref seq db shm key
				attachShm(seqShmKey, seqShmptr, seqShmName);
			}
			if (argc == 7)
			{
				if (argv[5][0] != "-" || argv[5][1] != "s")
				{
					std::cout << "Please specify refseq flag before refseq argument\n";
					return 1;
				}
				if (argv[5][0] == "-" || argv[5][1] == "s")
				{
					seqShmKey = 1235; //ref seq db shm key
					attachShm(seqShmKey, seqShmptr, seqShmName);
				}
			}
		}
		if (argc == 6 || argc==8) //a genomic region file must be included here
		{
			std::string arg4(argv[3]);
			std::string arg6(argv[5]);
			std::string arg8(argv[7]);
			if (arg4.find(".bed") == std::string::npos && arg6.find(".bed") == std::string::npos)
			{
				if (argc == 8)
				{
					if (arg8.find(".bed") == std::string::npos)
					{
						std::cout << "The region bed file must be 4th, 6th or 8th argument \n";
						return 1;
					}
					else
					{
						if (argv[3][0] == "-" && argv[3][1] == "s")
						{
							seqShmKey = 1235; //ref seq db shm key
							attachShm(seqShmKey, seqShmptr, seqShmName);
							if (argv[6][0] == "-" && argv[6][1] == "a")
							{
								annotShmKey = 1236; //ref seq db shm key
								attachShm(annotShmKey, annotShmptr, annotShmName);
							}
						}
						if (argv[3][0] == "-" && argv[3][1] == "a")
						{
							annotShmKey = 1236; //ref seq db shm key
							attachShm(annotShmKey, annotShmptr, annotShmName);
							if (argv[6][0] == "-" && argv[6][1] == "s")
							{
								seqShmKey = 1235; //ref seq db shm key
								attachShm(seqShmKey, seqShmptr, seqShmName);
							}
						}
					}
				}
				else
				{
					std::cout << "The region bed file must be 4th, 6th or 8th argument \n";
					return 1;
				}
				
			}
			if (arg4.find(".bed") != std::string::npos)
			{
				if (argv[4][0] == "-" && argv[4][1] == "s")
				{
					seqShmKey = 1235; //ref seq db shm key
					attachShm(seqShmKey, seqShmptr, seqShmName);
					if (argc == 8 && argv[6][0] == "-" && argv[6][1] == "a")
					{
						annotShmKey = 1236;
						attachShm(annotShmKey, annotShmptr, annotShmName);
					}
					if (argc == 8 && argv[6][0] != "-" && argv[6][1] != "s")
					{
						std::cout << "Your 6th argument must be a flag to annotation database\n";
					}
				}
				if (argv[4][0] == "-" && argv[4][1] == "a")
				{
					annotShmKey = 1236; //annot db shm key
					attachShm(annotShmKey, annotShmptr, annotShmName);
					if (argc == 8 && argv[6][0] == "-" && argv[6][1] == "s")
					{
						attachShm(seqShmKey, seqShmptr, seqShmName);
					}
					if (argc == 8 && argv[6][0] != "-" && argv[6][1] != "s")
					{
						std::cout << "Your 6th argument must be a flag to ref sequence database\n";
					}
				}
				
			}
		}
		

		//Then confirm it's truly written to (not empty)

	}
		
	