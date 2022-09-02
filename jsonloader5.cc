#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <ctype.h>	/*requried for tolower*/
#include <csignal>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <pthread.h> 
#include <semaphore.h> 
#include <unistd.h>
#include <fcntl.h>           /* For O_* constants */
#include <sys/stat.h>
#include <sys/mman.h>
#include <errno.h>


#define SEM_FILE_PATH "/newsem"

std::locale loc;	//a locale object encapsulates a set of features that are culture-specific, 
					//to enhance international portability. A default-constructed locale object is a global locale.
void lowerCase(std::string str)
{
	int str_length = str.length();
	for (std::string::size_type i = 0; i < str_length; ++i)
	{
		std::tolower(str[i]);
	}
}

std::vector<std::string> split(char const *s, char delimiter)
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


int main(int argc, char *argv[]) {

	std::cout << "Loading databases...please wait" << "\n";
	if (argc != 3)
	{
		std::cout << "usage: " << argv[0] << "-[v,r,g] for vcf, reference, annotation(gff/gtf) files " << "<json file>" << "\n";
	}

	key_t shmKey;
	std::string fileTypeIndicator(argv[1]);
	if (fileTypeIndicator[1] == 'v')
	{
		shmKey = 1234;
	}
	if (fileTypeIndicator[1] == 'r')
	{
		shmKey = 1235;
	}
	if (fileTypeIndicator[1] == 'g')
	{
		shmKey = 1236;
	}

	std::cout << "shmKey: " << shmKey << "\n";
	
	std::ifstream jsonInput(argv[2], std::ios::binary); 
	jsonInput.seekg(0, jsonInput.end); //moves file pointer (needle) from beginning to last byte
	int jsonLength = jsonInput.tellg(); 
	jsonInput.seekg(0, jsonInput.beg); //move needle back to begin of file cos file would be eventually read.
	
	//create the semaphore "/mysems1w" in locked state.
	//this is important prior to creating, attaching and writing to shared memory.
	//semaphore'll be unlocked below after write process.
	//other (read) processes will be required to check if semaphore is unlocked b4 use
	//O_CREAT | O_EXCL ensures error is returned if a semaphore with given name already exists.
	sem_t *mysemp;
	//this semaphore has value 0.
	mysemp = sem_open(SEM_FILE_PATH, O_CREAT, 0666, 0); //0666 means read/write by user and by grp and others

	if (mysemp == (sem_t *)-1)
	{
		perror("sem_open() \"mysemp\" failed ");
	}
	if (errno == EACCES)
	{
		perror("access to semaphore denied");
		exit(EXIT_FAILURE);
	}
	if (errno == EEXIST)
	{
		perror("Semaphore already exists");
		exit(EXIT_FAILURE);
	}
	if (errno == EINVAL)
	{
		perror("sem_open() not supported for given name");
		exit(EXIT_FAILURE);
	}
	if (errno == ENFILE)
	{
		perror("Too many semaphores currently open");
		exit(EXIT_FAILURE);
	}
	if (errno == ENOENT)
	{
		perror("O_CREAT not set and named semaphore doesn't exist");
		exit(EXIT_FAILURE);
	}
	std::cout << "seamphore opened" << "\n";

	size_t shmSz = (sizeof(char) * jsonLength) + 1000000; //extra 1MB added as caution against segfault
															 
	std::cout << "Shared memory size: " << shmSz << "\n";
	int shmid = shmget(shmKey, shmSz, 0644 | IPC_CREAT);//0644 means read/write by user and only read by grp and others
	shmError(shmid);
	std::cout << "shared memory opened\n";
	
	char *shmPtr, *shmWriter;
	//attach shm to processes data space
	shmPtr = (char*)shmat(shmid, 0, 0);
	if (shmPtr == (char *)-1)
	{
		perror("shmat(): \"shmPtr\" failed");
		exit(EXIT_FAILURE);
	}
	std::cout << "shared memory attached\n";
	shmWriter = shmPtr; /*'shmWriter' also now points to begin of shm segment.
						Will be moved along the segment to fill it with data but 'shmPtr'
						stays at the start*/
	/*Now read in json data on disk char by char into shm*/
	char c;
	while (jsonInput >> std::noskipws >> c)
	{
		*shmWriter++ = c;
	}

	//now that shared memory has been created, attached and written to, unlock semaphore
	//to allow reader processes access to shared mem.
	int sempost;
	sempost = sem_post(mysemp); //sem_post increments semaphore to 1.
	if (errno == EINVAL)
	{
		perror("semaphore is not valid\n");
		exit(EXIT_FAILURE);
	}
	else if (errno == EOVERFLOW)
	{
		perror("max semaphore value exceeded\n");
		exit(EXIT_FAILURE);
	}
	else if (sempost == 0)
	{
		printf("sem_post() succeeded!\n");
	}
	return 0;
}