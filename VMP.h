#ifndef _VMP_INC_
#define _VMP_INC_

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <boost/algorithm/string.hpp>
#include <array>
#include "platform.h"

using namespace std;

struct VirtualMachine
{
	int RAMReq;
	int CPUReq;
	int pmAssigned;
	int expiry;
	int id;
	int listPos;
	
	float resourceBalance;
	bool smallCPU = false;
	bool smallRAM = false;

	bool operator==(const VirtualMachine & other)
	{
		if(this->id == other.id)
		{
			return true;
		}

		else
			return false;
	}
};

struct PhysicalMachine
{
	int ID;
	int RAMCapacity;
	int CPUCapacity;
	int CPUUsage = 0;
	int RAMUsage= 0;
};

class VMP
{
	public:
	VirtualMachine* vmList;
	PhysicalMachine* pmList;
	int* CPUCapacity;
	int* RAMCapacity;
	int numPMs;
	int numVMs;
	int numVMsPadded;
	int numPMPadded;
	vector<int> PMCPU;
	vector<int> PMRAM;
	vector<int> subNums;
	bool multipleServerTypes;
	float* CPUReq;
	float* RAMReq;
	int CPUReqTotal;
	int RAMReqTotal;
	int aRAM;
	int aCPU;
	int bRAM;
	int bCPU;

	int CPUCapAv = 0;
	int RAMCapAv = 0;
	int* IDs;
	int count = 0;
	int highestCPU = 0;
	int highestRAM = 0;


	void readVMP(string filename)
	{
		string line;
		ifstream VMPFile (filename);
		std::cout << filename << "\n";
		pmList = NULL;

		

		
		if (VMPFile.is_open())
		{
			int lineNum = 0;
			
			while (getline(VMPFile, line))
			{


				if (lineNum == 1)
				{
					std::stringstream numss(line);
					while (numss.good())
					{
						std::string substr;
						getline(numss,substr,',');
						subNums.push_back(stoi(substr));
					}
				}
				if (lineNum == 2)
				{
					std::stringstream cpuss(line);
					while (cpuss.good())
					{
						std::string substr;
						getline(cpuss,substr,',');
						PMCPU.push_back(stoi(substr));
					}

					std::reverse(PMCPU.begin(), PMCPU.end());
				

				}
				else if (lineNum == 3)
				{
					std::stringstream RAMss(line);
					while (RAMss.good())
					{
						std::string substr;
						getline(RAMss,substr,',');
						PMRAM.push_back(stoi(substr));
					}
					std::reverse(PMRAM.begin(), PMRAM.end());

				}

				if (lineNum == 4)
				{
					numVMs = stoi(line);
					numVMsPadded = numVMs;
					if ( numVMsPadded%_VECSIZE )
						numVMsPadded = _VECSIZE * (numVMsPadded/_VECSIZE + 1 );

					CPUReq = (float*)ALLOC(numVMs * sizeof( float ) );
					RAMReq = (float*)ALLOC(numVMs * sizeof( float ) );
					CPUReqTotal = 0;
					RAMReqTotal = 0;
					IDs = (int*)ALLOC(numVMs * sizeof(int));
					vmList = new VirtualMachine[numVMs];
				}

				else if (lineNum > 4)
				{
					vector<string> splitString;
					VirtualMachine vm;

					boost::split(splitString, line, boost::is_any_of(" "));
					vm.CPUReq = stoi(splitString[0]);
					vm.RAMReq = stoi(splitString[1]);
					if(vm.CPUReq > highestCPU)
					{
						highestCPU = vm.CPUReq;
					}

					if(vm.RAMReq > highestRAM)
					{
						highestRAM = vm.RAMReq;
					}
					CPUReq[lineNum-5] = stoi(splitString[0]);
					RAMReq[lineNum-5] = stoi(splitString[1]);
					CPUReqTotal += stoi(splitString[0]);
					RAMReqTotal += stoi(splitString[1]);
					vm.expiry = stoi(splitString[2]);
					vm.pmAssigned = -1;
					vm.id = lineNum-5;
					IDs[lineNum-5] = lineNum-5;
					vm.resourceBalance = (float)(vm.CPUReq -vm.RAMReq)/200;
					
					vmList[count] = vm;

					count++;

				}

				
				lineNum++;

			}

		}

		
		numPMs = numVMs;
		numPMPadded = numPMs;
		if (numPMPadded %_VECSIZE)
			numPMPadded = _VECSIZE * (numPMPadded /_VECSIZE +1 );
		pmList = new PhysicalMachine[numPMs];
		std::vector<PhysicalMachine> pmVec;
		//RAMCapacity = new int[numPMs];
		//CPUCapacity = new int[numPMs];
		RAMCapacity = (int*)ALLOC(numPMPadded*sizeof(int));
		CPUCapacity = (int*)ALLOC(numPMPadded*sizeof(int));


		

		for(int i = 0; i < numVMs; i++)
		{
			if(vmList[i].CPUReq <= highestCPU/20)
			{
				vmList[i].smallCPU = true;
			}

			if(vmList[i].RAMReq <= highestRAM/20)
			{
				vmList[i].smallRAM = true;
			}
		}
		
		int count = 0;

		if(subNums.size() > 1)
		{
			multipleServerTypes = true;
			aRAM = PMRAM.at(1);
			aCPU = PMCPU.at(1);
			bRAM = PMRAM.at(0);
			bCPU = PMCPU.at(0);
		}

		for(auto i : subNums)
		{
			for(int j = 0; j < i; j++)
			{
				PhysicalMachine pm;

			
				pm.RAMCapacity = PMRAM.back();
				pm.CPUCapacity = PMCPU.back();
				CPUCapAv += PMCPU.back();
				RAMCapAv += PMRAM.back();
			
				pm.ID = count;
				pmVec.push_back(pm);
				//pmList[count] = pm;
				RAMCapacity[count] = PMRAM.back();
				CPUCapacity[count] = PMCPU.back();
				count++;

			}

			PMRAM.pop_back();
			PMCPU.pop_back();
			
		}

		int n = sizeof(pmList) / sizeof(pmList[0]);
		std::sort (pmVec.begin(), pmVec.end(), compFunc);
		std::sort (CPUCapacity, CPUCapacity + numPMs, greater<int>());
		std::sort (RAMCapacity, RAMCapacity + numPMs, greater<int>());


		for(int i = 0; i < numPMs; i++)
		{
			pmList[i] = pmVec.at(i);
		}
		CPUCapAv = ceil(CPUCapAv/numPMs);
		RAMCapAv = ceil(RAMCapAv/numPMs);
		int optimum;

		if(multipleServerTypes == true)
		{
			optimum = (numPMs/10) + std::max(((CPUReqTotal - bCPU * (numPMs/10))/aCPU), (RAMReqTotal - bRAM * (numPMs/10))/aRAM);
		}

		else
		{
			optimum =  std::max(CPUReqTotal/CPUCapAv, RAMReqTotal/RAMCapAv);

		}
		std::cout << "Optimal PM num: " << optimum << "\n";

	}

	static bool compFunc(const PhysicalMachine pm1, const PhysicalMachine pm2)
	{
		return (pm1.CPUCapacity > pm2.CPUCapacity);
	}
};
#endif
