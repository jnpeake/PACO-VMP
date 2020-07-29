//#include "VMP.h"
#include "antsystem.h"
#include "ranluxgen.h"
#include <iostream>
#include <array>
#include <iterator>
#include <algorithm>

#ifndef EMULATE
#include <immintrin.h>
#endif


void AntSystem::Init ( int nAnts, VMP *vmp, int seed, int _alpha, int _beta, float _rho)
{
    m_pVMP = vmp;
    m_nAnts = nAnts;
    alpha = _alpha;
    beta = _beta;
    rho = _rho;
    globalLowestPMUsage = 100000;
    globalBestSolution = new SolutionObj[m_pVMP->numPMs];
    iterBestSolution = new SolutionObj[m_pVMP->numPMs];
    m_pher = (float**)ALLOC((m_pVMP->numVMs*sizeof(float*)));
    m_weights = (float**)ALLOC((m_pVMP->numVMs*sizeof(float*)));
    numServers = m_pVMP->numPMs;
    timers = new Timers(5);
    int paddedVMs = m_pVMP->numVMs;

	if ( paddedVMs%_VECSIZE )
		paddedVMs = _VECSIZE * (paddedVMs/_VECSIZE + 1 );

    for(int i = 0; i < m_pVMP->numVMs; i++)
    {
        m_pher[i] = (float*)ALLOC((paddedVMs*sizeof(float)));
        m_weights[i] = (float*)ALLOC((paddedVMs*sizeof(float)));
        memset(m_pher[i],0,paddedVMs*sizeof(float));
        memset(m_weights[i],0,paddedVMs*sizeof(float));
    }

    m_pAnts = (Ant*)ALLOC((m_nAnts * sizeof(Ant)));
    rlgen.init(seed);

    for(int i = 0; i < m_nAnts;i++)
    {
        int ALIGN(seeds[_VECSIZE]);
        for(int j = 0; j < _VECSIZE; j++)
        {
            seeds[j] = rlgen.irand(0xFFFF);
        }
        m_pAnts[i].Init(this,seeds);
    }
}

void AntSystem::Reset(void)
{
    pher0 = 1.0f/numServers;

    for(int i = 0; i < m_pVMP->numVMs; i++)
    {
        for(int j = 0; j < m_pVMP->numVMs; j++)
        {
            m_pher[i][j] = pher0;
            float usageTotal = m_pVMP->vmList[j].RAMReq + m_pVMP->vmList[j].CPUReq;
            m_weights[i][j] = 0;
        }
    }

    float a = exp( log(0.05) / (float)m_pVMP->numVMs );
	int n = m_pVMP->numVMs;
    mmasConst = (1.0f - a ) / ( ((float)n + 1.0f) * a * 0.5f);
}

void AntSystem::DoTours(void)
{
    iterLowestPMUsage = 1000000;
    iterHighestUtilisation = -1;
    SolutionObj *bestSolution;
    
#pragma omp parallel for
    for(int i = 0; i < m_nAnts; i++)
    {
        m_pAnts[i].ConstructTour();
    }

    for(int i = 0; i < m_nAnts; i++)
    {
        if(m_pAnts[i].usedPMs < iterLowestPMUsage)
        {
            iterLowestPMUsage = m_pAnts[i].usedPMs;
            bestSolution = m_pAnts[i].solution;
            iterHighestUtilisation = m_pAnts[i].solUtilisation;
            bestAnt = m_pAnts[i];      
        }

        else if (m_pAnts[i].usedPMs == iterLowestPMUsage && m_pAnts[i].solUtilisation > iterHighestUtilisation)
        {
            iterLowestPMUsage = m_pAnts[i].usedPMs;
            bestSolution = m_pAnts[i].solution;
            iterHighestUtilisation = m_pAnts[i].solUtilisation;
            bestAnt = m_pAnts[i];
        }
    }

    SaveSolution(bestSolution, iterBestSolution);

    if(bestAnt.feasibleSolution == true)
    {  
        if(iterLowestPMUsage < globalLowestPMUsage)
        {
            globalLowestPMUsage = iterLowestPMUsage;
            globalHighestUtilisation = iterHighestUtilisation;
            lastImprovedIter = currIteration;
            SaveSolution(bestSolution, globalBestSolution);
            std::cout<<"New Global Best \n";
            timers->StopTimer(3);
            lastImprovedTime = timers->GetTimer(3);

            std::cout << "Iteration " << currIteration << "| PMs Utilised: "<< globalLowestPMUsage <<  " | Total Time: " << timers -> GetTimer(3) << "\n";
            timers->StartTimer(3);    
        }

        else if (iterLowestPMUsage == globalLowestPMUsage && iterHighestUtilisation > globalHighestUtilisation)
        {
            globalLowestPMUsage = iterLowestPMUsage;
            globalHighestUtilisation = iterHighestUtilisation;
            lastImprovedIter = currIteration;
            SaveSolution(bestSolution, globalBestSolution);
            std::cout<<"New Global Best \n";
            timers->StopTimer(3);
            lastImprovedTime = timers->GetTimer(3);
            std::cout << "Iteration " << currIteration <<"| PMs Utilised: "<< globalLowestPMUsage <<  " | Total Time: " << timers -> GetTimer(3)<< "\n";
            timers->StartTimer(3);
        }
        numServers = bestAnt.usedPMs -1;
    }

    else
    {
        LocalSearch(bestSolution);
    }
}

void AntSystem::Deposit(void)
{

    Vector ones;
    ones.set1(1);

    for(int i =0; i < numServers; i++)
    {  
        if(iterBestSolution[i].vm_list.size())
        {
            for(int j = 0; j < iterBestSolution[i].vm_list.size()-1; j++)
            {
                for(int k = 0; k < iterBestSolution[i].vm_list.size()-1; k++)
                {
                    float currPher = m_pher[iterBestSolution[i].vm_list[j].listPos][iterBestSolution[i].vm_list[k].listPos];
                    float newPher = currPher * (1 + (1/(CalculateSolutionPowerUsage(iterBestSolution))/365));
                    m_pher[iterBestSolution[i].vm_list[j].listPos][iterBestSolution[i].vm_list[k].listPos] = newPher;
                    m_pher[iterBestSolution[i].vm_list[k].listPos][iterBestSolution[i].vm_list[j].listPos] = newPher;
                }
            }
        }
    }

    float evapFac = 1-rho;
    Vector vEvap;
	vEvap.set1(evapFac);

    float pherMax = 1.0f/(rho * globalLowestPMUsage );
	float pherMin = pherMax * mmasConst;

    for(int i = 0; i < m_pVMP -> numVMs; i++)
    {
        for (int j = 0; j < m_pVMP -> numVMs/_VECSIZE; j++)
        { 	
            Vector pher;
            Vector CPUReq;
            CPUReq.load(m_pVMP->CPUReq+ j*_VECSIZE);
            pher.load(m_pher[i] + j*_VECSIZE);
            pher = pher * vEvap;
            pher.vecMax(pherMax);
            pher.vecMin(pherMin);
            store(m_pher[i] + j*_VECSIZE, pher);
        }
    }
}

void AntSystem::Iterate(void)
{
    timers->StartTimer(0);
    DoTours();
    timers->StopTimer(0);
    timers->StartTimer(1);
    timers->StopTimer(1);
    timers->StartTimer(2);
    Deposit();
    timers->StopTimer(2);
}

void AntSystem::FFDTour(void)
{
    int numPMs = numServers;
    std::vector<int> unusedVMs(numPMs);
	std::iota (std::begin(unusedVMs), std::end(unusedVMs), 0);

    PhysicalMachine* PMList = m_pVMP->pmList;
    VirtualMachine* VMList = m_pVMP->vmList;

    int* cpuUsage = new int[numPMs];
    int* ramUsage = new int[numPMs];
    int* pmAllocation = new int[numPMs];

    for(int i = 0; i < numServers; i++)
    {
        memset(cpuUsage,0,numServers*sizeof(int));
        memset(ramUsage,0,numServers*sizeof(int));
    }

    SolutionObj *solution = new SolutionObj[numServers];
    
    for (int i = 0; i < m_pVMP->numVMs;i++)
    {
        for(int j = 0; j < numServers; j++)
        {
            if(VMList[i].CPUReq < (PMList[j].CPUCapacity-cpuUsage[j]) && VMList[i].RAMReq < (PMList[j].RAMCapacity -ramUsage[j]))
            {
                pmAllocation[VMList[i].id] = j;
                cpuUsage[j] += VMList[i].CPUReq;
                ramUsage[j] += VMList[i].RAMReq;
                break;
            }
        }        
    }

    for(int i = 0; i < numServers; i++)
    {
        solution[i].pm_id = PMList[i].ID;
        solution[i].cpu_usage = cpuUsage[i];
        solution[i].ram_usage = ramUsage[i];
        unusedVMs.erase(std::remove(unusedVMs.begin(), unusedVMs.end(), i), unusedVMs.end());
        for(int j = 0; j < m_pVMP->numVMs; j++)
        {
            if(pmAllocation[j] == i)
            {
                solution[i].vm_list.push_back(VMList[j]);
            }
        }
    }

    int usedServers = 0;

    for(int i = 0; i < numServers; i++)
    {
        if(solution[i].vm_list.size() > 0)
        {
            usedServers++;
        }
    }
    
    float solUtilisation;
    for(int i = 0; i < numPMs; i++)
    {
        if(solution[i].vm_list.size() > 0 && solution[i].vm_list.size() < 10000)
            solUtilisation += ((m_pVMP->pmList[i].CPUCapacity - solution[i].cpu_usage / m_pVMP->pmList[i].CPUCapacity)/(m_pVMP->pmList[i].RAMCapacity - solution[i].ram_usage / m_pVMP->pmList[i].RAMCapacity));
    }

    SaveSolution(solution,globalBestSolution);
    globalLowestPMUsage = usedServers;
    globalHighestUtilisation = solUtilisation;
    delete[] solution;
    delete[] cpuUsage;
    delete[] ramUsage;
    delete[] pmAllocation;

    std::cout << "\ngreedy finished \n";
    std::cout << usedServers << "\n";
    for ( auto i = unusedVMs.begin(); i != unusedVMs.end(); i++ ) {
    	std::cout << *i << std::endl;
	}    
}

void AntSystem::Solve(int maxIterations)
{   
    Reset();
    timers->StartTimer(3);
    timers->StartTimer(4);
    for(int i = 0; i < maxIterations; i++)
    {
        int randInt = rlgen.irand(0xFFFF);
        std::srand(randInt);
        std::random_shuffle(m_pVMP->vmList, m_pVMP->vmList + m_pVMP->numVMs);

        for(int j = 0; j < m_pVMP->numVMs; j++)
        {
            m_pVMP->vmList[j].listPos = j;
        }

        std::srand(randInt);
        std::random_shuffle(m_pVMP->CPUReq, m_pVMP->CPUReq + m_pVMP->numVMs);

        std::srand(randInt);
        std::random_shuffle(m_pVMP->RAMReq, m_pVMP->RAMReq + m_pVMP->numVMs);

        std::srand(randInt);
        std::random_shuffle(m_pVMP->IDs, m_pVMP->IDs + m_pVMP->numVMs);

        for(int j = 0; j < m_pVMP->numVMs; j++)
        {
            std::srand(randInt);
            std::random_shuffle(m_pher[j], m_pher[j] + m_pVMP->numVMs);
        }

        std::srand(randInt);
        std::random_shuffle(m_pher, m_pher + m_pVMP->numVMs);

        Iterate();
        timers-> ClearTimer(0);
        timers-> ClearTimer(1);
        timers-> ClearTimer(2);
    }

    timers->StopTimer(3);
    std::cout << "\nCOMPLETE :"<<"| PMs Utilised: "<< globalLowestPMUsage <<" | Num VMs: " << m_pVMP->numVMs <<" | Last Improved Iter: "<< lastImprovedIter <<" | Last Improved Time: " <<lastImprovedTime <<  " | Total Time: " << timers -> GetTimer(3) <<"\n";

    for(int i = 0; i < m_pVMP->numPMs; i++)
    {
        if(globalBestSolution[i].vm_list.size() > 0)
        {
            std::cout << "\n PM: " << globalBestSolution[i].pm_id << ", CPU: " << globalBestSolution[i].cpu_usage << ", RAM: " << globalBestSolution[i].ram_usage << ", VMs:";

                for(int j = 0; j < globalBestSolution[i].vm_list.size(); j++)
                {
                    std::cout << globalBestSolution[i].vm_list[j].id << ", ";
                }
        }

    }
    std::cout << std::flush;
}

float AntSystem::CalculateSolutionPowerUsage(SolutionObj *solution)
{
    float usageTotal = 0;
    for(int i = 0; i < numServers; i++)
    {
        if(solution[i].cpu_usage > 0)
        {
            usageTotal+= (218-171) * (solution[i].cpu_usage/m_pVMP->pmList[i].CPUCapacity) +171;
        }
    }
    return usageTotal;
}

float AntSystem::CalculateSolutionPowerUsage(std::list<SolutionObj> solution)
{
    float usageTotal = 0;
    for(SolutionObj currSolution : solution)
    {
        if(currSolution.cpu_usage > 0)
        {
            usageTotal+= (218-171) * (currSolution.cpu_usage/m_pVMP->pmList[currSolution.pm_id].CPUCapacity) +171;
        }
    }
    return usageTotal;
}


float AntSystem::CalculateSolutionResourceWastage(SolutionObj *solution)
{
    float wastageTotal = 0;
    for(int i = 0; i < numServers; i++)
    {
        float cpuUsageRatio = solution[i].cpu_usage/m_pVMP->pmList[solution[i].pm_id].CPUCapacity;
        float RAMUsageRatio = solution[i].ram_usage/m_pVMP->pmList[solution[i].pm_id].RAMCapacity;

        wastageTotal = abs(cpuUsageRatio - RAMUsageRatio);

    }
    return wastageTotal;
}

float AntSystem::CalculateSolutionResourceWastage(std::list<SolutionObj> solution)
{
    float wastageTotal = 0;
    for(SolutionObj currSolution: solution)
    {
        float cpuUsageRatio = currSolution.cpu_usage/m_pVMP->pmList[currSolution.pm_id].CPUCapacity;
        float RAMUsageRatio = currSolution.ram_usage/m_pVMP->pmList[currSolution.pm_id].RAMCapacity;

        if(cpuUsageRatio <= RAMUsageRatio)
        {
            wastageTotal += RAMUsageRatio - cpuUsageRatio;
        }

        else
        {
            wastageTotal += cpuUsageRatio - RAMUsageRatio;
        }
    }
    return wastageTotal;
}

void AntSystem::LocalSearch(SolutionObj* solution)
{
   std::vector<int> resolvedPMs;
   for(int pm : bestAnt.overloadedPMs)
   {
        std::sort(bestAnt.solution[pm].vm_list.begin(), bestAnt.solution[pm].vm_list.end(),[](VirtualMachine const & a, VirtualMachine const & b) -> bool{ return a.resourceBalance > b.resourceBalance; } );
        bool swapFound = false;
        for(int i = 0 ; i < numServers; i++)
        {
            if(swapFound == true)
            {
                break;
            }
            if(bestAnt.solution[i].vm_list.size() > 0 && i != pm )
            {
                std::sort(bestAnt.solution[i].vm_list.begin(), bestAnt.solution[i].vm_list.end(),[](VirtualMachine const & a, VirtualMachine const & b) -> bool{ return a.resourceBalance < b.resourceBalance; } );

                for(int j = 0; j < bestAnt.solution[pm].vm_list.size(); j++)
                {
                    if(swapFound == true)
                    {
                        break;
                    }
                    for(int k = 0; k < bestAnt.solution[i].vm_list.size(); k++)
                    {
                        float tempOverloadedCPU = bestAnt.solution[pm].cpu_usage - bestAnt.solution[pm].vm_list[j].CPUReq + bestAnt.solution[i].vm_list[k].CPUReq;
                        float tempOverloadedRAM = bestAnt.solution[pm].ram_usage - bestAnt.solution[pm].vm_list[j].RAMReq + bestAnt.solution[i].vm_list[k].RAMReq;
                        float tempSwappedCPU = bestAnt.solution[i].cpu_usage - bestAnt.solution[i].vm_list[k].CPUReq + bestAnt.solution[pm].vm_list[j].CPUReq;
                        float tempSwappedRAM = bestAnt.solution[i].ram_usage - bestAnt.solution[i].vm_list[k].RAMReq + bestAnt.solution[pm].vm_list[j].RAMReq;

                        if(tempOverloadedCPU <= m_pVMP->pmList[pm].CPUCapacity && tempOverloadedRAM <= m_pVMP->pmList[pm].RAMCapacity)
                        {
                            if(tempSwappedCPU <= m_pVMP->pmList[i].CPUCapacity && tempSwappedRAM <= m_pVMP->pmList[i].RAMCapacity)
                            {
                                swapFound = true;
                                bestAnt.solution[pm].vm_list.push_back( bestAnt.solution[i].vm_list[k]);
                                bestAnt.solution[i].vm_list.push_back( bestAnt.solution[pm].vm_list[j]);

                                bestAnt.solution[pm].ram_usage = tempOverloadedRAM;
                                bestAnt.solution[i].ram_usage = tempSwappedRAM;

                                bestAnt.solution[pm].cpu_usage = tempOverloadedCPU;
                                bestAnt.solution[i].cpu_usage = tempSwappedCPU;

                                bestAnt.solution[pm].cpu_usage_ratio = bestAnt.solution[pm].cpu_usage / m_pVMP->pmList[pm].CPUCapacity;
                                bestAnt.solution[i].cpu_usage_ratio = bestAnt.solution[i].cpu_usage / m_pVMP->pmList[i].CPUCapacity;

                                bestAnt.solution[pm].ram_usage_ratio = bestAnt.solution[pm].ram_usage / m_pVMP->pmList[pm].RAMCapacity;
                                bestAnt.solution[i].ram_usage_ratio = bestAnt.solution[i].ram_usage / m_pVMP->pmList[i].RAMCapacity;

                                bestAnt.solution[pm].vm_list.erase(bestAnt.solution[pm].vm_list.begin()+j);
                                bestAnt.solution[i].vm_list.erase(bestAnt.solution[i].vm_list.begin()+k);

                                bestAnt.solution[pm].resourceBalance = bestAnt.solution[pm].cpu_usage_ratio - bestAnt.solution[pm].ram_usage_ratio;
                                bestAnt.solution[i].resourceBalance = bestAnt.solution[i].cpu_usage_ratio - bestAnt.solution[i].ram_usage_ratio;

                                resolvedPMs.push_back(pm);

                                break;
                            }
                        }
                    }
                }
            }
        }
   }

   for(int pm : bestAnt.overloadedPMs)
   {
       if(std::find(resolvedPMs.begin(), resolvedPMs.end(), pm) == resolvedPMs.end()) 
       {      
            for(int i = 0; i < bestAnt.solution[pm].vm_list.size(); i++)
            {
                for(int j = 0 ; j < numServers; j++)
                {
                    float tempOverloadedCPU = bestAnt.solution[pm].cpu_usage - bestAnt.solution[pm].vm_list[i].CPUReq;
                    float tempOverloadedRAM = bestAnt.solution[pm].ram_usage - bestAnt.solution[pm].vm_list[i].RAMReq;

                    float tempMigratedCPU = bestAnt.solution[j].cpu_usage + bestAnt.solution[pm].vm_list[i].CPUReq;
                    float tempMigratedRAM = bestAnt.solution[j].ram_usage + bestAnt.solution[pm].vm_list[i].RAMReq;

                    if(tempOverloadedCPU <= m_pVMP->pmList[pm].CPUCapacity && tempOverloadedRAM <= m_pVMP->pmList[pm].RAMCapacity)
                    {
                        if(tempMigratedCPU <= m_pVMP->pmList[j].CPUCapacity && tempMigratedRAM <= m_pVMP->pmList[j].RAMCapacity)
                        {
                            bestAnt.solution[j].vm_list.push_back( bestAnt.solution[pm].vm_list[i]);

                            bestAnt.solution[pm].ram_usage = tempOverloadedRAM;
                            bestAnt.solution[j].ram_usage = tempMigratedRAM;

                            bestAnt.solution[pm].cpu_usage = tempOverloadedCPU;
                            bestAnt.solution[j].cpu_usage = tempMigratedCPU;

                            bestAnt.solution[pm].cpu_usage_ratio = bestAnt.solution[pm].cpu_usage / m_pVMP->pmList[pm].CPUCapacity;
                            bestAnt.solution[j].cpu_usage_ratio = bestAnt.solution[j].cpu_usage / m_pVMP->pmList[j].CPUCapacity;

                            bestAnt.solution[pm].ram_usage_ratio = bestAnt.solution[pm].ram_usage / m_pVMP->pmList[pm].RAMCapacity;
                            bestAnt.solution[j].ram_usage_ratio = bestAnt.solution[j].ram_usage / m_pVMP->pmList[j].RAMCapacity;

                            bestAnt.solution[pm].vm_list.erase(bestAnt.solution[pm].vm_list.begin()+i);

                            bestAnt.solution[pm].resourceBalance = bestAnt.solution[pm].cpu_usage_ratio - bestAnt.solution[pm].ram_usage_ratio;
                            bestAnt.solution[j].resourceBalance = bestAnt.solution[j].cpu_usage_ratio - bestAnt.solution[j].ram_usage_ratio;

                            resolvedPMs.push_back(pm);

                            break;
                        }
                    }
                }
            }
        } 
    }
    bestAnt.solUtilisation = 0;
    for(int i = 0; i < numServers+1; i++)
    {
        if(bestAnt.solution[i].vm_list.size() > 0 && bestAnt.solution[i].vm_list.size() < 10000)
        {
            bestAnt.solUtilisation += ((m_pVMP->pmList[i].CPUCapacity - bestAnt.solution[i].cpu_usage / m_pVMP->pmList[i].CPUCapacity)/(m_pVMP->pmList[i].RAMCapacity - bestAnt.solution[i].ram_usage / m_pVMP->pmList[i].RAMCapacity));
        }
    }

    if(bestAnt.usedPMs < iterLowestPMUsage)
    {
        iterLowestPMUsage =bestAnt.usedPMs;
        iterHighestUtilisation = bestAnt.solUtilisation;
        SaveSolution(bestAnt.solution, iterBestSolution);
    }

    else if (bestAnt.usedPMs == iterLowestPMUsage && bestAnt.solUtilisation > iterHighestUtilisation)
    {
        iterLowestPMUsage = bestAnt.usedPMs;
        iterHighestUtilisation = bestAnt.solUtilisation;
        SaveSolution(bestAnt.solution, iterBestSolution);
    }

    if(bestAnt.overloadedPMs.size() == resolvedPMs.size())
    {
        int tempUsedPMs = 0;

        for(int i = 0; i < numServers; i++)
        {
            if(bestAnt.solution[i].vm_list.size() > 0)
            {
                tempUsedPMs++;
            }
        }

        numServers = tempUsedPMs-1;
        bestAnt.usedPMs = tempUsedPMs;
        std::cout << std::flush << "\n";
        iterLowestPMUsage = bestAnt.usedPMs;
        lastImprovedIter = currIteration;
        SaveSolution(bestAnt.solution, iterBestSolution);
        std::cout<<"Successful LS\n";
        timers->StopTimer(3);
        lastImprovedTime = timers->GetTimer(3);
        std::cout << "Iteration " << currIteration <<"| PMs Utilised: "<< globalLowestPMUsage <<" | Num VMs: " << m_pVMP->numVMs <<  " | Total Time: " << timers -> GetTimer(3);
        timers->StartTimer(3);

        if(iterLowestPMUsage<globalLowestPMUsage)
        {
            SaveSolution(bestAnt.solution, globalBestSolution);
            globalHighestUtilisation = iterHighestUtilisation =  bestAnt.solUtilisation;
            globalLowestPMUsage = bestAnt.usedPMs;
        }
   }   
}

void AntSystem::SaveSolution(SolutionObj* solution, SolutionObj* bestSolution)
{
    for(int i = 0; i < numServers; i++)
    {
        bestSolution[i].cpu_usage = solution[i].cpu_usage;
        bestSolution[i].ram_usage = solution[i].ram_usage;
        bestSolution[i].pm_id = solution[i].pm_id;
        bestSolution[i].cpu_usage_ratio = solution[i].cpu_usage_ratio;
        bestSolution[i].ram_usage_ratio = solution[i].ram_usage_ratio;
        bestSolution[i].resourceBalance = solution[i].resourceBalance;
        bestSolution[i].vm_list = solution[i].vm_list;
    }
}


