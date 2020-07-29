#include "ant.h"
#include "antsystem.h"
#include <cstdio> 
#include <cstdlib>
#include <iostream>
#include <algorithm>



void Ant::Init(AntSystem *as, int *seeds)
{
    m_as = as;
    std::random_device rd;

    rlgen.init(seeds[0]);

    seedVecRandom(rC0,rC1,factor,seeds,rSeed);

    numVMPadded = as->m_pVMP->numVMs;
    if (numVMPadded %_VECSIZE)
		numVMPadded = _VECSIZE * (numVMPadded /_VECSIZE +1 );

    tabu = (float*)ALLOC(numVMPadded*sizeof(float));

    
    ones.set1(1.0f);
    firstIt = false;

}

void Ant::ConstructTour()
{
    overloadedPMs.clear();
    availableVMs.clear();
    solUtilisation = 0;
    dominated = false;
    feasibleSolution = false;
    usedPMs = 0;
    usedVMs = 0;
    VMP* m_pVMP = m_as->GetVMP();
    int numPMs = m_as->numServers;

    memset(tabu, 0, numVMPadded*sizeof(float));

    Vector randoms = vecRandom(rC0, rC1, factor, rSeed);
    ALIGN(float r[_VECSIZE]);
    store(r,randoms);

    if(firstIt == true)
    {
        delete[] solution;
        
    }

    firstIt = false;
    solution = new SolutionObj[numPMs];

    
    std::fill_n(tabu,m_pVMP->numVMs, 1);
    int numPMsPadded = numPMs;
    if ( numPMs%_VECSIZE )
                numPMsPadded = _VECSIZE * (numPMs/_VECSIZE + 1 );

    for(int i = numPMs; i < numPMsPadded; i++)
    {
        tabu[i] = -1;
    }

    VecRoulette(m_as->m_weights, m_as->m_pher, tabu, numPMs, m_pVMP);    
}

void Ant::VecRoulette(float **weights, float **m_pher, float *tabu, int numPMs, VMP *m_pVMP)
{

    ALIGN(float indexSeed[_VECSIZE]) = { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f};
	ALIGN(float indexStep[_VECSIZE]) = { 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f};
	ALIGN(float minusOnes[_VECSIZE]) = { -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f};
    ALIGN(float pherList[_VECSIZE]);
    int numPMsPadded = numPMs;
    if ( numPMs%_VECSIZE )
                numPMsPadded = _VECSIZE * (numPMs/_VECSIZE + 1 );
                
    pherAvg = (float**) ALLOC ((numVMPadded * sizeof(float*)));

    for(int i = 0; i < numVMPadded; i++)
    {
        pherAvg[i] = (float*) ALLOC ((numPMsPadded * sizeof(float)));
    }

    ALIGN(float PMCUsage[numPMsPadded]);
    ALIGN(float PMRUsage[numPMsPadded]);

    ALIGN(float PMCCap[numPMsPadded]);
    ALIGN(float PMRCap[numPMsPadded]);


    float* allocatedVMs;
    allocatedVMs = (float*)ALLOC((numVMPadded*sizeof(float)));
    

    memset(PMCUsage,0,numPMsPadded*sizeof(float));
    memset(PMRUsage,0,numPMsPadded*sizeof(float));
    memset(PMCCap,0,numPMsPadded*sizeof(float));
    memset(PMRCap,0,numPMsPadded*sizeof(float));
    memset(allocatedVMs,0,numVMPadded*sizeof(float));
    for(int i = 0; i < numPMs; i++)
    {
        PMCCap[i] = m_pVMP->CPUCapacity[i];
        PMRCap[i] = m_pVMP->RAMCapacity[i];
    }

    for(int i = 0; i < m_pVMP->numVMs; i++)
    {
        std::fill(pherAvg[i], pherAvg[i] + numPMs, m_as->pher0);   
        for(int j = numPMs; j < numPMsPadded; j++)
        {
            pherAvg[i][j] = -10000000;
        }
    }

    for(int i = m_pVMP->numVMs; i < numVMPadded; i++)
    {
        for(int j = 0; j < numPMsPadded; j++)
        {
            pherAvg[i][j]=-10000000;
        }
    }
    Vector minusOne;
    Vector delta16;
    minusOne.load(minusOnes);
    delta16.load (indexStep);

    std::list<int> currIDs;
    std::vector<int> usedServers;
    Timers* timers;

    for(int i = 0; i < m_pVMP->numVMs; i++)
    {
        Vector runningIndex;
        Vector currentPMs = minusOne;
        Vector currentWeights = minusOne;
        runningIndex.load(indexSeed);

        for(int j = 0; j < numPMsPadded/_VECSIZE; j++)
        {
            Vector pher;
            Vector weights;
            Vector CPUCapacity;
            Vector RAMCapacity;
            Vector CPUUsage;
            Vector RAMUsage;
            Vector CPUReq;
            Vector RAMReq;
            Vector combinedWeight;
            Vector nextPMs = runningIndex;
            Vector tabuMask;
            tabuMask.load(tabu + j * _VECSIZE);

            float highestVal = 0.0f;
            int highestIndex = -1; 
            Vector randoms = vecRandom(rC0, rC1, factor, rSeed);
            CPUCapacity.load(PMCCap + j*_VECSIZE);
            RAMCapacity.load(PMRCap + j*_VECSIZE);
            CPUUsage.load(PMCUsage + j*_VECSIZE);
            RAMUsage.load(PMRUsage + j*_VECSIZE);
            Vector CPUMax = CPUCapacity - CPUUsage;
            Vector RAMMax = RAMCapacity - RAMUsage;
            CPUReq.set1(m_pVMP->vmList[i].CPUReq);
            RAMReq.set1(m_pVMP->vmList[i].RAMReq);
            Vector maxRAMMask = ltMask(RAMReq, RAMMax);
            Vector maxCPUMask = ltMask(CPUReq, CPUMax);

            pher.load(pherAvg[i] + j*_VECSIZE);
            weights = VecCalcHeuristic(CPUUsage,RAMUsage,CPUCapacity,RAMCapacity, m_pVMP->vmList[i].CPUReq, m_pVMP->vmList[i].RAMReq);
            Vector origPher = pher;
            Vector origWeight = weights;
            for(int k = 1; k < m_as->alpha; k++)
            {
                pher = pher*origPher;
            }

            for(int k = 1; k < m_as->beta; k++)
            {
                weights = weights*origWeight;
            }
            combinedWeight = pher * weights * randoms;
            combinedWeight = mask_mov(combinedWeight, maxCPUMask, minusOne);
            combinedWeight = mask_mov(combinedWeight, maxRAMMask, minusOne);
            combinedWeight = mask_mov(minusOne, tabuMask, combinedWeight);

            maxLocStep(currentWeights, currentPMs, combinedWeight, nextPMs);
            runningIndex = runningIndex + delta16;

            for(int k = 0; k < _VECSIZE; k++)
            {
                if(currentWeights.AVXVec[k] > 0 || currentWeights.AVXVec[k] == -1)
                {
                    break;
                }
            }
        }

        int reduced = reduceMax(currentWeights, currentPMs);

        if (reduced < 0) {
            reduced = -1;
        }

        if (reduced == -1)
        {
            availableVMs.push_back(i);
        }

        else
        {
            usedVMs++;
            PMCUsage[reduced] += m_pVMP->CPUReq[i];
            PMRUsage[reduced] += m_pVMP->RAMReq[i];
            solution[reduced].pm_id = m_pVMP->pmList[reduced].ID;
            solution[reduced].cpu_usage = PMCUsage[reduced];
            solution[reduced].ram_usage = PMRUsage[reduced];
            solution[reduced].cpu_usage_ratio = (float)PMCUsage[reduced]/m_pVMP->pmList[reduced].CPUCapacity;
            solution[reduced].ram_usage_ratio = (float)PMRUsage[reduced]/m_pVMP->pmList[reduced].RAMCapacity;
            solution[reduced].resourceBalance = solution[reduced].cpu_usage_ratio - solution[reduced].ram_usage_ratio;
            solution[reduced].vm_list.push_back(m_pVMP->vmList[i]);
            usedServers.push_back(reduced);
            allocatedVMs[reduced]++;
    
            if(allocatedVMs[reduced] == 1)
            {
                usedPMs++;
                for(int j = 0; j < m_pVMP->numVMs; j++)
                {
                    pherAvg[j][reduced] = m_pher[i][j];
                }
            }

            else
            {
                for(int j = 0; j < numVMPadded/_VECSIZE; j++)
                {
                    Vector currPher;
                    Vector newPher;
                    Vector assignedVMsVec;
                    ALIGN(float currPherArr[_VECSIZE]);
                    ALIGN(float newPherArr[_VECSIZE]);
                    ALIGN(float avPherArr[_VECSIZE]);

                    for(int k = 0; k < _VECSIZE; k++)
                    {
                        if((j*_VECSIZE) + k < m_pVMP->numVMs)
                        {
                            currPherArr[k] = pherAvg[(j*_VECSIZE) + k][reduced];
                            newPherArr[k] = m_pher[(j*_VECSIZE) + k][i];
                        }

                        else
                        {
                            
                            currPherArr[k] = -1;
                            newPherArr[k] = -1;
                        }
                        
                    }
                    assignedVMsVec.set1(allocatedVMs[reduced]);
                    currPher.load(currPherArr);
                    newPher.load(newPherArr);

                    currPher = currPher+newPher;
                    currPher = currPher/assignedVMsVec;

                    store(avPherArr, currPher);

                    for(int k = 0; k < _VECSIZE; k++)
                    {
                        pherAvg[(j*_VECSIZE)+k][reduced] = avPherArr[k];
                    }

                }
            }
            
        }
    }

    if(usedVMs == m_pVMP->numVMs)
        {
            feasibleSolution = true;
        }

        else
        {
            for(int vm : availableVMs)
            {
                float lowestOverload = 1e20f;
                int lowestOverloadPM = -1;

                for(int i = 0; i < m_as->numServers; i++)
                {
                    float cpuOverload = abs((m_pVMP->pmList[i].CPUCapacity - solution[i].cpu_usage -  m_pVMP->vmList[vm].CPUReq)/m_pVMP->pmList[i].CPUCapacity);
                    float ramOverload = abs((m_pVMP->pmList[i].RAMCapacity - solution[i].ram_usage -  m_pVMP->vmList[vm].RAMReq)/m_pVMP->pmList[i].RAMCapacity);
                    float totalOverload = cpuOverload + ramOverload;
                    if( totalOverload < lowestOverload)
                    {
                        lowestOverload = totalOverload;
                        lowestOverloadPM = i;
                    }

                }

                overloadedPMs.push_back(lowestOverloadPM);
                overload+= lowestOverload;
                solution[lowestOverloadPM].vm_list.push_back(m_pVMP->vmList[vm]);

                PMCUsage[lowestOverloadPM] += m_pVMP->CPUReq[vm];
                PMRUsage[lowestOverloadPM] += m_pVMP->RAMReq[vm];
                solution[lowestOverloadPM].cpu_usage = solution[lowestOverloadPM].cpu_usage + m_pVMP->vmList[vm].CPUReq;
                solution[lowestOverloadPM].ram_usage = solution[lowestOverloadPM].ram_usage + m_pVMP->vmList[vm].RAMReq;
                solution[lowestOverloadPM].cpu_usage_ratio = (float)solution[lowestOverloadPM].cpu_usage/m_pVMP->pmList[lowestOverloadPM].CPUCapacity;
                solution[lowestOverloadPM].ram_usage_ratio = (float)solution[lowestOverloadPM].ram_usage/m_pVMP->pmList[lowestOverloadPM].RAMCapacity;
                solution[lowestOverloadPM].resourceBalance = solution[lowestOverloadPM].cpu_usage_ratio - solution[lowestOverloadPM].ram_usage_ratio;


            }
        }
  
        FREE (allocatedVMs);
        for(int i = 0; i < numVMPadded; i++)
        {
           FREE(pherAvg[i]);
        }

        FREE(pherAvg);

       
    for(int i = 0; i < numPMs; i++)
    {
        if(solution[i].vm_list.size() > 0 && solution[i].vm_list.size() < 10000)
            solUtilisation += ((m_pVMP->pmList[i].CPUCapacity - solution[i].cpu_usage / m_pVMP->pmList[i].CPUCapacity)/(m_pVMP->pmList[i].RAMCapacity - solution[i].ram_usage / m_pVMP->pmList[i].RAMCapacity));
    }
    

}

Vector Ant::VecCalcHeuristic(Vector CPUUsage, Vector RAMUsage, Vector CPUCap, Vector RAMCap, int vmCPUReq, int vmRAMReq)
{
    Vector ones;
    ones.set1(1);

    Vector twos;
    twos.set1(2);

    Vector CPUUsageRatio, RAMUsageRatio;
    
    Vector vmCPU, vmRAM;
    vmCPU.set1(vmCPUReq);
    vmRAM.set1(vmRAMReq);
   
    CPUUsageRatio = (CPUCap - CPUUsage - vmCPU)/CPUCap;
    RAMUsageRatio = (RAMCap - RAMUsage - vmRAM)/RAMCap;

    Vector numerator, denominator;
    numerator = ones - abs(CPUUsageRatio - RAMUsageRatio);
    denominator = abs(CPUUsageRatio) + abs(RAMUsageRatio) + ones;

    return numerator/denominator;

}

float Ant::CalculateSolutionResourceWastage(SolutionObj *solution)
{

    float wastageTotal = 0;
    for(int i = 0; i < m_as->numServers; i++)
    {
        float cpuUsageRatio = solution[i].cpu_usage/m_as->GetVMP()->pmList[i].CPUCapacity;
        float RAMUsageRatio = solution[i].ram_usage/m_as->GetVMP()->pmList[i].RAMCapacity;

        wastageTotal = abs (cpuUsageRatio - RAMUsageRatio);

    }
    return wastageTotal;
}


