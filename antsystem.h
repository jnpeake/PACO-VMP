#ifndef _ANTSYSTEM_INC_
#define _ANTSYSTEM_INC_

#include "VMP.h"
#include "ant.h"
#include "ranluxgen.h"
#include "timers.h"
#include "platform.h"
#include <vector>


struct SolutionObj{
    int pm_id;
    float cpu_usage;
    float ram_usage;
    float cpu_usage_ratio;
    float ram_usage_ratio;
    float resourceBalance;
    std::vector<VirtualMachine> vm_list;
};


class Ant;

class AntSystem
{
    friend class Ant;

public:
   
   
    int m_nAnts; //number of ants
    int iterLowestPMUsage;
    int globalLowestPMUsage;
    int iterationBestAnt;
    int iterBestRWAnt;
    int iterBestPUAnt;
    int numServers;
    int currIteration;
    int alpha;
    int beta;
    int lastImprovedIter;
    float lastImprovedTime;
    float epsilon;
    float rho;
    float pher0;
    float mmasConst;
    float iterHighestUtilisation;
    float globalHighestUtilisation;
    float evapFac;
    float bestPU;
    float (**m_pher);
    float (**m_weights);

    Ant bestAnt;
    Ant globalBestAnt;
    SolutionObj* globalBestSolution;
    SolutionObj* iterBestSolution;
    Timers* timers;
    RanluxGenerator rlgen;
    VMP *GetVMP( void ) { return m_pVMP; }
    Ant *m_pAnts; //list of all ants
    VMP *m_pVMP; //Virtual Machine Placement Problem

    void Init (int nAnts, VMP *vmp, int seed, int alpha, int beta, float rho);
    void Reset(void);
    void DoTours(void);
    void Deposit(void);
    void Iterate(void);
    void FFDTour(void);
    float CalculatePowerUsage(PhysicalMachine pm);
    float CalculateSolutionPowerUsage(SolutionObj *solution);
    float CalculateSolutionPowerUsage(std::list<SolutionObj> solution);
    float CalculateSolutionResourceWastage(SolutionObj *solution);
    float CalculateSolutionResourceWastage(std::list<SolutionObj> solution);
    void Solve(int maxIterations);
    void LocalSearch(SolutionObj* solution);
    void SaveSolution(SolutionObj *solution, SolutionObj* bestSolution);
};
#endif