#pragma once
#ifndef _ANT_INC_
#define _ANT_INC_
#include "vector.h"
#include "platform.h"
#include "VMP.h"
#include "ranluxgen.h"

#include <random>

class Vector;

class AntSystem;

struct SolutionObj;

class Ant
{
public:
    int numVMPadded;
    int usedPMs;
    int usedVMs;
    float overload;
    std::vector<int> availableVMs;
    std::vector<int> overloadedPMs;
    float pherList[8];
    float solUtilisation;
    float** pherAvg;
    bool firstIt;

    
    AntSystem *m_as;
    //unsigned int rSeed[8];
    VMP *vmp;
    SolutionObj* solution;
    Vector VecCalculateResourceWastage(PhysicalMachine pm, Vector vmRAMReq, Vector vmCPUReq);
    Vector VecCalcHeuristic(Vector CPUUsage, Vector RAMUsage, Vector CPUCap, Vector RAMCap, int vmCPUReq, int vmRAMReq);
    float CalculateSolutionResourceWastage(SolutionObj *solution);
    RanluxGenerator rlgen;
    std::vector<pair<int,int>> assignmentList;



    Vector rSeed;
	Vector rC0;
	Vector rC1;
	Vector factor;
	Vector ones; // constant
    float *tabu; 
    int* pmAssigned;
    bool dominated;
    bool feasibleSolution;

    void VecRoulette(float **weights, float **m_pher_RW, float *tabu, int numPMs, VMP *m_pVMP);
    void ConstructTour();
    void Init(AntSystem *as, int *seeds);
    void VecCheck(int chr);




};
#endif