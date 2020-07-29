#include "VMP.h"
#include "antsystem.h"
#include <cstdlib>
#include <iostream>

int main( int argc, char *argv[] )
{
    Timers* timers = new Timers(1);

    VMP *vmp = new VMP();
    AntSystem *as = new AntSystem();
    

    if(argc < 7)
    {
        std::perror("ERROR: Not enough arguments (VMP, numAnts, numIterations, seed, alpha, beta, rho)");
        
    }    

    else
    {
        std::string filename(argv[1]);
        vmp->readVMP("Instances/" + filename);
        int numAnts;
        sscanf( argv[2], "%d", &numAnts);
        int numIterations;
        sscanf( argv[3], "%d", &numIterations);
        int seed;
        sscanf( argv[4], "%d", &seed);
        int alpha;
        sscanf( argv[5], "%d", &alpha);
        int beta;
        sscanf( argv[6], "%d", &beta);
        float rho;
        sscanf( argv[7], "%f", &rho);

        timers->StartTimer(0);
        as->Init(numAnts,vmp,seed, alpha, beta, rho);
        as->FFDTour();
        as->Solve(numIterations);
        timers->StopTimer(0);
        free(vmp);
       
        std::cout << "\nRuntime: " << timers->GetTimer(0) << " PMs: " << as->globalLowestPMUsage;
        free(as);
        free(timers);
    }
    
    
    
    
}