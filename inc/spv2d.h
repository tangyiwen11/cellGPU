//spv.h
#ifndef SPV_H
#define SPV_H


#include "std_include.h"
#include <stdio.h>
#include <cmath>
#include <random>
#include <sys/time.h>
#include "cuda_runtime.h"
#include "curand.h"
#include "curand_kernel.h"
#include "vector_types.h"
#include "vector_functions.h"

using namespace std;

#include "Matrix.h"
#include "cu_functions.h"

#include "DelaunayMD.h"

class SPV2D : public DelaunayMD
    {
    protected:
        Dscalar Dr;
        Dscalar v0;

        Dscalar gamma; //value of inter-cell surface tension
        bool useTension;
        bool particleExclusions;

        GPUArray<Dscalar2> AreaPeriPreferences;//(A0,P0) for each cell
        GPUArray<Dscalar2> AreaPeri;//(current A,P) for each cell
        GPUArray<Dscalar2> Moduli;//(KA,KP)
        GPUArray<Dscalar2> Motility;//(v0,Dr) for each cell

        GPUArray<Dscalar> cellDirectors_initial;// for testing
        GPUArray<Dscalar2> displacements;

//        curandState *devStates;
        GPUArray<curandState> devStates;

        //delSet.data[n_idx(nn,i)] are four consecutive delaunay neighbors, orientationally ordered, of point i (for use in computing forces on GPU)
        GPUArray<int4> delSets;
        //delOther.data[n_idx(nn,i)] contains the index of the "other" delaunay neighbor. i.e., the mutual neighbor of delSet.data[n_idx(nn,i)].y and delSet.data[n_idx(nn,i)].z that isn't point i
        GPUArray<int> delOther;

        //arrays indexed by (nn, pidx) of that particle and neighbor number's voronoi vertex (in order) and, and the last and next voro points
        GPUArray<Dscalar2> VoroCur;
        GPUArray<Dscalar4> VoroLastNext;

        //interactions are computed "per voronoi vertex"...forceSets are summed up to get total force on a particle
        GPUArray<Dscalar2> forceSets;

    public:
        int Timestep;
        int sortPeriod;
        bool spatialSortThisStep;
        Dscalar deltaT;
        GPUArray<int> CellType;
        GPUArray<Dscalar> cellDirectors;
        GPUArray<Dscalar2> forces;

        //"exclusiosn" zero out the force on a cell...the external force needed to do this is stored in external_forces
        GPUArray<Dscalar2> external_forces;
        GPUArray<int> exclusions;

        ~SPV2D()
            {
            };
        //initialize with random positions in a square box
        SPV2D(int n);
        //additionally set all cells to have uniform target A_0 and P_0 parameters
        SPV2D(int n, Dscalar A0, Dscalar P0);

        //initialize DelaunayMD, and set random orientations for cell directors
        void Initialize(int n);

        //set and get
        void setDeltaT(Dscalar dt){deltaT = dt;};
        ///the following set uniform motilities...for individual choices use setCellMotility
        void setv0Dr(Dscalar v0new,Dscalar drnew);

        void setTension(Dscalar g){gamma = g;};
        void setUseTension(bool u){useTension = u;};


        void setSortPeriod(int sp){sortPeriod = sp;};
        void setCellPreferencesUniform(Dscalar A0, Dscalar P0);
        void setModuliUniform(Dscalar KA, Dscalar KP);

        void setCellTypeUniform(int i);
        void setCellType(vector<int> &types);

        void setCellMotility(vector<Dscalar> &v0s,vector<Dscalar> &drs);

        //sets particles within an ellipse to type 0, outside to type 1. frac is fraction of area for the ellipse to take up, aspectRatio is (r_x/r_y)
        void setCellTypeEllipse(Dscalar frac, Dscalar aspectRatio);
        //sets particles within a strip (surface normal to x) to type 0, other particles to type 1. Fraction is the area of strip occupied by the system
        void setCellTypeStrip(Dscalar frac);

        void setCurandStates(int i);

        //set exclusions...if a particle is excluded (ex[idx]=1) then its force is zeroed out (the external force to do this is stored in excluded_forces) and its motility is set to zero
        void setExclusions(vector<int> &exes);

        //utility
        void getDelSets(int i);
        void allDelSets();
        void centerCells();
        void spatialSorting();

        //cell-dynamics related functions
        void performTimestep();
        void performTimestepCPU();
        void performTimestepGPU();


        //CPU functions
        void computeGeometryCPU();
        void computeSPVForceCPU(int i);
        void computeSPVForceWithTensionsCPU(int i,bool verbose = false);
        void calculateDispCPU();


        //GPU functions
        void DisplacePointsAndRotate();
        void computeGeometryGPU();
        void computeSPVForceSetsGPU();
        void computeSPVForceSetsWithTensionsGPU();
        void sumForceSets();
        void sumForceSetsWithExclusions();



        //


        //testing functions...
        void reportCellInfo();
        void reportForces();
        void reportDirectors();
        void meanForce();
        void meanArea();
        Dscalar reportq();
        void deltaAngle();

        Dscalar triangletiming, forcetiming;
    };





#endif
