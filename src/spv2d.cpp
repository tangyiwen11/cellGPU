using namespace std;

//definitions needed for DelaunayLoc, voroguppy namespace, Triangle, Triangle, and all GPU functions, respectively

#define EPSILON 1e-12
#define dbl float
#define REAL double
#define ANSI_DECLARATIONS
#define ENABLE_CUDA

#define PI 3.14159265358979323846

#include "spv2d.h"

SPV2D::SPV2D(int n)
    {
    printf("Initializing %i cells with random positions in a square box\n",n);
    Initialize(n);
    };

SPV2D::SPV2D(int n,float A0, float P0)
    {
    printf("Initializing %i cells with random positions in a square box\n",n);
    Initialize(n);
    setCellPreferencesUniform(A0,P0);
    };

void SPV2D::Initialize(int n)
    {
    setDeltaT(0.01);
    setDr(1.);
    initialize(n);
    forces.resize(n);
    AreaPeri.resize(n);
    cellDirectors.resize(n);
    ArrayHandle<float2> h_cd(cellDirectors,access_location::host, access_mode::overwrite);
    int randmax = 100000000;
    for (int ii = 0; ii < N; ++ii)
        {
        float theta = PI/(float)(randmax)* (float)(rand()%randmax);
        h_cd.data[ii].x = cos(theta);
        h_cd.data[ii].y = sin(theta);
        };
    };

void SPV2D::setCellPreferencesUniform(float A0, float P0)
    {
    AreaPeriPreferences.resize(N);
    ArrayHandle<float2> h_p(AreaPeriPreferences,access_location::host,access_mode::overwrite);
    for (int ii = 0; ii < N; ++ii)
        {
        h_p.data[ii].x = A0;
        h_p.data[ii].y = P0;
        };
    };

void SPV2D::computeSPVForces()
    {


    };

void SPV2D::performTimestep()
    {
    computeSPVForces();
    //vector of displacements is forces*timestep + v0's*timestep
    GPUArray<float2> ds; ds.resize(N);
    repel(ds,1e-3);

    movePoints(ds);
    //need to re-write a new movepoints that moves points and rotates cell directors
    testAndRepairTriangulation();

    };


void SPV2D::computeSPVForceCPU(int i)
    {
    printf("cell %i: \n",i);
    //for testing these routines...
    vector <int> test;
    DelaunayCell celltest;
    delLoc.triangulatePoint(i, test,celltest);

    //read in all the data we'll need
    ArrayHandle<float2> h_p(points,access_location::host,access_mode::read);
    ArrayHandle<float2> h_f(forces,access_location::host,access_mode::readwrite);
    ArrayHandle<float2> h_AP(AreaPeri,access_location::host,access_mode::readwrite);
    ArrayHandle<float2> h_APpref(AreaPeriPreferences,access_location::host,access_mode::read);

    ArrayHandle<int> h_nn(neigh_num,access_location::host,access_mode::read);
    ArrayHandle<int> h_n(neighs,access_location::host,access_mode::read);

    //get Delaunay neighbors of the cell
    int neigh = h_nn.data[i];
    vector<int> ns(neigh);
    for (int nn = 0; nn < neigh; ++nn)
        ns[nn]=h_n.data[n_idx(nn,i)];


    //compute base set of voronoi points
    vector<float2> voro(neigh);
    float2 circumcent;
    float2 origin; origin.x = 0.; origin.y=0.;
    float2 nnext,nlast;
    float2 nnextp,nlastp;
    float2 pi = h_p.data[i];

    nlastp = h_p.data[ns[ns.size()-1]];
    Box.minDist(nlastp,pi,nlast);
    for (int nn = 0; nn < neigh;++nn)
        {
        nnextp = h_p.data[ns[nn]];
        Box.minDist(nnextp,pi,nnext);
        Circumcenter(origin,nlast,nnext,circumcent);
        voro[nn] = circumcent;
        nlast=nnext;
        };

    float2 vlast,vnext;
    //compute Area and perimeter
    float Varea = 0.0;
    float Vperi = 0.0;
    vlast = voro[neigh-1];
    for (int nn = 0; nn < neigh; ++nn)
        {
        vnext=voro[nn];
        Varea += TriangleArea(vlast,vnext);
        float dx = vlast.x-vnext.x;
        float dy = vlast.y-vnext.y;
        Vperi += sqrt(dx*dx+dy*dy);
        vlast=vnext;
        };
    h_AP.data[i].x = Varea;
    h_AP.data[i].y = Vperi;

    };



void SPV2D::meanArea()
    {
    ArrayHandle<float2> h_AP(AreaPeri,access_location::host,access_mode::read);
    float Am = 0.0;
    for (int i = 0; i < N; ++i)
        Am += h_AP.data[i].x/N;
    printf("Mean area = %f\n" ,Am);

    };

