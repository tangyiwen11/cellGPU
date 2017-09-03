#define ENABLE_CUDA

#include "avm2d.h"
#include "avm2d.cuh"
#include "voronoi2d.h"
/*! \file avm2d.cpp */

/*!
\param n number of CELLS to initialize
\param A0 set uniform preferred area for all cells
\param P0 set uniform preferred perimeter for all cells
\param reprod should the simulation be reproducible (i.e. call a RNG with a fixed seed)
\param runSPVToInitialize the default constructor has the cells start as a Voronoi tesselation of
a random point set. Set this flag to true to relax this initial configuration via the Voronoi2D class
\post Initialize(n,runSPVToInitialize) is called, setCellPreferencesUniform(A0,P0), and
setModuliUniform(1.0,1.0)
*/
AVM2D::AVM2D(int n,Dscalar A0, Dscalar P0,bool reprod,bool runSPVToInitialize) :
    T1Threshold(0.01)
    {
    printf("Initializing %i cells with random positions as an initially Delaunay configuration in a square box... \n",n);
    Reproducible = reprod;
    GPUcompute=true;
    Initialize(n,runSPVToInitialize);
    setCellPreferencesUniform(A0,P0);
    setModuliUniform(1.0,1.0);
    setCellTypeUniform(0);
    };

/*!
Take care of all class initialization functions, this involves setting arrays to the right size, etc.
*/
void AVM2D::Initialize(int n,bool spvInitialize)
    {
    //set number of cells, and a square box
    Ncells=n;
    cellPositions.resize(Ncells);
    //put cells in box randomly
    setCellPositionsRandomly();
    //derive the vertices from a voronoi tesselation
    setCellsVoronoiTesselation(spvInitialize);
    setCellDirectorsRandomly();

    Timestep = 0;
    setDeltaT(0.01);
    setT1Threshold(0.01);

    //initializes per-cell lists
    AreaPeri.resize(Ncells);
    AreaPeriPreferences.resize(Ncells);
    initializeCellSorting();

    //initializes per-vertex lists
    vertexForces.resize(Nvertices);
    displacements.resize(Nvertices);
    initializeVertexSorting();
    initializeEdgeFlipLists();
    //initialize per-triple-vertex lists
    vertexForceSets.resize(3*Nvertices);
    voroCur.resize(3*Nvertices);
    voroLastNext.resize(3*Nvertices);

    growCellVertexListAssist.resize(1);
    ArrayHandle<int> h_grow(growCellVertexListAssist,access_location::host,access_mode::overwrite);
    h_grow.data[0]=0;
    };

/*!
A function of convenience.... initialize cell positions and vertices by constructing the Delaunay
triangulation of the current cell positions. If you want something more regular, run the Voronoi mode for a few
timesteps to smooth out the random point set first.
\param spvInitialize only use if the initial cell positions are to be random, and you want to make the points more uniform
\post After this is called, all topology data structures are initialized
*/
void AVM2D::setCellsVoronoiTesselation(bool spvInitialize)
    {
    ArrayHandle<Dscalar2> h_p(cellPositions,access_location::host,access_mode::readwrite);
    //use the Voronoi class to relax the initial configuration just a bit?
    if(spvInitialize)
        {
        EOMPtr spp = make_shared<selfPropelledParticleDynamics>(Ncells);

        ForcePtr spv = make_shared<Voronoi2D>(Ncells,1.0,3.8,Reproducible);
        spv->setCellPreferencesUniform(1.0,3.8);
        spv->setv0Dr(.1,1.0);

        SimulationPtr sim = make_shared<Simulation>();
        sim->setConfiguration(spv);
        sim->addUpdater(spp,spv);
        sim->setIntegrationTimestep(0.1);
        sim->setCPUOperation(true);
        spp->setDeltaT(0.1);
        sim->setReproducible(true);

        for (int ii = 0; ii < 100;++ii)
            sim->performTimestep();
        ArrayHandle<Dscalar2> h_pp(spv->cellPositions,access_location::host,access_mode::read);
        for (int ii = 0; ii < Ncells; ++ii)
            h_p.data[ii] = h_pp.data[ii];
        };

    //call CGAL to get Delaunay triangulation
    vector<pair<Point,int> > Psnew(Ncells);
    for (int ii = 0; ii < Ncells; ++ii)
        {
        Psnew[ii]=make_pair(Point(h_p.data[ii].x,h_p.data[ii].y),ii);
        };
    Dscalar b11,b12,b21,b22;
    Box.getBoxDims(b11,b12,b21,b22);
    Iso_rectangle domain(0.0,0.0,b11,b22);
    PDT T(Psnew.begin(),Psnew.end(),domain);
    T.convert_to_1_sheeted_covering();

    //set number of vertices
    Nvertices = 2*Ncells;
    vertexPositions.resize(Nvertices);
    ArrayHandle<Dscalar2> h_v(vertexPositions,access_location::host,access_mode::overwrite);

    //first, ask CGAL for the circumcenter of the face, and add it to the list of vertices, and make a map between the iterator and the vertex idx
    map<PDT::Face_handle,int> faceToVoroIdx;
    int idx = 0;
    for(PDT::Face_iterator fit = T.faces_begin(); fit != T.faces_end(); ++fit)
        {
        PDT::Point p(T.dual(fit));
        h_v.data[idx].x = p.x();
        h_v.data[idx].y = p.y();
        faceToVoroIdx[fit] = idx;
        idx +=1;
        };
    //create a list of what vertices are connected to what vertices,
    //and what cells each vertex is part of
    vertexNeighbors.resize(3*Nvertices);
    vertexCellNeighbors.resize(3*Nvertices);
    ArrayHandle<int> h_vn(vertexNeighbors,access_location::host,access_mode::overwrite);
    ArrayHandle<int> h_vcn(vertexCellNeighbors,access_location::host,access_mode::overwrite);
    for(PDT::Face_iterator fit = T.faces_begin(); fit != T.faces_end(); ++fit)
        {
        int vidx = faceToVoroIdx[fit];
        for(int ff =0; ff<3; ++ff)
            {
            PDT::Face_handle neighFace = fit->neighbor(ff);
            int vnidx = faceToVoroIdx[neighFace];
            h_vn.data[3*vidx+ff] = vnidx;
            h_vcn.data[3*vidx+ff] = fit->vertex(ff)->info();
            };
        };

    //now create a list of what vertices are associated with each cell
    //first get the maximum number of vertices for a cell, and the number of vertices per cell
    cellVertexNum.resize(Ncells);
    ArrayHandle<int> h_cvn(cellVertexNum,access_location::host,access_mode::overwrite);
    vertexMax = 0;
    int nnum = 0;
    for(PDT::Vertex_iterator vit = T.vertices_begin(); vit != T.vertices_end(); ++vit)
        {
        PDT::Vertex_circulator vc(vit);
        int base = vc ->info();
        int neighs = 1;
        ++vc;
        while(vc->info() != base)
            {
            neighs += 1;
            ++vc;
            };
        h_cvn.data[vit->info()] = neighs;
        if (neighs > vertexMax) vertexMax = neighs;
        nnum += neighs;
        };
    vertexMax += 2;
    cout << "Total number of neighs = " << nnum << endl;
    cellVertices.resize(vertexMax*Ncells);
    n_idx = Index2D(vertexMax,Ncells);

    //now use face circulators and the map to get the vertices associated with each cell
    ArrayHandle<int> h_cv(cellVertices,access_location::host, access_mode::overwrite);
    for(PDT::Vertex_iterator vit = T.vertices_begin(); vit != T.vertices_end(); ++vit)
        {
        int cellIdx = vit->info();
        PDT::Face_circulator fc(vit);
        int fidx = 0;
        for (int ff = 0; ff < h_cvn.data[vit->info()]; ++ff)
            {
            h_cv.data[n_idx(fidx,cellIdx)] = faceToVoroIdx[fc];
            ++fidx;
            ++fc;
            };
        };
   };

/*!
 *When sortPeriod < 0 this routine does not get called
 \post vertices are re-ordered according to a Hilbert sorting scheme, cells are reordered according
 to what vertices they are near, and all data structures are updated
 */
void AVM2D::spatialSorting()
    {
    //the avm class doesn't need to change any other unusual data structures at the moment
    spatiallySortVerticesAndCellActivity();
    };

/*!
Returns the quadratic energy functional:
E = \sum_{cells} K_A(A_i-A_i,0)^2 + K_P(P_i-P_i,0)^2
*/
Dscalar AVM2D::computeEnergy()
    {
    ArrayHandle<Dscalar2> h_AP(AreaPeri,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_APP(AreaPeriPreferences,access_location::host,access_mode::read);
    Energy = 0.0;
    for (int nn = 0; nn  < Ncells; ++nn)
        {
        Energy += KA * (h_AP.data[nn].x-h_APP.data[nn].x)*(h_AP.data[nn].x-h_APP.data[nn].x);
        Energy += KP * (h_AP.data[nn].y-h_APP.data[nn].y)*(h_AP.data[nn].y-h_APP.data[nn].y);
        };

    return Energy;
    };

/*!
compute the geometry and the forces and the vertices, on either the GPU or CPU as determined by
flags
*/
void AVM2D::computeForces()
    {
    //compute the current area and perimeter of every cell
    //use this information to compute the net force on the vertices
    if(GPUcompute)
        {
        computeGeometryGPU();
        computeForcesGPU();
        }
    else
        {
        computeGeometryCPU();
        computeForcesCPU();
        };
    };

/*!
enforce and update topology of vertex wiring on either the GPU or CPU
*/
void AVM2D::enforceTopology()
    {
    if(GPUcompute)
        {
        //see if vertex motion leads to T1 transitions...ONLY allow one transition per vertex and
        //per cell per timestep
        testAndPerformT1TransitionsGPU();
        }
    else
        {
        //see if vertex motion leads to T1 transitions
        testAndPerformT1TransitionsCPU();
        };
    };

/*!
Use the data pre-computed in the geometry routine to rapidly compute the net force on each vertex
*/
void AVM2D::computeForcesCPU()
    {
    ArrayHandle<int> h_vcn(vertexCellNeighbors,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_vc(voroCur,access_location::host,access_mode::read);
    ArrayHandle<Dscalar4> h_vln(voroLastNext,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_AP(AreaPeri,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_APpref(AreaPeriPreferences,access_location::host,access_mode::read);

    ArrayHandle<Dscalar2> h_fs(vertexForceSets,access_location::host, access_mode::overwrite);
    ArrayHandle<Dscalar2> h_f(vertexForces,access_location::host, access_mode::overwrite);

    //first, compute the contribution to the force on each vertex from each of its three cells
    Dscalar2 vlast,vcur,vnext;
    Dscalar2 dEdv;
    Dscalar Adiff, Pdiff;
    for(int fsidx = 0; fsidx < Nvertices*3; ++fsidx)
        {
        int cellIdx = h_vcn.data[fsidx];
        Dscalar Adiff = KA*(h_AP.data[cellIdx].x - h_APpref.data[cellIdx].x);
        Dscalar Pdiff = KP*(h_AP.data[cellIdx].y - h_APpref.data[cellIdx].y);
        vcur = h_vc.data[fsidx];
        vlast.x = h_vln.data[fsidx].x;  vlast.y = h_vln.data[fsidx].y;
        vnext.x = h_vln.data[fsidx].z;  vnext.y = h_vln.data[fsidx].w;

        //computeForceSetAVM is defined in inc/functions.h
        computeForceSetAVM(vcur,vlast,vnext,Adiff,Pdiff,dEdv);


        h_fs.data[fsidx].x = dEdv.x;
        h_fs.data[fsidx].y = dEdv.y;
        };

    //now sum these up to get the force on each vertex
    for (int v = 0; v < Nvertices; ++v)
        {
        Dscalar2 ftemp = make_Dscalar2(0.0,0.0);
        for (int ff = 0; ff < 3; ++ff)
            {
            ftemp.x += h_fs.data[3*v+ff].x;
            ftemp.y += h_fs.data[3*v+ff].y;
            };
        h_f.data[v] = ftemp;
        };
    };

/*!
A utility function for the CPU routine. Given two vertex indices representing an edge that will undergo
a T1 transition, return in the pass-by-reference variables a helpful representation of the cells in the T1
and the vertices to be re-wired...see the comments in "testAndPerformT1TransitionsCPU" for what that representation is
*/
void AVM2D::getCellVertexSetForT1(int vertex1, int vertex2, int4 &cellSet, int4 &vertexSet, bool &growList)
    {
    int cell1,cell2,cell3,ctest;
    int vlast, vcur, vnext, cneigh;
    ArrayHandle<int> h_cv(cellVertices,access_location::host, access_mode::read);
    ArrayHandle<int> h_cvn(cellVertexNum,access_location::host,access_mode::read);
    ArrayHandle<int> h_vcn(vertexCellNeighbors,access_location::host,access_mode::read);
    cell1 = h_vcn.data[3*vertex1];
    cell2 = h_vcn.data[3*vertex1+1];
    cell3 = h_vcn.data[3*vertex1+2];
    //cell_l doesn't contain vertex 1, so it is the cell neighbor of vertex 2 we haven't found yet
    for (int ff = 0; ff < 3; ++ff)
        {
        ctest = h_vcn.data[3*vertex2+ff];
        if(ctest != cell1 && ctest != cell2 && ctest != cell3)
            cellSet.w=ctest;
        };
    //find vertices "c" and "d"
    cneigh = h_cvn.data[cellSet.w];
    vlast = h_cv.data[ n_idx(cneigh-2,cellSet.w) ];
    vcur = h_cv.data[ n_idx(cneigh-1,cellSet.w) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = h_cv.data[n_idx(cn,cell1)];
        if(vcur == vertex2) break;
        vlast = vcur;
        vcur = vnext;
        };

    //classify cell1
    cneigh = h_cvn.data[cell1];
    vlast = h_cv.data[ n_idx(cneigh-2,cell1) ];
    vcur = h_cv.data[ n_idx(cneigh-1,cell1) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = h_cv.data[n_idx(cn,cell1)];
        if(vcur == vertex1) break;
        vlast = vcur;
        vcur = vnext;
        };
    if(vlast == vertex2)
        cellSet.x = cell1;
    else if(vnext == vertex2)
        cellSet.z = cell1;
    else
        {
        cellSet.y = cell1;
        };

    //classify cell2
    cneigh = h_cvn.data[cell2];
    vlast = h_cv.data[ n_idx(cneigh-2,cell2) ];
    vcur = h_cv.data[ n_idx(cneigh-1,cell2) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = h_cv.data[n_idx(cn,cell2)];
        if(vcur == vertex1) break;
        vlast = vcur;
        vcur = vnext;
        };
    if(vlast == vertex2)
        cellSet.x = cell2;
    else if(vnext == vertex2)
        cellSet.z = cell2;
    else
        {
        cellSet.y = cell2;
        };

    //classify cell3
    cneigh = h_cvn.data[cell3];
    vlast = h_cv.data[ n_idx(cneigh-2,cell3) ];
    vcur = h_cv.data[ n_idx(cneigh-1,cell3) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = h_cv.data[n_idx(cn,cell3)];
        if(vcur == vertex1) break;
        vlast = vcur;
        vcur = vnext;
        };
    if(vlast == vertex2)
        cellSet.x = cell3;
    else if(vnext == vertex2)
        cellSet.z = cell3;
    else
        {
        cellSet.y = cell3;
        };

    //get the vertexSet by examining cells j and l
    cneigh = h_cvn.data[cellSet.y];
    vlast = h_cv.data[ n_idx(cneigh-2,cellSet.y) ];
    vcur = h_cv.data[ n_idx(cneigh-1,cellSet.y) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = h_cv.data[n_idx(cn,cellSet.y)];
        if(vcur == vertex1) break;
        vlast = vcur;
        vcur = vnext;
        };
    vertexSet.x=vlast;
    vertexSet.y=vnext;
    cneigh = h_cvn.data[cellSet.w];
    vlast = h_cv.data[ n_idx(cneigh-2,cellSet.w) ];
    vcur = h_cv.data[ n_idx(cneigh-1,cellSet.w) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = h_cv.data[n_idx(cn,cellSet.w)];
        if(vcur == vertex2) break;
        vlast = vcur;
        vcur = vnext;
        };
    vertexSet.w=vlast;
    vertexSet.z=vnext;

    //Does the cell-vertex-neighbor data structure need to be bigger?
    if(h_cvn.data[cellSet.x] == vertexMax || h_cvn.data[cellSet.z] == vertexMax)
        growList = true;
    };

/*!
Test whether a T1 needs to be performed on any edge by simply checking if the edge length is beneath a threshold.
This function also performs the transition and maintains the auxiliary data structures
 */
void AVM2D::testAndPerformT1TransitionsCPU()
    {
    ArrayHandle<Dscalar2> h_v(vertexPositions,access_location::host,access_mode::readwrite);
    ArrayHandle<int> h_vn(vertexNeighbors,access_location::host,access_mode::readwrite);
    ArrayHandle<int> h_cvn(cellVertexNum,access_location::host,access_mode::readwrite);
    ArrayHandle<int> h_cv(cellVertices,access_location::host,access_mode::readwrite);
    ArrayHandle<int> h_vcn(vertexCellNeighbors,access_location::host,access_mode::readwrite);

    Dscalar2 edge;
    //first, scan through the list for any T1 transitions...
    int vertex2;
    //keep track of whether vertexMax needs to be increased
    int vMax = vertexMax;
    /*
     IF v1 is above v2, the following is the convention (otherwise flip CW and CCW)
     cell i: contains both vertex 1 and vertex 2, in CW order
     cell j: contains only vertex 1
     cell k: contains both vertex 1 and vertex 2, in CCW order
     cell l: contains only vertex 2
     */
    int4 cellSet;
    /*
    vertexSet (a,b,c,d) have those indices in which before the transition
    cell i has CCW vertices: ..., c, v2, v1, a, ...
    and
    cell k has CCW vertices: ..., b,v1,v2,d, ...
    */
    int4 vertexSet;
    Dscalar2 v1,v2;
    for (int vertex1 = 0; vertex1 < Nvertices; ++vertex1)
        {
        v1 = h_v.data[vertex1];
        //look at vertexNeighbors list for neighbors of vertex, compute edge length
        for (int vv = 0; vv < 3; ++vv)
            {
            vertex2 = h_vn.data[3*vertex1+vv];
            //only look at each pair once
            if(vertex1 < vertex2)
                {
                v2 = h_v.data[vertex2];
                Box.minDist(v1,v2,edge);
                if(norm(edge) < T1Threshold)
                    {
                    bool growCellVertexList = false;
                    getCellVertexSetForT1(vertex1,vertex2,cellSet,vertexSet,growCellVertexList);
                    //Does the cell-vertex-neighbor data structure need to be bigger?
                    if(growCellVertexList)
                        {
                        vMax +=1;
                        growCellVerticesList(vMax);
                        h_cv = ArrayHandle<int>(cellVertices,access_location::host,access_mode::readwrite);
                        };

                    //Rotate the vertices in the edge and set them at twice their original distance
                    Dscalar2 midpoint;
                    midpoint.x = v2.x + 0.5*edge.x;
                    midpoint.y = v2.y + 0.5*edge.y;

                    v1.x = midpoint.x-edge.y;
                    v1.y = midpoint.y+edge.x;
                    v2.x = midpoint.x+edge.y;
                    v2.y = midpoint.y-edge.x;
                    Box.putInBoxReal(v1);
                    Box.putInBoxReal(v2);
                    h_v.data[vertex1] = v1;
                    h_v.data[vertex2] = v2;

                    //re-wire the cells and vertices
                    //start with the vertex-vertex and vertex-cell  neighbors
                    for (int vert = 0; vert < 3; ++vert)
                        {
                        //vertex-cell neighbors
                        if(h_vcn.data[3*vertex1+vert] == cellSet.z)
                            h_vcn.data[3*vertex1+vert] = cellSet.w;
                        if(h_vcn.data[3*vertex2+vert] == cellSet.x)
                            h_vcn.data[3*vertex2+vert] = cellSet.y;
                        //vertex-vertex neighbors
                        if(h_vn.data[3*vertexSet.y+vert] == vertex1)
                            h_vn.data[3*vertexSet.y+vert] = vertex2;
                        if(h_vn.data[3*vertexSet.z+vert] == vertex2)
                            h_vn.data[3*vertexSet.z+vert] = vertex1;
                        if(h_vn.data[3*vertex1+vert] == vertexSet.y)
                            h_vn.data[3*vertex1+vert] = vertexSet.z;
                        if(h_vn.data[3*vertex2+vert] == vertexSet.z)
                            h_vn.data[3*vertex2+vert] = vertexSet.y;
                        };
                    //now rewire the cells
                    //cell i loses v2 as a neighbor
                    int cneigh = h_cvn.data[cellSet.x];
                    int cidx = 0;
                    for (int cc = 0; cc < cneigh-1; ++cc)
                        {
                        if(h_cv.data[n_idx(cc,cellSet.x)] == vertex2)
                            cidx +=1;
                        h_cv.data[n_idx(cc,cellSet.x)] = h_cv.data[n_idx(cidx,cellSet.x)];
                        cidx +=1;
                        };
                    h_cvn.data[cellSet.x] -= 1;

                    //cell j gains v2 in between v1 and b
                    cneigh = h_cvn.data[cellSet.y];
                    vector<int> cvcopy1(cneigh+1);
                    cidx = 0;
                    for (int cc = 0; cc < cneigh; ++cc)
                        {
                        int cellIndex = h_cv.data[n_idx(cc,cellSet.y)];
                        cvcopy1[cidx] = cellIndex;
                        cidx +=1;
                        if(cellIndex == vertex1)
                            {
                            cvcopy1[cidx] = vertex2;
                            cidx +=1;
                            };
                        };
                    for (int cc = 0; cc < cneigh+1; ++cc)
                        h_cv.data[n_idx(cc,cellSet.y)] = cvcopy1[cc];
                    h_cvn.data[cellSet.y] += 1;

                    //cell k loses v1 as a neighbor
                    cneigh = h_cvn.data[cellSet.z];
                    cidx = 0;
                    for (int cc = 0; cc < cneigh-1; ++cc)
                        {
                        if(h_cv.data[n_idx(cc,cellSet.z)] == vertex1)
                            cidx +=1;
                        h_cv.data[n_idx(cc,cellSet.z)] = h_cv.data[n_idx(cidx,cellSet.z)];
                        cidx +=1;
                        };
                    h_cvn.data[cellSet.z] -= 1;

                    //cell l gains v1 in between v2 and c
                    cneigh = h_cvn.data[cellSet.w];
                    vector<int> cvcopy2(cneigh+1);
                    cidx = 0;
                    for (int cc = 0; cc < cneigh; ++cc)
                        {
                        int cellIndex = h_cv.data[n_idx(cc,cellSet.w)];
                        cvcopy2[cidx] = cellIndex;
                        cidx +=1;
                        if(cellIndex == vertex2)
                            {
                            cvcopy2[cidx] = vertex1;
                            cidx +=1;
                            };
                        };
                    for (int cc = 0; cc < cneigh+1; ++cc)
                        h_cv.data[n_idx(cc,cellSet.w)] = cvcopy2[cc];
                    h_cvn.data[cellSet.w] = cneigh + 1;

                    };//end condition that a T1 transition should occur
                };
            };//end loop over vertex2
        };//end loop over vertices
    };

/*!
call kernels to (1) do force sets calculation, then (2) add them up
*/
void AVM2D::computeForcesGPU()
    {
    ArrayHandle<int> d_vcn(vertexCellNeighbors,access_location::device,access_mode::read);
    ArrayHandle<Dscalar2> d_vc(voroCur,access_location::device,access_mode::read);
    ArrayHandle<Dscalar4> d_vln(voroLastNext,access_location::device,access_mode::read);
    ArrayHandle<Dscalar2> d_AP(AreaPeri,access_location::device,access_mode::read);
    ArrayHandle<Dscalar2> d_APpref(AreaPeriPreferences,access_location::device,access_mode::read);
    ArrayHandle<Dscalar2> d_fs(vertexForceSets,access_location::device, access_mode::overwrite);
    ArrayHandle<Dscalar2> d_f(vertexForces,access_location::device, access_mode::overwrite);

    int nForceSets = voroCur.getNumElements();
    gpu_avm_force_sets(
                    d_vcn.data,
                    d_vc.data,
                    d_vln.data,
                    d_AP.data,
                    d_APpref.data,
                    d_fs.data,
                    nForceSets,
                    KA,
                    KP
                    );

    gpu_avm_sum_force_sets(
                    d_fs.data,
                    d_f.data,
                    Nvertices);
    };

/*!
perform whatever check is desired for T1 transtions (here just a "is the edge too short")
and detect whether the edge needs to grow. If so, grow it!
*/
void AVM2D::testEdgesForT1GPU()
    {
        {//provide scope for array handles
        ArrayHandle<Dscalar2> d_v(vertexPositions,access_location::device,access_mode::read);
        ArrayHandle<int> d_vn(vertexNeighbors,access_location::device,access_mode::read);
        ArrayHandle<int> d_vflip(vertexEdgeFlips,access_location::device,access_mode::overwrite);
        ArrayHandle<int> d_cvn(cellVertexNum,access_location::device,access_mode::read);
        ArrayHandle<int> d_cv(cellVertices,access_location::device,access_mode::read);
        ArrayHandle<int> d_vcn(vertexCellNeighbors,access_location::device,access_mode::read);
        ArrayHandle<int> d_grow(growCellVertexListAssist,access_location::device,access_mode::readwrite);

        //first, test every edge, and check if the cellVertices list needs to be grown
        gpu_avm_test_edges_for_T1(d_v.data,
                              d_vn.data,
                              d_vflip.data,
                              d_vcn.data,
                              d_cvn.data,
                              d_cv.data,
                              Box,
                              T1Threshold,
                              Nvertices,
                              vertexMax,
                              d_grow.data,
                              n_idx);
        }
    ArrayHandle<int> h_grow(growCellVertexListAssist,access_location::host,access_mode::readwrite);
    if(h_grow.data[0] ==1)
        {
        h_grow.data[0]=0;
        growCellVerticesList(vertexMax+1);
        };
    };

/*!
  Iterate through the vertexEdgeFlips list, selecting at most one T1 transition per cell to be done
  on each iteration, until all necessary T1 events have bee performed.
 */
void AVM2D::flipEdgesGPU()
    {
    bool keepFlipping = true;
    //By construction, this loop must always run at least twice...save one of the memory transfers
    int iterations = 0;
    while(keepFlipping)
        {
            {//provide scope for ArrayHandles
            ArrayHandle<Dscalar2> d_v(vertexPositions,access_location::device,access_mode::readwrite);
            ArrayHandle<int> d_vn(vertexNeighbors,access_location::device,access_mode::readwrite);
            ArrayHandle<int> d_vflip(vertexEdgeFlips,access_location::device,access_mode::readwrite);
            ArrayHandle<int> d_vflipcur(vertexEdgeFlipsCurrent,access_location::device,access_mode::readwrite);
            ArrayHandle<int> d_cvn(cellVertexNum,access_location::device,access_mode::readwrite);
            ArrayHandle<int> d_cv(cellVertices,access_location::device,access_mode::readwrite);
            ArrayHandle<int> d_vcn(vertexCellNeighbors,access_location::device,access_mode::readwrite);
            ArrayHandle<int> d_ffe(finishedFlippingEdges,access_location::device,access_mode::readwrite);

            gpu_avm_flip_edges(d_vflip.data,
                               d_vflipcur.data,
                               d_v.data,
                               d_vn.data,
                               d_vcn.data,
                               d_cvn.data,
                               d_cv.data,
                               d_ffe.data,
                               T1Threshold,
                               Box,
                               n_idx,
                               Nvertices,
                               Ncells);
            }; //scope for arrayhandles
        ArrayHandle<int> h_ffe(finishedFlippingEdges,access_location::host,access_mode::readwrite);
        if(h_ffe.data[0]==0)
            keepFlipping = false;
        h_ffe.data[0]=0;
        iterations += 1;
        };
    };

/*!
Because the cellVertexList might need to grow, it's convenient to break this into two parts
*/
void AVM2D::testAndPerformT1TransitionsGPU()
    {
    testEdgesForT1GPU();
    flipEdgesGPU();
    };
