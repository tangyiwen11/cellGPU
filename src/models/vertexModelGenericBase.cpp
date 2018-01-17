#define ENABLE_CUDA

#include "vertexModelGenericBase.h"
#include "DelaunayCGAL.h" 
#include "vertexModelGenericBase.cuh"
/*! \file vertexModelGenericBase.cpp */

/*!
Take care of all base class initialization functions, this involves setting arrays to the right size, etc.
*/
void vertexModelGenericBase::initializeVertexGenericModelBase(int n)
    {
    //call initializer chain...sets Ncells = n
    //Note that this, by default, uses a random voronoi construction with every vertex being three-fold coordinated
    initializeSimpleVertexModelBase(n);
    vertexCoordinationMaximum = 3;
    printf("vertices = %i\n", Nvertices);
    vertexNeighborIndexer = Index2D(vertexCoordinationMaximum,Nvertices);
    vector<int> vertexCoordination(Nvertices,3);
    fillGPUArrayWithVector(vertexCoordination,vertexNeighborNum);
    fillGPUArrayWithVector(vertexCoordination,vertexCellNeighborNum);
    //initializeEdgeFlipLists();//NOTE THAT when this is written it must be added to remove cells and any other topology-changing functions

    //growCellVertexListAssist.resize(1);
    //ArrayHandle<int> h_grow(growCellVertexListAssist,access_location::host,access_mode::overwrite);
    //h_grow.data[0]=0;
    };

/*!
When a transition increases the maximum coordination of the vertices in the system,
call this function first to copy over the vertexNeighbors and vertexCellNeighbors structures
and resize the forceSets, voroCur, and voroLastNext structures
 */
void vertexModelGenericBase::growVertexNeighborLists(int newCoordinationMax)
    {
    cout << "maximum vertex coordination grew from " << vertexCoordinationMaximum << " to " << newCoordinationMax << endl;
    vertexCoordinationMaximum = newCoordinationMax;

    Index2D old_indexer = vertexNeighborIndexer;
    vertexNeighborIndexer = Index2D(vertexCoordinationMaximum,Nvertices);

    copyReIndexed2DGPUArray(vertexNeighborNum,old_indexer,vertexNeighborIndexer);
    copyReIndexed2DGPUArray(vertexNeighbors,old_indexer,vertexNeighborIndexer);
    copyReIndexed2DGPUArray(vertexCellNeighborNum,old_indexer,vertexNeighborIndexer);
    copyReIndexed2DGPUArray(vertexCellNeighbors,old_indexer,vertexNeighborIndexer);

    //resize per-vertex-coordination lists
    vertexForceSets.resize(vertexCoordinationMaximum*Nvertices);
    voroCur.resize(vertexCoordinationMaximum*Nvertices);
    voroLastNext.resize(vertexCoordinationMaximum*Nvertices);
    };

/*!
Compute the area and perimeter of every cell on the CPU
*/
void vertexModelGenericBase::computeGeometryCPU()
    {
    ArrayHandle<Dscalar2> h_v(vertexPositions,access_location::host,access_mode::read);
    ArrayHandle<int> h_cvn(cellVertexNum,access_location::host,access_mode::read);
    ArrayHandle<int> h_cv(cellVertices,access_location::host,access_mode::read);
    ArrayHandle<int> h_vcn(vertexCellNeighbors,access_location::host,access_mode::read);
    ArrayHandle<Dscalar2> h_vc(voroCur,access_location::host,access_mode::readwrite);
    ArrayHandle<Dscalar4> h_vln(voroLastNext,access_location::host,access_mode::readwrite);
    ArrayHandle<Dscalar2> h_AP(AreaPeri,access_location::host,access_mode::readwrite);
    ArrayHandle<int> h_vcnn(vertexCellNeighborNum,access_location::host,access_mode::read);

    //compute the geometry for each cell
    for (int i = 0; i < Ncells; ++i)
        {
        int neighs = h_cvn.data[i];
//      Define the vertices of a cell relative to some (any) of its verties to take care of periodic boundaries
        Dscalar2 cellPos = h_v.data[h_cv.data[cellNeighborIndexer(neighs-2,i)]];
        Dscalar2 vlast, vcur,vnext;
        Dscalar Varea = 0.0;
        Dscalar Vperi = 0.0;
        //compute the vertex position relative to the cell position
        vlast.x=0.;vlast.y=0.0;
        int vidx = h_cv.data[cellNeighborIndexer(neighs-1,i)];
        Box->minDist(h_v.data[vidx],cellPos,vcur);
        for (int nn = 0; nn < neighs; ++nn)
            {
            //for easy force calculation, save the current, last, and next vertex position in the approprate spot.
            int forceSetIdx= -1;
            int vertexCoordination = h_vcnn.data[vidx];
            for (int ff = 0; ff < vertexCoordination; ++ff)
                if(h_vcn.data[vertexNeighborIndexer(ff,vidx)] == i)
                    forceSetIdx = vertexNeighborIndexer(ff,vidx);
            vidx = h_cv.data[cellNeighborIndexer(nn,i)];
            Box->minDist(h_v.data[vidx],cellPos,vnext);
            //contribution to cell's area is
            // 0.5* (vcur.x+vnext.x)*(vnext.y-vcur.y)
            Varea += SignedPolygonAreaPart(vcur,vnext);
            Dscalar dx = vcur.x-vnext.x;
            Dscalar dy = vcur.y-vnext.y;
            Vperi += sqrt(dx*dx+dy*dy);
            //save vertex positions in a convenient form
            h_vc.data[forceSetIdx] = vcur;
            h_vln.data[forceSetIdx] = make_Dscalar4(vlast.x,vlast.y,vnext.x,vnext.y);
            //advance the loop
            vlast = vcur;
            vcur = vnext;
            };
        h_AP.data[i].x = Varea;
        h_AP.data[i].y = Vperi;
        };
    };  

/*!
Remove any cell whose index is in the list. Also discard any vertices which are now no longer part
of a cell. This then requires some laborious -- but straightforward -- reindexing and resizing of
GPUArrays
*/
void vertexModelGenericBase::removeCells(vector<int> cellIndices)
    {
    //figure out the cells and vertices left
    sort(cellIndices.begin(), cellIndices.end());
    vector<int> verticesToRemove;
    vector<int> newVertexNeighborNumber;
    copyGPUArrayData(vertexNeighborNum,newVertexNeighborNumber);
    
    {//array handle scope
    ArrayHandle<int> cvn(cellVertexNum);
    ArrayHandle<int> cv(cellVertices);
    for (int cc = 0; cc <cellIndices.size();++cc)
        {
        int cIdx = cellIndices[cc];
        int vNeighs = cvn.data[cellIndices[cc]];
        for (int vv = 0; vv < vNeighs; ++vv)
            {
            int vidx = cv.data[cellNeighborIndexer(vv,cIdx)];
            newVertexNeighborNumber[vidx] -= 1;
            if (newVertexNeighborNumber[vidx] <= 0)
                verticesToRemove.push_back(vidx);
            };
        };
    }//array handle scope
    sort(verticesToRemove.begin(),verticesToRemove.end());
    //Now, create a map of (old index, new index) for cells and vertices, with new index = -1 if it is to be removed
    //so, e.g., cellMap[100] is the new index of what used to be cell 100, or is -1 if it was to be removed
    //cellMapInverse is a vector of length new Ncells. cellMapInverse[10] gives the index of the cell that will be mapped to the 10th cell
    createIndexMapAndInverse(cellMap,cellMapInverse,cellIndices,Ncells);
    createIndexMapAndInverse(vertexMap,vertexMapInverse,verticesToRemove,Nvertices);

    //Great... now we need to correct basically all of the data structures
    int NvOld = Nvertices;
    int NcOld = Ncells;
    Nvertices = NvOld - verticesToRemove.size();
    Ncells = NcOld - cellIndices.size();

    //first, let's handle the arrays that don't require any logic
    removeGPUArrayElement(vertexPositions,verticesToRemove);
    removeGPUArrayElement(vertexMasses,verticesToRemove);
    removeGPUArrayElement(vertexForces,verticesToRemove);
    removeGPUArrayElement(vertexVelocities,verticesToRemove);
    removeGPUArrayElement(displacements,verticesToRemove);

    removeGPUArrayElement(cellDirectors,cellIndices);
    removeGPUArrayElement(Motility,cellIndices);
    removeGPUArrayElement(AreaPeri,cellIndices);
    removeGPUArrayElement(AreaPeriPreferences,cellIndices);
    removeGPUArrayElement(cellForces,cellIndices);
    removeGPUArrayElement(cellMasses,cellIndices);
    removeGPUArrayElement(cellVelocities,cellIndices);
    removeGPUArrayElement(Moduli,cellIndices);
    removeGPUArrayElement(cellType,cellIndices);
    removeGPUArrayElement(cellPositions,cellIndices);

    //let's handle the vertex-based structures
    vertexNeighborIndexer = Index2D(vertexCoordinationMaximum,Nvertices);
    Index2D oldVNI = Index2D(vertexCoordinationMaximum,NvOld);
    GPUArray<int> newVertexNeighbors, newVertexNeighborNum, newVertexCellNeighbors, newVertexCellNeighborNum;
    newVertexNeighbors.resize(vertexCoordinationMaximum*Nvertices);
    newVertexCellNeighbors.resize(vertexCoordinationMaximum*Nvertices);
    newVertexNeighborNum.resize(Nvertices);
    newVertexCellNeighborNum.resize(Nvertices);
    {//scope for array handles
    ArrayHandle<int> vn(newVertexNeighbors); ArrayHandle<int> vnOld(vertexNeighbors);
    ArrayHandle<int> vcn(newVertexCellNeighbors); ArrayHandle<int> vcnOld(vertexCellNeighbors);
    ArrayHandle<int> vnn(newVertexNeighborNum); ArrayHandle<int> vnnOld(vertexNeighborNum);
    ArrayHandle<int> vcnn(newVertexCellNeighborNum); ArrayHandle<int> vcnnOld(vertexCellNeighborNum);
    for (int vv = 0; vv < Nvertices; ++vv)
        {
        int oldVidx = vertexMapInverse[vv];
        //first, do the vertex neighbors of new vertex vv
        int neighs = vnnOld.data[oldVidx];
        int neighCur = 0;
        for (int nn = 0; nn < neighs; ++nn)
            {
            int neighbor = vnOld.data[oldVNI(nn,oldVidx)];
            if(vertexMap[neighbor] != -1)
                {
                vn.data[vertexNeighborIndexer(neighCur,vv)] = vertexMap[neighbor];
                neighCur +=1;
                };
            };
        vnn.data[vv] = neighCur;
        //now do the cells that vertex vv borders
        neighCur = 0;
        neighs = vcnnOld.data[oldVidx];
        for (int nn = 0; nn < neighs; ++nn)
            {
            int cellNeighbor = vcnOld.data[oldVNI(nn,oldVidx)];
            if(cellMap[cellNeighbor] != -1)
                {
                vcn.data[vertexNeighborIndexer(neighCur,vv)] = cellMap[cellNeighbor];
                neighCur += 1;
                };
            };
            vcnn.data[vv] = neighCur;
        };
    }//scope for array handles
    vertexNeighbors.swap(newVertexNeighbors);
    vertexCellNeighbors.swap(newVertexCellNeighbors);
    vertexNeighborNum.swap(newVertexNeighborNum);
    vertexCellNeighborNum.swap(newVertexCellNeighborNum);

    vertexForceSets.resize(vertexCoordinationMaximum*Nvertices);
    voroCur.resize(vertexCoordinationMaximum*Nvertices);
    voroLastNext.resize(vertexCoordinationMaximum*Nvertices);

    
    //Great! Now the cell-based lists

    cellNeighborIndexer = Index2D(vertexMax,Ncells);  
    Index2D oldCNI = Index2D(vertexMax,NcOld);
    GPUArray<int> newCellVertices,newCellVertexNum;
    newCellVertices.resize(vertexMax*Ncells);
    newCellVertexNum.resize(Ncells);
    {//scope for array handles
    ArrayHandle<int> cv(newCellVertices); ArrayHandle<int> cvOld(cellVertices);
    ArrayHandle<int> cvn(newCellVertexNum); ArrayHandle<int> cvnOld(cellVertexNum);
    for (int cc = 0; cc < Ncells; ++cc)
        {
        int oldCidx = cellMapInverse[cc];
        int neighs = cvnOld.data[oldCidx];
        int neighCur = 0;
        for (int nn = 0; nn < neighs; ++nn)
            {
            int neighbor = cvOld.data[oldCNI(nn,oldCidx)];
            if(vertexMap[neighbor] != -1)
                {
                cv.data[cellNeighborIndexer(neighCur,cc)] = vertexMap[neighbor];
                neighCur +=1;
                };
            };
        cvn.data[cc] = neighCur;
        };
    }//scope for array handles

    cellVertexNum.swap(newCellVertexNum);
    cellVertices.swap(newCellVertices);

    //handle the spatial sorting arrays
    vector<int> newitt, newtti, newidxtotag, newtagtoidx;
    newitt.resize(Ncells);newtti.resize(Ncells);newidxtotag.resize(Ncells);newtagtoidx.resize(Ncells);
    for (int cc = 0; cc < Ncells; ++cc)
        {
        int oldCidx = cellMapInverse[cc];
        newitt[cc] = cellMap[itt[oldCidx]];
        newtti[cc] = cellMap[tti[oldCidx]];
        newidxtotag[cc] = cellMap[idxToTag[oldCidx]];
        newtagtoidx[cc] = cellMap[tagToIdx[oldCidx]];
        };
    itt=newitt;
    tti=newtti;
    idxToTag=newidxtotag;
    tagToIdx=newtagtoidx;

    vector<int> newittv, newttiv, newidxtotagv, newtagtoidxv;
    newittv.resize(Nvertices);newttiv.resize(Nvertices);newidxtotagv.resize(Nvertices);newtagtoidxv.resize(Nvertices);
    for (int vv = 0; vv < Nvertices; ++vv)
        {
        int oldVidx = vertexMapInverse[vv];
        newittv[vv] = vertexMap[ittVertex[oldVidx]];
        newttiv[vv] = vertexMap[ttiVertex[oldVidx]];
        newidxtotagv[vv] = vertexMap[idxToTagVertex[oldVidx]];
        newtagtoidxv[vv] = vertexMap[tagToIdxVertex[oldVidx]];
        };
    ittVertex=newittv;
    ttiVertex=newttiv;
    idxToTagVertex=newidxtotagv;
    tagToIdxVertex=newtagtoidxv;
    
    //edgeFlipLists
    computeGeometry();
    };

/*!
Take a vector of vertex indices and merge them into a single vertex, handling all of the data structure
updates necessary.
Of note, this will leave a single vertex at the position of the lowest-indexed vertex in the set of
vertices that is to be merged.
*/
void vertexModelGenericBase::mergeVertices(vector<int> verticesToMerge)
    {
    //for convenience, let's decide to keep the lowest-indexed vertex...for super-convenience, we'll sort the vector
    sort(verticesToMerge.begin(),verticesToMerge.end());

    int vertexIdx = verticesToMerge[0];
    //create vertexMap and vertexMapIndex
    vector<int> verticesToRemove(verticesToMerge.begin()+1,verticesToMerge.end());
    createIndexMapAndInverse(vertexMap,vertexMapInverse,verticesToRemove,Nvertices);

    //create a list of the new vertex and cell neighbors of the merged vertex... for convenience, this list contains the future vertex indexing, not the current one
    vector<int> vNeighbors, cNeighbors;
    int vNeighNum = 0;
    int cNeighNum = 0;
    vNeighbors.reserve(vertexCoordinationMaximum*verticesToMerge.size());
    cNeighbors.reserve(vertexCoordinationMaximum*verticesToMerge.size());
    {//array handle scope
    ArrayHandle<int> vn(vertexNeighbors);
    ArrayHandle<int> vnn(vertexNeighborNum);
    ArrayHandle<int> vcn(vertexCellNeighbors);
    ArrayHandle<int> vcnn(vertexCellNeighborNum);
    for (int vv = 0; vv < verticesToMerge.size(); ++vv)
        {
        int vidx = verticesToMerge[vv];
        int cneigh = vcnn.data[vidx];
        for (int nn = 0; nn < cneigh; ++nn)
            {
            int neighbor = vcn.data[vertexNeighborIndexer(nn,vidx)];
            cNeighbors.push_back(neighbor);
            cNeighNum +=1;
            };
        int vneigh = vnn.data[vidx];
        for (int nn = 0; nn < vneigh; ++nn)
            {
            int neighbor = vn.data[vertexNeighborIndexer(nn,vidx)];
            if (vertexMap[neighbor] != -1)
                {
                vNeighbors.push_back(vertexMap[neighbor]);
                vNeighNum +=1;
                }
            };
        };
    };
    removeDuplicateVectorElements(vNeighbors);
    removeDuplicateVectorElements(cNeighbors);

    if(vNeighNum > vertexCoordinationMaximum)
        {
        vertexCoordinationMaximum = vNeighNum;
        growVertexNeighborLists(vertexCoordinationMaximum);
        };


    //resize data structures as necessary...
    //note that cellNeighborIndexer will not need to change, but verteNeighborIndexer might

    //the new, combined vertex will be neighbors with all of the cells of the original vertex

    //each cell that neighbored any of the merged vertices is neighbors only with the remaining vertex

    //update 
    


    };

/*!
cell death, as opposed to cell removal, removes a cell but maintains the confluent nature of the
tissue. This is accomplished by subsequent calls to removeCell and then mergeVertices.
*/
void vertexModelGenericBase::cellDeath(int cellIndex)
    {
    //get a list of vertices that make up the cell
    //call removeCell
    //use vertexMap and vertexMapInverse to figure out the new vertex labels
    //merge those vertices together
    };
