
/*!
Trigger a cell division event, which involves some laborious re-indexing of various data structures.
This version uses a heavy-handed approach, hitting the cell positions with a full, global retriangulation
(rather than just updating the targeted cell and its neighbor with the local topology repair routines).
If a simulation is needed where the cell division rate is very rapid, this should be improved.
The idea of the division is that a targeted cell will divide normal to an axis specified by the
angle, theta, passed to the function. The final state cell positions are placed along the axis at a
distance away from the initial cell position set by a multiplicative factor (<1) of the in-routine determined
maximum distance in the cell along that axis.
parameters[0] = the index of the cell to undergo a division event
dParams[0] = an angle, theta
dParams[1] = a fraction of maximum separation of the new cell positions along the axis of the cell specified by theta
This function is meant to be called before the start of a new timestep. It should be immediately followed by a computeGeometry call.
\post the new cell is the final indexed entry of the various data structures (e.g.,
cellPositions[new number of cells - 1])
*/
void voronoiModelBase::cellDivision(const vector<int> &parameters, const vector<double> &dParams)
    {
    int cellIdx = parameters[0];
    int cdcellIdx = parameters[1];//(yw)
    double theta = dParams[0];
    double separationFraction = dParams[1];

    //First let's get the geometry of the cell in a convenient reference frame
    //computeGeometry has not yet been called, so need to find the voro positions
    vector<double2> voro;
    voro.reserve(10);
    int neigh;
    double2 initialCellPosition;
    {//arrayHandle scope
    ArrayHandle<double2> h_p(cellPositions,access_location::host,access_mode::read);
    initialCellPosition = h_p.data[cellIdx];
    ArrayHandle<int> h_nn(neighborNum,access_location::host,access_mode::read);
    ArrayHandle<int> h_n(neighbors,access_location::host,access_mode::read);
    neigh = h_nn.data[cellIdx];
    vector<int> ns(neigh);
    for (int nn = 0; nn < neigh; ++nn)
            ns[nn]=h_n.data[n_idx(nn,cellIdx)];
    double2 circumcent;
    double2 nnextp,nlastp;
    double2 pi = h_p.data[cellIdx];
    double2 rij, rik;
    nlastp = h_p.data[ns[ns.size()-1]];
    Box->minDist(nlastp,pi,rij);
    for (int nn = 0; nn < neigh;++nn)
        {
        nnextp = h_p.data[ns[nn]];
        Box->minDist(nnextp,pi,rik);
        Circumcenter(rij,rik,circumcent);
        voro.push_back(circumcent);
        rij=rik;
        }
    };//arrayHandle scope

    //find where the line emanating from the polygons intersects the edges
    double2 c,v1,v2,Int1,Int2;
    bool firstIntFound = false;
    v1 = voro[neigh-1];
    double2 ray; ray.x = Cos(theta); ray.y = Sin(theta);
    c.x = - Sin(theta);
    c.y = Cos(theta);
    double2 p; p.x =0.0; p.y=0.;
    for (int nn = 0; nn < neigh; ++nn)
        {
        v2=voro[nn];
        double2 a; a.x = p.x-v2.x; a.y = p.y-v2.y;
        double2 b; b.x = v1.x - v2.x; b.y = v1.y-v2.y;
        double t2 = -1.;
        if (dot(b,c) != 0)
            t2 = dot(a,c)/dot(b,c);
        if (t2 >= 0. && t2 <= 1.0)
            {
            if (firstIntFound)
                {
                Int2 = v2+t2*(v1-v2);
                }
            else
                {
                Int1 = v2+t2*(v1-v2);
                firstIntFound = true;
                };
            };

        v1=v2;
        };
    double maxSeparation = max(norm(p-Int1),norm(p-Int2));
    double2 newCellPos1 = initialCellPosition + separationFraction*maxSeparation*ray;
    double2 newCellPos2 = initialCellPosition - separationFraction*maxSeparation*ray;
    Box->putInBoxReal(newCellPos1);
    Box->putInBoxReal(newCellPos2);

    //This call updates many of the base data structres, but (among other things) does not actually
    //set the new cell position
    Simple2DActiveCell::cellDivision(parameters);
    {
    ArrayHandle<double2> cp(cellPositions);
    cp.data[cellIdx] = newCellPos1;
    cp.data[cdcellIdx] = newCellPos2;//(yw)cp.data[Ncells-1] = newCellPos2;
    }
    resizeAndReset();
    };


//////////////////////////
  
  /*!
This function supports cellDivisions, updating data structures in Simple2DActiveCell
This function will first call Simple2DCell's routine, and then
grows the cellDirectors and Motility arrays, and assign the new cell
(the last element of those arrays) the values of the cell given by parameters[0]
Note that dParams does nothing
 */
void Simple2DActiveCell::cellDivision(const vector<int> &parameters, const vector<double> &dParams)
    {
    //The Simple2DCell routine will increment Ncells by one, and then update other data structures
    Simple2DCell::cellDivision(parameters);
    int cellIdx = parameters[0];
    int cdcellIdx = parameters[1];//(yw)
    //growGPUArray(cellDirectors,1);
    //growGPUArray(Motility,1);
    noise.Reproducible = Reproducible;
        {//arrayhandle scope
        ArrayHandle<double2> h_mot(Motility); h_mot.data[cdcellIdx] = h_mot.data[cellIdx];
        ArrayHandle<double> h_cd(cellDirectors); h_cd.data[cdcellIdx] = noise.getRealUniform(0.,2*PI);
        ArrayHandle<double2> h_v(cellVelocities);
        h_v.data[cdcellIdx].x = h_mot.data[cdcellIdx].x*cos(h_cd.data[cdcellIdx]);
        h_v.data[cdcellIdx].y = h_mot.data[cdcellIdx].x*sin(h_cd.data[cdcellIdx]);
        };
    }; 


/////////////////////////////////////

/*!
This function supports cellDivisions, updating data structures in Simple2DCell
This function will grow the cell lists by 1 and assign the new cell
(the last element of those arrays) the values of the cell given by parameters[0]
Note that dParams does nothing by default, but allows more general virtual functions to be defined
downstream (used in the Voronoi branch)
 */
void Simple2DCell::cellDivision(const vector<int> &parameters, const vector<double> &dParams)
    {
    forcesUpToDate=false;
    //Ncells += 1;
    n_idx = Index2D(vertexMax,Ncells);
    int cellIdx = parameters[0];
    int cdcellIdx = parameters[1];

    //additions to the spatial sorting vectors...
    //itt.push_back(Ncells-1);
    //tti.push_back(Ncells-1);
    //tagToIdx.push_back(Ncells-1);
    //idxToTag.push_back(Ncells-1);

    //AreaPeri will have its values updated in a geometry routine... just change the length
    //AreaPeri.resize(Ncells);

    //use the copy and grow mechanism where we need to actually set values
    //growGPUArray(AreaPeriPreferences,1); //(nc)
    //growGPUArray(Moduli,1);
    //growGPUArray(cellMasses,1);
    //growGPUArray(cellVelocities,1);
    //growGPUArray(cellType,1);
    //growGPUArray(cellPositions,1);

        {//arrayhandle scope
        ArrayHandle<double2> h_APP(AreaPeriPreferences); h_APP.data[cdcellIdx] = h_APP.data[cellIdx];
        ArrayHandle<double2> h_Mod(Moduli); h_Mod.data[cdcellIdx] = h_Mod.data[cellIdx];
        ArrayHandle<int> h_ct(cellType); h_ct.data[cdcellIdx] = h_ct.data[cellIdx];
        ArrayHandle<double> h_cm(cellMasses);  h_cm.data[cdcellIdx] = h_cm.data[cellIdx];
        ArrayHandle<double2> h_v(cellVelocities); h_v.data[cdcellIdx] = make_double2(0.0,0.0);
        };
    };


