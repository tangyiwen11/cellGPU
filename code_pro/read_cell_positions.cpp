// read cell positions from txt

std::string line;
std::ifstream myFile("XY_npts_2430_Hex.txt");
//Ncells = newCellPositionx.size();
//if(cellPositions.getNumElements() != Ncells) cellPositions.resize(Ncells);
ArrayHandle<Dscalar2> h_p(cellPositions,access_location::host,access_mode::overwrite);
for (int ii = 0; ii < Ncells; ++ii)
    {
    getline(myFile, line);
    std::istringstream lineStream(line);
    double first, second;
    lineStream >> first >> second;
    printf("cx_%.4f_cy_%.4f\n", first,second);
    h_p.data[ii].x = first;
    h_p.data[ii].y = second;
    };




////////////////////////////

void Simple2DCell::setCellPositionsRandomly()
    {
    cellPositions.resize(Ncells);
    Dscalar boxsize = sqrt((Dscalar)Ncells);
    Box->setSquare(boxsize,boxsize);
    noise.Reproducible = Reproducible;

    std::string line;
    std::ifstream myFile("XY_npts_2430_Hex.txt");
    //Ncells = newCellPositionx.size();
    //if(cellPositions.getNumElements() != Ncells) cellPositions.resize(Ncells);

    ArrayHandle<Dscalar2> h_p(cellPositions,access_location::host,access_mode::overwrite);
    for (int ii = 0; ii < Ncells; ++ii)
        {
        getline(myFile, line);
        std::istringstream lineStream(line);
        double first, second;
        lineStream >> first >> second;
        Dscalar x = noise.getRealUniform(0.0,boxsize);
        Dscalar y = noise.getRealUniform(0.0,boxsize);
        printf("cx_%.4f_cy_%.4f\n", first,second);
        h_p.data[ii].x = first;
        h_p.data[ii].y = second;
        };
    };

///////////////////////////////


void Simple2DCell::setCellPositions(vector<Dscalar2> newCellPositions)
    {
    Ncells = newCellPositions.size();
    Dscalar boxsize = sqrt((Dscalar)Ncells);
    Box->setSquare(boxsize,boxsize);
    noise.Reproducible = Reproducible;
    if(cellPositions.getNumElements() != Ncells) cellPositions.resize(Ncells);
    ArrayHandle<Dscalar2> h_p(cellPositions,access_location::host,access_mode::overwrite);
    for (int ii = 0; ii < Ncells; ++ii)
        h_p.data[ii] = newCellPositions[ii];
    }




void Simple2DCell::initializeSimple2DCell(int n)
    {
    Ncells = n;
    Nvertices = 2*Ncells;

    //setting cell positions randomly also auto-generates a square box with L = sqrt(Ncells)
    //setCellPositionsRandomly();
    noise.Reproducible = Reproducible;
    vector<Dscalar2> positions(Ncells);
    Dscalar boxsize = sqrt((Dscalar)Ncells);
    Box->setSquare(boxsize,boxsize);
    noise.Reproducible = Reproducible;
    for (int ii = 0; ii < Ncells; ++ii)
        {
        Dscalar x = noise.getRealUniform(0.0,boxsize);
        Dscalar y = noise.getRealUniform(0.0,boxsize);
        positions[ii].x = x;
        positions[ii].y = y;
        };
    setCellPositions(positions);
    AreaPeri.resize(Ncells);
    cellForces.resize(Ncells);
    setCellPreferencesUniform(1.0,3.8);
    setModuliUniform(1.0,1.0);
    setCellTypeUniform(0);
    cellMasses.resize(Ncells);
    cellVelocities.resize(Ncells);
    vector<Dscalar> masses(Ncells,1.0);
    fillGPUArrayWithVector(masses,cellMasses);
    vector<Dscalar2> velocities(Ncells,make_Dscalar2(0.0,0.0));
    fillGPUArrayWithVector(velocities,cellVelocities);

    vertexForces.resize(Nvertices);
    };
