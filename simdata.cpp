#include "simdata.hpp"

#define SPECIAL_CELL = 999

SimData::SimData(string inputstream)
{  
  pStdElement = new StandardElements();      

  instream = inputstream;
  
  outstream = inputstream;
  stringstream streamword(inputstream);
  streamword.imbue(locale(locale(), new wordtokens()));
  istream_iterator<string> begin(streamword);
  istream_iterator<string> end;
  vector<string> vwords(begin, end);
  outstream = vwords[0];
  nNodes = 0;
  nGhostCells = 0;
  nCells = 0;  
  
  dNotNumber = -999.999;
  dNeimanForce = 1000.0;
  
  double SigmaH = 286.67;
  double SigmaV = 500.0;
  
  vPointPass.push_back(0);
  vPointPass.push_back(2);
  vPointPass.push_back(1);
  
  vPointCoord.resize(3, 0.0);
          
  vvPlate.resize(3);
  for (int i = 0; i < 3; i++) vvPlate[i].resize(3, 0.0);
      
  vsPhysicalBoundary.resize(11);
  // X
  vsPhysicalBoundary[0].ntype = 1;
  vsPhysicalBoundary[0].nmark = -1111111;  
  vsPhysicalBoundary[0].vCondition.push_back(0);
  vsPhysicalBoundary[0].vCondition.push_back(dNotNumber);
  vsPhysicalBoundary[0].vCondition.push_back(dNotNumber);
  
  vsPhysicalBoundary[1].ntype = 2;
  vsPhysicalBoundary[1].nmark = -1111112;  
  //vsPhysicalBoundary[1].vCondition.push_back(0.0);
  //vsPhysicalBoundary[1].vCondition.push_back(dNotNumber);
  //vsPhysicalBoundary[1].vCondition.push_back(dNotNumber);
  vsPhysicalBoundary[1].vCondition.push_back(300.0); // Bar
  
  // Y
  vsPhysicalBoundary[2].ntype = 2;
  vsPhysicalBoundary[2].nmark = -2222221;  
  //vsPhysicalBoundary[2].vCondition.push_back(dNotNumber);
  //vsPhysicalBoundary[2].vCondition.push_back(0.0);
  //vsPhysicalBoundary[2].vCondition.push_back(dNotNumber);
  vsPhysicalBoundary[2].vCondition.push_back(250.0); // Bar
  
  vsPhysicalBoundary[3].ntype = 2;
  vsPhysicalBoundary[3].nmark = -2222222;  
  //vsPhysicalBoundary[3].vCondition.push_back(dNotNumber);
  //vsPhysicalBoundary[3].vCondition.push_back(0.0);
  //vsPhysicalBoundary[3].vCondition.push_back(dNotNumber);
  vsPhysicalBoundary[3].vCondition.push_back(250.0); // Bar
  
  // Z
  vsPhysicalBoundary[4].ntype = 1;
  vsPhysicalBoundary[4].nmark = -3333331;  
  vsPhysicalBoundary[4].vCondition.push_back(dNotNumber); // Bar
  vsPhysicalBoundary[4].vCondition.push_back(dNotNumber); // Bar
  vsPhysicalBoundary[4].vCondition.push_back(0.0); // Bar
  
  vsPhysicalBoundary[5].ntype = 1;
  vsPhysicalBoundary[5].nmark = -3333332;
  vsPhysicalBoundary[5].vCondition.push_back(dNotNumber); // Bar
  vsPhysicalBoundary[5].vCondition.push_back(dNotNumber); // Bar
  vsPhysicalBoundary[5].vCondition.push_back(0.0); // Bar
  //vsPhysicalBoundary[5].vCondition.push_back(500.0); // Bar
  
  //wells 
  nWells = 1;
  vsWell.resize(nWells);
  
  // well 1
  vsWell[0].vWellCoordinate.clear();
  vsWell[0].vWellCoordinate.push_back(0.1); // x
  vsWell[0].vWellCoordinate.push_back(0.0); // y
  vsWell[0].vWellCoordinate.push_back(-100.0); // z0
  vsWell[0].vWellCoordinate.push_back(100.0); // z1  
  vsWell[0].Type = "WCONINJE";
  vsWell[0].radius_poisk = 1; // m

  // Kirill's renumbering
  pRenum = new renum(); 
}

SimData::~SimData()
{
}

void SimData::defineRockProperties()
{  
  vsCellRockProps.resize(nCells);
      
  for(int icell = 0; icell < nCells; icell++)
  {    
    vsCellRockProps[icell].zmf.assign(5,0.0);
    vsCellRockProps[icell].stress.assign(6,0.0);
    vsCellRockProps[icell].biot_plas = 0.0;
    
    
    
    if( vsCellCustom[icell].nMarker == 9999999 )
    {
      // Reservoir
      vsCellRockProps[icell].density = 2500.0;  // rock density kg/m3    
      vsCellRockProps[icell].heat_capacity = 2000; // rock heat capacity (same units as for liquids )
      vsCellRockProps[icell].thc = 150;            // rock thermal conductivity (same units as for liquids )
      vsCellRockProps[icell].thc_x = vsCellRockProps[icell].thc;
      vsCellRockProps[icell].thc_y = vsCellRockProps[icell].thc_x;
      vsCellRockProps[icell].thc_z = vsCellRockProps[icell].thc_x;

      vsCellRockProps[icell].poro = 0.3; // initial porosity
      vsCellRockProps[icell].perm = 25; // permeability mD
      
      vsCellRockProps[icell].perm_x = vsCellRockProps[icell].perm;
      vsCellRockProps[icell].perm_y = vsCellRockProps[icell].perm;
      vsCellRockProps[icell].perm_z = vsCellRockProps[icell].perm;
      
      vsCellRockProps[icell].biot = 1;      // Rock biot modulus (stress calsulation). Default 1.0.
      vsCellRockProps[icell].biot_flow = 0; // same as previous  (porosity calculation). Default = vsCellRockProps[icell].biot.
      
      vsCellRockProps[icell].young = 5.0;   // Young modulus *e10 Pa (default 1.0)
      vsCellRockProps[icell].poisson = 0.25; // Poisson ration (default 0.25 and < 0.45)
      
      vsCellRockProps[icell].cohesion = 1e10; // Plasticity cohesion *e10 Pa
      vsCellRockProps[icell].friction = 0;  // Plasticity friction angle < 90
      vsCellRockProps[icell].dilation = 0;  // Plasticity dilation angle < 90    
      
      vsCellRockProps[icell].temp = 343.0;       // initial temperature (same as the simulator input, K)
      
      // define pressure as a function of depth where 'fabs ( vsCellCustom[icell].vCenter[2] )' is depth
      vsCellRockProps[icell].pressure = 180; // initial pressure Bar
      
      //vsCellRockProps[icell].volmult = 1.0;            
      vsCellRockProps[icell].ref_pres = vsCellRockProps[icell].pressure; // reference pressure for porosity calculation      
      vsCellRockProps[icell].ref_temp = vsCellRockProps[icell].temp;     // reference temperature. If (ref_temp != temp) then we include initial thermal strain
      
      vsCellRockProps[icell].thermal_expansion = 2e-5; // rock thermal expansion (1/K)
      vsCellRockProps[icell].pore_thermal_expansion = 3.0 * vsCellRockProps[icell].thermal_expansion * (vsCellRockProps[icell].biot - vsCellRockProps[icell].poro); // pores thermal expansion (formula is exact). You may change values.
      vsCellRockProps[icell].poron = 0; // permability calculation k=k0*pow(phi/phi_o,poron)
      vsCellRockProps[icell].model = 0; // Model type 0 - Elasticity, 1 - Drucker Prager (see Manual)
      
      // initial composition - no water
      vsCellRockProps[icell].zmf[1] = 0.5;
      vsCellRockProps[icell].zmf[0] = 0.5;

      // general comments
      // you may add new variables, such as swat, pbub, and etc.
      // add output of these variables into femout.cpp
    }
      else
    {
      // Any other part
		  vsCellRockProps[icell].density = 2500.0;  // rock density kg/m3    
		  vsCellRockProps[icell].heat_capacity = 2000; // rock heat capacity (same units as for liquids )
		  vsCellRockProps[icell].thc = 150;            // rock thermal conductivity (same units as for liquids )
		  vsCellRockProps[icell].thc_x = vsCellRockProps[icell].thc;
		  vsCellRockProps[icell].thc_y = vsCellRockProps[icell].thc_x;
		  vsCellRockProps[icell].thc_z = vsCellRockProps[icell].thc_x;

		  vsCellRockProps[icell].poro = 0.3; // initial porosity
		  vsCellRockProps[icell].perm = 25; // permeability mD

		  vsCellRockProps[icell].perm_x = vsCellRockProps[icell].perm;
		  vsCellRockProps[icell].perm_y = vsCellRockProps[icell].perm;
		  vsCellRockProps[icell].perm_z = vsCellRockProps[icell].perm;

		  vsCellRockProps[icell].biot = 1;      // Rock biot modulus (stress calsulation). Default 1.0.
		  vsCellRockProps[icell].biot_flow = 0; // same as previous  (porosity calculation). Default = vsCellRockProps[icell].biot.

		  vsCellRockProps[icell].young = 5.0;   // Young modulus *e10 Pa (default 1.0)
		  vsCellRockProps[icell].poisson = 0.25; // Poisson ration (default 0.25 and < 0.45)

		  vsCellRockProps[icell].cohesion = 1e10; // Plasticity cohesion *e10 Pa
		  vsCellRockProps[icell].friction = 0;  // Plasticity friction angle < 90
		  vsCellRockProps[icell].dilation = 0;  // Plasticity dilation angle < 90    

		  vsCellRockProps[icell].temp = 343.0;       // initial temperature (same as the simulator input, K)

		  // define pressure as a function of depth where 'fabs ( vsCellCustom[icell].vCenter[2] )' is depth
		  vsCellRockProps[icell].pressure = 180; // initial pressure Bar

		  //vsCellRockProps[icell].volmult = 1.0;
		  vsCellRockProps[icell].ref_pres = vsCellRockProps[icell].pressure; // reference pressure for porosity calculation      
		  vsCellRockProps[icell].ref_temp = vsCellRockProps[icell].temp;     // reference temperature. If (ref_temp != temp) then we include initial thermal strain

		  vsCellRockProps[icell].thermal_expansion = 2e-5; // rock thermal expansion (1/K)
		  vsCellRockProps[icell].pore_thermal_expansion = 3.0 * vsCellRockProps[icell].thermal_expansion * (vsCellRockProps[icell].biot - vsCellRockProps[icell].poro); // pores thermal expansion (formula is exact). You may change values.
		  vsCellRockProps[icell].poron = 0; // permability calculation k=k0*pow(phi/phi_o,poron)
		  vsCellRockProps[icell].model = 0; // Model type 0 - Elasticity, 1 - Drucker Prager (see Manual)

		  // initial composition - no water
		  vsCellRockProps[icell].zmf[1] = 0.5;
		  vsCellRockProps[icell].zmf[0] = 0.5;
    }

      // intial stress
      double sz = 0;      
	  double sy = 0;
	  double sx = 0;
      double sxy = 0;      
      vsCellRockProps[icell].stress[0] = -sx; 
      vsCellRockProps[icell].stress[1] = -sy; 
      vsCellRockProps[icell].stress[2] = -sz; 
      vsCellRockProps[icell].stress[5] = sxy;
            
      }
}

void SimData::defineBoundaryAquifersEmil()
{

	// clean up all porosity values
	for (int ic = 0; ic < nCells; ic++)
	{
		vsCellRockProps[ic].volmult = 1.0;
	}

	// loop over all cells and assign the aquifers on the radius
	double radius_aquifer = 48;
	for (int icell = 0; icell < nCells; icell++)
	{
		double distance = sqrt(vsCellCustom[icell].vCenter[0] * vsCellCustom[icell].vCenter[0] +
			vsCellCustom[icell].vCenter[1] * vsCellCustom[icell].vCenter[1]);
		if (distance >= radius_aquifer) vsCellRockProps[icell].volmult = 1e5;
	}

	//// loop over all faces and find boundary faces
	//for (int iface = 0; iface < nFaces; iface++)
	//{
	//		// criterion 1 (right boundary -1111112)
	//		// X or Y are dominant then acitvate an aquifer
	//		if (vsFaceCustom[iface].nMarker == -1111112)
	//		{
	//			for (int ic = 0; ic < vsFaceCustom[iface].nNeighbors; ++ic)
	//			{
	//				int icell = vsFaceCustom[iface].vNeighbors[ic];
	//				vsCellRockProps[icell].volmult = 1e5;
	//			}
	//		}

	//		// criterion 2 (upper boundary -2222222)
	//		// X or Y are dominant then acitvate an aquifer
	//		if (vsFaceCustom[iface].nMarker == -2222222)
	//		{
	//			for (int ic = 0; ic < vsFaceCustom[iface].nNeighbors; ++ic)
	//			{
	//				int icell = vsFaceCustom[iface].vNeighbors[ic];
	//				vsCellRockProps[icell].volmult = 1e5;
	//			}
	//		}

	//		// criterion 3 (bottom boundary -2222221)
	//		// X or Y are dominant then acitvate an aquifer
	//		if (vsFaceCustom[iface].nMarker == -2222221)
	//		{
	//			for (int ic = 0; ic < vsFaceCustom[iface].nNeighbors; ++ic)
	//			{
	//				int icell = vsFaceCustom[iface].vNeighbors[ic];
	//				vsCellRockProps[icell].volmult = 1e5;
	//			}
	//		}

	//}
}

void SimData::createDoublePorosityModel()
{
  vCellsHasGhost.assign(nCells,false);
  return;
  
  for(int icell = 0; icell < nCells; icell++)
  {
    if(vsCellCustom[icell].nMarker == 9999992)
    {
      vCellsHasGhost[icell] = true;
    }
  }
  
  for ( int icell = 0; icell < nCells; icell++ )
  {
    if ( vsCellCustom[icell].nMarker == 9999991 )
    {
      for ( int ic = 0; ic < vsCellCustom[icell].nNeighbors; ++ic )
      {	
        if ( vsCellCustom[vsCellCustom[icell].vNeighbors[ic]].nMarker == 9999992)
          vCellsHasGhost[icell] = true;
      }
    }
  }
  
  for(int icell = 0; icell < nCells; icell++)
  {
    if(vsCellCustom[icell].nMarker > 9999992)
      vCellsHasGhost[icell] = false;
  }

  nGhostCells = 0;
  for ( int icell = 0; icell < nCells; icell++ )
  {
    if( vCellsHasGhost[icell] )
      nGhostCells++;
  }  
  
  dGhostAperture = 5e-3;    // m
  dJointDX = 30;            // m
  dJointDY = 60;            // m
  dGhostConductivity = 10;  // mD.m
  
  dGhostFraction = (dGhostAperture * dJointDX + dGhostAperture * dJointDY) / (dJointDX * dJointDY); // m3/m3
  dGhostPermeability = dGhostConductivity / dGhostAperture;                                         // mD
  dGhost2MasterTransGeomPart = 2.0 * dJointDX / (dJointDY / 4) + 2.0 * dJointDY / (dJointDX / 4);   // m/m
}

void SimData::readGmshFile()
{
  double scale = 1e0;
  Gelement Element3D, Element2D;
  
  fstream inputfile;
  inputfile.open (instream.c_str(), fstream::in);
  
  std::string inputline;  
       
  getline (inputfile, inputline); 
  
  std::stringstream streamline(stringstream::in | stringstream::out);
  streamline << inputline;  
  streamline.imbue(std::locale(std::locale(), new tokens()));
  std::istream_iterator<std::string> begin;
  std::istream_iterator<std::string> end;
  std::vector<std::string> vstrings;

  copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );
  
    while(vstrings[0] != "$Nodes")
    {
      getline(inputfile, inputline);
      cout << inputline << endl;
      streamline.clear(); vstrings.clear();
      streamline << inputline;
      streamline.imbue(std::locale(std::locale(), new tokens()));  
      copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );
    }
    getline(inputfile, inputline);
    streamline.clear(); vstrings.clear();
    streamline << inputline;
    streamline.imbue(std::locale(std::locale(), new tokens()));  
    copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );

    nNodes = atoi(vstrings[0].c_str());
    
    vvVrtxCoords.resize(nNodes);
    for(int i = 0; i < nNodes; i++ )
    {
      getline(inputfile, inputline);
      streamline.clear(); vstrings.clear();
      streamline << inputline;
      streamline.imbue(std::locale(std::locale(), new tokens()));  
      copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );

      for (int j = 1; j < 4; j++)
      {
        vvVrtxCoords[i].push_back( atof(vstrings[j].c_str()) / scale);
      }
    }
    nNodes = vvVrtxCoords.size();
    
    maxVrtxCoordsX = -1e10;
    maxVrtxCoordsY = -1e10;
    maxVrtxCoordsZ = -1e10;
    minVrtxCoordsX = 1e10;
    minVrtxCoordsY = 1e10;
    minVrtxCoordsZ = 1e10;
    
    for(int i = 0; i < nNodes; i++ )
    {
      if(vvVrtxCoords[i][0] > maxVrtxCoordsX ) maxVrtxCoordsX = vvVrtxCoords[i][0];
      if(vvVrtxCoords[i][1] > maxVrtxCoordsY ) maxVrtxCoordsY = vvVrtxCoords[i][1];
      if(vvVrtxCoords[i][2] > maxVrtxCoordsZ ) maxVrtxCoordsZ = vvVrtxCoords[i][2];
      
      if(vvVrtxCoords[i][0] < minVrtxCoordsX ) minVrtxCoordsX = vvVrtxCoords[i][0];
      if(vvVrtxCoords[i][1] < minVrtxCoordsY ) minVrtxCoordsY = vvVrtxCoords[i][1];
      if(vvVrtxCoords[i][2] < minVrtxCoordsZ ) minVrtxCoordsZ = vvVrtxCoords[i][2];
    }
cout << "deb 2" << endl;    
    
    while(vstrings[0] != "$Elements")
    {
     getline(inputfile, inputline);
     streamline.clear(); vstrings.clear();
     streamline << inputline;
     streamline.imbue(std::locale(std::locale(), new tokens()));  
     copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );
    }
cout << "deb 3" << endl;    
      
    getline(inputfile, inputline);
    streamline.clear(); vstrings.clear();
    streamline << inputline;
    streamline.imbue(std::locale(std::locale(), new tokens()));  
    copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );
        
    int n = atoi(vstrings[0].c_str());
	
	vector<int>face1(4);
	vector<int>face2(4);
	vector<int>faceX1(8);
	vector<int>faceX2(8);

	set<vector<int> >VerticesSorted;
	std::pair<std::set<vector<int>>::iterator, bool> ret;
        
    nCells = 0; nFaces = 0;
    for(int i = 0; i < n; i++ )
    {
      getline(inputfile, inputline);
      streamline.clear(); vstrings.clear();
      streamline << inputline;
      streamline.imbue(std::locale(std::locale(), new tokens()));  
      copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );

      int elem_type = atoi(vstrings[1].c_str());
      int ntags = atoi(vstrings[2].c_str());
      int node_list_run = 0;
      int node_list_end = 0;
      
      Element3D.nNeighbors = 0;
      Element3D.nVertices = 0;
      Element3D.vCenter.clear();
      Element3D.vNeighbors.clear();
      Element3D.vNormal.clear();
      Element3D.vVertices.clear();
      Element3D.vVerticesNewnum.clear();
      Element3D.vVerticesSorted.clear();
                  
      Element2D.nNeighbors = 0;
      Element2D.nVertices = 0;
      Element2D.vCenter.clear();
      Element2D.vNeighbors.clear();
      Element2D.vNormal.clear();
      Element2D.vVertices.clear();
      Element2D.vVerticesNewnum.clear();
      Element2D.vVerticesSorted.clear();
   
switch (elem_type)
{
case 2:
	node_list_run = 3 + ntags;
	node_list_end = node_list_run + 3;
	for (int j = node_list_run; j < node_list_end; j++) Element2D.vVertices.push_back(atoi(vstrings[j].c_str()));
	Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
	Element2D.vtkIndex = 5;
	Element2D.formIndex = TRGLE3;
	Element2D.nMarker = atoi(vstrings[3].c_str());
	Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
	vsFaceCustom.push_back(Element2D);
	sort(Element2D.vVertices.begin(), Element2D.vVertices.end());
	VerticesSorted.insert(Element2D.vVertices);
	nFaces++;
	break;

case 9:
	node_list_run = 3 + ntags;
	node_list_end = node_list_run + 6;
	for (int j = node_list_run; j < node_list_end; j++) Element2D.vVertices.push_back(atoi(vstrings[j].c_str()));            
	Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());	
	Element2D.vtkIndex = 22;
	Element2D.formIndex = TRGLE6;
	Element2D.nMarker = atoi(vstrings[3].c_str());            
	Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
	vsFaceCustom.push_back(Element2D);
	sort(Element2D.vVertices.begin(), Element2D.vVertices.end());
	VerticesSorted.insert(Element2D.vVertices);
	nFaces++;
	break;

case 3:
	node_list_run = 3 + ntags;
	node_list_end = node_list_run + 4;
	for (int j = node_list_run; j < node_list_end; j++) Element2D.vVertices.push_back(atoi(vstrings[j].c_str()));
	Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
	Element2D.vtkIndex = 9;
	Element2D.formIndex = QUAD4;
	Element2D.nMarker = atoi(vstrings[3].c_str());
	Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
	vsFaceCustom.push_back(Element2D);
	sort(Element2D.vVertices.begin(), Element2D.vVertices.end());
	VerticesSorted.insert(Element2D.vVertices);
	nFaces++;
	break;

case 16:
	node_list_run = 3 + ntags;
	node_list_end = node_list_run + 8;
	for (int j = node_list_run; j < node_list_end; j++) Element2D.vVertices.push_back(atoi(vstrings[j].c_str()));  
	Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
	Element2D.vtkIndex = 23;
	Element2D.formIndex = QUAD8;
	Element2D.nMarker = atoi(vstrings[3].c_str());
	Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
	vsFaceCustom.push_back( Element2D );
	sort(Element2D.vVertices.begin(), Element2D.vVertices.end());
	VerticesSorted.insert(Element2D.vVertices);
	nFaces++;
	break;

case 4:
	node_list_run = 3 + ntags;
	node_list_end = node_list_run + 4;
	for (int j = node_list_run; j < node_list_end; j++) Element3D.vVertices.push_back(atoi(vstrings[j].c_str()));
	Element3D.nMarkerGMSH = atoi(vstrings[4].c_str());
	Element3D.vtkIndex = 10;
	Element3D.formIndex = TETRA4;
	Element3D.nMarker = atoi(vstrings[3].c_str());
	vsCellCustom.push_back(Element3D);
	nCells++;
	break;

case 5:
	node_list_run = 3 + ntags;
	node_list_end = node_list_run + 8;
	for (int j = node_list_run; j < node_list_end; j++) Element3D.vVertices.push_back(atoi(vstrings[j].c_str()));
	Element3D.nMarkerGMSH = atoi(vstrings[4].c_str());
	Element3D.vtkIndex = 12;
	Element3D.formIndex = PRISM8;
	Element3D.nMarker = atoi(vstrings[3].c_str());
	vsCellCustom.push_back(Element3D);
	nCells++;

	if (Element3D.nMarker > 0 && Element3D.nMarker < 1111110)
	{
		face1[0] = Element3D.vVertices[0];
		face1[1] = Element3D.vVertices[1];
		face1[2] = Element3D.vVertices[5];
		face1[3] = Element3D.vVertices[4];

		face2 = face1;

		sort(face1.begin(), face1.end());
		ret = VerticesSorted.insert(face1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = face2;
			Element2D.vtkIndex = 9;
			Element2D.formIndex = QUAD4;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}

		face1[0] = Element3D.vVertices[1];
		face1[1] = Element3D.vVertices[2];
		face1[2] = Element3D.vVertices[6];
		face1[3] = Element3D.vVertices[5];

		face2 = face1;

		sort(face1.begin(), face1.end());
		ret = VerticesSorted.insert(face1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = face2;
			Element2D.vtkIndex = 9;
			Element2D.formIndex = QUAD4;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}

		face1[0] = Element3D.vVertices[2];
		face1[1] = Element3D.vVertices[3];
		face1[2] = Element3D.vVertices[7];
		face1[3] = Element3D.vVertices[6];

		face2 = face1;

		sort(face1.begin(), face1.end());
		ret = VerticesSorted.insert(face1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = face2;
			Element2D.vtkIndex = 9;
			Element2D.formIndex = QUAD4;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}

		face1[0] = Element3D.vVertices[3];
		face1[1] = Element3D.vVertices[0];
		face1[2] = Element3D.vVertices[4];
		face1[3] = Element3D.vVertices[7];

		face2 = face1;

		sort(face1.begin(), face1.end());
		ret = VerticesSorted.insert(face1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = face2;
			Element2D.vtkIndex = 9;
			Element2D.formIndex = QUAD4;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}
	}
	break;

case 6:
	node_list_run = 3 + ntags;
	node_list_end = node_list_run + 6;
	for (int j = node_list_run; j < node_list_end; j++) Element3D.vVertices.push_back(atoi(vstrings[j].c_str()));
	Element3D.nMarkerGMSH = atoi(vstrings[4].c_str());
	Element3D.vtkIndex = 13;
	Element3D.formIndex = PRISM6;
	Element3D.nMarker = atoi(vstrings[3].c_str());
	vsCellCustom.push_back(Element3D);
	nCells++;

	if (Element3D.nMarker > 0 && Element3D.nMarker < 1111110)
	{

		face1[0] = Element3D.vVertices[0];
		face1[1] = Element3D.vVertices[1];
		face1[2] = Element3D.vVertices[4];
		face1[3] = Element3D.vVertices[3];

		face2 = face1;

		sort(face1.begin(), face1.end());
		ret = VerticesSorted.insert(face1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = face2;
			Element2D.vtkIndex = 9;
			Element2D.formIndex = QUAD4;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}

		face1[0] = Element3D.vVertices[1];
		face1[1] = Element3D.vVertices[2];
		face1[2] = Element3D.vVertices[5];
		face1[3] = Element3D.vVertices[4];

		face2 = face1;

		sort(face1.begin(), face1.end());
		ret = VerticesSorted.insert(face1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = face2;
			Element2D.vtkIndex = 9;
			Element2D.formIndex = QUAD4;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}

		face1[0] = Element3D.vVertices[0];
		face1[1] = Element3D.vVertices[2];
		face1[2] = Element3D.vVertices[5];
		face1[3] = Element3D.vVertices[3];

		face2 = face1;

		sort(face1.begin(), face1.end());
		ret = VerticesSorted.insert(face1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = face2;
			Element2D.vtkIndex = 9;
			Element2D.formIndex = QUAD4;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}
	}
	break;


case 18:
	node_list_run = 3 + ntags;
	node_list_end = node_list_run + 15;
	for (int j = node_list_run; j < node_list_end; j++) Element3D.vVertices.push_back(atoi(vstrings[j].c_str()));  
	Element3D.nMarkerGMSH = atoi(vstrings[4].c_str());	
	Element3D.vtkIndex = 26;
	Element3D.formIndex = PRISM15;
	Element3D.nMarker = atoi(vstrings[3].c_str());
	vsCellCustom.push_back( Element3D );
	nCells++;
if (Element3D.nMarker > 0 && Element3D.nMarker < 1111110)
	{

		faceX1[0] = Element3D.vVertices[0];
		faceX1[1] = Element3D.vVertices[3];
		faceX1[2] = Element3D.vVertices[5];
		faceX1[3] = Element3D.vVertices[2];
		faceX1[4] = Element3D.vVertices[8];
		faceX1[5] = Element3D.vVertices[13];
		faceX1[6] = Element3D.vVertices[11];
		faceX1[7] = Element3D.vVertices[7];

		faceX2 = faceX1;

		sort(faceX1.begin(), faceX1.end());
		ret = VerticesSorted.insert(faceX1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = faceX2;
			Element2D.vtkIndex = 23;
			Element2D.formIndex = QUAD8;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}

		faceX1[0] = Element3D.vVertices[1];
		faceX1[1] = Element3D.vVertices[4];
		faceX1[2] = Element3D.vVertices[3];
		faceX1[3] = Element3D.vVertices[0];
		faceX1[4] = Element3D.vVertices[10];
		faceX1[5] = Element3D.vVertices[12];
		faceX1[6] = Element3D.vVertices[8];
		faceX1[7] = Element3D.vVertices[6];

		faceX2 = faceX1;

		sort(faceX1.begin(), faceX1.end());
		ret = VerticesSorted.insert(faceX1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = faceX2;
			Element2D.vtkIndex = 23;
			Element2D.formIndex = QUAD8;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}

		faceX1[0] = Element3D.vVertices[1];
		faceX1[1] = Element3D.vVertices[2];
		faceX1[2] = Element3D.vVertices[5];
		faceX1[3] = Element3D.vVertices[4];
		faceX1[4] = Element3D.vVertices[9];
		faceX1[5] = Element3D.vVertices[11];
		faceX1[6] = Element3D.vVertices[14];
		faceX1[7] = Element3D.vVertices[10];

		faceX2 = faceX1;

		sort(faceX1.begin(), faceX1.end());
		ret = VerticesSorted.insert(faceX1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = faceX2;
			Element2D.vtkIndex = 23;
			Element2D.formIndex = QUAD8;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}

}
	break;

case 11:
	node_list_run = 3 + ntags;
	node_list_end = node_list_run + 10;
	for (int j = node_list_run; j < node_list_end; j++) Element3D.vVertices.push_back(atoi(vstrings[j].c_str()));

	// because different numbering (vtk vs gmsh)
	std::swap(Element3D.vVertices[8], Element3D.vVertices[9]);
	Element3D.nMarkerGMSH = atoi(vstrings[4].c_str());
	Element3D.vtkIndex = 24;
	Element3D.formIndex = TETRA10;
	Element3D.nMarker = atoi(vstrings[3].c_str());
	vsCellCustom.push_back(Element3D);
	nCells++;
	break;


case 17:
	node_list_run = 3 + ntags;
	node_list_end = node_list_run + 20;
	for (int j = node_list_run; j < node_list_end; j++) Element3D.vVertices.push_back(atoi(vstrings[j].c_str()));
	Element3D.nMarkerGMSH = atoi(vstrings[4].c_str());
	Element3D.vtkIndex = 25;
	Element3D.formIndex = PRISM20;
	Element3D.nMarker = atoi(vstrings[3].c_str());
	vsCellCustom.push_back(Element3D);
	nCells++;

	if (Element3D.nMarker > 0 && Element3D.nMarker < 1111110)
	{
	//1	
	  faceX1[0] = Element3D.vVertices[0];
		faceX1[1] = Element3D.vVertices[4];
		faceX1[2] = Element3D.vVertices[7];
		faceX1[3] = Element3D.vVertices[3];
		faceX1[4] = Element3D.vVertices[10];
		faceX1[5] = Element3D.vVertices[17];
		faceX1[6] = Element3D.vVertices[15];
		faceX1[7] = Element3D.vVertices[9];

		faceX2 = faceX1;

		sort(faceX1.begin(), faceX1.end());
		ret = VerticesSorted.insert(faceX1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = faceX2;
			Element2D.vtkIndex = 23;
			Element2D.formIndex = QUAD8;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}

	//2
	faceX1[0] = Element3D.vVertices[4];
		faceX1[1] = Element3D.vVertices[5];
		faceX1[2] = Element3D.vVertices[6];
		faceX1[3] = Element3D.vVertices[7];
		faceX1[4] = Element3D.vVertices[16];
		faceX1[5] = Element3D.vVertices[18];
		faceX1[6] = Element3D.vVertices[19];
		faceX1[7] = Element3D.vVertices[17];

		faceX2 = faceX1;

		sort(faceX1.begin(), faceX1.end());
		ret = VerticesSorted.insert(faceX1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = faceX2;
			Element2D.vtkIndex = 23;
			Element2D.formIndex = QUAD8;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}

	//3
		faceX1[0] = Element3D.vVertices[5];
		faceX1[1] = Element3D.vVertices[1];
		faceX1[2] = Element3D.vVertices[2];
		faceX1[3] = Element3D.vVertices[6];
		faceX1[4] = Element3D.vVertices[12];
		faceX1[5] = Element3D.vVertices[11];
		faceX1[6] = Element3D.vVertices[14];
		faceX1[7] = Element3D.vVertices[18];

		faceX2 = faceX1;

		sort(faceX1.begin(), faceX1.end());
		ret = VerticesSorted.insert(faceX1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = faceX2;
			Element2D.vtkIndex = 23;
			Element2D.formIndex = QUAD8;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}

	//4	
		faceX1[0] = Element3D.vVertices[1];
		faceX1[1] = Element3D.vVertices[0];
		faceX1[2] = Element3D.vVertices[3];
		faceX1[3] = Element3D.vVertices[2];
		faceX1[4] = Element3D.vVertices[8];
		faceX1[5] = Element3D.vVertices[9];
		faceX1[6] = Element3D.vVertices[13];
		faceX1[7] = Element3D.vVertices[11];

		faceX2 = faceX1;

		sort(faceX1.begin(), faceX1.end());
		ret = VerticesSorted.insert(faceX1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = faceX2;
			Element2D.vtkIndex = 23;
			Element2D.formIndex = QUAD8;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}
	//5	
		faceX1[0] = Element3D.vVertices[7];
		faceX1[1] = Element3D.vVertices[6];
		faceX1[2] = Element3D.vVertices[2];
		faceX1[3] = Element3D.vVertices[3];
		faceX1[4] = Element3D.vVertices[19];
		faceX1[5] = Element3D.vVertices[14];
		faceX1[6] = Element3D.vVertices[13];
		faceX1[7] = Element3D.vVertices[15];

		faceX2 = faceX1;

		sort(faceX1.begin(), faceX1.end());
		ret = VerticesSorted.insert(faceX1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = faceX2;
			Element2D.vtkIndex = 23;
			Element2D.formIndex = QUAD8;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}
	//6	
		faceX1[0] = Element3D.vVertices[1];
		faceX1[1] = Element3D.vVertices[0];
		faceX1[2] = Element3D.vVertices[4];
		faceX1[3] = Element3D.vVertices[5];
		faceX1[4] = Element3D.vVertices[8];
		faceX1[5] = Element3D.vVertices[10];
		faceX1[6] = Element3D.vVertices[16];
		faceX1[7] = Element3D.vVertices[12];

		faceX2 = faceX1;

		sort(faceX1.begin(), faceX1.end());
		ret = VerticesSorted.insert(faceX1);
		if (ret.second)
		{
			Element2D.nMarkerGMSH = atoi(vstrings[4].c_str());
			Element2D.vVertices = faceX2;
			Element2D.vtkIndex = 23;
			Element2D.formIndex = QUAD8;
			Element2D.nMarker = atoi(vstrings[3].c_str());
			Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
			vsFaceCustom.push_back(Element2D);
			nFaces++;
		}
	}
	break;

default:
	cout << "Element type " << elem_type << endl;
	cout << "Wrong element type. Supported: {2, 3, 4, 5, 6, 9, 11, 16, 17, 18}\n";
	exit(-1);
	break;
}
    }
    // c++ rule
    for(int icell = 0; icell < nCells; icell++)
    {
      vsCellCustom[icell].nVertices = vsCellCustom[icell].vVertices.size();
      for(int ivrtx = 0; ivrtx < vsCellCustom[icell].nVertices; ivrtx++ )
      {
        vsCellCustom[icell].vVertices[ivrtx] = vsCellCustom[icell].vVertices[ivrtx] - 1;
      }
    }
    
    for(int iface = 0; iface < nFaces; iface++)
    {
      vsFaceCustom[iface].nVertices = vsFaceCustom[iface].vVertices.size();
      for(int ivrtx = 0; ivrtx < vsFaceCustom[iface].nVertices; ivrtx++ )
      {
        vsFaceCustom[iface].vVertices[ivrtx] = vsFaceCustom[iface].vVertices[ivrtx] - 1;
      }
    }
    cout << "GMSH vrtxs " << nNodes << endl;
    cout << "GMSH cells " << nCells << endl;
    cout << "GMSH faces " << nFaces << endl;
}

void SimData::extractInternalFaces()
{
  vector<int> vLocalPolygonVertices;
  set<int> setLocalPolygonVertices;
  stringstream vertices_stream;
  
  set<string>  setIdenticalPolygons;
  pair<set<string>::iterator, bool> pair_itstring_bool;
  
  // At first we should write all input polygons
  for(int iface = 0; iface < nFaces; iface++ )    
  {
    setLocalPolygonVertices.clear();
    for(int ivrtx = 0; ivrtx < vsFaceCustom[iface].nVertices; ivrtx++ )
    {
      setLocalPolygonVertices.insert ( vsFaceCustom[iface].vVertices[ivrtx] );
    }
    set<int>::iterator it_set;
    vertices_stream.str ( "" ); 
    for(it_set = setLocalPolygonVertices.begin(); it_set != setLocalPolygonVertices.end(); ++it_set) vertices_stream << *it_set;
    //try write into set
    pair_itstring_bool = setIdenticalPolygons.insert ( vertices_stream.str() );    
  }

  // we dont know how many polygons
  for ( int icell = 0; icell < nCells; icell++ )
  {
    int job_percent = int ( ( 100. * icell ) / ( nCells ) );
    cout << "\r    " << job_percent << "%";

    // loop by all element faces
    for ( int iface = 0; iface < pStdElement->elementProps[ vsCellCustom[icell].formIndex ].facesInElement; iface++ )
    {
      vLocalPolygonVertices.clear();
      setLocalPolygonVertices.clear();
      //num of nodes on face
      int nnodes_ = pStdElement->elementProps[ vsCellCustom[icell].formIndex].vvFacesNodes[iface].size();

      //write nodes into vector
      for ( int inode = 0; inode < nnodes_; inode++ )
      {
        int local_node = pStdElement->elementProps[ vsCellCustom[icell].formIndex].vvFacesNodes[iface][inode];
        int global_node = vsCellCustom[icell].vVertices[local_node];

        vLocalPolygonVertices.push_back ( global_node );
        setLocalPolygonVertices.insert ( global_node );
      }

      set<int>::iterator it_set;
      vertices_stream.str ( "" );

      for ( it_set = setLocalPolygonVertices.begin(); it_set != setLocalPolygonVertices.end(); ++it_set )
        vertices_stream << *it_set;

      //try write into set
      pair_itstring_bool = setIdenticalPolygons.insert ( vertices_stream.str() );

      //element will be written if we see it first time
      Gelement temporaryElement;
      if ( pair_itstring_bool.second == true )
      {
        // extend vFaceCustom
        temporaryElement.fluidElement = -1;
        temporaryElement.nNeighbors = 0;
        temporaryElement.nVertices = vLocalPolygonVertices.size();
        temporaryElement.vVertices.resize ( temporaryElement.nVertices );
        temporaryElement.vVertices = vLocalPolygonVertices;

        if ( temporaryElement.nVertices == 3 )
        {
          temporaryElement.vtkIndex = 5;
          temporaryElement.formIndex = TRGLE3;
        }

        if ( temporaryElement.nVertices == 6 )
        {
          temporaryElement.vtkIndex = 22;
          temporaryElement.formIndex = TRGLE6;
	  int j = 0;
	  for(unsigned int i = 0; i < temporaryElement.nVertices; i=i+2, j++)
	    temporaryElement.vVertices[j] = vLocalPolygonVertices[i];
	  for(unsigned int i = 1; i < temporaryElement.nVertices; i=i+2, j++)
	    temporaryElement.vVertices[j] = vLocalPolygonVertices[i];
        }

        if ( temporaryElement.nVertices == 4 )
        {
          temporaryElement.vtkIndex = 9;
          temporaryElement.formIndex = QUAD4;
        }

        if ( temporaryElement.nVertices == 8 )
        {
          temporaryElement.vtkIndex = 23;
          temporaryElement.formIndex = QUAD8;
	  int j = 0;
	  for(unsigned int i = 0; i < temporaryElement.nVertices; i=i+2, j++)
	    temporaryElement.vVertices[j] = vLocalPolygonVertices[i];
	  for(unsigned int i = 1; i < temporaryElement.nVertices; i=i+2, j++)
	    temporaryElement.vVertices[j] = vLocalPolygonVertices[i];
        }
        temporaryElement.nMarker = 0;
        vsFaceCustom.push_back ( temporaryElement );            
      }
    }
  }
  
  nFaces = vsFaceCustom.size();  
  
}

void SimData::convertGmsh2Sim()
{
  int counter_ = 0;  
  for(int iface = 0; iface < nFaces; iface++)
  {    
    methodElementCenter(iface, vsFaceCustom);    
    if(vsFaceCustom[iface].nMarker != 0) methodFaceNormalVector(iface, vsFaceCustom);              
  }
  methodChangeFacesNormalVector();
  
  counter_ = 0;
  for(int icell = 0; icell < nCells; icell++)
  {
    methodElementCenter(icell, vsCellCustom);    
  }
      
  cout << "\t sort verticies" << endl;
  for(int iface = 0; iface < nFaces; iface++)
  {    
    vsFaceCustom[iface].vVerticesSorted.resize( vsFaceCustom[iface].nVertices );
    vsFaceCustom[iface].vVerticesSorted = vsFaceCustom[iface].vVertices;
    sort(vsFaceCustom[iface].vVerticesSorted.begin(), vsFaceCustom[iface].vVerticesSorted.end());
  }
  
  for(int icell = 0; icell < nCells; icell++)
  {    
    vsCellCustom[icell].vVerticesSorted.resize( vsCellCustom[icell].nVertices );
    vsCellCustom[icell].vVerticesSorted = vsCellCustom[icell].vVertices;
    sort(vsCellCustom[icell].vVerticesSorted.begin(), vsCellCustom[icell].vVerticesSorted.end());
  }  
  
  cout << "\t find face neighbor cells / cell neighbor faces (slow)" << endl;
  vsetPolyhedronPolygon.resize(nCells);
  vsetPolygonPolyhedron.resize(nFaces);
  for ( int iface = 0; iface < nFaces; iface++ )
  {
    int job_percent = int ( ( 100. * iface ) / ( nFaces ) );
    cout << "\r    " << job_percent << "%";      
    for ( int icell = 0; icell < nCells; icell++ )
    {
      if ( includes ( vsCellCustom[icell].vVerticesSorted.begin(), vsCellCustom[icell].vVerticesSorted.end(),
                      vsFaceCustom[iface].vVerticesSorted.begin(), vsFaceCustom[iface].vVerticesSorted.end() ) )
      {
        vsetPolyhedronPolygon[icell].insert ( iface );
        vsetPolygonPolyhedron[iface].insert ( icell );
      }
    }
  }
  
  for ( int iface = 0; iface < nFaces; iface++ )
  {
    if(vsetPolygonPolyhedron[iface].size() == 0)
    {
      cout << endl << "Polygon " << iface << " (" << nFaces << ") has no connected polyhedrons!" << endl;
      exit(0);      
    }
  }

  for ( int ic = 0; ic < nCells; ic++ )
  {
    if(vsetPolyhedronPolygon[ic].size() == 0)
    {
      cout << endl << "Polyhedron " << ic << " has no connected polygons!" << endl;
      exit(0);      
    }
  }
  
  cout << "\t find cell neighbor cells " << endl;
  vector<set<int> > vsetCellNeighborCells;
  vsetCellNeighborCells.resize(nCells);
  for ( int icell = 0; icell < nCells; icell++ )
  {
    set<int>::iterator it_set, it_set_f;

    for ( it_set = vsetPolyhedronPolygon[icell].begin(); it_set != vsetPolyhedronPolygon[icell].end(); ++it_set )
    {
      int nface = *it_set;

      for ( it_set_f = vsetPolygonPolyhedron[nface].begin(); it_set_f != vsetPolygonPolyhedron[nface].end(); ++it_set_f )
      {
        int ncell = *it_set_f;
        if(ncell != icell) vsetCellNeighborCells[icell].insert(ncell);
      }
    }

  }
  
  for ( int icell = 0; icell < nCells; icell++ )
  {
    vsCellCustom[icell].vNeighbors.clear();
    set<int>::iterator it_set;
    for ( it_set = vsetCellNeighborCells[icell].begin(); it_set != vsetCellNeighborCells[icell].end(); ++it_set )
    {
      vsCellCustom[icell].vNeighbors.push_back( *it_set );
    }

    vsCellCustom[icell].nNeighbors = vsCellCustom[icell].vNeighbors.size();
  } 
  
  
  counter_ = 0;
  for ( int icell = 0; icell < nCells; icell++ )
  {

    for ( int i = 0; i < vsCellCustom[icell].nVertices; i++ )
    {      
      vsCellCustom[icell].vVerticesNewnum.push_back(counter_);
      counter_++;
    }
  }  
  
  cout << "\t choise support cell for internal facets " << endl;
  vector<double> vDatumDistance;
  for ( int iface = 0; iface < nFaces; iface++ )
  {
    vsFaceCustom[iface].vNeighbors.clear();

    set<int>::iterator it_set;

    for ( it_set = vsetPolygonPolyhedron[iface].begin(); it_set != vsetPolygonPolyhedron[iface].end(); ++it_set )
      vsFaceCustom[iface].vNeighbors.push_back ( *it_set );
    
    vsFaceCustom[iface].nNeighbors = vsFaceCustom[iface].vNeighbors.size();
    
    if( vsFaceCustom[iface].nNeighbors == 0 )
    {
      cout << "Face " << iface << " has no connected cells!" << endl;
      exit(0);
    }
    
    /// choise support cell for internal facets
    if ( vsFaceCustom[iface].nNeighbors == 2 && vsFaceCustom[iface].nMarker > 0 )
    {
      double cosa = 0;
      for ( int idx = 0; idx < 3; idx++ )
        cosa += (vsCellCustom[ vsFaceCustom[iface].vNeighbors[0] ].vCenter[idx] - vsFaceCustom[iface].vCenter[idx]) * vsFaceCustom[iface].vNormal[idx];
       
      if(cosa > 0) swap(vsFaceCustom[iface].vNeighbors[0], vsFaceCustom[iface].vNeighbors[1]);
    }
  }
  
}

void SimData::methodElementCenter(int nelem, vector<Gelement> &vsElement)
{
  int nodes_in_elem = vsElement[nelem].nVertices;

  vsElement[nelem].vCenter.resize(3, 0.0);

  for (int inodes = 0; inodes < nodes_in_elem; inodes++)
  {
    int gl_node_num = vsElement[nelem].vVertices[inodes];
    for (int idx = 0; idx < 3; idx++)
    {
      vsElement[nelem].vCenter[idx] += vvVrtxCoords[gl_node_num][idx] / nodes_in_elem;
    }
  }
  
  vsElement[nelem].center_distance = 0.0;  
  for (int idx = 0; idx < 3; idx++)
  {
    vsElement[nelem].center_distance += vsElement[nelem].vCenter[idx] * vsElement[nelem].vCenter[idx];
  }
  vsElement[nelem].center_distance = sqrt(vsElement[nelem].center_distance);  
  
  double distance, buf;
  int node_master, node_slave;
  switch ( vsElement[nelem].formIndex )
  {
  case TETRA4:
    distance = 0.0;
    node_master = vsElement[nelem].vVertices[0];
    for ( int inode = 1; inode < nodes_in_elem; ++inode )
    {
      node_slave = vsElement[nelem].vVertices[inode];
      buf = 0.0;
      for ( int idx = 0; idx < 3; idx++ )
      {
        buf += ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] ) * ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] );
      }
      distance += sqrt ( buf );

    }
    vsElement[nelem].thickness = distance / (nodes_in_elem-1);
    break;    
  case TETRA10:
    distance = 0.0;
    node_master = vsElement[nelem].vVertices[0];
    for ( int inode = 1; inode < nodes_in_elem; ++inode )
    {
      node_slave = vsElement[nelem].vVertices[inode];
      buf = 0.0;
      for ( int idx = 0; idx < 3; idx++ )
      {
        buf += ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] ) * ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] );
      }
      distance += sqrt ( buf );

    }
    vsElement[nelem].thickness = distance / (nodes_in_elem-1);
    break;   
  case PRISM6:
    distance = 0.0;
    for ( int inode = 0; inode < nodes_in_elem / 2; ++inode )
    {      
      node_master = vsElement[nelem].vVertices[inode];
      node_slave = vsElement[nelem].vVertices[inode + nodes_in_elem / 2];
      buf = 0.0;
      for ( int idx = 0; idx < 3; idx++ )
      {
        buf += ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] ) * ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] );
      }
      distance += sqrt ( buf );
    }
    vsElement[nelem].thickness = distance / (nodes_in_elem / 2);
    break;     
  case PRISM8:
    distance = 0.0;
    for ( int inode = 0; inode < nodes_in_elem / 2; ++inode )
    {      
      node_master = vsElement[nelem].vVertices[inode];
      node_slave = vsElement[nelem].vVertices[inode + nodes_in_elem / 2];
      buf = 0.0;
      for ( int idx = 0; idx < 3; idx++ )
      {
        buf += ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] ) * ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] );
      }
      distance += sqrt ( buf );
    }
    vsElement[nelem].thickness = distance / (nodes_in_elem / 2);
    break;     
  default:
    break;
  }
}


void SimData::methodFaceNormalVector(int nelem, vector<Gelement> &vsElement)
{  

  for (int i = 0; i < 3; i++) vvPlate[i].assign(3, 0.0);

  for (int inodes = 0; inodes < 3; inodes++)
  {
    int gl_node_num = vsElement[nelem].vVertices[vPointPass[inodes]];
    for (int idx = 0; idx < 3; idx++)
    {
      vvPlate[inodes][idx] = vvVrtxCoords[gl_node_num][idx];
    }
  }
  // calculate basis
  vPointCoord[0] = vvPlate[1][1] * vvPlate[2][2] - vvPlate[1][2] * vvPlate[2][1] +
                   vvPlate[2][1] * vvPlate[0][2] - vvPlate[0][1] * vvPlate[2][2] +
                   vvPlate[0][1] * vvPlate[1][2] - vvPlate[1][1] * vvPlate[0][2];

  vPointCoord[1] = vvPlate[0][0] * vvPlate[2][2] - vvPlate[0][0] * vvPlate[1][2] +
                   vvPlate[1][0] * vvPlate[0][2] - vvPlate[1][0] * vvPlate[2][2] +
                   vvPlate[2][0] * vvPlate[1][2] - vvPlate[2][0] * vvPlate[0][2];

  vPointCoord[2] = -vvPlate[1][0] * vvPlate[0][1] + vvPlate[0][0] * vvPlate[1][1] -
                   vvPlate[0][0] * vvPlate[2][1] + vvPlate[1][0] * vvPlate[2][1] +
                   vvPlate[2][0] * vvPlate[0][1] - vvPlate[2][0] * vvPlate[1][1];

  //calculate length
  double length = 0.0;
  for (int i = 0; i < 3; i++) length += vPointCoord[i] * vPointCoord[i];

  length = sqrt(length);

  // calculate normal basis
  vsElement[nelem].vNormal.resize(3, 0.0);

  for (int i = 0; i < 3; i++) vsElement[nelem].vNormal[i] = vPointCoord[i] / length;

  
}

void SimData::methodChangeFacesNormalVector()
{  
  vector<int> vVertices(256,0);
  pair<set<int>::iterator, bool> pairIterBool;
  vector<double> vDatumNormal;
  vDatumNormal.resize(3, 0.0);
  
  vector<int> vFacevVertices;
  
  nExternalBoundaryFaces = 0;
  nInternalBoundaryFaces = 0;
  
  for(int iface = 0; iface < nFaces; iface++)
  {
    if( vsFaceCustom[iface].nMarker < 0)
    {
      setIdenticalExternalMarker.insert( vsFaceCustom[iface].nMarker );
      nExternalBoundaryFaces++;
    }
    
    
    if( vsFaceCustom[iface].nMarker > 0) 
    {
      pairIterBool = setIdenticalInternalMarker.insert( vsFaceCustom[iface].nMarker );
      double a = 1e-3;
      vIdenticalInternalFacetPerm.push_back ( a*a/12/1e-15 );
      vIdenticalInternalFacetAperture.push_back ( 0.001 );
      vIdenticalInternalFacetFFpermMult.push_back ( 1.0 );
      
      nInternalBoundaryFaces++;
    }
  }
  set<int>::iterator itintset;
  for (itintset = setIdenticalInternalMarker.begin(); itintset != setIdenticalInternalMarker.end(); ++itintset)
  {
    // datum normal vector
    for(int iface = 0; iface < nFaces; iface++)
    {
      if( vsFaceCustom[iface].nMarker == *itintset)
      {
        vDatumNormal = vsFaceCustom[iface].vNormal;
        break;
      }
    }

    for(int iface = 0; iface < nFaces; iface++)
    {
      if( vsFaceCustom[iface].nMarker == *itintset)
      {
        double cosa = 0;
        for(int idx = 0; idx < 3; idx++) cosa += vDatumNormal[idx] * vsFaceCustom[iface].vNormal[idx];

        // non collinear vector. change verticies order
        if(cosa < 0.0)
        {
          if(vsFaceCustom[iface].formIndex == QUAD8)
          {
           reverse(vsFaceCustom[iface].vVertices.begin(),vsFaceCustom[iface].vVertices.begin()+4);
           for(unsigned int i = 4; i < 8; ++i)
             vVertices[i] = vsFaceCustom[iface].vVertices[i];
           vsFaceCustom[iface].vVertices[4] = vVertices[6];
           vsFaceCustom[iface].vVertices[5] = vVertices[5];
           vsFaceCustom[iface].vVertices[6] = vVertices[4];
           vsFaceCustom[iface].vVertices[7] = vVertices[7];
          }
          else if(vsFaceCustom[iface].formIndex == TRGLE6)
          {
            reverse(vsFaceCustom[iface].vVertices.begin(), vsFaceCustom[iface].vVertices.begin() + 3);
            for(unsigned int i = 3; i < 6; ++i)
              vVertices[i] = vsFaceCustom[iface].vVertices[i];
            vsFaceCustom[iface].vVertices[3] = vVertices[4];
            vsFaceCustom[iface].vVertices[4] = vVertices[3];
            vsFaceCustom[iface].vVertices[5] = vVertices[5];
          }
          else
          {
           reverse(vsFaceCustom[iface].vVertices.begin(),vsFaceCustom[iface].vVertices.end());
          }          

          for(int idx = 0; idx < 3; idx++) vsFaceCustom[iface].vNormal[idx] *= -1.0;
        }
      }
    }
  }    
}

void SimData::methodRandomRockProperties()
{
 for (int i = 0; i < nCells; i++)
 {   
   double x = vsCellCustom[i].vCenter[0];
   double y = vsCellCustom[i].vCenter[1];
   
   double lognormal_value = createLognormalDistribution(0.01, 0.001);
   double value = lognormal_value * 
                  1.0/(1.0+exp((x-200) / 5)) * 1.0/(1.0+exp(-(x-20) / 5)) * 
                  1.0/(1.0+exp((y-200) / 5)) * 1.0/(1.0+exp(-(y-20) / 5)); 
   
   vsCellRockProps[i].perm = value;

   lognormal_value = createLognormalDistribution(2.0, 0.0001);
   value = lognormal_value;   
   vsCellRockProps[i].young = value;
 }
}

double SimData::createLognormalDistribution(double E, double S)
{
  vector<double> vVars;
  double value = 0.0;
  
  for(int i = 0; i < 12; i++)
  {
    vVars.push_back( (rand() % 100 + 1) / 100.);    
  }
  
  for(int i = 0; i < 12; i++)
  {
    value += vVars[i];
  }
  value -= 6.0;
  
  return(E + S * value);
}


void SimData::splitInternalFaces()
{
  int counter_;  
  cout << endl << "Create set of stick vertices (slow)" << endl;  
  vector<set<int> > vsetGlueVerticies;
  counter_ = 0;
  for ( int icell = 0; icell < nCells; icell++ )
  {
    counter_ += vsCellCustom[icell].nVertices;
  }
  vsetGlueVerticies.resize(counter_);
  for(int i = 0; i < counter_; i++) vsetGlueVerticies[i].insert(i);
    
  vector<double> vVerticesPair; vVerticesPair.resize(2, 0);
  
  for ( int iface = 0; iface < nFaces; iface++ )
  {
    int job_percent = int ( ( 100. * iface ) / ( nFaces ) );
    cout << "\r    " << job_percent << "%";
    /// non phisical face            
    if ( vsFaceCustom[ iface ].nMarker == 0 && vsFaceCustom[ iface ].nNeighbors == 2 )
    {
      vector<int>::iterator it_face_vrtx;

      /// loop by all face vertices
      for ( it_face_vrtx = vsFaceCustom[ iface ].vVertices.begin();  it_face_vrtx != vsFaceCustom[ iface ].vVertices.end(); ++it_face_vrtx)
      {
        /// loop by neighbor cells
        int n_polyhedron = 0;
        set<int>::iterator it_polyhedron;
        for ( it_polyhedron = vsetPolygonPolyhedron[iface].begin(); it_polyhedron != vsetPolygonPolyhedron[iface].end(); ++it_polyhedron, n_polyhedron++)
        {
          /// loop by all cell vertices
          int ncell = *it_polyhedron;
          int ivrtx = 0;

          vector<int>::iterator it_cell_vrtx;
          for ( it_cell_vrtx = vsCellCustom[ncell].vVertices.begin() ; it_cell_vrtx < vsCellCustom[ncell].vVertices.end(); ++it_cell_vrtx, ivrtx++ )
          {
            if ( *it_cell_vrtx == *it_face_vrtx )
              break;
          }          
          vVerticesPair[n_polyhedron] = vsCellCustom[ncell].vVerticesNewnum[ivrtx];
        }
        vsetGlueVerticies[ vVerticesPair[0] ].insert ( vVerticesPair[1] );
        vsetGlueVerticies[ vVerticesPair[1] ].insert ( vVerticesPair[0] );
      }
    }
  }
  
  cout << endl << "Distinguish authenic vertices (might be very slow)" << endl;  
  int n_possible_verticies = vsetGlueVerticies.size();
  vector<int> v_buf_storage;
  int total_size_old = 0;
  int total_size_new = 1;
  int cycle = 0;
  while ( total_size_old != total_size_new )
  {
    total_size_old = total_size_new;
    total_size_new = 0;
    cout << endl << "cycle   :" << cycle << "\t \t Hopefully < 10"; cycle++;
    
    for ( int ivrtx = vsCellCustom[0].nVertices; ivrtx < n_possible_verticies; ivrtx++ )
    {
      int job_percent = int ( ( 100. * ivrtx ) / ( n_possible_verticies ) );
      cout << "\r    " << job_percent << "%";
      v_buf_storage.clear();
      set<int>::iterator it_set;

      for ( it_set = vsetGlueVerticies[ivrtx].begin(); it_set != vsetGlueVerticies[ivrtx].end(); ++it_set )
      {
        set<int>::iterator it_set_down;

        for ( it_set_down = vsetGlueVerticies[ *it_set ].begin(); it_set_down != vsetGlueVerticies[ *it_set ].end(); ++it_set_down )
        {
          v_buf_storage.push_back ( *it_set_down );
        }
      }

      vector<int>::iterator it_vec;

      for ( it_vec = v_buf_storage.begin(); it_vec != v_buf_storage.end(); it_vec++ )
      {
        vsetGlueVerticies[ivrtx].insert ( *it_vec );
      }
      total_size_new += vsetGlueVerticies[ivrtx].size();
    }
  }
      
  cout << endl << "Renumber vector of stick vertices" << endl;
  vector<int> vRenumVerticies;
  vRenumVerticies.resize(n_possible_verticies);
  
  /// take first cell
  for(int ivrtx = 0; ivrtx < vsCellCustom[0].nVertices; ivrtx++)
  {
    vRenumVerticies[ivrtx] = ivrtx;
  }
  
  counter_ = vsCellCustom[0].nVertices;
  for(int ivrtx = vsCellCustom[0].nVertices; ivrtx < n_possible_verticies; ivrtx++)
  {
    if( *vsetGlueVerticies[ivrtx].begin() == ivrtx )
    {
      // this vertex is not stick
      vRenumVerticies[ivrtx] = counter_;
      counter_++;
    }
    else
    {
      // this vertex is stick
      // set is sorted by c++ defaults, and we take minimum value
      vRenumVerticies[ivrtx] = vRenumVerticies[ *vsetGlueVerticies[ivrtx].begin() ];
    }
  }
  
  for(int icell = 0; icell < nCells; icell++)
  {
    for(int ivrtx = 0; ivrtx < vsCellCustom[icell].nVertices; ivrtx++)
    {
      vsCellCustom[icell].vVerticesNewnum[ ivrtx ] = vRenumVerticies[ vsCellCustom[icell].vVerticesNewnum[ ivrtx ] ];
    }
  }
  int nTotalNodes = counter_;
  
  cout << "\t check renumbering consistency " << endl;
  vector<double> vStatus;
  vStatus.resize(counter_, false);
  for ( int icell = 0; icell < nCells; icell++ )
  {
    for ( int ivrtx = 0; ivrtx < vsCellCustom[icell].nVertices; ivrtx++ )
      vStatus[ vsCellCustom[icell].vVerticesNewnum[ivrtx] ] = true;
  }

  for ( int i = 0; i < counter_; i++ )
  {
    if ( vStatus[i] == false )
      cout << i << endl;
  }
    
  cout << "\t change face numbering" << endl;
  for(int iface = 0; iface < nFaces; iface++)
  {
    vsFaceCustom[iface].vVerticesNewnum.resize( vsFaceCustom[iface].nVertices, -1);
    // we take always [0] - support cell
    int icell = vsFaceCustom[iface].vNeighbors[0];
    for(int ivrtx = 0; ivrtx < vsFaceCustom[iface].nVertices; ivrtx++)
    {
      for(int inode = 0; inode < vsCellCustom[icell].nVertices; inode++)
      {
        if( vsFaceCustom[iface].vVertices[ivrtx] == vsCellCustom[icell].vVertices[inode] )
        {
          vsFaceCustom[iface].vVerticesNewnum[ivrtx] = vsCellCustom[icell].vVerticesNewnum[inode];
          break;
        }
      }
    }
  }

  cout << "\t create atoms list" << endl;
  set<int>::iterator itintset;
  vector<int> vTempAtoms; 
  vTempAtoms.resize( nTotalNodes, -999 );
  for (itintset = setIdenticalInternalMarker.begin(); itintset != setIdenticalInternalMarker.end(); ++itintset)
  {
    for(int iface = 0; iface < nFaces; iface++)
    {
      int ncell = vsFaceCustom[iface].vNeighbors[0];   //EMIL
      if( vsFaceCustom[iface].nMarker == *itintset)
      {
        for(int ivrtx = 0; ivrtx < vsFaceCustom[iface].nVertices; ivrtx++)
        {
          for(int inode = 0; inode < vsCellCustom[ncell].nVertices; inode++)
          {
            if( vsFaceCustom[iface].vVertices[ivrtx] == vsCellCustom[ncell].vVertices[inode] )
            {
              vTempAtoms[ vsCellCustom[ncell].vVerticesNewnum[inode] ] = vsFaceCustom[iface].vVerticesNewnum[ivrtx];
            }
          }
        }
      }
    }
  }
        
  cout << endl << "Change coordanates vector" << endl;
  vector<set<int> > vsetOld2New;
  vsetOld2New.resize(nNodes);
  vector<vector<double> > vvNewCoordinates;
  vvNewCoordinates.resize( nTotalNodes, vector<double>(3,0) );
  for(int icell = 0; icell < nCells; icell++)
  {
    vector<int>::iterator it_old, it_new;
    it_new = vsCellCustom[icell].vVerticesNewnum.begin();
    for(it_old = vsCellCustom[icell].vVertices.begin(); it_old != vsCellCustom[icell].vVertices.end(); ++it_old, ++it_new)
    {
      vvNewCoordinates[ *it_new ] = vvVrtxCoords[ *it_old ];
      vsetOld2New[*it_old].insert(*it_new);
    }
  }
  
  vvVrtxCoords.resize(nTotalNodes, vector<double>(3,0));
  for(int ivrtx = 0; ivrtx < nTotalNodes; ivrtx++)
  {
    vvVrtxCoords[ivrtx] = vvNewCoordinates[ivrtx];
  }

  cout << endl << "Unify previous and splitted data" << endl;
  cout << "Verticies : " << nNodes << "\t \t After splitting : " << nTotalNodes << endl;
  nNodes = nTotalNodes;
  for(int icell = 0; icell < nCells; icell++)
  {
    vsCellCustom[icell].vVertices = vsCellCustom[icell].vVerticesNewnum;
  }
  
  for(int iface = 0; iface < nFaces; iface++)
  {
    vsFaceCustom[iface].vVertices = vsFaceCustom[iface].vVerticesNewnum;
  }

  vvAtoms.resize(nNodes, vector<int>(2,0) );
  nAtoms = 0;
  for(int iatom = 0; iatom < nNodes; iatom++)
  {
    if(vTempAtoms[iatom] >= 0) 
    {
      vvAtoms[nAtoms][0] = iatom;
      vvAtoms[nAtoms][1] = vTempAtoms[iatom];
      nAtoms++;
    }
  }
  
  // recast boundary conditions
  cout << "recast boundary conditions stage 1" << endl;
  vector<int> vDirichletNode_;  
  vector<int> vDirichletNodeIdx_;  
  vector<double> vDirichletNodeVal_;    
  for ( int i = 0; i < nDirichletNodes; ++i )
  {    
    for(set<int>::iterator it = vsetOld2New[vDirichletNode[i]].begin(); it != vsetOld2New[vDirichletNode[i]].end(); ++it)
    {
      vDirichletNode_.push_back(*it);
      vDirichletNodeIdx_.push_back(vDirichletNodeIdx[i]);
      vDirichletNodeVal_.push_back(vDirichletNodeVal[i]);
    }
  }
  
  cout << "recast boundary conditions stage 2" << endl;
  nDirichletNodes = vDirichletNode_.size();
  vDirichletNode.assign(nDirichletNodes,0);
  vDirichletNodeIdx.assign(nDirichletNodes,0);
  vDirichletNodeVal.assign(nDirichletNodes,0);
  for ( int i = 0; i < nDirichletNodes; ++i )
  {
    vDirichletNode[i] = vDirichletNode_[i];
    vDirichletNodeIdx[i] = vDirichletNodeIdx_[i];
    vDirichletNodeVal[i] = vDirichletNodeVal_[i];
  }
  
  //@EMIL We could not provide a smart renumbering for a Emil's meshes
  return;
  cout << endl << "Smart verticies renumbering" << endl;
  vector<int> vNodesID;
  vector<int> vIA;
  vector<int> vJA;
  vector<int> vRCM;
  
  vIA.push_back(0);
  for(int ic = 0; ic < nCells; ++ic)
  {
    int n = vsCellCustom[ic].vVertices.size();
    for(int iv = 0; iv < n; ++iv)
    {
      vJA.push_back( vsCellCustom[ic].vVertices[iv] );
    }
    n += vIA[ic];
    vIA.push_back(n);
  }
  
  vRCM.assign(nNodes,-1);
  pRenum->convert(nCells, nNodes, vIA, vJA, vRCM);
  //pRenum->convertEmil( vRCM, vvVrtxCoords);

  for(int icell = 0; icell < nCells; icell++)
  {
    for(int iv = 0; iv < vsCellCustom[icell].vVertices.size(); iv++)
      vsCellCustom[icell].vVertices[iv] = vRCM[vsCellCustom[icell].vVertices[iv]];
  }
 
  for(int iface = 0; iface < nFaces; iface++)
  {
    for(int iv = 0; iv < vsFaceCustom[iface].vVertices.size(); iv++)
      vsFaceCustom[iface].vVertices[iv] = vRCM[vsFaceCustom[iface].vVertices[iv]];
  }
  
  for(int iatom = 0; iatom < nNodes; iatom++)
  {
    if(vTempAtoms[iatom] >= 0) 
    {
      vvAtoms[nAtoms][0] = vRCM[vvAtoms[nAtoms][0]];
      vvAtoms[nAtoms][1] = vRCM[vvAtoms[nAtoms][1]];
    }
  }
  
  for(int ivrtx = 0; ivrtx < nTotalNodes; ivrtx++)
    vvNewCoordinates[vRCM[ivrtx]] = vvVrtxCoords[ivrtx];
  
  for(int ivrtx = 0; ivrtx < nTotalNodes; ivrtx++)
    vvVrtxCoords[ivrtx] = vvNewCoordinates[ivrtx];  
}
 

void SimData::handleConnections()
{  
  // we always start to count fractures
  int counter = 0;
  for ( int iface = 0; iface < nFaces; iface++ )
  {
    if ( vsFaceCustom[iface].nMarker > 0 && vsFaceCustom[iface].nMarker < 1111110 )
    {
      vsFaceCustom[iface].fluidElement = counter; counter++;
    }

    if ( vsFaceCustom[iface].nMarker < 0 )
      vsFaceCustom[iface].fluidElement = -1;
  }
  
  // then we count cells
  for(int icell = 0; icell < nCells; icell++)
  {
    vsCellCustom[icell].fluidElement = counter; counter++;
  }
}

#include <random>
void SimData::definePhysicalFacets()
{
  int nbnd = vsPhysicalBoundary.size();
  
  nNeumannFaces = 0;
  nDirichletNodes = 0;
  int nDirichletFaces = 0;
   
  int nfacets = 0;
  int nfluid = 0;
  vsPhysicalFacet.resize(nFaces);
  vector<set<int>> vDirichletNodesSorted;
  vDirichletNodesSorted.resize(3);
  pair<set<int>::iterator, bool> ret;
  for(int iface = 0; iface < nFaces; iface++)
  {
    if( vsFaceCustom[iface].nMarker < 0)
    {
      for(int i = 0; i < nbnd; i++)
      {
        if( vsFaceCustom[iface].nMarker == vsPhysicalBoundary[i].nmark)
        {
          vsPhysicalFacet[nfacets].nface = nfacets;
          vsPhysicalFacet[nfacets].ntype = vsPhysicalBoundary[i].ntype;
          vsPhysicalFacet[nfacets].nmark = vsPhysicalBoundary[i].nmark;
          vsPhysicalFacet[nfacets].vCondition.resize( vsPhysicalBoundary[i].vCondition.size() );          
          vsPhysicalFacet[nfacets].vCondition = vsPhysicalBoundary[i].vCondition; 
          vsPhysicalFacet[nfacets].nfluid = -1;
	  
	  if(vsPhysicalBoundary[i].ntype == 1)
	  {
	  for(unsigned int iv = 0; iv < vsFaceCustom[iface].nVertices; ++iv)
	  {
	    for(unsigned int k = 0; k < 3; ++k)
	    {
	      if(vsPhysicalFacet[nfacets].vCondition[k] != dNotNumber)
	      {
			  ret = vDirichletNodesSorted[k].insert(vsFaceCustom[iface].vVertices[iv]);
			  if (ret.second == true)
			  {
				  vDirichletNode.push_back(vsFaceCustom[iface].vVertices[iv]);
				  vDirichletNodeIdx.push_back(k);
				  vDirichletNodeVal.push_back(vsPhysicalFacet[nfacets].vCondition[k]);
				  nDirichletNodes++;
			  }
		  }
	    }
	  }
	  }
	  
          if(vsPhysicalBoundary[i].ntype == 1)
          {
            nDirichletFaces++;
          }
          else
          {
            nNeumannFaces++;
          }
          nfacets++;
        }
      }
    }
    
    if ( vsFaceCustom[iface].nMarker > 0 && vsFaceCustom[iface].nMarker < 1111110 )
    {
      vsPhysicalFacet[nfacets].nface = nfacets;
      vsPhysicalFacet[nfacets].ntype = 0;
      vsPhysicalFacet[nfacets].nmark = vsFaceCustom[iface].nMarker;
      vsPhysicalFacet[nfacets].nfluid = nfluid;
      nfacets++;
      nfluid++;
    }
  }
  nPhysicalFacets = nfacets; 
  
  std::default_random_engine generator;
  std::default_random_engine generator2;
  std::normal_distribution<double> distribution(100.0, 20.0);
  std::normal_distribution<double> distribution2(10.0, 2.0); //mD
  for(int iface = 0; iface < nFaces; iface++)
  {    
    vsFaceCustom[iface].aperture = 1e-5; 
    vsFaceCustom[iface].conductivity = 2.5e-004;
    //@EMIL
    vsFaceCustom[iface].active_segment = 0; // passive   
	if (vsFaceCustom[iface].nMarker == 1)
	{
		vsFaceCustom[iface].active_segment = 1;    // active
		vsFaceCustom[iface].conductivity = 2.5e-004;
	}      
  }
}

void SimData::defineStressAndDispOnBoundary()
{
  PhysicalFace pfFace;
  int inputBC = vvsBCIn.size();
  
  vvsBCOut.resize( inputBC );
  for(int ibnd = 0; ibnd < inputBC; ibnd++)
  {
    int nfacets = -1;
    for ( int iface = 0; iface < nFaces; iface++ )
    {
      if ( vsFaceCustom[iface].nMarker != 0 ) nfacets++;
      
      if ( vsFaceCustom[iface].nMarker < 0 )
      {
        for(int i = 0; i < vvsBCIn[ibnd].size(); i++ )
        {
          if( vsFaceCustom[iface].nMarker == vvsBCIn[ibnd][i].nmark )
          {
            pfFace.nface = nfacets;
            pfFace.ntype = vvsBCIn[ibnd][i].ntype;
            pfFace.nmark = vvsBCIn[ibnd][i].nmark;
            pfFace.nfluid = -1;
            pfFace.axle = vvsBCIn[ibnd][i].axle;
            
            pfFace.vCondition.clear(); 
            pfFace.vCondition.push_back( vvsBCIn[ibnd][i].vCondition[0] );                                  
            vvsBCOut[ibnd].push_back( pfFace );
          }
        }
      }
    }    
  }
}


void SimData::createSimpleWells()
{
  vector<double> vCenter(3,0);
  /// choise support cell for internal facets
  for ( int iwell = 0; iwell < nWells; iwell++ )
  {
    int icell = 0;
    for ( int iface = 0; iface < nFaces; iface++ )
    {
      vCenter[0] = vsFaceCustom[iface].vCenter[0];
      vCenter[1] = vsFaceCustom[iface].vCenter[1];
      vCenter[2] = vsFaceCustom[iface].vCenter[2];
      
	  if (vsFaceCustom[iface].nMarker > 0 )
      {
	if ( abs ( vCenter[0] - vsWell[iwell].vWellCoordinate[0] ) < vsWell[iwell].radius_poisk && abs ( vCenter[1] - vsWell[iwell].vWellCoordinate[1] ) < vsWell[iwell].radius_poisk )
        {
          vsWell[iwell].vID.push_back ( icell );
          vsWell[iwell].vWi.push_back ( 100.0 );
        }
        ++icell;
      }
    }
  }
  /*
  // MATRIX PART
  for ( int iwell = 0; iwell < nWells; iwell++ )
  {
    int icell = 0;
    for ( int ic = 0; ic < nCells; ic++ )
    {
      vCenter[0] = vsCellCustom[ic].vCenter[0];
      vCenter[1] = vsCellCustom[ic].vCenter[1];
      vCenter[2] = vsCellCustom[ic].vCenter[2];
      if ( abs ( vCenter[0] - vsWell[iwell].vWellCoordinate[0] ) < vsWell[iwell].radius_poisk && abs ( vCenter[1] - vsWell[iwell].vWellCoordinate[1] ) < vsWell[iwell].radius_poisk )
      {
        vsWell[iwell].vID.push_back ( ic );
        vsWell[iwell].vWi.push_back ( 1.0 );
      }
    }
  }
  */
  for ( int iwell = 0; iwell < nWells; iwell++ )
  {
    vsWell[iwell].datum = 1e16;
    for(int i = 0; i < vsWell[iwell].vID.size(); ++i)
    {
      int icell = vsWell[iwell].vID.size();
      if( fabs(vsCellCustom[icell].vCenter[2]) < vsWell[iwell].datum )
	vsWell[iwell].datum = abs(vsCellCustom[icell].vCenter[2]);
    }
  } 
}

