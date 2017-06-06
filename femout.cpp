#include "femout.hpp"

OutputData::OutputData(SimData * pSimData)
{
  pSim = pSimData;
}

OutputData::~OutputData()
{
}

void OutputData::writeGeomechDataNewKeywords()
{
	string outstring = pSim->outstream + ".geomech";
	ofstream geomechfile;
	outstring = "fl_dimens.txt";
	geomechfile.open(outstring.c_str());
	geomechfile << "DIMENS" << endl;
	geomechfile << pSim->nInternalBoundaryFaces + pSim->nCells + pSim->nGhostCells << "\t" << 1 << "\t" << 1 << " /" << endl;
	geomechfile.close();

	cout << pSim->nInternalBoundaryFaces + pSim->nCells << "\t" << pSim->nGhostCells << endl;

	// GEOMETRY
	outstring = "gm_geometry.txt";
	geomechfile.open(outstring.c_str());
	geomechfile << "GMDIMS" << endl;
	geomechfile << pSim->nNodes << "\t" << pSim->nCells << "\t" << pSim->nFaces;
	geomechfile << "/" << endl << endl;

	geomechfile.precision(12);
	cout << "write all coordinates\n";
	geomechfile << "GMNODE_COORDS" << endl;
	for (int i = 0; i < pSim->nNodes; i++)
	{
		geomechfile << pSim->vvVrtxCoords[i][0] << "\t" << pSim->vvVrtxCoords[i][1] << "\t" << pSim->vvVrtxCoords[i][2] << "\n";
	}
	geomechfile << "/" << endl << endl;

	cout << "write all elements\n";
	geomechfile << "GMCELL_NODES" << endl;
	for (int i = 0; i < pSim->nCells; i++)
	{
		geomechfile << pSim->vsCellCustom[i].nVertices << "\t";
		if (pSim->vsCellCustom[i].vtkIndex == 25)
		{
			// super wierd element 25
			for (int j = 0; j < 8; j++)
				geomechfile << pSim->vsCellCustom[i].vVertices[j] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[8] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[11] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[13] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[9] + 1 << "\t";

			geomechfile << pSim->vsCellCustom[i].vVertices[16] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[18] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[19] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[17] + 1 << "\t";

			geomechfile << pSim->vsCellCustom[i].vVertices[10] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[12] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[14] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[15] + 1 << "\t";
			geomechfile << endl;
		}
		else if (pSim->vsCellCustom[i].vtkIndex == 26)
		{
			for (int j = 0; j < 6; j++)
				geomechfile << pSim->vsCellCustom[i].vVertices[j] + 1 << "\t";

			geomechfile << pSim->vsCellCustom[i].vVertices[6] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[9] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[7] + 1 << "\t";

			geomechfile << pSim->vsCellCustom[i].vVertices[12] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[14] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[13] + 1 << "\t";

			geomechfile << pSim->vsCellCustom[i].vVertices[8] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[10] + 1 << "\t";
			geomechfile << pSim->vsCellCustom[i].vVertices[11] + 1 << "\t";
			geomechfile << endl;
		}
		else
		{
			for (int j = 0; j < pSim->vsCellCustom[i].nVertices; j++)
				geomechfile << pSim->vsCellCustom[i].vVertices[j] + 1 << "\t";
			geomechfile << endl;
		}
	}
	geomechfile << "/" << endl << endl;


	geomechfile << "GMCELL_TYPE" << endl;
	for (int i = 0; i < pSim->nCells; i++)
		geomechfile << pSim->vsCellCustom[i].vtkIndex << endl;
	geomechfile << "/" << endl << endl;

	geomechfile << "GMCELL_TO_FLOWCELLS" << endl;
	int i_ = 0;
	for (int i = 0; i < pSim->nCells; i++)
	{
		geomechfile << 1 << "\t" << pSim->vsCellCustom[i].fluidElement + 1 << endl;
	}
	geomechfile << "/\n\n";

	cout << "write all faces\n";
	geomechfile << "GMFACE_NODES\n";
	for (int i = 0; i < pSim->nFaces; i++)
	{
		geomechfile << pSim->vsFaceCustom[i].nVertices << "\t";

		for (int j = 0; j < pSim->vsFaceCustom[i].nVertices; j++) geomechfile << pSim->vsFaceCustom[i].vVertices[j] + 1 << "\t";

		geomechfile << endl;
	}
	geomechfile << "/" << endl << endl;

	geomechfile << "GMFACE_TYPE" << endl;
	for (int i = 0; i < pSim->nFaces; i++)
		geomechfile << pSim->vsFaceCustom[i].vtkIndex << endl;
	geomechfile << "/" << endl << endl;

	geomechfile << "GMFACE_GMCELLS" << endl;
	for (int i = 0; i < pSim->nFaces; i++)
	{
		geomechfile << pSim->vsFaceCustom[i].nNeighbors << "\t";

		for (int j = 0; j < pSim->vsFaceCustom[i].nNeighbors; j++)
			geomechfile << pSim->vsFaceCustom[i].vNeighbors[j] + 1 << "\t";

		geomechfile << endl;
	}
	geomechfile << "/" << endl << endl;
	geomechfile.close();

	outstring = "gm_model.txt";
	geomechfile.open(outstring.c_str());
	geomechfile << "GMCELL_MODEL" << endl;
	for (int ib = 0; ib < pSim->nCells; ++ib)
	{
		geomechfile << pSim->vsCellRockProps[ib].model << endl;
	}
	geomechfile << "/" << endl << endl;
	geomechfile.close();

	outstring = "gm_capacity.txt";
	geomechfile.open(outstring.c_str());
	geomechfile << "GMCELL_HEAT_CAPACITY" << endl;
	for (int ib = 0; ib < pSim->nCells; ++ib)
	{
		geomechfile << pSim->vsCellRockProps[ib].heat_capacity << endl;
	}
	geomechfile << "/" << endl << endl;
	geomechfile.close();

	outstring = "gm_density.txt";
	geomechfile.open(outstring.c_str());

	geomechfile << "GMCELL_DENSITY\n";
	int j = 0;
	for (int i = 0; i < pSim->nCells; i++, j++)
	{
		geomechfile << pSim->vsCellRockProps[i].density << "\t";
		if (j == 10)
		{
			j = 0;
			geomechfile << endl;
		}
	}
	geomechfile << "\n/\n\n";
	geomechfile.close();

	outstring = "gm_biot.txt";
	geomechfile.open(outstring.c_str());
	geomechfile << "GMCELL_BIOT\n";
	j = 0;
	for (int i = 0; i < pSim->nCells; i++, j++)
	{
		geomechfile << pSim->vsCellRockProps[i].biot << "\t";
		if (j == 10)
		{
			j = 0;
			geomechfile << endl;
		}
	}
	geomechfile << "\n/\n\n";

	geomechfile << "GMCELL_BIOT_FLOW\n";
	j = 0;
	for (int i = 0; i < pSim->nCells; i++, j++)
	{
		geomechfile << pSim->vsCellRockProps[i].biot_flow << "\t";
		if (j == 10)
		{
			j = 0;
			geomechfile << endl;
		}
	}
	geomechfile << "\n/\n\n";
	geomechfile.close();

	outstring = "gm_elastic.txt";
	geomechfile.open(outstring.c_str());
	geomechfile << "GMCELL_YOUNG\n";
	j = 0;
	for (int i = 0; i < pSim->nCells; i++, j++)
	{
		geomechfile << pSim->vsCellRockProps[i].young << "\t";
		if (j == 10)
		{
			j = 0;
			geomechfile << endl;
		}
	}
	geomechfile << "\n/\n\n";

	geomechfile << "GMCELL_POISSON\n";
	j = 0;
	for (int i = 0; i < pSim->nCells; i++, j++)
	{
		geomechfile << pSim->vsCellRockProps[i].poisson << "\t";
		if (j == 8)
		{
			j = 0;
			geomechfile << endl;
		}
	}
	geomechfile << "\n/\n\n";
	geomechfile.close();

	outstring = "gm_thermal.txt";
	geomechfile.open(outstring.c_str());

	geomechfile << "GMCELL_THERMAL_EXPANSION\n";
	j = 0;
	for (int i = 0; i < pSim->nCells; i++, j++)
	{
		geomechfile << pSim->vsCellRockProps[i].thermal_expansion << "\t";
		if (j == 10)
		{
			j = 0;
			geomechfile << endl;
		}
	}
	geomechfile << "\n/\n\n";

	geomechfile << "GMCELL_PORE_THERMAL_EXPANSION\n";
	j = 0;
	for (int i = 0; i < pSim->nCells; i++, j++)
	{
		geomechfile << pSim->vsCellRockProps[i].pore_thermal_expansion << "\t";
		if (j == 10)
		{
			j = 0;
			geomechfile << endl;
		}
	}
	geomechfile << "\n/\n\n";
	geomechfile.close();


	outstring = "gm_perm_update.txt";
	geomechfile.open(outstring.c_str());

	geomechfile << "GMCELL_PERM_PORO_POW_N\n";
	j = 0;
	for (int i = 0; i < pSim->nCells; i++, j++)
	{
		geomechfile << pSim->vsCellRockProps[i].poron << "\t";
		if (j == 10)
		{
			j = 0;
			geomechfile << endl;
		}
	}
	geomechfile << "\n/\n\n";
	geomechfile.close();


	outstring = "gm_plastic.txt";
	geomechfile.open(outstring.c_str());

	geomechfile << "GMCELL_COHESION\n";
	j = 0;
	for (int i = 0; i < pSim->nCells; i++, j++)
	{
		geomechfile << pSim->vsCellRockProps[i].cohesion << "\t";
		if (j == 8)
		{
			j = 0;
			geomechfile << endl;
		}
	}
	geomechfile << "\n/\n\n";

	geomechfile << "GMCELL_FRICTION\n";
	j = 0;
	for (int i = 0; i < pSim->nCells; i++, j++)
	{
		geomechfile << pSim->vsCellRockProps[i].friction << "\t";
		if (j == 8)
		{
			j = 0;
			geomechfile << endl;
		}
	}
	geomechfile << "\n/\n\n";

	geomechfile << "GMCELL_DILATION\n";
	j = 0;
	for (int i = 0; i < pSim->nCells; i++, j++)
	{
		geomechfile << pSim->vsCellRockProps[i].dilation << "\t";
		if (j == 8)
		{
			j = 0;
			geomechfile << endl;
		}
	}
	geomechfile << "\n/\n\n";

	geomechfile.close();


	outstring = "gm_reference.txt";
	geomechfile.open(outstring.c_str());

	geomechfile << "GMREF_TEMPERATURE" << endl;
	for (int ib = 0; ib < pSim->nCells; ++ib)
	{
		geomechfile << pSim->vsCellRockProps[ib].ref_temp << endl;
	}
	geomechfile << "/" << endl << endl;

	geomechfile << "GMREF_PRESSURE" << endl;
	for (int ib = 0; ib < pSim->nCells; ++ib)
	{
		geomechfile << pSim->vsCellRockProps[ib].ref_pres << endl;
	}
	geomechfile << "/" << endl << endl;

	geomechfile << "/" << endl << endl;
	geomechfile.close();


	outstring = "gm_init.txt";
	geomechfile.open(outstring.c_str());
	geomechfile << "GMCELL_INITSTRESS" << endl;
	for (int ib = 0; ib < pSim->nCells; ++ib)
	{
		for (int i = 0; i < pSim->vsCellRockProps[ib].stress.size(); ++i)
			geomechfile << pSim->vsCellRockProps[ib].stress[i] << "\t";
		geomechfile << endl;
	}
	geomechfile << "/" << endl << endl;
	geomechfile.close();


	outstring = "gm_flow_stress.txt";
	geomechfile.open(outstring.c_str());
	geomechfile << "GMREF_FLOW_STRESS" << endl;
	for (int ib = 0; ib < pSim->nCells; ++ib)
	{
		for (int i = 0; i < pSim->vsCellRockProps[ib].stress.size(); ++i)
			geomechfile << pSim->vsCellRockProps[ib].stress[i] << "\t";
		geomechfile << endl;
	}
	geomechfile << "/" << endl << endl;
	geomechfile.close();

	if (pSim->vvVrtxDisp.size() > 0)
	{
		outstring = "gm_init_disp.txt";
		geomechfile.open(outstring.c_str());
		geomechfile << "GMINITDISP" << endl;
		for (int i = 0; i < pSim->nNodes; i++)
		{
			geomechfile << pSim->vvVrtxDisp[i][0] << "\t" << pSim->vvVrtxDisp[i][1] << "\t" << pSim->vvVrtxDisp[i][2] << "\n";
		}
		geomechfile << "/" << endl << endl;
	}
	geomechfile.close();


	if (pSim->vConstraintVertex.size() > 0)
	{
		outstring = "gm_constraint.txt";
		geomechfile.open(outstring.c_str());
		geomechfile << "GMEQCONSTRAINT" << endl;
		for (int i = 1; i < pSim->vConstraintVertex.size(); ++i)
		{
			geomechfile << "2" << "\t" << pSim->vConstraintVertex[i - 1] + 1 << "\t" << pSim->vConstraintVertex[i] + 1 << "\t";
			geomechfile << "2 2 1 -1 1.e10" << endl;
		}
		geomechfile << "/" << endl << endl;
	}
	geomechfile.close();

	// write fractured faces
	outstring = "gm_bcond.txt";
	geomechfile.open(outstring.c_str());

	if (pSim->nNeumannFaces > 0)
	{
		cout << "write all Neumann faces\n";
		geomechfile << "GMFACE_TRACTION_N" << endl;

		for (int i = 0; i < pSim->nPhysicalFacets; i++)
		{
			if (pSim->vsPhysicalFacet[i].ntype == 2)
			{
				geomechfile << pSim->vsPhysicalFacet[i].nface + 1 << "\t";
				geomechfile << pSim->vsPhysicalFacet[i].vCondition[0] << endl;
			}
		}
		geomechfile << "/\n\n";
	}

	/*
	if ( pSim->nDirichletFaces > 0 )
	{
	cout << "write all Dirichlet faces\n";

	geomechfile << "GMNODE_BCDISPX\n";
	for ( int i = 0; i < pSim->nPhysicalFacets; i++ )
	{
	if ( pSim->vsPhysicalFacet[i].ntype == 1 )
	{
	int n = pSim->vsPhysicalFacet[i].vCondition.size();

	int j = 0;
	if ( pSim->vsPhysicalFacet[i].vCondition[j] != pSim->dNotNumber )
	{
	int nvrtx = pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nVertices;
	for ( int ivrtx = 0; ivrtx < nvrtx; ++ivrtx )
	{
	geomechfile <<  pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].vVertices[ivrtx] + 1 << "\t";
	geomechfile << pSim->vsPhysicalFacet[i].vCondition[j] << endl;
	}
	}
	}
	}
	geomechfile << "/\n\n";

	geomechfile << "GMNODE_BCDISPY\n";
	for ( int i = 0; i < pSim->nPhysicalFacets; i++ )
	{
	if ( pSim->vsPhysicalFacet[i].ntype == 1 )
	{
	int n = pSim->vsPhysicalFacet[i].vCondition.size();

	int j = 1;
	if ( pSim->vsPhysicalFacet[i].vCondition[j] != pSim->dNotNumber )
	{
	int nvrtx = pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nVertices;
	for ( int ivrtx = 0; ivrtx < nvrtx; ++ivrtx )
	{
	geomechfile <<  pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].vVertices[ivrtx] + 1 << "\t";
	geomechfile << pSim->vsPhysicalFacet[i].vCondition[j] << endl;
	}
	}
	}
	}
	geomechfile << "/\n\n";

	geomechfile << "GMNODE_BCDISPZ\n";
	for ( int i = 0; i < pSim->nPhysicalFacets; i++ )
	{
	if ( pSim->vsPhysicalFacet[i].ntype == 1 )
	{
	int n = pSim->vsPhysicalFacet[i].vCondition.size();

	int j = 2;
	if ( pSim->vsPhysicalFacet[i].vCondition[j] != pSim->dNotNumber )
	{
	int nvrtx = pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nVertices;
	for ( int ivrtx = 0; ivrtx < nvrtx; ++ivrtx )
	{
	geomechfile <<  pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].vVertices[ivrtx] + 1 << "\t";
	geomechfile << pSim->vsPhysicalFacet[i].vCondition[j] << endl;
	}
	}
	}
	}
	geomechfile << "/\n\n";
	geomechfile.close();
	}
	*/

	if (pSim->nDirichletNodes > 0)
	{
		cout << "write all Dirichlet nodes\n";

		geomechfile << "GMNODE_BCDISPX\n";
		for (int i = 0; i < pSim->nDirichletNodes; i++)
		{
			if (pSim->vDirichletNodeIdx[i] == 0)
			{
				geomechfile << pSim->vDirichletNode[i] + 1 << "\t";
				geomechfile << pSim->vDirichletNodeVal[i] << endl;

			}
		}
		geomechfile << "/\n\n";

		geomechfile << "GMNODE_BCDISPY\n";
		for (int i = 0; i < pSim->nDirichletNodes; i++)
		{
			if (pSim->vDirichletNodeIdx[i] == 1)
			{
				geomechfile << pSim->vDirichletNode[i] + 1 << "\t";
				geomechfile << pSim->vDirichletNodeVal[i] << endl;
			}
		}
		geomechfile << "/\n\n";

		geomechfile << "GMNODE_BCDISPZ\n";
		for (int i = 0; i < pSim->nDirichletNodes; i++)
		{
			if (pSim->vDirichletNodeIdx[i] == 2)
			{
				geomechfile << pSim->vDirichletNode[i] + 1 << "\t";
				geomechfile << pSim->vDirichletNodeVal[i] << endl;
			}
		}
		geomechfile << "/\n\n";

		geomechfile.close();
	}

	outstring = "gm_fractures.txt";
	geomechfile.open(outstring.c_str());
	set<int>::iterator itsetint;

	int counter = 0;
	if (pSim->nInternalBoundaryFaces > 0)
	{
		int nFractures_ = 0;
		cout << "write all fractured faces\n";
		nFractures_ = 0;
		geomechfile << "GMFACE_FRACTURE_ACTIVE_SEGMENT" << endl;
		for (itsetint = pSim->setIdenticalInternalMarker.begin(); itsetint != pSim->setIdenticalInternalMarker.end(); itsetint++, nFractures_++)
		{
			for (int i = 0; i < pSim->nPhysicalFacets; i++)
			{
				if (pSim->vsPhysicalFacet[i].nmark == *itsetint)
				{
					geomechfile << pSim->vsFaceCustom[pSim->vsPhysicalFacet[i].nface].active_segment << endl;
				}
			}
		}
		geomechfile << "/" << endl << endl;

		geomechfile << "GMFACE_FRACTURE_TO_FLOWCELL\n";
		for (itsetint = pSim->setIdenticalInternalMarker.begin(); itsetint != pSim->setIdenticalInternalMarker.end(); itsetint++, nFractures_++)
		{
			for (int i = 0; i < pSim->nPhysicalFacets; i++)
			{
				if (pSim->vsPhysicalFacet[i].nmark == *itsetint)
				{
					geomechfile << pSim->vsPhysicalFacet[i].nface + 1 << "\t";
					geomechfile << pSim->vsPhysicalFacet[i].nfluid + 1 << endl;
					if (pSim->vsFaceCustom[pSim->vsPhysicalFacet[i].nface].nNeighbors != 2)
					{
						cout << "Fracture interface # " << nFractures_ << endl;
						cout << "Global interface   # " << pSim->vsPhysicalFacet[i].nface << endl;
						cout << "Number of neighbors  " << pSim->vsFaceCustom[pSim->vsPhysicalFacet[i].nface].nNeighbors << endl;
						cout << "Wrong msh file. Mesh verticies are not connected on fracture interface" << endl;
						exit(0);
					}
				}
			}
		}
		geomechfile << "/" << endl << endl;

		nFractures_ = 0;
		geomechfile << "GMFACE_FRACTURE_CONDUCTIVITY" << endl;
		for (itsetint = pSim->setIdenticalInternalMarker.begin(); itsetint != pSim->setIdenticalInternalMarker.end(); itsetint++, nFractures_++)
		{
			for (int i = 0; i < pSim->nPhysicalFacets; i++)
			{
				if (pSim->vsPhysicalFacet[i].nmark == *itsetint)
					geomechfile << pSim->vsFaceCustom[pSim->vsPhysicalFacet[i].nface].conductivity << endl;
			}
		}
		geomechfile << "/" << endl << endl;

		nFractures_ = 0;
		geomechfile << "GMFACE_FRACTURE_REGION\n";
		for (itsetint = pSim->setIdenticalInternalMarker.begin(); itsetint != pSim->setIdenticalInternalMarker.end(); itsetint++, nFractures_++)
		{
			for (int i = 0; i < pSim->nPhysicalFacets; i++)
			{
				if (pSim->vsPhysicalFacet[i].nmark == *itsetint)
				{
					geomechfile << nFractures_ + 1 << endl;
				}
			}
		}
		geomechfile << "/" << endl << endl;

		nFractures_ = 0;
		geomechfile << "GMFACE_FRACTURE_GROUP\n";
		for (itsetint = pSim->setIdenticalInternalMarker.begin(); itsetint != pSim->setIdenticalInternalMarker.end(); itsetint++, nFractures_++)
		{
			for (int i = 0; i < pSim->nPhysicalFacets; i++)
			{
				if (pSim->vsPhysicalFacet[i].nmark == *itsetint)
				{
					geomechfile << nFractures_ + 1 << endl;
				}
			}
		}
		geomechfile << "/" << endl << endl;
		geomechfile.close();
	}
	geomechfile.close();

	outstring = "fl_pres.txt";
	geomechfile.open(outstring.c_str());

	// fractures first
	geomechfile << "PRESSURE" << endl;
	for (int iface = 0; iface < pSim->nFaces; iface++)
	{
		// internal suface
		if (pSim->vsFaceCustom[iface].nMarker > 0)
		{
			int n1_ = pSim->vsFaceCustom[iface].vNeighbors[0];
			int n2_ = pSim->vsFaceCustom[iface].vNeighbors[1];
			geomechfile << min(pSim->vsCellRockProps[n1_].pressure, pSim->vsCellRockProps[n2_].pressure) << endl;
		}
	}
	// matrix
	for (int ib = 0; ib < pSim->nCells; ++ib)
		geomechfile << pSim->vsCellRockProps[ib].pressure << endl;
	//ghost
	for (int ib = 0; ib < pSim->nCells; ++ib)
	{
		if (pSim->vCellsHasGhost[ib])
			geomechfile << pSim->vsCellRockProps[ib].pressure << endl;
	}

	geomechfile << "\n/\n\n";
	geomechfile.close();

	outstring = "fl_temp.txt";
	geomechfile.open(outstring.c_str());

	// fractures first
	geomechfile << "RTEMP" << endl;
	for (int iface = 0; iface < pSim->nFaces; iface++)
	{
		// internal suface
		if (pSim->vsFaceCustom[iface].nMarker > 0)
		{
			int n1_ = pSim->vsFaceCustom[iface].vNeighbors[0];
			int n2_ = pSim->vsFaceCustom[iface].vNeighbors[1];
			geomechfile << min(pSim->vsCellRockProps[n1_].temp, pSim->vsCellRockProps[n2_].temp) << endl;
		}
	}
	// matrix
	for (int ib = 0; ib < pSim->nCells; ++ib)
		geomechfile << pSim->vsCellRockProps[ib].temp << endl;
	//ghost
	for (int ib = 0; ib < pSim->nCells; ++ib)
	{
		if (pSim->vCellsHasGhost[ib])
			geomechfile << pSim->vsCellRockProps[ib].temp << endl;
	}

	geomechfile << "\n/\n\n";
	geomechfile.close();

	outstring = "fl_zmf.txt";
	geomechfile.open(outstring.c_str());

	geomechfile << "ZMF\n";
	for (int k = 0; k < 5; ++k)
	{
		for (int iface = 0; iface < pSim->nFaces; iface++)
		{
			// internal suface
			if (pSim->vsFaceCustom[iface].nMarker > 0)
			{
				int n1_ = pSim->vsFaceCustom[iface].vNeighbors[0];
				int n2_ = pSim->vsFaceCustom[iface].vNeighbors[1];
				geomechfile << min(pSim->vsCellRockProps[n1_].zmf[k], pSim->vsCellRockProps[n2_].zmf[k]) << endl;
			}
		}

		for (int ib = 0; ib < pSim->nCells; ++ib)
			geomechfile << pSim->vsCellRockProps[ib].zmf[k] << endl;
		//ghost
		for (int ib = 0; ib < pSim->nCells; ++ib)
		{
			if (pSim->vCellsHasGhost[ib])
				geomechfile << pSim->vsCellRockProps[ib].zmf[k] << endl;
		}
	}

	geomechfile << "\n/\n\n";
	geomechfile.close();

	outstring = "gm_fracprops.txt";
	geomechfile.open(outstring.c_str());

	geomechfile << "GMCONTACT_THERMAL_EXPANSION" << endl;
	geomechfile << "0.0" << endl;
	geomechfile << "/" << endl << endl;

	geomechfile << "GMCONTACT_NORMAL_PROPS\n";
	geomechfile << "-0.01	0.041190632086458	1.02471437925187	3.237789418e3	3.237789418e3	2.118894709" << endl;
	geomechfile << "0	0.041190632086458	1.02471437925187	3.237789418	3.237789418	2.118894709" << endl;
	geomechfile << "0	50.3728278581928	31.2236967149157	2.1726124	2.1726124	1.5863062" << endl;
	geomechfile << "0	119.320276113133	72.5921656678797	1.6837264	1.6837264	1.3418632" << endl;
	geomechfile << "0	188.267724368073	113.960634620844	1.4272834	1.4272834	1.2136417" << endl;
	geomechfile << "0	257.215172623013	155.329103573808	1.26216874	1.26216874	1.13108437" << endl;
	geomechfile << "0	326.162620877953	196.697572526772	1.14439306	1.14439306	1.07219653" << endl;
	geomechfile << "0	395.110069132893	238.066041479736	1.05491054	1.05491054	1.02745527" << endl;
	geomechfile << "0	464.057517387833	279.4345104327	0.9833746	0.9833746	0.9916873" << endl;
	geomechfile << "0	533.004965642773	320.802979385664	0.92438524	0.92438524	0.96219262" << endl;
	geomechfile << "0	601.952413897713	362.171448338628	0.87479436	0.87479436	0.93739718" << endl;
	geomechfile << "0	670.899862152653	403.539917291592	0.83220196	0.83220196	0.91610098" << endl;
	geomechfile << "0	739.847310407594	444.908386244556	0.79420804	0.79420804	0.89710402" << endl;
	geomechfile << "0	808.794758662534	486.27685519752	0.76051336	0.76051336	0.88025668" << endl;
	geomechfile << "0	877.742206917474	527.645324150484	0.72991792	0.72991792	0.86495896" << endl;
	geomechfile << "0	939.79491034692	564.876946208152	0.70534172	0.70534172	0.85267086" << endl;

	geomechfile << "/" << endl << endl;
	geomechfile.close();


	outstring = "fl_satnum.txt";
	geomechfile.open(outstring.c_str());

	geomechfile << "SATNUM" << endl;
	for (int iface = 0; iface < pSim->nFaces; iface++)
	{
		if (pSim->vsFaceCustom[iface].nMarker > 0)
			geomechfile << 1 << endl;
	}

	for (int ib = 0; ib < pSim->nCells; ++ib)
		geomechfile << 0 << endl;
	//ghost
	for (int ib = 0; ib < pSim->nCells; ++ib)
	{
		if (pSim->vCellsHasGhost[ib])
			geomechfile << 1 << endl;
	}

	geomechfile << "/" << endl;
	geomechfile.close();


	outstring = "fl_thc.txt";
	geomechfile.open(outstring.c_str());

	geomechfile << "THCROCK" << endl;
	for (int iface = 0; iface < pSim->nFaces; iface++)
	{
		// internal suface
		if (pSim->vsFaceCustom[iface].nMarker > 0)
		{
			int n1_ = pSim->vsFaceCustom[iface].vNeighbors[0];
			int n2_ = pSim->vsFaceCustom[iface].vNeighbors[1];
			geomechfile << min(pSim->vsCellRockProps[n1_].thc, pSim->vsCellRockProps[n2_].thc) << endl;
		}
	}
	for (int ib = 0; ib < pSim->nCells; ++ib)
		geomechfile << pSim->vsCellRockProps[ib].thc << endl;
	// ghost
	for (int ib = 0; ib < pSim->nCells; ++ib)
	{
		if (pSim->vCellsHasGhost[ib])
			geomechfile << pSim->vsCellRockProps[ib].thc << endl;
	}

	geomechfile << "/" << endl;
	geomechfile.close();


	cout << "write all wells\n";
	if (pSim->nWells > 0)
	{
		outstring = "fl_wells.txt";
		geomechfile.open(outstring.c_str());
		geomechfile << "WELSPECS\n";

		for (int iw = 0; iw < pSim->nWells; iw++)
		{
			if (pSim->vsWell[iw].vID.size() > 0)
				geomechfile << "W" << iw << " 1* " << pSim->vsWell[iw].vID[0] + 1 << " 1 " << pSim->vsWell[iw].datum << " WATER /" << endl;
		}

		geomechfile << "/\n\n";

		geomechfile << "COMPDAT\n";
		for (int iw = 0; iw < pSim->nWells; iw++)
		{
			for (int i = 0; i < pSim->vsWell[iw].vID.size(); ++i)
				geomechfile << "W" << iw << "\t" << pSim->vsWell[iw].vID[i] + 1 << " 1 1 1 OPEN 1* " << pSim->vsWell[iw].vWi[i] << " 4* Z/" << endl;
		}
		geomechfile << "/\n\n";
		geomechfile.close();
	}

}
