#pragma once 
#include "Utils/stdafx.h"
#include "MeshBuffer.h"
#include "CGAL/CGAL_Mesh.h"
#include "CGAL/CGAL_Macros.h"
#include <iostream>
using namespace std;

class Parser
{
public:
	Parser() {}
	~Parser() {}
	static bool loadOBJ(const char* filename, MeshBuffer* mb, MeshBuffer* smb);
	static bool loadVectorField(const char* filename, Mesh& cgalMesh);
	static void setMeshAdditionalData(Mesh& cgalMesh, MeshBuffer& m1);
};