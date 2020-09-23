#ifndef MESHFROMFILEHEADERDEF
#define MESHFROMFILEHEADERDEF

#include "Mesh.hpp"
#include <armadillo>
#include <string>

class MeshFromFile : public Mesh
{
public:
	// Constructors for reading in mesh from a file
	MeshFromFile(const std::string fileName);
};

#endif