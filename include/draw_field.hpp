//C++ Text strings and functions for strings
#include<string>
#include<iostream>

//Cross platform loading and saving of files, compatible with all operating systems!
#include<filesystem>
#include<fstream>

//Fixed size integers
#include<cstdint>

//GLM, openGL Mathematics library, originally created to do the matrix and vector calculations needed for OpenGL powered 3D programs, it is also a decent vector and matrix mathematics library.
#include<glm/glm.hpp>


using vec = glm::dvec3;//Default to glm::douple precision 3D vector, by default glm uses single precision for comptatibility with the GLSL shading language (Single precision is good enough for computer-graphics application where speed is more important than accuracy)

#include"field.hpp"
#include"constants.hpp"
#include"particle.hpp"

using namespace std;

void draw_field(vector<vec> origins,composite_field& Fields, ofstream& OUT, bool B=true);
