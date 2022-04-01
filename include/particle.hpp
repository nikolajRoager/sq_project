#pragma once

/*THIS FILE

A charged particle

*/

#include<fstream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>



//GLM, openGL Mathematics library, originally created to do the matrix and vector calculations needed for OpenGL powered 3D programs, it is also a decent vector and matrix mathematics library.
#include<glm/glm.hpp>


using vec = glm::dvec3;//Default to glm::douple precision 3D vector, by default glm uses single precision for comptatibility with the GLSL shading language (Single precision is good enough for computer-graphics application where speed is more important than accuracy)

typedef boost::array< double, 6 > state_type;//(position,velocity), need to be as doubles as odeint need to use > operator


#include"field.hpp"
#include"constants.hpp"

//A particle with a single mass, single charge and some position and momentum, as of yet, not build in magnetic dipole moment (though that would be cool)

//In this case, the default copy and move constructors are both good enough, and not required
class particle
{
private:
    double mass;
    double inv_mass;//1/mass, I here some people insist that divisions are evil and must be avoided at all cost, but my tests indicate that the performance is not *that* much worse than multiplication, still, might as well avoid divisions when I can
    vector<vec> positions;//In c++ context, vector means a dynamic sized list, in this case a list of position vectors, this does not need to include every single point calculated, just a representative sample.
    vector<vec> velocities;
    vector<double> timestamps;
    double charge;

    vec v0;
    vec pos0;

    string name;

    static double print_interval;//how far (in simulation time units) is there between subsequent prints of data

    double p_print=0;//When was data last saved last print
    //Compatible with boost/odeint, save the current data




public:

    particle(vec pos0, vec v0, double m, double q, string name);

    //Run the entire simulation for this particle from t0=0, using odeint library (RK45)
    void calculate(const composite_field& Fields, double T, double dt);

    void calculate_euler(const composite_field& Fields, double T, double dt);

    void calculate_RK4(const composite_field& Fields, double T, double dt);

    void calculate_RKDP45(const composite_field& Fields, double T, double dt);

    void calculate_3rdparty_RK4(const composite_field& Fields, double T, double dt);

    void analytical_solenoid(const composite_field& Fields, double T, double dt);

    //For saving, return the ID of the particle position/velocity/time vectors
    const char* positionID() const;
    const char* velocityID() const;
    const char* timeID() const;

    const string& getName() const {return name;}

    static void set_print_interval(double interval);
    void dump_raw(std::ofstream& txt_output) const;
    void write_final(std::ofstream& txt_output) const;
    size_t get_data_size() const {return positions.size();}
};
