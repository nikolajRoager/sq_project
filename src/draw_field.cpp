//C++ Text strings and functions for strings
#include<string>
#include<iostream>

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

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


void draw_field(vector<vec> origins,composite_field& Fields, ofstream& OUT, bool B=true)
{
    vector<vector< vec> > out_positions = vector<vector< vec> >(origins.size());
    size_t n_steps_printed = 128;


    //Use a simple ODE to draw the path of the fields
    auto ODE = [&Fields, B]( const boost::array< double, 3 > pos, boost::array< double, 3 >& dpos , const double t )
    {
        //Extract position and velocity from data
        vec v_pos = (vec(pos[0],pos[1],pos[2]));

        vec v_dpos = (B? Fields.get_Bfield(v_pos,t) : Fields.get_Efield(v_pos,t));
        if (v_dpos.x != 0 || v_dpos.y != 0 || v_dpos.z != 0)
            v_dpos = glm::normalize(v_dpos);
        dpos[0] = v_dpos.x;
        dpos[1] = v_dpos.y;
        dpos[2] = v_dpos.z;
    };

    size_t p_print_i =0;
    double T=50;




    for (size_t j = 0; j < origins.size(); ++j)
    {
        out_positions[j]=vector< vec>(n_steps_printed);
        out_positions[j][0]=origins[j];
        auto save_step = [&p_print_i,T, &out_positions, j , &n_steps_printed]( const boost::array< double, 3 >& Data , const double t )
        {

            size_t i = n_steps_printed*t/T;
            if (i==n_steps_printed)
            {
                --i;
            }
            if (i>p_print_i)
            {
                out_positions[j][i] = vec(Data[0],Data[1],Data[2]);
                p_print_i=i;

            }
        };

        p_print_i=0;
        boost::array< double, 3 > POS0 = {origins[j].x,origins[j].y,origins[j].z};

        size_t steps = integrate_const(
            runge_kutta4<boost::array< double, 3 > >(),
            ODE,
            POS0 ,
            0.0 ,
            T ,
            T/(10*n_steps_printed),
            save_step
        );


    }




    for (size_t i = 0; i < n_steps_printed; ++i)
    {
        for (size_t j = 0; j < origins.size(); ++j)
        {

            OUT<<out_positions[j][i].x<<"\t"<<out_positions[j][i].y<<"\t"<<out_positions[j][i].z;
            if (j+1 == origins.size())
                OUT<<endl;
            else
                OUT<<"\t";
        }
    }
}
