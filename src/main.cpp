/*THIS FILE

main function, main simulation loop

*/

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
#include"draw_field.hpp"

//This is too long to write explicitly, shorten these down a little
using namespace std;
namespace fs = std::filesystem;

//I am not quite sure what type of integer makes most sense to use, so lets just define a signed int type here



//The main loop, I try to keep anything graphics related out of this.
int main(int argc, char* argv[])
{
    if (argc!=3)
    {
        cout<<"Usage "<<argv[0]<<" setupfile.txt outpath"<<endl;
        return 1;
    }

    fs::path outpath = fs::path(string(argv[2]));

    if (!fs::exists(outpath))
    {
        cout<<"Creating output directory"<<endl;
        fs::create_directory(outpath);
    }
    else if (!fs::is_directory(outpath))
    {
        cout<<argv[2]<<" is a file, must be a directory or not exist"<<endl;
        return 1;
    }


    bool alt_save_B_field = false;
    bool alt_save_E_field = false;
    vector<vec> alt_save_B_pos;
    vector<vec> alt_save_E_pos;
    string alt_save_B_name;
    string alt_save_E_name;

    bool save_B_field=false;
    bool save_E_field=false;
    double save_field_min;
    double save_field_max;
    size_t save_field_steps;

    composite_field Field;
    vector<particle> particles;

    double T=0;
    vector<double> T_custom;
    double dt=0;
    bool has_extrafile=false;
    string Extra_file="";
    bool txt_output=false;

    uint8_t engine_type = 0;//0=Default boost::numeric::odeint adaptive stepsize RK45, 1=Euler fixed step, 2=RK4 fixed step 3, RK4 from build in library

    //Load setup information in this limited scope
    {
        cout<<"Trying to open setup file"<<endl;
        fs::path setup_path = fs::path(argv[1]);
        ifstream setup_file(setup_path );

        if (!setup_file.is_open())
        {
            cout<<"Setup file "<<setup_path.string()<<" could not be opened"<<endl;
            return 1;
        }


        //Simple, and not very safe plain-text setup loading, I assume the correct format is used and this may crash horribly if it is not. It should be pretty clear what happens here
        string this_line;
        uint16_t line_count=0;


        while (getline(setup_file,this_line))//Stream the setup file, line by line
        {

            stringstream ss(this_line);//now stream each of the elements on the line

            string input;

            //Read the first command on this line (if the line is not blank) this eats any blank space or tabs before the command
            if (ss>>input)
            {
                if(input[0]=='#')
                {
                    //Comment, skip
                }
                else if (input.compare("engine")==0)
                {
                    if (!(ss>>input))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }
                    if (input.compare("RK45lib")==0)
                        engine_type =0;
                    else if (input.compare("euler")==0)
                        engine_type =1;
                    else if (input.compare("RK4")==0)
                        engine_type =2;
                    else if (input.compare("RK4lib")==0)
                        engine_type =3;
                    else if (input.compare("RK45")==0)
                        engine_type =4;
                    else
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" not valid engine"<<endl;
                        cout<<"\""<<input<<"\""<<endl;

                        return 1;
                    }

                }
                else if (input.compare("distance_unit")==0)
                {
                    if (!(ss>>input))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                    }

                    constants::set_name_distance(input);
                }
                else if (input.compare("time_unit")==0)
                {
                    if (!(ss>>input))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                    }

                    constants::set_name_time(input);
                }
                else if (input.compare("charge_unit")==0)
                {
                    if (!(ss>>input))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                    }

                    constants::set_name_charge(input);
                }
                else if (input.compare("mass_unit")==0)
                {
                    if (!(ss>>input))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                    }

                    constants::set_name_mass(input);
                }
                else if (input.compare("force_unit")==0)
                {
                    if (!(ss>>input))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                    }

                    constants::set_name_force(input);
                }
                else if (input.compare("Epotential_unit")==0)
                {
                    if (!(ss>>input))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                    }

                    constants::set_name_Epotential(input);
                }
                else if (input.compare("eps0")==0)
                {
                    double eps0;
                    if (!(ss>>eps0))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                    }


                    constants::set_eps0(eps0);
                }
                else if (input.compare("mu0")==0)
                {
                    double mu0;
                    if (!(ss>>mu0))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                    }


                    constants::set_mu0(mu0);
                }
                else if (input.compare("T")==0)
                {
                    if (!(ss>>T))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                    }
                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting argument to floating point number"<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }
                }
                else if (input.compare("save_dt")==0)
                {
                    double Dt=0.1;
                    if (!(ss>>Dt))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting argument to float"<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }

                    particle::set_print_interval(Dt);
                }
                else if (input.compare("dt")==0)
                {
                    if (!(ss>>dt))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting argument to float"<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }

                    if (dt<0)
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error"<<endl;
                        cout<<"Backwards time travel is not allowed"<<endl;
                        return 1;
                    }

                }
                else if (input.compare("field")==0)
                {
                    bool B = false;

                    if(!(ss>>input))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 1st argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                    if (input.compare("B")==0)
                        B=true;
                    else if (input.compare("E")==0)
                        B=false;
                    else
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error"<<endl;
                        cout<<"Unknown field, must be B or E, got "<<input<<endl;
                    }


                    if(!(ss>>input))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 2nd argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }
                    if (input.compare("constant")==0)
                    {
                        vec strength;
                        if(!(ss>>strength.x))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 3rd argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 3rd argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>strength.y))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 4th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 4th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>strength.z))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 5th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 5th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }

                        Field+=new constant_field(B,strength);

                    }
                    else if (input.compare("Earth")==0)
                    {
                        double B0, Re;
                        if(!(ss>>Re))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 3rd argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 3rd argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }

                        if(!(ss>>B0))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 4th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 4th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }


                        Field+=new Earth_dipole2(Re,B0);

                    }
                    else if (input.compare("dipole")==0)
                    {
                        double mu0per4pi;
                        vec m;
                        if(!(ss>>mu0per4pi))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 3rd argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 3rd argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>m.x))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 4th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 4th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>m.y))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 5th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 5th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }

                        if(!(ss>>m.z))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 6th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 6th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }

                        Field+=new dipole(m,mu0per4pi);

                    }
                    else if (input.compare("solenoid")==0)
                    {
                        vec origin;
                        vec dir;
                        unsigned int N;
                        double I;
                        double R;
                        if(!(ss>>origin.x))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 3rd argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 3rd argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>origin.y))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 4th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 4th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>origin.z))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 5th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 5th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>dir.x))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 6th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 6th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>dir.y))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 7th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 7th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }

                        if(!(ss>>dir.z))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 8th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 8th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if(!(ss>>N))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 9th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 9th argument to unsigned integer"<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }

                        if(!(ss>>I))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 10th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 10th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if(!(ss>>R))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 11th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 11th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        Field+=new solenoid(origin,dir,N,I,R);

                    }
                    else if (input.compare("solenoid1")==0)
                    {
                        vec origin;
                        vec dir;
                        double B,R;
                        if(!(ss>>origin.x))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 3rd argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 3rd argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>origin.y))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 4th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 4th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>origin.z))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 5th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 5th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>dir.x))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 6th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 6th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>dir.y))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 7th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 7th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }

                        if(!(ss>>dir.z))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 8th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 8th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if(!(ss>>B))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 9th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 9th argument to floating point number"<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }


                        if(!(ss>>R))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 10th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 10th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        Field+=new solenoid(origin,dir,B,R);

                    }
                    else if (input.compare("cyclotron_gab")==0)
                    {
                        double R,gab,B,omega,phase;
                        if(!(ss>>R))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 3rd argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 3rd argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>gab))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 4th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 4th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>B))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 5th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 5th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>omega))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 5th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 5th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>phase))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 5th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 5th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        Field+=new cyclotron_gab(R,gab,B,omega,phase);
                    }
                    else if (input.compare("torus")==0)
                    {
                        double B,R,r;
                        if(!(ss>>R))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 3rd argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 3rd argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>r))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 4th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 4th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        if(!(ss>>B))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 5th argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                        if (ss.fail())
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 5th argument to floating point number "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;

                        }
                        Field+=new torus(R,r,B);
                    }
                    else
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error "<<endl;
                        cout<<"Unsupported field type"<<endl;
                        return 1;
                    }


                }
                else if (input.compare("printfield")==0)
                {
                    if(!(ss>>input))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }
                    if (input.compare("B")==0)
                        save_B_field=true;
                    else if (input.compare("E")==0)
                        save_E_field=true;
                    else
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" argument should be B or E"<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                }
                else if (input.compare("save_field")==0)
                {
                    if(!(ss>>input))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }
                    if (input.compare("B")==0)
                        alt_save_B_field=true;
                    else if (input.compare("E")==0)
                        alt_save_E_field=true;
                    else
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" argument should be B or E"<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                    if(alt_save_B_field)
                    {
                        if(!(ss>>alt_save_B_name))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 2nd argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }
                    }
                    else
                        if(!(ss>>alt_save_E_name))
                        {
                            cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 2nd argument "<<endl;
                            cout<<"\""<<this_line<<"\""<<endl;
                            return 1;
                        }

                    vec V;
                    while (ss>>V.x && ss>>V.y && ss>>V.z)
                    {
                        if(alt_save_B_field)
                            alt_save_B_pos.push_back(V);
                        else
                            alt_save_E_pos.push_back(V);

                    }


                }
                else if (input.compare("field_res")==0)
                {

                    if(!(ss>>save_field_min))
                    {

                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 1st argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 1st argument to int "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }
                    else if(!(ss>>save_field_max))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 2nd argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }
                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 2nd argument to int "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }
                    else if(!(ss>>save_field_steps))
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 3rd argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }
                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 3rd argument to int "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }
                    if (save_field_steps<=0)//I could have used an unsigned integer here, but then I would have to see a lot of compiler warnings, which I don't want
                    {

                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" 3rd argument must be positive"<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }




                }
                else if (input.compare("extra")==0)//A file with extra graphical setup, will be copied to the output
                {
                    ss>>Extra_file;

                    if (!fs::exists(Extra_file))
                    {

                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" File not found "<<Extra_file<<endl;
                        return 1;
                    }
                    has_extrafile=true;
                }
                else if (input.compare("write_txt")==0)
                {
                    txt_output=true;
                }
                else if (input.compare("particle")==0)//Spawn a particle and make it move
                {

                    double mass, q;
                    vec pos,V;
                    if(!(ss>>mass))
                    {

                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 1st argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 1st argument to floating point number"<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }

                    if(!(ss>>q))
                    {

                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 2nd argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 2nd argument to floating point number "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }
                    if(!(ss>>pos.x))
                    {

                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 3rd argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 3rd argument to floating point number "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }
                    if(!(ss>>pos.y))
                    {

                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 4th argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 4th argument to floating point number "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }
                    if(!(ss>>pos.z))
                    {

                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 5th argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 5th argument to floating point number "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }
                    if(!(ss>>V.x))
                    {

                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 6th argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 6th argument to floating point number "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }
                    if(!(ss>>V.y))
                    {

                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 7th argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 7th argument to floating point number "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }
                    if(!(ss>>V.z))
                    {

                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" missing 8th argument "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;
                    }

                    if (ss.fail())
                    {
                        cout<<'('<<setup_path.string()<<") line "<<line_count<<" error on casting 8th argument to floating point number "<<endl;
                        cout<<"\""<<this_line<<"\""<<endl;
                        return 1;

                    }

                    //Optionally, this particle can be simulated for a different time
                    double thisT = T;


                    if (!(ss >> thisT))
                        thisT=T;//Default time

                    T_custom.push_back (thisT);


                    //Optinally pick a different name
                    string name_input;
                    string name="";

                    while(ss>>name_input)
                    {
                        name = name+name_input+" ";
                    }

                    if (name.length()==0)
                        name = "particle "+to_string(particles.size());

                    particles.push_back(particle(pos,V,mass,q,name));
                }
                else
                {
                    cout<<'('<<setup_path.string()<<") line "<<line_count<<" Not valid "<<endl;
                    cout<<"\""<<this_line<<"\""<<endl;
                    return 1;
                }


            }//Ends, read the first word in every line



            ++line_count;//This makes error messages actually helpful
        }//Ends, loop all lines in the setup file

        setup_file.close();
        cout<<"setup loaded"<<endl;
    }//Ends limited scope



    //Check if the setup is sane
    if (save_field_min>save_field_max)
    {
        cout<<"Can't save field from "<<save_field_min<<" to "<<save_field_max<<endl;
        return 1;
    }

    //Tell us what we run
    cout<<"===Simulation setup==="<<endl;
    cout<<"From t=0 to (default) t="<<T<<' '<<constants::get_name_time()<<endl;
    cout<<"Δt="<<dt<<' '<<constants::get_name_time()<<endl;

    cout<<"----Fields used----"<<endl;
    Field.log_info();

    //Before we run, lets save the electrostatic setup, it does not change
    if (save_B_field || save_B_field )
    {
        //std::array must have output size set at compiletime, hence I use a vector here, even though I do not intend to resize it
        //I want to save this in whatever format is easiest and fastest to load in python afterwards, and matplotlib vector-plot (quiveR) wants the coordinates separate
        vector<double> output_data_Bx;
        vector<double> output_data_By;
        vector<double> output_data_Bz;

        vector<double> output_data_Ex;
        vector<double> output_data_Ey;
        vector<double> output_data_Ez;


        size_t save_field_steps2=save_field_steps*save_field_steps;
        size_t save_field_steps3 =save_field_steps2*save_field_steps;

        //the numpy meshgrid function is bad, I can not get it to line up with the data both in the 3D and 2D case, so, I simply save the grid location for display
        vector<double> output_data_refX=vector<double>(save_field_steps3);
        vector<double> output_data_refY=vector<double>(save_field_steps3);
        vector<double> output_data_refZ=vector<double>(save_field_steps3);


        if (save_B_field)
        {
            cout<<"B";
            output_data_Bx=vector<double>(save_field_steps3);
            output_data_By=vector<double>(save_field_steps3);
            output_data_Bz=vector<double>(save_field_steps3);
            if (save_E_field)
                cout <<" & ";

        }
        if (save_E_field)
        {
            output_data_Ex=vector<double>(save_field_steps3);
            output_data_Ey=vector<double>(save_field_steps3);
            output_data_Ez=vector<double>(save_field_steps3);
            cout<<"E";
        }
        cout<<" field values "<<endl;

        //Loop and save everything
        double D = (save_field_max-save_field_min)/double(save_field_steps-1);
        for (size_t k = 0; k<save_field_steps; ++k)
            for (size_t j = 0; j<save_field_steps; ++j)
                for (size_t i = 0; i<save_field_steps; ++i)
                {

                    vec pos = vec(D*i+save_field_min,D*j+save_field_min,D*k+save_field_min);


                    //Pythons meshgrid function is bad, I can not get it to line up with the data I save in both the 2D and 3D plot at the same time, So I just save the position to be used in the plot outright
                    output_data_refX[k+j*save_field_steps+i*save_field_steps2]=pos.x;
                    output_data_refY[k+j*save_field_steps+i*save_field_steps2]=pos.y;
                    output_data_refZ[k+j*save_field_steps+i*save_field_steps2]=pos.z;

                    //Save this in whatever order is most optimal for us to load it in python (we need to cast the x,y,z indices to wherever they would end up with the 3D meshgrid function
                    if (save_B_field)
                    {
                        vec B = Field.get_Bfield(pos);
                        output_data_Bx[k+j*save_field_steps+i*save_field_steps2]=B.x;
                        output_data_By[k+j*save_field_steps+i*save_field_steps2]=B.y;
                        output_data_Bz[k+j*save_field_steps+i*save_field_steps2]=B.z;
                    }
                    if (save_E_field)
                    {
                        vec E = Field.get_Efield(pos);
                        output_data_Ex[k+j*save_field_steps+i*save_field_steps2]=E.x;
                        output_data_Ey[k+j*save_field_steps+i*save_field_steps2]=E.y;
                        output_data_Ez[k+j*save_field_steps+i*save_field_steps2]=E.z;

                    }
                }

        cout<<endl;


        ofstream B_file(outpath/"pos_ref.bin" ,std::ios::binary );
        if (B_file.is_open())
        {
            //I save this in the order it is most practical to read this in python
            B_file.write((const char*) &(output_data_refX[0]),sizeof(double)*save_field_steps3);
            B_file.write((const char*) &(output_data_refY[0]),sizeof(double)*save_field_steps3);
            B_file.write((const char*) &(output_data_refZ[0]),sizeof(double)*save_field_steps3);
            B_file.close();
        }
        else
        {
            save_B_field=false;
            cout<<"Grid position reference could not be saved, file "<<(outpath/"pos_ref.bin").string()<<" could not be opened"<<endl;//Just carry on, maybe the rest of the things can be saved
        }

        if (save_B_field)
        {
            ofstream B_file(outpath/"B.bin" ,std::ios::binary );

            if (B_file.is_open())
            {
                //I save this in the order it is most practical to read this in python
                B_file.write((const char*) &(output_data_Bx[0]),sizeof(double)*save_field_steps3);
                B_file.write((const char*) &(output_data_By[0]),sizeof(double)*save_field_steps3);
                B_file.write((const char*) &(output_data_Bz[0]),sizeof(double)*save_field_steps3);
                B_file.close();
            }
            else
            {
                save_B_field=false;
                cout<<"B field could not be saved, file "<<(outpath/"B.bin").string()<<" could not be opened"<<endl;//Just carry on, maybe the rest of the things can be saved
            }
        }
        if (save_E_field)
        {
            ofstream E_file(outpath/"E.bin" ,std::ios::binary );

            if (E_file.is_open())
            {
                //Save in the order it makes most sense to load in
                E_file.write((const char*) &(output_data_Ex[0]),sizeof(double)*save_field_steps3);
                E_file.write((const char*) &(output_data_Ey[0]),sizeof(double)*save_field_steps3);
                E_file.write((const char*) &(output_data_Ez[0]),sizeof(double)*save_field_steps3);
                E_file.close();
            }
            else
            {
                save_E_field=false;
                cout<<"E field could not be saved, file "<<(outpath/"E.bin").string()<<" could not be opened"<<endl;//Just carry on, maybe the rest of the things can be saved
            }

        }

    }

    //If we need to use the alternative field saving function
    if (alt_save_B_field)
    {
        ofstream B_file(outpath/alt_save_B_name);
        if (B_file.is_open())
        {

            draw_field(alt_save_B_pos,Field, B_file, true);
            B_file.close();
        }
    }

    if (alt_save_E_field)
    {
        ofstream E_file(outpath/alt_save_E_name);
        if (E_file.is_open())
        {

            draw_field(alt_save_E_pos,Field, E_file, false);
            E_file.close();
        }

    }

    cout<<"Running simulation"<<endl;
    //Run the simulation of all particles



    for (size_t i = 0; i< particles.size(); ++i)
    {
        cout<<(i+1)<<"/"<<particles.size()<<"..."<<flush;//calculate prints how many steps were used and include endl;
        particle& P = particles[i];

        if (engine_type ==0)
            P.calculate(Field,T_custom[i],dt);
        else if (engine_type ==1)
            P.calculate_euler(Field,T_custom[i],dt);
        else if (engine_type ==2)
            P.calculate_RK4(Field,T_custom[i],dt);
        else if (engine_type ==3)
            P.calculate_3rdparty_RK4(Field,T_custom[i],dt);
        else if (engine_type ==4)
            P.calculate_RKDP45(Field,T_custom[i],dt);

        //Time crunch analytical solenoid, sorry
        //P.analytical_solenoid(Field,T_custom[i],dt);
    }

    //Save a template file telling the python display program what to load, csv is easy to do without requiring a third party csv library
    ofstream data_setup_file(outpath/"setup.csv" );

    if (!(data_setup_file.is_open()))
    {
        cout<<"Could not open "<<(outpath/"setup.csv").string()<<endl;
        return 1;
    }
    if (save_E_field)
        data_setup_file<<"display_E_field , True"<<endl;
    else
        data_setup_file<<"display_E_field , False"<<endl;

    if (save_B_field)
        data_setup_file<<"display_B_field , True"<<endl;
    else
        data_setup_file<<"display_B_field , False"<<endl;


    data_setup_file<<"field_min , "<<save_field_min<<endl;
    data_setup_file<<"field_max , "<<save_field_max<<endl;
    data_setup_file<<"field_steps , "<<save_field_steps<<endl;
    data_setup_file<<"T, "<<T<<endl;
    data_setup_file<<"dt, "<<dt<<endl;

//I needed to get some data asap, so I hardcoded this, YES HARDCODED, I know it is really bad, but sorry, uncomment this to only do one particle but get NSTEPS timestep sizes
//#define CRUNCH_DT

#ifndef CRUNCH_DT
    data_setup_file<<"particles , "<<particles.size()<<endl;
    for (size_t i = 0; i <particles.size(); ++i)
    {
        data_setup_file<<"name"<<" , "<<particles[i].getName()<<endl;
        data_setup_file<<"timesteps , "<<particles[i].get_data_size()<<endl;
    }
#else

uint8_t NSTEPS = 128;

//Time crunch hardcoded cheat for testing many dt
    data_setup_file<<"particles , "<<NSTEPS<<endl;//Yep, I am in that much of a hurry
    for (size_t i = 0; i <NSTEPS ; ++i)
    {
        data_setup_file<<"name"<<" , h="<<dt*i<<constants::get_name_time()<<endl;
        data_setup_file<<"timesteps , "<<particles[0].get_data_size()<<endl;
    }
#endif

    //Now save all the particles
    data_setup_file.close();

    cout<<"Saving"<<endl;


#ifndef CRUNCH_DT
    for (unsigned short i = 0; i < particles.size(); ++i )
    {

    unsigned short j =0;
    particle& P = particles[i];
#else

    ofstream final_data (outpath/("Final.tsv"));
//Time crunch hardcoded cheat for testing many dt, change i to j
    double dt_step = dt;
    for (unsigned short j = 0; j < NSTEPS ; ++j )
    {
        dt+=dt_step;
        unsigned short i =0;

        //This did not fit into the way it was set up, I wanted to keep it separate, but if I am just using one particle, this is the fastest way to implement this
        cout<<"Rerun"<<(j+1)<<"/"<<NSTEPS<<"..."<<flush;//calculate prints how many steps were used and include endl;
        particle& P = particles[i];


        if (engine_type ==0)
            P.calculate(Field,T_custom[i],dt);
        else if (engine_type ==1)
            P.calculate_euler(Field,T_custom[i],dt);
        else if (engine_type ==2)
            P.calculate_RK4(Field,T_custom[i],dt);
        else if (engine_type ==3)
            P.calculate_3rdparty_RK4(Field,T_custom[i],dt);


#endif


        ofstream particle_path_file(outpath/("particles"+to_string(j+i)+".bin"),std::ios::binary );

        if (!(particle_path_file.is_open()))
        {
            cout<<"Could not open "<<(outpath/("particles"+to_string(j+i)+".bin")).string()<<endl;
            return 1;
        }


        particle_path_file.write(P.positionID(),sizeof(vec)*P.get_data_size());


        particle_path_file.close();

        //raw text output for gnuplot (which does not handle binary well)
        if (txt_output)
        {
            ofstream txt_file(outpath/("particles"+to_string(j+i)+".tsv"));

            P.dump_raw(txt_file);

            txt_file.close();


        }

#ifdef CRUNCH_DT
        final_data <<dt<<'\t';
        P.write_final(final_data );
#endif

    }

#ifdef CRUNCH_DT
    final_data.close();
#endif
    if(has_extrafile)
    {
        fs::copy(Extra_file,outpath/"extra.txt",fs::copy_options::update_existing);
    }

    cout<<"Program quits"<<endl;
    return 0;
}
