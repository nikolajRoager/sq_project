#include"constants.hpp"

#define PI 3.141592653589793

namespace constants
{
    //Names of the units used
    string name_Epotential = "V";
    string name_distance = "m";
    string name_time = "s";
    string name_charge = "C";
    string name_mass = "kg";
    string name_force = "N";

    void set_name_Epotential(string N){name_Epotential = N;}
    void set_name_distance(string N)  {name_distance = N;}
    void set_name_charge(string N)    {name_charge = N;}
    void set_name_time(string N)      {name_time = N;}
    void set_name_mass(string N)      {name_mass = N;}
    void set_name_force(string N)      {name_force = N;}

    string get_name_current() {return "("+name_charge+" "+name_time+"^{-1})";}
    string get_name_B_field() {return get_name_current()+" "+name_distance+"^{-1}";}
    string get_name_E_field()  {return name_Epotential+" "+name_distance+"^{-1}";}
    string get_name_mass()     {return name_mass;}
    string get_name_time()     {return name_time;}
    string get_name_inv_time()     {return name_time+"^{-1}";}
    string get_name_distance() {return name_distance;}
    string get_name_charge()   {return name_charge;}
    string get_name_force()   {return name_force;}

    //Default, SI units
    double mu0 = 4*PI*1e-7;
    void set_mu0(double val)
    {

        mu0=val;
    }
    double get_mu0(){return mu0;}


    double eps0 =  8.8541878128e-12;
    void set_eps0(double val)
    {
        eps0=val;
    }
    double get_eps0(){return eps0;}
}
