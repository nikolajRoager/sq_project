#pragma once
/*This file
Some physical constants
*/

#include<string>

using namespace std;

namespace constants
{

    void set_name_Epotential(string N);
    void set_name_distance(string N);
    void set_name_charge(string N);
    void set_name_time(string N);
    void set_name_mass(string N);
    void set_name_force(string N);

    string get_name_B_field();
    string get_name_E_field();
    string get_name_mass();
    string get_name_distance();
    string get_name_time();
    string get_name_inv_time();
    string get_name_charge();
    string get_name_force();
    string get_name_current();

    void set_mu0(double val);
    double get_mu0();

    void set_eps0(double val);
    double get_eps0();
}
