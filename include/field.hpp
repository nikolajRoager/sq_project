
#pragma once
/*
THIS FILE

A field of some sort, this class is polymorphic and can be specified as a dipole, constant field or similar

*/


//GLM, openGL Mathematics library, originally created to do the matrix and vector calculations needed for OpenGL powered 3D programs, it is also a decent vector and matrix mathematics librar.
#include<glm/glm.hpp>
#include<vector>



using vec = glm::dvec3;//Default to glm::douple precision 3D vector


using namespace std;

class field
{
protected:
    bool is_B_field=true;//Is this a B or E-field
public:
    field(bool B){is_B_field=B;}
    virtual ~field(){};
    virtual vec get_field(vec pos, double t=0) const=0;//Get the field strength and direction at this amplitude, this is a pure virtual class

    //For debugging, print the information about the field to the console
    virtual void log_info() const=0;
    virtual bool is_B() const {return is_B_field;}
};

class constant_field : public field
{
private:
    vec my_field;

public:

    constant_field (bool B,vec F);
    ~constant_field(){};
    virtual vec get_field(vec pos, double t=0) const;
    virtual void log_info() const;
};



//Dipole approximation of the Earth magnetosphere, 2nd attempt
class Earth_dipole2 : public field
{
private:
    double Re;
    double B0;

    vec North;

    double factor;//All the constants involved in a dipole field, combined with the dipole moment size to get B0 at Re at the pole

public:

    Earth_dipole2(double _Re, double _B0):field(true)
    {
        Re=_Re;
        B0=_B0;


        //Need size of this to be B0
        //(factor/(Re*Re*Re))*(3.0*(glm::dot(North,nr)*nr)-North);//Equation 5.89 in Griffiths, Introduction to electrodynamics

        //B0 = factor*|1/(Re*Re*Re))|*|(3.0*(glm::dot(North,nr)*nr)-North)|;
        North=vec(0,0,1);//Magnetic north, really the south pole
        //B0 =factor* 2.0*|North|/(Re*Re*Re);
        factor = B0*(Re*Re*Re)/2.0;

    }
    ~Earth_dipole2(){};

    virtual vec get_field(vec pos, double t=0) const;

    //The dipole fields are, as of yet, hardcoded magnetic fields
    virtual bool is_B() const {return true;}
    virtual void log_info() const;
};



//Dipole approximation of the Earth magnetosphere
class Earth_dipole : public field
{
private:
    double Re;
    double B0;

    vec North;

public:

    Earth_dipole(double _Re, double _B0):field(true)
    {
        Re=_Re;
        B0=_B0;
        North=vec(0,0,1);
    }
    ~Earth_dipole(){};

    virtual vec get_field(vec pos, double t=0) const;

    //The dipole fields are, as of yet, hardcoded magnetic fields
    virtual bool is_B() const {return true;}
    virtual void log_info() const;
};




class dipole : public field
{
private:
    vec dipole_moment;
    double factor;//Pre calculate the factor mu0/4pi to work in whatever units are used internally

public:

    dipole(vec m,double mu0Per4pi=1e-6):field(true)//Vacuum permeability in SI units, ok, it is not exactly μ0 = 4π*10^-6 N/A^2 after some recent reshuffling of the SI system, but for all practical purposes it is.
    {
        dipole_moment = m;
        factor = mu0Per4pi;
    }
    ~dipole(){};

    virtual vec get_field(vec pos, double t=0) const;

    //The dipole fields are, as of yet, hardcoded magnetic fields
    virtual bool is_B() const {return true;}
    virtual void log_info() const;
};


class solenoid: public field
{
private:
    vec dir;
    vec origin;
    double R;

    double I;//Use negative if the winding is the other way around
    unsigned int N;
    vec inside_field;

public:

    solenoid(vec O,vec D, unsigned int _N,double _I,double R);
    solenoid(vec O,vec D, double B, double R);//Skip calculations, just give the field
    ~solenoid(){};

    virtual vec get_field(vec pos, double t=0) const;

    virtual bool is_B() const {return true;}
    virtual void log_info() const;
};


class torus: public field
{
private:
    double Rminor;
    double Rmajor;
    double center_field;

public:

    torus(double R, double r, double B);
    ~torus(){};

    virtual vec get_field(vec pos,double t=0) const;

    virtual bool is_B() const {return true;}
    virtual void log_info() const;
};



//Paul trap, simplified
class Paul_trap: public field
{
private:
    double R;//Distance from center to charged "rods"
    double omega;//Frequency

    double Qmax;
    double kQmax;

public:


    Paul_trap(double R, double omega, double Qmax);
    ~Paul_trap(){};

    virtual vec get_field(vec pos,double t=0) const;

    virtual bool is_B() const {return false;}
    virtual void log_info() const;

};


//The gab between the plates of a cyclotron, the field varies sinusoidally
class cyclotron_gab: public field
{
private:
    double R;//The plates are typically D-shaped, what is the max radius after which the field will just be pictured as stopping
    double omega;//Frequency
    double phase_offset;

    //The gab is always placed around the x-axis, it is this wide
    double gab_size;
    double half_gab_size;
    double Emax;

public:


    cyclotron_gab(double R, double gab_size, double Emax, double omega, double phase_offset=0);
    ~cyclotron_gab(){};

    virtual vec get_field(vec pos,double t=0) const;

    virtual bool is_B() const {return false;}
    virtual void log_info() const;

};



//A class for handling multiple E and B fields added together
class composite_field
{
private:
    //As field is polymorphic, this will include both B an E fields, and may include multiple different types
    vector<field*> fields;

public:

    //A little sloppy design, I assume whoever called this made sure the fields are not nullptr
    composite_field(vector<field*> Fields = vector<field*>())
    {
        fields = Fields;
    }
    ~composite_field()
    {
        //Delete all the field pointers from memory
        for (field* F : fields)
            if (F!=nullptr)
                delete F;
    }

    void operator+=(field* new_field)
    {
        if (new_field !=nullptr)
            fields.push_back(new_field);
    }

    vec get_Bfield(vec pos, double t=0) const;

    vec get_Efield(vec pos, double t=0) const;
    virtual void log_info() const;
};
