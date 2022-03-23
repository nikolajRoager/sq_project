#include"field.hpp"
#include"constants.hpp"

#include<iostream>
#define PI 3.141592653589793

using namespace std;

constant_field::constant_field (bool B,vec F):field(B)
{
    my_field = F;
}

vec constant_field::get_field(vec pos,double t) const
{
    return my_field;
}
void constant_field::log_info() const
{
    cout<<"Constant "<<(is_B_field ? "B" : "E")<<" field ["<<my_field.x<<' '<<my_field.y<<' '<<my_field.z<< "] " <<(is_B_field  ? constants::get_name_B_field() : constants::get_name_E_field())<<endl;
}

vec dipole::get_field(vec pos,double t) const
{
    double r = sqrt(glm::dot(pos,pos));
    if (r==0)
        return vec(0);//Will break at origin, just ignore this

    vec nr = glm::normalize(pos);
    return (factor/(r*r*r))*(3.0*(glm::dot(dipole_moment,nr)*nr)-dipole_moment);//Equation 5.89 in Griffiths, Introduction to electrodynamics
}

void dipole::log_info() const
{
    cout<<"Magnetic dipole at origin with magnetic moment vector ["<<dipole_moment.x<<' '<<dipole_moment.y<<' '<<dipole_moment.z<< "] " <<constants::get_name_charge()<<" "<<constants::get_name_time()<<"^-1 "<<constants::get_name_distance()<<"^2 " <<endl;
}

vec Earth_dipole2::get_field(vec pos,double t) const
{
    double r = sqrt(glm::dot(pos,pos));
    if (r==0)
        return vec(0);//Will break at origin, just ignore this

    vec nr = glm::normalize(pos);
    return (factor/(r*r*r))*(3.0*(glm::dot(North,nr)*nr)-North);//Equation 5.89 in Griffiths, Introduction to electrodynamics
}

void Earth_dipole2::log_info() const
{
    cout<<"Magnetic dipole Earth approximation at origin with factor "<<factor <<endl;
    cout<<"Check  "<<get_field(vec(0,0,Re)).z<<" "<<get_field(vec(0,0,10*Re)).z<<endl;

    cout<<"Check  "<<get_field(vec(Re,0,0)).z<<" "<<get_field(vec(Re*10,0,0)).z<<endl;
}


vec Earth_dipole::get_field(vec pos,double t) const
{
    double r = sqrt(glm::dot(pos,pos));
    if (r==0)
        return vec(0);

    vec nr = glm::normalize(pos);

    double costheta= nr.z;
    double sintheta = sqrt(nr.x*nr.x+nr.y*nr.y);

    double Br = -2*B0* pow(Re/r,3)*costheta;
    double Btheta =  -B0* pow(Re/r,3)*sintheta;

    vec ntheta = normalize(vec(-nr.z*nr.x/sintheta,-nr.z*nr.y/sintheta,sintheta));

    return Br*nr+Btheta*ntheta;


}

void Earth_dipole::log_info() const
{
    cout<<"Earth dipole model at origin with Earh radius "<<Re <<constants::get_name_distance()<<" Reference feild strength"<<constants::get_name_B_field() <<endl;

    cout<<"Check  "<<get_field(vec(0,0,Re)).z<<" "<<get_field(vec(0,0,10*Re)).z<<endl;
}

solenoid::solenoid(vec O,vec D, unsigned int _N,double _I,double _R):field(true)
{
    dir = glm::normalize(D);
    I=_I;
    N=_N;
    origin = O;
    R=_R;
    inside_field =constants::get_mu0()*I*N*dir;
}

solenoid::solenoid(vec O,vec D, double B,double _R):field(true)
{
    dir = glm::normalize(D);
    I=-1;//Undefined
    N=-1;//Undefined
    origin = O;
    R=_R;
    inside_field =B*dir;
}



vec solenoid::get_field(vec pos,double t) const
{

    //If c is the distance to the origin, we want to find the distance to the center line
    //We want to find r = c sin(theta) = c sqrt(1-cos(theta)^2)
    double c = sqrt(glm::dot(pos-origin,pos-origin));

    //If pos == origin EXACTLY (within floating point precision) then we would get a zero division, but then we know we are at 0 distance anyway. It is not unthinkable that I would ask for exactly origin
    if (c==0)
        return inside_field;

    double costheta = glm::dot((pos-origin)/c,dir);
    double r = c*sqrt(1-costheta*costheta);

    return ((r<R) ? inside_field : vec(0));
}



void solenoid::log_info() const
{
    if (N != (unsigned int) -1)
        cout<<"Solenoid with direction vector ["<<dir.x<<' '<<dir.y<<' '<<dir.z<< "] going though ["<<origin.x<<' '<<origin.y<<' '<<origin.z<< "] " <<constants::get_name_distance()<<" with "<<N<<" turns per "<<constants::get_name_distance()<<" carrying "<<I<<constants::get_name_current()<<" with radius "<<R<<constants::get_name_distance()<<", where also μ0 = "<<constants::get_mu0()<<" "<<constants::get_name_force()<<" "<<constants::get_name_current()<<"^{-2}, resulting in field strengths "<< (constants::get_mu0()*I*N)<<" "<<constants::get_name_B_field()<<endl;
    else//If N is flagged as undefined, the simulation only knows the field strength
        cout<<"Solenoid with direction vector ["<<dir.x<<' '<<dir.y<<' '<<dir.z<< "] going though ["<<origin.x<<' '<<origin.y<<' '<<origin.z<< "] " <<constants::get_name_distance()<<" with field ["<<inside_field.x<<' '<<inside_field.y<<' '<<inside_field.z<< "] "<<constants::get_name_B_field()<<endl;


}

torus::torus(double R, double r, double B):field(true)
{
    Rminor =r;
    Rmajor =R;
    center_field=B;
}

void torus::log_info() const
{
    cout<<"Torus with major radius "<<Rmajor<<" " <<constants::get_name_distance()<< " minor radius "<<Rminor<<" " <<constants::get_name_distance()<<" with field "<<center_field<< " "<<constants::get_name_B_field()<<" at r = "<<Rmajor<<endl;


}

vec torus::get_field(vec pos,double t) const
{
    //get distance to origin axis
    double r = sqrt(pos.x*pos.x+pos.y*pos.y);

    if (pos.z> Rminor || pos.z<-Rminor)
        return vec(0);

    //What range of r is inside at this hide
    double this_r_minor= sqrt(Rminor*Rminor-pos.z*pos.z);

    if (r>Rmajor-this_r_minor && r<Rmajor+this_r_minor)
        return center_field*normalize(vec(-pos.y,pos.x,0));//*(Rmajor/r);
    else
        return vec(0);

}



Paul_trap::Paul_trap(double _R, double _omega, double _Qmax):field(false)
{
    R=_R;
    Qmax=_Qmax;
    omega=_omega;
    kQmax=Qmax/(4*PI*constants::get_eps0());
}

vec Paul_trap::get_field(vec pos,double t) const
{
    vec r1 = -pos+vec( R,0,0);
    vec r2 = -pos+vec(-R,0,0);
    vec r3 = -pos+vec(0, R,0);
    vec r4 = -pos+vec(0,-R,0);
    vec E1 = normalize(r1)*kQmax/(dot(r1,r1));
    vec E2 = normalize(r2)*kQmax/(dot(r2,r2));
    vec E3 = -normalize(r3)*kQmax/(dot(r3,r3));
    vec E4 = -normalize(r4)*kQmax/(dot(r4,r4));

    return E1+E2+E3+E4;
}

void Paul_trap::log_info() const
{
    cout<<"Paul trap, distance to poles "<<R<<" " <<constants::get_name_distance()<< " Oscillating field "<<Qmax<< " "<<constants::get_name_charge()<<" with frequency of oscillation "<<omega<<" "<<constants::get_name_inv_time()<<endl;
}


cyclotron_gab::cyclotron_gab(double _R, double _gab_size, double _Emax, double _omega, double _phase_offset):field(false)
{
    R=_R;
    gab_size=_gab_size;
    Emax=_Emax;
    omega=_omega;
    phase_offset=_phase_offset;
    half_gab_size = 0.5*gab_size;
}

vec cyclotron_gab::get_field(vec pos,double t) const
{
    //If we are somewhere between the two D shaped plates, apply a constant oscillating field (not exactly true near the edges, but good enough)
    if (pos.y>-half_gab_size && pos.y<half_gab_size)
        if (sqrt(pos.x*pos.x+pos.y*pos.y)<R)
            return vec(0,Emax*cos(omega*t),0);
    return vec(0);
}

void cyclotron_gab::log_info() const
{

    cout<<"Gab in a "<<R<<" " <<constants::get_name_distance()<< " cyclotron with size "<<gab_size<<" " <<constants::get_name_distance()<<" with max field "<<Emax<< " "<<constants::get_name_E_field()<<" with frequency of oscillation "<<omega<<" "<<constants::get_name_inv_time()<<endl;
}


vec composite_field::get_Bfield(vec pos,double t) const
{
    vec out=vec(0);
    for (field* F : fields)
        if (F!=nullptr)
            if (F->is_B())
                out+=F->get_field(pos,t);
    return out;
}

vec composite_field::get_Efield(vec pos,double t) const
{
    vec out=vec(0);
    for (field* F : fields)
        if (F!=nullptr)
            if (!(F->is_B()))
                out+=F->get_field(pos,t);
    return out;
}

void composite_field::log_info() const
{
    cout<<"Field made up off:"<<endl;
    for (field* F : fields)
        if (F!=nullptr)
            F->log_info();
}



