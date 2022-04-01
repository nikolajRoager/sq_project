#include"particle.hpp"
//GLM, openGL Mathematics library, originally created to do the matrix and vector calculations needed for OpenGL powered 3D programs, it is also a decent vector and matrix mathematics library.
#include<glm/glm.hpp>
#include<cstdint>

using namespace boost::numeric::odeint;

using uint = uint32_t;
using vec = glm::dvec3;//Default to glm::douple precision 3D vector, by default glm uses single precision for comptatibility with the GLSL shading language (Single precision is good enough for computer-graphics application where speed is more important than accuracy)

#include"field.hpp"
#include"constants.hpp"

#include<iostream>


double particle::print_interval=0.1;

void particle::set_print_interval(double interval)
{
    print_interval=interval;
}

void particle::dump_raw(std::ofstream& txt_output) const
{
    cout<<"Saving "<<positions.size()<<" data points"<<endl;
    for (size_t i =0; i < positions.size(); ++i)
    {
        txt_output<<timestamps[i]<<'\t';
        txt_output<<positions[i].x<<'\t';
        txt_output<<positions[i].y<<'\t';
        txt_output<<positions[i].z<<'\t';
        txt_output<<velocities[i].x<<'\t';
        txt_output<<velocities[i].y<<'\t';
        txt_output<<velocities[i].z;

        if (i+1<positions.size())
            txt_output<<'\t'<< (timestamps[i+1]-timestamps[i]);
        txt_output<<endl;

    }
}

void particle::write_final(std::ofstream& txt_output) const
{
        size_t i = positions.size()-1;
        txt_output<<positions[i].x<<'\t';
        txt_output<<positions[i].y<<'\t';
        txt_output<<positions[i].z<<'\t';
        txt_output<<velocities[i].x<<'\t';
        txt_output<<velocities[i].y<<'\t';
        txt_output<<velocities[i].z<<endl;
}


particle::particle(vec _pos0, vec _v0, double m, double q, string _name)
{
    v0=_v0;
    pos0=_pos0;
    positions = vector<vec>();
    velocities=vector<vec>();
    timestamps=vector<double>();
    charge=q;
    mass = m;
    inv_mass = 1/mass;
    name = _name;

}

//Time Crunch, assume the particle is rotating around the y=z=0 axis
void particle::analytical_solenoid(const composite_field& Fields, double T, double dt)
{
    double Bx    = Fields.get_Bfield(vec(0,0,0),0).x;
    double vperp = sqrt(v0.y*v0.y+v0.z*v0.z);

    double omega =  -charge*Bx*inv_mass;

    double R     =  -vperp*mass/(charge*Bx);



    double theta0=0;


    size_t time_res = T/print_interval;


    positions  = vector<vec>   (time_res);
    velocities = vector<vec>   (time_res);
    timestamps = vector<double>(time_res);

    for (size_t i = 0; i < time_res; ++i)
    {
        double time = i*print_interval;
        positions[i] =vec(pos0.x+time*v0.x, R*cos(omega*time),R*sin(omega*time));
        velocities[i]=vec(v0.x, vperp*sin(omega*time),vperp*cos(omega*time));
        timestamps[i]=time;
    }

}

void particle::calculate_RK4(const composite_field& Fields, double T, double dt)
{
    cout<<"Running Runge Kutta 4 "<<endl;
    state_type Data0 = {pos0.x,pos0.y,pos0.z,v0.x,v0.y,v0.z};


    double Charge = charge;
    double Inv_mass= inv_mass;

    //Declare the differential equation we want to solve
    auto ODE = [Inv_mass,Charge,&Fields]( const state_type Data , state_type &dDatadt , const double t )
    {
        //Extract position and velocity from data
        vec pos = vec(Data[0],Data[1],Data[2]);
        vec velocity = vec(Data[3],Data[4],Data[5]);

        //Get current force
        vec dpdt = Charge*(Fields.get_Efield(pos,t)+glm::cross(velocity,Fields.get_Bfield(pos,t)));
        //force -> acceleration
        vec dV = dpdt*Inv_mass;

        //Save acceleration as derivative of velocity
        dDatadt[3] = dV.x;
        dDatadt[4] = dV.y;
        dDatadt[5] = dV.z;
        //velocity as derivative of position
        dDatadt[0] = Data[3];
        dDatadt[1] = Data[4];
        dDatadt[2] = Data[5];
    };



    vector<vec> position_out;
    position_out.push_back(pos0);

    vector<vec> velocity_out;
    velocity_out.push_back(v0);

    vector<double> time_out;
    time_out.push_back(0);


    double p_print = 0;
    double Print_interval = print_interval;


    auto save_step = [&p_print,Print_interval,&position_out,&velocity_out,&time_out]( const state_type& Data , const double t )
    {   //Data[0]=position, Data[1]=velocity
        if (t>p_print+Print_interval)
        {
            position_out.push_back(vec(Data[0],Data[1],Data[2]));//I do not know how many points there will be, so reserving is not possible, hopefully there will not be too many resizes, //The print inteval ensures that there will not be crazy much data to save
            velocity_out.push_back(vec(Data[3],Data[4],Data[5]));
            time_out.push_back(t);

            p_print=t;
        }
    };




    state_type Data = Data0;
    state_type temp=Data0;
    state_type K1,K2,K3,K4;
    size_t time_res = T/dt;
    for (size_t i = 1; i < time_res; ++i)
    {
        double t=i*dt;

        //substep 1
        ODE(Data,K1,t);
        for (uint i = 0; i<Data.size(); ++i)
            temp[i]=Data[i]+dt*K1[i]/2;

        //substep 2
        ODE(temp,K2,t+dt/2);
        for (uint i = 0; i<Data.size(); ++i)
            temp[i]=Data[i]+dt*K2[i]/2;

        //substep 3
        ODE(temp,K3,t+dt/2);
        for (uint i = 0; i<Data.size(); ++i)
            temp[i]=Data[i]+dt*K3[i];

        //substep 4
        ODE(temp,K4,t+dt);

        //Read data

        for (uint i = 0; i<Data.size(); ++i)
            Data[i]+=dt*(K1[i]+2.0*K2[i]+2.0*K3[i]+K4[i])/6.0;

        save_step( Data , i*dt );
    }

    cout<<"Solved in "<<time_res<<" steps (t=0 "<<constants::get_name_time()<<" to t="<<T<<' '<<constants::get_name_time()<<')'<<endl;
    //std::move only moves the pointer to the data, so it is much faster than the default copy assignment, which would be used if I just used =, ... in practice there is a good chance the compiler is smart enough to use move here instead, but being explicit never hurts
    positions = std::move(position_out);
    velocities = std::move(velocity_out);
    timestamps = std::move(time_out);

}
void particle::calculate_euler(const composite_field& Fields, double T, double dt)
{
    cout<<dt<<" ... ";
    state_type Data0 = {pos0.x,pos0.y,pos0.z,v0.x,v0.y,v0.z};


    double Charge = charge;
    double Inv_mass= inv_mass;

    //Declare the differential equation we want to solve
    auto ODE = [Inv_mass,Charge,&Fields]( const state_type Data , state_type &dDatadt , const double t )
    {
        //Extract position and velocity from data
        vec pos = vec(Data[0],Data[1],Data[2]);
        vec velocity = vec(Data[3],Data[4],Data[5]);

        //Get current force
        vec dpdt = Charge*(Fields.get_Efield(pos,t)+glm::cross(velocity,Fields.get_Bfield(pos,t)));
        //force -> acceleration
        vec dV = dpdt*Inv_mass;

        //Save acceleration as derivative of velocity
        dDatadt[3] = dV.x;
        dDatadt[4] = dV.y;
        dDatadt[5] = dV.z;
        //velocity as derivative of position
        dDatadt[0] = Data[3];
        dDatadt[1] = Data[4];
        dDatadt[2] = Data[5];
    };



    vector<vec> position_out;
    position_out.push_back(pos0);

    vector<vec> velocity_out;
    velocity_out.push_back(v0);

    vector<double> time_out;
    time_out.push_back(0);


    double p_print = 0;
    double Print_interval = print_interval;


    auto save_step = [&p_print,Print_interval,&position_out,&velocity_out,&time_out]( const state_type& Data , const double t )
    {   //Data[0]=position, Data[1]=velocity
        if (t>p_print+Print_interval)
        {
            position_out.push_back(vec(Data[0],Data[1],Data[2]));//I do not know how many points there will be, so reserving is not possible, hopefully there will not be too many resizes, //The print inteval ensures that there will not be crazy much data to save
            velocity_out.push_back(vec(Data[3],Data[4],Data[5]));
            time_out.push_back(t);

            p_print=t;
        }
    };




    state_type Data = Data0;
    state_type dDatadt;
    size_t time_res = T/dt;
    for (size_t i = 1; i < time_res; ++i)
    {
        double t=i*dt;
        ODE(Data,dDatadt,t);
        //Euler time evolution
        //Data +=dt*dDatadt; 1 variable
        for (uint i = 0; i<Data.size(); ++i)
            Data[i]+=dt*dDatadt[i];
        save_step( Data , i*dt );
    }

    cout<<"Solved in "<<time_res<<" steps (t=0 "<<constants::get_name_time()<<" to t="<<T<<' '<<constants::get_name_time()<<')'<<endl;
    //std::move only moves the pointer to the data, so it is much faster than the default copy assignment, which would be used if I just used =, ... in practice there is a good chance the compiler is smart enough to use move here instead, but being explicit never hurts
    positions = std::move(position_out);
    velocities = std::move(velocity_out);
    timestamps = std::move(time_out);

}

void particle::calculate_3rdparty_RK4(const composite_field& Fields, double T, double dt)
{




    state_type Data0 = {pos0.x,pos0.y,pos0.z,v0.x,v0.y,v0.z};

    //p_print=0;
    double Charge = charge;
    double Inv_mass= inv_mass;

    //Declare the differential equation we want to solve
    auto ODE = [Inv_mass,Charge,&Fields]( const state_type Data , state_type &dDatadt , const double t )
    {
        //Extract position and velocity from data
        vec pos = vec(Data[0],Data[1],Data[2]);
        vec velocity = vec(Data[3],Data[4],Data[5]);

        //Get current force
        vec dpdt = Charge*(Fields.get_Efield(pos,t)+glm::cross(velocity,Fields.get_Bfield(pos,t)));
        //force -> acceleration
        vec dV = dpdt*Inv_mass;

        //Save acceleration as derivative of velocity
        dDatadt[3] = dV.x;
        dDatadt[4] = dV.y;
        dDatadt[5] = dV.z;
        //velocity as derivative of position
        dDatadt[0] = Data[3];
        dDatadt[1] = Data[4];
        dDatadt[2] = Data[5];
    };


    vector<vec> position_out;
    position_out.push_back(pos0);

    vector<vec> velocity_out;
    velocity_out.push_back(v0);

    vector<double> time_out;
    time_out.push_back(0);
    double p_print = 0;
    double Print_interval = print_interval;


    auto save_step = [&p_print,Print_interval,&position_out,&velocity_out,&time_out]( const state_type& Data , const double t )
    {   //Data[0]=position, Data[1]=velocity
        if (t>p_print+Print_interval)
        {
            position_out.push_back(vec(Data[0],Data[1],Data[2]));//I do not know how many points there will be, so reserving is not possible, hopefully there will not be too many resizes, //The print inteval ensures that there will not be crazy much data to save
            velocity_out.push_back(vec(Data[3],Data[4],Data[5]));
            time_out.push_back(t);

            p_print=t;
        }
    };

    size_t steps = integrate_const(
        runge_kutta4< state_type >(),
        ODE,   //Lorentz-force
        Data0 ,//{pos0,v0}
        0.0 ,  //t0=0
        T ,    //max time
        dt,//length of each step
        save_step //User defined save data function
);


    cout<<"Solved in "<<steps<<" steps (t=0 "<<constants::get_name_time()<<" to t="<<T<<' '<<constants::get_name_time()<<')'<<endl;

    //std::move only moves the pointer to the data, so it is much faster than the default copy assignment, which would be used if I just used =, ... in practice there is a good chance the compiler is smart enough to use move here instead, but being explicit never hurts
    positions = std::move(position_out);
    velocities = std::move(velocity_out);
    timestamps = std::move(time_out);

}

void particle::calculate_RKDP45(const composite_field& Fields, double T, double dt)
{
    double h_max = print_interval;
    state_type Data0 = {pos0.x,pos0.y,pos0.z,v0.x,v0.y,v0.z};

    //p_print=0;
    double Charge = charge;
    double Inv_mass= inv_mass;

    //Declare the differential equation we want to solve
    auto ODE = [Inv_mass,Charge,&Fields]( const state_type Data , state_type &dDatadt , const double t )
    {
        //Extract position and velocity from data
        vec pos = vec(Data[0],Data[1],Data[2]);
        vec velocity = vec(Data[3],Data[4],Data[5]);

        //Get current force
        vec dpdt = Charge*(Fields.get_Efield(pos,t)+glm::cross(velocity,Fields.get_Bfield(pos,t)));
        //force -> acceleration
        vec dV = dpdt*Inv_mass;

        //Save acceleration as derivative of velocity
        dDatadt[3] = dV.x;
        dDatadt[4] = dV.y;
        dDatadt[5] = dV.z;
        //velocity as derivative of position
        dDatadt[0] = Data[3];
        dDatadt[1] = Data[4];
        dDatadt[2] = Data[5];
    };


    vector<vec> position_out;
    position_out.push_back(pos0);

    vector<vec> velocity_out;
    velocity_out.push_back(v0);

    vector<double> time_out;
    time_out.push_back(0);
    double p_print = 0;
    double Print_interval = print_interval;


    auto save_step = [&p_print,Print_interval,&position_out,&velocity_out,&time_out]( const state_type& Data , const double t )
    {   //Data[0]=position, Data[1]=velocity
        if (t>p_print+Print_interval)
        {
            position_out.push_back(vec(Data[0],Data[1],Data[2]));//I do not know how many points there will be, so reserving is not possible, hopefully there will not be too many resizes, //The print inteval ensures that there will not be crazy much data to save
            velocity_out.push_back(vec(Data[3],Data[4],Data[5]));
            time_out.push_back(t);

            p_print=t;
        }
    };

    cout<<"Using Custom RK DP 45"<<flush;

    size_t steps = 0;


    double t=0,h=dt;

    //Relative and absolute error on all the data
    double relErr=10e-7;
    double absErr=10e-7;


    state_type Data = Data0;
    state_type temp=Data0;
    state_type K1,K2,K3,K4,K5,K6,K7;



    //Very first step requires one more step than the rest
    ODE(Data,K1,t);
    size_t rejections = 0;
    while (t<T)
    {
        bool reject = false;

        do
        {

            //Do not overshoot
            if(t+h>T)
                h =T-t;

            //substep 1, was  calculated that last substep of last step
            //ODE(Data,K1,t);


            //substep 2
            for (uint i = 0; i<Data.size(); ++i)
                temp[i]=Data[i]+h*K1[i]/5.0;
            ODE(temp,K2,t+h/5.0);


            //substep 3
            for (uint i = 0; i<Data.size(); ++i)
                temp[i]=Data[i]+h*K1[i]*3/40.0+h*K2[i]*9/40.0;
            ODE(temp,K3,t+h*3/10.0);

            //substep 4
            for (uint i = 0; i<Data.size(); ++i)
                temp[i]=Data[i]+h*K1[i]*44/45.0-h*K2[i]*56/15.0+h*K3[i]*32/9.0;
            ODE(temp,K4,t+h*4/5.0);


            //substep 5
            for (uint i = 0; i<Data.size(); ++i)
                temp[i]=Data[i]+h*K1[i]*19372/6561.0-h*K2[i]*25360/2187.0+h*K3[i]*64448/6561.0-h*K4[i]*212/729.0;
            ODE(temp,K5,t+h*8/9.0);//   19372/6561  -        25360/2187           64448/6561  -        212/729

            //substep 6
            for (uint i = 0; i<Data.size(); ++i)
                temp[i]=Data[i]+h*K1[i]*9017/3168.0-h*K2[i]*355/33.0+h*K3[i]*46732/5247.0+h*K4[i]*49/176.0-h*K5[i]*5103/18656.0;
            ODE(temp,K6,t+h);//         9017/3168  -        355/33           46732/5247           49/176  -      5103/18656

            //substep 7, and ALSO the 5'th order estimate
            state_type Data_5 = Data;
            for (uint i = 0; i<Data.size(); ++i)
                Data_5[i]=Data[i]+h*K1[i]*35/384.0+h*K3[i]*500/1113.0+h*K4[i]*125/192.0-h*K5[i]*2187/6784.0+h*K6[i]*11/84.0;
            ODE(Data_5,K7,t+h);//         35/384  0        500/1113           125/192  −        2187/6784           11/84

            state_type Data_4 = Data;//
            //for (uint i = 0; i<Data.size(); ++i)
            //    Data_5[i]+=h*((35/384.0)*K1[i]+(500/1113.0)*K3[i]+(125/192.0)*K4[i]-(2187/6784.0)*K5[i]+(11/84.0)*K6[i]);
            for (uint i = 0; i<Data.size(); ++i)
                Data_4[i]+=h*((5179/57600.0)*K1[i]+(7571/16695.0)*K3[i]+(393/640.0)*K4[i]-(92097/339200.0)*K5[i]+(187/2100.0)*K6[i]+(1/40.0)*K7[i]);

            //Estimate the error
            state_type Error;
            for (uint i = 0; i<Data.size(); ++i)
                Error[i] = abs(Data_5[i]-Data_4[i]);


            double h_new=h;
            reject = false;
            bool first_run = true;
            for (uint i = 0; i<Data.size() ; ++i)
            {
                double delta = sqrt(h)*(relErr*abs(Data[i])+absErr);
                if (Error[i]>delta)
                {
                    reject = true;
                }

                if (Error[i]!=0)
                {
                    double my_h_new = 0.9*h*pow(delta/Error[i],0.25);//What step size would fix THIS

                    //If this is the worst parameter so far (or the first this loop), resize to fit that
                    if (my_h_new<h_new || first_run )
                    {
                        h_new =my_h_new;

                        if (h_new>h_max && h_max !=0)//Don't overstep, if h_max is 0, no limit is used
                        {
                            h_new =h_max;
                        }

                        first_run = false;
                    }
                }
            }
            if (!reject)
            {
                t+=h;
                Data=Data_5;
            }
            else
                ++rejections;

            h=h_new;

        }
        while(reject);

        //If we KEPT the previous step, then K7 is the same as K1 for the next step
        K1=K7;

        ++steps;
        save_step( Data , t);

    }

    cout<<" Solved in "<<steps<<" steps and "<<rejections<<" rejections (t=0 "<<constants::get_name_time()<<" to t="<<T<<' '<<constants::get_name_time()<<')'<<endl;


    //std::move only moves the pointer to the data, so it is much faster than the default copy assignment, which would be used if I just used =, ... in practice there is a good chance the compiler is smart enough to use move here instead, but being explicit never hurts
    positions = std::move(position_out);
    velocities = std::move(velocity_out);
    timestamps = std::move(time_out);

}


void particle::calculate(const composite_field& Fields, double T, double dt)
{




    state_type Data0 = {pos0.x,pos0.y,pos0.z,v0.x,v0.y,v0.z};

    //p_print=0;
    double Charge = charge;
    double Inv_mass= inv_mass;

    //Declare the differential equation we want to solve
    auto ODE = [Inv_mass,Charge,&Fields]( const state_type Data , state_type &dDatadt , const double t )
    {
        //Extract position and velocity from data
        vec pos = vec(Data[0],Data[1],Data[2]);
        vec velocity = vec(Data[3],Data[4],Data[5]);

        //Get current force
        vec dpdt = Charge*(Fields.get_Efield(pos,t)+glm::cross(velocity,Fields.get_Bfield(pos,t)));
        //force -> acceleration
        vec dV = dpdt*Inv_mass;

        //Save acceleration as derivative of velocity
        dDatadt[3] = dV.x;
        dDatadt[4] = dV.y;
        dDatadt[5] = dV.z;
        //velocity as derivative of position
        dDatadt[0] = Data[3];
        dDatadt[1] = Data[4];
        dDatadt[2] = Data[5];
    };


    vector<vec> position_out;
    position_out.push_back(pos0);

    vector<vec> velocity_out;
    velocity_out.push_back(v0);

    vector<double> time_out;
    time_out.push_back(0);
    double p_print = 0;
    double Print_interval = print_interval;


    auto save_step = [&p_print,Print_interval,&position_out,&velocity_out,&time_out]( const state_type& Data , const double t )
    {   //Data[0]=position, Data[1]=velocity
        if (t>p_print+Print_interval)
        {
            position_out.push_back(vec(Data[0],Data[1],Data[2]));//I do not know how many points there will be, so reserving is not possible, hopefully there will not be too many resizes, //The print inteval ensures that there will not be crazy much data to save
            velocity_out.push_back(vec(Data[3],Data[4],Data[5]));
            time_out.push_back(t);

            p_print=t;
        }
    };

    cout<<"Using odeint "<<flush;

    //All these things live in the namespace boost::numeric::odeint
    size_t steps = integrate_adaptive(
        //runge_kutta_dopri5< state_type >(),
        make_controlled( 1e-7 , 1e-7 , runge_kutta_dopri5< state_type >() ) ,//Create stepper
        ODE,   //Lorentz-force
        Data0 ,//{pos0,v0}
        0.0 ,  //t0=0
        T ,    //max time
        dt ,   //Initial dt, odeint will change if it needs to
        save_step //User defined save data function
    );

/*
size_t steps = integrate_const(
    runge_kutta4< state_type >(),
    ODE,   //Lorentz-force
    Data0 ,//{pos0,v0}
    0.0 ,  //t0=0
    T ,    //max time
    dt,//length of each step
    save_step //User defined save data function
);
*/

    cout<<"Solved in "<<steps<<" steps (t=0 "<<constants::get_name_time()<<" to t="<<T<<' '<<constants::get_name_time()<<')'<<endl;

    //std::move only moves the pointer to the data, so it is much faster than the default copy assignment, which would be used if I just used =, ... in practice there is a good chance the compiler is smart enough to use move here instead, but being explicit never hurts
    positions = std::move(position_out);
    velocities = std::move(velocity_out);
    timestamps = std::move(time_out);

}


//For saving, return the ID of the particles position
const char* particle::positionID() const {return (const char*) &(positions[0]);}
