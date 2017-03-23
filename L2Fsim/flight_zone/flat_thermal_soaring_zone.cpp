#include <L2Fsim/flight_zone/flat_thermal_soaring_zone.hpp>

namespace L2Fsim {

using namespace std;

/*---------------------------------------------------------------------------------------
 ---------------------------            CONSTRUCTOR             -------------------------
 --------------------------------------------------------------------------------------*/

///Constructor. Create an empty zone defined by its dimension, wind, end time and height of thermals.
flat_thermal_soaring_zone::flat_thermal_soaring_zone(double Tend,int miX,int maX,int miY,int maY,
                                                     int miZ,int maZ,double wx,double wy,double Zi)
{
    tstart=0.0;
    tend=Tend;
    minX=miX;
    maxX=maX;
    minY=miY;
    maxY=maX;
    minZ=miZ;
    maxZ=maZ;
    windx = wx;
    windy = wy;
    zi=Zi;
}

///Constructor. Read a saved thermal scenario in order to play an identical simulation of the flight zone.
/**
    \param filename a file including a thermal scenario, saved thanks to the saveConfig function.
*/
flat_thermal_soaring_zone::flat_thermal_soaring_zone(string filename)
{
    ifstream file(filename);
    if(!file){cerr << "Unable to open the file ("<< filename<<") read by FZ !" << endl;}

    string  line;
    vector<double> data;
    vector<std::vector<double> > datavector;
    int readline=1;
    while(getline(file, line))
    {
        istringstream stream(line);
        if(line[0]=='#')
            continue;
        switch(readline)
        {
            case 1: stream>>minX>>maxX;break;
            case 2: stream>>minY>>maxY;break;
            case 3: stream>>minZ>>maxZ;break;
            case 4: stream>>tstart>>tend;break;
            case 5: stream>>windx>>windy;break;
            case 6: stream>>zi;break;
        }

        if(readline>6)
        {
            double th[]={0,0.0,0,0,0,0.0,0.0,0,0.0};
            stream>>th[0]>>th[1]>>th[2]>>th[3]>>th[4]>>th[5]>>th[6]>>th[7]>>th[8];
            data.push_back(th[0]);
            data.push_back(th[1]);
            data.push_back(th[2]);
            data.push_back(th[3]);
            data.push_back(th[4]);
            data.push_back(th[5]);
            data.push_back(th[6]);
            data.push_back(th[7]);
            data.push_back(th[8]);
            datavector.push_back(data);
            data.clear();
        }
        readline++;
    }
    file.close();
    for(auto data: datavector)
    {
        //std_thermal(mod,tB,XC0,YC0,ZC0,Zi,wstar,Lifetime,Ksi,read):
        std_thermal* newTH = new std_thermal(data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],1);
        newTH->setwind(windx,windy);
        thermals.push_back(newTH);
    }
    cout<<"Configurations read file : "<<filename<<" ; and thermals are created accordingly."<<endl;
}

/*---------------------------------------------------------------------------------------
 ---------------------------          PRIVATE METHODS           -------------------------
 --------------------------------------------------------------------------------------*/
///Define the maximum number of thermals in the zone.
int flat_thermal_soaring_zone::nbMaxThermals()
{
    double ziAverage = 1200;
    return floor(0.6*(maxX-minX)*(maxY-minY)/(dmin*ziAverage));
}

///Define the center of a new thermal, dependent on positions of other thermals.
/**
    \param center a list of the thermal centers in the simulation.
    \param t the simulated time.
*/
bool flat_thermal_soaring_zone::createThermalCenter(vector<double>& center, double t)
{
    bool newCenterIsValid = false;

    // find the farthest point from each thermals
    double farthestDist=0.;
    double farthestPointFromThermalsX=0.;
    double farthestPointFromThermalsY=0.;
    double dist;
    bool farthestCenterDontConflitWithAnotherCenter;
    for(int x=minX; x<=maxX;x+=20){
         for(int y=minY; y<=maxY;y+=20){
             dist=0.;
             farthestCenterDontConflitWithAnotherCenter=1;
             for(auto& th:thermals)
             {
                 if (th->isAlive(t)){
                 dist+=th->distToUpdraftCenter(x,y,0.);
                 bool farthestCenterDontConflitWithThisThermal=
                 sqrt( (x-th->getCenter()[0])*(x-th->getCenter()[0])
                      +(y-th->getCenter()[1])*(y-th->getCenter()[1]))
                 > 2*dmin;
                 farthestCenterDontConflitWithAnotherCenter=farthestCenterDontConflitWithAnotherCenter*farthestCenterDontConflitWithThisThermal;
                 }
             }
             if (farthestDist<dist && farthestCenterDontConflitWithAnotherCenter)
             {
                 farthestDist=dist;
                 farthestPointFromThermalsX=x;
                 farthestPointFromThermalsY=y;
             }
         }
    }

    // we're going to create the center according to a normal law around the farthestPoint we founded
    // nevertheless, it needs to be in the windfield and not in conflict with other thermals :
    double newCenterX=minX-1.;
    double newCenterY=minY-1.;
    while (!newCenterIsValid)
    {
        newCenterIsValid = true;

        newCenterX = farthestPointFromThermalsX+normalLaw()*dmin*2.5;
        newCenterY = farthestPointFromThermalsY+normalLaw()*dmin*2.5;

        // center is in the field
        while ( newCenterX<minX+dmin || newCenterX>maxX-dmin || newCenterY<minY+dmin || newCenterY>maxY-dmin)
        {
            newCenterX = farthestPointFromThermalsX+normalLaw()*dmin*2.5;
            newCenterY = farthestPointFromThermalsY+normalLaw()*dmin*2.5;
        }

        // Check it is valid for all the other centers
        for(auto& th: thermals)
        {
            if(th->isAlive(t))
            {
                bool newCenterIsValidForThisThermal =
                sqrt((newCenterX-th->getCenter()[0])*(newCenterX-th->getCenter()[0])
                     +(newCenterY-th->getCenter()[1])*(newCenterY-th->getCenter()[1]))
                > 2.0*dmin;
                newCenterIsValid=newCenterIsValid*newCenterIsValidForThisThermal;
            }
        }

    }
    center.push_back(newCenterX);
    center.push_back(newCenterY);
    center.push_back(0.);

    // return a boolean
    return(newCenterIsValid);
}

///Count the number of thermals alive at time t.
/**
    \param t the simulated time.
*/
int flat_thermal_soaring_zone::numberAliveAtTime(double t)
{
    int count=0;
    for(int i=0; i<thermals.size();i++)
    {
        thermal* th = thermals[i];
        if(th->isAlive(t))
        {
            count++;
        }
    }
    return count;
}


double flat_thermal_soaring_zone::environSink(double z, double t)
{
    double areaTh=0.0, massFlowTh=0.0, w_e;
    double avgUpdraft,radius,swd;
    double regX=maxX-minX;
    double regY=maxY-minY;
    for(auto th: thermals)
    {

        if (th->isAlive(t))
        {
            swd=((0.5 < (z / th->getzi())) && ((z / th->getzi()) < 0.9))?2.5 * (z/th->getzi() - 0.5):0;
            avgUpdraft=th->getw_star() * pow((z / th->getzi()),1.0/3.0) * (1.0 - 1.1*z/th->getzi())*(1-swd)*th->timeCoeff(t);
            radius = 0.102 * pow((z / th->getzi()),1/3) * (1 - 0.25*z/th->getzi()) * th->getzi();
            massFlowTh += avgUpdraft*PI*radius*radius;
            areaTh += PI*radius*radius;
        }

    }

    w_e = - massFlowTh / (regX*regY - areaTh);

    if(w_e>0.0)
    {
        w_e=0.0;
    }

    return(w_e);
}

/*---------------------------------------------------------------------------------------
 ---------------------------           PUBLIC METHODS           -------------------------
 --------------------------------------------------------------------------------------*/


/// Computes the wind vector w, at point (x,y,z), at time t.
/**
    \param x a horizontal coordinate of the system (x,y,z).
    \param y a horizontal coordinate of the system (x,y,z).
    \param z the vertical coordinate of the system (x,y,z).
    \param t the time.
    \param w the wind vector: windx, windy, windz.
*/
flat_thermal_soaring_zone& flat_thermal_soaring_zone::wind
    (double x, double y, double z, double t, vector<double> &w)
{
    vector<double>(3, 0.).swap(w);
    w[0]=windx;w[1]=windy;w[2]=0.;
    double  w_e=environSink(z,t);

    for(auto& therm: thermals)
    {
        if (therm->isAlive(t))
        {
            therm->wind(x,y,z,t,w);
        }
    }

    // compute environmental sink rate
    if (thermals.size()!=0) {
        if(thermals[0]->getModel()==1)
        {
            w[2] = w[2] + w_e;
        }
    }

    return *this;
}

///Simulate a scenario of thermals.
/**
    \param deltaT the interval of thermal actualization.
    \param model the chosen model for the simulation.
*/
void flat_thermal_soaring_zone::createScenario(double deltaT,int model)
{
    cout << "--> Create scenario" << endl;
    int nb_Alive_t,total=0;
    int nbThermalmax=nbMaxThermals();

    for (double t=tstart; t<tend ; t=t+deltaT)
    {
        cout << "t = " << t << endl;
        nb_Alive_t=numberAliveAtTime(t);

        while(nb_Alive_t<nbThermalmax)
        {
            vector<double> center;

            if(createThermalCenter(center,t))
            {
                //std_thermal(mod,tB,XC0,YC0,ZC0,Zi,wstar,Lifetime,Ksi,read):
                std_thermal* newTH = new std_thermal(model,t,center[0],center[1],center[2],zi,0.,0,0.,0);
                newTH->setwind(windx,windy);
                thermals.push_back(newTH);
                total++;
            }
            nb_Alive_t=numberAliveAtTime(t);
        }
    }
    cout << endl;
}

///Write the whole wind data for the visualization of a zslice in a file.
/**
    \param deltaT the interval of thermal actualization.
    \param deltax the definition of the mesh precision onto x-axis.
    \param deltay the definition of the mesh precision onto y-axis.
    \param zslice the height of the windfield you want to write.
    \param filename the name of the file you want to write data.
*/
void flat_thermal_soaring_zone::writeScenario
    (double deltaT, double deltax, double deltay, double zslice,string filename)
{
    cout << "--> Write scenario" << endl;
    ofstream file;
    file.open (filename);

    file << "t x y z updraft"<<endl;

    for (double t=tstart; t<tend ; t=t+deltaT) {
        cout << "t = " << t << endl;
        vector<double> w;
        for (int x=minX;x<maxX;x=x+deltax)
        {
            for (int y=minY;y<maxY;y=y+deltay)
            {
                //for (int z=300;z<1250;z=z+300)
                //{
                    //cout << t << "  " << x << "  " << y << endl;
                    this->wind(x,y,zslice,t,w);
                    file << t << " " ;
                    file << x << " " ;
                    file << y << " " ;
                    file << zslice << " " ;
                    file << w[2];
                    file << endl;
                //}
            }
        }
    }
    file.close();
    cout << endl;
}

///Save a thermal scenario in order to play it again in an other simulation.
/**
    \param filename the name of the .txt file you want to save the scenario.
*/
void flat_thermal_soaring_zone::saveConfig(string filename)
{
    cout << "--> Save config"<< endl;
    ofstream myfile;
    myfile.open (filename);
    myfile<<"# Definition of the domain :\n# MinX and MaxX"<<endl
    <<minX<<" "<<maxX<<endl
    <<"# MinY and MaxY"<<endl
    <<minY<<" "<<maxY<<endl
    <<"# MinZ and MaxZ"<<endl
    <<minZ<<" "<<maxZ<<endl
    <<"# tstart and tend"<<endl
    <<tstart<<" "<<tend<<endl
    <<"# windx and windy"<<endl
    <<windx<<" "<<windy<<endl
    <<"# zi "<<endl
    <<zi<<endl;

    myfile<<"# The list of Thermals created :"<<endl
    <<"# model   tBirth   CentreX   CentreY   CentreZ   zi   w_star   lifetime   ksi"<<endl;

    for(auto th: thermals)
    {
        myfile<<th->getModel()<< " " <<th->gettBirth() << " " << th->getCenter()[0]<<" "<<th->getCenter()[1]<<" "<<th->getCenter()[2]<<" "<<th->getzi()<< " " << th->getw_star() << " " << th->getlifeTime()  << " " << th->getksi()<< endl;
    }
    myfile.close();
}

///Save a thermal scenario in order to play it again in an other simulation.
/**
    \param filename the name of the .csv file you want to save the scenario.
*/
void flat_thermal_soaring_zone::saveConfigToCSV(string filename)
{
    ofstream myfile;
    myfile.open (filename);
    myfile<<"model tBirth CentreX CentreY CentreZ zi w_star lifetime ksi"<<endl;

    for(auto th: thermals)
    {
        myfile<<th->getModel()<< " " <<th->gettBirth() << " " << th->getCenter()[0]<<" "<<th->getCenter()[1]<<" "<<th->getCenter()[2]<<" "<<th->getzi()<< " " << th->getw_star() << " " << th->getlifeTime()  << " " << th->getksi()<< endl;
    }
    myfile.close();
}
}
