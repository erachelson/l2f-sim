#ifndef L2FSIM_THERMAL_HPP_
#define L2FSIM_THERMAL_HPP_

#include <vector>

namespace L2Fsim {

class thermal {
    
    /*--------------------------------------
     ------------- Attributes --------------
     -------------------------------------*/
    
protected:
    
    double windx,windy;     // wind implying the drift of the thermal [m/s]
    
    /*--------------------------------------
     --------------- Methods ---------------
     -------------------------------------*/
    
public:
        
    // get
    virtual double getw_star() =0;
    virtual double gettBirth() =0;
    virtual double getlifeTime() =0;
    virtual double getzi() =0;
    virtual int getModel() =0;
    virtual double getksi() =0;
    virtual std::vector<double> getCenter() =0;
    
    // method
    virtual bool isAlive(double t) =0;
    virtual double timeCoeff(double t) =0;
    virtual double distToUpdraftCenter(double x, double y, double z) =0;
    
    // set
    void setwind(double wx,double wy)
    {
        windx=wx;
        windy=wy;
    }
    
    // wind
	virtual thermal& wind(double x,double y,double z,double t,std::vector<double> &w) =0;
};

}

#endif
