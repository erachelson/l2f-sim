#ifndef L2FSIM_BEELERGLIDER_HPP_
#define L2FSIM_BEELERGLIDER_HPP_

#include <L2Fsim/aircraft/aircraft.hpp>
#include <L2Fsim/flight_zone/flight_zone.hpp>
#include <L2Fsim/utils/quaternion.hpp>
#include <L2Fsim/utils/utils.hpp>
#include <vector>
#include <cmath>

namespace L2Fsim {
    
    /// Beeler's glider model.
    /**
     * Equations derived from:
     * Beeler, Moerder and Cox. A Flight Dynamics Model for a Small Glider
     * in Ambient Winds. NASA/TM-2003-212665. 2003.
     * Model validity:
     * - wingspan in [1.5m; 3.5m]
     * - aspect ratio in [6; 20]
     * - weight in [0.22kg; 5.45kg]
     * Notations:
     * (x,y,z) glider position in an Earth-based coordinate system
     * V ground speed of the glider
     * gamma elevation angle
     * khi azimuth angle
     * sigma bank angle
     * alpha angle of attack
     * beta sideslip angle
     */
    
    class BeelerGlider : public aircraft {
    public:
        
        /* attributes */
        
        /** Aircraft state: (x,y,z,V,gamma,khi) */
        std::vector<double> state;
        /** Aircraft mass */
        double mass;
        /** Wingspan */
        double wingspan;
        /** Aspect ratio */
        double aspect_ratio;
        
        /** Aspect ratio of vertical tail */
        double ARv = 0.5*aspect_ratio;
        /** Fuselage area */
        double Sf = 86 * (2.54/100) * (2.54/100);
        /** Fuselage moment arm length */
        double lt = 0.28*(2.54/100);
        /** Horizontal tail volume ratio */
        double Vh = 0.4;
        /** Vertical tail volume ratio */
        double Vv = 0.02;
        /** Mean aerodynamic cord */
        double c_ = 1.03*wingspan/aspect_ratio;
        /** Wing surface area */
        double S = wingspan*wingspan/aspect_ratio;
        /** Horizontal tail surface */
        double St = Vh *c_*S/lt;
        /** Vertical tail surface */
        double Sv = Vv*wingspan*S/lt;
        
        /** Oswald efficiency number */
        double e = 0.95;
        /** Reynolds number */
        double Re = 150000;
        /** Lift curve slope */
        double a0 = 0.1*(180/3.14);
        /** Zero point */
        double alpha0 = -2.5*(3.14/180);
        /** Minimum wing profile drag */
        double Cd0 = 0.01;
        double Cdl = 0.05;
        /** Minimum lift */
        double Clmin = 0.4;
        double Cl_alpha = a0/(1+a0/(3.14*e*aspect_ratio));
        double Cc_beta = a0 / (1+a0*(3.14*ARv)) * (Sv / S);
        
        /* methods */
        
        /// Constructor
        /**
         \param m the aircraft's mass.
         \param ws the aircraft's wingspan.
         \param ar the aircraft's aspect ratio.
         */
        BeelerGlider(double m=1.35, double ws=2., double ar=16.)
        : mass(m),
        wingspan(ws),
        aspect_ratio(ar) {
            
        }
        
        /// Computes and updates the aircraft's state vector.
        /**
         \param st the aircraft's state vector.
         \param cd the aircraft's command vector: alpha, beta, sigma.
         \param fz the flight zone.
         \param t the time.
         \param statedot the derivatives vector: x_dot, y_dot, z_dot, V_dot, khi_dot, gamma_dot.
         */
        const BeelerGlider& state_dynamics(const std::vector<double> &st,
                                           const std::vector<double> &cd,
                                           flight_zone &fz,
                                           double t,
                                           std::vector<double> &statedot) const {
            double lift =0;
            double drag =0;
            double sideforce=0;
            
            //printf("Ready to caculate aero force\n" );
            this->calcAeroForces(st, cd, fz, t, lift, drag, sideforce);
            std::vector<double>(6, 0.).swap(statedot);
            
            //printf("Calculation finishing\n" );
            const double &x = st.at(0);
            const double &y = st.at(1);
            const double &z = st.at(2);
            const double &V = st.at(3);
            const double &gamma = st.at(4);
            const double &khi = st.at(5);
            const double &alpha = cd.at(0);
            const double &beta = cd.at(1);
            const double &sigma = cd.at(2);
            
            statedot.at(0) = V * cos(gamma) * cos(khi);
            statedot.at(1) = V * cos(gamma) * sin(khi);
            statedot.at(2) = V * sin(gamma);
            statedot.at(3) = - drag / mass - 9.81 * sin(gamma);
            statedot.at(4) = (lift * cos(sigma) + sideforce * sin(sigma)) / (mass * V) - 9.81 * cos(gamma) / V;
            statedot.at(5) = (lift * sin(sigma) - sideforce * cos(sigma)) / (mass * V * cos(gamma));
            
            return *this;
        }
        
        /// Gets the observation vector.
        /**
         \param observation the observation vector.
         */
        BeelerGlider& observation(std::vector<double> &observation)  {
            observation = {state.at(2),state.at(3),state.at(4),state.at(5)};
            return *this;
        }
        
        /// Gets the aircraft's state vector.
        /**
         \param state the vector to copy the aircraft's state in.
         */
        BeelerGlider& get_state(std::vector<double> &st)  {
            st = state;
            return *this;
        }
        
        /// Sets the aircraft's state vector.
        /**
         \param state the vector to be copied to the aircraft's state.
         */
        BeelerGlider& set_state(const std::vector<double> &st) {
            state = st;
            return *this;
        }
        
        /// Checks if the state vector contains values that are out of the model's range of validity.
        /**
         \param state the aircraft's state vector.
         \param command the aircraft's command vector: alpha, beta, sigma.
         */
        BeelerGlider& is_in_model(std::vector<double> &st,std::vector<double> &cd) {
            
            double limit_gamma_angle = 0.3;
            double limit_alphagamma_angle = 0.3;
            
            if(st.at(2) < 0)
            {
                printf("The aircraft is below 0 (altitude)\n");
                exit(-1);
            }
            if(st.at(4) > limit_gamma_angle){
                printf("The aircraft have an inclinaison (gamma) > 17 degree\n");
                exit(-1);
            }
            if(st.at(4) < -limit_gamma_angle){
                printf("The aircraft have an inclinaison (gamma) < -17 degree\n");
                exit(-1);
            }
            
            if(st.at(4) + cd.at(1) < -limit_alphagamma_angle){
                printf("The aircraft have an bank (gamma + alpha) < -17 degree\n");
                exit(-1);
            }
            if(st.at(4) + cd.at(1) > limit_alphagamma_angle){
                printf("The aircraft have an bank (gamma + alpha) > 17 degree\n");
                exit(-1);
            }
            
            return *this;
        }
        
    protected:
        
        /// Computes lift, drag and sideforce.
        /**
         \param st the aircraft's state vector.
         \param cd the aircraft's command vector: alpha, beta, sigma.
         \param fz the flight zone.
         \param t the time.
         \param lift the lift.
         \param drag the drag.
         \param statedot the derivatives vector: x_dot, y_dot, z_dot, V_dot, khi_dot, gamma_dot.
         */
        void calcAeroForces(const std::vector<double> &st,
                            const std::vector<double> &cd,
                            flight_zone &fz,
                            double t,
                            double &lift,
                            double &drag,
                            double &sideforce) const {
            
            /** Retrieve aircraft's state */
            const double &x = st.at(0);
            const double &y = st.at(1);
            const double &z = st.at(2);
            const double &V = st.at(3);
            const double &gamma = st.at(4);
            const double &khi = st.at(5);
            const double &alpha = cd.at(0);
            const double &beta = cd.at(1);
            const double &sigma = cd.at(2);
            
            /** Relative wind */
            std::vector<double> w(3);
            
            /** Wind relative velocity */
            std::vector<double> V_w(3);
            std::vector<double> X_w(3);
            
            /** Wind relative angles */
            double alpha_w, beta_w, gamma_w, khi_w, sigma_w;
            
            /** Aerodynamic force coefficients with wind */
            double Cc_w, Cl_w, Cd_w;
            
            /** Dynamic pressure */
            double q;
            
            /** Rotation matrices and quaternions */
            quaternion rviq; // quaternion version of the Euler rotation sequence R_VI
            
            std::vector<double> rbv(9); // rotation from velocity frame to body frame R_BV
            std::vector<double> rbv1(9); // rotation of alpha
            std::vector<double> rbv2(9); // rotation of -beta
            quaternion rbvq; //quaternion version R_BV
            quaternion rbv1q;
            quaternion rbv2q;
            
            std::vector<double> m(9); // M matrix, used to retrieve wind relative angles
            std::vector<double> m11(9);
            std::vector<double> m12(9);
            quaternion mq; // quaternion version of M
            quaternion m11q;
            quaternion m12q;
            
            /** Calc relative wind */
            fz.wind(x, y, z, t, w);
            
            V_w.at(0) = V * cos(gamma) * cos(khi) - w.at(0);
            V_w.at(1) = V * cos(gamma) * sin(khi) - w.at(1);
            V_w.at(2) = V * sin(gamma) - w.at(2);
            
            X_w.at(0) = V_w.at(0) / sqrt(V_w.at(0)*V_w.at(0) + V_w.at(1)*V_w.at(1) + V_w.at(2)*V_w.at(2));
            X_w.at(1) = V_w.at(1) / sqrt(V_w.at(0)*V_w.at(0) + V_w.at(1)*V_w.at(1) + V_w.at(2)*V_w.at(2));
            X_w.at(2) = V_w.at(2) / sqrt(V_w.at(0)*V_w.at(0) + V_w.at(1)*V_w.at(1) + V_w.at(2)*V_w.at(2));
            
            /** Calc of gamma_w and khi_w */
            gamma_w = asin(X_w.at(2));
            
            if(X_w.at(0) / cos(gamma_w) > 1) {
                khi_w = 0;
            }else if(X_w.at(0) / cos(gamma_w) < -1){
                khi_w = M_1_PI;
            }else{
                khi_w = sgn(X_w.at(1) / cos(gamma_w)) * acos(X_w.at(0) / cos(gamma_w));
            }
            
            /** Calc of alpha_w, beta_w and sigma_w :
             Use of rotation matrices and quaternions */
            
            rviq.fromEuler(khi, gamma, sigma);
            
            rbv1.at(0) = cos(alpha);
            rbv1.at(2) = sin(alpha);
            rbv1.at(4) = 1;
            rbv1.at(6) = -sin(alpha);
            rbv1.at(8) = cos(alpha);
            rbv1.at(1) = rbv1.at(3) = rbv1.at(5) = rbv1.at(7) = 0;
            
            rbv2.at(0) = cos(beta);
            rbv2.at(1) = sin(beta);
            rbv2.at(3) = -sin(beta);
            rbv2.at(4) = cos(beta);
            rbv2.at(8) = 1;
            rbv2.at(2) = rbv2.at(5) = rbv2.at(6) = rbv2.at(7) = 0;
            
            rbv1q.fromRotationMatrix(rbv1);
            rbv2q.fromRotationMatrix(rbv2);
            rbvq = rbv1q;
            rbvq.multRight(rbv2q);
            
            // On doit commencer par les lignes, il me semble
            m11.at(0) = cos(gamma_w);
            m11.at(2) = -sin(gamma_w);
            m11.at(4) = 1;
            m11.at(6) = sin(gamma_w);
            m11.at(8) = cos(gamma_w);
            m11.at(1) = m11.at(3) = m11.at(5) = m11.at(7) = 0;
            
            m12.at(0) = cos(khi_w);
            m12.at(1) = sin(khi_w);
            m12.at(3) = -sin(khi_w);
            m12.at(4) = cos(khi_w);
            m12.at(8) = 1;
            m12.at(2) = m12.at(5) = m12.at(6) = m12.at(7) = 0;
            
            m11q.fromRotationMatrix(m11);
            m12q.fromRotationMatrix(m12);
            mq = m11q;
            mq.multRight(m12q);
            mq.multRight(rviq);
            mq.multRight(rbvq);
            mq.toRotationMatrix(m);
            
            alpha_w = asin(m.at(2));
            
            if(m.at(8) / cos(alpha_w) > 1) {
                sigma_w = 0;
            }else if(m.at(8) / cos(alpha_w) < -1){
                sigma_w = M_1_PI;
            }else{
                sigma_w = sgn(-m.at(5) / cos(alpha_w)) * acos(m.at(8) / cos(alpha_w));
            }
            
            if(m.at(0) / cos(alpha_w) > 1) {
                beta_w = 0;
            }else if(m.at(0) / cos(alpha_w) < -1){
                beta_w = M_1_PI;
            }else{
                beta_w = sgn(m.at(1) / cos(alpha_w)) * acos(m.at(0) / cos(alpha_w));
            }
            
            /** Calc of the aerodynamic force coefficients with wind */
            Cc_w = Cc_beta * beta_w;
            Cl_w = Cl_alpha * (alpha_w - alpha0);
            Cd_w = Cd0 + Cdl * (Cl_w - Clmin) * (Cl_w - Clmin) + Cl_w * Cl_w / (3.14*e*aspect_ratio)
            + Cc_w * Cc_w / (3.14*e*aspect_ratio) * (S / Sv);
            
            /** Calc of the dynamic pressure */
            q = 0.5 * 1.225 * (V_w.at(0)*V_w.at(0) + V_w.at(1)*V_w.at(1) + V_w.at(2)*V_w.at(2));
            
            /** Calc of the aerodynamic forces */
            drag = q * S * Cd_w;
            sideforce = q * S * Cc_w;
            lift = q * S * Cl_w;
            
            //printf("ok\n");
            
        }
    };
    
}

#endif
