#ifndef TRACTOR_PHYSICS_H
#define TRACTOR_PHYSICS_H

#include <ode/ode.h>

/*-------------Algrithm constant-------------*/
#define SPEED_THRES 2
#define MS_TO_S_FRAC (1000.0f)
#define TIME_STEP ((current_time - last_time) / MS_TO_S_FRAC)    //s
#define GRAVITY_CONSTANT 9.81

// Soil constant
#define DEFORMATION_COHESION_MODULUS   0.0f
#define DEFORMATION_FRICTIONAL_MODULUS 410400.0f  //4.104*10^5
#define SINKAGE_EXP 0.8f
#define CONE_INDEX_CONSTANT 1000000.0    //If known, just set it here
#define FRICTION_FRACTOR 1000

// Vehivle specification
#define TOTAL_GRAVITY      11575.8f    //N  
#define WHEEL_LOAD_APPROX  2500.0f     //Approximation of the wheel load
#define FRONT_WHEEL_LOAD_INIT    2400.0f
#define REAR_WHEEL_LOAD_INIT     2600.0f
#define FRONT_WHEEL_RADIUS 0.38f
#define FRONT_WHEEL_WIDTH  0.19f
#define REAR_WHEEL_RADIUS  0.6f
#define REAR_WHEEL_WIDTH   0.37f
#define THETA_CONSTANT     0.523598f   // pi/6

/*--------------------------------------------*/

#define EMMITER_RECV_BUF_LENGTH 10
#define EMMITER_BUF_DATA_TYPE   double

enum{X = 0, Y = 1, Z = 2};

typedef enum{
  FRONT_RIGHT_ID = 0,
  FRONT_LEFT_ID,
  REAR_RIGHT_ID,
  REAR_LEFT_ID
}wheel_id_e;

typedef struct {
  double thr_angular_vel;  //theoretical angular velocity
  double thr_angular_acc;  //theoretical angulat acceleration
  double torque_integral;
  double friction_coef_mid;
  double last_friction_coef_mid;

}calc_mid_value_t;

typedef struct {
  double driving_force;
  double thoery_slip_ratio;
  double slip_ratio;
  double wheel_numeric;
  double towed_force;
  double torque;
  double drawbar_pull;
  double trct_effc;         //tractive efficiency
  double actual_distance;   //actual distance the tractor travels taking slip ratio into consideration
  double inertia;           //inertia around y-axis for each wheel
  double friction_cff;
  double friction;
  
}mechanics_t;

typedef struct{
  wheel_id_e id;      //ID for wheel, can be used to distinguish each wheel when store wheels in vector

  double theta;       //dependent on load 100N
  double angular_vel; //angular velocity
  double horiz_vel;   //horizontal velcity
  double wheel_w;     //wheel width
  double wheel_r;     //wheel radius
  double encode;      //value from motor encoder
  double load;        //normal force form the ground
  double steering_angle;  //only for front wheels

  dBodyID wheel_body;
  dGeomID wheel_geom;
  dJointFeedback joint_feed_back;
  dJointID joint;
  dMass mass;

  mechanics_t mechanics;
  calc_mid_value_t calc_mid_value;  //used for mid value in caculation for each wheel 
}wheel_t;

typedef struct{
    wheel_t frontRightWheel;
    wheel_t frontLeftWheel;
    wheel_t rearRightWheel;
    wheel_t rearLeftWheel;

}wheel_qt_data_t;

const wheel_qt_data_t *get_qt_data(); 

#endif