#ifndef TRACTOR_PHYSICS_H
#define TRACTOR_PHYSICS_H

#include <ode/ode.h>

/*-------------Algrithm constant-------------*/
#define SPEED_THRES 4.5
#define MS_TO_S_FRAC (1000.0f)
#define TIME_STEP ((current_time - last_time) / MS_TO_S_FRAC)    //s
#define GRAVITY_CONSTANT 9.81

// Soil constant
#define DEFORMATION_COHESION_MODULUS   0.0f // previously 0.0f
#define DEFORMATION_FRICTIONAL_MODULUS 410400.0f  //4.104*10^5
#define SINKAGE_EXP 0.8f
#define CONE_INDEX_CONSTANT 1000000.0  //Previously 1000000.0 If known, just set it here
#define FRICTION_FRACTOR 1000

// Vehicle specification
#define TOTAL_GRAVITY      11575.8f    //N  
#define WHEEL_LOAD_APPROX  2500.0f     //Approximation of the wheel load
#define FRONT_WHEEL_LOAD_INIT    2400.0f
#define REAR_WHEEL_LOAD_INIT     2450.0f
#define FRONT_WHEEL_RADIUS 0.38f
#define FRONT_WHEEL_WIDTH  0.19f
#define REAR_WHEEL_RADIUS  0.6f
#define REAR_WHEEL_WIDTH   0.37f
#define THETA_CONSTANT     0.523598f   // pi/6

/*--------------------------------------------*/

#define EMMITER_RECV_BUF_LENGTH 10
#define EMMITER_BUF_DATA_TYPE   double

//Webots coordinate
//Green -  Y
//Red   -  X
//Blue  -  Z
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
  double driving_force;     //theoretical driving force
  double thoery_slip_ratio; //theoretical slip ratio
  double slip_ratio;        //real slip ratio in simulation
  double wheel_numeric;     
  double towed_force;       //resistance in Wismer-Luth model
  double torque;            //motor torque
  double drawbar_pull;      //friction_force-towed_force
  double trct_effc;         //tractive efficiency
  double actual_distance;   //actual distance the tractor travels taking slip ratio into consideration
  double inertia;           //inertia around y-axis for each wheel
  double friction_coef;     //friction coefficient
  
}mechanics_t;

typedef struct{
  wheel_id_e id;      //ID for wheel, use for debug

  double theta;       //dependent on load 100N
  double angular_vel; //wheel angular velocity
  double horiz_vel;   //wheel horizontal velcity
  double wheel_w;     //wheel width
  double wheel_r;     //wheel radius
  double encode;      //value from motor encoder
  double load;        //normal force form the ground
  double steering_angle;  //only for front wheels
  double friction_force;  //friction force from ground

  dBodyID wheel_body;     //body id used in ODE
  dGeomID wheel_geom;     //geom id used in ODE
  dJointFeedback joint_feed_back;  //used to get force from ground
  dJointID joint;         //wheel-ground collision joint
  dMass mass;             //check out dMass structure in ODE manual

  mechanics_t mechanics;
  calc_mid_value_t calc_mid_value;  //used for medium value in caculation for each wheel 
}wheel_t;

#endif
