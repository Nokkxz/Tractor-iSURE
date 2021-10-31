/*
 * File: Bekker_model.cpp
 * Date: 2019.8.19
 * Description: Physics plugin for Tractor's soil-tire interactions.
 * Author: Yonglin Jing
 * Email: 11712605@mail.sustech.edu.cn
 * Edited by: Ranbao Deng
 * Email: 11712002@mail.sustech.edu.cn
 *
 * AgJunction Inc, USA
 * University of Notre Dame, USA
 * Southern University of Science and Technology, PRC
 */ 

/*
 * Note: This plugin will become operational only after it was compiled and associated with the current world (.wbt).
 * To associate this plugin with the world follow these steps:
 *  1. In the Scene Tree, expand the "WorldInfo" node and select its "physics" field
 *  2. Then hit the [Select] button at the bottom of the Scene Tree
 *  3. In the list choose the name of this plugin (same as this file without the extention)
 *  4. Then save the .wbt by hitting the "Save" button in the toolbar of the 3D view
 *  5. Then reload the world: the plugin should now load and execute with the current simulation
 */
 
/*
 * Major Problem for this program:
 * The original formulas of traction_force, resisting_moment, and vertical_force can not be applied on the tractor.
 * Traction force is now applied only on the cab body instead of the wheels
 * Resisting moment is not applied
 * Vertical force might be good actuually?
 */
 
 /*
 * It is noted that the Bekker and WR models were developed
 * under certain assumptions, where a cylindrical wheel moves on
 * a flat and horizontal soil surface under steady-state operations.
 */
 

#include <ode/ode.h>
#include <plugins/physics.h>
#include <math.h>

//Green -  Y
//Red   -  X
//Blue  -  Z

//Using Wismer Luth model

/*
 * What do each of the soil variables mean?
 
 * DEFORMATION_COHESION_MODULUS : k_c [N/m^(n+1)] -- empirical pressure-sinkage parameter for Bekker 
 * higher = higher normal pressure on tire = more motion resistance

 * DEFORMATION_FRICTIONAL_MODULUS : k_theta [N/m^(n+2)] -- empirical pressure-sinkage parameter for Bekker 
 * higher = higher normal pressure on tire = more motion resistance

 * SINKAGE_EXP : n [N/A] -- empirical pressure-sinkage parameter for Bekker 
 * higher = more sinkage = more motion resistance (want smaller sinkage exp, so its not completely linear, but intuitively makes sense)

 * COHESION_STRESS : c [Pa] -- Cohesion is the bond that cements particles of the material together, regardless of the normal pressure between the particles. 
 * high cohesion = stronger soil = less slip
 
 * INTERNAl_SHEARING_RESISTANCE : phi [rad] -- angle of shearing in the soil (20 degrees)
 * higher = higher shear stress
 
 * SHEAR_DEFORMATION_MODULUS : K_d [m] -- soil deformation from shear stress (usually tangent to wheel surface)
 * higher = higher shear stress = lower traction ??
 
 * FRICTION_FACTOR : ?? not used anywhere
 
 * NDV = number of divisions for for loops in motion resistance stuff
 
 * MU = test coefficient of friction - this is probably not the correct way to do it, but we couldn't figure it out otherwise
*/

/*-------------Soil constants-------------*/
// Uncomment whichever soil type you want to use. 
// If different soil parameters are needed for different wheels, they should be put in the struct instead.

//#define SNOW
#define SANDY_LOAM_MICHIGAN
//#define BEKKER_ARTICLE_SOIL
//#define HEAVY_CLAY
//#define PARKING_LOT
// you have to change the set torque for reasonable results (for horizontal test)

//record the data or not
#define RECORD

#ifdef RECORD
  FILE * file = NULL;
  // file can be renamed here to fit MATLAB code, but also must be changed in Bekker
  char * path = "file_circle.txt";
#endif

//add lateral force functions into it
#define ADD_LATERAL_FORCE

#ifdef PARKING_LOT
  // USE TORQUE OF 70 NM FOR HORIZONTAL TEST
  // SET RADIUS OF 12.5 M FOR CIRCLE TEST
  // Parking lot (for reference)
  #define DEFORMATION_COHESION_MODULUS   0.0f
  #define DEFORMATION_FRICTIONAL_MODULUS 112797000.0f
  #define SINKAGE_EXP 0.0f
  #define COHESION_STRESS 0.0f
  #define INTERNAL_SHEARING_RESISTANCE 0.0f
  #define SHEAR_DEFORMATION_MODULUS 0.0f //estimated
  #define FRICTION_FACTOR 1000
  #define Ndv 10
  #define MU dInfinity
#endif  
#ifdef SANDY_LOAM_MICHIGAN
  // USE TORQUE OF 192 NM AND MU = 502.25 FOR HORIZONTAL TEST
  // SET RADIUS OF 12.5 M FOR CIRCLE TEST
  // Soil constant: sandy loam michigan, 
  #define DEFORMATION_COHESION_MODULUS   52230.0f
  #define DEFORMATION_FRICTIONAL_MODULUS 1127970.0f
  #define SINKAGE_EXP 0.9f
  #define COHESION_STRESS 4830.0f
  #define INTERNAL_SHEARING_RESISTANCE 0.3489f
  #define SHEAR_DEFORMATION_MODULUS 0.02f //estimated
  #define FRICTION_FACTOR 1000
  #define Ndv 10
  #define MU 502.25
  //#define MU dInfinity
#endif  
#ifdef HEAVY_CLAY
  // USE TORQUE OF 70 NM AND MU = 227 FOR HORIZONTAL TEST
  // SET RADIUS OF 12.5 M FOR CIRCLE TEST
  // Heavy clay (hard soil)
  #define DEFORMATION_COHESION_MODULUS   12700.0f
  #define DEFORMATION_FRICTIONAL_MODULUS 1555950.0f
  #define SINKAGE_EXP 0.13f
  #define COHESION_STRESS 68950.0f
  // (34 degrees)
  #define INTERNAL_SHEARING_RESISTANCE 0.5934f
  #define SHEAR_DEFORMATION_MODULUS 0.02f //estimated
  #define FRICTION_FACTOR 1000
  #define Ndv 10
  //#define MU 227
  #define MU dInfinity
#endif  
#ifdef BEKKER_ARTICLE_SOIL
  // USE TORQUE OF 250 NM AND MU = 654.75 FOR HORIZONTAL TEST
  // SET RADIUS OF 12.5 M FOR CIRCLE TEST
  // Bekker soil (from article)
  #define DEFORMATION_COHESION_MODULUS   0.0f
  #define DEFORMATION_FRICTIONAL_MODULUS 410400.0f
  #define SINKAGE_EXP 0.8f
  #define COHESION_STRESS 234.0f
  // 30 deg
  #define INTERNAL_SHEARING_RESISTANCE 0.5236f
  #define SHEAR_DEFORMATION_MODULUS 0.013f //estimated
  #define FRICTION_FACTOR 1000
  #define Ndv 10
  //#define MU 654.75
  #define MU dInfinity
#endif   
#ifdef SNOW
  // This causes super extreme cases that are probably just too unrealistic
  // Snow (really light)
  #define DEFORMATION_COHESION_MODULUS   10550.0f
  #define DEFORMATION_FRICTIONAL_MODULUS 66080.0f
  #define SINKAGE_EXP 1.44f
  #define COHESION_STRESS 6000.0f
  // (20.7  degrees)
  #define INTERNAL_SHEARING_RESISTANCE 0.3613f
  #define SHEAR_DEFORMATION_MODULUS 0.02f //estimated
  #define FRICTION_FACTOR 1000
  #define Ndv 10
  #define MU 5000
#endif 

      
// Vehicle specification

#define VEHICLE_GRAVITY    10000.0f //N  Modified by NERanger 20190807
#define TOTAL_MASS         940.0f   //kg 
#define MOMENT_OF_INERTIA  2557.88f   //kg*m^2 I assumed this was around the y-axis
#define WHEEL_LOAD_APPROX  2500.0f  //Approximation of the wheel load
#define WHEEL_LOAD_INIT    2500.0f
#define FRONT_WHEEL_RADIUS 0.38f
#define FRONT_WHEEL_WIDTH  0.19f
#define REAR_WHEEL_RADIUS  0.6f
#define REAR_WHEEL_WIDTH   0.37f


// Vehicle specification (to test against MATLAB Bekker)
/*
#define VEHICLE_GRAVITY    400.0f //N  Modified by NERanger 20190807
#define TOTAL_MASS         40.0f   //kg 
#define MOMENT_OF_INERTIA  2557.88f   //kg*m^2 I assumed this was around the y-axis
#define WHEEL_LOAD_APPROX  100.0f  //Approximation of the wheel load
#define WHEEL_LOAD_INIT    100.0f
#define FRONT_WHEEL_RADIUS 0.15f
#define FRONT_WHEEL_WIDTH  0.11f
#define REAR_WHEEL_RADIUS  0.15f
#define REAR_WHEEL_WIDTH   0.1f
*/

typedef enum{
  FRONT_RIGHT_ID = 0,
  FRONT_LEFT_ID,
  REAR_RIGHT_ID,
  REAR_LEFT_ID
}wheel_id_e;

typedef struct {
//the following parameters can be found in the week 3 slides. more details are shown there.
  double sinkage;                //sinkage of the wheel based on radial coodinates.
  double normal_stress;          //based on radial coodinates.
  double shear_stress;           //same as above.
  double shear_displacement;     //same as above.
  double vertical_force;         //Verification for load. Theoretically, this should equal load.
  double max_sinkage;            //sinkage at bottom, where theta is zero.
  double theta_1;                //exit angle of the wheel
  double traction_force;
  double motion_resistance;
  double motion_resistance_2;    //used to try to solve the integration force problem
  double motion_resistance_fast; //a transfromed equation was found to describe motion resistance without integration.
  double resisting_moment;
  double driving_force;          //Wismer-Luth model. used to compare results.
  double theory_slip_ratio;      //theoretical value of slip ratio. calculated by a transformed formula using Wismer-Luth model.
  double slip_ratio;             //definition formula of slip ratio.
  double wheel_numeric;          //Wismer-Luth model. used to compare results.
  double torque;                 //torque given by the engine.
  double drawbar_pull;
  double trct_effc;   //tractive efficiency, Wismer-Luth model. used to compare results.
  double lateral_accel;  // lateral acceleration
  double angular_accel_yaxis; // angular acceleration
}mechanics_t;

typedef struct{
  wheel_id_e id;      //ID for wheel, can be used to distinguish each wheel when store wheels in vector
  double angular_vel; //angular velocity
  double horiz_vel;   //horizontal velcity
  double wheel_w;     //wheel width
  double wheel_r;     //wheel radius
  double load;        //normal force form the ground
  double steering_angle;  //only for front wheels

  dBodyID wheel_body;
  dGeomID wheel_geom;
  dJointFeedback joint_feed_back;
  dJointID joint;

  mechanics_t mechanics;
}wheel_t;

//std::vector<wheel_t> wheels;
//static double friction_frac = 0.0;

static wheel_t front_right_w, front_left_w, rear_right_w, rear_left_w;

static dBodyID tractorBody = NULL;
static dBodyID floorBody = NULL;
static dGeomID floorGeom = NULL;

static pthread_mutex_t mutex;  // needed to run with multi-threaded version of ODE

inline double sinkage_calc(double r,   //wheel radius
                           double b,   //wheel width
                           double L,   //load, vertical force of the wheel
                           double t,   //theta,as location
                           double kc,  //cohesion modulus of deformation
                           double se,  //sinkage exponent
                           double kp){ //frictional modulus of deformation
	double zm = 0.0;
	zm = pow((3 * L / (b * (3 - se) * (kc / b + kp) * pow(2 * r, 0.5))), (2 / (2 * se + 1)));
	
	return (r * (cos(t) - 1) + zm);
}

inline double cone_index_calc(double r,   //wheel radius
	double kc,  //cohesion modulus of deformation
	double se,  //sinkage exponent
	double kp,  //frictional modulus of deformation
	double b,   //wheel width
	double L) { //wheel load
	double zetam = 0.0;
	zetam = sinkage_calc(r,b,L,0,kc,se,kp);

	return 1.625 * (kc / (se + 1) * (pow((zetam + 1.5), (se + 1)) - pow(zetam, (se + 1))) + 0.517 * kp * (pow((zetam + 1.5), \
		(se + 2)) / (se + 1) / (se + 2) + pow(zetam, (se + 2)) / (se + 2) - (zetam + 1.5) * pow(zetam, (se + 1)) / (se + 1)));
}


inline double wheel_numeric_calc(double b,  //wheel width 
	double r,  //wheel radius
	double L) { //wheel load
	return cone_index_calc(r, DEFORMATION_COHESION_MODULUS, SINKAGE_EXP, DEFORMATION_FRICTIONAL_MODULUS, b, L) * b * 2 * r / L;
}

inline double torque_calc(double r,   //wheel radius
	double L,   //wheel load
	double wn,  //wheel numeric
	double sr) { //slip ratio

	return r * L * 0.75 * (1 - exp(-0.3 * wn * sr));
}


inline double normal_stress_calc(double b,   //wheel width
                                 double z,   //sinkage
                                 double kc,  //cohesion modulus of deformation
                                 double se,  //sinkage exponent
                                 double kp){ //frictional modulus of deformation
	        	

      if(z<0)
        z=0;

	return (kc / b + kp) * pow(z, se);
	
}

inline double shear_stress_calc(double c,    //cohesion stress
                                double phi,  //internal shearing resistance
	                            double tau,  //normal stress distribution
	                            double jd,   //shear displacement
	                            double kd){ //shear deformation modulus
	return ((c + tau * tan(phi)) * (1 - exp(-jd / kd)));
}

inline double slip_ratio_calc(double r,  //wheel radius 
                              double w,  //wheel angular velocity
                              double vx){  //wheel horizontal velocity
	return (1 - (vx / (r * w)));
}

inline double theory_slip_ratio_calc(double wn,    //wheel numeric
	double t,     //torque
	double r,     //wheel radius
	double L) {   //wheel load

	return (-10 / (3 * wn)) * log(1 - (4 * t) / (3 * r * L));
}

inline double shear_displacement_calc(double r,   //wheel radius 
                                      double t,   //theta,as location
                                      double t1,  //maximum sinkage
                                      double sr){ //slip ratio
	return r * ((t1 - t) - (1 - sr) * (sin(t1) - sin(t)));
}

inline double motion_resistance_calc_fast(double b,    //wheel radius
                                          double zm,   //maximum sinkage
                                          double kc,   //cohesion modulus of deformation
                                          double se,   //sinkage exponent
                                          double kp) { //frictional modulus of deformation
	return b * ((kc / b + kp) * pow(zm, (se + 1)) / (se + 1));
}

inline double motion_resistance_calc(double r,   //wheel radius
	double b,   //wheel width
	double L,   //load
	double kc,  //cohesion modulus of deformation
	double se,  //sinkage exponent
	double kp,  //frictional modulus of deformation
	double t1) {//theta_1
	double result = 0;
	double delta = t1 / Ndv;
	
	for (double i = 0 + delta; i < t1; i += delta)
	{
		result += normal_stress_calc(b, sinkage_calc(r, b, L, 0 + i, kc, se, kp), kc, se, kp) * sin(0 + i) * delta;
	}
	
	return (result * r * b);
}

inline double resisting_moment_calc(double r,   //wheel radius
	                                double b,   //wheel width
	                                double L,   //load
	                                double kc,  //cohesion modulus of deformation
	                                double se,  //sinkage exponent
	                                double kp,  //frictional modulus of deformation
	                                double c,   //cohesion stress
      				           double phi,  //internal shearing resistance
				           double kd,   //shear deformation modulus
  					    double w,    //angular velocity
					    double vx,   //horizontal velocity
				           double t1) { //theta_1
	double result = 0;
	double delta = t1 / Ndv;
	double sr = slip_ratio_calc(r, w, vx);
	for (double i = 0 + delta; i < t1; i += delta)
	{
		result += shear_stress_calc(c, phi, normal_stress_calc(b, sinkage_calc(r, b, L, 0 + i, kc, se, kp), kc, se, kp), shear_displacement_calc(r, 0 + i, t1, sr), kd) * delta;
	}
	return result*r*r*b;
}

inline double traction_force_calc(double r,   //wheel radius
	double b,   //wheel width
	double L,   //load
	double kc,  //cohesion modulus of deformation
	double se,  //sinkage exponent
	double c,   //cohesion stress
	double phi, //internal shearing resistance
	double kd,  //shear deformation modulus
	double kp,  //frictional modulus of deformation
	double w,   //angular velocity
	double vx,  //horizontal velocity
	double t1,  //theta_1
	double sr) {//slip ratio
	
	double result = 0;
	double delta = t1 / Ndv;
	for (double i = 0 + delta; i < t1; i += delta)
	{
		result += shear_stress_calc(c, phi, normal_stress_calc(b, sinkage_calc(r, b, L, 0 + i, kc, se, kp), kc, se, kp), shear_displacement_calc(r, 0 + i, t1, sr), kd) * cos(0 + i) * delta;
	}
	return result * r * b;
}

inline double vertical_force_calc(double r,   //wheel radius
	double b,   //wheel width
	double L,   //load
	double kc,  //cohesion modulus of deformation
	double se,  //sinkage exponent
	double c,   //cohesion stress
	double phi, //internal shearing resistance
	double kd,  //shear deformation modulus
	double kp,  //frictional modulus of deformation
	double w,   //angular velocity
	double vx,  //horizontal velocity
	double t1,  //theta_1
	double sr) {//slip ratio
	
	double result = 0;
	double delta = t1 / Ndv;
	for (double i = 0 + delta; i < t1; i += delta)
	{
		result += (normal_stress_calc(b, sinkage_calc(r, b, L, 0 + i, kc, se, kp), kc, se, kp) * cos(0 + i) + shear_stress_calc(c, phi, normal_stress_calc(b, sinkage_calc(r, b, L, 0 + i, kc, se, kp), kc, se, kp), shear_displacement_calc(r, 0 + i, t1, sr), kd) * sin(0 + i)) * delta;
	}
	return result * r * b;
}

inline double drawbar_pull_calc(double ft,   //traction force
                                double rc){ //motion resistance
  return ft-rc;
}

inline double tractive_efficiency_calc(double dp,  //drawbar pull
                                       double vx,  //horizontal velocity
                                       double t,   //torque
                                       double w){  //angular velocity
  return dp*vx/t/w;
}

inline double wheel_driving_force_calc(double t,     //motor torque
                                       double r){    //wheel radius
  return (t*r);
}

inline double lateral_accel_calc(double m, // mass of tractor
                                 double Caf, // front tire cornering stiffness coefficient 
                                 double Carear, // rear tire cornering stiffness coefficient 
                                 double u, // z (forward) velocity
                                 double v, // x (side) velocity
                                 double a, // distance from front wheels to CoM
                                 double b, // distance from rear wheels to CoM
                                 double y, // angular velocity around y-axis
                                 double delta){ // slip angle -- angle between the current direction of travel and the direction the front wheel faces                         
  
  double v_prime;  
  v_prime = (-Caf * (atan2(v+(a*y), u) - delta) - Carear * atan2(v-(b*y), u) - (m*u*y)) / m;
  return v_prime;
}

inline double angular_accel_yaxis_calc(double I, // moment of inertia around CoM
                                       double Caf, // front tire cornering stiffness coefficient 
                                       double Carear, // rear tire cornering stiffness coefficient 
                                       double u, // z (forward) velocity
                                       double v, // x (side) velocity
                                       double a, // distance from front wheels to CoM
                                       double b, // distance from rear wheels to CoM
                                       double y, // angular velocity around y-axis
                                       double delta){ // slip angle -- angle between the current direction of travel and the direction the front wheel faces
                     
  double y_prime;                   
  y_prime = (-a*Caf*(atan2(v+(a*y), u) - delta) + b*Carear*atan2(v-(b*y), u)) / I;
  return y_prime;
}


void wheel_parameter_init() {
  front_right_w.id      = FRONT_RIGHT_ID;
  front_right_w.wheel_r = FRONT_WHEEL_RADIUS;
  front_right_w.wheel_w = FRONT_WHEEL_WIDTH;
  //front_right_w.load    = WHEEL_LOAD_INIT;

  front_left_w.id       = FRONT_LEFT_ID;
  front_left_w.wheel_r  = FRONT_WHEEL_RADIUS;
  front_left_w.wheel_w  = FRONT_WHEEL_WIDTH;
  //front_left_w.load     = WHEEL_LOAD_INIT;
  
  rear_right_w.id       = REAR_RIGHT_ID;
  rear_right_w.wheel_r  = REAR_WHEEL_RADIUS;
  rear_right_w.wheel_w  = REAR_WHEEL_WIDTH;
  //rear_right_w.load     = WHEEL_LOAD_INIT;
  
  rear_left_w.id        = REAR_LEFT_ID;
  rear_left_w.wheel_r   = REAR_WHEEL_RADIUS;
  rear_left_w.wheel_w   = REAR_WHEEL_WIDTH;
  //rear_left_w.load      = WHEEL_LOAD_INIT;


   /*
   * TFRW - Tractor Front Right Wheel
   * TFLW - Tractor Front Left  Wheel
   * TRRW - Tractor Rear  Right Wheel
   * TRLW - Tractor Rear  Left  Wheel
   */

  front_right_w.wheel_body = dWebotsGetBodyFromDEF("TFRW");
  front_right_w.wheel_geom = dWebotsGetGeomFromDEF("TFRW");

  if (front_right_w.wheel_body == NULL){
    dWebotsConsolePrintf("!!! error : could not get frontRightWheelBody.\r\n");
  }else if (front_right_w.wheel_geom == NULL){
    dWebotsConsolePrintf("!!! error : could not get frontRightWheelGeom.\r\n");
  }else {
    dWebotsConsolePrintf("frontRightWheel Body&Geom found !\r\n");
  }

  front_left_w.wheel_body = dWebotsGetBodyFromDEF("TFLW");
  front_left_w.wheel_geom = dWebotsGetGeomFromDEF("TFLW");

  if (front_left_w.wheel_body == NULL){
    dWebotsConsolePrintf("!!! error : could not get frontLeftWheelBody.\r\n");
  }else if(front_left_w.wheel_geom == NULL){
    dWebotsConsolePrintf("!!! error : could not get frontLeftWheelGeom.\r\n");
  }else {
    dWebotsConsolePrintf("frontLeftWheel Body&Geom found !\r\n");
  }
  
  rear_right_w.wheel_body = dWebotsGetBodyFromDEF("TRRW");
  rear_right_w.wheel_geom = dWebotsGetGeomFromDEF("TRRW");

  if (rear_right_w.wheel_body == NULL){
    dWebotsConsolePrintf("!!! error : could not get rearRightWheelBody.\r\n");
  }else if(rear_right_w.wheel_geom == NULL){
    dWebotsConsolePrintf("!!! error : could not get rearRightWheelGeom.\r\n");
  }else {
    dWebotsConsolePrintf("rearRightWheel Body&Geom found !\r\n");
  }
  
  rear_left_w.wheel_body = dWebotsGetBodyFromDEF("TRLW");
  rear_left_w.wheel_geom = dWebotsGetGeomFromDEF("TRLW");

  if (rear_left_w.wheel_body == NULL){
    dWebotsConsolePrintf("!!! error : could not get rearLeftWheelBody.\r\n");
  }else if(rear_left_w.wheel_geom == NULL){
    dWebotsConsolePrintf("!!! error : could not get rearLeftWheelGeom.\r\n");
  }else {
    dWebotsConsolePrintf("rearLeftWheel Body&Geom found !\r\n");
  }
}

void wheel_mechanics_calc(wheel_t* wheel){
  //Note that calculation order should not be changed

  wheel->mechanics.slip_ratio = slip_ratio_calc(wheel->wheel_r, 
                                                wheel->angular_vel, 
                                                wheel->horiz_vel);

  wheel->mechanics.driving_force = wheel_driving_force_calc(wheel->mechanics.torque,
	  wheel->wheel_r);

  wheel->mechanics.wheel_numeric = wheel_numeric_calc(wheel->wheel_w,
	  wheel->wheel_r,
	  wheel->load);

  wheel->mechanics.theory_slip_ratio = theory_slip_ratio_calc(wheel->mechanics.wheel_numeric,
	  wheel->mechanics.torque,
	  wheel->wheel_r,
	  wheel->load);
  
  wheel->mechanics.max_sinkage = sinkage_calc(wheel->wheel_r, 
	  wheel->wheel_w, 
	  wheel->load, 0, 
	  DEFORMATION_COHESION_MODULUS, 
	  SINKAGE_EXP, 
	  DEFORMATION_FRICTIONAL_MODULUS);

  wheel->mechanics.theta_1 = acos(1 - wheel->mechanics.max_sinkage / wheel->wheel_r);

  wheel->mechanics.motion_resistance_fast = motion_resistance_calc_fast(wheel->wheel_w, 
	  sinkage_calc(wheel->wheel_r, 
		  wheel->wheel_w, 
		  wheel->load, 
		  0, 
		  DEFORMATION_COHESION_MODULUS, 
		  SINKAGE_EXP, 
		  DEFORMATION_FRICTIONAL_MODULUS), 
	  DEFORMATION_COHESION_MODULUS, 
	  SINKAGE_EXP, 
	  DEFORMATION_FRICTIONAL_MODULUS);

  wheel->mechanics.motion_resistance = motion_resistance_calc(wheel->wheel_r,
	  wheel->wheel_w,
	  wheel->load,
	  DEFORMATION_COHESION_MODULUS,
	  SINKAGE_EXP,
	  DEFORMATION_FRICTIONAL_MODULUS,
	  wheel->mechanics.theta_1);
	  
      //dWebotsConsolePrintf("motion resistance: \t %f\n", wheel->mechanics.motion_resistance);


  wheel->mechanics.resisting_moment = resisting_moment_calc(wheel->wheel_r,
	  wheel->wheel_w,
	  wheel->load,
	  DEFORMATION_COHESION_MODULUS,
	  SINKAGE_EXP,
	  DEFORMATION_FRICTIONAL_MODULUS,
	  COHESION_STRESS,
	  INTERNAL_SHEARING_RESISTANCE,
	  SHEAR_DEFORMATION_MODULUS,
	  wheel->angular_vel,
	  wheel->horiz_vel,
	  wheel->mechanics.theta_1);
  
  wheel->mechanics.traction_force = traction_force_calc(wheel->wheel_r,
	  wheel->wheel_w,
	  wheel->load,
	  DEFORMATION_COHESION_MODULUS,
	  SINKAGE_EXP, COHESION_STRESS, INTERNAL_SHEARING_RESISTANCE, SHEAR_DEFORMATION_MODULUS, DEFORMATION_FRICTIONAL_MODULUS,
	  wheel->angular_vel,
	  wheel->horiz_vel,
	  wheel->mechanics.theta_1,
	  wheel->mechanics.theory_slip_ratio);

  wheel->mechanics.vertical_force = vertical_force_calc(wheel->wheel_r,
	  wheel->wheel_w,
	  wheel->load,
	  DEFORMATION_COHESION_MODULUS,
	  SINKAGE_EXP,
	  COHESION_STRESS,
	  INTERNAL_SHEARING_RESISTANCE,
	  SHEAR_DEFORMATION_MODULUS,
	  DEFORMATION_FRICTIONAL_MODULUS,
	  wheel->angular_vel,
	  wheel->horiz_vel,
	  wheel->mechanics.theta_1,
	  wheel->mechanics.theory_slip_ratio);

  //Terminate motion resistance in case that motion resistance will cause tractor move backwards
  /*
  // This was commented out so you could go into reverse if you wanted to
  if(wheel->horiz_vel < 0.01){
    wheel->mechanics.motion_resistance_fast = 0;
	wheel->mechanics.motion_resistance = 0;
    dWebotsConsolePrintf("motion resistance terminated\n");
  }
  */
  wheel->mechanics.drawbar_pull = drawbar_pull_calc(wheel->mechanics.traction_force, 
                                                    wheel->mechanics.motion_resistance);

  wheel->mechanics.trct_effc = tractive_efficiency_calc(wheel->mechanics.drawbar_pull, 
                                                        wheel->horiz_vel, 
                                                        wheel->mechanics.torque, 
                                                        wheel->angular_vel);

       // this eq is only for front wheels not back wheels I believe	
	const dReal *tmp_world;
	dVector3 side_vel;
       //tmp_world = dBodyGetLinearVel(front_right_w.wheel_body); // do we use tractor direction or wheels?
       tmp_world = dBodyGetLinearVel(tractorBody);
       dBodyVectorFromWorld(tractorBody, tmp_world[0], tmp_world[1], tmp_world[2], side_vel);
	
	const dReal * tmp_angular;
	tmp_angular = dBodyGetAngularVel(tractorBody);	

   // I am pretty sure these variables were entered correctly, but are definitely worth checkng over
  wheel->mechanics.lateral_accel = lateral_accel_calc(TOTAL_MASS, //kg -- mass of tractor  + 4 wheels (from PROTO)
                                                      0.15*wheel->load, //16000, //N/rad -- estimated from article (COULD BE SO WRONG)
                                                      0.15*wheel->load, //16000, //N/rad -- estimated from article (COULD BE SO WRONG)
                                                      wheel->horiz_vel,
                                                      side_vel[0], // x (side )velocity -- needs to be fixed?
                                                      0.861931, // 1.807 - 0.945069 m
                                                      0.870069,// |0.075 - 0.945069| m
                                                      tmp_angular[1], // angular velocity around y-axis   
                                                      wheel->steering_angle);
       
  wheel->mechanics.angular_accel_yaxis = angular_accel_yaxis_calc(MOMENT_OF_INERTIA, //kg -- mass of tractor  + 4 wheels (from PROTO)
                                                      0.15*wheel->load, //16000, //N/rad -- estimated from article (COULD BE SO WRONG)
                                                      0.15*wheel->load, //16000, //N/rad -- estimated from article (COULD BE SO WRONG)
                                                      wheel->horiz_vel,
                                                      side_vel[0], // x (side )velocity -- needs to be fixed?
                                                      0.861931, // 1.807 - 0.945069 m
                                                      0.870069,// |0.075 - 0.945069| m
                                                      tmp_angular[1], // angular velocity around y-axis   
                                                      wheel->steering_angle);
       //dWebotsConsolePrintf("side_vel: \t %f \n", side_vel[0] );
                      
}

double wheel_mu_calc(wheel_t* wheel){
  dWebotsConsolePrintf("Order: RL, FR, FL, RR\n");
  dWebotsConsolePrintf("mu: \t %f \n", abs(wheel->mechanics.torque / wheel->wheel_r - wheel->mechanics.motion_resistance) / wheel->load);
  dWebotsConsolePrintf("load: \t %f \n", wheel->load);

  return (wheel->mechanics.motion_resistance / wheel->load);
  // mu = F/N
}

void wheel_steering_angle_update(double* recv_buf) {
	front_right_w.steering_angle = recv_buf[8];
	front_left_w.steering_angle = recv_buf[9];
}

void wheel_load_update(){
  front_right_w.load = front_right_w.joint_feed_back.f1[1];
  front_left_w.load  = front_left_w.joint_feed_back.f1[1];
  rear_right_w.load  = rear_right_w.joint_feed_back.f1[1];
  rear_left_w.load   = rear_left_w.joint_feed_back.f1[1];
}

void wheel_force_unit_vector(dVector3 result) {
	dVector3 tmp, tmp1;
	dBodyVectorToWorld(tractorBody, 0, sin(0.125), cos(0.125), tmp);
	dBodyVectorToWorld(tractorBody, 0, sin(0.125), cos(0.125), tmp1);

	tmp[0] = tmp1[0] / sqrt(tmp1[2] * tmp1[2] + tmp1[1] * tmp1[1] + tmp1[0] * tmp1[0]);
	tmp[1] = tmp1[1] / sqrt(tmp1[2] * tmp1[2] + tmp1[1] * tmp1[1] + tmp1[0] * tmp1[0]);
	tmp[2] = tmp1[2] / sqrt(tmp1[2] * tmp1[2] + tmp1[1] * tmp1[1] + tmp1[0] * tmp1[0]);
	result[0] = tmp[0];
	result[1] = tmp[1];
	result[2] = tmp[2];
	//dWebotsConsolePrintf("tractor direction: %f,%f,%f\n",result[0], result[1], result[2]);	
}

void steer_unit_vector(double t, dVector3 result) {
	dVector3 tmp, tmp1;
	//front_right_w.wheel_body
	dBodyVectorToWorld(tractorBody, cos(0.125) * sin(t), sin(0.125), cos(0.125) * cos(t), tmp);
	dBodyVectorToWorld(tractorBody, cos(0.125) * sin(t), sin(0.125), cos(0.125) * cos(t), tmp1);

	tmp[0] = tmp1[0] / sqrt(tmp1[2] * tmp1[2] + tmp1[1] * tmp1[1] + tmp1[0] * tmp1[0]);
	tmp[1] = tmp1[1] / sqrt(tmp1[2] * tmp1[2] + tmp1[1] * tmp1[1] + tmp1[0] * tmp1[0]);
	tmp[2] = tmp1[2] / sqrt(tmp1[2] * tmp1[2] + tmp1[1] * tmp1[1] + tmp1[0] * tmp1[0]);
	result[0] = tmp[0];
	result[1] = tmp[1];
	result[2] = tmp[2];
	//dWebotsConsolePrintf("steer direction: %f,%f,%f\n",result[0], result[1], result[2]);
}

void wheel_motor_torque_update(double* torque_feedback){

  front_left_w.mechanics.torque  = torque_feedback[0];
  front_right_w.mechanics.torque = torque_feedback[1];
  rear_left_w.mechanics.torque   = torque_feedback[2];
  rear_right_w.mechanics.torque  = torque_feedback[3];

  //friction_frac = torque_feedback[0];

  //dWebotsConsolePrintf("motor torque update: %f %f %f %f", torque_feedback[0], torque_feedback[1], torque_feedback[2], torque_feedback[3]);
}

void wheel_angular_vel_update(){
  const dReal* tmp;

  //coordinate transform required if consider steering
  // uses wheels coordinate frame
  tmp = dBodyGetAngularVel(front_right_w.wheel_body);
  front_right_w.angular_vel = sqrt(tmp[0]*tmp[0] + tmp[2]*tmp[2]);
  tmp = dBodyGetAngularVel(front_left_w.wheel_body);
  front_left_w.angular_vel = sqrt(tmp[0]*tmp[0] + tmp[2]*tmp[2]);
  tmp = dBodyGetAngularVel(rear_right_w.wheel_body);
  rear_right_w.angular_vel = sqrt(tmp[0]*tmp[0] + tmp[2]*tmp[2]);
  tmp = dBodyGetAngularVel(rear_left_w.wheel_body);
  rear_left_w.angular_vel = sqrt(tmp[0]*tmp[0] + tmp[2]*tmp[2]);
}

void wheel_linear_vel_update(){
  const dReal* tmp;

  //coordinate transform required if consider steering
  tmp = dBodyGetLinearVel(front_right_w.wheel_body);
  front_right_w.horiz_vel = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
  tmp = dBodyGetLinearVel(front_left_w.wheel_body);
  front_left_w.horiz_vel = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
  tmp = dBodyGetLinearVel(rear_right_w.wheel_body);
  rear_right_w.horiz_vel = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
  tmp = dBodyGetLinearVel(rear_left_w.wheel_body);
  rear_left_w.horiz_vel = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
}

void print_output(){
/*
	dWebotsConsolePrintf("Bekker-Front_Right(%d): slip ratio:%f, thry slip ratio:%f, wheel numeric:%f, tractive eficiency:%f, drawbar pull:%f, motion resistance:%f, torque:%f, w:%f, vx:%f, traction force:%f\nresisting moment:%f, vertical force:%f, load:%f, motion resistance fast:%f, max sinkage:%f, max shear stress:%f, max normal stress:%f, max shear displacement:%f\ntheta_1:%f, fdir:[%f, %f, %f], wheel numeric:%f\n",
		Ndv, front_right_w.mechanics.slip_ratio, front_right_w.mechanics.theory_slip_ratio, front_right_w.mechanics.wheel_numeric, front_right_w.mechanics.trct_effc,
		front_right_w.mechanics.drawbar_pull, front_right_w.mechanics.motion_resistance, front_right_w.mechanics.torque,
		front_right_w.angular_vel, front_right_w.horiz_vel, front_right_w.mechanics.traction_force, front_right_w.mechanics.resisting_moment,
		front_right_w.mechanics.vertical_force,
		front_right_w.load,
		front_right_w.mechanics.motion_resistance_fast,
		front_right_w.mechanics.max_sinkage,
		shear_stress_calc(COHESION_STRESS, INTERNAL_SHEARING_RESISTANCE, normal_stress_calc(front_right_w.wheel_w, front_right_w.mechanics.max_sinkage, DEFORMATION_COHESION_MODULUS, SINKAGE_EXP, DEFORMATION_FRICTIONAL_MODULUS), shear_displacement_calc(front_right_w.wheel_r, 0, front_right_w.mechanics.max_sinkage, front_right_w.mechanics.theory_slip_ratio), SHEAR_DEFORMATION_MODULUS),
		normal_stress_calc(front_right_w.wheel_w, front_right_w.mechanics.max_sinkage, DEFORMATION_COHESION_MODULUS, SINKAGE_EXP, DEFORMATION_FRICTIONAL_MODULUS),
		shear_displacement_calc(front_right_w.wheel_r, 0, front_right_w.mechanics.max_sinkage, front_right_w.mechanics.theory_slip_ratio),
		front_right_w.mechanics.theta_1,
		front_right_w.joint_feed_back.f1[0] + front_left_w.joint_feed_back.f1[0] + rear_right_w.joint_feed_back.f1[0] + rear_left_w.joint_feed_back.f1[0],
		front_right_w.joint_feed_back.f1[1] + front_left_w.joint_feed_back.f1[1] + rear_right_w.joint_feed_back.f1[1] + rear_left_w.joint_feed_back.f1[1],
		front_right_w.joint_feed_back.f1[2] + front_left_w.joint_feed_back.f1[2] + rear_right_w.joint_feed_back.f1[2] + rear_left_w.joint_feed_back.f1[2],
		front_right_w.mechanics.wheel_numeric);
	*/
  
  
  //dWebotsConsolePrintf("FR angular_vel: \t %f \n", front_right_w.angular_vel);
  //dWebotsConsolePrintf("FR slip_ratio: \t %f\n", front_right_w.mechanics.slip_ratio);
  //dWebotsConsolePrintf("BR angular_vel: \t %f \n", rear_right_w.angular_vel);
  //dWebotsConsolePrintf("BR slip_ratio: \t %f\n", rear_right_w.mechanics.slip_ratio);
  //dWebotsConsolePrintf("theory_sr: \t %f\n",front_right_w.mechanics.theory_slip_ratio);
  //dWebotsConsolePrintf("torque (Nm): \t %f\n",front_right_w.mechanics.torque);
  //dWebotsConsolePrintf("lateral accel: \t %f\n",front_right_w.mechanics.lateral_accel);
  //dWebotsConsolePrintf("angular accel: \t %f\n",front_right_w.mechanics.angular_accel_yaxis);
  //dWebotsConsolePrintf("motion resistance: \t %f\n",front_right_w.mechanics.motion_resistance);
  //dWebotsConsolePrintf("normal stress: \t %f\n",normal_stress_calc(front_right_w.wheel_w, front_right_w.mechanics.max_sinkage, DEFORMATION_COHESION_MODULUS, SINKAGE_EXP, DEFORMATION_FRICTIONAL_MODULUS));
  //dWebotsConsolePrintf("sinkage: \t %f\n",front_right_w.mechanics.max_sinkage);
  //dWebotsConsolePrintf("traction force: \t %f\n",front_right_w.mechanics.traction_force);
  //dWebotsConsolePrintf("resisting moment: \t %f\n",front_right_w.mechanics.resisting_moment);
  //dWebotsConsolePrintf("vertical load: \t %f\n",front_right_w.mechanics.vertical_force);

}

void webots_physics_init() {
  pthread_mutex_init(&mutex, NULL);
  /*
   * Get ODE object from the .wbt model, e.g.
   *   dBodyID body1 = dWebotsGetBodyFromDEF("MY_ROBOT");
   *   dBodyID body2 = dWebotsGetBodyFromDEF("MY_SOLID");
   *   dGeomID geom2 = dWebotsGetGeomFromDEF("MY_SOLID");
   * If an object is not found in the .wbt world, the function returns NULL.
   * Your code should correcly handle the NULL cases because otherwise a segmentation fault will crash Webots.
   *
   * This function is also often used to add joints to the simulation, e.g.
   *   dWorldID world = dBodyGetWorld(body1);
   *   pthread_mutex_lock(&mutex);
   *   dJointID joint = dJointCreateBall(world, 0);
   *   dJointAttach(joint, body1, body2);
   *   pthread_mutex_unlock(&mutex);
   *   ...
   */

  floorBody = dWebotsGetBodyFromDEF("FLOOR");
  floorGeom = dWebotsGetGeomFromDEF("FLOOR");
  
  tractorBody = dWebotsGetBodyFromDEF("TRACTOR");

  if(tractorBody == NULL){
    dWebotsConsolePrintf("!!! error : could not get boxBody.\r\n");
  }
  else {
    dWebotsConsolePrintf("tractor Body&Geom found !\r\n");
  }
  
  //if(floorBody == NULL){
  //  dWebotsConsolePrintf("!!! error : could not get floorBody.\r\n");
  //}
  //else if(floorGeom == NULL){
  //  dWebotsConsolePrintf("!!! error : could not get floorGeom.\r\n");
  //}else {
  //  dWebotsConsolePrintf("floor Body&Geom found !\r\n");
  //}

  //Init wheel id, radius, width, body, geom
  wheel_parameter_init();

  //wheels.push_back(front_right_w);
  //wheels.push_back(front_left_w);
  //wheels.push_back(rear_right_w);
  //wheels.push_back(rear_left_w);

#ifdef RECORD
  file = fopen(path, "w");
#endif
}

int tick = 0; //time step
void webots_physics_step() {
  tick++;
  /*
   * Do here what needs to be done at every time step, e.g. add forces to bodies
   *   dBodyAddForce(body1, f[0], f[1], f[2]);
   *   ...
   */
  
  double* wheel_torque_feedback;
  int size;
  wheel_torque_feedback = (double *)dWebotsReceive(&size);
  if (size != 10 * sizeof(double)) {
    dWebotsConsolePrintf("invalid receive buffer length\n");
    return;
  }else{
    wheel_motor_torque_update(wheel_torque_feedback);
    wheel_steering_angle_update(wheel_torque_feedback);
  }
  
  wheel_angular_vel_update();
  wheel_linear_vel_update();

  wheel_mechanics_calc(&front_right_w);
  wheel_mechanics_calc(&front_left_w);
  wheel_mechanics_calc(&rear_right_w);
  wheel_mechanics_calc(&rear_left_w);

  pthread_mutex_lock(&mutex);
  //const dReal* tmp = dBodyGetLinearVel(front_right_w.wheel_body);
  //dWebotsConsolePrintf("get linear vel: %f %f %f\r\n", tmp[0], tmp[1], tmp[2]);

  //FORCE
  //applying forces at the correct points and directions.
  //unit vectors in the direction of the tractor, and unit vectors in the direction of the front wheels.
  dVector3 rear_unit_vector;
  dVector3 front_unit_vector;

  //unit vectors initiation, steering angle of the two front wheels are considered the same
  wheel_force_unit_vector(rear_unit_vector);
  steer_unit_vector(front_right_w.steering_angle, front_unit_vector);

  //initiation of the motion resistance, traction force in vector form
  //dVector3 rear_traction_force;
  //initiation of the motion resistance, traction force on the front wheels in vector form.
  //dVector3 front_traction_force;
  
  // correct the direction of the resistance
  dVector3 tractor_vec;
  dBodyVectorToWorld(tractorBody, 0, 0, 1, tractor_vec); // get vector[0,0,1] in world frame
  // dWebotsConsolePrintf("tractor_vec: \t %f \t %f \t %f",tractor_vec[0],tractor_vec[1],tractor_vec[2]);

  const double* tractor_vel;
  tractor_vel = dBodyGetLinearVel(tractorBody); // get linearVel in world frame
  // dWebotsConsolePrintf("tractor_vel: \t %f \t %f \t %f",tractor_vel[0],tractor_vel[1],tractor_vel[2]);

  double lin_vel = tractor_vec[2]*tractor_vel[2]+tractor_vec[0]*tractor_vel[0]; // calculate linearVel in tractor frame
  dWebotsConsolePrintf("lin_vel: \t %f", lin_vel);


  if (lin_vel < 0){
    front_right_w.mechanics.motion_resistance*=-1; // reverse the resistance if the tractor is moving backward
    front_left_w.mechanics.motion_resistance*=-1; // reverse the resistance if the tractor is moving backward
    rear_right_w.mechanics.motion_resistance*=-1; // reverse the resistance if the tractor is moving backward
    rear_left_w.mechanics.motion_resistance*=-1; // reverse the resistance if the tractor is moving backward
  }
  
  
  //calculation of the forces. Just doubles the torque of wheels?
  // It was probably meant for when the friction coefficient was = 0
  /*
  front_traction_force[0] = front_unit_vector[0] * front_right_w.mechanics.torque / front_right_w.wheel_r;
  front_traction_force[1] = front_unit_vector[1] * front_right_w.mechanics.torque / front_right_w.wheel_r;
  front_traction_force[2] = front_unit_vector[2] * front_right_w.mechanics.torque / front_right_w.wheel_r;

  rear_traction_force[0] = rear_unit_vector[0] * rear_right_w.mechanics.torque / rear_right_w.wheel_r;
  rear_traction_force[1] = rear_unit_vector[1] * rear_right_w.mechanics.torque / rear_right_w.wheel_r;
  rear_traction_force[2] = rear_unit_vector[2] * rear_right_w.mechanics.torque / rear_right_w.wheel_r;

  //Get the global coordinates of the wheels' bottom points. the numbers are specified for this tractor.
  dVector3 fr_t, fl_t, rr_t, rl_t;
  dBodyGetRelPointPos(tractorBody, -0.7, -0.376, 1.807, fr_t);
  dBodyGetRelPointPos(tractorBody, 0.7, -0.376, 1.807, fl_t);
  dBodyGetRelPointPos(tractorBody, -0.7, -0.595, 0.075, rr_t);
  dBodyGetRelPointPos(tractorBody, 0.7, -0.595, 0.075, rl_t);
 
  //apply traction force at the contact point of the wheels.
  dBodyAddForceAtPos(front_right_w.wheel_body, front_traction_force[0], front_traction_force[1], front_traction_force[2], fr_t[0], fr_t[1], fr_t[2]);
  dBodyAddForceAtPos(front_left_w.wheel_body, front_traction_force[0], front_traction_force[1], front_traction_force[2], fl_t[0], fl_t[1], fl_t[2]);
  dBodyAddForceAtPos(rear_right_w.wheel_body, rear_traction_force[0], rear_traction_force[1], rear_traction_force[2], rr_t[0], rr_t[1], rr_t[2]);
  dBodyAddForceAtPos(rear_left_w.wheel_body, rear_traction_force[0], rear_traction_force[1], rear_traction_force[2], rl_t[0], rl_t[1], rl_t[2]);

  dWebotsConsolePrintf("front traction force: \t %f", front_traction_force[2]);
  dWebotsConsolePrintf("rear traction force: \t %f", rear_traction_force[2]);
*/

  //apply motion resistance at the center of the mass of the tractor. 
  // applying it to Wheels doesnt work because I think the coordinate frames spin
 double total_motion_resistance = -front_right_w.mechanics.motion_resistance-front_left_w.mechanics.motion_resistance-rear_right_w.mechanics.motion_resistance-rear_left_w.mechanics.motion_resistance;
 dBodyAddRelForce(tractorBody, 0, 0, total_motion_resistance);
 dWebotsConsolePrintf("motion resistance: \t %f", total_motion_resistance);

/*
  dWebotsConsolePrintf("motion resistance: \t %f", -front_right_w.mechanics.motion_resistance);
  dWebotsConsolePrintf("motion resistance: \t %f", -front_left_w.mechanics.motion_resistance);
  dWebotsConsolePrintf("motion resistance: \t %f", -rear_right_w.mechanics.motion_resistance);
  dWebotsConsolePrintf("motion resistance: \t %f", -rear_left_w.mechanics.motion_resistance);
*/

#ifdef ADD_LATERAL_FORCE
  // Problem: I followed this all the way thorugh and it mostly works correctly, 
  // Why does it keep coming back to where it started??
  // but does not work with motion resistance on crcle testing
  if (front_right_w.steering_angle != 0){
    // This adds lateral force and torque required for more accurate turning motion (understeer)
    dWebotsConsolePrintf("applied lateral force: \t %f", TOTAL_MASS*front_right_w.mechanics.lateral_accel);
    dBodyAddRelForce(tractorBody, TOTAL_MASS*front_right_w.mechanics.lateral_accel, 0, 0);
    
    dWebotsConsolePrintf("applied torque: \t %f", MOMENT_OF_INERTIA*front_right_w.mechanics.angular_accel_yaxis);
    //dWebotsConsolePrintf("angular torque: \t %f", MOMENT_OF_INERTIA*front_right_w.mechanics.angular_accel_yaxis);
    dBodyAddRelTorque(tractorBody, 0, MOMENT_OF_INERTIA*front_right_w.mechanics.angular_accel_yaxis, 0);
  }
#endif
  
  const dReal* tractor_force = dBodyGetForce(tractorBody);
  const dReal* tractor_torque = dBodyGetTorque(tractorBody);
  // This is useful if you apply forces on the wheel, but we dont
  const dReal* FR_wheel_force = dBodyGetForce(front_right_w.wheel_body);
  const dReal* FL_wheel_force = dBodyGetForce(front_left_w.wheel_body);
  const dReal* RR_wheel_force = dBodyGetForce(rear_right_w.wheel_body);
  const dReal* RL_wheel_force = dBodyGetForce(rear_left_w.wheel_body);
  
  dWebotsConsolePrintf("x tractor force: \t %f", tractor_force[0]);
  dWebotsConsolePrintf("z tractor force: \t %f", tractor_force[2]);
  //dWebotsConsolePrintf("lateral force: \t %f", TOTAL_MASS*4.0*4.0/10.9 + tractor_force[0]);
  dWebotsConsolePrintf("tractor torque: \t %f", tractor_torque[1]);
  //dWebotsConsolePrintf("FR Wheel force: \t %f", front_right_w.mechanics.torque / front_right_w.wheel_r - FR_wheel_force[2]);
  //dWebotsConsolePrintf("FL Wheel force: \t %f", front_left_w.mechanics.torque / front_left_w.wheel_r - FL_wheel_force[2]);
  //dWebotsConsolePrintf("RR Wheel force: \t %f", rear_right_w.mechanics.torque / rear_right_w.wheel_r - RR_wheel_force[2]);
  //dWebotsConsolePrintf("RL Wheel force: \t %f", rear_left_w.mechanics.torque / rear_left_w.wheel_r - RL_wheel_force[2]);
  
  double sum_of_z_forces = sqrt(tractor_force[0]*tractor_force[0] + tractor_force[2]*tractor_force[2]) 
                                            + front_right_w.mechanics.torque / front_right_w.wheel_r - FR_wheel_force[2] 
                                            + front_left_w.mechanics.torque / front_left_w.wheel_r - FL_wheel_force[2]
                                            + rear_right_w.mechanics.torque / rear_right_w.wheel_r - RR_wheel_force[2]
                                            + rear_left_w.mechanics.torque / rear_left_w.wheel_r - RL_wheel_force[2];
  dWebotsConsolePrintf("Total forward force: \t %f", sum_of_z_forces);
  
  pthread_mutex_unlock(&mutex);

  print_output();

  //dBodyAddRelForce(frontRightWheelBody, f[0], f[1], f[2]);
  //dBodyAddRelForce(frontLeftWheelBody, f[0], f[1], f[2]);

#ifdef RECORD
  fprintf(file, "%d\n", tick);
  if(tick>=100){
    fclose(file);
  }
#endif

}

void webots_physics_step_end(){

  wheel_load_update();
  //dWebotsConsolePrintf("Load:%f %f %f %f\n", front_right_w.joint_feed_back.f1[1], front_left_w.joint_feed_back.f1[1],
  //                                           rear_right_w.joint_feed_back.f1[1], rear_left_w.joint_feed_back.f1[1]);
  //dWebotsConsolePrintf("F1:%f %f %f  F2:%f %f %f\n", front_right_w.joint_feed_back.f1[0], front_right_w.joint_feed_back.f1[1], front_right_w.joint_feed_back.f1[2],
  //                                                   front_right_w.joint_feed_back.f2[0], front_right_w.joint_feed_back.f2[1], front_right_w.joint_feed_back.f2[2]);
  
  //pthread_mutex_lock(&mutex);
  //pthread_mutex_unlock(&mutex);
}

int webots_physics_collide(dGeomID g1, dGeomID g2) {
  /*
   * This function needs to be implemented if you want to overide Webots collision detection.
   * It must return 1 if the collision was handled and 0 otherwise.
   * Note that contact joints should be added to the contact_joint_group which can change over the time, e.g.
   *   n = dCollide(g1, g2, MAX_CONTACTS, &contact[0].geom, sizeof(dContact));
   *   dJointGroupID contact_joint_group = dWebotsGetContactJointGroup();
   *   dWorldID world = dBodyGetWorld(body1);
   *   ...
   *   pthread_mutex_lock(&mutex);
   *   dJointCreateContact(world, contact_joint_group, &contact[i])
   *   dJointAttach(contact_joint, body1, body2);
   *   pthread_mutex_unlock(&mutex);
   *   ...
   */
  dContact contact;
  dBodyID body;
  dWorldID world;
  dJointGroupID contact_joint_group;
  int contact_num;


  if (dAreGeomsSame(g1, front_right_w.wheel_geom) || dAreGeomsSame(g2, front_right_w.wheel_geom)){
    if((dAreGeomsSame(g1,front_right_w.wheel_geom) && dAreGeomsSame(g2,floorGeom)) || (dAreGeomsSame(g1,floorGeom) && dAreGeomsSame(g2,front_right_w.wheel_geom))){

      contact_num = dCollide(g1, g2, 1, &contact.geom, sizeof(dContact));
      if(contact_num == 0)
        return 1;

      body = dGeomGetBody(g1);
      world = dBodyGetWorld(body);  
      contact_joint_group = dWebotsGetContactJointGroup();

      contact.surface.mode =  dContactBounce | dContactSoftCFM | dContactSoftERP;
      contact.surface.soft_cfm = 0.0;
      contact.surface.soft_erp = 0.2;
      contact.surface.bounce = 0.0;
      //contact.surface.mu = dInfinity;
      contact.surface.mu = MU;
      //contact.surface.mu = wheel_mu_calc(&front_right_w);

      pthread_mutex_lock(&mutex);
      front_right_w.joint = dJointCreateContact(world, contact_joint_group, &contact);
      dJointAttach(front_right_w.joint, front_right_w.wheel_body, NULL);
      dJointSetFeedback(front_right_w.joint, &(front_right_w.joint_feed_back));
      pthread_mutex_unlock(&mutex);
      return 1;
    }
  }

  if (dAreGeomsSame(g1, front_left_w.wheel_geom) || dAreGeomsSame(g2, front_left_w.wheel_geom)){
    if((dAreGeomsSame(g1,front_left_w.wheel_geom) && dAreGeomsSame(g2,floorGeom)) || (dAreGeomsSame(g1,floorGeom) && dAreGeomsSame(g2,front_left_w.wheel_geom))){

      contact_num = dCollide(g1, g2, 1, &contact.geom, sizeof(dContact));
      if(contact_num == 0)
        return 1;

      body = dGeomGetBody(g1);
      world = dBodyGetWorld(body);  
      contact_joint_group = dWebotsGetContactJointGroup();

      contact.surface.mode =  dContactBounce | dContactSoftCFM | dContactSoftERP;
      contact.surface.soft_cfm = 0.0;
      contact.surface.soft_erp = 0.2;
      contact.surface.bounce = 0.0;
      //contact.surface.mu = dInfinity;
      contact.surface.mu = MU;
      //contact.surface.mu = wheel_mu_calc(&front_left_w);

      pthread_mutex_lock(&mutex);
      front_left_w.joint = dJointCreateContact(world, contact_joint_group, &contact);
      dJointAttach(front_left_w.joint, front_left_w.wheel_body, NULL);
      dJointSetFeedback(front_left_w.joint, &(front_left_w.joint_feed_back));
      pthread_mutex_unlock(&mutex);
      return 1;
    }
  }

  if (dAreGeomsSame(g1, rear_right_w.wheel_geom) || dAreGeomsSame(g2, rear_right_w.wheel_geom)){
    if((dAreGeomsSame(g1,rear_right_w.wheel_geom) && dAreGeomsSame(g2,floorGeom)) || (dAreGeomsSame(g1,floorGeom) && dAreGeomsSame(g2,rear_right_w.wheel_geom))){

      contact_num = dCollide(g1, g2, 1, &contact.geom, sizeof(dContact));
      if(contact_num == 0)
        return 1;

      body = dGeomGetBody(g1);
      world = dBodyGetWorld(body);  
      contact_joint_group = dWebotsGetContactJointGroup();

      contact.surface.mode =  dContactBounce | dContactSoftCFM | dContactSoftERP;
      contact.surface.soft_cfm = 0.0;
      contact.surface.soft_erp = 0.2;
      contact.surface.bounce = 0.0;
      //contact.surface.mu = dInfinity;
      contact.surface.mu = MU;
      //contact.surface.mu = wheel_mu_calc(&rear_right_w);

      pthread_mutex_lock(&mutex);
      rear_right_w.joint = dJointCreateContact(world, contact_joint_group, &contact);
      dJointAttach(rear_right_w.joint, rear_right_w.wheel_body, NULL);
      dJointSetFeedback(rear_right_w.joint, &(rear_right_w.joint_feed_back));
      pthread_mutex_unlock(&mutex);
      return 1;
    }
  }

  if (dAreGeomsSame(g1, rear_left_w.wheel_geom) || dAreGeomsSame(g2, rear_left_w.wheel_geom)){
    if((dAreGeomsSame(g1,rear_left_w.wheel_geom) && dAreGeomsSame(g2,floorGeom)) || (dAreGeomsSame(g1,floorGeom) && dAreGeomsSame(g2,rear_left_w.wheel_geom))){

      contact_num = dCollide(g1, g2, 1, &contact.geom, sizeof(dContact));
      if(contact_num == 0)
        return 1;

      body = dGeomGetBody(g1);
      world = dBodyGetWorld(body);  
      contact_joint_group = dWebotsGetContactJointGroup();

      contact.surface.mode =  dContactBounce | dContactSoftCFM | dContactSoftERP;
      contact.surface.soft_cfm = 0.0;
      contact.surface.soft_erp = 0.2;
      contact.surface.bounce = 0.0;
      //contact.surface.mu = dInfinity;
      contact.surface.mu = MU;
      //contact.surface.mu = wheel_mu_calc(&rear_left_w);
      
      //dWebotsConsolePrintf("mr: \t %f", rear_left_w.mechanics.motion_resistance );
      //dWebotsConsolePrintf("theta: \t %f", rear_left_w.mechanics.theta_1);
      //dWebotsConsolePrintf("load: \t %f", rear_left_w.load);
      //dWebotsConsolePrintf("mu: \t %f", contact.surface.mu);


      pthread_mutex_lock(&mutex);
      rear_left_w.joint = dJointCreateContact(world, contact_joint_group, &contact);
      dJointAttach(rear_left_w.joint, rear_left_w.wheel_body, NULL);
      dJointSetFeedback(rear_left_w.joint, &(rear_left_w.joint_feed_back));
      pthread_mutex_unlock(&mutex);
      return 1;
    }
  }

  return 0;
}

void webots_physics_cleanup() {
  /*
   * Here you need to free any memory you allocated in above, close files, etc.
   * You do not need to free any ODE object, they will be freed by Webots.
   */
  pthread_mutex_destroy(&mutex);

#ifdef RECORD
  fclose(file);
#endif
}