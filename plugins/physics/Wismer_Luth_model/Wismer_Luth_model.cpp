/*
 * File: Wismer_luth_model.cpp
 * Date: 2019.8.1
 * Description: physics plugin for Webots simulation 
 * Author: Yonglin Jing 
 * E-mail: 11712605@mail.sustech.edu.cn 
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

#include <plugins/physics.h>
#include <math.h>
#include "Wismer_luth_model.h"

//Using Wismer Luth model

static double last_time = 0.0;
static double current_time = 0.0;

static wheel_t front_right_w, front_left_w, rear_right_w, rear_left_w;
static double cone_index = 0.0;

static dBodyID floorBody = NULL;
static dGeomID floorGeom = NULL;

static dBodyID tractorBody = NULL;

static pthread_mutex_t mutex;  // needed to run with multi-threaded version of ODE

//recommand using constant cone index instead of the calculation function
//Note that cone index is not decided by radius and theta1
double cone_index_calc(double r,   //wheel radius
                       double t,   //theta1, dependent on load, or use pre-calculated THETA_CONSTANT
                       double kc,  //cohesion modulus of deformation
                       double se,  //sinkage exponent
                       double kp){ //frictional modulus of deformation
  double zetam = 0.0;
  zetam = r*(1-cos(t));

  return 1.625*(kc/(se+1)*(pow((zetam+1.5),(se+1))-pow(zetam, (se+1)))+0.517*kp*(pow((zetam+1.5),\
         (se+2))/(se+1)/(se+2)+pow(zetam, (se+2))/(se+2)-(zetam+1.5)*pow(zetam, (se+1))/(se+1)));
}

inline double slip_ratio_calc(double r,  //wheel radius 
                              double w,  //wheel angular velocity
                              double vx){  //wheel horizontal velocity
  return 1-(vx/(r*w));
}

double theory_slip_ratio_calc(double t,     //torque
                              double r,     //wheel radius
                              double I,     //inertia of wheel
                              double tf,    //towed force
                              double fc,    //friction coefficient
                              double vx,    //horizontal velocity
                              double L,     //wheel load
                              double& omega,  //theoretical angular velocity
                              double& alpha,  //theoretical angular acceleration
                              wheel_id_e id){   

  //static double omega = 0.0, alpha = 0.0;
  double  m = 0.0, threshold = 0.0, s = 0.0;

  m = L/GRAVITY_CONSTANT;
  threshold = fc*L*(I/(m*r) + r);

  if(t < threshold && (omega <= vx/r)){
    omega = vx/r;
    alpha = 0.0;
  }else{
    alpha = (t - fc*L*r)/I;
    omega += alpha * TIME_STEP;       
    if(id == FRONT_RIGHT_ID)
      dWebotsConsolePrintf("!-------");
  }
    
  s = 1 - (vx/(omega*r));

  if(id == FRONT_RIGHT_ID)
    dWebotsConsolePrintf("sr:%f tsr:%f vx:%f w:%f omega:%f alpha:%f torque:%f thres:%f Load:%f tmp:%f error_persent:%f\n", 
        front_right_w.mechanics.slip_ratio, s, vx, front_right_w.angular_vel,
        omega, alpha, t, threshold, L, fc*L*r, (omega - front_right_w.angular_vel)/front_right_w.angular_vel);
  
  return s;
  
  //return (-10/(3*wn))*log(1-(4*t)/(3*r*L));
}

double friction_coefficient_calc(double t,    //torque
                                double r,     //wheel radius
                                double I,     //inertia of wheel
                                double vx,    //horizontal velocity
                                double L,     //wheel load
                                double wn,    //wheel numeric
                                double& ti,     //torque integral
                                double& u_mid,  //midium value used for friction coefficient calculation 
                                double& last_u_mid,  //midium value used for friction coefficient calculation
                                wheel_id_e id ){     //wheel id for debug

  double u = 0.0;

  ti += t * TIME_STEP;
  
  u_mid = ((-vx*I)/(1+(10*log(1-(4*t)/(3*r*L)))/(3*wn)) + r*ti)/(L*r*r);

  u = ((u_mid-last_u_mid)/TIME_STEP);
  //u = ((u_mid)/(current_time/1000));

  last_u_mid = u_mid;

  return u;
}

inline double wheel_numeric_calc(double b,  //wheel width 
                                 double r,  //wheel radius
                                 double L){ //wheel load
  return cone_index*b*2*r/L;
}

inline double towed_force_calc(double L,   //wheel load
                               double wn){ //wheel numeric
  return L*(1.2/wn+0.04);
}

inline double torque_calc(double r,   //wheel radius
                          double L,   //wheel load
                          double wn,  //wheel numeric
                          double sr){ //slip ratio
  
  return r*L*0.75*(1-exp(-0.3*wn*sr));
}

inline double drawbar_pull_calc(double ff,  //friction force
                                double tf){ //towed force
  return (ff-tf > 0) ? (ff-tf) : 0;
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

//A cheating way to handle the effect of slip ratio out of Webots
//Which means there will no slip effect in Webots
inline double actual_distance_calc(double e,        //wheel encoder data
                                   double r,        //wheel radius
                                   double sr){      //slip ratio
  return e*r*sr;
}

void wheel_parameter_init() {
  front_right_w.id      = FRONT_RIGHT_ID;
  front_right_w.wheel_r = FRONT_WHEEL_RADIUS;
  front_right_w.wheel_w = FRONT_WHEEL_WIDTH;
  front_right_w.load    = FRONT_WHEEL_LOAD_INIT;

  front_left_w.id       = FRONT_LEFT_ID;
  front_left_w.wheel_r  = FRONT_WHEEL_RADIUS;
  front_left_w.wheel_w  = FRONT_WHEEL_WIDTH;
  front_left_w.load     = FRONT_WHEEL_LOAD_INIT;
  
  rear_right_w.id       = REAR_RIGHT_ID;
  rear_right_w.wheel_r  = REAR_WHEEL_RADIUS;
  rear_right_w.wheel_w  = REAR_WHEEL_WIDTH;
  rear_right_w.load     = REAR_WHEEL_LOAD_INIT;
  
  rear_left_w.id        = REAR_LEFT_ID;
  rear_left_w.wheel_r   = REAR_WHEEL_RADIUS;
  rear_left_w.wheel_w   = REAR_WHEEL_WIDTH;
  rear_left_w.load      = REAR_WHEEL_LOAD_INIT;

  


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

  dBodyGetMass(front_right_w.wheel_body, &front_right_w.mass);
  dBodyGetMass(front_left_w.wheel_body, &front_left_w.mass);
  dBodyGetMass(rear_right_w.wheel_body, &rear_right_w.mass);
  dBodyGetMass(rear_left_w.wheel_body, &rear_left_w.mass);

  front_right_w.mechanics.inertia = front_right_w.mass.I[5];
  front_left_w.mechanics.inertia = front_left_w.mass.I[5];
  rear_right_w.mechanics.inertia = rear_right_w.mass.I[5];
  rear_left_w.mechanics.inertia = rear_left_w.mass.I[5];
  
}

void wheel_mechanics_calc(wheel_t* wheel){
  //Be careful when you change the caculation order

  wheel->mechanics.slip_ratio = slip_ratio_calc(wheel->wheel_r, 
                                                wheel->angular_vel, 
                                                wheel->horiz_vel);
  
  wheel->mechanics.driving_force = wheel_driving_force_calc(wheel->mechanics.torque, 
                                                            wheel->wheel_r);

  wheel->mechanics.wheel_numeric = wheel_numeric_calc(wheel->wheel_w, 
                                                      wheel->wheel_r, 
                                                      wheel->load);

  wheel->mechanics.towed_force = towed_force_calc(wheel->load, 
                                                  wheel->mechanics.wheel_numeric);

  //Terminate towed force in case that towed force will cause tractor move backwards
  if(wheel->horiz_vel < 0.1){
    wheel->mechanics.towed_force = 0;
    dWebotsConsolePrintf("towed forced terminated\n");
  }

  wheel->mechanics.friction_coef = friction_coefficient_calc(wheel->mechanics.torque, 
                                                             wheel->wheel_r,
                                                             wheel->mechanics.inertia, 
                                                             wheel->horiz_vel, 
                                                             wheel->load, 
                                                             wheel->mechanics.wheel_numeric, 
                                                             wheel->calc_mid_value.torque_integral, 
                                                             wheel->calc_mid_value.friction_coef_mid, 
                                                             wheel->calc_mid_value.last_friction_coef_mid, 
                                                             wheel->id);

  wheel->mechanics.drawbar_pull = drawbar_pull_calc(wheel->friction_force,
                                                    wheel->mechanics.towed_force);

  wheel->mechanics.trct_effc = tractive_efficiency_calc(wheel->mechanics.drawbar_pull, 
                                                        wheel->horiz_vel, 
                                                        wheel->mechanics.torque, 
                                                        wheel->angular_vel);

  wheel->mechanics.actual_distance = actual_distance_calc(wheel->encode, 
                                                          wheel->wheel_r, 
                                                          wheel->mechanics.thoery_slip_ratio);
}

double wheel_mu_calc(wheel_t* wheel){
  return (wheel->mechanics.towed_force / wheel->load);
}

void wheel_load_update(){
  front_right_w.load = front_right_w.joint_feed_back.f1[Y];
  front_left_w.load  = front_left_w.joint_feed_back.f1[Y];
  rear_right_w.load  = rear_right_w.joint_feed_back.f1[Y];
  rear_left_w.load   = rear_left_w.joint_feed_back.f1[Y];
}

void wheel_friction_force_update(){
  front_right_w.friction_force = front_right_w.joint_feed_back.f1[Z];
  front_left_w.friction_force  = front_left_w.joint_feed_back.f1[Z];
  rear_right_w.friction_force  = rear_right_w.joint_feed_back.f1[Z];
  rear_left_w.friction_force   = rear_left_w.joint_feed_back.f1[Z];
}

void wheel_motor_update(double* recv_buf){

  front_left_w.mechanics.torque  = recv_buf[0];
  front_right_w.mechanics.torque = recv_buf[1];
  rear_left_w.mechanics.torque   = recv_buf[2];
  rear_right_w.mechanics.torque  = recv_buf[3];

  front_left_w.encode = recv_buf[4];
  front_right_w.encode = recv_buf[5];
  rear_left_w.encode = recv_buf[6];
  rear_right_w.encode = recv_buf[7];

  //friction_frac = feedback[0];

  //dWebotsConsolePrintf("encoder:%f, %f, %f, %f\n", feedback[4], feedback[5], feedback[6], feedback[7]);
}

void wheel_angular_vel_update(){
  const dReal* tmp;

  //only consider one axis 
  //coordinate transform required if consider steering
  tmp = dBodyGetAngularVel(front_right_w.wheel_body);
  front_right_w.angular_vel = tmp[X];
  tmp = dBodyGetAngularVel(front_left_w.wheel_body);
  front_left_w.angular_vel = tmp[X];
  tmp = dBodyGetAngularVel(rear_right_w.wheel_body);
  rear_right_w.angular_vel = tmp[X];
  tmp = dBodyGetAngularVel(rear_left_w.wheel_body);
  rear_left_w.angular_vel = tmp[X];
}

void wheel_linear_vel_update(){
  const dReal* tmp;

  //only consider one axis 
  //coordinate transform required if consider steering
  tmp = dBodyGetLinearVel(front_right_w.wheel_body);
  front_right_w.horiz_vel = tmp[Z];
  tmp = dBodyGetLinearVel(front_left_w.wheel_body);
  front_left_w.horiz_vel = tmp[Z];
  tmp = dBodyGetLinearVel(rear_right_w.wheel_body);
  rear_right_w.horiz_vel = tmp[Z];
  tmp = dBodyGetLinearVel(rear_left_w.wheel_body);
  rear_left_w.horiz_vel = tmp[Z];

}

void wheel_steering_angle_update(double* recv_buf){
  front_right_w.steering_angle = recv_buf[8];
  front_left_w.steering_angle = recv_buf[9];
}

void print_output(){
  dWebotsConsolePrintf("Front Right-slip ratio:%f, Rear Right-slip ratio:%f, thry slip ratio:%f, wheel numeric:%f, tractive eficiency:%f, drawbar pull:%f, towed force:%f, torque:%f, w:%f, vx:%f, \n", //frac:%f\n",
                      front_right_w.mechanics.slip_ratio, rear_right_w.mechanics.slip_ratio, front_right_w.mechanics.thoery_slip_ratio, front_right_w.mechanics.wheel_numeric, front_right_w.mechanics.trct_effc,
                       front_right_w.mechanics.drawbar_pull, front_right_w.mechanics.towed_force, front_right_w.mechanics.torque, 
                        front_right_w.angular_vel, front_right_w.horiz_vel); //, front_right_w.mechanics.friction_frac);
  
  //dWebotsConsolePrintf("FrontRight: %.3f FrontLeft:%.3f RearRight:%.3f RearLeft:%.3f FrontSlipRatio:%.3f RearSlipRatio:%.3f\n", 
  //                      front_right_w.mechanics.friction_frac, front_left_w.mechanics.friction_frac,
  //                      rear_right_w.mechanics.friction_frac, rear_left_w.mechanics.friction_frac,
  //                      front_right_w.mechanics.slip_ratio, rear_right_w.mechanics.slip_ratio);

  //dWebotsConsolePrintf("FrontRight: ti:%f u_mid:%f error:%f fric:%f vx:%f slip:%f\n", front_right_w.calc_mid_value.torque_integral, 
  //                                                                front_right_w.calc_mid_value.friction_coef_mid,
  //                                                                front_right_w.calc_mid_value.friction_coef_mid - front_right_w.calc_mid_value.last_friction_coef_mid,
  //                                                                front_right_w.mechanics.friction_frac,
  //                                                                front_right_w.horiz_vel,
  //                                                                front_right_w.mechanics.slip_ratio);

  //dWebotsConsolePrintf("time step:%f\n", current_time - last_time);

  //dWebotsConsolePrintf("Load:%f %f %f vel:%f torque:%f %f\n", front_right_w.joint_feed_back.f1[0] + front_left_w.joint_feed_back.f1[0] + rear_right_w.joint_feed_back.f1[0] + rear_left_w.joint_feed_back.f1[0], 
  //                                                         front_right_w.joint_feed_back.f1[1] + front_left_w.joint_feed_back.f1[1] + rear_right_w.joint_feed_back.f1[1] + rear_left_w.joint_feed_back.f1[1], 
  //                                                         front_right_w.joint_feed_back.f1[2] + front_left_w.joint_feed_back.f1[2] + rear_right_w.joint_feed_back.f1[2] + rear_left_w.joint_feed_back.f1[2], 
  //                                                         front_right_w.horiz_vel, front_right_w.mechanics.torque * 4, front_right_w.mechanics.torque + front_left_w.mechanics.torque + rear_right_w.mechanics.torque + rear_left_w.mechanics.torque);

  //dWebotsConsolePrintf("BOX: vel:%f, load:%f, fractor:%f, tow:%f\n", box_vel, box_load, friction_coef, box_tow);
  //dWebotsConsolePrintf("[%f %f %f]\n[%f %f %f]\n[%f %f %f]\n -----------------\n", 
  //                      rear_right_w.mass.I[0], rear_right_w.mass.I[1], rear_right_w.mass.I[2],
  //                      rear_right_w.mass.I[4], rear_right_w.mass.I[5], rear_right_w.mass.I[6],
  //                      rear_right_w.mass.I[8], rear_right_w.mass.I[9], rear_right_w.mass.I[10]);
  //dWebotsConsolePrintf("FR-vx:%f torque:%f friction:%f towedForce:%f drawbarPull:%f uL:%f slipRatio:%f fricForce:%f Q/r:%f Load:%f\nFL-vx:%f torque:%f friction:%f towedForce:%f drawbarPull:%f uL:%f slipRatio:%f fricForce:%f Q/r:%f Load:%f\nRR-vx:%f torque:%f friction:%f towedForce:%f drawbarPull:%f uL:%f slipRatio:%f fricForce:%f Q/r:%f Load:%f\nRL-vx:%f torque:%f friction:%f towedForce:%f drawbarPull:%f uL:%f slipRatio:%f fricForce:%f Q/r:%f Load:%f",
  //                      front_right_w.horiz_vel, front_right_w.mechanics.torque, front_right_w.mechanics.friction_coef, front_right_w.mechanics.towed_force, front_right_w.mechanics.drawbar_pull, front_right_w.mechanics.friction_coef * front_right_w.load, front_right_w.mechanics.slip_ratio, front_right_w.friction_force, front_right_w.mechanics.torque/front_right_w.wheel_r, front_right_w.load,
  //                      front_left_w.horiz_vel, front_left_w.mechanics.torque, front_left_w.mechanics.friction_coef, front_left_w.mechanics.towed_force, front_left_w.mechanics.drawbar_pull, front_left_w.mechanics.friction_coef * front_left_w.load, front_left_w.mechanics.slip_ratio, front_left_w.friction_force, front_left_w.mechanics.torque/front_left_w.wheel_r, front_left_w.load,
  //                      rear_right_w.horiz_vel, rear_right_w.mechanics.torque, rear_right_w.mechanics.friction_coef, rear_right_w.mechanics.towed_force, rear_right_w.mechanics.drawbar_pull, rear_right_w.mechanics.friction_coef * rear_right_w.load, rear_right_w.mechanics.slip_ratio, rear_right_w.friction_force, rear_right_w.mechanics.torque/rear_right_w.wheel_r, rear_right_w.load,
  //                     rear_left_w.horiz_vel, rear_left_w.mechanics.torque, rear_left_w.mechanics.friction_coef, rear_left_w.mechanics.towed_force, rear_left_w.mechanics.drawbar_pull, rear_left_w.mechanics.friction_coef * rear_left_w.load, rear_left_w.mechanics.slip_ratio, rear_left_w.friction_force, rear_left_w.mechanics.torque/rear_left_w.wheel_r, rear_left_w.load);
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
  
  if(floorBody == NULL){
    dWebotsConsolePrintf("!!! error : could not get floorBody.\r\n");
  }
  else if(floorGeom == NULL){
    dWebotsConsolePrintf("!!! error : could not get floorGeom.\r\n");
  }else {
    dWebotsConsolePrintf("floor Body&Geom found !\r\n");
  }

  tractorBody = dWebotsGetBodyFromDEF("TRACTOR");

  if(tractorBody == NULL){
    dWebotsConsolePrintf("!!! error : could not get boxBody.\r\n");
  }
  else {
    dWebotsConsolePrintf("tractor Body&Geom found !\r\n");
  }

  //If cone index is known then there is no need to use the calculate function
  cone_index = CONE_INDEX_CONSTANT;

  //Init wheel id, radius, width, body, geom
  wheel_parameter_init();

}

void webots_physics_step() {
  /*
   * Do here what needs to be done at every time step, e.g. add forces to bodies
   *   dBodyAddForce(body1, f[0], f[1], f[2]);
   *   ...
   */
  
  current_time = dWebotsGetTime();

  double* recv_buf;
  int size;
  recv_buf = (double *)dWebotsReceive(&size);
  if (size != EMMITER_RECV_BUF_LENGTH * sizeof(EMMITER_BUF_DATA_TYPE)) {
    dWebotsConsolePrintf("invalid receive buffer length\n");
    return;
  }else{
    wheel_motor_update(recv_buf);
    wheel_steering_angle_update(recv_buf);
  }

  wheel_angular_vel_update();
  wheel_linear_vel_update();

  wheel_mechanics_calc(&front_right_w);
  wheel_mechanics_calc(&front_left_w);
  wheel_mechanics_calc(&rear_right_w);
  wheel_mechanics_calc(&rear_left_w);

  pthread_mutex_lock(&mutex);

  dBodyAddForce(front_right_w.wheel_body, 0, 0, -front_right_w.mechanics.towed_force);
  dBodyAddForce(front_left_w.wheel_body, 0, 0, -front_left_w.mechanics.towed_force);
  dBodyAddForce(rear_right_w.wheel_body, 0, 0, -rear_right_w.mechanics.towed_force);
  dBodyAddForce(rear_left_w.wheel_body, 0, 0, -rear_left_w.mechanics.towed_force);

  pthread_mutex_unlock(&mutex);

  print_output();

  last_time = current_time;

}

void webots_physics_step_end(){

  //wheel_load_update();  //uncomment this to get dynamic load
  wheel_friction_force_update();

  //total_load = front_right_w.joint_feed_back.f1[1] + front_left_w.joint_feed_back.f1[1] + rear_right_w.joint_feed_back.f1[1] + rear_left_w.joint_feed_back.f1[1];
  //error = total_load-TOTAL_GRAVITY;
  //error_persent = (error/TOTAL_GRAVITY)*100;
  //dWebotsConsolePrintf("Load:%f %f %f %f, Total:%f, Error:%f, Error persentage:%f%\n", front_right_w.joint_feed_back.f1[1], front_left_w.joint_feed_back.f1[1],
  //                                           rear_right_w.joint_feed_back.f1[1], rear_left_w.joint_feed_back.f1[1], total_load, error, error_persent);
  //dWebotsConsolePrintf("F1:%f %f %f  F2:%f %f %f\n", front_right_w.joint_feed_back.f1[0], front_right_w.joint_feed_back.f1[1], front_right_w.joint_feed_back.f1[2],
  //                                                  front_right_w.joint_feed_back.f2[0], front_right_w.joint_feed_back.f2[1], front_right_w.joint_feed_back.f2[2]);
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

      contact.surface.mode =  dContactBounce | dContactSoftCFM | dContactSoftERP | dContactApprox1;
      contact.surface.soft_cfm = 0.0;
      contact.surface.soft_erp = 0.2;
      contact.surface.bounce = 0.0;
      contact.surface.bounce_vel = 0.0;
      contact.surface.mu = front_right_w.mechanics.friction_coef;
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

      contact.surface.mode =  dContactBounce | dContactSoftCFM | dContactSoftERP | dContactApprox1;
      contact.surface.soft_cfm = 0.0;
      contact.surface.soft_erp = 0.2;
      contact.surface.bounce = 0.0;
      contact.surface.bounce_vel = 0.0;
      contact.surface.mu = front_left_w.mechanics.friction_coef;
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

      contact.surface.mode =  dContactBounce | dContactSoftCFM | dContactSoftERP | dContactApprox1;
      contact.surface.soft_cfm = 0.0;
      contact.surface.soft_erp = 0.2;
      contact.surface.bounce = 0.0;
      contact.surface.bounce_vel = 0.0;
      contact.surface.mu = rear_right_w.mechanics.friction_coef;
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

      contact.surface.mode =  dContactBounce | dContactSoftCFM | dContactSoftERP | dContactApprox1;
      contact.surface.soft_cfm = 0.0;
      contact.surface.soft_erp = 0.2;
      contact.surface.bounce = 0.0;
      contact.surface.bounce_vel = 0.0;
      contact.surface.mu = rear_left_w.mechanics.friction_coef;
      //contact.surface.mu = wheel_mu_calc(&rear_left_w);

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
}
