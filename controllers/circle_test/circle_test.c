/*
 * File:          circle_test.c
 * Date:
 * Description:
 * Author:
 * Modifications:
 */

/*
 * FIND THE RADIUS BY TESTIING IT WITH THE PARKING LOT 
 * RADIUS = MAX_LEFT / 2 
 * RADIUS = 12 m right now
 * MOVE TRACTOR 10 m BACK AND RADIUS LENGTH TO THE RIGHT 
 */
#include <stdio.h>
#include <string.h>
#include <webots/camera.h>
#include <webots/gps.h>
#include <webots/keyboard.h>
#include <webots/led.h>
#include <webots/lidar.h>
#include <webots/motor.h>
#include <webots/robot.h>
#include <webots/emitter.h>
#include <webots/touch_sensor.h>
#include <webots/position_sensor.h>
#include <math.h>

// to be used as array indices
enum { X, Y, Z };

// This needs to be changed if the .wbt model changes
#define FRONT_WHEEL_RADIUS 0.38
#define REAR_WHEEL_RADIUS 0.6

//comment SPEED_CONTROL to use torque control
#define SPEED_CONTROL

//uncomment to recrod
#define RECORD

double emmiter_send_buf[10];
double friction_frac = 400.0;

#ifdef RECORD
//create file to export data
FILE *file_circle = NULL;
char* file_path = "file_circle.txt"; // file path
#endif


// devices
WbDeviceTag left_steer, right_steer;
WbDeviceTag left_front_wheel, right_front_wheel;
WbDeviceTag left_rear_wheel, right_rear_wheel;
WbDeviceTag left_front_encoder, right_front_encoder;
WbDeviceTag left_rear_encoder, right_rear_encoder;

// lights
WbDeviceTag left_flasher, right_flasher, tail_lights;
WbDeviceTag work_head_lights, road_head_lights;

// camera
WbDeviceTag camera;
int camera_width = -1;
int camera_height = -1;
double camera_fov = -1.0;

// SICK laser
WbDeviceTag sick;
int sick_width = -1;
double sick_range = -1.0;
double sick_fov = -1.0;

// GPS
WbDeviceTag gps;
double gps_coords[3];
double gps_speed = 0.0;

//Touch sensor
//WbDeviceTag touch_sensor;
//WbTouchSensorType touch_type;
//const double* wheel_force_3d;

//Emitter
WbDeviceTag emitter;

// misc variables
int time_step = -1;
double speed = 0.0;
double steering_angle = 0.0;
double manual_steering = 0.0;
double torque = 0.0; //Added by NERanger 20190731

//Debug parameters 
//Added by NERanger 20190731
double  left_front_wheel_torque_feedback    = 0, 
        right_front_wheel_torque_feedback   = 0,
        left_rear_wheel_torque_feedback     = 0,
        right_rear_wheel_torque_feedback    = 0,

        left_front_wheel_velocity_feedback  = 0,
        right_front_wheel_velocity_feedback = 0,
        left_rear_wheel_velocity_feedback   = 0,
        right_rear_wheel_velocity_feedback  = 0;


void update_emitter_send_info(){
  #ifdef SPEED_CONTROL
    emmiter_send_buf[0] = wb_motor_get_torque_feedback(left_front_wheel);
    emmiter_send_buf[1] = wb_motor_get_torque_feedback(right_front_wheel);
    emmiter_send_buf[2] = wb_motor_get_torque_feedback(left_rear_wheel);
    emmiter_send_buf[3] = wb_motor_get_torque_feedback(right_rear_wheel);
  #else
  //emmiter_send_buf[0] = friction_frac;

    emmiter_send_buf[0] = torque;
    emmiter_send_buf[1] = torque;
    emmiter_send_buf[2] = torque;
    emmiter_send_buf[3] = torque;
  #endif
  
  emmiter_send_buf[4] = wb_position_sensor_get_value(left_front_encoder);
  emmiter_send_buf[5] = wb_position_sensor_get_value(right_front_encoder);
  emmiter_send_buf[6] = wb_position_sensor_get_value(left_rear_encoder);
  emmiter_send_buf[7] = wb_position_sensor_get_value(right_rear_encoder);

  emmiter_send_buf[8] = steering_angle;
  emmiter_send_buf[9] = steering_angle;
}

void blink_lights() {
  int on = (int)wb_robot_get_time() % 2;
  wb_led_set(left_flasher, on);
  wb_led_set(right_flasher, on);
  wb_led_set(tail_lights, on);
  wb_led_set(work_head_lights, on);
  wb_led_set(road_head_lights, on);
}

void print_help() {
  printf("CIRCLE TEST");
}

// set target speed
void set_speed(double kmh) {

  double front_ang_vel = kmh * 1000.0 / 3600.0 / FRONT_WHEEL_RADIUS;
  double rear_ang_vel = kmh * 1000.0 / 3600.0 / REAR_WHEEL_RADIUS;

  // set motor rotation speed
  wb_motor_set_velocity(left_front_wheel, front_ang_vel);
  wb_motor_set_velocity(right_front_wheel, front_ang_vel);
  wb_motor_set_velocity(left_rear_wheel, rear_ang_vel);
  wb_motor_set_velocity(right_rear_wheel, rear_ang_vel);
}

void set_torque(double torq) {

  if (torq > 200.0)
    torq = 200.0;

  torque = torq;

  printf("setting torque to %g N*m\n", torq);

  wb_motor_set_torque(left_front_wheel,torq);
  wb_motor_set_torque(right_front_wheel,torq);
  wb_motor_set_torque(left_rear_wheel,torq);
  wb_motor_set_torque(right_rear_wheel,torq);
}

double target_speed = 6.0; //kmh
void compute_gps_speed() {
  const double *coords = wb_gps_get_values(gps);
  double vel[3] = {coords[X] - gps_coords[X], coords[Y] - gps_coords[Y], coords[Z] - gps_coords[Z]};
  double dist = sqrt(vel[X] * vel[X] + vel[Y] * vel[Y] + vel[Z] * vel[Z]);

  // store into global variables
  gps_speed = dist / time_step * 3600.0;
  memcpy(gps_coords, coords, sizeof(gps_coords));

  printf("target speed: %f km/h\n", target_speed);
  printf("current speed: %g km/h\n", gps_speed);
}


void check_keyboard() {
  int key = wb_keyboard_get_key();
  switch (key) {
    case ' ':
      #ifdef SPEED_CONTROL
        set_speed(0.0);
      #else
        set_torque(0.0);
      #endif  
  }
}

double max_radius = 0;
void calc_max_radius(double current_radius){
  if(current_radius > max_radius)
    max_radius = current_radius;
}

double min_radius = 100; // this should change if radius is changed!!
void calc_min_radius(double current_radius){
  if(current_radius < min_radius)
    min_radius = current_radius;
}

bool circle_started = false;
bool end_test = true;
bool start_record = false;

int end_counter = 0;
double original_radius = 0;
double max_error = 0;
double time = 0;
double radian = 0;
bool flip_signs = false;

void circle_test() {
  
  const double *coords = wb_gps_get_values(gps);
  
  //Print if you want to see these / debug
  //printf("Z-coords: %f \n", coords[Z]);
  //printf("X-coords: %f \n", coords[X]);
  
  // coords[Z] > -0.870069 <-- to start turninig at CoM?
  if (coords[Z] > 0 && circle_started == false){ //starts the circle
    steering_angle = -0.2;
    wb_motor_set_position(left_steer, steering_angle);
    wb_motor_set_position(right_steer, steering_angle);
    
    original_radius = sqrt(coords[X]*coords[X] + coords[Z]*coords[Z]);

    circle_started = true;
  }
  
  if(circle_started){
    double current_radius = sqrt(coords[X]*coords[X] + coords[Z]*coords[Z]);
    calc_max_radius(current_radius);
    calc_min_radius(current_radius);
    
    printf("Original Radius [m]: %f \n", original_radius);
    printf("Current Radius: %f \n", current_radius);
    printf("Max radius: %f \n", max_radius);
    printf("Min radius: %f \n", min_radius);
    
    /*
    double error_percentage = (original_radius - current_radius)/original_radius*100;
    printf("error: %f \n", error_percentage);
    if (error_percentage > fabs(max_error))
      max_error = error_percentage;
    printf("max error: %f \n", max_error);
    */
    #ifdef RECORD
      fprintf(file_circle, "%f\t%f\t%f\t%f\n", time, current_radius, coords[X], coords[Z]); // input time, radius, x and z into the file
      printf("recording...\n");
    #endif
    
    //radian = gps_speed/current_radius*time;
    time += time_step/1000.0;
    
    // ends when end counter is = 3 (one full circle)
    if (coords[Z] > 0 && flip_signs == false){
      end_counter++;     
      flip_signs = !flip_signs;
    }else if(coords[Z] < 0 && flip_signs == true){
      end_counter++;     
      flip_signs = !flip_signs;
    }
    
      
    printf("end counter: %i\n", end_counter);
       
      if(end_counter == 3){
        set_speed(0.0);
        #ifdef RECORD
          fclose(file_circle); // close the file after full circle
        #endif
        }
  }
}

int main(int argc, char **argv) {
  wb_robot_init();

  time_step = (int)wb_robot_get_basic_time_step();

  // find wheels
  left_front_wheel = wb_robot_get_device("left_front_wheel");
  right_front_wheel = wb_robot_get_device("right_front_wheel");
  left_rear_wheel = wb_robot_get_device("left_rear_wheel");
  right_rear_wheel = wb_robot_get_device("right_rear_wheel");

#ifdef RECORD
  //if(circle_started) {
    //file_circle = fopen(file_path, "w");
    file_circle = fopen(file_path, "w");
  //}
#endif
  
#ifdef SPEED_CONTROL
  wb_motor_set_position(left_front_wheel, INFINITY);
  wb_motor_set_position(right_front_wheel, INFINITY);
  wb_motor_set_position(left_rear_wheel, INFINITY);
  wb_motor_set_position(right_rear_wheel, INFINITY); 
  wb_motor_enable_torque_feedback(left_front_wheel, time_step);
  wb_motor_enable_torque_feedback(right_front_wheel, time_step);
  wb_motor_enable_torque_feedback(left_rear_wheel, time_step);
  wb_motor_enable_torque_feedback(right_rear_wheel, time_step);
#endif

  //find encoders
  left_front_encoder = wb_motor_get_position_sensor(left_front_wheel);
  right_front_encoder = wb_motor_get_position_sensor(right_front_wheel);
  left_rear_encoder = wb_motor_get_position_sensor(left_rear_wheel);
  right_rear_encoder = wb_motor_get_position_sensor(right_rear_wheel);

  wb_position_sensor_enable(left_front_encoder, time_step);
  wb_position_sensor_enable(right_front_encoder, time_step);
  wb_position_sensor_enable(left_rear_encoder, time_step);
  wb_position_sensor_enable(right_rear_encoder, time_step);

  // get steering motors
  left_steer = wb_robot_get_device("left_steer");
  right_steer = wb_robot_get_device("right_steer");

  // camera device
  camera = wb_robot_get_device("camera");
  wb_camera_enable(camera, time_step);
  camera_width = wb_camera_get_width(camera);
  camera_height = wb_camera_get_height(camera);
  camera_fov = wb_camera_get_fov(camera);

  // SICK sensor
  sick = wb_robot_get_device("Sick LMS 291");
  wb_lidar_enable(sick, time_step);
  sick_width = wb_lidar_get_horizontal_resolution(sick);
  sick_range = wb_lidar_get_max_range(sick);
  sick_fov = wb_lidar_get_fov(sick);

  // initialize gps
  gps = wb_robot_get_device("gps");
  wb_gps_enable(gps, time_step);

  // find lights
  left_flasher = wb_robot_get_device("left_flasher");
  right_flasher = wb_robot_get_device("right_flasher");
  tail_lights = wb_robot_get_device("tail_lights");
  work_head_lights = wb_robot_get_device("work_head_lights");
  road_head_lights = wb_robot_get_device("road_head_lights");

  //touch sensor
  //touch_sensor = wb_robot_get_device("front left wheel");
  //wb_touch_sensor_enable(touch_sensor, time_step);
  //if(!touch_sensor)
  //  printf("Touch sensor not found !\n");

  //emitter
  emitter = wb_robot_get_device("emitter");
  if (!emitter)
    printf("!!! emitter is not available.\n");

  // start engine

#ifdef SPEED_CONTROL
  set_speed(target_speed);  // km/h
#else
  set_torque(40.0);
#endif

  print_help();
  

  // allow to switch to manual control
  wb_keyboard_enable(time_step);

  // main loop
  while (wb_robot_step(time_step) != -1) {
    
    // run circle test
    circle_test();
    
    check_keyboard();

    update_emitter_send_info();
    wb_emitter_send(emitter, emmiter_send_buf, sizeof(emmiter_send_buf));

    // read sensors
    //const unsigned char *camera_image = wb_camera_get_image(camera);
    //const float *sick_data = wb_lidar_get_range_image(sick);

    // update stuff
    compute_gps_speed();
    blink_lights();
    
  }

  wb_robot_cleanup();
  
  return 0;  // ignored
}
