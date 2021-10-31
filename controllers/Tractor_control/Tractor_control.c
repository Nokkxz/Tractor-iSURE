/*
 * Copyright 1996-2019 Cyberbotics Ltd.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * Description:  Boomer tractor example
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

// to be used as array indices
enum { X, Y, Z };

// This needs to be changed if the .wbt model changes
#define FRONT_WHEEL_RADIUS 0.38
#define REAR_WHEEL_RADIUS 0.6

//comment SPEED_CONTROL to use torque control
//#define SPEED_CONTROL

double emmiter_send_buf[10];
double friction_frac = 400.0;

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
  printf("You can drive this vehicle!\n");
  printf("Select the 3D window and then use the cursor keys to:\n");
  printf("[LEFT]/[RIGHT] - steer\n");
  printf("[UP]/[DOWN] - accelerate/slow down\n");
}

// set target speed
void set_speed(double kmh) {
  if (kmh > 30.0)
    kmh = 30.0;

  speed = kmh;

  printf("setting speed to %g km/h\n", kmh);

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

// positive: turn right, negative: turn left
void set_steering_angle(double wheel_angle) {
  steering_angle = wheel_angle;
  wb_motor_set_position(left_steer, steering_angle);
  wb_motor_set_position(right_steer, steering_angle);
}

void change_manual_steer_angle(double inc) {
  double new_manual_steering = manual_steering + inc;
  printf("steer %f\n", new_manual_steering);
  if (new_manual_steering <= 0.94 && new_manual_steering >= -0.94) {
    manual_steering = new_manual_steering;
    set_steering_angle(manual_steering);
  }

  if (manual_steering == 0)
    printf("going straight\n");
  else
    printf("turning %.2f rad (%s)\n", steering_angle, steering_angle < 0 ? "left" : "right");
}

void check_keyboard() {
  int key = wb_keyboard_get_key();
  switch (key) {
    case WB_KEYBOARD_UP:
      #ifdef SPEED_CONTROL
        set_speed(speed + 1.0);
      #else
        set_torque(torque + 1.0);
      #endif
        //friction_frac++;
      break;
    case WB_KEYBOARD_DOWN:
      #ifdef SPEED_CONTROL
        set_speed(speed - 1.0);
      #else
        set_torque(torque - 1.0);
      #endif  
        //friction_frac--;
      break;
    case ' ':
      #ifdef SPEED_CONTROL
        set_speed(0.0);
      #else
        set_torque(0.0);
      #endif  
      break;
    case WB_KEYBOARD_RIGHT:
      change_manual_steer_angle(+0.02);
      break;
    case WB_KEYBOARD_LEFT:
      change_manual_steer_angle(-0.02);
      break;
    case 'C': // center the steeriing
      change_manual_steer_angle(-manual_steering);
      break;
  }
}

void compute_gps_speed() {
  const double *coords = wb_gps_get_values(gps);
  double vel[3] = {coords[X] - gps_coords[X], coords[Y] - gps_coords[Y], coords[Z] - gps_coords[Z]};
  double dist = sqrt(vel[X] * vel[X] + vel[Y] * vel[Y] + vel[Z] * vel[Z]);

  // store into global variables
  gps_speed = dist / time_step * 3600.0;
  memcpy(gps_coords, coords, sizeof(gps_coords));

  //printf("current speed: %g km/h\n", gps_speed);
}

//void get_wheel_load(){
//  wheel_force_3d = wb_touch_sensor_get_values(touch_sensor);
//  touch_type = wb_touch_sensor_get_type(touch_sensor);
//  printf("Type: %d Force %f, %f, %f\n", touch_type, wheel_force_3d[0], wheel_force_3d[1], wheel_force_3d[2]);
//}

//void print_speed_torque() {
//  printf("left_front_wheel-Current torque: %g ,Traget velocity: %g rad/s\n", left_front_wheel_torque_feedback, left_front_wheel_velocity_feedback);
//}

int main(int argc, char **argv) {
  wb_robot_init();

  time_step = (int)wb_robot_get_basic_time_step();

  // find wheels
  left_front_wheel = wb_robot_get_device("left_front_wheel");
  right_front_wheel = wb_robot_get_device("right_front_wheel");
  left_rear_wheel = wb_robot_get_device("left_rear_wheel");
  right_rear_wheel = wb_robot_get_device("right_rear_wheel");

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
  set_speed(0.0);  // km/h
#else
  set_torque(0.0);
#endif

  print_help();

  // allow to switch to manual control
  wb_keyboard_enable(time_step);

  // main loop
  while (wb_robot_step(time_step) != -1) {
    // get user input
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
