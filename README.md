# AgJunction-Project-Webots

Webots project for AgJunction, by University of Notre Dame and Southern University of Science and Technology

## Version

* Webots Version: R2019b

## Webots Source Code

<https://github.com/omichel/webots>

## Documentation

* [Webots User Guide](https://www.cyberbotics.com/doc/guide/foreword)

* [Webots Reference Manual](https://cyberbotics.com/doc/reference/thanks)

## Get Webots and Preparation

You can get the download link on their offical website: <https://cyberbotics.com>

* [Installing instruction](https://www.cyberbotics.com/doc/guide/installing-webots)

* [Start-up guide](https://www.cyberbotics.com/doc/guide/getting-started-with-webots)

* [development environments setup](https://www.cyberbotics.com/doc/guide/development-environments)

* [Tutorials](https://www.cyberbotics.com/doc/guide/tutorials)

* [Physics plugin](https://cyberbotics.com/doc/reference/physics-plugin)

> If you want to build Webots from source code or contribute to Webots, check out the [repository](https://github.com/omichel/webots) and [wiki](https://github.com/omichel/webots/wiki#installation-of-the-webots-development-environment)

## Setup physics plugin

This part is under the assumption that you have read the tutorials and know the basic operations in Webots.

1. Create a new project in Webots

    > Wizards -> new project directory

2. Download this project
3. Close Webots and copy all the directories into the project folder newly created and override everything.
4. Restart Webots and the simulation is expected to run normally

    > Be sure to hit the "Reload World" button every time you want to reset the simulation. Otherwise the physics plugin will not be reload and you can get unexpected outcome.

5. you can choose different physics plugin here

    > WorldInfo -> physics

6. you can choose different tractor controller plugin here

    > tractor -> controller

**Note that you need to recompile the controller/physics plugin after editing it.**

## Use the controller

In the project we provide a tractor controller that can be used with any of the three physics plugin, the control method is as following:

1. Select the 3D window
2. Use the cursor keys to control

    * [LEFT]/[RIGHT] - steer
    * [UP]/[DOWN] - increase/decrease torque or target speed

> Note that if you are using speed control, then the motor torque is computed automatically by Webots and it is not PID control. the algrithm within needs further exploration. Webots offical documentation about it can be found [here](https://cyberbotics.com/doc/reference/motor#velocity-control).

To change the controller mode between torque control and speed control, simply comment/uncomment the macro definition in the source code and recompile

```c
//comment SPEED_CONTROL to use torque control
#define SPEED_CONTROL
```

## Physics plugin

The project provides 3 implementations for the theoretical models from literature review. **More detailed information about the models and implementation can be found in our literature review summary**

### Wismer_Luth_model

This plugin is based on Wismer-Luth model. Wismer-Luth model is the choice for our first implementation because of the following advantages:

* Easy math form (friendly to programming and computation)
* Better approximation for rolling resistance.

But there are several problems with the implementation

* When using dynamic load, the friction coefficient calculated from slip ratio is very unstable because differentiation is sensitive to noisy data. Use a constant approximation load can make the output stable. But the load is no longer dynamic and accurate. A 5-point derivative or a gaussian filter might help with it (but we did not validate it).

* Wismer-Luth model is under the assumption that the tractor wheel is running at a constant speed. However, torque control is hard to result in a constant speed and speed control in Webots will result in unstable motor torque. Besides, acceleration is not considered but is actually in simulation.

* slip ratio calculated by the model is lower than expected. According to the literature, normally there should be a slip ratio around 10%. But the outcome is lower than 0.001

> Note that this plugin do not support steering and moving backward, which can result in unexpected outcome. This is mainly bacause the angular velocity and horizontal velocity are gotton from z-axis in world coornidate. So is the applied force. Transform coordinate from world to the tractor or the wheel should solve this problem.

### direct_force_method

This plugin is also based on Wismer-Luth model. The main difference is that this approach removes the friction given by the Webots and replace it with our own ones, which partly get rid of the unknown calculation process within Webots.

Another improvement we made in this version is that we can apply forces on the right spot with a realistic direction. This was achieved by using some solid geometry and vector transformations.

### Bekker_Model

Bekker model is an advanced model compared to Wismer-Luth

* It has a larger range of applicable situations. And many other models are based on this one. The accuracy and feasibility of this model had been proved by many other researchers.

* And another advantage of Bekker’s model is that we have found some of the real soil parameters for this model. While for Wismer-Luth model, we have few. Wismer-Luth model is supposed to have a better prediction on the soil the original test were based on, but Bekker’s model could handle more types of soil.

* And the math form is more complicated and realistic. Wismer-Luth should be applied under constant speed, but Bekker’s model didn’t have this assumption.

The problems are as the following

* The numerical results of those variables are realistic. However, Webots would crash to desktop when these forces are applied on the tractor.the problem might be that there’s too much calculations. Further optimization and better algorithms could help solve this problem. However, the physics plugin is written in C/C++, the lack of proficiency of our team in this programming language might be another cause.

## File Content

### controllers

controller program for the tractor in Webots.

More information about controller programming can be found [here](https://cyberbotics.com/doc/guide/controller-programming).

### libraries

Nothing yet, just a project folder generated by Webots

### plugins

* **physics**: physics plugin for Webots. This is a vital part for the project, check out [here](https://cyberbotics.com/doc/reference/physics-plugin) to know more about physics plugin.

* **remote_control**: remote-control controller plugin used to remote control the real robot. Nothing for now. You can get more info [here](https://www.cyberbotics.com/doc/guide/transfer-to-your-own-robot).

* **robot_windows**:  robot window controller plugin used to display the robot window to help data visualization. Nothing for now. You can get more info [here](https://cyberbotics.com/doc/automobile/robot-window).

### proto

PROTO file lists the fields of the PROTO and defines how these fields impact the underlying object which is defined using base nodes and/or PROTO nodes.

For more information about PROTO, click [here](https://cyberbotics.com/doc/reference/proto).

### worlds

Webots world files.

## Author

### Yonglin Jing (Dale)

Main author for the physics plugin program structure and Wismer_Luth_model.

* University: Southern University of Science and Technology
* GitHub page: <https://github.com/NERanger>
* E-mail: 11712605@mail.sustech.edu.cn

### Ranbao Deng (Rambo)

Edit plugin algrithm to get direct_force_method and Bekker_Model

* University: Southern University of Science and Technology
* Email: 11712002@mail.sustech.edu.cn

## Get Support

* You can post the question on [Webots Discord Community](https://discordapp.com/invite/nTWbN9m), some Webots stuff will help you.
* E-mail the author.

> If you have any question regarding program structure and Wismer_Luth_model you can e-mail Dale. Questions about algrithm in direct_force_method and Bekker_Model are recommanded to be sent to Rambo.
