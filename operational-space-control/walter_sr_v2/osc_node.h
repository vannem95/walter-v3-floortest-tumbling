#pragma once

#include <vector>
#include <string>
#include <mutex>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <filesystem>
#include <set>
#include <array>
#include <memory>

#include "rclcpp/rclcpp.hpp"
#include "osc_2_in_interface/msg/osc_mujoco_state.hpp"
#include "std_msgs/msg/float64_multi_array.hpp"

#include "walter_msgs/msg/command.hpp"
#include "walter_msgs/msg/motor_command.hpp"
#include "walter_msgs/msg/wheel_motor_command.hpp"

// We keep containers/aliases if State struct is defined there, 
// otherwise we can just define a simple local struct.
#include "operational-space-control/walter_sr_v2/aliases.h"
#include "operational-space-control/walter_sr_v2/containers.h"

using namespace operational_space_controller::containers;
using namespace operational_space_controller::aliases;

using rclcpp::Node;
using osc_2_in_interface::msg::OSCMujocoState;
using walter_msgs::msg::Command;
using walter_msgs::msg::MotorCommand;

class OSCNode : public Node {
public:
    // Constructor still takes xml_path just to match your existing launch file API, 
    // even though it is completely ignored internally now.
    OSCNode(const std::string& xml_path);
    ~OSCNode();

    void stop_robot(); 

private:
    // ROS 2 Callbacks
    void state_callback(const OSCMujocoState::SharedPtr msg);
    void timer_callback();

    void publish_torque_command(bool safety_override_active_local, 
                                std::chrono::time_point<std::chrono::high_resolution_clock> state_read_time_local,
                                const std::vector<double>& torques);

    // ROS 2 members
    rclcpp::Subscription<OSCMujocoState>::SharedPtr state_subscriber_;
    rclcpp::Publisher<Command>::SharedPtr torque_publisher_;
    rclcpp::Publisher<std_msgs::msg::Float64MultiArray>::SharedPtr data_pub_;    
    
    std_msgs::msg::Float64MultiArray data_msg_; 
    rclcpp::TimerBase::SharedPtr timer_;

    // Mutexes for thread-safe access
    std::mutex state_mutex_;

    // --- Persistent Control History & Status ---
    double last_time_ = 0.0;
    bool safety_override_active_ = false;
    bool is_state_received_ = false;    

    // --- Shared State Variable ---
    State state_;

    // --- Time storage for latency calculation ---
    std::chrono::time_point<std::chrono::high_resolution_clock> state_read_time_;

    double time_wait_for_execution_ms_ = 0.0;        
    double time_mujoco_update_ms_ = 0.0;
};