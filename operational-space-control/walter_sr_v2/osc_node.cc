#include "operational-space-control/walter_sr_v2/osc_node.h"
#include <chrono>
#include <cmath>
#include <algorithm>

using namespace std::chrono_literals;

const int NUM_MOTORS = 8;
const int CONTACT_SITES = 8; 

OSCNode::OSCNode(const std::string& xml_path)
    : Node("osc_node")
{
    // --- ROS 2 communication setup ---
    state_subscriber_ = this->create_subscription<OSCMujocoState>(
        "/state_estimator/state", 1, std::bind(&OSCNode::state_callback, this, std::placeholders::_1));
    
    torque_publisher_ = this->create_publisher<Command>("walter/command", 1);
    
    auto data_qos = rclcpp::QoS(rclcpp::KeepLast(10));
    data_qos.reliability(RMW_QOS_POLICY_RELIABILITY_BEST_EFFORT);
    data_pub_ = this->create_publisher<std_msgs::msg::Float64MultiArray>("walter/data", data_qos);

    // Keep MATLAB layout perfectly intact (54 Elements)
    data_msg_.layout.dim.resize(4);
    data_msg_.layout.dim[0].label = "Telemetry Data Array (54 elements)"; 
    data_msg_.data.reserve(54); 

    timer_ = this->create_wall_timer(std::chrono::microseconds(5000), std::bind(&OSCNode::timer_callback, this));

    rclcpp::on_shutdown([this]() {
        RCLCPP_WARN(this->get_logger(), "Shutdown signal received. Attempting to stop robot...");
        this->stop_robot();
    });
    
    RCLCPP_INFO(this->get_logger(), "Pure Joint-Space PD Node Initialized. MuJoCo Completely Bypassed.");
}

OSCNode::~OSCNode() {
    // MuJoCo pointers completely removed. Nothing to clean up.
}

void OSCNode::state_callback(const OSCMujocoState::SharedPtr msg) {
    std::lock_guard<std::mutex> lock(state_mutex_);
    for (size_t i = 0; i < NUM_MOTORS; ++i) {
        state_.motor_position(i) = static_cast<double>(msg->motor_position[i]);
        state_.motor_velocity(i) = static_cast<double>(msg->motor_velocity[i]);
        state_.torque_estimate(i) = static_cast<double>(msg->torque_estimate[i]);
    }
    state_read_time_ = std::chrono::high_resolution_clock::now();    
    
    for (size_t i = 0; i < 4; ++i) state_.body_rotation(i) = static_cast<double>(msg->body_rotation[i]);
    for (size_t i = 0; i < 3; ++i) {
        state_.linear_body_velocity(i) = static_cast<double>(msg->linear_body_velocity[i]);
        state_.angular_body_velocity(i) = static_cast<double>(msg->angular_body_velocity[i]);
    }
    for (size_t i = 0; i < CONTACT_SITES; ++i) state_.contact_mask(i) = static_cast<double>(msg->contact_mask[i]);
    is_state_received_ = true;    
}

void OSCNode::timer_callback() {
    State local_state; 
    bool local_safety_override_active;
    std::chrono::time_point<std::chrono::high_resolution_clock> local_state_read_time;
    double current_time = this->now().seconds();

    auto t_start_execution = std::chrono::high_resolution_clock::now(); 

    {
        std::lock_guard<std::mutex> lock_state(state_mutex_);
        local_state = state_; 
        local_safety_override_active = safety_override_active_;
        local_state_read_time = state_read_time_;
        if (!is_state_received_) return; 
    } 

    time_wait_for_execution_ms_ = std::chrono::duration<double, std::milli>(t_start_execution - local_state_read_time).count();    

    // --- STATIC PERSISTENT VARIABLES ---
    static double shin_pos_0_initial = 0.0;
    static double shin_pos_1_initial = 0.0;
    static double shin_pos_2_initial = 0.0;
    static double shin_pos_3_initial = 0.0;    
    static double gait_start_time = 0.0;

    // Check for first call 
    if (last_time_ == 0.0) {
        gait_start_time = current_time;

        // Capture initial positions of the 4 knees (Odds: 1, 3, 5, 7)
        shin_pos_0_initial  = local_state.motor_position(1);
        shin_pos_1_initial  = local_state.motor_position(3);
        shin_pos_2_initial  = local_state.motor_position(5);
        shin_pos_3_initial  = local_state.motor_position(7);

        RCLCPP_INFO(this->get_logger(), "Starting Pure PD Tumbling Trajectory...");
        last_time_ = current_time;
        return; 
    }

    // --- SAFETY LIMITS ---
    bool limit_hit = local_safety_override_active; 
    if (!local_safety_override_active) {
        const double HIP_MIN_RAD = -0.8; // Expanded safely to allow deep robot tucks
        const double HIP_MAX_RAD = 1.8;
        
        for (size_t i : {0, 2, 4, 6}) { // Hips only
            double pos = local_state.motor_position(i);
            if (pos < HIP_MIN_RAD || pos > HIP_MAX_RAD) {
                limit_hit = true;
                RCLCPP_WARN_ONCE(this->get_logger(), "Hip limit hit on motor index %zu (val: %.2f). Overriding.", i, pos);
                break; 
            }
        }
        if (limit_hit) { 
            std::lock_guard<std::mutex> lock_state(state_mutex_);
            safety_override_active_ = true;
            local_safety_override_active = true; 
        }
    }

    std::vector<double> commanded_torques(NUM_MOTORS, 0.0);
    double shin_vel_target = 0.0;
    double pos_offset = 0.0;

    // --- PURE JOINT SPACE PD CONTROLLER ---
    if (!local_safety_override_active) {
        
        // Trajectory Generation
        double elapsed_t = current_time - gait_start_time;
        double MAX_SHIN_VEL = 0.1; // Flipping speed
        double RAMP_TIME = 1.0;    

        if (elapsed_t < RAMP_TIME) {
            shin_vel_target = MAX_SHIN_VEL * (elapsed_t / RAMP_TIME);
            pos_offset = 0.5 * (MAX_SHIN_VEL / RAMP_TIME) * (elapsed_t * elapsed_t);
        } else {
            shin_vel_target = MAX_SHIN_VEL;
            double distance_during_ramp = 0.5 * MAX_SHIN_VEL * RAMP_TIME;
            pos_offset = distance_during_ramp + MAX_SHIN_VEL * (elapsed_t - RAMP_TIME);
        }

        // --- PD Gains ---
        double target_hip_pos = 0.0; // Straight down relative to chassis
        double target_hip_vel = 0.0;
        
        double hip_kp = 10.0; // Stiff enough to hold the body
        double hip_kd = 1.0;
        
        double shin_kp = 10.0; // Aggressive tracking for the spin
        double shin_kd = 1.0;
        
        for (int i = 0; i < NUM_MOTORS; ++i) {
            bool is_hip = (i % 2 == 0);
            double current_pos = local_state.motor_position(i);
            double current_vel = local_state.motor_velocity(i);
            
            if (is_hip) {
                // Hip Control: Hold 0.0
                commanded_torques[i] = hip_kp * (target_hip_pos - current_pos) + hip_kd * (target_hip_vel - current_vel);
            } else {
                // Shin Control: Sweep backwards
                double initial_shin = (i==1) ? shin_pos_0_initial : (i==3) ? shin_pos_1_initial : (i==5) ? shin_pos_2_initial : shin_pos_3_initial;
                double target_shin_pos = initial_shin + pos_offset;
                commanded_torques[i] = shin_kp * (target_shin_pos - current_pos) + shin_kd * (shin_vel_target - current_vel);
            }
        }

        // ==============================================================================
        // --- TELEMETRY DUMP (Strictly padded to 54 elements for MATLAB) ---
        // ==============================================================================
        data_msg_.data.clear();
        
        // [0-9] Hip Z targets and velocities (Unused)
        for(int i=0; i < 10; i++) data_msg_.data.push_back(0.0); 

        // [10-13] Target Shin Positions
        data_msg_.data.push_back(shin_pos_0_initial + pos_offset);
        data_msg_.data.push_back(shin_pos_1_initial + pos_offset);
        data_msg_.data.push_back(shin_pos_2_initial + pos_offset);
        data_msg_.data.push_back(shin_pos_3_initial + pos_offset);

        // [14-17] Actual Shin Positions (Odd indices)
        data_msg_.data.push_back(local_state.motor_position(1));
        data_msg_.data.push_back(local_state.motor_position(3));
        data_msg_.data.push_back(local_state.motor_position(5));
        data_msg_.data.push_back(local_state.motor_position(7));

        // [18] Target Shin Velocity
        data_msg_.data.push_back(shin_vel_target);

        // [19-22] Actual Shin Velocities (Odd indices)
        data_msg_.data.push_back(local_state.motor_velocity(1));
        data_msg_.data.push_back(local_state.motor_velocity(3));
        data_msg_.data.push_back(local_state.motor_velocity(5));
        data_msg_.data.push_back(local_state.motor_velocity(7));

        // [23-25] Body X, Y, Z (Unused without MuJoCo)
        data_msg_.data.push_back(0.0); 
        data_msg_.data.push_back(0.0); 
        data_msg_.data.push_back(0.0); 
        
        // [26-33] Contact Mask (8 elements)
        for (int i = 0; i < 8; ++i) data_msg_.data.push_back(local_state.contact_mask(i));
        
        // [34-37] Solver Health (Unused)
        data_msg_.data.push_back(0.0); // OSQP Exit Code
        data_msg_.data.push_back(0.0); // Time CasADi
        data_msg_.data.push_back(0.0); // Time OSQP
        data_msg_.data.push_back(0.0); // QP Obj

        // [38-45] Commanded Torques (RLH, RLK, RRH, RRK, FLH, FLK, FRH, FRK)
        for (int i = 0; i < NUM_MOTORS; ++i) data_msg_.data.push_back(commanded_torques[i]);
        
        // [46-53] FZ Ground Reaction Forces (Unused)
        for (int i = 0; i < CONTACT_SITES; ++i) data_msg_.data.push_back(0.0); 
        
        // Safety Check: Ensure exactly 54 elements
        assert(data_msg_.data.size() == 54 && "Data array size mismatch!");

        data_pub_->publish(data_msg_);
    }

    if (!local_safety_override_active) {
        std::lock_guard<std::mutex> lock(state_mutex_);
        if (safety_override_active_) local_safety_override_active = true; 
    }    

    last_time_ = current_time;
    publish_torque_command(local_safety_override_active, local_state_read_time, commanded_torques); 
}

void OSCNode::publish_torque_command(bool safety_override_active_local, 
                                     std::chrono::time_point<std::chrono::high_resolution_clock> state_read_time_local,
                                     const std::vector<double>& torques) 
{
    const std::set<std::string> reversed_joints_ = {
        "rear_left_hip", "rear_left_knee", "front_left_hip", "front_left_knee"};
    const std::array<std::string, NUM_MOTORS> MOTOR_NAMES = {
        "rear_left_hip", "rear_left_knee", "rear_right_hip", "rear_right_knee",
        "front_left_hip", "front_left_knee", "front_right_hip", "front_right_knee"};
    
    const double MAX_TORQUE = 8.0; // Safe Joint Limit
    
    auto command_msg = std::make_unique<Command>(); 
    command_msg->master_gain = 1.0; 
    command_msg->motor_commands.resize(NUM_MOTORS);

    if (safety_override_active_local) {
        command_msg->high_level_control_mode = 2; 
        for (size_t i = 0; i < NUM_MOTORS; ++i) {
            command_msg->motor_commands[i].name = MOTOR_NAMES[i];
            command_msg->motor_commands[i].control_mode = 2; // VELOCITY_CONTROL
            command_msg->motor_commands[i].position_setpoint = 0.0; 
            command_msg->motor_commands[i].velocity_setpoint = 0.0;
            command_msg->motor_commands[i].feedforward_torque = 0.0; 
            command_msg->motor_commands[i].kp = 0.0; 
            command_msg->motor_commands[i].kd = 0.0;
            command_msg->motor_commands[i].input_mode = 1;   
            command_msg->motor_commands[i].enable = true; 
        }
    } else {
        command_msg->high_level_control_mode = 2;
        
        std::stringstream ss;
        ss << "PD Torques: [ ";
        for (size_t i = 0; i < NUM_MOTORS; ++i) {
            ss << std::fixed << std::setprecision(2) << torques[i] << " ";
        }
        ss << "]";
        RCLCPP_INFO_THROTTLE(this->get_logger(), *this->get_clock(), 500, "%s", ss.str().c_str());

        for (size_t i = 0; i < NUM_MOTORS; ++i) {
            double final_torque = torques[i];
            
            // Apply sign flip for URDF mirroring
            if (reversed_joints_.count(MOTOR_NAMES[i])) final_torque *= -1.0;
            
            // Hard clamp before network sending
            final_torque = std::clamp(final_torque, -MAX_TORQUE, MAX_TORQUE);

            command_msg->motor_commands[i].name = MOTOR_NAMES[i];
            command_msg->motor_commands[i].control_mode = 1; // TORQUE_CONTROL
            command_msg->motor_commands[i].feedforward_torque = final_torque; 
            command_msg->motor_commands[i].position_setpoint = 0.0;
            command_msg->motor_commands[i].velocity_setpoint = 0.0; 
            command_msg->motor_commands[i].kp = 0.0; 
            command_msg->motor_commands[i].kd = 0.0;
            command_msg->motor_commands[i].input_mode = 1;   
            command_msg->motor_commands[i].enable = true; 
        }
    }
    
    std::chrono::high_resolution_clock::time_point torque_ready_time = std::chrono::high_resolution_clock::now();
    auto latency = std::chrono::duration_cast<std::chrono::microseconds>(torque_ready_time - state_read_time_local);
    double latency_ms = static_cast<double>(latency.count()) / 1000.0;

    torque_publisher_->publish(std::move(command_msg));

    if (safety_override_active_local) {
        RCLCPP_WARN_THROTTLE(this->get_logger(), *this->get_clock(), 1000, "Safety Override Active. Latency: %.3f ms", latency_ms);
    }
}

void OSCNode::stop_robot() {
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        safety_override_active_ = true;        
    }    
    
    auto command_msg = std::make_unique<Command>(); 
    command_msg->master_gain = 1.0; 
    command_msg->motor_commands.resize(NUM_MOTORS);
    command_msg->high_level_control_mode = 2; 

    const std::array<std::string, NUM_MOTORS> MOTOR_NAMES = {
        "rear_left_hip", "rear_left_knee", "rear_right_hip", "rear_right_knee",
        "front_left_hip", "front_left_knee", "front_right_hip", "front_right_knee"};

    for (size_t i = 0; i < NUM_MOTORS; ++i) {
        command_msg->motor_commands[i].name = MOTOR_NAMES[i];
        command_msg->motor_commands[i].control_mode = 2; 
        command_msg->motor_commands[i].position_setpoint = 0.0;
        command_msg->motor_commands[i].velocity_setpoint = 0.0;
        command_msg->motor_commands[i].feedforward_torque = 0.0; 
        command_msg->motor_commands[i].kp = 0.0; 
        command_msg->motor_commands[i].kd = 0.0;
        command_msg->motor_commands[i].input_mode = 1;   
        command_msg->motor_commands[i].enable = true; 
    }

    if (torque_publisher_) {
        torque_publisher_->publish(std::move(command_msg));
        RCLCPP_INFO(this->get_logger(), ">>> SAFETY STOP COMMAND SENT <<<");
    }
}