#pragma once

#include <filesystem>
#include <vector>
#include <string>
#include <mutex>
#include <chrono>

#include <sstream>
#include <iomanip>

#include "absl/status/status.h"
#include "absl/log/absl_check.h"

#include "mujoco/mujoco.h"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "osqp++.h"
#include "osqp.h"
#include "GLFW/glfw3.h"

#include "rclcpp/rclcpp.hpp"
#include "osc_2_in_interface/msg/osc_mujoco_state.hpp"
#include "osc_2_in_interface/msg/osc_torque_command.hpp"
#include "osc_2_in_interface/msg/osc_taskspace_targets.hpp"

#include "std_msgs/msg/float64_multi_array.hpp"

// #include "osc_2_in_interface/msg/command.hpp"
// #include "osc_2_in_interface/msg/motor_command.hpp"
// #include "osc_2_in_interface/msg/wheel_motor_command.hpp"

#include "walter_msgs/msg/command.hpp"
#include "walter_msgs/msg/motor_command.hpp"
#include "walter_msgs/msg/wheel_motor_command.hpp"


#include "operational-space-control/walter_sr_v2/utilities.h"
#include "operational-space-control/utilities.h"
#include "operational-space-control/walter_sr_v2/autogen/autogen_functions.h"
#include "operational-space-control/walter_sr_v2/autogen/autogen_defines.h"
#include "operational-space-control/walter_sr_v2/aliases.h"
#include "operational-space-control/walter_sr_v2/constants.h"
#include "operational-space-control/walter_sr_v2/containers.h"


using namespace operational_space_controller::constants;
using namespace operational_space_controller::containers;
using namespace operational_space_controller::aliases;
using namespace osqp;

using rclcpp::Node;
using osc_2_in_interface::msg::OSCMujocoState;
using osc_2_in_interface::msg::OSCTorqueCommand;
using osc_2_in_interface::msg::OSCTaskspaceTargets;

using walter_msgs::msg::Command;
using walter_msgs::msg::MotorCommand;
using walter_msgs::msg::WheelMotorCommand;


// namespace {
//     FunctionOperations Aeq_ops{
//         .incref=Aeq_incref, .checkout=Aeq_checkout, .eval=Aeq, .release=Aeq_release, .decref=Aeq_decref};
//     FunctionOperations beq_ops{
//         .incref=beq_incref, .checkout=beq_checkout, .eval=beq, .release=beq_release, .decref=beq_decref};
//     FunctionOperations Aineq_ops{
//         .incref=Aineq_incref, .checkout=Aineq_checkout, .eval=Aineq, .release=Aineq_release, .decref=Aineq_decref};
//     FunctionOperations bineq_ops{
//         .incref=bineq_incref, .checkout=bineq_checkout, .eval=bineq, .release=bineq_release, .decref=bineq_decref};
//     FunctionOperations H_ops{
//         .incref=H_incref, .checkout=H_checkout, .eval=H, .release=H_release, .decref=H_decref};
//     FunctionOperations f_ops{
//         .incref=f_incref, .checkout=f_checkout, .eval=f, .release=f_release, .decref=f_decref};
//     using AeqParams = FunctionParams<Aeq_SZ_ARG, Aeq_SZ_RES, Aeq_SZ_IW, Aeq_SZ_W, optimization::Aeq_rows, optimization::Aeq_cols, optimization::Aeq_sz, 4>;
//     using beqParams = FunctionParams<beq_SZ_ARG, beq_SZ_RES, beq_SZ_IW, beq_SZ_W, optimization::beq_sz, 1, optimization::beq_sz, 4>;
//     using AineqParams = FunctionParams<Aineq_SZ_ARG, Aineq_SZ_RES, Aineq_SZ_IW, Aineq_SZ_W, optimization::Aineq_rows, optimization::Aineq_cols, optimization::Aineq_sz, 1>;
//     using bineqParams = FunctionParams<bineq_SZ_ARG, bineq_SZ_RES, bineq_SZ_IW, bineq_SZ_W, optimization::bineq_sz, 1, optimization::bineq_sz, 1>;
//     using HParams = FunctionParams<H_SZ_ARG, H_SZ_RES, H_SZ_IW, H_SZ_W, optimization::H_rows, optimization::H_cols, optimization::H_sz, 4>;
//     using fParams = FunctionParams<f_SZ_ARG, f_SZ_RES, f_SZ_IW, f_SZ_W, optimization::f_sz, 1, optimization::f_sz, 4>;
// }

class OSCNode : public Node {
public:
    OSCNode(const std::string& xml_path);
    ~OSCNode();

    void stop_robot(); // New helper to send zero commands                                    


private:
    // ROS 2 Callbacks
    void state_callback(const OSCMujocoState::SharedPtr msg);
    void timer_callback();

    // Internal helper functions (REFACTORED SIGNATURES)
    // Now takes the state to eliminate reliance on the shared state_ member.
    void update_mj_data(const State& current_state);
    void update_osc_data();
    // void update_optimization_data();
    void update_optimization_data(const Eigen::Vector<double, model::contact_site_ids_size>& force_limits);

    absl::Status set_up_optimization(const Vector<model::contact_site_ids_size>& force_limits);
    // Now takes contact mask data to update bounds
    absl::Status update_optimization(const Vector<model::contact_site_ids_size>& force_limits); 
    // void solve_optimization();
    bool solve_optimization();
    void reset_optimization();
    // Now takes the safety status
    void publish_torque_command(bool safety_override_active_local, 
                                    std::chrono::time_point<std::chrono::high_resolution_clock> state_read_time_local);
    // Helper functions from the original main
    
    template <typename T>
    bool contains(const std::vector<T>& vec, const T& value);
    std::vector<int> getSiteIdsOnSameBodyAsGeom(const mjModel* m, int geom_id);
    std::vector<int> getBinaryRepresentation_std_find(const std::vector<int>& A, const std::vector<int>& B);

    // ROS 2 members
    rclcpp::Subscription<OSCMujocoState>::SharedPtr state_subscriber_;
    // rclcpp::Subscription<OSCTaskspaceTargets>::SharedPtr taskspace_targets_subscriber_;
    // rclcpp::Publisher<OSCTorqueCommand>::SharedPtr torque_publisher_;
    rclcpp::Publisher<Command>::SharedPtr torque_publisher_;

    // Inside OSCNode class private members:
    rclcpp::Publisher<std_msgs::msg::Float64MultiArray>::SharedPtr data_pub_;    
    std_msgs::msg::Float64MultiArray data_msg_; // Persistent message

    rclcpp::TimerBase::SharedPtr timer_;

    // Mutexes for thread-safe access
    std::mutex state_mutex_;
    // std::mutex taskspace_targets_mutex_;


    
    // --- Persistent PD Control History & References ---
    double last_time_ = 0.0;
    
    bool safety_override_active_ = false;

    // Initial joint angles (for position tracking reference)
    double initial_tl_angular_position_ = 0.0;
    double initial_tr_angular_position_ = 0.0;
    double initial_hl_angular_position_ = 0.0;
    double initial_hr_angular_position_ = 0.0;
    double initial_tlh_angular_position_ = 0.0;
    double initial_trh_angular_position_ = 0.0;
    double initial_hlh_angular_position_ = 0.0;
    double initial_hrh_angular_position_ = 0.0;



    // Shared Variables (Inputs/Outputs)
    State state_;
    Matrix<model::site_ids_size, 6> taskspace_targets_;
    Vector<model::nu_size> torque_command_;

    Vector<3> initial_position_;    

    // Mujoco Variables
    mjModel* mj_model_;
    mjData* mj_data_;
    std::string xml_path_;
    std::vector<std::string> sites_;
    std::vector<int> site_ids_;
    std::vector<int> wheel_sites_mujoco_ = {3, 4, 7, 8, 11, 12, 15, 16};
    std::vector<std::string> bodies_;
    std::vector<std::string> noncontact_sites_;
    std::vector<std::string> contact_sites_;
    std::vector<int> noncontact_site_ids_;
    std::vector<int> contact_site_ids_;
    std::vector<int> body_ids_;
    Matrix<model::site_ids_size, 3> points_;
    static constexpr bool is_fixed_based = false;

    // OSQP Solver, settings, and matrices
    OsqpInstance instance_;
    OsqpSolver solver_;
    OsqpSettings settings_;
    OsqpExitCode exit_code_;
    Vector<optimization::design_vector_size> solution_;
    Vector<optimization::constraint_matrix_rows> dual_solution_;
    Vector<optimization::design_vector_size> design_vector_;
    double infinity_ = OSQP_INFTY;
    OSCData osc_data_;
    OptimizationData opt_data_;
    float big_number_ = 1e4;

    // Constraints
    MatrixColMajor<optimization::design_vector_size, optimization::design_vector_size> Abox_;
    Vector<optimization::dv_size> dv_lb_;
    Vector<optimization::dv_size> dv_ub_;
    Vector<model::nu_size> u_lb_;
    Vector<model::nu_size> u_ub_;
    Vector<optimization::z_size> z_lb_;
    Vector<optimization::z_size> z_ub_;
    Vector<optimization::bineq_sz> bineq_lb_;

    // Time storage for latency calculation
    std::chrono::time_point<std::chrono::high_resolution_clock> state_read_time_;
    std::chrono::time_point<std::chrono::high_resolution_clock> torque_ready_time_;

    // Storage for the motor position detected at state read time
    Vector<model::nu_size> last_detected_motor_position_;    
    bool is_state_received_ = false;    
    double time_mujoco_update_ms_ = 0.0;
    double time_casadi_update_ms_ = 0.0;
    double time_osqp_solve_ms_ = 0.0;
    double time_wait_for_execution_ms_ = 0.0;        
};