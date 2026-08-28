#pragma once

#include <filesystem>
#include <vector>
#include <string>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <iostream>
#include <cassert>

#include <typeinfo>

#include "absl/status/status.h"
#include "absl/log/absl_check.h"

#include "mujoco/mujoco.h"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "osqp++.h"
#include "osqp.h"

#include "operational-space-control/walter_sr/utilities.h"
#include "operational-space-control/utilities.h"

#include "operational-space-control/walter_sr/autogen/autogen_functions.h"
#include "operational-space-control/walter_sr/autogen/autogen_defines.h"

#include "operational-space-control/walter_sr/aliases.h"
#include "operational-space-control/walter_sr/constants.h"
#include "operational-space-control/walter_sr/containers.h"

#include "rclcpp/rclcpp.hpp"
#include "osc_2_in/msg/osc_mujoco_state.hpp"
#include "osc_2_in/msg/osc_torque_command.hpp"


using namespace operational_space_controller::constants;
using namespace operational_space_controller::containers;
using namespace operational_space_controller::aliases;
using namespace osqp;

// ROS 2 namespaces
using osc_2_in::msg::OSCMujocoState;
using osc_2_in::msg::OSCTorqueCommand;

// Anonymous Namespace for shorthand constants:
namespace {
    // Map Casadi Functions to FunctionOperations Struct:
    FunctionOperations Aeq_ops{
        .incref=Aeq_incref,
        .checkout=Aeq_checkout,
        .eval=Aeq,
        .release=Aeq_release,
        .decref=Aeq_decref
    };

    FunctionOperations beq_ops{
        .incref=beq_incref,
        .checkout=beq_checkout,
        .eval=beq,
        .release=beq_release,
        .decref=beq_decref
    };

    FunctionOperations Aineq_ops{
        .incref=Aineq_incref,
        .checkout=Aineq_checkout,
        .eval=Aineq,
        .release=Aineq_release,
        .decref=Aineq_decref
    };

    FunctionOperations bineq_ops{
        .incref=bineq_incref,
        .checkout=bineq_checkout,
        .eval=bineq,
        .release=bineq_release,
        .decref=bineq_decref
    };

    FunctionOperations H_ops{
        .incref=H_incref,
        .checkout=H_checkout,
        .eval=H,
        .release=H_release,
        .decref=H_decref
    };

    FunctionOperations f_ops{
        .incref=f_incref,
        .checkout=f_checkout,
        .eval=f,
        .release=f_release,
        .decref=f_decref
    };

    // Casadi Functions
    using AeqParams = 
        FunctionParams<Aeq_SZ_ARG, Aeq_SZ_RES, Aeq_SZ_IW, Aeq_SZ_W, optimization::Aeq_rows, optimization::Aeq_cols, optimization::Aeq_sz, 4>;
    using beqParams =
        FunctionParams<beq_SZ_ARG, beq_SZ_RES, beq_SZ_IW, beq_SZ_W, optimization::beq_sz, 1, optimization::beq_sz, 4>;
    using AineqParams =
        FunctionParams<Aineq_SZ_ARG, Aineq_SZ_RES, Aineq_SZ_IW, Aineq_SZ_W, optimization::Aineq_rows, optimization::Aineq_cols, optimization::Aineq_sz, 1>;
    using bineqParams =
        FunctionParams<bineq_SZ_ARG, bineq_SZ_RES, bineq_SZ_IW, bineq_SZ_W, optimization::bineq_sz, 1, optimization::bineq_sz, 1>;
    using HParams =
        FunctionParams<H_SZ_ARG, H_SZ_RES, H_SZ_IW, H_SZ_W, optimization::H_rows, optimization::H_cols, optimization::H_sz, 4>;
    using fParams =
        FunctionParams<f_SZ_ARG, f_SZ_RES, f_SZ_IW, f_SZ_W, optimization::f_sz, 1, optimization::f_sz, 4>;
}

//TODO(jeh15): Refactor all voids with absl::Status
class OSCNode : public rclcpp::Node {
    public:
        // Constructor to set up all ROS 2 communication and controller initialization.
        OSCNode(const std::string& xml_path) : Node("osc_node"), xml_path(xml_path) 

        // --- Initialization of variables using initializer list ---
        solution_(Vector<optimization::design_vector_size>::Zero()),
        dual_solution_(Vector<optimization::constraint_matrix_rows>::Zero()),
        design_vector_(Vector<optimization::design_vector_size>::Zero()),
        infinity_(OSQP_INFTY),
        big_number_(1e4),
        Abox_(MatrixColMajor<optimization::design_vector_size, optimization::design_vector_size>::Identity()),
        dv_lb_(Vector<optimization::dv_size>::Constant(-infinity_)),
        dv_ub_(Vector<optimization::dv_size>::Constant(infinity_)),
        u_lb_({
            -1000, -1000,
            -1000, -1000,
            -1000, -1000,
            -1000, -1000
        }),
        u_ub_({
            1000, 1000,
            1000, 1000,
            1000, 1000,
            1000, 1000
        }),
        z_lb_({
            -infinity_, -infinity_, 0.0,
            -infinity_, -infinity_, 0.0,
            -infinity_, -infinity_, 0.0,
            -infinity_, -infinity_, 0.0,
            -infinity_, -infinity_, 0.0,
            -infinity_, -infinity_, 0.0,
            -infinity_, -infinity_, 0.0,
            -infinity_, -infinity_, 0.0
        }),
        z_ub_({
            infinity_, infinity_, big_number_,
            infinity_, infinity_, big_number_,
            infinity_, infinity_, big_number_,
            infinity_, infinity_, big_number_,
            infinity_, infinity_, big_number_,
            infinity_, infinity_, big_number_,
            infinity_, infinity_, big_number_,
            infinity_, infinity_, big_number_
        }),
        bineq_lb_(Vector<optimization::bineq_sz>::Constant(-infinity_))        

        {
            
            // --- Initialization logic from your `initialize` function ---
            char error[1000];
            mj_model = mj_loadXML(xml_path.c_str(), nullptr, error, 1000);
            if (!mj_model) {
                printf("%s\n", error);
                // In a ROS 2 node, it's best to throw an exception or log a fatal error
                // if a critical resource like the model cannot be loaded.
                RCLCPP_FATAL(this->get_logger(), "Failed to load Mujoco Model: %s", error);
                throw std::runtime_error("Failed to load Mujoco Model.");
            }
            mj_model->opt.timestep = 0.002;
            mj_data = mj_makeData(mj_model);
            
            // Populate the site and body ID vectors

            for(const std::string_view& site : model::site_list) {
                std::string site_str = std::string(site);
                int id = mj_name2id(mj_model, mjOBJ_SITE, site_str.data());
                assert(id != -1 && "Site not found in model.");
                sites.push_back(site_str);
                site_ids.push_back(id);
            }
            for(const std::string_view& site : model::noncontact_site_list) {
                std::string site_str = std::string(site);
                int id = mj_name2id(mj_model, mjOBJ_SITE, site_str.data());
                assert(id != -1 && "Site not found in model.");
                noncontact_sites.push_back(site_str);
                noncontact_site_ids.push_back(id);
            }
            for(const std::string_view& site : model::contact_site_list) {
                std::string site_str = std::string(site);
                int id = mj_name2id(mj_model, mjOBJ_SITE, site_str.data());
                assert(id != -1 && "Site not found in model.");
                contact_sites.push_back(site_str);
                contact_site_ids.push_back(id);
            }
            for(const std::string_view& body : model::body_list) {
                std::string body_str = std::string(body);
                int id = mj_name2id(mj_model, mjOBJ_BODY, body_str.data());
                assert(id != -1 && "Body not found in model.");
                bodies.push_back(body_str);
                body_ids.push_back(id);
            }
            
            assert(site_ids.size() == body_ids.size() && "Number of Sites and Bodies must be equal.");
            
            // --- Initialization logic from your `initialize_optimization` function ---
            
            // You'll need to pass an initial state from somewhere. 
            // Create an instance of the message.
            osc_2_in::msg::OSCMujocoState initial_state_msg;

            // Set all fields to zeros. The default constructor usually handles this,
            // but explicit assignment is safer for clarity and consistency.
            initial_state_msg.motor_position.assign(model::nu_size, 0.0);
            initial_state_msg.motor_velocity.assign(model::nu_size, 0.0);
            initial_state_msg.torque_estimate.assign(model::nu_size, 0.0);
            initial_state_msg.body_rotation.assign(4, 0.0);
            initial_state_msg.linear_body_velocity.assign(3, 0.0);
            initial_state_msg.angular_body_velocity.assign(3, 0.0);
            initial_state_msg.contact_mask.assign(model::contact_site_ids_size, 0.0);

            // Initialize optimization ------------------
            update_mj_data(initial_state_msg); // Assuming you have an `update_mj_data` function
            
            // Set up optimization
            absl::Status result = set_up_optimization();
            if (!result.ok()) {
                RCLCPP_FATAL(this->get_logger(), "Failed to initialize optimization: %s", result.message().data());
                throw std::runtime_error("Failed to initialize optimization.");
            }
            // ------------------

            // --- Set up ROS 2 communication (Subscribers, Publishers, Timers) ---
            state_subscriber_ = this->create_subscription<osc_2_in::msg::OSCMujocoState>(
                "osc/mujoco_state", 10,
                std::bind(&OSCNode::state_callback, this, std::placeholders::_1)
            );

            
            torque_publisher_ = this->create_publisher<osc_2_in::msg::OSCTorqueCommand>(
                "osc/torque_command", 10
            );

            // Replace your manual thread with a ROS 2 timer
            timer_ = this->create_wall_timer(
                std::chrono::microseconds(2000), // control_rate_us
                std::bind(&OSCNode::timer_callback, this)
            );
        }

    private:

        // Shared Variables (Inputs/Outputs)
        State state_;
        Matrix<model::site_ids_size, 6> taskspace_targets_;
        Vector<model::nu_size> torque_command_;

        // Mujoco Variables
        mjModel* mj_model_;
        mjData* mj_data_;
        std::string xml_path_;
        std::vector<std::string> sites_;
        std::vector<std::string> bodies_;
        std::vector<std::string> noncontact_sites_;
        std::vector<std::string> contact_sites_;
        std::vector<int> site_ids_;
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
        const double infinity_ = OSQP_INFTY;
        OSCData osc_data_;
        OptimizationData opt_data_;
        const float big_number_ = 1e4;

        // Constraints (as member variables or defined within a function)
        MatrixColMajor<optimization::design_vector_size, optimization::design_vector_size> Abox_;
        Vector<optimization::dv_size> dv_lb_;
        Vector<optimization::dv_size> dv_ub_;
        Vector<model::nu_size> u_lb_;
        Vector<model::nu_size> u_ub_;
        Vector<optimization::z_size> z_lb_;
        Vector<optimization::z_size> z_ub_;
        Vector<optimization::bineq_sz> bineq_lb_;



        absl::Status set_up_optimization() {
                // Initialize the Optimization: (Everything should be Column Major for OSQP)
                // Get initial data from initial state:
                update_osc_data();
                update_optimization_data();

                // Concatenate Constraint Matrix:
                MatrixColMajor<optimization::constraint_matrix_rows, optimization::constraint_matrix_cols> A;
                A << opt_data.Aeq, opt_data.Aineq, Abox;
                // Calculate Bounds:
                Vector<optimization::bounds_size> lb;
                Vector<optimization::bounds_size> ub;
                Vector<optimization::z_size> z_lb_masked = z_lb;
                Vector<optimization::z_size> z_ub_masked = z_ub;
                for(int i = 0; i < model::contact_site_ids_size; i++) {
                    z_lb_masked(Eigen::seqN(3 * i, 3)) *= state.contact_mask(i);
                    z_ub_masked(Eigen::seqN(3 * i, 3)) *= state.contact_mask(i);
                }
                lb << opt_data.beq, bineq_lb, dv_lb, u_lb, z_lb_masked;
                ub << opt_data.beq, opt_data.bineq, dv_ub, u_ub, z_ub_masked;
                
                // Initialize Sparse Matrix:
                Eigen::SparseMatrix<double> sparse_H = opt_data.H.sparseView();
                Eigen::SparseMatrix<double> sparse_A = A.sparseView();
                sparse_H.makeCompressed();
                sparse_A.makeCompressed();

                // Initalize OSQP workspace:
                instance.objective_matrix = sparse_H;
                instance.objective_vector = opt_data.f;
                instance.constraint_matrix = sparse_A;
                instance.lower_bounds = lb;
                instance.upper_bounds = ub;
                
                // Check initialization:
                absl::Status result = solver.Init(instance, settings);
                return result;
            }        


        void update_mj_data(const State& state_in) {
            Vector<model::nq_size> qpos = Vector<model::nq_size>::Zero();
            Vector<model::nv_size> qvel = Vector<model::nv_size>::Zero();

            if constexpr (is_fixed_based) {
                qpos << state_in.motor_position;
                qvel << state_in.motor_velocity;
            } else {
                const Vector<3> zero_vector = {0.0, 0.0, 0.0};
                qpos << zero_vector, state_in.body_rotation, state_in.motor_position;
                qvel << state_in.linear_body_velocity, state_in.angular_body_velocity, state_in.motor_velocity;
            }

            // Update Mujoco Data members of the node class
            mj_data_->qpos = qpos.data();
            mj_data_->qvel = qvel.data();

            // Runs steps: 2-12, 12-18:
            mj_fwdPosition(mj_model_, mj_data_);
            mj_fwdVelocity(mj_model_, mj_data_);

            // Update Points:
            points_ = Eigen::Map<Matrix<model::site_ids_size, 3>>(mj_data_->site_xpos)(site_ids_, Eigen::placeholders::all);
        }



        void update_osc_data() {
            // Mass Matrix:
            Matrix<model::nv_size, model::nv_size> mass_matrix = 
                Matrix<model::nv_size, model::nv_size>::Zero();
            mj_fullM(mj_model_, mass_matrix.data(), mj_data_->qM);

            // Coriolis Matrix:
            Vector<model::nv_size> coriolis_matrix = 
                Eigen::Map<Vector<model::nv_size>>(mj_data_->qfrc_bias);

            // Generalized Positions and Velocities:
            Vector<model::nq_size> generalized_positions = 
                Eigen::Map<Vector<model::nq_size>>(mj_data_->qpos);
            Vector<model::nv_size> generalized_velocities = 
                Eigen::Map<Vector<model::nv_size>>(mj_data_->qvel);

            // Jacobian Calculation:
            Matrix<optimization::p_size, model::nv_size> jacobian_translation = 
                Matrix<optimization::p_size, model::nv_size>::Zero();
            Matrix<optimization::r_size, model::nv_size> jacobian_rotation = 
                Matrix<optimization::r_size, model::nv_size>::Zero();
            Matrix<optimization::p_size, model::nv_size> jacobian_dot_translation = 
                Matrix<optimization::p_size, model::nv_size>::Zero();
            Matrix<optimization::r_size, model::nv_size> jacobian_dot_rotation = 
                Matrix<optimization::r_size, model::nv_size>::Zero();
            
            for (int i = 0; i < model::body_ids_size; i++) {
                Matrix<3, model::nv_size> jacp = Matrix<3, model::nv_size>::Zero();
                Matrix<3, model::nv_size> jacr = Matrix<3, model::nv_size>::Zero();
                Matrix<3, model::nv_size> jacp_dot = Matrix<3, model::nv_size>::Zero();
                Matrix<3, model::nv_size> jacr_dot = Matrix<3, model::nv_size>::Zero();

                mj_jac(mj_model_, mj_data_, jacp.data(), jacr.data(), points_.row(i).data(), body_ids_[i]);
                mj_jacDot(mj_model_, mj_data_, jacp_dot.data(), jacr_dot.data(), points_.row(i).data(), body_ids_[i]);

                int row_offset = i * 3;
                for(int row_idx = 0; row_idx < 3; row_idx++) {
                    for(int col_idx = 0; col_idx < model::nv_size; col_idx++) {
                        jacobian_translation(row_idx + row_offset, col_idx) = jacp(row_idx, col_idx);
                        jacobian_rotation(row_idx + row_offset, col_idx) = jacr(row_idx, col_idx);
                        jacobian_dot_translation(row_idx + row_offset, col_idx) = jacp_dot(row_idx, col_idx);
                        jacobian_dot_rotation(row_idx + row_offset, col_idx) = jacr_dot(row_idx, col_idx);
                    }
                }
            }

            Matrix<optimization::s_size, model::nv_size> taskspace_jacobian = Matrix<optimization::s_size, model::nv_size>::Zero();
            Matrix<optimization::s_size, model::nv_size> jacobian_dot = Matrix<optimization::s_size, model::nv_size>::Zero();
            taskspace_jacobian << jacobian_translation, jacobian_rotation;
            jacobian_dot << jacobian_dot_translation, jacobian_dot_rotation;

            Vector<optimization::s_size> taskspace_bias = Vector<optimization::s_size>::Zero();
            taskspace_bias = jacobian_dot * generalized_velocities;

            Matrix<model::nv_size, optimization::z_size> contact_jacobian = 
                Matrix<model::nv_size, optimization::z_size>::Zero();

            contact_jacobian = jacobian_translation(
                Eigen::seq(Eigen::placeholders::end - Eigen::fix<optimization::z_size>, Eigen::placeholders::last),
                Eigen::placeholders::all
            ).transpose();

            // Assign to OSCData member variable:
            osc_data_.mass_matrix = mass_matrix;
            osc_data_.coriolis_matrix = coriolis_matrix;
            osc_data_.contact_jacobian = contact_jacobian;
            osc_data_.taskspace_jacobian = taskspace_jacobian;
            osc_data_.taskspace_bias = taskspace_bias;
            osc_data_.previous_q = generalized_positions;
            osc_data_.previous_qd = generalized_velocities;
        }        


        void update_optimization_data() {
            // Convert OSCData to Column Major for Casadi Functions:
            auto mass_matrix = matrix_utils::transformMatrix<double, model::nv_size, model::nv_size, matrix_utils::ColumnMajor>(osc_data_.mass_matrix.data());
            auto coriolis_matrix = matrix_utils::transformMatrix<double, model::nv_size, 1, matrix_utils::ColumnMajor>(osc_data_.coriolis_matrix.data());
            auto contact_jacobian = matrix_utils::transformMatrix<double, model::nv_size, optimization::z_size, matrix_utils::ColumnMajor>(osc_data_.contact_jacobian.data());
            auto taskspace_jacobian = matrix_utils::transformMatrix<double, optimization::s_size, model::nv_size, matrix_utils::ColumnMajor>(osc_data_.taskspace_jacobian.data());
            auto taskspace_bias = matrix_utils::transformMatrix<double, optimization::s_size, 1, matrix_utils::ColumnMajor>(osc_data_.taskspace_bias.data());
            
            // The taskspace_targets are now a member variable
            auto desired_taskspace_ddx = matrix_utils::transformMatrix<double, model::site_ids_size, 6, matrix_utils::ColumnMajor>(taskspace_targets_.data());
            
            // Evaluate Casadi Functions:
            auto Aeq_matrix = evaluate_function<AeqParams>(Aeq_ops, {design_vector_.data(), mass_matrix.data(), coriolis_matrix.data(), contact_jacobian.data()});
            auto beq_matrix = evaluate_function<beqParams>(beq_ops, {design_vector_.data(), mass_matrix.data(), coriolis_matrix.data(), contact_jacobian.data()});
            auto Aineq_matrix = evaluate_function<AineqParams>(Aineq_ops, {design_vector_.data()});
            auto bineq_matrix = evaluate_function<bineqParams>(bineq_ops, {design_vector_.data()});
            auto H_matrix = evaluate_function<HParams>(H_ops, {design_vector_.data(), desired_taskspace_ddx.data(), taskspace_jacobian.data(), taskspace_bias.data()});
            auto f_matrix = evaluate_function<fParams>(f_ops, {design_vector_.data(), desired_taskspace_ddx.data(), taskspace_jacobian.data(), taskspace_bias.data()});

            // Assign to OptimizationData member variable:
            opt_data_.H = H_matrix;
            opt_data_.f = f_matrix;
            opt_data_.Aeq = Aeq_matrix;
            opt_data_.beq = beq_matrix;
            opt_data_.Aineq = Aineq_matrix;
            opt_data_.bineq = bineq_matrix;
        }        

        
        absl::Status update_optimization() {
            // Concatenate Constraint Matrix:
            MatrixColMajor<optimization::constraint_matrix_rows, optimization::constraint_matrix_cols> A;
            A << opt_data_.Aeq, opt_data_.Aineq, Abox_;
            // Calculate Bounds:
            Vector<optimization::bounds_size> lb;
            Vector<optimization::bounds_size> ub;
            Vector<optimization::z_size> z_lb_masked = z_lb_;
            Vector<optimization::z_size> z_ub_masked = z_ub_;
            
            // Use a lock to safely read the contact mask from the shared state
            {
                std::lock_guard<std::mutex> lock(state_mutex_);
                for(int i = 0; i < model::contact_site_ids_size; i++) {
                    z_lb_masked(Eigen::seqN(3 * i, 3)) *= state_.contact_mask(i);
                    z_ub_masked(Eigen::seqN(3 * i, 3)) *= state_.contact_mask(i);
                }
            }
            
            lb << opt_data_.beq, bineq_lb_, dv_lb_, u_lb_, z_lb_masked;
            ub << opt_data_.beq, opt_data_.bineq, dv_ub_, u_ub_, z_ub_masked;
            
            // Initialize Sparse Matrix:
            Eigen::SparseMatrix<double> sparse_H = opt_data_.H.sparseView();
            Eigen::SparseMatrix<double> sparse_A = A.sparseView();
            sparse_H.makeCompressed();
            sparse_A.makeCompressed();

            // Check if sparsity changed:
            absl::Status result;
            auto sparsity_check = solver_.UpdateObjectiveAndConstraintMatrices(sparse_H, sparse_A);
            if(sparsity_check.ok()) {
                // Update Internal OSQP workspace:
                result.Update(solver_.SetObjectiveVector(opt_data_.f));
                result.Update(solver_.SetBounds(lb, ub));
            } else {
                // Reinitialize OSQP workspace:
                instance_.objective_matrix = sparse_H;
                instance_.objective_vector = opt_data_.f;
                instance_.constraint_matrix = sparse_A;
                instance_.lower_bounds = lb;
                instance_.upper_bounds = ub;
                
                result.Update(solver_.Init(instance_, settings_));
                
                // Set warmstart:
                result.Update(solver_.SetWarmStart(solution_, dual_solution_));
            }

            return result;
        }


        void solve_optimization() {
            // Solve the Optimization:
            exit_code_ = solver_.Solve();
            solution_ = solver_.primal_solution();
            dual_solution_ = solver_.dual_solution();
        }


        void reset_optimization() {
            // Set Warm Start to Zero:
            Vector<optimization::constraint_matrix_cols> primal_vector = Vector<optimization::constraint_matrix_cols>::Zero();
            Vector<optimization::constraint_matrix_rows> dual_vector = Vector<optimization::constraint_matrix_rows>::Zero();
            std::ignore = solver_.SetWarmStart(primal_vector, dual_vector);
        }




        // Inside your OSCNode.cpp file
        void state_callback(const osc_2_in::msg::OSCMujocoState::SharedPtr msg) {
            // Use a lock guard to ensure thread-safe access to the state variable.
            std::lock_guard<std::mutex> lock(state_mutex_);

            // Convert the data from the ROS 2 message to your internal State struct.
            state_.motor_position = Eigen::Map<Vector<model::nu_size>>(msg->motor_position.data());
            state_.motor_velocity = Eigen::Map<Vector<model::nu_size>>(msg->motor_velocity.data());
            state_.torque_estimate = Eigen::Map<Vector<model::nu_size>>(msg->torque_estimate.data());
            state_.body_rotation = Eigen::Map<Vector<4>>(msg->body_rotation.data());
            state_.linear_body_velocity = Eigen::Map<Vector<3>>(msg->linear_body_velocity.data());
            state_.angular_body_velocity = Eigen::Map<Vector<3>>(msg->angular_body_velocity.data());
            state_.contact_mask = Eigen::Map<Vector<model::contact_site_ids_size>>(msg->contact_mask.data());
        }

   

        void timer_callback() {
            // Lock the state mutex to get the latest robot state
            std::lock_guard<std::mutex> lock(state_mutex_);
            
            // 1. Locally compute task-space targets
            // The `taskspace_targets` variable is now local to this function.
            TaskspaceTargets taskspace_targets = Matrix<model::site_ids_size, 6>::Zero();
            // Insert your task-space target computation logic here.
            // For example, this could be a simple function call or a complex calculation.
            
            // 2. Update the controller with the latest state and the newly computed targets
            // These internal functions will use the `state_` and `taskspace_targets` variables.
            update_mj_data(state_);
            update_osc_data(taskspace_targets);
            update_optimization_data();

            // 3. Update and solve the optimization
            std::ignore = update_optimization();
            solve_optimization();

            // 4. Get the torque command from the solution and publish it
            publish_torque_command();
        }



        // New private function to handle publishing
        void publish_torque_command() {
            // `solution_` is a member variable, updated by `solve_optimization()`
            Vector<model::nu_size> torque_command = solution_(Eigen::seqN(optimization::dv_idx, optimization::u_size));

            auto torque_msg = std::make_unique<osc_2_in::msg::OSCTorqueCommand>();
            torque_msg->torque_command.assign(torque_command.data(), torque_command.data() + torque_command.size());
            torque_publisher_->publish(std::move(torque_msg));
        }        




        // ROS 2 members
        rclcpp::Subscription<osc_2_in::msg::OSCMujocoState>::SharedPtr state_subscriber_;
        rclcpp::Publisher<osc_2_in::msg::OSCTorqueCommand>::SharedPtr torque_publisher_;
        rclcpp::TimerBase::SharedPtr timer_;

        // Controller
        std::unique_ptr<OperationalSpaceController> controller_;
    };
