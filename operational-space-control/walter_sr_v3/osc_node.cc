#include "operational-space-control/walter_sr_v3/osc_node.h"

// Your anonymous namespace with Casadi functions goes here.
namespace {
    FunctionOperations Aeq_ops{
        .incref=Aeq_incref, .checkout=Aeq_checkout, .eval=Aeq, .release=Aeq_release, .decref=Aeq_decref};
    FunctionOperations beq_ops{
        .incref=beq_incref, .checkout=beq_checkout, .eval=beq, .release=beq_release, .decref=beq_decref};
    FunctionOperations Aineq_ops{
        .incref=Aineq_incref, .checkout=Aineq_checkout, .eval=Aineq, .release=Aineq_release, .decref=Aineq_decref};
    FunctionOperations bineq_ops{
        .incref=bineq_incref, .checkout=bineq_checkout, .eval=bineq, .release=bineq_release, .decref=bineq_decref};
    FunctionOperations H_ops{
        .incref=H_incref, .checkout=H_checkout, .eval=H, .release=H_release, .decref=H_decref};
    FunctionOperations f_ops{
        .incref=f_incref, .checkout=f_checkout, .eval=f, .release=f_release, .decref=f_decref};
    using AeqParams = FunctionParams<Aeq_SZ_ARG, Aeq_SZ_RES, Aeq_SZ_IW, Aeq_SZ_W, optimization::Aeq_rows, optimization::Aeq_cols, optimization::Aeq_sz, 4>;
    using beqParams = FunctionParams<beq_SZ_ARG, beq_SZ_RES, beq_SZ_IW, beq_SZ_W, optimization::beq_sz, 1, optimization::beq_sz, 4>;
    using AineqParams = FunctionParams<Aineq_SZ_ARG, Aineq_SZ_RES, Aineq_SZ_IW, Aineq_SZ_W, optimization::Aineq_rows, optimization::Aineq_cols, optimization::Aineq_sz, 2>;
    using bineqParams = FunctionParams<bineq_SZ_ARG, bineq_SZ_RES, bineq_SZ_IW, bineq_SZ_W, optimization::bineq_sz, 1, optimization::bineq_sz, 2>;
    using HParams = FunctionParams<H_SZ_ARG, H_SZ_RES, H_SZ_IW, H_SZ_W, optimization::H_rows, optimization::H_cols, optimization::H_sz, 4>;
    using fParams = FunctionParams<f_SZ_ARG, f_SZ_RES, f_SZ_IW, f_SZ_W, optimization::f_sz, 1, optimization::f_sz, 4>;

    // Helper function definitions
    template <typename T>
    bool contains(const std::vector<T>& vec, const T& value) {
        return std::find(vec.begin(), vec.end(), value) != vec.end();
    }
    
    std::vector<int> getSiteIdsOnSameBodyAsGeom(const mjModel* m, int geom_id) {
        std::vector<int> associated_site_ids;
        if (geom_id < 0 || geom_id >= m->ngeom) {
            std::cerr << "Error: Invalid geom ID: " << geom_id << std::endl;
            return associated_site_ids;
        }
        int geom_body_id = m->geom_bodyid[geom_id];
        for (int i = 0; i < m->nsite; ++i) {
            if (m->site_bodyid[i] == geom_body_id) {
                associated_site_ids.push_back(i);
            }
        }
        return associated_site_ids;
    }
    
    std::vector<int> getBinaryRepresentation_std_find(const std::vector<int>& A, const std::vector<int>& B) {
        std::vector<int> C;
        C.reserve(B.size());
        for (int b_element : B) {
            auto it = std::find(A.begin(), A.end(), b_element);
            C.push_back((it != A.end()) ? 1 : 0);
        }
        return C;
    }

    double get_propeller_leg_height(
            const Eigen::Quaterniond& body_quat, 
            double q_hip, double q_knee, 
            double q_hip_offset, double q_knee_offset,
            double L_thigh, double L_shin, double R_wheel) 
        {
            double theta_hip_calibrated  = q_hip + q_hip_offset;
            double theta_knee_calibrated = q_knee + q_knee_offset;

            Eigen::Vector3d v_thigh(0, 0, -L_thigh);
            double y_offset = -0.04675; 
            
            Eigen::Vector3d v_wheel_A( L_shin, y_offset, 0); 
            Eigen::Vector3d v_wheel_B(-L_shin, y_offset, 0);

            double global_shin = theta_hip_calibrated + theta_knee_calibrated;

            Eigen::AngleAxisd rot_thigh(theta_hip_calibrated, Eigen::Vector3d::UnitY());
            Eigen::AngleAxisd rot_shin(global_shin,           Eigen::Vector3d::UnitY());

            Eigen::Vector3d p_knee = rot_thigh * v_thigh;
            Eigen::Vector3d p_wheel_A = p_knee + (rot_shin * v_wheel_A);
            Eigen::Vector3d p_wheel_B = p_knee + (rot_shin * v_wheel_B);

            // Eigen::Vector3d g_A = body_quat * p_wheel_A;
            // Eigen::Vector3d g_B = body_quat * p_wheel_B;

            // return std::max(-g_A.z(), -g_B.z()) + R_wheel;
            return std::max(-p_wheel_A.z(), -p_wheel_B.z()) + R_wheel;            
        }

    double get_propeller_leg_height_velocity(
        double q_hip, double q_knee, 
        double dq_hip, double dq_knee,
        double q_hip_offset, double q_knee_offset,
        double L_thigh, double L_shin) 
    {
        double theta_hip  = q_hip + q_hip_offset;
        double theta_knee = q_knee + q_knee_offset;
        double global_shin = theta_hip + theta_knee;

        // Chain rule inner derivatives
        double d_theta_hip = dq_hip; 
        double d_global_shin = dq_hip + dq_knee;

        double sin_global_shin = std::sin(global_shin);
        double cos_global_shin = std::cos(global_shin);
        
        // Derivative of absolute value |x| is sign(x) * dx
        double sign_sin = (sin_global_shin >= 0.0) ? 1.0 : -1.0;

        // Analytical velocity calculation
        double v_height = -L_thigh * std::sin(theta_hip) * d_theta_hip 
                        + L_shin * sign_sin * cos_global_shin * d_global_shin;
        
        return v_height;
    }    

}

// Full constructor implementation
OSCNode::OSCNode(const std::string& xml_path)
    : Node("osc_node"),
      xml_path_(xml_path),
      solution_(Vector<optimization::design_vector_size>::Zero()),
      dual_solution_(Vector<optimization::constraint_matrix_rows>::Zero()),
      design_vector_(Vector<optimization::design_vector_size>::Zero()),
      infinity_(OSQP_INFTY),
      big_number_(1e4),
      Abox_(MatrixColMajor<optimization::design_vector_size, optimization::design_vector_size>::Identity()),
      dv_lb_(Vector<optimization::dv_size>::Constant(-infinity_)),
      dv_ub_(Vector<optimization::dv_size>::Constant(infinity_)),
      //   u_lb_({-20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0}),
      //   u_ub_({20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0}),
      // Hips get 12.0 N.m, Knees get 25.0 N.m
      u_lb_({
            -11.9, -23.8,
            -11.9, -23.8,
            -11.9, -23.8,
            -11.9, -23.8
        }),
      u_ub_({
            11.9, 23.8,
            11.9, 23.8,
            11.9, 23.8,
            11.9, 23.8
        }),       
      z_lb_({
          -infinity_, -infinity_, 0.0, -infinity_, -infinity_, 0.0, -infinity_, -infinity_, 0.0, -infinity_, -infinity_, 0.0,
          -infinity_, -infinity_, 0.0, -infinity_, -infinity_, 0.0, -infinity_, -infinity_, 0.0, -infinity_, -infinity_, 0.0}),
      z_ub_({
          infinity_, infinity_, big_number_, infinity_, infinity_, big_number_, infinity_, infinity_, big_number_, infinity_, infinity_, big_number_,
          infinity_, infinity_, big_number_, infinity_, infinity_, big_number_, infinity_, infinity_, big_number_, infinity_, infinity_, big_number_}),
      bineq_lb_(Vector<optimization::bineq_sz>::Constant(-infinity_))

{
    // --- Mujoco initialization ---
    char error[1000];
    mj_model_ = mj_loadXML(xml_path_.c_str(), nullptr, error, 1000);
    if (!mj_model_) {
        RCLCPP_FATAL(this->get_logger(), "Failed to load Mujoco Model: %s", error);
        throw std::runtime_error("Failed to load Mujoco Model.");
    }
    mj_model_->opt.timestep = 0.002;
    mj_data_ = mj_makeData(mj_model_);

    // check if keyframe is correct======================================================================
    mj_resetDataKeyframe(mj_model_, mj_data_, 11); // 
    mj_forward(mj_model_, mj_data_); // Compute initial kinematics
    
    
    // Thighs: 0, 2, 4, 6 in the 8-DOF motor array.
    // They correspond to indices 7, 9, 11, 13 in the full mj_data_->qpos array.
    // initial_tlh_angular_position_ = mj_data_->qpos[7 + 0]; // Index 7 (Motor 0)
    // initial_tlh_angular_position_ = mj_data_->qpos[7 + 0]; // Index 7 (Motor 0)
    // initial_trh_angular_position_ = mj_data_->qpos[7 + 2]; // Index 9 (Motor 2)
    // initial_hlh_angular_position_ = mj_data_->qpos[7 + 4]; // Index 11 (Motor 4)
    // initial_hrh_angular_position_ = mj_data_->qpos[7 + 6]; // Index 13 (Motor 6)

    // // Shins: 1, 3, 5, 7 in the 8-DOF motor array.
    // // They correspond to indices 8, 10, 12, 14 in the full mj_data_->qpos array.
    // initial_tl_angular_position_ = mj_data_->qpos[7 + 1]; // Index 8 (Motor 1)
    // initial_tr_angular_position_ = mj_data_->qpos[7 + 3]; // Index 10 (Motor 3)
    // initial_hl_angular_position_ = mj_data_->qpos[7 + 5]; // Index 12 (Motor 5)
    // initial_hr_angular_position_ = mj_data_->qpos[7 + 7]; // Index 14 (Motor 7)

    
    
    
    // Populate the site and body ID vectors
    for (const std::string_view& site : model::site_list) {
        std::string site_str = std::string(site);
        int id = mj_name2id(mj_model_, mjOBJ_SITE, site_str.data());
        assert(id != -1 && "Site not found in model.");
        sites_.push_back(site_str);
        site_ids_.push_back(id);
    }
    for (const std::string_view& site : model::noncontact_site_list) {
        std::string site_str = std::string(site);
        int id = mj_name2id(mj_model_, mjOBJ_SITE, site_str.data());
        assert(id != -1 && "Site not found in model.");
        noncontact_sites_.push_back(site_str);
        noncontact_site_ids_.push_back(id);
    }
    for (const std::string_view& site : model::contact_site_list) {
        std::string site_str = std::string(site);
        int id = mj_name2id(mj_model_, mjOBJ_SITE, site_str.data());
        assert(id != -1 && "Site not found in model.");
        contact_sites_.push_back(site_str);
        contact_site_ids_.push_back(id);
    }
    for (const std::string_view& body : model::body_list) {
        std::string body_str = std::string(body);
        int id = mj_name2id(mj_model_, mjOBJ_BODY, body_str.data());
        assert(id != -1 && "Body not found in model.");
        bodies_.push_back(body_str);
        body_ids_.push_back(id);
    }

    // SAFETY CHECK: Verify Site-Body Alignment
    for (size_t i = 0; i < site_ids_.size(); i++) {
        int site_id = site_ids_[i];
        int body_id_from_list = body_ids_[i]; // The one from your generated list
        int body_id_real = mj_model_->site_bodyid[site_id]; // The truth from MuJoCo

        if (body_id_from_list != body_id_real) {
            RCLCPP_FATAL(this->get_logger(), 
                "MISMATCH at index %zu! Site '%s' is on Body %d, but List says Body %d", 
                i, sites_[i].c_str(), body_id_real, body_id_from_list);
            throw std::runtime_error("Kinematic Chain Mismatch");
        }
    }
    
    assert(site_ids_.size() == body_ids_.size() && "Number of Sites and Bodies must be equal.");

    // --- Optimization Initialization ---
    // Create an initial state message to use for setup.
    // OSCMujocoState initial_state_msg;
    // initial_state_msg.motor_position.assign(model::nu_size, 0.0);
    // initial_state_msg.motor_velocity.assign(model::nu_size, 0.0);
    // initial_state_msg.torque_estimate.assign(model::nu_size, 0.0);
    // initial_state_msg.body_rotation.assign(4, 0.0);
    // initial_state_msg.linear_body_velocity.assign(3, 0.0);
    // initial_state_msg.angular_body_velocity.assign(3, 0.0);
    // initial_state_msg.contact_mask.assign(model::contact_site_ids_size, 0.0);

    // std::fill(initial_state_msg.motor_position.begin(), initial_state_msg.motor_position.end(), 0.0f);
    // std::fill(initial_state_msg.motor_velocity.begin(), initial_state_msg.motor_velocity.end(), 0.0f);
    // std::fill(initial_state_msg.torque_estimate.begin(), initial_state_msg.torque_estimate.end(), 0.0f);
    // std::fill(initial_state_msg.body_rotation.begin(), initial_state_msg.body_rotation.end(), 0.0f);
    // std::fill(initial_state_msg.linear_body_velocity.begin(), initial_state_msg.linear_body_velocity.end(), 0.0f);
    // std::fill(initial_state_msg.angular_body_velocity.begin(), initial_state_msg.angular_body_velocity.end(), 0.0f);
    // std::fill(initial_state_msg.contact_mask.begin(), initial_state_msg.contact_mask.end(), false);
    
    
    // state_callback(std::make_shared<OSCMujocoState>(initial_state_msg));
    
    // update_mj_data();

    // --- Optimization Initialization ---
    // Instead of using a dummy ROS message to set state_ to zero, 
    // populate state_ with the actual initial Keyframe 5 data from mj_data_.
    
    // 1. Populate state_.motor_position from mj_data_->qpos
    //    Motor positions start at index 7 in the floating-base qpos array (3-pos + 4-quat).
    
    // Assuming model::nu_size is 8:
    // for (size_t i = 0; i < model::nu_size; ++i) {
    //     // qpos index = 7 (base pos/quat end) + i (motor index)
    //     state_.motor_position(i) = mj_data_->qpos[7 + i];
    //     // Ensure other essential fields are also non-zero if needed, 
    //     // e.g., base rotation:
    //     if (i < 4) {
    //         state_.body_rotation(i) = mj_data_->qpos[3 + i];
    //     }
    // }
    // // You can clear velocities and torques as they should start at zero.
    // state_.motor_velocity.setZero();
    // state_.linear_body_velocity.setZero();
    // state_.angular_body_velocity.setZero();
    // state_.torque_estimate.setZero();



    Vector<model::nq_size> qpos = Eigen::Map<Vector<model::nq_size>>(mj_data_->qpos);
    // initial_position_ = qpos(Eigen::seqN(0, 3));    
    
    Vector<model::contact_site_ids_size> initial_contact_mask = Vector<model::contact_site_ids_size>::Constant(1.0);

    absl::Status result = set_up_optimization(initial_contact_mask);
    if (!result.ok()) {
        RCLCPP_FATAL(this->get_logger(), "Failed to initialize optimization: %s", result.message().data());
        throw std::runtime_error("Failed to initialize optimization.");
    }    
    
    // --- ROS 2 communication setup ---
    state_subscriber_ = this->create_subscription<OSCMujocoState>(
        "/state_estimator/state", 1, std::bind(&OSCNode::state_callback, this, std::placeholders::_1));
    // taskspace_targets_subscriber_ = this->create_subscription<OSCTaskspaceTargets>(
    //     "osc/taskspace_targets", 10, std::bind(&OSCNode::taskspace_targets_callback, this, std::placeholders::_1));
    torque_publisher_ = this->create_publisher<Command>("walter/command", 1);
    // torque_publisher_ = this->create_publisher<OSCTorqueCommand>("walter/command", 10);
    // New: 5000 microseconds (5 ms = 200 Hz)
    
    //===========================================================================
    auto data_qos = rclcpp::QoS(rclcpp::KeepLast(10));
    data_qos.reliability(RMW_QOS_POLICY_RELIABILITY_BEST_EFFORT);
    data_qos.durability(RMW_QOS_POLICY_DURABILITY_VOLATILE);

    data_pub_ = this->create_publisher<std_msgs::msg::Float64MultiArray>(
        "walter/data", data_qos);

    // Initialize the telemetry message layout ONCE
    // int num_sites = model::site_ids_size;
    // int num_dof = 6; 
    
    data_msg_.layout.dim.resize(4);
    

    // Add your metadata note here!

    data_msg_.layout.dim[0].label = "target_hip_z, hip_z_tl, hip_z_tr, hip_z_hl, hip_z_hr, target_hip_z_vel, hip_zv_tl, hip_zv_tr, hip_zv_hl, hip_zv_hr, shin_pos_tl_target, shin_pos_tr_target, shin_pos_hl_target, shin_pos_hr_target, shin_pos_tl, shin_pos_tr, shin_pos_hl, shin_pos_hr, shin_vel_target, shin_vel_tl, shin_vel_tr, shin_vel_hl, shin_vel_hr, body_x, body_y, body_z, contact_tlf, contact_tlr, contact_trf, contact_trr, contact_hlf, contact_hlr, contact_hrf, contact_hrr, osqp_exit, time_casadi, time_osqp, qp_obj, tau_rlh, tau_rlk, tau_rrh, tau_rrk, tau_flh, tau_flk, tau_frh, tau_frk, fz_tlf, fz_tlr, fz_trf, fz_trr, fz_hlf, fz_hlr, fz_hrf, fz_hrr"; 
    
    data_msg_.data.reserve(54); // Exactly 54 elements now    
    // // Reserve memory so push_back is zero-overhead
    // // data_msg_.data.reserve(num_sites * num_dof);        
    // data_msg_.data.reserve(1);        
    //===========================================================================

    timer_ = this->create_wall_timer(std::chrono::microseconds(5000), std::bind(&OSCNode::timer_callback, this));

    rclcpp::on_shutdown([this]() {
        RCLCPP_WARN(this->get_logger(), "Shutdown signal received. Attempting to stop robot...");
        this->stop_robot();
    });

}

OSCNode::~OSCNode() {
    mj_deleteData(mj_data_);
    mj_deleteModel(mj_model_);
}

// ===============================================================================================================
// Full implementation of all methods
void OSCNode::state_callback(const OSCMujocoState::SharedPtr msg) {
    std::lock_guard<std::mutex> lock(state_mutex_);
    // Manually copy and cast each member to the correct double type
    for (size_t i = 0; i < model::nu_size; ++i) {
        state_.motor_position(i) = static_cast<double>(msg->motor_position[i]);
        state_.motor_velocity(i) = static_cast<double>(msg->motor_velocity[i]);
        state_.torque_estimate(i) = static_cast<double>(msg->torque_estimate[i]);

        // CAPTURE DETECTED POSITION
        last_detected_motor_position_(i) = state_.motor_position(i);        
        
    }

    // CAPTURE STATE READ TIME
    state_read_time_ = std::chrono::high_resolution_clock::now();    
    
    for (size_t i = 0; i < 4; ++i) {
        state_.body_rotation(i) = static_cast<double>(msg->body_rotation[i]);
    }

    for (size_t i = 0; i < 3; ++i) {
        state_.linear_body_velocity(i) = static_cast<double>(msg->linear_body_velocity[i]);
        state_.angular_body_velocity(i) = static_cast<double>(msg->angular_body_velocity[i]);
    }

    for (size_t i = 0; i < model::contact_site_ids_size; ++i) {
        state_.contact_mask(i) = static_cast<double>(msg->contact_mask[i]);
    }
    is_state_received_ = true;    
}



// ===============================================================================================================
void OSCNode::timer_callback() {

    // Initial variables
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
        if (!is_state_received_) {
            RCLCPP_WARN_THROTTLE(this->get_logger(), *this->get_clock(), 1000, "Waiting for initial state message...");
            return; 
        }
    }
    time_wait_for_execution_ms_ = std::chrono::duration<double, std::milli>(t_start_execution - local_state_read_time).count();    
    double dt = current_time - last_time_;
    if (dt == 0.0) dt = 0.0001;
    
    // initial hip height
    static double hip_z_tl_initial = 0.14;
    static double hip_z_tr_initial = 0.14;
    static double hip_z_hl_initial = 0.14;
    static double hip_z_hr_initial = 0.14;

    // initial shin position
    static double shin_pos_tl_initial = 0.0;
    static double shin_pos_tr_initial = 0.0;
    static double shin_pos_hl_initial = 0.0;
    static double shin_pos_hr_initial = 0.0;    

    // initial contact switch info
    static double gait_start_time = 0.0;
    static Eigen::Vector<double, model::contact_site_ids_size> contact_start_times = Eigen::Vector<double, model::contact_site_ids_size>::Constant(-100.0);
    static Eigen::Vector<double, model::contact_site_ids_size> prev_contact_mask = Eigen::Vector<double, model::contact_site_ids_size>::Zero();
    const double soft_switch_max_force = 770.0; 
    const double soft_switch_ramp_time = 0.5;    
    

    Eigen::Quaterniond q_body(local_state.body_rotation(0), local_state.body_rotation(1), 
                            local_state.body_rotation(2), local_state.body_rotation(3));

    const double L_THIGH = 0.1016;  // hip to knee
    const double L_SHIN  = 0.08255; // knee to wheel center
    const double R_WHEEL = 0.0635;  // Radius of the wheel
    double hip_z_tl = get_propeller_leg_height(q_body, local_state.motor_position(0), local_state.motor_position(1), 0, 0, L_THIGH, L_SHIN, R_WHEEL);
    double hip_z_tr = get_propeller_leg_height(q_body, local_state.motor_position(2), local_state.motor_position(3), 0, 0, L_THIGH, L_SHIN, R_WHEEL);
    double hip_z_hl = get_propeller_leg_height(q_body, local_state.motor_position(4), local_state.motor_position(5), 0, 0, L_THIGH, L_SHIN, R_WHEEL);
    double hip_z_hr = get_propeller_leg_height(q_body, local_state.motor_position(6), local_state.motor_position(7), 0, 0, L_THIGH, L_SHIN, R_WHEEL);


    // Zero time step
    if (last_time_ == 0.0) {
        update_mj_data(local_state); // Pass local state

        // Set initial variables
        hip_z_tl_initial = hip_z_tl;
        hip_z_tr_initial = hip_z_tr;
        hip_z_hl_initial = hip_z_hl;
        hip_z_hr_initial = hip_z_hr;
        gait_start_time = current_time;
        shin_pos_tl_initial  =  local_state.motor_position(1);
        shin_pos_tr_initial  =  local_state.motor_position(3);
        shin_pos_hl_initial  =  local_state.motor_position(5);
        shin_pos_hr_initial  =  local_state.motor_position(7);

        
        // Build the OSQP matrix
        update_osc_data();

        // Initial contact force limits
        Eigen::Vector<double, model::contact_site_ids_size> initial_force_limits;
        for (int i=0; i < model::contact_site_ids_size; ++i) {
            // If starting on the ground, give it the full 770N limit instantly
            initial_force_limits(i) = local_state.contact_mask(i) > 0.5 ? soft_switch_max_force : 0.0;
        }
        prev_contact_mask = local_state.contact_mask;

        update_optimization_data(initial_force_limits);
        absl::Status result = set_up_optimization(initial_force_limits);

        if (!result.ok()) {
            RCLCPP_FATAL(this->get_logger(), "Failed to initialize optimization: %s", result.message().data());
            std::lock_guard<std::mutex> lock_state(state_mutex_);
            safety_override_active_ = true;
        } else {
            RCLCPP_INFO(this->get_logger(), "OSQP Optimizer successfully initialized with first robot state.");
        }

        last_time_ = current_time;
        return; 
    }

    
    // Hip height and hip angle safety checks
    bool limit_hit = local_safety_override_active;     
    if (!local_safety_override_active) {
        const double THIGH_LIMIT = 1.95;
        const double MIN_HIP_HEIGHT = 0.11; // 0.13 meters
        // 1. Check Hip Heights
        if (hip_z_tl < MIN_HIP_HEIGHT || hip_z_tr < MIN_HIP_HEIGHT || 
            hip_z_hl < MIN_HIP_HEIGHT || hip_z_hr < MIN_HIP_HEIGHT) {
            limit_hit = true;
            RCLCPP_WARN_ONCE(this->get_logger(), "Hip height limit (< %.2fm) hit. Overriding control.", MIN_HIP_HEIGHT);
        }
        
        // 2. Check Hip (Thigh) Joints (0, 2, 4, 6) using local_state.motor_position
        if (!limit_hit) {
            // Check Thighs (0, 2, 4, 6) using local_state.motor_position
            for (size_t i : {0, 2, 4, 6}) {
                if (std::abs(local_state.motor_position(i)) >= THIGH_LIMIT) {
                    limit_hit = true;
                    RCLCPP_WARN_ONCE(this->get_logger(), "Absolute THIGH limit (%.2f rad) hit on motor index %zu. Overriding control.", THIGH_LIMIT, i);
                    break; 
                }
            }
        }

        // If a new limit was hit, set the SHARED safety flag to true (permanently)
        if (limit_hit) { 
            std::lock_guard<std::mutex> lock_state(state_mutex_);
            safety_override_active_ = true;
            local_safety_override_active = true; // Update local for subsequent steps
        }
    }


    if (!local_safety_override_active) {


        // ===============================================================
        // --- HARDCODED TIME-BASED CONTACT MASK AND TUMBLE TARGETS ---
        // ===============================================================
        double elapsed_t = current_time - gait_start_time;
        
        // 1. Constant Tumbling Velocity (Matches Sim: 1.5 rad/s)
        double shin_rot_vel = 1.5; 
        double shin_vel_target = shin_rot_vel;

        // 2. Open-Loop Position Targets
        double shin_pos_tl_target = shin_pos_tl_initial + shin_rot_vel * elapsed_t;
        double shin_pos_tr_target = shin_pos_tr_initial + shin_rot_vel * elapsed_t;
        double shin_pos_hl_target = shin_pos_hl_initial + shin_rot_vel * elapsed_t;
        double shin_pos_hr_target = shin_pos_hr_initial + shin_rot_vel * elapsed_t;

        // 3. Timing Parameters (Calibrated from Simulation)
        const double x_front_period = 1.42; 
        const double y_back_period  = 1.42;
        const double z_front_start  = 1.30;
        const double w_back_start   = 0.638;
        
        const bool initial_front_F = true;
        const bool initial_back_F  = true;

        // 4. Calculate Front/Head Flips
        bool front_F = initial_front_F;
        if (elapsed_t >= z_front_start) {
            int flips = static_cast<int>((elapsed_t - z_front_start) / x_front_period) + 1;
            if (flips % 2 != 0) front_F = !initial_front_F;
        }

        // 5. Calculate Back/Torso Flips
        bool back_F = initial_back_F;
        if (elapsed_t >= w_back_start) {
            int flips = static_cast<int>((elapsed_t - w_back_start) / y_back_period) + 1;
            if (flips % 2 != 0) back_F = !initial_back_F;
        }

        // 6. Apply to Contact Mask
        // Back/Torso: indices 0 (F), 1 (R), 2 (F), 3 (R)
        local_state.contact_mask(0) = back_F ? 1.0 : 0.0;
        local_state.contact_mask(1) = back_F ? 0.0 : 1.0;
        local_state.contact_mask(2) = back_F ? 1.0 : 0.0;
        local_state.contact_mask(3) = back_F ? 0.0 : 1.0;

        // Front/Head: indices 4 (F), 5 (R), 6 (F), 7 (R)
        local_state.contact_mask(4) = front_F ? 1.0 : 0.0;
        local_state.contact_mask(5) = front_F ? 0.0 : 1.0;
        local_state.contact_mask(6) = front_F ? 1.0 : 0.0;
        local_state.contact_mask(7) = front_F ? 0.0 : 1.0;



        // Calculate contact force limit
        Eigen::Vector<double, model::contact_site_ids_size> current_force_limits;
        for(int i=0; i < model::contact_site_ids_size; ++i) {
            bool is_contact = (local_state.contact_mask(i) > 0.5);
            bool was_contact = (prev_contact_mask[i] > 0.5);

            if (is_contact && !was_contact) contact_start_times[i] = current_time;

            double limit = 0.0;
            if (is_contact) {
                double duration = current_time - contact_start_times[i];
                double ratio = std::clamp(duration / soft_switch_ramp_time, 0.0, 1.0);
                limit = ratio * soft_switch_max_force;
            }
            current_force_limits[i] = limit;
        }
        prev_contact_mask = local_state.contact_mask;


        // 1. Update Mujoco Data for Kinematics 
        auto t_start_kinematics = std::chrono::high_resolution_clock::now();        
        update_mj_data(local_state); 
        
        // ===============================================================
        // Z control targets
        // ===============================================================
        // thigh - (kp - 600.0 — kd - 45.0)
        // double thigh_z_kp = 1300.0; double thigh_z_kv = 72.0;
        double thigh_z_kp = 1150.0; double thigh_z_kv = 68.0;
        // double thigh_z_kp = 2200.0; double thigh_z_kv = 90.0;

        // Use instantaneous motor velocities to calculate exact Z velocity
        double hip_zv_tl = get_propeller_leg_height_velocity(
            local_state.motor_position(0), local_state.motor_position(1),
            local_state.motor_velocity(0), local_state.motor_velocity(1),
            0, 0, L_THIGH, L_SHIN);

        double hip_zv_tr = get_propeller_leg_height_velocity(
            local_state.motor_position(2), local_state.motor_position(3),
            local_state.motor_velocity(2), local_state.motor_velocity(3),
            0, 0, L_THIGH, L_SHIN);

        double hip_zv_hl = get_propeller_leg_height_velocity(
            local_state.motor_position(4), local_state.motor_position(5),
            local_state.motor_velocity(4), local_state.motor_velocity(5),
            0, 0, L_THIGH, L_SHIN);

        double hip_zv_hr = get_propeller_leg_height_velocity(
            local_state.motor_position(6), local_state.motor_position(7),
            local_state.motor_velocity(6), local_state.motor_velocity(7),
            0, 0, L_THIGH, L_SHIN);


        // ===============================================================
        // HIP HEIGHT TARGET
        // ===============================================================
        const double BASE_HEIGHT = hip_z_tl_initial; // Anchor around the starting height
        double target_hip_z = BASE_HEIGHT;
        double target_hip_z_vel = 0.0;
        // ===============================================================        


        // ===============================================================
        // SHIN VELOCITY AND POSITION TARGET
        // ===============================================================        
        double shin_pos_tl  =  local_state.motor_position(1);
        double shin_pos_tr  =  local_state.motor_position(3);
        double shin_pos_hl  =  local_state.motor_position(5);
        double shin_pos_hr  =  local_state.motor_position(7);

        double shin_vel_tl  =  local_state.motor_velocity(1);
        double shin_vel_tr  =  local_state.motor_velocity(3);
        double shin_vel_hl  =  local_state.motor_velocity(5);
        double shin_vel_hr  =  local_state.motor_velocity(7);

        double shin_kp = 100.0; 
        double shin_kv = 5.0;
        // ===============================================================        



        // Shin DDQ Commands (using local_state)
        double tl_ddq_cmd  = shin_kp * (shin_pos_tl_target - shin_pos_tl) + shin_kv * (shin_vel_target - shin_vel_tl);
        double tr_ddq_cmd  = shin_kp * (shin_pos_tr_target - shin_pos_tr) + shin_kv * (shin_vel_target - shin_vel_tr);
        double hl_ddq_cmd  = shin_kp * (shin_pos_hl_target - shin_pos_hl) + shin_kv * (shin_vel_target - shin_vel_hl);
        double hr_ddq_cmd  = shin_kp * (shin_pos_hr_target - shin_pos_hr) + shin_kv * (shin_vel_target - shin_vel_hr);

        // Thigh DDQ Commands (using local_state) -- Hip height
        double tl_hip_z_ddq_cmd = thigh_z_kp * (target_hip_z - hip_z_tl) + thigh_z_kv * (target_hip_z_vel - hip_zv_tl);
        double tr_hip_z_ddq_cmd = thigh_z_kp * (target_hip_z - hip_z_tr) + thigh_z_kv * (target_hip_z_vel - hip_zv_tr);
        double hl_hip_z_ddq_cmd = thigh_z_kp * (target_hip_z - hip_z_hl) + thigh_z_kv * (target_hip_z_vel - hip_zv_hl);
        double hr_hip_z_ddq_cmd = thigh_z_kp * (target_hip_z - hip_z_hr) + thigh_z_kv * (target_hip_z_vel - hip_zv_hr);

        // Populate Taskspace Targets Matrix 
        taskspace_targets_.setZero(); 
        taskspace_targets_.row(1)(4) = tl_ddq_cmd; taskspace_targets_.row(2)(4) = tr_ddq_cmd;
        taskspace_targets_.row(3)(4) = hl_ddq_cmd; taskspace_targets_.row(4)(4) = hr_ddq_cmd;
        taskspace_targets_.row(5)(2) = tl_hip_z_ddq_cmd; taskspace_targets_.row(6)(2) = tr_hip_z_ddq_cmd;
        taskspace_targets_.row(7)(2) = hl_hip_z_ddq_cmd; taskspace_targets_.row(8)(2) = hr_hip_z_ddq_cmd;



        data_msg_.data.clear();

        data_msg_.data.push_back(target_hip_z);

        data_msg_.data.push_back(hip_z_tl);
        data_msg_.data.push_back(hip_z_tr);
        data_msg_.data.push_back(hip_z_hl);
        data_msg_.data.push_back(hip_z_hr);

        data_msg_.data.push_back(target_hip_z_vel);

        data_msg_.data.push_back(hip_zv_tl);
        data_msg_.data.push_back(hip_zv_tr);
        data_msg_.data.push_back(hip_zv_hl);
        data_msg_.data.push_back(hip_zv_hr);


        data_msg_.data.push_back(shin_pos_tl_target);
        data_msg_.data.push_back(shin_pos_tr_target);
        data_msg_.data.push_back(shin_pos_hl_target);
        data_msg_.data.push_back(shin_pos_hr_target);

        data_msg_.data.push_back(shin_pos_tl);
        data_msg_.data.push_back(shin_pos_tr);
        data_msg_.data.push_back(shin_pos_hl);
        data_msg_.data.push_back(shin_pos_hr);

        data_msg_.data.push_back(shin_vel_target);

        data_msg_.data.push_back(shin_vel_tl);
        data_msg_.data.push_back(shin_vel_tr);
        data_msg_.data.push_back(shin_vel_hl);
        data_msg_.data.push_back(shin_vel_hr);

        data_msg_.data.push_back(mj_data_->qpos[0]); // body_x
        data_msg_.data.push_back(mj_data_->qpos[1]); // body_y
        data_msg_.data.push_back(mj_data_->qpos[2]); // body_z
        
        // --- RECORD CONTACT MASK ---
        // (Indices: 0=TLF, 1=TLR, 2=TRF, 3=TRR, 4=HLF, 5=HLR, 6=HRF, 7=HRR)
        for (int i = 0; i < 8; ++i) {
            data_msg_.data.push_back(local_state.contact_mask(i));
        }        
        // 3. Publish
        // --- SOLVER HEALTH ---
        data_msg_.data.push_back(static_cast<double>(exit_code_));
        data_msg_.data.push_back(time_casadi_update_ms_);
        data_msg_.data.push_back(time_osqp_solve_ms_);
        
        // Calculate QP Objective (0.5 * x'Hx + f'x)
        double quadratic_term = 0.5 * solution_.transpose() * opt_data_.H * solution_;
        double linear_term = opt_data_.f.transpose() * solution_;
        data_msg_.data.push_back(quadratic_term + linear_term);

        // --- TORQUES ---
        Vector<model::nu_size> osc_torque = solution_(Eigen::seqN(optimization::dv_idx, optimization::u_size));
        for (int i = 0; i < model::nu_size; ++i) {
            data_msg_.data.push_back(osc_torque(i));
        }

        // --- PREDICTED CONTACT FORCES (Z-Axis Only) ---
        // Z forces are every 3rd element starting at u_idx + 2
        for (int i = 0; i < model::contact_site_ids_size; ++i) {
            data_msg_.data.push_back(solution_(optimization::z_idx + (3 * i) + 2));
        }
        
        // 3. Publish
        data_pub_->publish(data_msg_);
        // ==============================================================================


        
        // Solve Optimization
        update_osc_data();
        
        // --- TIMING POINT B: END MUJOCO/KINEMATICS, START CASADI/OSQP DATA ---        
        auto t_start_casadi = std::chrono::high_resolution_clock::now();        
        
        update_optimization_data(current_force_limits);
        std::ignore = update_optimization(current_force_limits); 

        // --- TIMING POINT C: END CASADI/OSQP DATA, START SOLVE ---
        auto t_start_solve = std::chrono::high_resolution_clock::now();


        // OLD: solve_optimization();
        // NEW: Check result
        bool solver_success = solve_optimization();
        
        if (!solver_success) {
            // CRITICAL: Trigger safety override immediately if math fails
            std::lock_guard<std::mutex> lock(state_mutex_);
            safety_override_active_ = true;
            local_safety_override_active = true; 
            RCLCPP_ERROR(this->get_logger(), "Optimization failed! Engaging safety override.");
        }
                
        // --- TIMING POINT D: END SOLVE ---
        auto t_end_solve = std::chrono::high_resolution_clock::now();
        
        // Calculate and store internal execution times 
        time_mujoco_update_ms_ = std::chrono::duration<double, std::milli>(t_start_casadi - t_start_kinematics).count();
        time_casadi_update_ms_ = std::chrono::duration<double, std::milli>(t_start_solve - t_start_casadi).count();
        time_osqp_solve_ms_ = std::chrono::duration<double, std::milli>(t_end_solve - t_start_solve).count();                
    }

    if (!local_safety_override_active) {
        std::lock_guard<std::mutex> lock(state_mutex_); // Quick lock
        if (safety_override_active_) {
            // Aha! The flag changed from false to true while we were calculating.
            // Update our local copy so we send the safety command instead.
            local_safety_override_active = true; 
        }
    }    

    last_time_ = current_time;
    // Publish using the determined safety status and the captured timestamp
    publish_torque_command(local_safety_override_active, local_state_read_time); 
}
// ---------------------------------------------------------------------------------------------------------




void OSCNode::update_mj_data(const State& current_state) {
    // -------------------------------------------------------------
    // 1. DIRECT KINEMATICS (Orientation & Joints)
    // -------------------------------------------------------------
    mj_data_->qpos[3] = current_state.body_rotation(0); // w
    mj_data_->qpos[4] = current_state.body_rotation(1); // x
    mj_data_->qpos[5] = current_state.body_rotation(2); // y
    mj_data_->qpos[6] = current_state.body_rotation(3); // z

    if (current_state.body_rotation.norm() < 1e-6) {
        mj_data_->qpos[3] = 1.0;
    }

    for (int i = 0; i < model::nu_size; ++i) {
        mj_data_->qpos[7 + i] = current_state.motor_position(i);
        mj_data_->qvel[6 + i] = current_state.motor_velocity(i);
    }

    // -------------------------------------------------------------
    // 2. TRANSLATION (Invariant XY + Filtered Smooth Z)
    // -------------------------------------------------------------
    // Keep X and Y centered at 0. Dynamics and torques are invariant to XY translation.
    mj_data_->qpos[0] = 0.0;
    mj_data_->qpos[1] = 0.0;

    // Probe FK to get instantaneous foot heights relative to base
    mj_fwdPosition(mj_model_, mj_data_);

    double lowest_contact_z = 100.0;
    for (int id : contact_site_ids_) {
        double z = mj_data_->site_xpos[3 * id + 2];
        if (z < lowest_contact_z) {
            lowest_contact_z = z;
        }
    }

    // Kinematic ground-referenced height
    double raw_base_z = -lowest_contact_z + 0.0635;

    // Low-pass / Complementary filter on Z to prevent discontinuous jumps on contact changes
    static double filtered_base_z = 0.20;
    const double alpha_z = 0.05; // 5% correction per tick at 200 Hz (smooth ~0.1s time constant)
    filtered_base_z = (1.0 - alpha_z) * filtered_base_z + alpha_z * raw_base_z;
    
    mj_data_->qpos[2] = std::clamp(filtered_base_z, 0.10, 0.50);

    // -------------------------------------------------------------
    // 3. VELOCITIES (Fused IMU + Kinematics)
    // -------------------------------------------------------------
    // Angular rates directly from IMU Gyro
    mj_data_->qvel[3] = current_state.angular_body_velocity(0);
    mj_data_->qvel[4] = current_state.angular_body_velocity(1);
    mj_data_->qvel[5] = current_state.angular_body_velocity(2);

    // Linear velocity from EKF / estimator (or fused contact kinematics)
    mj_data_->qvel[0] = current_state.linear_body_velocity(0);
    mj_data_->qvel[1] = current_state.linear_body_velocity(1);
    mj_data_->qvel[2] = current_state.linear_body_velocity(2);

    // -------------------------------------------------------------
    // 4. FULL EVALUATION FOR CASADI MATRICES
    // -------------------------------------------------------------
    mj_fwdPosition(mj_model_, mj_data_);
    mj_fwdVelocity(mj_model_, mj_data_);

    // Map site positions for operational space tasks
    points_ = Eigen::Map<Matrix<model::site_ids_size, 3>>(
        mj_data_->site_xpos)(site_ids_, Eigen::placeholders::all);
}



// ===============================================================================================================
void OSCNode::update_osc_data() {
    Matrix<model::nv_size, model::nv_size> mass_matrix = Matrix<model::nv_size, model::nv_size>::Zero();
    mj_fullM(mj_model_, mass_matrix.data(), mj_data_->qM);
    Vector<model::nv_size> coriolis_matrix = Eigen::Map<Vector<model::nv_size>>(mj_data_->qfrc_bias);

    // ===============================
    // --- NEW: Add Passive Force Compensation ---
    // Mapping qfrc_passive (Damping/Springs)
    Vector<model::nv_size> passive_forces = Eigen::Map<Vector<model::nv_size>>(mj_data_->qfrc_passive);
    
    // Total non-actuated forces = Bias - Passive
    // This ensures the solver accounts for the 'help' or 'resistance' from joint damping
    coriolis_matrix -= passive_forces;
    // ===============================

    Vector<model::nq_size> generalized_positions = Eigen::Map<Vector<model::nq_size>>(mj_data_->qpos);
    Vector<model::nv_size> generalized_velocities = Eigen::Map<Vector<model::nv_size>>(mj_data_->qvel);

    Matrix<optimization::p_size, model::nv_size> jacobian_translation = Matrix<optimization::p_size, model::nv_size>::Zero();
    Matrix<optimization::r_size, model::nv_size> jacobian_rotation = Matrix<optimization::r_size, model::nv_size>::Zero();
    Matrix<optimization::p_size, model::nv_size> jacobian_dot_translation = Matrix<optimization::p_size, model::nv_size>::Zero();
    Matrix<optimization::r_size, model::nv_size> jacobian_dot_rotation = Matrix<optimization::r_size, model::nv_size>::Zero();
    
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
    Matrix<model::nv_size, optimization::z_size> contact_jacobian = Matrix<model::nv_size, optimization::z_size>::Zero();
    contact_jacobian = jacobian_translation(Eigen::seq(Eigen::placeholders::end - Eigen::fix<optimization::z_size>, Eigen::placeholders::last), Eigen::placeholders::all).transpose();

    osc_data_.mass_matrix = mass_matrix;
    osc_data_.coriolis_matrix = coriolis_matrix;
    osc_data_.contact_jacobian = contact_jacobian;
    osc_data_.taskspace_jacobian = taskspace_jacobian;
    osc_data_.taskspace_bias = taskspace_bias;
    osc_data_.previous_q = generalized_positions;
    osc_data_.previous_qd = generalized_velocities;
}

// ===============================================================================================================
void OSCNode::update_optimization_data(const Eigen::Vector<double, model::contact_site_ids_size>& force_limits) {
    auto mass_matrix = matrix_utils::transformMatrix<double, model::nv_size, model::nv_size, matrix_utils::ColumnMajor>(osc_data_.mass_matrix.data());
    auto coriolis_matrix = matrix_utils::transformMatrix<double, model::nv_size, 1, matrix_utils::ColumnMajor>(osc_data_.coriolis_matrix.data());
    auto contact_jacobian = matrix_utils::transformMatrix<double, model::nv_size, optimization::z_size, matrix_utils::ColumnMajor>(osc_data_.contact_jacobian.data());
    auto taskspace_jacobian = matrix_utils::transformMatrix<double, optimization::s_size, model::nv_size, matrix_utils::ColumnMajor>(osc_data_.taskspace_jacobian.data());
    auto taskspace_bias = matrix_utils::transformMatrix<double, optimization::s_size, 1, matrix_utils::ColumnMajor>(osc_data_.taskspace_bias.data());
    auto desired_taskspace_ddx = matrix_utils::transformMatrix<double, model::site_ids_size, 6, matrix_utils::ColumnMajor>(taskspace_targets_.data());
    
    auto Aeq_matrix = evaluate_function<AeqParams>(Aeq_ops, {design_vector_.data(), mass_matrix.data(), coriolis_matrix.data(), contact_jacobian.data()});
    auto beq_matrix = evaluate_function<beqParams>(beq_ops, {design_vector_.data(), mass_matrix.data(), coriolis_matrix.data(), contact_jacobian.data()});
    // Updates to integrate the required Casadi Inputs mapping into C++ arrays
    auto Aineq_matrix = evaluate_function<AineqParams>(Aineq_ops, {design_vector_.data(), const_cast<double*>(force_limits.data())});
    auto bineq_matrix = evaluate_function<bineqParams>(bineq_ops, {design_vector_.data(), const_cast<double*>(force_limits.data())});
    auto H_matrix = evaluate_function<HParams>(H_ops, {design_vector_.data(), desired_taskspace_ddx.data(), taskspace_jacobian.data(), taskspace_bias.data()});
    auto f_matrix = evaluate_function<fParams>(f_ops, {design_vector_.data(), desired_taskspace_ddx.data(), taskspace_jacobian.data(), taskspace_bias.data()});

    opt_data_.H = H_matrix;
    opt_data_.f = f_matrix;
    opt_data_.Aeq = Aeq_matrix;
    opt_data_.beq = beq_matrix;
    opt_data_.Aineq = Aineq_matrix;
    opt_data_.bineq = bineq_matrix;
}

// ===============================================================================================================
absl::Status OSCNode::set_up_optimization(const Vector<model::contact_site_ids_size>& force_limits) {
    MatrixColMajor<optimization::constraint_matrix_rows, optimization::constraint_matrix_cols> A;
    A << opt_data_.Aeq, opt_data_.Aineq, Abox_;
    Vector<optimization::bounds_size> lb;
    Vector<optimization::bounds_size> ub;
    Vector<optimization::z_size> z_lb_masked = z_lb_;
    Vector<optimization::z_size> z_ub_masked = z_ub_;
    
    for(int i = 0; i < model::contact_site_ids_size; i++) {
        if (force_limits(i) <= 0.001) {
            z_lb_masked(Eigen::seqN(3 * i, 3)).setZero();
            z_ub_masked(Eigen::seqN(3 * i, 3)).setZero();
        } else {
            z_ub_masked(3 * i + 2) = force_limits(i);
        }
    }

    lb << opt_data_.beq, bineq_lb_, dv_lb_, u_lb_, z_lb_masked;
    ub << opt_data_.beq, opt_data_.bineq, dv_ub_, u_ub_, z_ub_masked;
    
    Eigen::SparseMatrix<double> sparse_H = opt_data_.H.sparseView();
    Eigen::SparseMatrix<double> sparse_A = A.sparseView();
    sparse_H.makeCompressed();
    sparse_A.makeCompressed();

    instance_.objective_matrix = sparse_H;
    instance_.objective_vector = opt_data_.f;
    instance_.constraint_matrix = sparse_A;
    instance_.lower_bounds = lb;
    instance_.upper_bounds = ub;
    
    absl::Status result = solver_.Init(instance_, settings_);
    return result;
}

// ===============================================================================================================

absl::Status OSCNode::update_optimization(const Vector<model::contact_site_ids_size>& force_limits) {
    MatrixColMajor<optimization::constraint_matrix_rows, optimization::constraint_matrix_cols> A;
    A << opt_data_.Aeq, opt_data_.Aineq, Abox_;
    Vector<optimization::bounds_size> lb;
    Vector<optimization::bounds_size> ub;
    Vector<optimization::z_size> z_lb_masked = z_lb_;
    Vector<optimization::z_size> z_ub_masked = z_ub_;

    for(int i = 0; i < model::contact_site_ids_size; i++) {
        if (force_limits(i) <= 0.001) {
            z_lb_masked(Eigen::seqN(3 * i, 3)).setZero();
            z_ub_masked(Eigen::seqN(3 * i, 3)).setZero();
        } else {
            z_ub_masked(3 * i + 2) = force_limits(i);
        }

    }
    
    lb << opt_data_.beq, bineq_lb_, dv_lb_, u_lb_, z_lb_masked;
    ub << opt_data_.beq, opt_data_.bineq, dv_ub_, u_ub_, z_ub_masked;
    
    Eigen::SparseMatrix<double> sparse_H = opt_data_.H.sparseView();
    Eigen::SparseMatrix<double> sparse_A = A.sparseView();
    sparse_H.makeCompressed();
    sparse_A.makeCompressed();


    absl::Status result;
    auto sparsity_check = solver_.UpdateObjectiveAndConstraintMatrices(sparse_H, sparse_A);
    if(sparsity_check.ok()) {
        result.Update(solver_.SetObjectiveVector(opt_data_.f));
        result.Update(solver_.SetBounds(lb, ub));
    } else {
        instance_.objective_matrix = sparse_H;
        instance_.objective_vector = opt_data_.f;
        instance_.constraint_matrix = sparse_A;
        instance_.lower_bounds = lb;
        instance_.upper_bounds = ub;
        result.Update(solver_.Init(instance_, settings_));
        result.Update(solver_.SetWarmStart(solution_, dual_solution_));
    }

    return result;
}

// ===============================================================================================================
bool OSCNode::solve_optimization() {
    exit_code_ = solver_.Solve();
    
    if (exit_code_ == OsqpExitCode::kOptimal) {
        solution_ = solver_.primal_solution();
        dual_solution_ = solver_.dual_solution();
        return true;
    } else {
        // Clear the solution so we don't accidentally use old values
        solution_.setZero(); 
        
        RCLCPP_WARN(this->get_logger(), "OSQP Solve Failed. Exit Code: %d", static_cast<int>(exit_code_));
        return false;
    }
}

// ===============================================================================================================
void OSCNode::reset_optimization() {
    Vector<optimization::constraint_matrix_cols> primal_vector = Vector<optimization::constraint_matrix_cols>::Zero();
    Vector<optimization::constraint_matrix_rows> dual_vector = Vector<optimization::constraint_matrix_rows>::Zero();
    std::ignore = solver_.SetWarmStart(primal_vector, dual_vector);
}

// ===============================================================================================================

void OSCNode::publish_torque_command(bool safety_override_active_local, 
                                     std::chrono::time_point<std::chrono::high_resolution_clock> state_read_time_local) 
{
    // --- Constants ---
    const std::set<std::string> reversed_joints_ = {
        "rear_left_hip", "rear_left_knee", "front_left_hip", "front_left_knee"};
    const std::array<std::string, model::nu_size> MOTOR_NAMES = {
        "rear_left_hip", "rear_left_knee", "rear_right_hip", "rear_right_knee",
        "front_left_hip", "front_left_knee", "front_right_hip", "front_right_knee"};
    
    // const double MAX_TORQUE = 25.0;
    const int TORQUE_CONTROL_MODE = 1; 
    const int VELOCITY_CONTROL_MODE = 2; 

    // --- 1. Initialize Command Message ---
    auto command_msg = std::make_unique<Command>(); 
    command_msg->master_gain = 1.0; 
    command_msg->motor_commands.resize(model::nu_size);

    // --- 2. Determine Overall Mode and Populate Commands ---
    if (safety_override_active_local) {
        // SCENARIO A: PERMANENT SAFETY OVERRIDE
        command_msg->high_level_control_mode = 2; 
        
        for (size_t i = 0; i < model::nu_size; ++i) {
            command_msg->motor_commands[i].name = MOTOR_NAMES[i];
            
            // NOTE: HIGH GAIN HOLD IS THE RECOMMENDED SAFETY FIX (using Velocity mode is weak)
            command_msg->motor_commands[i].control_mode = VELOCITY_CONTROL_MODE;
            command_msg->motor_commands[i].position_setpoint = 0.0; 
            command_msg->motor_commands[i].velocity_setpoint = 0.0;
            command_msg->motor_commands[i].feedforward_torque = 0.0; 
            command_msg->motor_commands[i].kp = 0.0; 
            command_msg->motor_commands[i].kd = 0.0;
            command_msg->motor_commands[i].input_mode = 1;   
            command_msg->motor_commands[i].enable = true; 
        }

    } else {
        // SCENARIO B: NORMAL OPERATION 
        command_msg->high_level_control_mode = 2;
        
        // Vector<model::nu_size> osc_torque = solution_(Eigen::seqN(optimization::dv_idx, optimization::u_size));
        Vector<model::nu_size> osc_torque = solution_(Eigen::seqN(optimization::u_idx, optimization::u_size));
        

        std::stringstream ss;
        ss << "OSC Torques: [ ";
        for (size_t i = 0; i < model::nu_size; ++i) {
            ss << std::fixed << std::setprecision(2) << osc_torque(i) << " ";
        }
        ss << "]";
        RCLCPP_INFO(this->get_logger(), "%s", ss.str().c_str());
        

        for (size_t i = 0; i < model::nu_size; ++i) {
            double final_torque = osc_torque(i);
            
            if (reversed_joints_.count(MOTOR_NAMES[i])) {
                final_torque *= -1.0;
            }
            final_torque = std::clamp(final_torque, u_lb_(i), u_ub_(i));

            command_msg->motor_commands[i].name = MOTOR_NAMES[i];
            command_msg->motor_commands[i].control_mode = TORQUE_CONTROL_MODE;
            command_msg->motor_commands[i].feedforward_torque = static_cast<double>(final_torque); 
            // Zero out unused PD terms
            command_msg->motor_commands[i].position_setpoint = 0.0;
            command_msg->motor_commands[i].velocity_setpoint = 0.0; 
            command_msg->motor_commands[i].kp = 0.0; 
            command_msg->motor_commands[i].kd = 0.0;
            command_msg->motor_commands[i].input_mode = 1;   
            command_msg->motor_commands[i].enable = true; 
        }
    }
    
    // --- Time Logging ---
    std::chrono::high_resolution_clock::time_point torque_ready_time_local = std::chrono::high_resolution_clock::now();
    
    // Calculate the control loop latency
    auto latency = std::chrono::duration_cast<std::chrono::microseconds>(
        torque_ready_time_local - state_read_time_local
    );
    double latency_ms = static_cast<double>(latency.count()) / 1000.0;

    // --- 3. Publish ---
    torque_publisher_->publish(std::move(command_msg));

    if (!safety_override_active_local) {
        RCLCPP_INFO(this->get_logger(), 
            "Latency: %.3f ms | SolvCode: %d |**OS Wait: %.3f ms** | Kinematics: %.3f ms | CasADi: %.3f ms | OSQP Solve: %.3f ms | Total Internal: %.3f ms",
            latency_ms,
            static_cast<int>(exit_code_),
            time_wait_for_execution_ms_, // NEW: Shows delay before computation starts            
            time_mujoco_update_ms_, 
            time_casadi_update_ms_, 
            time_osqp_solve_ms_,
            time_mujoco_update_ms_ + time_casadi_update_ms_ + time_osqp_solve_ms_);
    } else {
        RCLCPP_WARN(this->get_logger(), "Safety Override Active. Latency: %.3f ms", latency_ms);
    }
}


void OSCNode::stop_robot() {

    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        safety_override_active_ = true;        
    }    
    
    // 1. Constants for Safety Braking
    const double SAFETY_KP = 0.0; // High stiffness to hold position
    const double SAFETY_KD = 0.0;  // High damping to kill velocity
    const int POSITION_CONTROL_MODE = 3; // Ensure this matches your driver's Enum
    const int VELOCITY_CONTROL_MODE = 2; 

    // 2. Prepare Command
    auto command_msg = std::make_unique<Command>(); 
    command_msg->master_gain = 1.0; 
    command_msg->motor_commands.resize(model::nu_size);
    command_msg->high_level_control_mode = 2; // Safety / Idle mode

    const std::array<std::string, model::nu_size> MOTOR_NAMES = {
        "rear_left_hip", "rear_left_knee", "rear_right_hip", "rear_right_knee",
        "front_left_hip", "front_left_knee", "front_right_hip", "front_right_knee"};

    // 4. Fill Command: "Freeze at current position"
    for (size_t i = 0; i < model::nu_size; ++i) {
        command_msg->motor_commands[i].name = MOTOR_NAMES[i];
        
        // Switch to Position Mode (or Velocity Mode with 0 target)
        command_msg->motor_commands[i].control_mode = VELOCITY_CONTROL_MODE; 
        
        // Target = Last known position (Freeze)
        command_msg->motor_commands[i].position_setpoint = 0.0;
        command_msg->motor_commands[i].velocity_setpoint = 0.0;
        command_msg->motor_commands[i].feedforward_torque = 0.0; // CRITICAL: Zero torque
        
        // Set gains high to resist movement
        command_msg->motor_commands[i].kp = SAFETY_KP; 
        command_msg->motor_commands[i].kd = SAFETY_KD;
        command_msg->motor_commands[i].input_mode = 1;   
        command_msg->motor_commands[i].enable = true; 
    }

    // 5. Publish Immediate Stop
    // We assume the publisher is still valid because on_shutdown runs before destruction
    if (torque_publisher_) {
        torque_publisher_->publish(std::move(command_msg));
        RCLCPP_INFO(this->get_logger(), ">>> SAFETY STOP COMMAND SENT <<<");
        
        // Optional: Sleep briefly to ensure message hits the network before process dies
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}