#include <csignal> // Required for signal handling
#include <thread>  // Required for sleep
#include "rclcpp/rclcpp.hpp"
#include "rules_cc/cc/runfiles/runfiles.h"
#include "operational-space-control/walter_sr_v2/osc_node.h"
#include "rclcpp/executors/multi_threaded_executor.hpp"

// Global pointer allows the static signal handler to access the executor
rclcpp::executors::MultiThreadedExecutor* g_executor = nullptr;

// Custom Signal Handler: Stops the executor loop but KEEPS ROS ALIVE
void signal_handler(int signum) {
    (void)signum;
    if (g_executor) {
        // This causes executor.spin() to return, but context remains valid
        g_executor->cancel(); 
    }
}

int main(int argc, char** argv) {
    // 1. Initialize ROS 2 WITHOUT default signal handlers
    // This prevents "Context Invalid" errors during shutdown
    rclcpp::InitOptions init_options;
    init_options.shutdown_on_signal = false; 
    rclcpp::init(argc, argv, init_options);

    // 2. Install our custom signal handler
    std::signal(SIGINT, signal_handler);

    std::string error;
    std::unique_ptr<rules_cc::cc::runfiles::Runfiles> runfiles(
        rules_cc::cc::runfiles::Runfiles::Create(argv[0], "walter-v2-floortest-z-control", &error));

    if (!error.empty()) {
        std::cerr << "Failed to create runfiles: " << error << std::endl;
        return 1;
    }

    std::filesystem::path model_path = 
        runfiles->Rlocation("mujoco-models~/models/walter_sr/scene_walter_sr_v2.xml");

    // 3. Instantiate Executor and set global pointer
    rclcpp::executors::MultiThreadedExecutor executor;
    g_executor = &executor; 

    try {
        auto node = std::make_shared<OSCNode>(model_path.string());
        executor.add_node(node);
        
        RCLCPP_INFO(rclcpp::get_logger("main"), "Starting MultiThreaded Executor (Priority 99)...");
        RCLCPP_INFO(rclcpp::get_logger("main"), "Press Ctrl+C to trigger Safety Stop.");
        
        // 4. Spin blocks here until Ctrl+C is pressed
        executor.spin(); 

        // =========================================================
        // 5. SAFETY SHUTDOWN SEQUENCE (Runs AFTER Ctrl+C)
        // =========================================================
        // The ROS context is STILL VALID here because we disabled auto-shutdown.
        
        RCLCPP_WARN(rclcpp::get_logger("main"), "Ctrl+C detected. Sending Safety Stop...");

        // Manually call the stop function (Must be public in osc_node.h!)
        node->stop_robot(); 

        // CRITICAL: Wait briefly to ensure the message physically leaves the network buffer
        // before we destroy the publisher context.
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        RCLCPP_INFO(rclcpp::get_logger("main"), "Safety command sent. Shutting down.");

    } catch (const std::exception& e) {
        RCLCPP_FATAL(rclcpp::get_logger("main"), "Exception caught: %s", e.what());
        executor.cancel(); 
        rclcpp::shutdown();        
        return 1;
    }

    // 6. Finally, safe to destroy context
    rclcpp::shutdown();
    return 0;
}