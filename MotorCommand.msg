# MotorCommand
# Motor Name: (Should match Topic name for the motor)
string name

# Command Data:
float64 position_setpoint
float64 velocity_setpoint
float64 feedforward_torque
float64 kp
float64 kd

# ODrive Command Modes:
uint8 control_mode
uint8 input_mode

# Enable Flag:
bool enable