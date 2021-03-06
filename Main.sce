m = mode();
mode(-1);
funcprot(0);

global min_angle;
global max_angle;
global robot;
global alpha_variable;
global beta_variable;
global gamma_variable;
global delta_variable;
global zeta_variable;
global first_level_params;
global second_level_params;

//robot = p560;

path = get_absolute_file_path("Main.sce");
exec(path + "Read_Files.sce");

// Read data.json to get robot params
[origin, destination, obstacles, joints, joints_origin] = Read_Robot(path + "data.json");

// Read param.json to get algo genetic params
[alpha_variable, beta_variable, gamma_variable, delta_variable, zeta_variable, first_level_params, second_level_params] = Read_Params(path + "params.json");

// To get min and max angles of each joints
joints = joints';
min_angle = joints(1,:);
max_angle = joints(2, :);

// trajectoire generation 
exec(path + "Two_Levels.sce");
f = Two_Levels(destination, obstacles, joints_origin);

//disp(theta_result);
