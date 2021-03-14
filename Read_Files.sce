funcprot(0);


// read json file
// 
// origin: a 1*3 matrix, the origin position of the end effector
// destination: a 1*3 matrix, the destination position of the end effector
// obstacle: a n*3 matrix, the position of the obstacles, there are n obstacle
function [origin, destination, obstacles, joints] = Read_Robot(file_name)
    try
        data = mgetl(file_name);
        mystruct = JSONParse(data);
        origin = mystruct.origin;
        destination = mystruct.destination;
        obstacles = mystruct.obstacles;
        joints = mystruct.joints;
    catch
        [error_message,error_number]=lasterror(%t);
        disp("There is an error in function Read_File, error message:" + error_message);
    end
endfunction

function [alpha_variable, beta_variable, gamma_variable, delta_variable, zeta_variable, first_level, second_level] = Read_Params(file_name)
    try
        data = mgetl(file_name);
        my_struct = JSONParse(data);
        alpha_variable = my_struct.alpha;
        beta_variable = my_struct.beta;
        gamma_variable = my_struct.gamma;
        delta_variable = my_struct.delta;
        zeta_variable = my_struct.zeta;
        first_level = my_struct.first_level;
        second_level = my_struct.second_level;
    catch
        [error_message,error_number]=lasterror(%t);
        disp("There is an error in function Read_Params, error message:" + error_message);
    end
endfunction
