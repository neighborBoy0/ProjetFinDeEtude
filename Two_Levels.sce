m = mode();
mode(-1);
path = get_absolute_file_path("Two_Levels.sce");
exec(path + "Functions.sce");

// geneation of trajectoire by bigenetic algorithm
//
// destination: The end point of the end joint.
// obstacle: The Cartesian space coordinates of the obstacle.
// joints_origin: The angle of the joint at the previous position.
function Two_Levels(destination, obstacle, joints_origin)
    Xc = [0 0 0];                 // found by the second level

    size_obstacle = size(obstacle);
    M = size_obstacle(1,1)        // M is the number of obstacles
    N = size(robot.links);        // N is the total number of control points
    
    CPosition = robot.links(N).r;           // the current position of the end effector
    FPosition = destination;                // the final position to be reached
    
    population_size_1 = first_level_params(1,1);       // Population size of the 1st level
    number_generations_1 = first_level_params(1,2);     // Number of generations of the 1st level
    mutation_rate_1 = first_level_params(1,3);        // Mutation rate of the 1st level
    crossover_probability_1 = first_level_params(1,4);// Crossover probability of the 1st level
    
    // disp("number of generations 1st level: ", number_generations_1);
    
    //global theta_result;
    global i;
    //theta_result = [];
    i = 0;
    
    try
        // add parameters to the 1st level
        ga_params_1 = init_param();
        ga_params_1 = add_param(ga_params_1, "dimension", 3);
        ga_params_1 = add_param(ga_params_1, 'minbound', CPosition-[0.1 0.1 0.1]);
        ga_params_1 = add_param(ga_params_1, 'maxbound', CPosition+[0.1 0.1 0.1]);
        ga_params_1 = add_param(ga_params_1, 'N', N);
        ga_params_1 = add_param(ga_params_1, 'joints_origin', joints_origin);
    catch
        [error_message,error_number]=lasterror(%t)
        disp("There is an error in function TwoLevel when add some params to the 1st algo genetic, error message:" + error_message);
    end
        
    try
        // the goal function of 1st level, and incoming parameters.
        myObjFun_1= list(firstLevel, FPosition, Xc, M, N, obstacle);
        
        path = get_absolute_file_path("main.sce");
        exec(path + "my_optim_moga.sce");
        
        [pop_opt_1, fobj_pop_opt_1, pop_init_1, fobj_pop_init_1, theta_results] = my_optim_moga(myObjFun_1, population_size_1, number_generations_1, mutation_rate_1, crossover_probability_1, %T, ga_params_1);
        [fmin_1,k_1] = min(fobj_pop_opt_1);        // the result of algo genetic
        xmin = pop_opt_1(k_1(1));
        disp(xmin);
        
        disp(theta_results);
        
    catch
        [error_message,error_number]=lasterror(%t)
        disp("There is an error in function TwoLevel when the 1st algo genetic runs, error message:" + error_message);
    end

endfunction

// The objectif function of the first level
//
// CPosition: the current position of the end effector
// FPosition: the final position of the end effector
// Xc: used to take into account the error found by both levels for the same position of the end effector. found by the 2nd level
// M: the number of obstacles
// N: the total number of control points
function result = firstLevel(CPosition, FPosition, Xc, M, N, obstacle)
    // 1st level in book: min(x)F(X, theta) = alpha * F1(X) + beta * (X - Xc) + gamma * F2(theta)
    result = alpha_variable * F1(CPosition, FPosition) + beta_variable * (CPosition - Xc) + gamma_variable * F2(M, N, obstacle);
endfunction

// The objectif function of the second level
// 
// Q: Angle de joint, c'est variable indépendant.
// N: The number of joints
// joints_origin: The angle of the joint at the previous position.
function result = secondLevel(Q, N, joints_origin)
    result = delta_variable * F3(Q) + zeta_variable * F4(Q, N, joints_origin);
endfunction
