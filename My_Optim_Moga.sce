// Mettre en œuvre un algorithme génétique
//
// ga_f: Une liste qui comprend le pointeur de la fonction cible et plusieurs paramètres ou la fonction objectif.
// pop_size: Le nombre de population.
// nb_generation: Le nombre de generation.
// p_mut: La probabilité de mutation.
// p_cross: La probabilité de croisement.
// Log: Si %T, appellera la fonction de sortie à la fin de chaque itération.
// param: Une liste de paramètres.
// 
// pop_opt: La population d'individus optimaux.
// fobj_pop_opt: L'ensemble des valeurs de fonctions multi-objectifs associées à pop_opt.
// pop_init: La population initiale d'individus.
// fobj_pop_init: L'ensemble des valeurs de fonctions multi-objectifs associées à pop_init.
// theta_results: La solution optimale calculée par l'algorithme génétique de second niveau est exécutée à chaque fois.

function [pop_opt, fobj_pop_opt, pop_init, fobj_pop_init, theta_results] = my_optim_moga(ga_f, pop_size, nb_generation, p_mut, p_cross, Log, param)

    [nargout, nargin] = argn();

    if ~isdef("param","local") then
        param = [];
    end

    [codage_func,err]    = get_param(param,"codage_func",coding_ga_identity);
    [init_func,err]      = get_param(param,"init_func",init_ga_default);
    [crossover_func,err] = get_param(param,"crossover_func",crossover_ga_default);
    [mutation_func,err]  = get_param(param,"mutation_func",mutation_ga_default);
    [selection_func,err] = get_param(param,"selection_func",selection_ga_elitist);
    [nb_couples,err]     = get_param(param,"nb_couples",100);
    [pressure,err]       = get_param(param,"pressure",0.05);
    [output_func, err]   = get_param(param, "output_func", output_moga_default);
    theta_results = []

    if ~isdef("ga_f","local") then
        error(gettext("optim_moga: ga_f is mandatory"));
    else
        if typeof(ga_f)=="list" then
            disp(ga_f(2:$));
            deff("y=_ga_f(x)","y=ga_f(1)(x, ga_f(2:$))");
        else
            deff("y=_ga_f(x)","y=ga_f(x)");
        end
    end

    if ~isdef("pop_size","local") then
        pop_size = 100;
    end
    if ~isdef("nb_generation","local") then
        nb_generation = 10;
    end
    if ~isdef("p_mut","local") then
        p_mut = 0.01;
    end
    if ~isdef("p_cross","local") then
        p_cross = 0.7;
    end
    if ~isdef("Log","local") then
        Log = %F;
    end

    // Initialization of the population
    if (Log) then
        printf(gettext("%s: Initialization of the population\n"),"optim_moga");
    end

    Pop = list();
    Pop = init_func(pop_size, param);

    if (nargout>=3) then
        pop_init = Pop;
    end

    // Code the individuals
    Pop = codage_func(Pop,"code",param);

    for i=1:length(Pop)
        MO_FObj_Pop(i,:) = _ga_f(Pop(i));
    end

    // Compute the domination rank
    for i=1:size(MO_FObj_Pop,1)
        Index = 0;
        for j=1:size(MO_FObj_Pop,1)
            Index = Index + double(and(MO_FObj_Pop(i,:)<=MO_FObj_Pop(j,:)) & or(MO_FObj_Pop(i,:)<MO_FObj_Pop(j,:)));
        end
        FObj_Pop(i) = - (Index + 1);
    end

    FObj_Pop_Max = max(FObj_Pop);
    FObj_Pop_Min = min(FObj_Pop);

    // Normalization of the efficiency

    Efficiency = (1 - pressure) * (FObj_Pop_Max - FObj_Pop) / max([FObj_Pop_Max - FObj_Pop_Min %eps]) + pressure;

    if (nargout==4) then
        fobj_pop_init = MO_FObj_Pop;
    else
        fobj_pop_init = -1;
    end

    // The genetic algorithm
    for i=1:nb_generation
        disp("first level")
        //
        // Selection
        //
        Indiv1 = list();
        Indiv2 = list();
        Wheel = cumsum(Efficiency);
        for j=1:nb_couples
            // Selection of the first individual in the couple
            Shoot = grand(1,1,"unf", 0, Wheel($));
            Index = find(Shoot <= Wheel, 1);
            Indiv1(j)           = Pop(Index);
            MO_FObj_Indiv1(j,:) = MO_FObj_Pop(Index,:);
            // Selection of the second individual in the couple
            Shoot = grand(1,1,"unf", 0, Wheel($));
            Index = find(Shoot <= Wheel, 1);
            Indiv2(j)           = Pop(Index);
            MO_FObj_Indiv2(j,:) = MO_FObj_Pop(Index,:);
        end
        //
        // Crossover
        //
        for j=1:nb_couples
            if (p_cross>grand(1,1,"def")) then
                [x1, x2] = crossover_func(Indiv1(j), Indiv2(j),param);
                Indiv1(j) = x1;
                Indiv2(j) = x2;
                ToCompute_I1(j) = %T;
                ToCompute_I2(j) = %T;
            else
                ToCompute_I1(j) = %F;
                ToCompute_I2(j) = %F;
            end
        end
        //
        // Mutation
        //
        for j=1:nb_couples
            if (p_mut>grand(1,1,"def")) then
                x1 = mutation_func(Indiv1(j),param);
                Indiv1(j) = x1;
                ToCompute_I1(j) = %T;
            end
            if (p_mut>grand(1,1,"def")) then
                x2 = mutation_func(Indiv2(j),param);
                Indiv2(j) = x2;
                ToCompute_I2(j) = %T;
            end
        end
        //
        // Computation of the objective functions
        //
        for j=1:length(Indiv1)
            if ToCompute_I1(j) then MO_FObj_Indiv1(j,:) = _ga_f(Indiv1(j)); end
            if ToCompute_I2(j) then MO_FObj_Indiv2(j,:) = _ga_f(Indiv2(j)); end
        end

        // Reinit ToCompute lists
        ToCompute_I1 = ToCompute_I1 & %F;
        ToCompute_I2 = ToCompute_I2 & %F;

        // Compute the domination rank
        for j=1:size(MO_FObj_Indiv1,1)
            // We compute the rank for Indiv1
            Index1 = 0; Index2 = 0; Index3 = 0;
            for k=1:size(MO_FObj_Indiv1,1)
                Index1 = Index1 + double(and(MO_FObj_Indiv1(j,:)<=MO_FObj_Indiv1(k,:)) & or(MO_FObj_Indiv1(j,:)<MO_FObj_Indiv1(k,:)));
                Index2 = Index2 + double(and(MO_FObj_Indiv1(j,:)<=MO_FObj_Indiv2(k,:)) & or(MO_FObj_Indiv1(j,:)<MO_FObj_Indiv2(k,:)));
            end
            for k=1:size(MO_FObj_Pop,1)
                Index3 = Index3 + double(and(MO_FObj_Indiv1(j,:)<=MO_FObj_Pop(k,:)) & or(MO_FObj_Indiv1(j,:)<MO_FObj_Pop(k,:)));
            end
            FObj_Indiv1(j) = - (Index1 + Index2 + Index3 + 1);

            // We compute the rank for Indiv2
            Index1 = 0; Index2 = 0; Index3 = 0;
            for k=1:size(MO_FObj_Indiv1,1)
                Index1 = Index1 + double(and(MO_FObj_Indiv2(j,:)<=MO_FObj_Indiv1(k,:)) & or(MO_FObj_Indiv2(j,:)<MO_FObj_Indiv1(k,:)));
                Index2 = Index2 + double(and(MO_FObj_Indiv2(j,:)<=MO_FObj_Indiv2(k,:)) & or(MO_FObj_Indiv2(j,:)<MO_FObj_Indiv2(k,:)));
            end
            for k=1:size(MO_FObj_Pop,1)
                Index3 = Index3 + double(and(MO_FObj_Indiv2(j,:)<=MO_FObj_Pop(k,:)) & or(MO_FObj_Indiv2(j,:)<MO_FObj_Pop(k,:)));
            end
            FObj_Indiv2(j) = - (Index1 + Index2 + Index3 + 1);
        end

        // We compute the rank for Pop
        for j=1:size(MO_FObj_Pop,1)
            Index1 = 0; Index2 = 0; Index3 = 0;
            for k=1:size(MO_FObj_Indiv1,1)
                Index1 = Index1 + double(and(MO_FObj_Pop(j,:)<=MO_FObj_Indiv1(k,:)) & or(MO_FObj_Pop(j,:)<MO_FObj_Indiv1(k,:)));
                Index2 = Index2 + double(and(MO_FObj_Pop(j,:)<=MO_FObj_Indiv2(k,:)) & or(MO_FObj_Pop(j,:)<MO_FObj_Indiv2(k,:)));
            end
            for k=1:size(MO_FObj_Pop,1)
                Index3 = Index3 + double(and(MO_FObj_Pop(j,:)<=MO_FObj_Pop(k,:)) & or(MO_FObj_Pop(j,:)<MO_FObj_Pop(k,:)));
            end
            FObj_Pop(j) = - (Index1 + Index2 + Index3 + 1);
        end

        //
        // Recombination
        //

        [Pop, FObj_Pop, Efficiency, MO_FObj_Pop] = selection_func(Pop, Indiv1, Indiv2, FObj_Pop, FObj_Indiv1, FObj_Indiv2, ...
                MO_FObj_Pop, MO_FObj_Indiv1, MO_FObj_Indiv2, param);
        if (Log) then
            stop = output_func(i, nb_generation, Pop, MO_FObj_Pop, param);
            if stop then
                break
            end
        end
        
        // 
        //Appelez l'algorithme génétique de deuxième niveau et calculez.
        // 
        population_size_2 = second_level_params(1,1);           // Population size of the 1st level
        number_generations_2 = second_level_params(1,2);        // Number of generations of the 1st level
        mutation_rate_2 = second_level_params(1,3);             // Mutation rate of the 1st level
        crossover_probability_2 = second_level_params(1,4);     // Crossover probability of the 1st level
        
        try
            N = get_param(param,"N",6);
            joints_origin = get_param(param, "joints_origin", [0 0 0 0 0 0]);
            disp(joints_origin);
        catch
            disp("There is an error in read param in 1st level to call 2nd level");
        end
        
        try
            ga_params_2 = init_param();
            ga_params_2 = add_param(ga_params_2, "dimension", 6);

            ga_params_2 = add_param(ga_params_2, 'minbound', min_angle);
            ga_params_2 = add_param(ga_params_2, 'maxbound', max_angle);
        catch
            [error_message,error_number]=lasterror(%t)
            disp("There is an error in function TwoLevel when add some params to the 2nd algo genetic, error message:" + error_message);
        end
        
        try
            myObjFun_2 = list(secondLevel, N, joints_origin);
            [pop_opt_2, fobj_pop_opt_2] = optim_moga(myObjFun_2, population_size_2, number_generations_2, mutation_rate_2, crossover_probability_2, %T, ga_params_2);
            [fmin_2, k_2] = min(fobj_pop_opt_2);        // the result of algo genetic
            theta_min = pop_opt_2(k_2);
            disp(theta_min);
            
            theta_results = [theta_results; theta_min];
            
        catch
            [error_message,error_number]=lasterror(%t)
            disp("There is an error in function TwoLevel when the 1st algo genetic runs, error message:" + error_message);
        end

    end
    pop_opt      = codage_func(Pop, 'decode', param);
    fobj_pop_opt = MO_FObj_Pop;
endfunction
