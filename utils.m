function func = utils
%FUNCTIONS Summary of this function goes here
%   Detailed explanation goes here
    func.change_topology_model= @change_topology_model;
    func.constraintHandle = @constraintHandle;
    func.fitnessEval_Scenario1_Case1 = @fitnessEval_Scenario1_Case1;
    func.custom_SE = @custom_SE;
    func.BDD = @BDD;
    func.fitnessEval_Scenario2 = @fitnessEval_Scenario2;
    func.BDA_cost = @BDA_cost;
    func.MOEAD = @MOEAD;
    func.PowerFlowEq= @PowerFlowEq;
    func.MOEAD_SCENARIO_2 = @MOEAD_SCENARIO_2;
    func.BDA_cost_Scenario_2=@BDA_cost_Scenario_2;
end



%% Function to change topology

function new_model = change_topology_model(model, z_digital)
    model.branch(:,11) = z_digital;
    new_model = model;
end


function [baseMVA, bus, gen, branch, success, et, z, z_est, error_sqrsum] = custom_SE_with_noise(model, z_measures)

    
    %The available measures are expressed as z = [PFrom; QFrom; PTo; QTo; Pg; Qg; slack_Vm]
    % Make the spread of the Gaussians be 20% of the a values
    sigmas = 0.1 * z_measures; % Also a column vector.
    % Create the noise values that we'll add to a.
    randomNoise = randn(length(z_measures), 1) .* sigmas;
    % Add noise to a to make an output column vector.
    z_analog = z_measures + 0.5*randomNoise;
    
    %% Test state estimation
    idx.idx_zPF = transpose(1:46);
    idx.idx_zPT = transpose(1:46);
    idx.idx_zPG = transpose(1:10);
    idx.idx_zVa = []; % transpose(1:39)
    idx.idx_zQF = transpose(1:46); %transpose(1:46);
    idx.idx_zQT = transpose(1:46);
    idx.idx_zQG = transpose(1:10);
    idx.idx_zVm = 31; %transpose(1:39);
     
    %% specify measurements
    measure.PF = z_analog(1:46); 
    measure.PT = z_analog(93:138); 
    measure.PG = z_analog(185:194); 
    measure.Va = []; 
    measure.QF = z_analog(47:92); 
    measure.QT = z_analog(139:184);
    measure.QG = z_analog(195:204);
    measure.Vm = z_analog(205);
    
    %% specify measurement variances
    sigma.sigma_PF = sigmas(1); 
    sigma.sigma_PT = sigmas(93);
    sigma.sigma_PG = sigmas(185); 
    sigma.sigma_Va = [];
    sigma.sigma_QF = sigmas(47); 
    sigma.sigma_QT = sigmas(139);
    sigma.sigma_QG = sigmas(195);
    sigma.sigma_Vm = sigmas(205);
    
    %% check input data integrity
    nbus = 39;
    [success, measure, idx, sigma] = checkDataIntegrity(measure, idx, sigma, nbus);
    if ~success
         error('State Estimation input data are not complete or sufficient!');
    end
         
    %% run state estimation
    type_initialguess = 2; % flat start
    [baseMVA, bus, gen, branch, success, et, z, z_est, error_sqrsum] = run_se(model, measure, idx, sigma, type_initialguess);
        
end


%% Functions and utilities for the NAA algorithm

function [newInd, newUserObj] = constraintHandle(ind, userObj)

    display('---------- CONSTRAINT EVALUATION --------------------')
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B,RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST,ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
    [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN,MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX,QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
    
    
    display('Individual');
    display(ind);

    %% Assignment fo the values given from the global variable UserObj
    tau = userObj.tau;
    sigmas = userObj.sigmas;
    nominal_model = userObj.nominal_model;
    z_analog = userObj.z_analog;
    z_digital = userObj.z_digital;
    omega_r_set = userObj.omega_r_set;
    omega_a_set = userObj.omega_a_set;
    omega_r_set_size = length(omega_r_set);
    omega_a_set_size = length(omega_a_set);
    a_limit = userObj.a_limit;
            
    %% Now we are taking the attack vector a from the individual and it has to pass the BDD Test 
    a = ind(omega_r_set_size+omega_a_set_size+1 : length(ind));
    a = transpose(a);
    for i=1:size(a)
            if (-a_limit<a(i)) && (a(i)<a_limit)
                a(i)=0;
            end 
    end
    display('a display');
    display(a);
    z_attack = a + z_analog;

    %% Take the current omega_plus and omega_minus sets 
    %% The ind is composed by |omega_minus||omega_plus||a vector|
    %% The indices corresponds to |omega_r||omega_a|
    omega_minus = ind(1:omega_r_set_size);
    omega_plus = ind(omega_r_set_size+1 : omega_a_set_size+omega_r_set_size);

    display('Omega_Plus');
    display(sum(omega_plus));
    display('Omega_Minus');
    display(sum(omega_minus));

    %% Test BDD

    success = BDD(z_attack, nominal_model, tau, sigmas);

    display('BDD result')
    display(success)

    if(success)
        newInd = ind;
        newInd(omega_r_set_size+omega_a_set_size+1 : length(ind)) = transpose(a);      
    else
        [a, success] = make_attack_vector(nominal_model, omega_r_set, omega_a_set, omega_plus, omega_minus, z_digital, z_analog, sigmas, a_limit);
        newInd(omega_r_set_size+omega_a_set_size+1 : length(ind)) = transpose(a);
        if(success == 0)
            newInd = -1;
        end
       
    end

    newUserObj = userObj;
    
end

%% NAA - Scenario 1 Case 1_Case_2
function fitness = fitnessEval_Scenario1_Case1(ind, userObj)

    display('------------ FITNESS EVALUATION --------------');
    
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B,RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST,ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
    [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN,MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX,QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
    
    display('Individual');
    display(ind);

    %% Assignment fo the values given from the global variable UserObj
    tau = userObj.tau;
    sigmas = userObj.sigmas;
    nominal_model = userObj.nominal_model;
    z_analog = userObj.z_analog;
    z_digital = userObj.z_digital;
    omega_r_set = userObj.omega_r_set;
    omega_a_set = userObj.omega_a_set;
    omega_r_set_size = length(omega_r_set);
    omega_a_set_size = length(omega_a_set);

    display(find(ind(1:omega_r_set_size) == 1));
    display(size(omega_a_set));

    %% The ind is composed by  |omega_minus||omega_plus||a vector|
    %% The indices corresponds to |omega_r||omega_a   |
    omega_minus_index = omega_r_set(find(ind(1:omega_r_set_size) == 1));
    omega_plus_index = omega_a_set(find(ind(omega_r_set_size+1 : omega_r_set_size + omega_a_set_size) == 1));

    %% We change the status of the transmission lines

    z_digital(omega_plus_index) = 1;
    z_digital(omega_minus_index) = 0;

    Pg_star = userObj.Pg_star;
    lmp_star = userObj.lmp_star;
    PFrom_star = userObj.PFrom_star;


    %% Creating the new model with the new topology 
    attack_model = utils().change_topology_model(nominal_model, z_digital);

    %% Solving the SCED model under attack
    display('RUN OPF ON ATTACK TOPOLOGY');
    opf_attack_model_res = runopf(attack_model);
    if (opf_attack_model_res.success == 0)
        fitness = +Inf;
        return;
    end
    
    %% Calculating the Power Generated and the other parameters "crossed" (in the system attacked")
    Pg_cross = opf_attack_model_res.gen(:,PG);
    lmp_cross = opf_attack_model_res.bus(:,LAM_P);
    PFrom_cross = opf_attack_model_res.branch(:,PF);
    
    %% We substitute the Power Generated in the nominal model
    nominal_model2 = nominal_model;
    nominal_model2.gen(:,PG) = Pg_cross;

    display('RUN POWER FLOW EQUATION WITH NOMINAL TOPOLOGY');
    [MVAbase, result.bus, result.gen, result.branch, success, et] = runpf(nominal_model2);
    
    %% Showing results
    display(result);
    PFrom_circled_times = result.branch(:,PF);
    Pg_circled_times = result.gen(:,PG);
    display('Old_PG')
    display(Pg_cross);
    display(PFrom_circled_times);
    display('New_PG');
    display(Pg_circled_times);

    %% Set fitness function the fitness function i relative to Scenario 1 - Case1_Case2
    gencost = nominal_model.gencost;
    % case 1 and case 2 ( two different fitness functions)
%     fitness = sum(totcost(gencost, Pg_circled_times)) - sum(totcost(gencost,Pg_star));
    fitness = sum(lmp_cross-lmp_star);
    eta = (fitness / sum(totcost(gencost,Pg_star) ) ) * 100;
    display(sum(totcost(gencost, Pg_circled_times)));
    display(sum(totcost(gencost,Pg_star)));
    display("ETA VALUE");
    display(eta);
    display(fitness);
    
    %% Since the NAA makes the min, to have te max is necessary to change the sign of the fitness function
    fitness = -fitness;

end

%% NAA - Scenario 2
function fitness = fitnessEval_Scenario2(ind, userObj)

    display('------------ FITNESS EVALUATION --------------');
    
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B,RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST,ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
    [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN,MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX,QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
    
    display('Individual');
    display(ind);

    %% Assignment fo the values given from the global variable UserObj
    nominal_model = userObj.nominal_model;
    z_digital = userObj.z_digital;
    omega_r_set = userObj.omega_r_set;
    omega_a_set = userObj.omega_a_set;
    omega_r_set_size = length(omega_r_set);
    omega_a_set_size = length(omega_a_set);

    display(find(ind(1:omega_r_set_size) == 1));
    display(size(omega_a_set));

    %% The ind is composed by  |omega_minus||omega_plus||a vector|
    %% The indices corresponds to |omega_r||omega_a   |
    omega_minus_index = omega_r_set(find(ind(1:omega_r_set_size) == 1));
    omega_plus_index = omega_a_set(find(ind(omega_r_set_size+1 : omega_r_set_size + omega_a_set_size) == 1));

    %% We change the status of the transmission lines
    z_digital(omega_plus_index) = 1;
    z_digital(omega_minus_index) = 0;
    
    %% We change the status of the transmission lines
    Pg_star = userObj.Pg_star;
    lmp_star = userObj.lmp_star;
    PFrom_star = userObj.PFrom_star;

    %% Creating the new model with the new topology 
    attack_model = utils().change_topology_model(nominal_model, z_digital);

    %% Solving the SCED model under attack
    display('RUN OPF ON ATTACK TOPOLOGY');
    opf_attack_model_res = runopf(attack_model);
    if (opf_attack_model_res.success == 0)
        fitness = +Inf;
        return;
    end

    %% Calculating the Power Generated and the other parameters "crossed" (in the system attacked")
    Pg_cross = opf_attack_model_res.gen(:,PG);
    lmp_cross = opf_attack_model_res.bus(:,LAM_P);
    PFrom_cross = opf_attack_model_res.branch(:,PF);
    
    %% We substitute the Power Generated in the nominal model
    nominal_model2 = nominal_model;
    nominal_model2.gen(:,PG) = Pg_cross;

    display('RUN POWER FLOW EQUATION WITHNOMINAL TOPOLOGY');
    [MVAbase, result.bus, result.gen, result.branch, success, et] = runpf(nominal_model2);

    %% Taking results
    display(result);
    PFrom_circled_times = result.branch(:,PF);
    Pg_circled_times = result.gen(:,PG);
    
    display('Old_PG')
    display(Pg_cross);
    display(PFrom_circled_times);
    

    display('New_PG');
    display(Pg_circled_times);

    %% Set fitness function the fitness function relative to Scenario 2 
    gencost = nominal_model.gencost;
    eta = ( (sum(totcost(gencost, Pg_circled_times)) - sum(totcost(gencost,Pg_star))) / sum(totcost(gencost,Pg_star) ) ) * 100;   
    fitness = sum(PFrom_circled_times - PFrom_star);
    display("ETA VALUE");
    display(eta);
    display(fitness);
    fitness = -fitness;

end

%% State Estimation
function [H, baseMVA, bus, gen, branch, success, et, z, z_est, error_sqrsum] = custom_SE(model, z_analog,sigmas)
    
    
    
    %% Test state estimation
    idx.idx_zPF = transpose(1:46);
    idx.idx_zPT = transpose(1:46);
    idx.idx_zPG = transpose(1:10);
    idx.idx_zVa = []; 
    idx.idx_zQF = transpose(1:46); 
    idx.idx_zQT = transpose(1:46);
    idx.idx_zQG = transpose(1:10);
    idx.idx_zVm = []; ;
     
    %% specify measurements
    branch_nb = 46;
    gen_nb = 10;
    measure.PF = z_analog(1:branch_nb); 
    measure.PT = z_analog(branch_nb+1:branch_nb*2); 
    measure.PG = z_analog(branch_nb*2+1:branch_nb*2+gen_nb); 
    measure.Va = []; 
    measure.QF = z_analog(branch_nb*2+gen_nb+1:branch_nb*3+gen_nb); 
    measure.QT = z_analog(branch_nb*3+gen_nb+1:branch_nb*4+gen_nb); 
    measure.QG = z_analog(branch_nb*4+gen_nb+1:branch_nb*4+gen_nb*2);
    measure.Vm = []; 
    
    %% specify measurement variances
    sigma.sigma_PF = sigmas(1);
    sigma.sigma_PT = sigmas(93);
    sigma.sigma_PG = sigmas(185);
    sigma.sigma_Va = [];
    sigma.sigma_QF = sigmas(47);
    sigma.sigma_QT = sigmas(139);
    sigma.sigma_QG = sigmas(195);
    sigma.sigma_Vm = [];
    
    %% check input data integrity
    nbus = 39;
    [success, measure, idx, sigma] = checkDataIntegrity(measure, idx, sigma, nbus);
    if ~success
         error('State Estimation input data are not complete or sufficient!');
    end
         
    %% run state estimation
    type_initialguess = 2; 
    display('DATA FOR STATE ESTIMATION');
    display(model);
    display(measure);
    [H, baseMVA, bus, gen, branch, success, et, z, z_est, error_sqrsum] = run_se(model, measure, idx, sigma, type_initialguess);
        
end

%% This function creates the attack vector a
function [a, success] = make_attack_vector(model, omega_r_set, omega_a_set, omega_plus, omega_minus, z_digital, z_analog, sigmas, a_limit)
        
        [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
        [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B,RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST,ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
        [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN,MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX,QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
        
        
        %% Take the index of digital measurmeents to modify
        omega_plus_index = omega_a_set(find(omega_plus == 1));
        omega_minus_index = omega_r_set(find(omega_minus == 1));

        
        %% Chenging the status of the transmission lines in the set of the removable and in the set of the ones not currently used
        z_digital(omega_plus_index) = 1;
        z_digital(omega_minus_index) = 0;

        %% New model with the new topology
        attack_model = utils().change_topology_model(model, z_digital);
        display('Print attack model branch status');
        display(attack_model.branch(:, BR_STATUS));
       

        %% The attacker, as the control center, solves the same opf model of the control center
        
        opf_attack_model_res = runopf(attack_model);
        
        % The results of the new opf are taken: this process is needed to make the state estimation with the attack model. 
        % Then it is possible to take the H matrix relative to the attack model

        %% Get generation levels
        Pg_attack_model = opf_attack_model_res.gen(:,PG);
        Qg_attack_model = opf_attack_model_res.gen(:,QG);
        
        %% Measurement vector 
        PFrom_attack_model = opf_attack_model_res.branch(:,PF);
        PTo_attack_model = opf_attack_model_res.branch(:,PT);
        QFrom_attack_model = opf_attack_model_res.branch(:,QF);
        QTo_attack_model = opf_attack_model_res.branch(:,QT);

        z_measures_attack_model = [PFrom_attack_model; PTo_attack_model; Pg_attack_model; QFrom_attack_model; QTo_attack_model; Qg_attack_model];
        
        %% Adding noise 
        % The available measures are expressed as z = [PFrom; QFrom; PTo; QTo; Pg; Qg; slack_Vm]
        % Make the spread of the Gaussians be 20% of the a values
        sigmas_attack_model = 0.2 * z_measures_attack_model; % Also a column vector.
        % Create the noise values that we'll add to a.
        randomNoise_attack_model = randn(length(z_measures_attack_model), 1) .* sigmas_attack_model;
        % Add noise to a to make an output column vector.
        z_analog_attack_model = z_measures_attack_model + randomNoise_attack_model;
        
        try
            [H_attack_model, baseMVA, bus_attack_model, gen_attack_model, branch_attack_model, success, et, z_attack_model, z_est_attack_model, error_sqrsum_attack_model] = utils().custom_SE(attack_model, z_analog_attack_model, sigmas_attack_model);
        catch ME

            display(ME.identifier);
            a = 0;
            success = 0;
            return
        end
        
        % The attacker makes its own estimate
        %% State preserving attack
        % The attacker makes the state  estimation of the normal operation state
        % The estimate of the atacker on the normal grid is made by:
        [H_normal_by_attacker, baseMVA, bus_by_attacker, gen_by_attacker, branch_by_attacker, success, et, z_by_attacker, z_est_by_attacker, error_sqrsum_normal] = utils().custom_SE(model, z_analog, sigmas);
        
        
        % Once obtained the H_attack_matrix and the H_normal_by_attacker matrix.
        % Therefore it is possible to compute the a vector as shown in 
        % [On Topology Attack of a Smart Grid:Undetectable Attacks and Countermeasures]
        
        a = ( ((H_attack_model-H_normal_by_attacker)*inv(H_normal_by_attacker' * H_normal_by_attacker))*H_normal_by_attacker' )*z_analog;
        
        % The a vector so obtained has many entries different
        % from 0. Caused by the presence of the measurement errors in z_analog and
        % since the H matrices are not obtained using the power flow equation but
        % using an iteration algorithm. An Euristic approach is to set to zero the
        % components of the a vector with a relatively small magnitude. 

        for i=1:size(a)
            if (-a_limit<a(i)) && (a(i)<a_limit)
                a(i)=0;
            end
        end
end

%% BDD test
function success = BDD(z_attack, model, tau, sigmas)
    
   

    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B,RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST,ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
    [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN,MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX,QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

    [Hl, baseMVA, bus, gen, branch, success, et, z, z_est, error_sqrsum] = utils().custom_SE(model, z_attack, sigmas);
    success = 0;
    if error_sqrsum < tau
        success = 1;

    end
end

function F= PowerFlowEq(v,t)         
         r12=0.01;
         r13=0.05;
         r23=0.05;
         g12=1/r12;
         g13=1/r13;
         g23=1/r23;
         b12=0.05;
         b13=0.025;
         b23=0.025;
         F=[-v(1)*v(2)*(g12*cos((t(1)-t(2)))+b12*sin((t(1)-t(2))))+g12*v(1)*v(1);
             +v(1)*v(2)*(b12*cos((t(1)-t(2)))-g12*sin((t(1)-t(2))))+b12*v(1)*v(1);
             -v(1)*v(3)*(g13*cos((t(1)-t(3)))+b13*sin((t(1)-t(3))))+g13*v(1)*v(1);
             +v(1)*v(3)*(b13*cos((t(1)-t(3)))-g13*sin((t(1)-t(3))))+b13*v(1)*v(1);
             -v(2)*v(3)*(g23*cos((t(2)-t(3)))+b23*sin((t(2)-t(3))))+g23*v(2)*v(2);
             +v(2)*v(3)*(b23*cos((t(2)-t(3)))-g23*sin((t(2)-t(3))))+b23*v(2)*v(2)
             ];

end

%% Cost function Scenario 1_Case 1_Case_2 ( For GA and DragonFly)
function o = BDA_cost(x)
    
    global userObj;
    ind = x;
    display('------------ FITNESS EVALUATION --------------');
    
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B,RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST,ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
    [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN,MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX,QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

    
    nominal_model = userObj.nominal_model;
    z_digital = userObj.z_digital;
    omega_r_set = userObj.omega_r_set;
    omega_a_set = userObj.omega_a_set;
    omega_r_set_size = length(omega_r_set);
    omega_a_set_size = length(omega_a_set);

    %% Take the index of digital measurmeents to modify
    display(find(ind(1:omega_r_set_size) == 1));
    display(size(omega_a_set));

    %% The ind is composed by         |omega_minus||omega_plus||a vector|
    %% The indices corresponds to     |omega_r    ||omega_a   |
    ind(omega_r_set_size+1 : omega_r_set_size + omega_a_set_size)
    omega_minus_index = omega_r_set(find(ind(1:omega_r_set_size) == 1));
    omega_plus_index = omega_a_set(find(ind(omega_r_set_size+1 : omega_r_set_size + omega_a_set_size) == 1));
    omega_minus = ind(1:omega_r_set_size);
    omega_plus = ind(omega_r_set_size+1 : omega_a_set_size+omega_r_set_size);
    display('Omega_Plus');
    display(sum(omega_plus));
    display('Omega_Minus');
    display(sum(omega_minus));   

    %%  We change the status of the transmission lines

    z_digital(omega_plus_index) = 1;
    z_digital(omega_minus_index) = 0;
    
    %% Taking results
    Pg_star = userObj.Pg_star;
    lmp_star = userObj.lmp_star;
    PFrom_star = userObj.PFrom_star;

    %% Creating the attack model
    attack_model = utils().change_topology_model(nominal_model, z_digital);

    %% Solving the SCED model under attack
    display('RUN OPF ON ATTACK TOPOLOGY');
    opf_attack_model_res = runopf(attack_model);
    if (opf_attack_model_res.success == 0)
        o = +Inf;
        return;
    end

    Pg_cross = opf_attack_model_res.gen(:,PG);
    lmp_cross = opf_attack_model_res.bus(:,LAM_P);
    PFrom_cross = opf_attack_model_res.branch(:,PF);
    
    %% Substitute the generationm power and run pf problem on the nominal model
    nominal_model2 = nominal_model;
    nominal_model2.gen(:,PG) = Pg_cross;

    display('RUN POWER FLOW EQUATION WITHNOMINAL TOPOLOGY');
    [MVAbase, result.bus, result.gen, result.branch, success, et] = runpf(nominal_model2);  

    %% Displaying results
    display(result);
    PFrom_circled_times = result.branch(:,PF);
    Pg_circled_times = result.gen(:,PG);
    
    display('Old_PG')
    display(Pg_cross);
    display(PFrom_circled_times);
    

    display('New_PG');
    display(Pg_circled_times);

    %% Setting the fitness value for Scenario 1- Case 1_Case_2
    gencost = nominal_model.gencost;
%     fitness = sum(totcost(gencost, Pg_circled_times)) - sum(totcost(gencost,Pg_star));
    fitness = sum(lmp_cross-lmp_star);
    eta = ( (sum(totcost(gencost, Pg_circled_times)) - sum(totcost(gencost,Pg_star))) / sum(totcost(gencost,Pg_star) ) ) * 100;
    display(sum(totcost(gencost, Pg_circled_times)));
    display(sum(totcost(gencost,Pg_star)));
    display("ETA VALUE");
    display(eta);
    display(fitness);
    
    fitness = -fitness;

    o = fitness;
end



%% Cost function Scenario_2( For GA and DragonFly)
function o = BDA_cost_Scenario_2(x)
    
    global userObj;
    ind = x;
    display('------------ FITNESS EVALUATION --------------');
    
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B,RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST,ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
    [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN,MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX,QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

    
    nominal_model = userObj.nominal_model;
    z_digital = userObj.z_digital;
    omega_r_set = userObj.omega_r_set;
    omega_a_set = userObj.omega_a_set;
    omega_r_set_size = length(omega_r_set);
    omega_a_set_size = length(omega_a_set);

    %% Take the index of digital measurmeents to modify
    display(find(ind(1:omega_r_set_size) == 1));
    display(size(omega_a_set));

    %% The ind is composed by         |omega_minus||omega_plus||a vector|
    %% The indices corresponds to     |omega_r    ||omega_a   |
    ind(omega_r_set_size+1 : omega_r_set_size + omega_a_set_size)
    omega_minus_index = omega_r_set(find(ind(1:omega_r_set_size) == 1));
    omega_plus_index = omega_a_set(find(ind(omega_r_set_size+1 : omega_r_set_size + omega_a_set_size) == 1));
    omega_minus = ind(1:omega_r_set_size);
    omega_plus = ind(omega_r_set_size+1 : omega_a_set_size+omega_r_set_size);
    display('Omega_Plus');
    display(sum(omega_plus));
    display('Omega_Minus');
    display(sum(omega_minus));   

    %%  We change the status of the transmission lines

    z_digital(omega_plus_index) = 1;
    z_digital(omega_minus_index) = 0;
    
    %% Taking results
    Pg_star = userObj.Pg_star;
    lmp_star = userObj.lmp_star;
    PFrom_star = userObj.PFrom_star;

    %% Creating the attack model
    attack_model = utils().change_topology_model(nominal_model, z_digital);

    %% Solving the SCED model under attack
    display('RUN OPF ON ATTACK TOPOLOGY');
    opf_attack_model_res = runopf(attack_model);
    if (opf_attack_model_res.success == 0)
        o = +Inf;
        return;
    end

    Pg_cross = opf_attack_model_res.gen(:,PG);
    lmp_cross = opf_attack_model_res.bus(:,LAM_P);
    PFrom_cross = opf_attack_model_res.branch(:,PF);
    
    %% Substitute the generationm power and run pf problem on the nominal model
    nominal_model2 = nominal_model;
    nominal_model2.gen(:,PG) = Pg_cross;

    display('RUN POWER FLOW EQUATION WITHNOMINAL TOPOLOGY');
    [MVAbase, result.bus, result.gen, result.branch, success, et] = runpf(nominal_model2);  

    %% Displaying results
    display(result);
    PFrom_circled_times = result.branch(:,PF);
    Pg_circled_times = result.gen(:,PG);
    
    display('Old_PG')
    display(Pg_cross);
    display(PFrom_circled_times);
    

    display('New_PG');
    display(Pg_circled_times);

    %% Setting the fitness value for Scenario 2
    gencost = nominal_model.gencost;
    fitness = norm(PFrom_circled_times - PFrom_star);
    eta = ( (sum(totcost(gencost, Pg_circled_times)) - sum(totcost(gencost,Pg_star))) / sum(totcost(gencost,Pg_star) ) ) * 100;
    display(sum(totcost(gencost, Pg_circled_times)));
    display(sum(totcost(gencost,Pg_star)));
    display("ETA VALUE");
    display(eta);
    display(fitness);
    
    fitness = -fitness;
    o = fitness;
end


%% MOEAD-(Solving two objective functions: one about the Case 1 andthe other about the Case 2)
function o=MOEAD(x)
    
    display(x);
    size(x)
    global userObj;
    ind = x;
    display('------------ FITNESS EVALUATION --------------');
    
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B,RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST,ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
    [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN,MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX,QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
    
    %% Assignment fo the values given from the global variable UserObj
    tau = userObj.tau;
    sigmas = userObj.sigmas;
    nominal_model = userObj.nominal_model;
    a_limit = userObj.a_limit;
    z_analog = userObj.z_analog;
    z_digital = userObj.z_digital;
    omega_r_set = userObj.omega_r_set;
    omega_a_set = userObj.omega_a_set;
    omega_r_set_size = length(omega_r_set);
    omega_a_set_size = length(omega_a_set);
   
    %% The ind is composed by         |omega_minus||omega_plus||a vector|
    %% The indices corresponds to     |omega_r    ||omega_a   |
    omega_minus_index = omega_r_set(find(ind(1:omega_r_set_size) == 1));
    omega_plus_index = omega_a_set(find(ind(omega_r_set_size+1 : omega_r_set_size + omega_a_set_size) == 1));

    omega_minus = ind(1:omega_r_set_size);
    omega_plus = ind(omega_r_set_size+1 : omega_a_set_size+omega_r_set_size);

    display('Omega_Plus');
    display(sum(omega_plus));
    display('Omega_Minus');
    display(sum(omega_minus));

    %% Changing The status of the transmission lines
    
    z_digital(omega_plus_index) = 1;
    z_digital(omega_minus_index) = 0;

    Pg_star = userObj.Pg_star;
    lmp_star = userObj.lmp_star;
    PFrom_star = userObj.PFrom_star;
        
    
    %% The new model with the nominal topology is
    display('Z Digital attack model');
    display(z_digital);
    attack_model = utils().change_topology_model(nominal_model, z_digital);
    display(nominal_model.branch(:,BR_STATUS));
    display(attack_model.branch(:,BR_STATUS));

    %% Solving the SCED model under attack
    display('RUN OPF ON ATTACK TOPOLOGY');
    opf_attack_model_res = runopf(attack_model);
    if (opf_attack_model_res.success == 0)
        o = [+Inf; +Inf];
        return;
    end

    Pg_cross = opf_attack_model_res.gen(:,PG);
    lmp_cross = opf_attack_model_res.bus(:,LAM_P);
    PFrom_cross = opf_attack_model_res.branch(:,PF);

    %% Substitute the generation power and run pf problem on the nominal model
    nominal_model2 = nominal_model;
    nominal_model2.gen(:,PG) = Pg_cross;

    %% Displaying results
    display('RUN POWER FLOW EQUATION WITHNOMINAL TOPOLOGY');
    [MVAbase, result.bus, result.gen, result.branch, success, et] = runpf(nominal_model2);
    display(result);
    PFrom_circled_times = result.branch(:,PF);
    Pg_circled_times = result.gen(:,PG);
    display('Old_PG')
    display(PFrom_circled_times);
    display('New_PG');
    display(Pg_circled_times);

    %% Exposing the fitness value , whcih in this case is composed from a vector (fitness1 & fitness2) due to the fatc that we are solving
    %% a multi-object optimization problem 
    gencost = nominal_model.gencost;
    eta = ( (sum(totcost(gencost, Pg_circled_times)) - sum(totcost(gencost,Pg_star))) / sum(totcost(gencost,Pg_star) ) ) * 100;
    fitness1 = -(sum(totcost(gencost, Pg_circled_times)) - sum(totcost(gencost,Pg_star)));
    fitness2 = -(sum(lmp_cross-lmp_star)); 
    display(sum(totcost(gencost, Pg_circled_times)));
    display(sum(totcost(gencost,Pg_star)));
    display("ETA VALUE");
    display(eta);
    display("fitness1");
    display(fitness1);
    display("Lmp_star & Lmp_cross");
    display(lmp_star);
    display(lmp_cross);
    display("fitness2")
    display(fitness2);
    o = [fitness1; fitness2]; 
    
end


%% Function with reproduce the 2nd scenario for the first fitness with randomically the Case 1 or Case2 for the second fitness value
function o=MOEAD_SCENARIO_2(x)
    
    display(x);
    size(x)
    global userObj;
    ind = x;
    display('------------ FITNESS EVALUATION --------------');
    
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B,RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST,ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
    [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN,MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX,QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
    
    %% Assignment fo the values given from the global variable UserObj
    tau = userObj.tau;
    sigmas = userObj.sigmas;
    nominal_model = userObj.nominal_model;
    a_limit = userObj.a_limit;
    z_analog = userObj.z_analog;
    z_digital = userObj.z_digital;
    omega_r_set = userObj.omega_r_set;
    omega_a_set = userObj.omega_a_set;
    omega_r_set_size = length(omega_r_set);
    omega_a_set_size = length(omega_a_set);
   
    %% The ind is composed by         |omega_minus||omega_plus||a vector|
    %% The indices corresponds to     |omega_r    ||omega_a   |
    omega_minus_index = omega_r_set(find(ind(1:omega_r_set_size) == 1));
    omega_plus_index = omega_a_set(find(ind(omega_r_set_size+1 : omega_r_set_size + omega_a_set_size) == 1));

    omega_minus = ind(1:omega_r_set_size);
    omega_plus = ind(omega_r_set_size+1 : omega_a_set_size+omega_r_set_size);

    display('Omega_Plus');
    display(sum(omega_plus));
    display('Omega_Minus');
    display(sum(omega_minus));

    %% Changing The status of the transmission lines
    
    z_digital(omega_plus_index) = 1;
    z_digital(omega_minus_index) = 0;

    Pg_star = userObj.Pg_star;
    lmp_star = userObj.lmp_star;
    PFrom_star = userObj.PFrom_star;
        
    
    %% The new model with the nominal topology is
    display('Z Digital attack model');
    display(z_digital);
    attack_model = utils().change_topology_model(nominal_model, z_digital);
    display(nominal_model.branch(:,BR_STATUS));
    display(attack_model.branch(:,BR_STATUS));

    %% Solving the SCED model under attack
    display('RUN OPF ON ATTACK TOPOLOGY');
    opf_attack_model_res = runopf(attack_model);
    if (opf_attack_model_res.success == 0)
        o = [+Inf; +Inf];
        return;
    end

    Pg_cross = opf_attack_model_res.gen(:,PG);
    lmp_cross = opf_attack_model_res.bus(:,LAM_P);
    PFrom_cross = opf_attack_model_res.branch(:,PF);

    %% Substitute the generation power and run pf problem on the nominal model
    nominal_model2 = nominal_model;
    nominal_model2.gen(:,PG) = Pg_cross;

    %% Displaying results
    display('RUN POWER FLOW EQUATION WITHNOMINAL TOPOLOGY');
    [MVAbase, result.bus, result.gen, result.branch, success, et] = runpf(nominal_model2);
    display(result);
    PFrom_circled_times = result.branch(:,PF);
    Pg_circled_times = result.gen(:,PG);
    display('Old_PG')
    display(PFrom_circled_times);
    display('New_PG');
    display(Pg_circled_times);

    %% Exposing the fitness value , whcih in this case is composed from a vector (fitness1 & fitness2) due to the fatc that we are solving
    %% a multi-object optimization problem 
    gencost = nominal_model.gencost;
    eta = ( (sum(totcost(gencost, Pg_circled_times)) - sum(totcost(gencost,Pg_star))) / sum(totcost(gencost,Pg_star) ) ) * 100;
    fitness = norm(PFrom_circled_times - PFrom_star);
    fitness2 = -(sum(lmp_cross-lmp_star)); 
    display("ETA VALUE");
    display(eta);
    display("fitness");
    display(fitness);
    display("fitness2");
    display(fitness2);
    o = [fitness; fitness2]; 
    
end




