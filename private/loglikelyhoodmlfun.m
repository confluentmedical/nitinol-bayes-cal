function p = loglikelyhoodmlfun(m, sim_mdl_arr, abq_ML_fitting_dic_val_matrix, err_poi_weight)
    %% Calculates log likelyhood of a parameter set m
    % m = input parameters, 
    % abq_ML_fitting_param_matrix = parameter matrix from simulation study (dictionary), 
    % abq_ML_fitting_sim_val_matrix = simulation value array from simulation study, 
    % abq_ML_fitting_dic_val_matrix) = DIC strain, load value array
    %
    % Find the parameters from abq_ML_fitting_param_matrix closest to m
    if(size(m, 1) > 1)
        m = reshape(m, 1, size(m, 1));
    end
    % Obtain simulation values from trained SVM model matrix
    n_rows = size(sim_mdl_arr, 1);
    n_cols = size(sim_mdl_arr, 2);
    sim_poi_matrix = zeros(n_rows, n_cols);
    for ii = 1:n_rows
        for jj = 1:n_cols
            sim_poi_matrix(ii, jj) = predict(sim_mdl_arr{ii, jj}, m);
        end
    end
    % Store total error between dic and simulation in abq_sim_dic_error 
    abq_sim_dic_error = 0;
    err_frame_weight  = [0.01 0.01 0.1 0.1 1.2 1.4 1.4 1.4 1.6 1.6 2 2 2 2 2 2 2 2 2 2 2 1.5]; % The weight for each POI.
    % Loop over each poi set
    for jj = 1:size(sim_poi_matrix, 1)
        % Loop over qoi values in each set
        for kk = 1:size(sim_poi_matrix, 2)
            err_kk = abs((abq_ML_fitting_dic_val_matrix(jj, kk) - sim_poi_matrix(jj, kk)) / abq_ML_fitting_dic_val_matrix(jj, kk) * err_frame_weight(jj) * err_poi_weight(kk))^2;
            if(isnan(err_kk))
                abq_sim_dic_error = abq_sim_dic_error + 0.0;
            else
                abq_sim_dic_error = abq_sim_dic_error + err_kk;
            end
        end
    end
    %
    abq_sim_dic_error = sqrt(abq_sim_dic_error);
    %
    p = -abq_sim_dic_error^2 / 2 / 0.01;
    if(isnan(p) || isinf(p))
        p = -100;
    end
end