function p = loglikelyhoodfun(m, abq_ML_fitting_param_matrix, abq_ML_fitting_sim_val_matrix, abq_ML_fitting_dic_val_matrix, err_poi_weight)
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
    sim_param_nn = (abq_ML_fitting_param_matrix - repmat(m, size(abq_ML_fitting_param_matrix, 1), 1)).^2;
    sim_param_nn = (sim_param_nn - mean(sim_param_nn, 1, 'omitnan')) / std(sim_param_nn, 1, 1, 'omitnan');
    [~, idx] = min(sum(sim_param_nn, 2));
    % [~, idx] = min(sum((abq_ML_fitting_param_matrix - repmat(m, size(abq_ML_fitting_param_matrix, 1), 1)).^2, 2));
    % Obtain simulation values correspoding to idx
    sim_poi_matrix = squeeze(abq_ML_fitting_sim_val_matrix(idx, :, :));
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