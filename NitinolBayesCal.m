% Calibration of Nitinol Superelastic Model Parameters: Data Analysis
%
% "A Probabilistic Approach with Built-in Uncertainty Quantification for 
% the Calibration of a Superelastic Constitutive Model from Full-field 
% Strain Data"
%
% https://github.com/confluentmedical/nitinol-bayes-cal
%
%    Copyright 2020 Harshad M. Paranjape
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.
%
%% Set appropriate process mode and run the script.
% 1 = Process data for a diamond tensile experiment and store QoI.
% 2 = Create a library of simulations.
% 3 = Batch-process simulation results and store QoI.
% 4 = Least-squares calibration.
% 5 = Calibration using Bayesian Inference.
% 6 = Fit a surrogate model to simulation QoI using Machine Learning.
% 7 = Calibration using Bayesian Inference and Machine Leaning.
process_mode = 1;
%
if(process_mode == 1)
    %% Process data for a diamond tensile experiment and store QoI.
    % Inputs
    dic_strain_data_file     = 'data/QoI_expt_DIC_strains.mat'; % DIC data processed by NCORR. See notes below for explanation of DIC inputs.
    dic_strain_data_frames   = 1:22; % Frames at which DIC data is available.
    %
    coord_dic_apex{1}        = [-1.659  0.185]; % Coordinate of apex manually picked from DIC plots by initially setting coord_dic_apex to 0, 0.
    coord_dic_apex{2}        = [-1.659  0.1494];
    coord_dic_apex{3}        = [-1.663  0.1672];
    coord_dic_apex{4}        = [-1.651  0.0959];
    coord_dic_apex{5}        = [-1.623  0.0484];
    coord_dic_apex{6}        = [-1.6   -0.005];
    coord_dic_apex{7}        = [-1.576 -0.1059];
    coord_dic_apex{8}        = [-1.558 -0.1237];
    coord_dic_apex{9}        = [-1.546 -0.1237];
    coord_dic_apex{10}       = [-1.523 -0.2069];
    coord_dic_apex{11}       = [-1.505 -0.2425];
    coord_dic_apex{12}       = [-1.505 -0.2959];
    coord_dic_apex{13}       = [-1.523 -0.2247];
    coord_dic_apex{14}       = [-1.528 -0.2069];
    coord_dic_apex{15}       = [-1.54  -0.1416];
    coord_dic_apex{16}       = [-1.546 -0.1772];
    coord_dic_apex{17}       = [-1.558 -0.1416];
    coord_dic_apex{18}       = [-1.564 -0.1237];
    coord_dic_apex{19}       = [-1.570 -0.1594];
    coord_dic_apex{20}       = [-1.582 -0.07031];
    coord_dic_apex{21}       = [-1.588 -0.08813];
    coord_dic_apex{22}       = [-1.594 -0.04656];
    %
    dic_roi_size             = [3.52 4.75]; % Second dimension is up-down in the as-mounted diamond
    dic_subset_radius        = 0.093; % DIC subset radius in mm (= 15 px)
    search_radius_factor     = [2, 2.5, 3.5, 2.1];
    %
    % Points of interest for FEA-DIC comparison. Points picked in the
    % downsampled DIC plot at zero load. Four points in the order:
    % bottom-left (intrados), bottom-right inside grip), top-right (outside
    % grip), and top-left (extrados)
    poi    = [0.59 0.23; ...
              2.8  1.02; ...
              2.34 1.5; ...
              0.25 0.46];
    % Load-disp data details from Instron during the data set
    dic_load_disp_file       = 'data/QoI_expt_load.csv'; % Name of the load-disp data file from Instron during DIC test
    dic_load_disp_file_max1  = 528; % Row number in the file for max disp. Note: row 3 is the first data row
    dic_load_disp_file_max2  = 570; % Row number in the file for max disp. after the hold.
    dic_load_disp_file_min   = 633; % Row number in the file for min disp. during the subcycle
    %
    qoi = 5;
    %
    show_debugging_plots = 0;
    %% Notes:
    % i. The DIC strain data is obtained at a specific frames during the
    % tensile loading experiments. It's a good practice to obtain DIC data
    % at specific time intervals (e.g., 1 second) and then pick the frames
    % for extracting data. The simulation library for calibration is
    % processed in such a way that the quantities of interest (strain and
    % load) are processed at the same frames during the tensile test.
    % ii. The DIC images are processed using NCORR software in this
    % example. However, any other DIC software can be used. The results of
    % NCORR data processing are saved in the NCORR GUI using the File ->
    % Save Data menu option.
    % iii. In this example shear strain (e_xy) component is used for
    % calibration.
    % iv. The load displacement data file is obtained directly from the
    % Instron load output. In the Instron test method, be sure to export
    % the raw data. The file has three columns: time (s), extension (mm), 
    % load (N).
    %
    % Matrix for reduced data. Size: num_dic_frames x num_quantities_of_interest
    abq_ML_fitting_dic_val_matrix = zeros(numel(dic_strain_data_frames), qoi);
    strain_mode_dic = 'exy'; % DO NOT ALTER THIS. THIS VARIABLE STORES THE STRAIN MODE USED FOR FITTING. SAVE IT. USEFUL FOR SANITY CHECK LATER.
    %
    %% Load disp-load file from DIC experiemnt (Instron raw data)
    data_cell = importdata(dic_load_disp_file);
    data_dic_load = zeros(size(data_cell, 1), 3);
    for jj = 7:size(data_cell, 1)
        tmp = strsplit(erase(data_cell{jj, :}, '"'), ',');
        for kk = 1:3
            % FInally store three data columns. time, disp, load
            data_dic_load(jj, kk) = str2num(tmp{kk});
        end
    end
    %% Load DIC strain data from NCORR
    load(dic_strain_data_file);
    %% Loop over DIC frames of ineterst
    coord_dic_apex2 = zeros(1, 2);
    initial_grip_pos = zeros(1, 2);
    for xx = 1:numel(dic_strain_data_frames)
        % Visualize DIC data
        % NCORR reports Green-Lagrangian shear strain. http://www.ncorr.com/index.php/dic-algorithms#5
        exy_dic = exy{dic_strain_data_frames(xx)};
        exy_dic(exy_dic == 0) = NaN;
        [coord_1_dic, coord_2_dic] = meshgrid(1:size(exy_dic, 2), 1:size(exy_dic, 1));
        coord_1_dic = (coord_1_dic - max(coord_1_dic(:))/2) / max(coord_1_dic(:)) * dic_roi_size(1);
        coord_2_dic = (coord_2_dic - max(coord_2_dic(:))/2) / max(coord_2_dic(:)) * dic_roi_size(2);
        % Delete NaN data
        coord_1_dic(isnan(exy_dic(:))) = [];
        coord_2_dic(isnan(exy_dic(:))) = [];
        exy_dic(isnan(exy_dic(:))) = [];
        % Move apex to 0, 0
        coord_1_dic = coord_1_dic - coord_dic_apex{xx}(1);
        coord_2_dic = coord_2_dic - coord_dic_apex{xx}(2);
        %
        [~, idx2] = max(coord_2_dic(:));
        coord_dic_apex2 = [coord_1_dic(idx2) coord_2_dic(idx2)];
        % Save the y location of the grip in order to track the verical
        % displacement later
        if(xx == 1)
            initial_grip_pos = coord_dic_apex2;
        end
        %
        exy_dic_minmax = [min(exy_dic(:), max(exy_dic(:)))]; % Min/max strain values
        %
        if(show_debugging_plots)
            % Vidualize DIC data
            figure('Position', [100 100 600 400]);
            scatter(coord_1_dic(:), coord_2_dic(:), 36, exy_dic(:), 'filled');
            hold on
            scatter(coord_dic_apex{xx}(1), coord_dic_apex{xx}(2), 49, 'r');
            scatter(coord_dic_apex2(1), coord_dic_apex2(2), 49, 'r');
            xlim([0 3.5]);
            ylim([0 2.5]);
            axis square;
            colorbar;
            title('DIC data, FEA data as circles')
        end
        %
        poi_xx = [poi(:, 1) (poi(:, 2) + poi(:, 1) / initial_grip_pos(1) * (coord_dic_apex2(2) - initial_grip_pos(2)))];
        %
        poi_near_pts_1 = rangesearch([coord_1_dic(:), coord_2_dic(:)], poi_xx(1, :), search_radius_factor(1) * dic_subset_radius);
        poi_near_pts_2 = rangesearch([coord_1_dic(:), coord_2_dic(:)], poi_xx(2, :), search_radius_factor(2) * dic_subset_radius);
        poi_near_pts_3 = rangesearch([coord_1_dic(:), coord_2_dic(:)], poi_xx(3, :), search_radius_factor(3) * dic_subset_radius);
        poi_near_pts_4 = rangesearch([coord_1_dic(:), coord_2_dic(:)], poi_xx(4, :), search_radius_factor(4) * dic_subset_radius);
        %
        poi_near_pts = cat(1, poi_near_pts_1, poi_near_pts_2, poi_near_pts_3, poi_near_pts_4);
        mean_dic_strain_poi = zeros(4, 1); % Mean strain in DIC at points near each poi above
        for pp = 1:4
            mean_dic_strain_poi(pp) = mean(exy_dic(poi_near_pts{pp}), 'omitnan');
            if(show_debugging_plots)
                % Show points near poi over which mean strain is calculated
                scatter(coord_1_dic(poi_near_pts{pp}), coord_2_dic(poi_near_pts{pp}), 16, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
            end
        end
        %% Process load qoi
        if(show_debugging_plots)
            figure;
            hold on;
            p1 = plot(data_dic_load(:, 2), data_dic_load(:, 3), 'LineWidth', 2, 'Color', 'k');
            %
            xlim([0 0.85]);
            ylim([0 50]);
            xlabel('Displacement (mm)');
            ylabel('Load (N)');
            title('DIC load-disp curve');
            set(gca, 'FontSize', 20);
        end
        % Find load corresponding to current frame in DIC load-disp
        % data
        if(xx <= 11)
            curr_load_dic_idx = round(dic_load_disp_file_max1 / 10 * (xx - 1) + 1);
        else
            curr_load_dic_idx = round(dic_load_disp_file_max2 + (dic_load_disp_file_min - dic_load_disp_file_max2) / 10 * (xx - 12) + 1);
        end
        curr_load_dic = data_dic_load(curr_load_dic_idx, 3);
        if(show_debugging_plots)
            s1 = scatter(data_dic_load(curr_load_dic_idx, 2), curr_load_dic, 36, 'r', 'LineWidth', 2);
        end
        %% Store strain and load qoi
        abq_ML_fitting_dic_val_matrix(xx, 1:4) = mean_dic_strain_poi; % Strain at poi
        abq_ML_fitting_dic_val_matrix(xx, 5) = curr_load_dic; % Load
    end
    % Save data
    save('private/abq_ML_fitting_dic_val_matrix', 'abq_ML_fitting_dic_val_matrix', 'strain_mode_dic');
elseif(process_mode == 2)
    %% Create a library of simulations.
    abq_inp_template = 'SimLib_ftgdia-0p8_0p4.template'; % Name of Abaqus inp template file. This should have VAR1, VAR2 placeholders for substitution.
    % Ranges of the material parameters to vary in the simulations.
    var_Ea           = [40000 80000]; % Austenite stiffness MPa
    var_Em           = [20000 50000]; % Martensite stiffness MPa
    var_etr          = [0.035 0.055]; % Transformation strain
    var_UPS          = [300 500]; % UPS MPa
    var_LPS          = [100 300]; % LPS MPa
    var_CPS          = [350 700]; % Compression plateau MPa
    %
    num_samples      = 1000; % Number of simulation files to create using latin hypercube sampling
    % Run the function to create a library simulations from a template
    % file. See notes below.
    NitinolBayesCalSimLib(abq_inp_template, var_Ea, var_Em, var_etr, var_UPS, var_LPS, var_CPS, num_samples);
    %% Notes:
    % i. The Abaqus inp files are created in private/SimLib directory.
    % Copy these input files to an appropriate location and run the
    % simulations. 
    % ii. Batch-process the output of simulations using the
    % included Python post-processing scripts fem-get_*.py. Copy the data 
    % files generated by the post-processing scripts back into
    % private/SimLib directory.
elseif(process_mode == 3)
    %% Batch-process simulation results and store QoI.
    % Inputs
    num_samples              = 1000; % Number of test runs
    %
    fea_strain_data_files    = 'private/SimLib/s_e_SimLib_ftgdia-0p8_0p4-%04d_%02d.data';
    fea_strain_data_frames   = 0:21;
    fea_load_disp_file       = 'private/SimLib/rf_u_SimLib_ftgdia-0p8_0p4-%04d.data';
    show_debugging_plots     = 0;
    show_debugging_text      = 0;
    %
    dic_roi_size             = [3.52 4.75]; % Second dimension is up-down in the as-mounted diamond
    dic_subset_radius        = 0.093; % DIC subset radius in mm (= 15 px)
    search_radius_factor     = [2, 2.5, 3.5, 2.1];
    %
    % Points of interest for FEA-DIC comparison. Points picked in the
    % downsampled DIC plot at zero load. Four points in the order:
    % bottom-left (intrados), bottom-right inside grip), top-right (outside
    % grip), and top-left (extrados)
    poi                      = [0.59 0.23; ...
                                2.8  1.02; ...
                                2.34 1.5; ...
                                0.25 0.46];
    % Number of quantities of interest that will be used for error
    % calculation later = number of strain poi + load = 5.
    qoi                      = 5;
    % Matrix to reduced data. Size: num_samples x num_fea_frames x num_quantities_of_interest
    abq_ML_fitting_sim_val_matrix = zeros(num_samples, numel(fea_strain_data_frames), qoi);
    strain_mode_fea = 'exy'; % DO NOT ALTER THIS. THIS VARIABLE STORES THE STRAIN MODE USED FOR FITTING. SAVE IT. USEFUL FOR SANITY CHECK LATER.
    %%
    % Loop over simulation results
    for ss = 1:num_samples
        try
            coord_fea_apex2 = zeros(1, 2); % Position of the top-right corner of the FEA domain
            initial_grip_pos = zeros(1, 2); % Initial position of the top-right corner of the FEA domain
            %% Loop over FEA frames
            for xx = 1:numel(fea_strain_data_frames)
                %% Strain analysis
                data_fea = importdata(sprintf(fea_strain_data_files, ss, fea_strain_data_frames(xx)));
                % Visualize FEA data
                coord_1_fea = data_fea(:, 1);
                coord_2_fea = data_fea(:, 2);
                % Find the apex of fea data and move to 0, 0
                [~, idx] = min(coord_1_fea(:));
                coord_fea_apex = [coord_1_fea(idx) coord_2_fea(idx)];
                coord_1_fea = coord_1_fea - coord_fea_apex(1);
                coord_2_fea = coord_2_fea - coord_fea_apex(2);
                % Save top-right corner position
                if(xx == 1)
                    initial_grip_pos = [max(coord_1_fea) max(coord_2_fea)];
                    coord_fea_apex2 = [max(coord_1_fea) max(coord_2_fea)];
                else
                    coord_fea_apex2 = [max(coord_1_fea) max(coord_2_fea)];
                end
                % s_e file columns: Coord-X, coord-Y, e11, e22, e12, e_princ
                % Abaqus always reports engineering shear strain exy+eyx. Abaqus Analysis User's Guide 1.2.2 Conventions. Thus, divide by 2.
                exy_fea = data_fea(:, 5) / 2;
                exy_fea_smooth = scatstat2(coord_1_fea, coord_2_fea, exy_fea, dic_subset_radius);
                %
                if(show_debugging_plots)
                    subplot(1, 2, 1);
                    hold on
                    scatter(coord_1_fea(:), coord_2_fea(:), 4, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
                    xlim([0 3.5]);
                    ylim([0 2.5]);
                    axis square;
                    colorbar;
                    axis off;
                    set(gca, 'FontSize', 20);
                end
                % Separately plot smoothed FEA data
                if(show_debugging_plots)
                    s2 = subplot(1, 2, 2);
                    hold on
                    scatter(coord_1_fea(:), coord_2_fea(:), 64, exy_fea_smooth(:), 'filled');
                    xlim([0 3.5]);
                    ylim([0 2.5]);
                    axis square;
                    colorbar;
                    title('FEA data, smoothed');
                    axis off;
                    set(gca, 'FontSize', 20);
                end
                % Calculate mean strains at POIs
                poi_xx = [poi(:, 1) (poi(:, 2) + poi(:, 1) / initial_grip_pos(1) * (coord_fea_apex2(2) - initial_grip_pos(2)))];
                %
                poi_near_pts_1 = rangesearch([coord_1_fea(:), coord_2_fea(:)], poi_xx(1, :), search_radius_factor(1) * dic_subset_radius);
                poi_near_pts_2 = rangesearch([coord_1_fea(:), coord_2_fea(:)], poi_xx(2, :), search_radius_factor(2) * dic_subset_radius);
                poi_near_pts_3 = rangesearch([coord_1_fea(:), coord_2_fea(:)], poi_xx(3, :), search_radius_factor(3) * dic_subset_radius);
                poi_near_pts_4 = rangesearch([coord_1_fea(:), coord_2_fea(:)], poi_xx(4, :), search_radius_factor(4) * dic_subset_radius);
                %
                poi_near_pts = cat(1, poi_near_pts_1, poi_near_pts_2, poi_near_pts_3, poi_near_pts_4);
                mean_fea_strain_poi = zeros(4, 1); % Mean strain in FEA at points near each poi above
                for pp = 1:4
                    mean_fea_strain_poi(pp) = mean(exy_fea_smooth(poi_near_pts{pp}), 'omitnan');
                    if(show_debugging_plots)
                        % Show points near poi over which mean strain is calculated
                        scatter(coord_1_fea(poi_near_pts{pp}), coord_2_fea(poi_near_pts{pp}), 16, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
                    end
                end
                %
                %% Process load
                data_fea_load = importdata(sprintf(fea_load_disp_file, ss));
                curr_load_fea = 2 * data_fea_load(xx, 2); % Multiply FEA load by 2 because quarter symmetry is modeled.
                if(show_debugging_plots)
                    s1 = plot(data_fea_load(:, 1), data_fea_load(:, 2), 'k');
                    s2 = scatter(data_fea_load(xx, 1), data_fea_load(xx, 2), 49, 'b', 'LineWidth', 2, 'Marker', '+');
                end
                %% Store strain and load qoi
                abq_ML_fitting_sim_val_matrix(ss, xx, 1:4) = mean_fea_strain_poi; % Strain
                abq_ML_fitting_sim_val_matrix(ss, xx, 5) = curr_load_fea; % Load
            end
        catch
            warning('Skipping dataset %04d', ss)
        end
    end
    %
    save('private/abq_ML_fitting_sim_val_matrix', 'abq_ML_fitting_sim_val_matrix', 'strain_mode_fea');
elseif(process_mode == 4)
    %% Least-squares calibration
    dic_val_matrix    = 'private/abq_ML_fitting_dic_val_matrix.mat'; % Name of MAT file with experimental (DIC) QOI
    sim_val_matrix    = 'private/abq_ML_fitting_sim_val_matrix.mat'; % Name of MAT file with simulation QOI
    err_frame_weight  = [0.01 0.01 0.1 0.1 1.2 1.4 1.4 1.4 1.6 1.6 2 2 2 2 2 2 2 2 2 2 2 1.5]; % The weight for each POI.
    err_poi_weight    = [0.25 0.25 0.25 0.25 1.00]; % Weight for error at each poi. This is useful to balance influence of strain vs. load etc.
    % Load data
    load(dic_val_matrix);
    load(sim_val_matrix);
    %
    n_sim_samples = size(abq_ML_fitting_sim_val_matrix, 1);
    n_sim_qoi = size(abq_ML_fitting_sim_val_matrix, 3);
    %
    abq_sim_dic_error_matrix = zeros(n_sim_samples, n_sim_qoi); % Matrix to store errors.
    % Loop over all simulation samples
    for ii = 1:n_sim_samples
        sim_poi_matrix = squeeze(abq_ML_fitting_sim_val_matrix(ii, :, :));
        if(sum(sim_poi_matrix(:, end)) == 0)
            abq_sim_dic_error_matrix(ii, :) = nan;
        else
            % Loop over each poi set
            for jj = 1:size(sim_poi_matrix, 1)
                % Loop over qoi values in each set
                for kk = 1:size(sim_poi_matrix, 2)
                    err_kk = abs((abq_ML_fitting_dic_val_matrix(jj, kk) - sim_poi_matrix(jj, kk)) / abq_ML_fitting_dic_val_matrix(jj, kk) * err_frame_weight(jj));
                    if(isnan(err_kk))
                        abq_sim_dic_error_matrix(ii, kk) = abq_sim_dic_error_matrix(ii, kk) + 0.0;
                    else
                        abq_sim_dic_error_matrix(ii, kk) = abq_sim_dic_error_matrix(ii, kk) + err_kk;
                    end
                end
            end
        end
    end
    % Standardize error for each qoi
    for kk = 1:size(abq_ML_fitting_sim_val_matrix, 3)
        abq_sim_dic_error_matrix(:, kk) = (abq_sim_dic_error_matrix(:, kk) - mean(abq_sim_dic_error_matrix(:, kk), 'omitnan')) / std(abq_sim_dic_error_matrix(:, kk), 'omitnan');
    end
    %
    save('abq_sim_dic_error_matrix', 'abq_sim_dic_error_matrix', 'strain_mode_dic');
    %% Plot strain and load data for verification
    [~, idx] = min(sum(abq_sim_dic_error_matrix .* repmat(err_poi_weight, n_sim_samples, 1), 2));
    for jj = 1:5
        figure; hold on
        for ii = 1:1000
            plot(squeeze(abq_ML_fitting_sim_val_matrix(ii, :, jj)), 'Color', [0.8 0.8 0.8]);
        end
        plot(abq_ML_fitting_dic_val_matrix(:, jj), 'LineWidth', 5)
        plot(squeeze(abq_ML_fitting_sim_val_matrix(idx, :, jj)), 'LineWidth', 5, 'LineStyle', '--');
    end
elseif(process_mode == 5)
    %% Calibration using Bayesian Inference.
    % Run process_mode = 1, 2, and 3 first to have all the matrices in
    % place.
    load('private/abq_ML_fitting_dic_val_matrix.mat'); % DIC data
    load('private/abq_ML_fitting_sim_val_matrix.mat'); % Simulation data
    load('private/SimLib/step_1_2_ML_fitting/abq_ML_fitting_param_matrix.mat'); % Parameters
    % Trim the data by deleting simulations that were unsuccessful
    abq_ML_fitting_param_matrix(sum(squeeze(abq_ML_fitting_sim_val_matrix(:, :, 5)), 2) == 0, :) = [];
    abq_ML_fitting_sim_val_matrix(sum(squeeze(abq_ML_fitting_sim_val_matrix(:, :, 5)), 2) == 0, :, :) = [];
    % Define parameter limits for prior distribution generation
    var_Ea           = [40000 80000]; % Austenite stiffness (min max) VAR1
    var_Em           = [20000 50000]; % Martensite stiffness VAR2
    var_etr          = [0.035 0.055]; % Transformation strain VAR3
    var_UPS          = [300 500]; % UPS VAR4A
    var_LPS          = [100 300]; % LPS VAR5A
    var_CPS          = [350 700]; % Compression plateau VAR6
    err_poi_weight   = [0.25 0.25 0.25 0.25 1.00]; % Weight for error at each poi. This is useful to balance influence of strain vs. load etc.
    logpriorfun = @(m) (m(1)>var_Ea(1))&&(m(1)<var_Ea(2)) && (m(2)>var_Em(1))&&(m(2)<var_Em(2)) && (m(3)>var_etr(1))&&(m(3)<var_etr(2)) && (m(4)>var_UPS(1))&&(m(4)<var_UPS(2)) && (m(5)>var_LPS(1))&&(m(5)<var_LPS(2)) && (m(6)>var_CPS(1))&&(m(6)<var_CPS(2)) ;
    % Also define shorthand log likelyhood function
    loglikefun = @(m) loglikelyhoodfun(m, abq_ML_fitting_param_matrix, abq_ML_fitting_sim_val_matrix, abq_ML_fitting_dic_val_matrix, err_poi_weight);
    % Initialize walkers
    % m_maxlike = fminsearch(@(m)-loglikefun(m) , [mean(var_Ea) mean(var_Em) mean(var_etr) mean(var_UPS) mean(var_LPS) mean(var_CPS)]);
    n_walkers = 100;
    minit = bsxfun(@plus, repmat([mean(var_Ea) mean(var_Em) mean(var_etr) mean(var_UPS) mean(var_LPS) mean(var_CPS)], n_walkers, 1), 0.1*randn(n_walkers, 6).*repmat([(var_Ea(2) - var_Ea(1)) (var_Em(2) - var_Em(1)) (var_etr(2) - var_etr(1)) (var_UPS(2) - var_UPS(1)) (var_LPS(2) - var_LPS(1)) (var_CPS(2) - var_CPS(1))], n_walkers, 1));
    %
    tic
    m = gwmcmc(minit', {logpriorfun loglikefun}, 1000000, 'ThinChain', 10, 'burnin', 0.2);
    toc
    % Save results (nn stands for nearrest neighbor since we used nearest neighbors in the likelyhood function)
    m_nn = m;
    save('abq_ML_fitting_fitted_params_nn', 'm_nn')
    %
    figure
    [C,lags,ESS]=eacorr(m);
    plot(lags,C,'.-',lags([1 end]),[0 0],'k');
    grid on
    xlabel('lags')
    ylabel('autocorrelation');
    text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
    title('Markov Chain Auto Correlation')
    %
    figure
    ecornerplot(m,'ks',true,'color',[.6 .35 .3])
    %
    m = m(:,:)';
    format short g;
    disp('Result from BI (median)')
    disp(median(m, 1))
    % Median: 49858        38749     0.044961       399.54       197.04       525.71
    % 5% credible interval
    disp('Result from BI (95% credible interval)')
    disp(prctile(m, 2.5, 1))
    disp(prctile(m, 97.5, 1))
    % Result from BI (95% credible interval)
    %         49163        22055     0.035488       304.97       104.38       359.08 
    %         75016        40258     0.054486       495.14       294.75       691.23
    %
    disp('Result from least squares')
    disp([49395        39228      0.03739       340.86       204.49        434.5])
elseif(process_mode == 6)
    %% Fit a surrogate model to simulation QoI using Machine Learning.
    load('private/abq_ML_fitting_sim_val_matrix.mat'); % Simulation data matrix
    load('private/SimLib/abq_ML_fitting_param_matrix.mat'); % Simulation parameter matrix
    %
    n_rows = size(abq_ML_fitting_sim_val_matrix, 2);
    n_cols = size(abq_ML_fitting_sim_val_matrix, 3);
    sim_mdl_arr = cell(n_rows, n_cols); % Cell array for storing SVM models. One model for each qoi.
    %
    % Trim the data by deleting simulations that were unsuccessful
    abq_ML_fitting_param_matrix(sum(squeeze(abq_ML_fitting_sim_val_matrix(:, :, 5)), 2) == 0, :) = [];
    abq_ML_fitting_sim_val_matrix(sum(squeeze(abq_ML_fitting_sim_val_matrix(:, :, 5)), 2) == 0, :, :) = [];
    % Loop over qoi's
    for ii = 1:n_rows
        for jj = 1:n_cols
            sim_mdl_arr{ii, jj} = fitrsvm(abq_ML_fitting_param_matrix, squeeze(abq_ML_fitting_sim_val_matrix(:, ii, jj)), ...
                                          'KernelFunction', 'gaussian', 'KernelScale', 'auto', 'Standardize', true, ...
                                          'OptimizeHyperparameters', 'auto');
            close all;
        end
    end
    save('private/abq_ML_fitting_sim_SVM_models', 'sim_mdl_arr');
elseif(process_mode == 7)
    %% Calibration using Bayesian Inference and Machine Leaning.
    load('private/abq_ML_fitting_dic_val_matrix.mat'); % DIC data
    load('private/abq_ML_fitting_sim_val_matrix.mat'); % Simulation data
    load('private/SimLib/abq_ML_fitting_param_matrix.mat'); % Parameters
    load('private/abq_ML_fitting_sim_SVM_models.mat'); % Fitted SVM models for simulation data
    % Trim the data by deleting simulations that were unsuccessful
    abq_ML_fitting_param_matrix(sum(squeeze(abq_ML_fitting_sim_val_matrix(:, :, 5)), 2) == 0, :) = [];
    abq_ML_fitting_sim_val_matrix(sum(squeeze(abq_ML_fitting_sim_val_matrix(:, :, 5)), 2) == 0, :, :) = [];
    % Define parameter limits for prior distribution generation
    var_Ea           = [40000 80000]; % Austenite stiffness (min max) VAR1
    var_Em           = [20000 50000]; % Martensite stiffness VAR2
    var_etr          = [0.035 0.055]; % Transformation strain VAR3
    var_UPS          = [300 500]; % UPS VAR4A
    var_LPS          = [100 300]; % LPS VAR5A
    var_CPS          = [350 700]; % Compression plateau VAR6
    err_poi_weight   = [0.25 0.25 0.25 0.25 1.00]; % Weight for error at each poi. This is useful to balance influence of strain vs. load etc.
    logpriorfun = @(m) (m(1)>var_Ea(1))&&(m(1)<var_Ea(2)) && (m(2)>var_Em(1))&&(m(2)<var_Em(2)) && (m(3)>var_etr(1))&&(m(3)<var_etr(2)) && (m(4)>var_UPS(1))&&(m(4)<var_UPS(2)) && (m(5)>var_LPS(1))&&(m(5)<var_LPS(2)) && (m(6)>var_CPS(1))&&(m(6)<var_CPS(2)) ;
    % Also define shorthand log likelyhood function
    loglikefun = @(m) loglikelyhoodmlfun(m, sim_mdl_arr, abq_ML_fitting_dic_val_matrix, err_poi_weight);
    % Initialize walkers
    % m_maxlike = fminsearch(@(m)-loglikefun(m) , [mean(var_Ea) mean(var_Em) mean(var_etr) mean(var_UPS) mean(var_LPS) mean(var_CPS)]);
    n_walkers = 100;
    minit = bsxfun(@plus, repmat([mean(var_Ea) mean(var_Em) mean(var_etr) mean(var_UPS) mean(var_LPS) mean(var_CPS)], n_walkers, 1), 0.1*randn(n_walkers, 6).*repmat([(var_Ea(2) - var_Ea(1)) (var_Em(2) - var_Em(1)) (var_etr(2) - var_etr(1)) (var_UPS(2) - var_UPS(1)) (var_LPS(2) - var_LPS(1)) (var_CPS(2) - var_CPS(1))], n_walkers, 1));
    %
    tic
    m = gwmcmc(minit', {logpriorfun loglikefun}, 100000, 'ThinChain', 10, 'burnin', 0.2);
    toc
    % Save results (ml stands for machine learning since we used SVM in this method)
    m_ml = m;
    save('abq_ML_fitting_fitted_params_ml', 'm_ml')
    %
    figure
    [C,lags,ESS]=eacorr(m);
    plot(lags,C,'.-',lags([1 end]),[0 0],'k');
    grid on
    xlabel('lags')
    ylabel('autocorrelation');
    text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
    title('Markov Chain Auto Correlation')
    %
    figure
    ecornerplot(m,'ks',true,'color',[.6 .35 .3])
    %
    m = m(:,:)';
    format short g;
    disp('Result from BI (median)')
    disp(median(m, 1))
    % Median:         52578        42849     0.038979       353.22       203.62       417.07
    % Std:           6157.7         5632     0.003738        15.75       9.2001       18.581
    % Calculate MAP
    m_MAP = zeros(1, 6);
    for ii = 1:6
        [N, e] = histcounts(m(:, ii), 100);
        [~, idxh] = max(N);
        m_MAP(ii) = e(idxh);
    end
    disp('Results from BI (MAP)')
    disp(m_MAP)
    % MAP: 27955        32585     0.051261       342.21        55.44       481.62
    % 95% credible interval
    disp('Result from BI (95% credible interval)')
    disp(prctile(m, 2.5, 1))
    disp(prctile(m, 97.5, 1))
    %    43201        27785     0.035188       326.68       187.09          384
    %    66201        49610     0.047803       388.73       221.84       452.35
    %
    disp('Result from least squares')
    disp([49395        39228      0.03739       340.86       204.49        434.5])
end
