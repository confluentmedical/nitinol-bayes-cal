function NitinolBayesCalSimLib(abq_inp_template, var_Ea, var_Em, var_etr, var_UPS, var_LPS, var_CPS, num_samples)
% Calibration of Nitinol Superelastic Model Parameters: Simulation Library
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

%% Inputs
% abq_inp_template: Name of Abaqus inp template file. This should have VAR1, VAR2 placeholders for substitution.
%
% var_Ea: Austenite stiffness (min max) VAR1
% var_Em: % Martensite stiffness VAR2
% var_etr: Transformation strain VAR3
% var_UPS: UPS VAR4A
% var_LPS: LPS VAR5A
% var_CPS: Compression plateau VAR6
%
% num_samples: Number of simulation files to create using latin hypercube sampling
%%
abq_inp_template_txt = fileread(abq_inp_template); % Read template file text
lhsamples = lhsdesign(num_samples, 6);
abq_ML_fitting_param_matrix = zeros(size(lhsamples));
[~, abq_inp_template_name, ~] = fileparts(abq_inp_template);
%
for ii = 1:num_samples
    % Make temp copy of inp file string
    abq_inp_template_txt_tmp = abq_inp_template_txt;
    % Define paarmeetr substitution strings
    var1_sample = var_Ea(1) + lhsamples(ii, 1) * (var_Ea(2) - var_Ea(1));
    var2_sample = var_Em(1) + lhsamples(ii, 2) * (var_Em(2) - var_Em(1));
    var3_sample = var_etr(1) + lhsamples(ii, 3) * (var_etr(2) - var_etr(1));
    var4A_sample = var_UPS(1) + lhsamples(ii, 4) * (var_UPS(2) - var_UPS(1));
    var4B_sample = var_UPS(1) + lhsamples(ii, 4) * (var_UPS(2) - var_UPS(1)) + 30;
    var5A_sample = var_LPS(1) + lhsamples(ii, 5) * (var_LPS(2) - var_LPS(1));
    var5B_sample = var_LPS(1) + lhsamples(ii, 5) * (var_LPS(2) - var_LPS(1)) - 30;
    var6_sample = var_CPS(1) + lhsamples(ii, 6) * (var_CPS(2) - var_CPS(1));
    %
    VAR1SUB = sprintf('%d', var1_sample);
    VAR2SUB = sprintf('%d', var2_sample);
    VAR3SUB = sprintf('%10.6f', var3_sample);
    VAR4ASUB = sprintf('%8.2f', var4A_sample);
    VAR4BSUB = sprintf('%8.2f', var4B_sample);
    VAR5ASUB = sprintf('%8.2f', var5A_sample);
    VAR5BSUB = sprintf('%8.2f', var5B_sample);
    VAR6SUB = sprintf('%8.2f', var6_sample);
    % Save parameters for later analysis
    abq_ML_fitting_param_matrix(ii, 1) = var1_sample;
    abq_ML_fitting_param_matrix(ii, 2) = var2_sample;
    abq_ML_fitting_param_matrix(ii, 3) = var3_sample;
    abq_ML_fitting_param_matrix(ii, 4) = var4A_sample;
    abq_ML_fitting_param_matrix(ii, 5) = var5A_sample;
    abq_ML_fitting_param_matrix(ii, 6) = var6_sample;
    % Replace substrings
    abq_inp_template_txt_tmp_rep = replace(abq_inp_template_txt_tmp, ...
                                           {'VAR1', 'VAR2', 'VAR3', 'VAR4A', 'VAR4B', 'VAR5A', 'VAR5B', 'VAR6'}, ...
                                           {VAR1SUB, VAR2SUB, VAR3SUB, VAR4ASUB, VAR4BSUB, VAR5ASUB, VAR5BSUB, VAR6SUB});
    % write inp file with name ftgdia_230120a-0p8_0p4-XXXX.inp IF
    % compression plateau is between UPS and 1.5 * UPS
    if(var6_sample > var4A_sample && var6_sample < 1.5*var4A_sample)
        abq_inp_outname = sprintf('%s-%04d.inp', abq_inp_template_name, ii);
        f = fopen(abq_inp_outname, 'w');
        fwrite(f, abq_inp_template_txt_tmp_rep);
        fclose(f);
    end
end
% Save parameter matrix
save('SimLib/abq_ML_fitting_param_matrix', 'abq_ML_fitting_param_matrix');
