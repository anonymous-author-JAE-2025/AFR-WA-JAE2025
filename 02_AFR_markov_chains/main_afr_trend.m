% main_afr_calculation.m

%% Clean environment %%
clc; clear;

%% Load parameters %%
% Load demographic parameters
data = load('data/demographic_parameters_trend.mat');  % Changed to trend data
survival = data.survival; % [15 ages × 49 years × 100 iterations × 2 sexes]
breeding = data.breeding; % Same structure

%% Initialize arrays to store results %%
n_iter = 100;
n_years = 49;  % Changed from n_cohorts to n_years
afr_values = zeros(n_iter, n_years, 2); % [iterations × years × sex]
R = [16 17]; % Reproductive states

%% Process each sex, iteration, and year %%
for sex = 1:2
    fprintf('Processing sex %d/2\n', sex);
    for iter = 1:n_iter
        if mod(iter, 10) == 0
            fprintf('Iteration %d/100\n', iter);
        end
        for year = 1:n_years  % Changed from cohort to year
            % Extract survival and breeding probabilities for this combination
            surv_probs = squeeze(survival(:, year, iter, sex));
            breed_probs = squeeze(breeding(:, year, iter, sex));
            
            % Construct matrices
            U = construct_U_matrix(surv_probs, breed_probs);
            F = construct_F_matrix();
            
            % Calculate AFR using interval_stats
            out = interval_stats(U, F, R);
            
            % Store mean AFR value
            afr_values(iter, year, sex) = out.eta1cond(1); % Take AFR from starting at age 1
        end
    end
end

%% Prepare output for R %%
% For trend data
trend_output = [];
for sex = 1:2
    sex_label = {'F', 'M'};
    for year = 1:n_years
        iterations = (1:n_iter)';
        years = repmat(year, n_iter, 1);
        sexes = repmat(sex_label(sex), n_iter, 1);
        afr_vals = afr_values(:, year, sex);
        
        trend_output = [trend_output; table(iterations, years, sexes, afr_vals, ...
            'VariableNames', {'Iteration', 'Year', 'Sex', 'AFR'})];
    end
end

% Write to CSV
writetable(trend_output, 'outputs/afr_trend_results.csv');

%% Save MATLAB workspace %%
save('workspace_afr_trend.mat', 'afr_values');

%% Display summary statistics %%
fprintf('\nSummary of AFR values:\n');
fprintf('Females: Mean = %.2f, SD = %.2f\n', ...
    mean(afr_values(:,:,1), 'all', 'omitnan'), ...
    std(afr_values(:,:,1), [], 'all', 'omitnan'));
fprintf('Males: Mean = %.2f, SD = %.2f\n', ...
    mean(afr_values(:,:,2), 'all', 'omitnan'), ...
    std(afr_values(:,:,2), [], 'all', 'omitnan'));

%% PLOT %%
% Calculate means and CIs for each year and sex
female_means = mean(afr_values(:,:,1));
male_means = mean(afr_values(:,:,2));

% Calculate 95% CI
female_ci = zeros(n_years, 2);
male_ci = zeros(n_years, 2);

for year = 1:n_years
    % Female CIs
    female_vals = afr_values(:,year,1);
    ci = prctile(female_vals, [2.5 97.5]);
    female_ci(year,:) = ci;
    
    % Male CIs
    male_vals = afr_values(:,year,2);
    ci = prctile(male_vals, [2.5 97.5]);
    male_ci(year,:) = ci;
end

% Create figure
figure('Position', [100 100 1000 600])
hold on

% Plot female data
f1 = fill([1:n_years, fliplr(1:n_years)], ...
    [female_ci(:,1)', fliplr(female_ci(:,2)')], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(1:n_years, female_means, 'r-', 'LineWidth', 2);

% Plot male data
f2 = fill([1:n_years, fliplr(1:n_years)], ...
    [male_ci(:,1)', fliplr(male_ci(:,2)')], ...
    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(1:n_years, male_means, 'b-', 'LineWidth', 2);

% Customize plot
xlabel('Year')
ylabel('Age at First Reproduction')
title('Age at First Reproduction over Time by Sex')
legend([f1, f2], {'Females', 'Males'}, 'Location', 'best')
grid on
ylim([min([female_ci(:); male_ci(:)])-0.5 max([female_ci(:); male_ci(:)])+0.5])
xlim([1 n_years])

% Add axis customization
ax = gca;
ax.Box = 'on';
ax.FontSize = 12;