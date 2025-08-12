% main_afr_calculation.m

%% Clean environment %%
clc; clear;

%% Load parameters %%
% Load demographic parameters
data = load('data/demographic_parameters_cohort.mat');
survival = data.survival;  % [15 ages × 49 cohorts × 100 iterations × 2 sexes]
breeding = data.breeding;  % Same structure

%% Initialize arrays to store results %%
n_iter = 100;
n_cohorts = 49;
afr_values = zeros(n_iter, n_cohorts, 2);  % [iterations × cohorts × sex]
R = [16 17];  % Reproductive states

%% Process each sex, iteration, and cohort %%
for sex = 1:2
    fprintf('Processing sex %d/2\n', sex);
    
    for iter = 1:n_iter
        if mod(iter, 10) == 0
            fprintf('Iteration %d/100\n', iter);
        end
        
        for cohort = 1:n_cohorts
            % Extract survival and breeding probabilities for this combination
            surv_probs = squeeze(survival(:, cohort, iter, sex));
            breed_probs = squeeze(breeding(:, cohort, iter, sex));
            
            % Construct matrices
            U = construct_U_matrix(surv_probs, breed_probs);
            F = construct_F_matrix();
            
            % Calculate AFR using interval_stats
            out = interval_stats(U, F, R);
            
            % Store mean AFR value
            afr_values(iter, cohort, sex) = out.eta1cond(1);  % Take AFR from starting at age 1
        end
    end
end

%% Prepare output for R %%
% For cohort data
cohort_output = [];
for sex = 1:2
    sex_label = {'F', 'M'};
    for cohort = 1:n_cohorts
        iterations = (1:n_iter)';
        cohorts = repmat(cohort, n_iter, 1);
        sexes = repmat(sex_label(sex), n_iter, 1);
        afr_vals = afr_values(:, cohort, sex);
        
        cohort_output = [cohort_output; table(iterations, cohorts, sexes, afr_vals, ...
            'VariableNames', {'Iteration', 'Cohort', 'Sex', 'AFR'})];
    end
end

% Write to CSV
writetable(cohort_output, 'outputs/afr_cohort_results.csv');

%% Save MATLAB workspace %%
save('workspace_afr_cohort.mat', 'afr_values');

%% Display summary statistics %%
fprintf('\nSummary of AFR values:\n');
fprintf('Females: Mean = %.2f, SD = %.2f\n', ...
    mean(afr_values(:,:,1), 'all', 'omitnan'), ...
    std(afr_values(:,:,1), [], 'all', 'omitnan'));
fprintf('Males: Mean = %.2f, SD = %.2f\n', ...
    mean(afr_values(:,:,2), 'all', 'omitnan'), ...
    std(afr_values(:,:,2), [], 'all', 'omitnan'));

%% PLOT %%

% Calculate means and CIs for each cohort and sex
n_cohorts = 49;
female_means = mean(afr_values(:,:,1));
male_means = mean(afr_values(:,:,2));

% Calculate 95% CI
female_ci = zeros(n_cohorts, 2);
male_ci = zeros(n_cohorts, 2);

for cohort = 1:n_cohorts
    % Female CIs
    female_vals = afr_values(:,cohort,1);
    ci = prctile(female_vals, [2.5 97.5]);
    female_ci(cohort,:) = ci;
    
    % Male CIs
    male_vals = afr_values(:,cohort,2);
    ci = prctile(male_vals, [2.5 97.5]);
    male_ci(cohort,:) = ci;
end

% Create figure
figure('Position', [100 100 1000 600])
hold on

% Plot female data
f1 = fill([1:n_cohorts, fliplr(1:n_cohorts)], ...
    [female_ci(:,1)', fliplr(female_ci(:,2)')], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(1:n_cohorts, female_means, 'r-', 'LineWidth', 2);

% Plot male data
f2 = fill([1:n_cohorts, fliplr(1:n_cohorts)], ...
    [male_ci(:,1)', fliplr(male_ci(:,2)')], ...
    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(1:n_cohorts, male_means, 'b-', 'LineWidth', 2);

% Customize plot
xlabel('Cohort')
ylabel('Age at First Reproduction')
title('Age at First Reproduction by Cohort and Sex')
legend([f1, f2], {'Females', 'Males'}, 'Location', 'best')
grid on
ylim([min([female_ci(:); male_ci(:)])-0.5 max([female_ci(:); male_ci(:)])+0.5])
xlim([1 n_cohorts])

% Add axis customization
ax = gca;
ax.Box = 'on';
ax.FontSize = 12;