% Initial parameters
% Some important info
pb = 16; % number of PB stages
repro_start_ageclass = 7; % First age class of adult reproductive states
max_ageclass = 7; 
adult_classes = max_ageclass - repro_start_ageclass +1; % number of adult age classes
nb_reprostates = 5; % number of adult reproductive states
tau = pb + nb_reprostates*adult_classes; % matrix dimensions
n_p = pb*3 + nb_reprostates*3*adult_classes; % number of parameters

% number of posteriors used
n_it = 1500;