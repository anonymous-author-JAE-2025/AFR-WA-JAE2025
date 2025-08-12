function U = construct_U_matrix(survival_probs, breeding_probs)
    n_states = 17;  % 15 age states + 2 breeding states
    U = zeros(n_states, n_states);
    
    % Survival transitions without breeding
    for age = 1:14  % Transitions up to age 14->15
        U(age+1, age) = survival_probs(age) * (1 - breeding_probs(age));
    end
    
    % Age 15 can stay in age 15 if survives and doesn't breed
    U(15, 15) = survival_probs(15) * (1 - breeding_probs(15));
    
    % Add transitions to breeding states from all living ages
    for age = 1:15
        U(16, age) = survival_probs(age) * breeding_probs(age) * 0.5;  % First breeding state
        U(17, age) = survival_probs(age) * breeding_probs(age) * 0.5;  % Second breeding state
    end
    
    % Make breeding states absorbing (self-loops with probability 1)
    U(16,16) = 1;
    U(17,17) = 1;
end
