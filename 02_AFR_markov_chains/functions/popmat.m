function out = popmat(theta)

% Build the population model for Wandering albatrosses - To estimate LRS and LEX

% First stage is fledgling, so PB1 = Fledged chicks of age 0

Max_age = 55;
max_age_class = 15; 

% Age first reproduction
age_first = 6; 
age_last = 15;

% offspring production: 
rho = 0.5;

% There are 6 states:
% State 1 = Pre-breeders (PB1-PB15); Age = 0-14 (14 stages)
n_PB = length(1:15);
% Stage 2 = Successful breeders (SB7); Age = 6-54, but group them from age 6+ (1 stage)
n_SB = length(1);
% Stage 3 = Failed breeders (FB7); Age = 6-54, but group them from age 6+ (1 stage)
n_FB = length(1);
% stage 4 = Post-successful breeders (PSB8); Age = 6-54, but group them from 6+ (1 stages); can only be reached from SB age 6 onward
n_PSB = length(1);
% stage 5 = Post-failed breeders (PFB8); Age = 6-54, but group them from 6+ (1 stages); can only be reached from FB age 6 onward
n_PFB = length(1);
% stage 6 = Non-breeders (NB9); Age = 6-54, but group them from 6+ (1 stages); can only be reached from PSB and PFB age 7 onward
n_NB = length(1);

% Pop matrix dimensions
s = n_PB + n_SB + n_FB + n_PSB + n_PFB + n_NB;

% Initiate the pop matrix A
A = zeros(s, s);
U = A; F = A;

% Parameters

%%%%%%%%%%%%%%%%%%%%%%%%
% Survival, sigma (phi in the NIMBLE model): 

% Pre-breeders
%   theta 1: survival probability PB1; Probability that a fledgling survives until age 1
%   theta 2: survival probability PB2 
%   theta 3: survival probability PB3
%   theta 4: survival probability PB4
%   theta 5: survival probability PB5
%   theta 6: survival probability PB6; survival probability of PB6 (age 5 years-old)
%   theta 7: survival probability PB7
%   theta 8: survival probability PB8
%   theta 9: survival probability PB9
%   theta 10: survival probability PB10
%   theta 11: survival probability PB11
%   theta 12: survival probability PB12
%   theta 13: survival probability PB13
%   theta 14: survival probability PB14
%   theta 15: survival probability PB15; survival probability of PB15 (age = 14)

% Successful breeders
%  theta 16: survival probability SB7; survival probability of SB 6 years old and older

% Failed breeders
%  theta 17: survival probability FB7; survival probability of FB 6 years old and older

% Post-Successful breeders
%  theta 18: survival probability PSB8; survival probability of PSB 7 years old and older

% Post-Failed breeders
%  theta 19: survival probability PFB8; survival probability of PFB 7 years old and older

% Non-breeders
%  theta 20: survival probability NB9; survival probability of NB 8 years old and older

%%%%%%%%%%%%%%%%%%%%%%%%
% Breeding probability, beta (psi in the NIMBLE model): 
% Pre-breeders
%  theta 21: breeding probability of PB1; probability that a PB1 (age 0) breeds from year t to t+1 
%  theta 22: breeding probability of PB2
%  theta 23: breeding probability of PB3
%  theta 24: breeding probability of PB4
%  theta 25: breeding probability of PB5
%  theta 26: breeding probability of PB6; probability that a PB6 (age 5) breeds from year t to t+1 
%  theta 27: breeding probability of PB7
%  theta 28: breeding probability of PB8
%  theta 29: breeding probability of PB9
%  theta 30: breeding probability of PB10
%  theta 31: breeding probability of PB11
%  theta 32: breeding probability of PB12
%  theta 33: breeding probability of PB13
%  theta 34: breeding probability of PB14
%  theta 35: breeding probability of PB15; probability that a PB15 (age 10) breeds from year t to t+1 

% Successful breeders
%  theta 36: breeding probability of SB7; probability that a SB 6 years old and older breeds from year t to t+1 

% Failed breeders
%  theta 37: breeding probability of FB7; probability that a FB 6 years old and older breeds from year t to t+1 

% Post-successful breeders
%  theta 38: breeding probability of PSB8; probability that a PSB age 7 and older breeds from year t to t+1 

% Post-failed breeders
%  theta 39: breeding probability of PFB8; probability that a PFB age 7 and older breeds from year t to t+1 

% Non-breeders
%  theta 40: breeding probability of NB9; probability that a NB age 8 and older breeds from year t to t+1 

%%%%%%%%%%%%%%%%%%%%%%%%
% Breeding success, gamma (rho in the NIMBLE model): 
% Pre-breeders
% theta 41: breeding success of PB1 (age 0); given having breed from age 0 to 1, probability that PB1 raises successfully a chick.
% theta 42: breeding success of PB2
% theta 43: breeding success of PB3
% theta 44: breeding success of PB4
% theta 45: breeding success of PB5
% theta 46: breeding success of PB6
% theta 47: breeding success of PB7
% theta 48: breeding success of PB8
% theta 49: breeding success of PB9
% theta 50: breeding success of PB10
% theta 51: breeding success of PB11
% theta 52: breeding success of PB12
% theta 53: breeding success of PB13
% theta 54: breeding success of PB14
% theta 55: breeding success of PB15

% Successful breeders
% theta 56: breeding success of SB7 (age 6); probability that a SB 6 years old and older breeds successfully 

% Failed breeders
% theta 57: breeding success of FB7 (age 6); probability that a FB 6 years old and older breeds successfully 

% Post successful breeders
% theta 58: breeding success of PSB8 (age 7); probability that a PSB 7 years old and older breeds successfully  

% Post failed breeders
% theta 59: breeding success of PFB8 (age 7); probability that a PFB 7 years old and older breeds successfully 

% Non breeders
% theta 60: breeding success of NB9 (age 8); probability that a NB 8 years old and older breeds successfully 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% TRANSITIONS
%%% Pre-breeders, PB, transitions 
A(2,1) = theta(1); % survival PB1 to PB2
A(3,2) = theta(2); % survival of PB2 to PB3
A(4,3) = theta(3); % survival of PB3 to PB4
A(5,4) = theta(4); % survival of PB4 to PB5
A(6,5) = theta(5); % survival of PB5 to PB6

%% Upon age 5, PB can transit to SB or FB or remain in PB
% Remaining a prebreeder, PB to PB
A(7,6) = theta(6)*(1-theta(26)); % survival of PB6 to PB7 * (1- probability of breeding as PB6 at age 5)
A(8,7) = theta(7)*(1-theta(27)); % survival of PB7 to PB8 * (1- probability of breeding at age 6)
A(9,8) = theta(8)*(1-theta(28)); % survival of PB8 to PB9 * (1- probability of breeding at age 7)
A(10,9) = theta(9)*(1-theta(29)); % survival of PB9 to PB10 * (1- probability of breeding at age 8)
A(11,10) = theta(10)*(1-theta(30)); % survival of PB10 to PB11 * (1- probability of breeding at age 9)
A(12,11) = theta(11)*(1-theta(31)); % survival of PB11 to PB12 * (1- probability of breeding at age 10)
A(13, 12) = theta(12)*(1-theta(32)); % survival of PB12 to PB13 * (1- probability of breeding at age 11)
A(14,13) = theta(13)*(1-theta(33)); % survival of PB13 to PB14 * (1- probability of breeding at age 12)
A(15, 14) = theta(14)*(1-theta(34)); % survival of PB14 to PB15 * (1- probability of breeding at age 13)
A(15,15) = theta(15)*(1-theta(35)); % survival of PB15 to PB15 * (1- probability of breeding at age 14)
% this last one could be removed

% Transit to successful breeder, PB to SB
A(16,6) = theta(6)*theta(26)*theta(46); % survival of PB6 * breeding probability of PB6 * Breeding success of PB6
A(16,7) = theta(7)*theta(27)*theta(47); % survival of PB7 * breeding probability of PB7 * Breeding success of PB7
A(16,8) = theta(8)*theta(28)*theta(48); % survival of PB8 * breeding probability of PB8 * Breeding success of PB8
A(16,9) = theta(9)*theta(29)*theta(49); % survival of PB9 * breeding probability of PB9 * Breeding success of PB9
A(16,10) = theta(10)*theta(30)*theta(50); % survival of PB10 * breeding probability of PB10 * Breeding success of PB10
A(16,11) = theta(11)*theta(31)*theta(51); % survival of PB11 * breeding probability of PB11 * Breeding success of PB11
A(16,12) = theta(12)*theta(32)*theta(52); % survival of PB12 * breeding probability of PB12 * Breeding success of PB12
A(16,13) = theta(13)*theta(33)*theta(53); % survival of PB13 * breeding probability of PB13 * Breeding success of PB13
A(16,14) = theta(14)*theta(34)*theta(54); % survival of PB14 * breeding probability of PB14 * Breeding success of PB14
A(16,15) = theta(15)*theta(35)*theta(55); % survival of PB15 * breeding probability of PB15 * Breeding success of PB15

% Transit to failed breeder, PB to FB
A(17,6) = theta(6)*theta(26)*(1-theta(46)); % survival of PB6 * breeding probability of PB6 * 1- Breeding success of PB6
A(17,7) = theta(7)*theta(27)*(1-theta(47)); % survival of PB7 * breeding probability of PB7 * 1- Breeding success of PB7
A(17,8) = theta(8)*theta(28)*(1-theta(48)); % survival of PB8 * breeding probability of PB8 * 1- Breeding success of PB8
A(17,9) = theta(9)*theta(29)*(1-theta(49)); % survival of PB9 * breeding probability of PB9 * 1- Breeding success of PB9
A(17,10) = theta(10)*theta(30)*(1-theta(50)); % survival of PB10 * breeding probability of PB10 * 1- Breeding success of PB10
A(17,11) = theta(11)*theta(31)*(1-theta(51)); % survival of PB11 * breeding probability of PB11 * 1- Breeding success of PB11
A(17,12) = theta(12)*theta(32)*(1-theta(52)); % survival of PB12 * breeding probability of PB12 * 1- Breeding success of PB12
A(17,13) = theta(13)*theta(33)*(1-theta(53)); % survival of PB13 * breeding probability of PB13 * 1- Breeding success of PB13
A(17,14) = theta(14)*theta(34)*(1-theta(54)); % survival of PB14 * breeding probability of PB14 * 1- Breeding success of PB14
A(17,15) = theta(15)*theta(35)*(1-theta(55)); % survival of PB15 * breeding probability of PB15 * 1- Breeding success of PB15


%% Successful breeders, SB, transitions
% Remaining successful breeder, SB to SB
A(16,16) = theta(16)*theta(36)*theta(56); % survival of SB7 * breeding probability of SB7 * breeding success of SB7

% Transit to failed breeder, SB to FB
A(17,16) = theta(16)*theta(36)*(1-theta(56)); % survival of SB7 * breeding probability of SB7 * 1- breeding success of SB7

% Transit to post-successful breeder, SB to PSB
A(18,16) = theta(16)*(1-theta(36)); % survival of SB7 * 1- breeding probability of SB7


%% Failed breeders, FB, transitions
% Remaining failed breeder, FB to FB
A(17,17) = theta(17)*theta(37)*(1-theta(57)); % survival of FB7 * breeding probability of FB7 * 1- breeding success of FB7

% Transit to successful breeder, FB to SB
A(16,17) = theta(17)*theta(37)*theta(57); % survival of FB7 * breeding probability of FB7 * breeding success of FB7

% Transit to post-failed breeder, FB to PFB
A(19,17) = theta(17)*(1-theta(37)); % survival of FB7 * 1 - breeding probability of FB7


%% Post-successful breeders, PSB, transitions
% Transit to successful breeder, PSB to SB
A(16,18) = theta(18)*theta(38)*theta(58); % survival of PSB7 * breeding probability of PSB7 * breeding success of PSB7; 

% Transit to failed breeder, PSB to FB
A(17,18) = theta(18)*theta(38)*(1-theta(58)); % survival of PSB8 * breeding probability of PSB8 * 1- breeding success of PSB8;

% Transit to non-breeder, PSB to NB
A(20,18) = theta(18)*(1-theta(38)); % survival of PSB8 * 1 - breeding probability of PSB8;


%% Post-failed breeders, PFB, transitions
% Transit to successful breeder, PFB to SB
A(16,19) = theta(19)*theta(39)*theta(59); % survival of PFB8 * breeding probability of PF8 * breeding success of PFB8

% Transit to failed breeder, PFB to FB
A(17,19) = theta(19)*theta(39)*(1-theta(59)); % survival of PFB8 * breeding probability of PF8 * 1- breeding success of PFB8

% Transit to non-breeder, PFB to NB
A(20,19) = theta(19)*(1-theta(39)); % survival of PFB8 * 1- breeding probability of PF8


%% Non-breeders, NB, transitions
% Transit to successful breeder, NB to SB
A(16,20) = theta(20)*theta(40)*theta(60); % survival of NB9 * breeding probability of NB9* breeding success of NB9

% Transit to failed breeder, NB to FB
A(17,20) = theta(20)*theta(40)*(1-theta(60)); % survival of NB9 * breeding probability of NB9 * 1- breeding success of NB9

% Remain a non-breeder, NB to NB
A(20,20) = theta(20)*(1-theta(40)); % survival of NB9 * 1- breeding probability of NB9


U = A;

%% Fertilities
% from pre-breeders
A(1,6) = theta(6)*theta(26)*theta(46)*rho; 
A(1,7) = theta(7)*theta(27)*theta(47)*rho; 
A(1,8) = theta(8)*theta(28)*theta(48)*rho; 
A(1,9) = theta(9)*theta(29)*theta(49)*rho; 
A(1,10) = theta(10)*theta(30)*theta(50)*rho; 
A(1,11) = theta(11)*theta(31)*theta(51)*rho; 
A(1,12) = theta(12)*theta(32)*theta(52)*rho;
A(1,13) = theta(13)*theta(33)*theta(53)*rho; 
A(1,14) = theta(14)*theta(34)*theta(54)*rho; 
A(1,15) = theta(15)*theta(35)*theta(55)*rho; 

% From successful breeders
A(1,16) = theta(16)*theta(36)*theta(56)*rho; 

% From failed breeders
A(1,17) = theta(17)*theta(37)*theta(57)*rho; 

% From Post-successful breeders
A(1,18) = theta(18)*theta(38)*theta(58)*rho;

% From Post-failed breeders
A(1,19) = theta(19)*theta(39)*theta(59)*rho; 

% From non-breeders
A(1,20) = theta(20)*theta(40)*theta(60)*rho;


F(1,:)=A(1,:);

% transition matrix

out.A=A;
out.F=F;
out.U=U;
out.theta=theta;

return

