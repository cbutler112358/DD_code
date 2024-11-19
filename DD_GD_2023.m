%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program for studying the interaction of density dependence and 
% population suppression. This is a newer script, created on 07/26/23, so
% consolidate all existing code and update things where appropriate.
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Capital letters in genotypes are wild-type, while lower-case denote
% transgenes. Pre-release equilibria are equivalent for all DD cases. 
% 
% Infected populations are typically ignored (ICs set equal to 0) since the
% focus is more on population dynamics, and epi performance can be inferred 
% by looking at R0. 
% 
% Note that a wild-type-only population can be simulated by using a 0
% release ratio for the case study simulations. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sandbox for general genetic load; individual replicate
% simulation parameters
params          = struct();
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.beta     = 0.5; 
tVec            = 0:(100*365); 

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
N               = 2*10^4;

params.lethality_type   = "LA";
params.lethality_case   = "BSL";

opts = odeset('RelTol',1e-8,'AbsTol',1e-9);

params.g        = 0.95;
% initial conditions (pre-release equil.)
WT_I0 = [(muY/m)*N, (muY/m)*N, (muY/muZ)*N, N];

[~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
% did the pop go extinct?
% plot(sum(xout,2), 'linewidth',1.5);
sum(xout(end,:))

%% sandbox for 2L UD; individual replicate
close all

% set parameters
params          = struct();
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.s_t      = 1.0;
N               = 2*(10^4); 

params.lethality_case   = "FSL";
params.lethality_type   = "LA";
params.supp             = "SS";


releaseType             = "MOR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

% [AABB_larvae_male, AABB_larvae_female    2
%  AABb_larvae_male, AABb_larvae_female    4
%  AAbb_larvae_male, AAbb_larvae_female    6
%  AaBB_larvae_male, AaBB_larvae_female    8
%  AaBb_larvae_male, AaBb_larvae_female    10
%  Aabb_larvae_male, Aabb_larvae_female    12
%  aaBB_larvae_male, aaBB_larvae_female    14
%  aaBb_larvae_male, aaBb_larvae_female    16
%  aabb_larvae_male, aabb_larvae_female    18
%  AABB_adult_males, AABb_adult_males      20
%  AAbb_adult_males, AaBB_adult_males      22
%  AaBb_adult_males, Aabb_adult_males      24
%  aaBB_adult_males, aaBb_adult_males      26
%  aabb_adult_males, AABB_adult_females    28
%  AABb_adult_females, AAbb_adult_females  30
%  AaBB_adult_females, AaBb_adult_females  32
%  Aabb_adult_females, aaBB_adult_females  34
%  aaBb_adult_females, aabb_adult_females] 36

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
beta            = 0.5;
alpha           = 0.0321302;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

wt_I0                   = zeros(1,36);    
wt_I0(28)               = N;
wt_I0(1)                = (1/2)*fInv_LAMBDA;
wt_I0(2)                = (1/2)*fInv_LAMBDA; 
wt_I0(19)               = (m/(2*muZ))*fInv_LAMBDA;

% release transgenic males
release_I0              = wt_I0;

tmax = 5*365;
tVec = 0:1:tmax; 
% study release thresholds up to 12
rFlag = 100;
%%% counter = 1;

% set s_a
params.s_a = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% GEN LOGISTIC (BETA = 0.5) %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.beta             = 0.9;
params.alpha            = 0.014543805956894;

% adjusted release based on releaseType
if releaseType == "MOR"
    release_I0(27) = rFlag*wt_I0(19);
elseif releaseType == "BSR"
    release_I0(27) = rFlag*wt_I0(19);
    release_I0(36) = rFlag*wt_I0(28);
elseif releaseType == "FOR"
    release_I0(36) = rFlag*wt_I0(28);
end

[~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts);

% relative female pop
femalePop = sum(xout(:,28:end),2)/(2*(10^4));
plot(femalePop)

% sum(xout(end,28:end))/(2*(10^4));

% check allelic frequency
ind1 = [3,4,7,8,20,22,29,31];
    % one copy of a or b
ind2 = [5,6,9,10,13,14,21,23,25,30,32,34];
    % two copies of a or b
ind3 = [11,12,15,16,24,26,33,35];
    % three copies of a or b
ind4 = [17,18,27,36];
    % four copies of a or b      

% is it above frequency of 5%?
fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
    4*sum(xout(end,ind4)))/(4*sum(xout(end,:)))
% is it unchanging?
freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
    4*sum(xout(:,ind4),2))./(4*sum(xout,2));
freq_steady = var(freq_steady((end-99):end)) < 0.01
thresh_invade = 0.05; 

%% sandbox for 2L UD

% set parameters
params          = struct();
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.s_t      = 1.0;
N               = 2*(10^4); 

params.lethality_case   = "BSL";
params.lethality_type   = "LA";
params.supp             = "SS";


releaseType             = "MOR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

% study invasion over range of fitness costs
s_aVec                  = 0.07; % 0:0.01:0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
beta            = 0.5;
alpha           = 0.0321302;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

% [AABB_larvae_male, AABB_larvae_female    2
%  AABb_larvae_male, AABb_larvae_female    4
%  AAbb_larvae_male, AAbb_larvae_female    6
%  AaBB_larvae_male, AaBB_larvae_female    8
%  AaBb_larvae_male, AaBb_larvae_female    10
%  Aabb_larvae_male, Aabb_larvae_female    12
%  aaBB_larvae_male, aaBB_larvae_female    14
%  aaBb_larvae_male, aaBb_larvae_female    16
%  aabb_larvae_male, aabb_larvae_female    18
%  AABB_adult_males, AABb_adult_males      20
%  AAbb_adult_males, AaBB_adult_males      22
%  AaBb_adult_males, Aabb_adult_males      24
%  aaBB_adult_males, aaBb_adult_males      26
%  aabb_adult_males, AABB_adult_females    28
%  AABb_adult_females, AAbb_adult_females  30
%  AaBB_adult_females, AaBb_adult_females  32
%  Aabb_adult_females, aaBB_adult_females  34
%  aaBb_adult_females, aabb_adult_females] 36

wt_I0                   = zeros(1,36);    
wt_I0(28)               = N;
wt_I0(1)                = (1/2)*fInv_LAMBDA;
wt_I0(2)                = (1/2)*fInv_LAMBDA; 
wt_I0(19)               = (m/(2*muZ))*fInv_LAMBDA;

% release transgenic males
release_I0              = wt_I0;

% matrices for storing everything
invasionThresh          = zeros(4,length(s_aVec)); 
controlEquil            = zeros(4,length(s_aVec)); 

tmax = 10*365;
tVec = 0:1:tmax; 
% study release thresholds up to 12
rFlagMax = 12;     
%%% counter = 1;

for i = 1:length(s_aVec) 
    % loading message
    fprintf("Running sim %.0f of %.0f: s = %0.2f.\n", i, length(s_aVec), ...
        s_aVec(i));

    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;
    %%% releaseMult = 1.5; 

    % set s_a
    params.s_a = s_aVec(i);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 0.5) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    params.beta             = 0.5;
    params.alpha            = 0.0321302;

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 

    % [AABB_larvae_male, AABB_larvae_female    2
    %  AABb_larvae_male, AABb_larvae_female    4
    %  AAbb_larvae_male, AAbb_larvae_female    6
    %  AaBB_larvae_male, AaBB_larvae_female    8
    %  AaBb_larvae_male, AaBb_larvae_female    10
    %  Aabb_larvae_male, Aabb_larvae_female    12
    %  aaBB_larvae_male, aaBB_larvae_female    14
    %  aaBb_larvae_male, aaBb_larvae_female    16
    %  aabb_larvae_male, aabb_larvae_female    18
    %  AABB_adult_males, AABb_adult_males      20
    %  AAbb_adult_males, AaBB_adult_males      22
    %  AaBb_adult_males, Aabb_adult_males      24
    %  aaBB_adult_males, aaBb_adult_males      26
    %  aabb_adult_males, AABB_adult_females    28
    %  AABb_adult_females, AAbb_adult_females  30
    %  AaBB_adult_females, AaBb_adult_females  32
    %  Aabb_adult_females, aaBB_adult_females  34
    %  aaBb_adult_females, aabb_adult_females] 36  

    % check allelic frequency
    ind1 = [3,4,7,8,20,22,29,31];
        % one copy of a or b
    ind2 = [5,6,9,10,13,14,21,23,25,30,32,34];
        % two copies of a or b
    ind3 = [11,12,15,16,24,26,33,35];
        % three copies of a or b
    ind4 = [17,18,27,36];
        % four copies of a or b        

    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
        4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
        4*sum(xout(:,ind4),2))./(4*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    %%%
    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
            4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
            4*sum(xout(:,ind4),2))./(4*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(27) = midFlag*wt_I0(19);
            elseif releaseType == "BSR"
                release_I0(27) = midFlag*wt_I0(19);
                release_I0(36) = midFlag*wt_I0(28);
            elseif releaseType == "FOR"
                release_I0(36) = midFlag*wt_I0(28);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                4*sum(xout(:,ind4),2))./(4*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                    4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                    4*sum(xout(:,ind4),2))./(4*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(1,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(1,i)   = sum(xout(end,28:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(1,i)   = nan; 
    end

end

%% Does early-acting mortality increase pop size? 
% Exactly what the posed question asks; Does an early-acting bi-sex lethal
% genetic load increase mosquito pop size? 
% Illustrates Hydra effect

% simulation parameters
params          = struct();
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
tVec            = 0:365; 

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
N               = 2*10^4;

params.lethality_type   = "EA";
params.lethality_case   = "FSL";

opts = odeset('RelTol',1e-8,'AbsTol',1e-9);

gVec = 0.0:0.05:1.0; 
% initial conditions (pre-release equil.)
WT_I0 = [(muY/m)*N, (muY/m)*N, (muY/muZ)*N, N];

params.lethality_case   = "BSL";

% vectors
genLog1_num = zeros(1,length(gVec));
genLog2_num = zeros(1,length(gVec));
genLog3_num = zeros(1,length(gVec));
log1_num    = zeros(1,length(gVec));

for i = 1:length(gVec)
    disp(i);
    params.g        = gVec(i);

    params.beta     = 0.5;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
    % store adult female population
    genLog1_num(i)  = sum(xout(end,3:4));

    params.beta     = 1.0;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
    genLog2_num(i)  = sum(xout(end,3:4));

    params.beta     = 1.5;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
    genLog3_num(i)  = sum(xout(end,3:4));

    params.beta     = 0.9;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_iii(t,x,params),tVec,WT_I0,opts); 
    log1_num(i)     = sum(xout(end,3:4));

    % [female_juveniles, male_juveniles, ...
    % [female_adults, male_adults]
end

% for a 1x2 plot...
subplot(1,2,1)

% plot(gVec, genLog1_an, '-k','linewidth',1.5)
% hold on
plot(gVec, genLog1_num/((1 + muY/muZ)*N), '-k','linewidth',1.5)
title('(a)')
hold on
% plot(gVec, genLog2_an, '-k','linewidth',1.5)
plot(gVec, genLog2_num/((1 + muY/muZ)*N), '--k','linewidth',1.5)
% plot(gVec, genLog3_an, '-k','linewidth',1.5)
plot(gVec, genLog3_num/((1 + muY/muZ)*N), '-.k','linewidth',1.5)
% plot(gVec, log1_an, '-k','linewidth',1.5)
plot(gVec, log1_num/((1 + muY/muZ)*N), ':k','linewidth',1.5)
% title('FSL','interpreter','latex')
xlabel('genetic load, $g$','interpreter','latex')
ylabel('relative total pop.','interpreter','latex')
set(gca,'fontsize',15)
ylim([0,1])


%% sandbox for 1L UD
close all 

% set parameters
params          = struct();
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.s_t      = 1;
N               = 2*(10^4); 

params.lethality_case   = "BSL";
params.lethality_type   = "LA";


releaseType             = "BSR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-10,'AbsTol',1e-10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
beta            = 0.5;
alpha           = 0.0321302;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

% [AA_larvae_males, Aa_larvae_males,        2
%  Ab_larvae_males, ab_larvae_males,        4
%  aa_larvae_males, bb_larvae_males,        6
%  AA_larvae_females, Aa_larvae_females,    8
%  Ab_larvae_females, ab_larvae_females,    10
%  aa_larvae_females, bb_larvae_females,    12
%  AA_adults_males, Aa_adults_males,        14
%  Ab_adults_males, ab_adults_males,        16
%  aa_adults_males, bb_adults_males,        18
%  AA_adults_females, Aa_adults_females,    20
%  Ab_adults_females, ab_adults_females,    22
%  aa_adults_females, bb_adults_females]    24

wt_I0                   = zeros(1,24);    
wt_I0(19)               = N;
wt_I0(7)                = (1/2)*fInv_LAMBDA;
wt_I0(1)                = (1/2)*fInv_LAMBDA; 
wt_I0(13)               = (m/(2*muZ))*fInv_LAMBDA;
% release transgenic males
release_I0              = wt_I0;

tmax = 10*365;
tVec = 0:1:tmax; 
% study release thresholds up to 12
rFlagMax = 12;
%%% counter = 1;

% left flag, right flag
lFlag = 0.01;
rFlag = rFlagMax;
%%% releaseMult = 1.5; 

% set s_a
params.s_a = 0.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% GEN LOGISTIC (BETA = 0.5) %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.beta             = 0.5;
params.alpha            = 0.0321302;

% adjusted release based on releaseType
if releaseType == "MOR"
    release_I0(16) = rFlag*wt_I0(13);
elseif releaseType == "BSR"
    release_I0(16) = rFlag*wt_I0(13);
    release_I0(22) = rFlag*wt_I0(19);
elseif releaseType == "FOR"
    release_I0(22) = rFlag*wt_I0(19);
end

[~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts);

% check allelic frequency
ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
% is it above frequency of 5%?
fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
% is it unchanging?
freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
freq_steady = var(freq_steady((end-99):end)) < 0.01;
thresh_invade = 0.05; 

if ((fixProp_tmp > thresh_invade) && freq_steady) 
    while ~(abs(lFlag - rFlag) < 10^(-2))        
        midFlag = (lFlag + rFlag)/2;
        disp(midFlag)
        % adjusted release based on releaseType
        if releaseType == "MOR"
            release_I0(16) = midFlag*wt_I0(13);
        elseif releaseType == "BSR"
            release_I0(16) = midFlag*wt_I0(13);
            release_I0(22) = midFlag*wt_I0(19);
        elseif releaseType == "FOR"
            release_I0(22) = midFlag*wt_I0(19);
        end
        
        [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 
        
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
        fprintf("Freq: %.2f\n", fixProp_tmp);
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = 2;
        while (~freq_steady)
            disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
            tVecTmp = 0:1:(increaseCount*tmax); 
            
            [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = increaseCount + 1;
        end
    
        % boolean variable -> 1 if "fixed," 0 otherwise
        if ((fixProp_tmp > thresh_invade) && freq_steady) 
            tmpBool = 1;
        else
            tmpBool = 0;
        end
        % no fixation, right flag is right flag
        % fixation, right flag is mid flag
        rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
        % no fixation, left flag is mid flag
        % fixation, left flag is left flag
        lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
    end % end of while flag loop
else
    rFlag = NaN;
end

% run a sim with rFlag
if releaseType == "MOR"
    release_I0(16) = rFlag*wt_I0(13);
elseif releaseType == "BSR"
    release_I0(16) = rFlag*wt_I0(13);
    release_I0(22) = rFlag*wt_I0(19);
elseif releaseType == "FOR"
    release_I0(22) = rFlag*wt_I0(19);
end

[~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts);
fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));

plot(tVec, xout);
fprintf("Release at ratio %.2f produces allele frequency %.2f\n", rFlag, fixProp_tmp);
fprintf("Control equil relative to pre-release: %.2f\n", sum(xout(end,19:end))/(2*(10^4)));

%% sandbox for 1L UD individual replicates
% set parameters
params          = struct();
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.s_t      = 0.9;
N               = 2*(10^4); 

params.lethality_case   = "BSL";
params.lethality_type   = "LA";


releaseType             = "MOR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-12,'AbsTol',1e-12,'NonNegative',1:24);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
beta            = 0.5;
alpha           = 0.0321302;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

wt_I0                   = zeros(1,24);    
wt_I0(19)               = N;
wt_I0(7)                = (1/2)*fInv_LAMBDA;
wt_I0(1)                = (1/2)*fInv_LAMBDA; 
wt_I0(13)               = (m/(2*muZ))*fInv_LAMBDA;
% release transgenic males
release_I0              = wt_I0;

tmax = 100; % 1*365;
tVec = 0:1:tmax; 
% study release thresholds up to 12
rFlagMax = 20;     

% left flag, right flag
lFlag = 0.01;
rFlag = rFlagMax;
%%% releaseMult = 1.5; 

% set s_a
params.s_a = 0.05;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% GEN LOGISTIC (BETA = 0.5) %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params.beta             = 0.5;
% params.alpha            = 0.0321302;

params.beta             = 0.9;
params.alpha            = 0.014543805956894; 

% adjusted release based on releaseType
if releaseType == "MOR"
    release_I0(16) = rFlag*wt_I0(13);
elseif releaseType == "BSR"
    release_I0(16) = rFlag*wt_I0(13);
    release_I0(22) = rFlag*wt_I0(19);
elseif releaseType == "FOR"
    release_I0(22) = rFlag*wt_I0(19);
end

% [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts);
[~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts); 


% [AA_larvae_males, Aa_larvae_males,        2
%  Ab_larvae_males, ab_larvae_males,        4
%  aa_larvae_males, bb_larvae_males,        6
%  AA_larvae_females, Aa_larvae_females,    8
%  Ab_larvae_females, ab_larvae_females,    10
%  aa_larvae_females, bb_larvae_females,    12
%  AA_adults_males, Aa_adults_males,        14
%  Ab_adults_males, ab_adults_males,        16
%  aa_adults_males, bb_adults_males,        18
%  AA_adults_females, Aa_adults_females,    20
%  Ab_adults_females, ab_adults_females,    22
%  aa_adults_females, bb_adults_females]    24

ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
% transgene frequency
allele_freq = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));

% plot(xout(:,16))
% plot(xout(:,19)/N)
sum(xout(end,:))

% plot(xout);
plot(allele_freq,'linewidth',1.5)
ylim([0,1])

%% sandbox for Khamis version of 1LUD, individual replicate

release         = 0.95;         
tVec            = 0:365;        

%%%%%%%%%%%%%%%%%%%%%
%%%% HOMOLOGOUS %%%%%
%%%%%%%%%%%%%%%%%%%%%

params                  = struct();
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 16;
params.muX      = 0.03;
params.muZ      = log(10/9);
params.muY      = log(10/9);
params.m        = 0.1;
params.s_t      = 0.9;
params.s_a      = 0.05;
N               = 2*(10^3); 

params.lethality_type   = "EA";
params.supp             = "SS";

releaseType             = "MOR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

beta            = 0.9;
alpha           = (m/(2*2*1000*muZ))*((exp(lambda*m/(2*muZ)-m-muX) - 1)^(1/beta));
fInv            = @(x) (1/alpha)*((exp(x)-1).^(1/beta));
fInv_LAMBDA     = fInv(LAMBDA); 

params.beta     = beta;
params.alpha    = alpha;                

wt_I0                   = zeros(1,24);    
wt_I0(19)               = N;
wt_I0(7)                = (1/2)*fInv_LAMBDA;
wt_I0(1)                = (1/2)*fInv_LAMBDA; 
wt_I0(13)               = (m/(2*muZ))*fInv_LAMBDA;

release_I0              = wt_I0;
release_ratio           = release/(1-release);
disp(release_ratio);

if releaseType == "MOR"
    release_I0(16) = release_ratio*wt_I0(13);
elseif releaseType == "BSR"
    release_I0(16) = release_ratio*wt_I0(13);
    release_I0(22) = release_ratio*wt_I0(19);
elseif releaseType == "FOR"
    release_I0(22) = release_ratio*wt_I0(19);
end

[~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts); 
% calculate drive success metric; in Khamis, measured total vector pop
% after 1 year
% outputMat(matInd,j) = sum(xout(end,19:end))/2000;

% drive frequency
ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
% transgene frequency
allele_freq = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));


plot(allele_freq,'linewidth',1.5)
ylim([0,1])

%% sandbox for 1LUD

subplot(2,1,1)

releaseVec      = 0.0:0.01:1.0;         
outputMat       = zeros(1,length(releaseVec));
tVec            = 0:365;        

%%%%%%%%%%%%%%%%%%%%%
%%%% HOMOLOGOUS %%%%%
%%%%%%%%%%%%%%%%%%%%%

matInd                  = 1;
params                  = struct();
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 16;
params.muX      = 0.03;
params.muZ      = log(10/9);
params.muY      = log(10/9);
params.m        = 0.1;
params.s_t      = 0.9;
params.s_a      = 0.05;
N               = 2*(10^3); 

params.lethality_type   = "EA";
params.supp             = "SS";

releaseType             = "BSR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

beta            = 0.9;
alpha           = (m/(2*2*1000*muZ))*((exp(lambda*m/(2*muZ)-m-muX) - 1)^(1/beta));
fInv            = @(x) (1/alpha)*((exp(x)-1).^(1/beta));
fInv_LAMBDA     = fInv(LAMBDA); 

params.beta     = beta;
params.alpha    = alpha;                

wt_I0                   = zeros(1,24);    
wt_I0(19)               = N;
wt_I0(7)                = (1/2)*fInv_LAMBDA;
wt_I0(1)                = (1/2)*fInv_LAMBDA; 
wt_I0(13)               = (m/(2*muZ))*fInv_LAMBDA;

for j = 1:length(releaseVec)
    disp(j);
    release_I0          = wt_I0;
    release_ratio       = releaseVec(j)/(1-releaseVec(j));

    if releaseType == "MOR"
        release_I0(16) = release_ratio*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = release_ratio*wt_I0(13);
        release_I0(22) = release_ratio*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = release_ratio*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts); 
    % calculate drive success metric
    outputMat(matInd,j) = sum(xout(end,19:end))/2000;

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24

end

plot(releaseVec, outputMat,'linewidth',1.5)
ylim([0,1])
yticks([0,0.2,0.4,0.6,0.8,1.0])
xticks([0,0.2,0.4,0.6,0.8,1.0])
grid on
title('(a)')
set(gca, 'fontsize',16)

%% sandbox for homing drive individual replicates

tVec            = 0:(10*365);        

release_ratio   = 0.01;    

% set parameters
params          = struct();
params.lambda   = 8;
params.CONV_EFF = 1; 
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.s        = 1.0;
params.h        = 0.0;
N               = 2*(10^4); 

params.lethality_type   = "LA";
params.lethality_case   = "BSL";

opts                    = odeset('RelTol',1e-9,'AbsTol',1e-10,'NonNegative',[1:12]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
beta            = 0.5;
alpha           = 0.0321302;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

% [AA_male_juveniles, Aa_male_juveniles,      2
%  aa_male_juveniles, AA_female_juveniles,    4
%  Aa_female_juveniles, aa_female_juveniles,  6
%  AA_adult_males, Aa_adult_males,            8
%  aa_adult_males, AA_adult_females,          10
%  Aa_adult_females, aa_adult_females]        12

wt_I0                   = zeros(1,12);    
wt_I0(10)               = N;
wt_I0(4)                = (1/2)*fInv_LAMBDA;
wt_I0(1)                = (1/2)*fInv_LAMBDA; 
wt_I0(7)                = (m/(2*muZ))*fInv_LAMBDA;
% release transgenic males
release_I0              = wt_I0;
release_I0(8)           = ceil(release_ratio*wt_I0(7));

params.beta             = 0.9;
params.alpha            = 0.014543805956894;

% see how the homing drive performs
[~,xout] = ode45(@(t,x) DDGD_iii_HOM(t,x,params),tVec,release_I0,opts); 
plot(xout)


%% STRESS TESTING to ensure that all alpha values are correct for each DD 
% functional form (included point checks in each function to ensure
% offspring genotype probabilities sum to 1, barring numerical error)

% population parameters
params = struct();
params.lambda   = 8;        % egg-laying rate
params.muX      = 0.029;    % density-independent mortality rate
params.muZ      = 0.28;     % male adult mortality
params.muY      = 0.1;      % female adult mortality
params.m        = 0.14;     % maturation rate
params.M        = 10000;    % total human pop.
params.N        = 20000;    % total female mosquito pop. 
% epi parameters
params.a        = 1.9;      % human-biting rate/female mosq.
params.b        = 0.75;     % prop. of bits resulting in infection (M->H)
params.c        = 0.75;     % prop. of bits resulting in infection (H->M)
params.gamma    = 1/3;      % recovery rate
% gene drive parameters
params.s_a      = 0.1;      % ambient transgene fitness cost
params.s_t      = 1;        % toxin transgene fitness cost
% ode inputs
opts            = odeset('RelTol',1e-8,'AbsTol',1e-9);
tmax            = 5*365;
dt              = 0.01;
tVec            = 0:dt:tmax;

% equilibrium pops
XStar = (params.muY/params.m)*params.N; 
XHatStar = XStar;
YHatStar = (params.muY/params.muZ)*params.N; % 7142.86;

% initial conditions for system (pre-release equilibrium)
IC = [XHatStar, XStar, ...                  % 2
        0, 0, ...                           % 4  
        0, 0, ...                           % 6  
        0, 0, ...                           % 8
        0, 0, ...                           % 10
        0, 0, ...                           % 12
        0, 0, ...                           % 14
        0, 0, ...                           % 16
        0, 0, ...                           % 18
        YHatStar, 0, ...                    % 20
        0, 0,  ...                          % 22  
        0, 0, ...                           % 24  
        0, 0, ...                           % 26
        0, params.N, ...                    % 28
        0, 0, ...                           % 30
        0, 0, ...                           % 32
        0, 0, ...                           % 34
        0, 0, ...                           % 36
        18535.13565, 8879.397446];          % 38 

% run a sim for the first case (gen log, beta = 0.5)
[t,xout1] = ode45(@(t,x) DDGD_ii(t,x,params,0.0321302,0.5,"BSL","SS","EA"),tVec,IC,opts); 
% run a sim for the second case (gen log, beta = 1.0)
[t,xout2] = ode45(@(t,x) DDGD_ii(t,x,params,1.90085*10^(-4),1,"BSL","SS","EA"),tVec,IC,opts); 
% run a sim for the third case (gen log, beta = 1.5)
[t,xout3] = ode45(@(t,x) DDGD_ii(t,x,params,1.12456*10^(-6),1.5,"BSL","SS","EA"),tVec,IC,opts); 
% run a sim for the fourth case (log)
[t,xout4] = ode45(@(t,x) DDGD_iii(t,x,params,0.014543805956894,0.9,"BSL","SS","EA"),tVec,IC,opts); 

% plot everything to ensure equilibria are correct
subplot(2,2,1)
plot(xout1)
subplot(2,2,2)
plot(xout2)
subplot(2,2,3)
plot(xout3)
subplot(2,2,4)
plot(xout4)


%% 1.1. restore time (days) vs. relative perturbation plot
% comparison of DD strength and time it takes for the system to return to
% equilibrium (<10^-4) following some perturbation to BOTH adult pops

% population parameters
params = struct();
params.lambda   = 8;        % egg-laying rate
params.muX      = 0.029;    % density-independent mortality rate
params.muZ      = 0.28;     % male adult mortality
params.muY      = 0.1;      % female adult mortality
params.m        = 0.14;     % maturation rate
params.M        = 10000;    % total human pop.
params.N        = 2*params.M;      % total female mosquito pop. 
% epi parameters
params.a        = 1.9;      % human-biting rate/female mosq.
params.b        = 0.75;     % prop. of bits resulting in infection (M->H)
params.c        = 0.75;     % prop. of bits resulting in infection (H->M)
params.gamma    = 1/3;      % recovery rate
% gene drive parameters
params.s_a      = 0.1;      % ambient transgene fitness cost
params.s_t      = 1;        % toxin transgene fitness cost
% params.init_ratio = 0;      % transgene:wild-type for initial release

% ode inputs
opts            = odeset('RelTol',1e-8,'AbsTol',1e-9);
tmax            = 5*365;
dt              = 0.01;
tVec            = 0:dt:tmax;
% threshold to determine when population has reached equilibrium 
conThresh = 10^(-4); 

% equilibrium pops
XStar = (params.muY/params.m)*params.N; 
XHatStar = XStar;
YHatStar = (params.muY/params.muZ)*params.N; % 7142.86;

% input/output key
% [AABB_larvae_male, AABB_larvae_female     2
%  AABb_larvae_male, AABb_larvae_female     4
%  AAbb_larvae_male, AAbb_larvae_female     6
%  AaBB_larvae_male, AaBB_larvae_female     8
%  AaBb_larvae_male, AaBb_larvae_female     10
%  Aabb_larvae_male, Aabb_larvae_female     12
%  aaBB_larvae_male, aaBB_larvae_female     14
%  aaBb_larvae_male, aaBb_larvae_female     16
%  aabb_larvae_male, aabb_larvae_female     18
%  AABB_adult_males, AABb_adult_males       20
%  AAbb_adult_males, AaBB_adult_males       22
%  AaBb_adult_males, Aabb_adult_males       24
%  aaBB_adult_males, aaBb_adult_males       26
%  aabb_adult_males, AABB_adult_females     28
%  AABb_adult_females, AAbb_adult_females   30
%  AaBB_adult_females, AaBb_adult_females   32
%  Aabb_adult_females, aaBB_adult_females   34
%  aaBb_adult_females, aabb_adult_females   36
%  infected_mosqs, infected_humans]         38

% test time to reach equilibrium following a 1% perturbation to the
% male and female mosquito populations

pVec = 0:0.0025:0.1;
N = params.N;
restoreVec_i = zeros(1,length(pVec)); % restorative force 
restoreVec_ii = zeros(1,length(pVec)); % restorative force 
restoreVec_iii = zeros(1,length(pVec)); % restorative force
restoreVec_iv = zeros(1,length(pVec)); % restorative force
restoreVec_i(1) = 0;
restoreVec_ii(1) = 0;
restoreVec_iii(1) = 0; 
restoreVec_iv(1) = 0; 
for i = 2:length(pVec)
    % run a simulation and find the time corresponding to when the male AND
    % female adult mosquito populations returned to equilibrium
    IC = [XHatStar, XStar, ...                  % 2
            0, 0, ...                           % 4  
            0, 0, ...                           % 6  
            0, 0, ...                           % 8
            0, 0, ...                           % 10
            0, 0, ...                           % 12
            0, 0, ...                           % 14
            0, 0, ...                           % 16
            0, 0, ...                           % 18
            YHatStar*(1+pVec(i)), 0, ...        % 20
            0, 0,  ...                          % 22  
            0, 0, ...                           % 24  
            0, 0, ...                           % 26
            0, N*(1+pVec(i)), ...               % 28
            0, 0, ...                           % 30
            0, 0, ...                           % 32
            0, 0, ...                           % 34
            0, 0, ...                           % 36
            18535.13565, 8879.397446];          % 38 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% first case (gen log, beta=0.5) %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    [t,xout] = ode45(@(t,x) DDGD_ii(t,x,params,0.0321302,0.5,"BSL","SS","EA"),tVec,IC,opts); 

    % when do the state variable amounts return to equilibrium? 
    XHat_per = flip(abs(xout(:,1) - xout(1,1))/xout(1,1)); 
    XHat_index = t(length(t)- find(XHat_per > conThresh,1) + 1);

    X_per = flip(abs(xout(:,2) - xout(1,2))/xout(1,2));
    X_index = t(length(t) - find(X_per > conThresh,1) + 1);

    YHat_per = abs(xout(:,19) - YHatStar)/YHatStar;
    YHat_index = t(find(YHat_per < conThresh, 1)); 

    Y_per = abs(xout(:,28) - N)/N; 
    Y_index = t(find(Y_per < conThresh, 1)); 

    if (isempty(Y_index) || isempty(YHat_index) || isempty(X_index) || isempty(XHat_index)) 
        restoreVec_i(i) = NaN; 
    else
        restoreVec_i(i) = max([Y_index, YHat_index, XHat_index, X_index]); 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% second case (gen log, beta=1) %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [t,xout] = ode45(@(t,x) DDGD_ii(t,x,params,1.90085*10^(-4),1,"BSL","SS","EA"),tVec,IC,opts); 

    % when do the state variable amounts return to equilibrium? 
    XHat_per = flip(abs(xout(:,1) - xout(1,1))/xout(1,1)); 
    XHat_index = t(length(t)- find(XHat_per > conThresh,1) + 1);

    X_per = flip(abs(xout(:,2) - xout(1,2))/xout(1,2));
    X_index = t(length(t) - find(X_per > conThresh,1) + 1);

    YHat_per = abs(xout(:,19) - YHatStar)/YHatStar;
    YHat_index = t(find(YHat_per < conThresh, 1)); 

    Y_per = abs(xout(:,28) - N)/N; 
    Y_index = t(find(Y_per < conThresh, 1)); 

    if (isempty(Y_index) || isempty(YHat_index) || isempty(X_index) || isempty(XHat_index)) 
        restoreVec_ii(i) = NaN; 
    else
        restoreVec_ii(i) = max([Y_index, YHat_index, XHat_index, X_index]); 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% third case (gen log, beta=1.5) %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [t,xout] = ode45(@(t,x) DDGD_ii(t,x,params,1.12456*10^(-6),1.5,"BSL","SS","EA"),tVec,IC,opts); 

    % when do the state variable amounts return to equilibrium? 
    XHat_per = flip(abs(xout(:,1) - xout(1,1))/xout(1,1)); 
    XHat_index = t(length(t)- find(XHat_per > conThresh,1) + 1);

    X_per = flip(abs(xout(:,2) - xout(1,2))/xout(1,2));
    X_index = t(length(t) - find(X_per > conThresh,1) + 1);

    YHat_per = abs(xout(:,19) - YHatStar)/YHatStar;
    YHat_index = t(find(YHat_per < conThresh, 1)); 

    Y_per = abs(xout(:,28) - N)/N; 
    Y_index = t(find(Y_per < conThresh, 1)); 

    if (isempty(Y_index) || isempty(YHat_index) || isempty(X_index) || isempty(XHat_index)) 
        restoreVec_iii(i) = NaN; 
    else
        restoreVec_iii(i) = max([Y_index, YHat_index, XHat_index, X_index]); 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% fourth case (logarithmic) %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [t,xout] = ode45(@(t,x) DDGD_iii(t,x,params,0.014543805956894,0.9,"BSL","SS","EA"),tVec,IC,opts); 

    % when do the state variable amounts return to equilibrium? 
    XHat_per = flip(abs(xout(:,1) - xout(1,1))/xout(1,1)); 
    XHat_index = t(length(t)- find(XHat_per > conThresh,1) + 1);

    X_per = flip(abs(xout(:,2) - xout(1,2))/xout(1,2));
    X_index = t(length(t) - find(X_per > conThresh,1) + 1);

    YHat_per = abs(xout(:,19) - YHatStar)/YHatStar;
    YHat_index = t(find(YHat_per < conThresh, 1)); 

    Y_per = abs(xout(:,28) - N)/N; 
    Y_index = t(find(Y_per < conThresh, 1)); 

    if (isempty(Y_index) || isempty(YHat_index) || isempty(X_index) || isempty(XHat_index)) 
        restoreVec_iv(i) = NaN; 
    else
        restoreVec_iv(i) = max([Y_index, YHat_index, XHat_index, X_index]); 
    end

    disp([i, restoreVec_i(i),restoreVec_ii(i),restoreVec_iii(i), restoreVec_iv(i)]); 

end

% subplot(1,2,1)
plot(pVec, restoreVec_i,'-k','Linewidth',1.5)
% title('(a)')
xlim([0,0.1]);
ylim([0,10^3]);
hold on
plot(pVec, restoreVec_ii,'--k','Linewidth',1.5)
plot(pVec, restoreVec_iii,'-.k','Linewidth',1.5)
plot(pVec, restoreVec_iv,':k','Linewidth',1.5)
% grid on
set(gca,'Fontsize',16,'YScale','log');
xlabel('relative size of perturbation','Interpreter','latex')
ylabel('return time (days)','Interpreter','latex')
% legend("generalized logistic $\beta = 0.5$","generalized logistic $\beta = 1$",...
%     "generalized logistic $\beta = 1.5$",...
%     "logarithmic","location","best",...
%     "interpreter","latex")
hold off

%% 1.2. recreating the results of Bellows and plotting our own findings; 
% the results of Bellows presented in his Fig. 6 claims that
% overcompensation results from his eqtn (2) if you integrate over a short
% time scale (t=1 in his case). However, the following shows this is not
% the case

% Figure plotted here is meant to be part (b) of the fig generated in 1.1.

N0Vec = 1:150;
ePop = zeros(6,length(N0Vec));
tVec = 0:0.001:5;
ode_options = odeset('RelTol',1e-6);
% max time to check if approximations hold for small values of t (Bellows
% used t=1)
maxTime = 1;

for j = 1:length(N0Vec)
    disp(j);
    IC = N0Vec(j);

    [~, xout5a] = ode45(@bellows_5,[0:0.01:maxTime],IC,ode_options, 0.01, 1);
    [~, xout5b] = ode45(@bellows_5,[0:0.01:maxTime],IC,ode_options, 0.01, 5);
    [~, xout5c] = ode45(@bellows_5,[0:0.01:maxTime],IC,ode_options, 0.02, 5);

    % integrated solutions
    ePop(1,j) = xout5a(end); 
    ePop(2,j) = xout5b(end); 
    ePop(3,j) = xout5c(end); 
    % Bellows' curves assuming fixed pop
    ePop(4,j) = IC*exp(-log(1+(0.01*IC)^1)); 
    ePop(5,j) = IC*exp(-log(1+(0.01*IC)^5)); 
    ePop(6,j) = IC*exp(-log(1+(0.02*IC)^5)); 
end


% subplot(1,2,2)
plot(N0Vec,ePop(1,:),'r--','Linewidth',1.5);
% title('(b)');
xlabel('initial number','Interpreter','latex'); 
ylabel('number surviving','Interpreter','latex'); 
set(gca,'Fontsize',16);
hold on
plot(N0Vec,ePop(2,:),'r-.','Linewidth',1.5);
plot(N0Vec,ePop(3,:),'r-','Linewidth',1.5);
plot(N0Vec,ePop(4,:),'k--','Linewidth',1.5);
plot(N0Vec,ePop(5,:),'k-.','Linewidth',1.5);
plot(N0Vec,ePop(6,:),'k-','Linewidth',1.5);

% do the approximations hold over small time scales?


%% 1.3. Confirming that analytic equilibria agree with numerical results 
% using the system with a general genetic load and with parameters
% befitting Ae. aegypti; conditions are tested both with (g > 0) and
% without (g = 0) a genetic load
params          = struct();
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
tVec            = 0:365; 

% for calculating equilibria
lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;

% [female_juveniles, male_juveniles, ...
% [female_adults, male_adults]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% generalized logistic case %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% beta = 0.5, EA

% parameter values
params.beta             = 0.5;
params.lethality_type   = "EA";
params.g                = 0.35;

beta                    = params.beta;
g                       = params.g;

% determine alpha
switch beta
case 0.5
    alpha = 0.0321302;
case 1
    alpha = 1.90085*10^(-4);
case 1.5
    alpha = 1.12456*10^(-6);
otherwise
    error('Error! Beta must be 0.5, 1, or 1.5.')
end

LAMBDA                  = (m*lambda/(2*muY))-muX-m;
LAMBDA_g                = (1-g)*(m*lambda/(2*muY))-muX-m;
fInv                    = @(x) (x/alpha).^(1/beta);

% all other parameters and individual sims...
params.lethality_case   = "BSL";
% initial condition
I0 = [(1/2)*fInv(LAMBDA_g), (1/2)*fInv(LAMBDA_g), (m/muY)*(1/2)*fInv(LAMBDA_g), ...
    (m/muZ)*(1/2)*fInv(LAMBDA_g)];
% I0 = [((1-g)/(2-g))*fInv(LAMBDA_g), (1/(2-g))*fInv(LAMBDA_g), (m/muY)*((1-g)/(2-g))*fInv(LAMBDA_g), ...
%     (m/muZ)*(1/(2-g))*fInv(LAMBDA_g)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,I0,opts); 

subplot(4,5,1)
plot(tVec, xout)
title('GL, EA, BSL, $\beta=0.5$','interpreter','latex')
ylim([0,max(max(xout))])

params.lethality_case   = "FSL";
I0 = [((1-g)/(2-g))*fInv(LAMBDA_g), (1/(2-g))*fInv(LAMBDA_g), (m/muY)*((1-g)/(2-g))*fInv(LAMBDA_g), ...
    (m/muZ)*(1/(2-g))*fInv(LAMBDA_g)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,I0,opts); 

subplot(4,5,2)
plot(tVec, xout)
title('GL, EA, FSL, $\beta=0.5$','interpreter','latex')
ylim([0,max(max(xout))])

params.lethality_case   = "MSL";
I0 = [(1/(2-g))*fInv(LAMBDA), ((1-g)/(2-g))*fInv(LAMBDA), (m/muY)*(1/(2-g))*fInv(LAMBDA), ...
    (m/muZ)*((1-g)/(2-g))*fInv(LAMBDA)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,I0,opts); 

subplot(4,5,3)
plot(tVec, xout)
title('GL, EA, MSL, $\beta=0.5$','interpreter','latex')
ylim([0,max(max(xout))])

%%% beta = 0.5, LA
params.lethality_type   = "LA";

params.lethality_case   = "BSL";
I0 = [(1/2)*fInv(LAMBDA_g), (1/2)*fInv(LAMBDA_g), (m/muY)*((1-g)/2)*fInv(LAMBDA_g), ...
    (m/muZ)*((1-g)/2)*fInv(LAMBDA_g)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,I0,opts); 

subplot(4,5,4)
plot(tVec, xout)
title('GL, LA, BSL, $\beta=0.5$','interpreter','latex')
ylim([0,max(max(xout))])

params.lethality_case   = "FSL";
I0 = [(1/2)*fInv(LAMBDA_g), (1/2)*fInv(LAMBDA_g), (m/muY)*((1-g)/2)*fInv(LAMBDA_g), ...
    (m/muZ)*(1/2)*fInv(LAMBDA_g)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,I0,opts); 

subplot(4,5,5)
plot(tVec, xout)
title('GL, LA, FSL, $\beta=0.5$','interpreter','latex')
ylim([0,max(max(xout))])

params.lethality_case   = "MSL";
I0 = [(1/2)*fInv(LAMBDA), (1/2)*fInv(LAMBDA), (m/(2*muY))*fInv(LAMBDA), ...
    (m/(2*muZ))*(1-g)*fInv(LAMBDA)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,I0,opts); 

subplot(4,5,6)
plot(tVec, xout)
title('GL, LA, MSL, $\beta=0.5$','interpreter','latex')
ylim([0,max(max(xout))])


%%% beta = 1.0, EA [can be changed to beta = 1.5 very easily...]

% parameter values
params.beta             = 1.0;
params.lethality_type   = "EA";
beta                    = params.beta;

% determine alpha
switch beta
case 0.5
    alpha = 0.0321302;
case 1
    alpha = 1.90085*10^(-4);
case 1.5
    alpha = 1.12456*10^(-6);
otherwise
    error('Error! Beta must be 0.5, 1, or 1.5.')
end

LAMBDA                  = (m*lambda/(2*muY))-muX-m;
LAMBDA_g                = (1-g)*(m*lambda/(2*muY))-muX-m;
fInv                    = @(x) (x/alpha).^(1/beta);

% all other parameters and individual sims...
params.lethality_case   = "BSL";
% initial condition
I0 = [(1/2)*fInv(LAMBDA_g), (1/2)*fInv(LAMBDA_g), (m/muY)*(1/2)*fInv(LAMBDA_g), ...
    (m/muZ)*(1/2)*fInv(LAMBDA_g)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,I0,opts); 

subplot(4,5,7)
plot(tVec, xout)
title('GL, EA, BSL, $\beta=1.0$','interpreter','latex')
ylim([0,max(max(xout))])

params.lethality_case   = "FSL";
I0 = [((1-g)/(2-g))*fInv(LAMBDA_g), (1/(2-g))*fInv(LAMBDA_g), (m/muY)*((1-g)/(2-g))*fInv(LAMBDA_g), ...
    (m/muZ)*(1/(2-g))*fInv(LAMBDA_g)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,I0,opts); 

subplot(4,5,8)
plot(tVec, xout)
title('GL, EA, FSL, $\beta=1.0$','interpreter','latex')
ylim([0,max(max(xout))])

params.lethality_case   = "MSL";
I0 = [(1/(2-g))*fInv(LAMBDA), ((1-g)/(2-g))*fInv(LAMBDA), (m/muY)*(1/(2-g))*fInv(LAMBDA), ...
    (m/muZ)*((1-g)/(2-g))*fInv(LAMBDA)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,I0,opts); 

subplot(4,5,9)
plot(tVec, xout)
title('GL, EA, MSL, $\beta=1.0$','interpreter','latex')
ylim([0,max(max(xout))])

%%% beta = 1.0, LA
params.lethality_type   = "LA";

params.lethality_case   = "BSL";
I0 = [(1/2)*fInv(LAMBDA_g), (1/2)*fInv(LAMBDA_g), (m/muY)*((1-g)/2)*fInv(LAMBDA_g), ...
    (m/muZ)*((1-g)/2)*fInv(LAMBDA_g)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,I0,opts); 

subplot(4,5,10)
plot(tVec, xout)
title('GL, LA, BSL, $\beta=1.0$','interpreter','latex')
ylim([0,max(max(xout))])

params.lethality_case   = "FSL";
I0 = [(1/2)*fInv(LAMBDA_g), (1/2)*fInv(LAMBDA_g), (m/muY)*((1-g)/2)*fInv(LAMBDA_g), ...
    (m/muZ)*(1/2)*fInv(LAMBDA_g)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,I0,opts); 

subplot(4,5,11)
plot(tVec, xout)
title('GL, LA, FSL, $\beta=1.0$','interpreter','latex')
ylim([0,max(max(xout))])

params.lethality_case   = "MSL";
I0 = [(1/2)*fInv(LAMBDA), (1/2)*fInv(LAMBDA), (m/(2*muY))*fInv(LAMBDA), ...
    (m/(2*muZ))*(1-g)*fInv(LAMBDA)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,I0,opts); 

subplot(4,5,12)
plot(tVec, xout)
title('GL, LA, MSL, $\beta=1.0$','interpreter','latex')
ylim([0,max(max(xout))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% logarithmic case %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% log case, EA

% parameter values
params.beta             = 0.9;
params.alpha            = 0.014543805956894;
params.lethality_type   = "EA";
params.g                = 0.35;

alpha                   = params.alpha;
beta                    = params.beta;
g                       = params.g;

LAMBDA                  = (m*lambda/(2*muY))-muX-m;
LAMBDA_g                = (1-g)*(m*lambda/(2*muY))-muX-m;
fInv                    = @(x) (1/alpha)*(exp(x) -1)^(1/beta);

% all other parameters and individual sims...
params.lethality_case   = "BSL";
I0 = [(1/2)*fInv(LAMBDA_g), (1/2)*fInv(LAMBDA_g), (m/muY)*(1/2)*fInv(LAMBDA_g), ...
    (m/muZ)*(1/2)*fInv(LAMBDA_g)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_iii(t,x,params),tVec,I0,opts); 

subplot(4,5,13)
plot(tVec, xout)
title('Log, EA, BSL','interpreter','latex')
ylim([0,max(max(xout))])

params.lethality_case   = "FSL";
I0 = [((1-g)/(2-g))*fInv(LAMBDA_g), (1/(2-g))*fInv(LAMBDA_g), (m/muY)*((1-g)/(2-g))*fInv(LAMBDA_g), ...
    (m/muZ)*(1/(2-g))*fInv(LAMBDA_g)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_iii(t,x,params),tVec,I0,opts); 

subplot(4,5,14)
plot(tVec, xout)
title('Log, EA, FSL','interpreter','latex')
ylim([0,max(max(xout))])

params.lethality_case   = "MSL";
I0 = [(1/(2-g))*fInv(LAMBDA), ((1-g)/(2-g))*fInv(LAMBDA), (m/muY)*(1/(2-g))*fInv(LAMBDA), ...
    (m/muZ)*((1-g)/(2-g))*fInv(LAMBDA)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_iii(t,x,params),tVec,I0,opts); 

subplot(4,5,15)
plot(tVec, xout)
title('Log, EA, MSL','interpreter','latex')
ylim([0,max(max(xout))])

%%% log case, LA
params.lethality_type   = "LA";

params.lethality_case   = "BSL";
I0 = [(1/2)*fInv(LAMBDA_g), (1/2)*fInv(LAMBDA_g), (m/muY)*((1-g)/2)*fInv(LAMBDA_g), ...
    (m/muZ)*((1-g)/2)*fInv(LAMBDA_g)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_iii(t,x,params),tVec,I0,opts); 

subplot(4,5,16)
plot(tVec, xout)
title('Log, LA, BSL','interpreter','latex')
ylim([0,max(max(xout))])

params.lethality_case   = "FSL";
I0 = [(1/2)*fInv(LAMBDA_g), (1/2)*fInv(LAMBDA_g), (m/muY)*((1-g)/2)*fInv(LAMBDA_g), ...
    (m/muZ)*(1/2)*fInv(LAMBDA_g)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_iii(t,x,params),tVec,I0,opts); 

subplot(4,5,17)
plot(tVec, xout)
title('Log, LA, FSL','interpreter','latex')
ylim([0,max(max(xout))])

params.lethality_case   = "MSL";
I0 = [(1/2)*fInv(LAMBDA), (1/2)*fInv(LAMBDA), (m/(2*muY))*fInv(LAMBDA), ...
    (m/(2*muZ))*(1-g)*fInv(LAMBDA)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
[~,xout] = ode45(@(t,x) vanilla_DDGD_iii(t,x,params),tVec,I0,opts); 

subplot(4,5,18)
plot(tVec, xout)
title('Log, LA, MSL','interpreter','latex')
ylim([0,max(max(xout))])

%% 1.4. Plots for LA system for FSL and BSL; generating plots comparing 
% different cases for the system with a general genetic load; both 
% numerical and analytical results are shown,then a final plot produced 
% based on these results. Only FSL is plotted since this is equivalent 
% to BSL (in terms of adult female pop. at control equilibrium), and MSL 
% does not reduce female control equil. 

% simulation parameters
params          = struct();
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
tVec            = 0:365; 

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
N               = 2*10^4;

params.lethality_type   = "LA";
params.lethality_case   = "BSL";

opts = odeset('RelTol',1e-8,'AbsTol',1e-9);

gVec = 0.0:0.05:1.0; 
% initial conditions (pre-release equil.)
WT_I0 = [(muY/m)*N, (muY/m)*N, (muY/muZ)*N, N];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% BI-SEX LETHAL (LA) %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % gen logistic, beta = 0.5
% beta        = 0.5; 
% alpha       = 0.0321302;
% fInv        = @(x) (x/alpha).^(1/beta);
% 
% tmp         = (1-gVec)*(m*lambda/(2*muY))-muX-m;
% genLog1_an  = (m/muY)*(1/2)*(1-gVec).*fInv(tmp);
% 
% % gen logistic, beta = 1.0
% beta        = 1.0; 
% alpha       = 1.90085*10^(-4);
% fInv        = @(x) (x/alpha).^(1/beta);
% 
% genLog2_an  = (m/muY)*(1/2)*(1-gVec).*fInv(tmp);
% 
% % gen logistic, beta = 1.5
% beta        = 1.5; 
% alpha       = 1.12456*10^(-6);
% fInv        = @(x) (x/alpha).^(1/beta);
% 
% genLog3_an  = (m/muY)*(1/2)*(1-gVec).*fInv(tmp);
% 
% % logarithmic
% beta        = 0.9; 
% alpha       = 0.014543805956894;
% fInv        = @(x) (1/alpha)*(exp(x) -1).^(1/beta);
% 
% log1_an     = (m/muY)*(1/2)*(1-gVec).*fInv(tmp);
% 
% % vectors
% genLog1_num = zeros(1,length(gVec));
% genLog2_num = zeros(1,length(gVec));
% genLog3_num = zeros(1,length(gVec));
% log1_num    = zeros(1,length(gVec));
% 
% % make sure LAMBDA_g or LAMBDA is greater than zero...
% LAMBDA_g    = (1-gVec)*(m*lambda/(2*muY))-muX-m;
% 
% genLog1_an(find(LAMBDA_g < 0))     = nan;
% genLog2_an(find(LAMBDA_g < 0))     = nan;
% genLog3_an(find(LAMBDA_g < 0))     = nan;
% log1_an(find(LAMBDA_g < 0))        = nan;
% 
% for i = 1:length(gVec)
%     disp(i);
%     params.g        = gVec(i);
% 
%     params.beta     = 0.5;
%     [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
%     % store adult female population
%     genLog1_num(i)  = xout(end,3);
% 
%     params.beta     = 1.0;
%     [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
%     genLog2_num(i)  = xout(end,3);    
% 
%     params.beta     = 1.5;
%     [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
%     genLog3_num(i)  = xout(end,3);
% 
%     params.beta     = 0.9;
%     [~,xout] = ode45(@(t,x) vanilla_DDGD_iii(t,x,params),tVec,WT_I0,opts); 
%     log1_num(i)     = xout(end,3);    
% 
%     % [female_juveniles, male_juveniles, ...
%     % [female_adults, male_adults]
% end
% 
% subplot(1,3,1)
% plot(gVec, genLog1_an, '-k','linewidth',1.5)
% hold on
% plot(gVec, genLog1_num, '--r','linewidth',1.5)
% plot(gVec, genLog2_an, '-k','linewidth',1.5)
% plot(gVec, genLog2_num, '--r','linewidth',1.5)
% plot(gVec, genLog3_an, '-k','linewidth',1.5)
% plot(gVec, genLog3_num, '--r','linewidth',1.5)
% plot(gVec, log1_an, '-k','linewidth',1.5)
% plot(gVec, log1_num, '--r','linewidth',1.5)
% title('BSL','interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% FEMALE-SEX LETHAL (LA) %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.lethality_case   = "FSL";

% gen logistic, beta = 0.5
beta        = 0.5; 
alpha       = 0.0321302;
fInv        = @(x) (x/alpha).^(1/beta);

tmp         = (1-gVec)*(m*lambda/(2*muY))-muX-m;
genLog1_an  = (m/muY)*(1/2)*(1-gVec).*fInv(tmp);

% gen logistic, beta = 1.0
beta        = 1.0; 
alpha       = 1.90085*10^(-4);
fInv        = @(x) (x/alpha).^(1/beta);

genLog2_an  = (m/muY)*(1/2)*(1-gVec).*fInv(tmp);

% gen logistic, beta = 1.5
beta        = 1.5; 
alpha       = 1.12456*10^(-6);
fInv        = @(x) (x/alpha).^(1/beta);

genLog3_an  = (m/muY)*(1/2)*(1-gVec).*fInv(tmp);

% logarithmic
beta        = 0.9; 
alpha       = 0.014543805956894;
fInv        = @(x) (1/alpha)*(exp(x) -1).^(1/beta);

log1_an     = (m/muY)*(1/2)*(1-gVec).*fInv(tmp);

% vectors
genLog1_num = zeros(1,length(gVec));
genLog2_num = zeros(1,length(gVec));
genLog3_num = zeros(1,length(gVec));
log1_num    = zeros(1,length(gVec));

% make sure LAMBDA_g or LAMBDA is greater than zero...
LAMBDA_g    = (1-gVec)*(m*lambda/(2*muY))-muX-m;

genLog1_an(find(LAMBDA_g < 0))     = nan;
genLog2_an(find(LAMBDA_g < 0))     = nan;
genLog3_an(find(LAMBDA_g < 0))     = nan;
log1_an(find(LAMBDA_g < 0))        = nan;

for i = 1:length(gVec)
    disp(i);
    params.g        = gVec(i);

    params.beta     = 0.5;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
    % store adult female population
    genLog1_num(i)  = xout(end,3);

    params.beta     = 1.0;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
    genLog2_num(i)  = xout(end,3);    

    params.beta     = 1.5;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
    genLog3_num(i)  = xout(end,3);

    params.beta     = 0.9;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_iii(t,x,params),tVec,WT_I0,opts); 
    log1_num(i)     = xout(end,3);    

    % [female_juveniles, male_juveniles, ...
    % [female_adults, male_adults]
end

% for a 1x2 plot... [this has been removed since we no longer show the
% inverse functions]
% subplot(1,2,1)

% plot(gVec, genLog1_an, '-k','linewidth',1.5)
% hold on
plot(gVec, genLog1_num/N, '-k','linewidth',1.5)
% title('(a)')
hold on
% plot(gVec, genLog2_an, '-k','linewidth',1.5)
plot(gVec, genLog2_num/N, '--k','linewidth',1.5)
% plot(gVec, genLog3_an, '-k','linewidth',1.5)
plot(gVec, genLog3_num/N, '-.k','linewidth',1.5)
% plot(gVec, log1_an, '-k','linewidth',1.5)
plot(gVec, log1_num/N, ':k','linewidth',1.5)
% title('FSL','interpreter','latex')
xlabel('control load, $g$','interpreter','latex')
ylabel('relative adult female population','interpreter','latex')
legend('gen. logistic, $\beta=0.5$', 'gen. logistic, $\beta=1.0$', ...
    'gen. logistic, $\beta=1.5$', 'logarithmic', ...
    'interpreter','latex');
set(gca,'fontsize',15)
ylim([0,1])

%% 1.5. Plotting the inverted functional forms of each per-capita DD 
% juvenile mortality considered to see how they compare... 

subplot(1,2,2)
xVec = 0:0.25:6;

% gen logistic, beta = 0.5
beta        = 0.5; 
alpha       = 0.0321302;
fInv        = @(x) (x/alpha).^(1/beta);

plot(xVec, fInv(xVec),'-k','linewidth',1.5);
title('(b)')
hold on

% gen logistic, beta = 1.0
beta        = 1.0; 
alpha       = 1.90085*10^(-4);
fInv        = @(x) (x/alpha).^(1/beta);
plot(xVec, fInv(xVec),'--k','linewidth',1.5);

% gen logistic, beta = 1.5
beta        = 1.5; 
alpha       = 1.12456*10^(-6);
fInv        = @(x) (x/alpha).^(1/beta);
plot(xVec, fInv(xVec),'-.k','linewidth',1.5);

% logarithmic
beta        = 0.9; 
alpha       = 0.014543805956894;
fInv        = @(x) (1/alpha)*(exp(x) -1).^(1/beta);
plot(xVec, fInv(xVec),':k','linewidth',1.5);

% add point for equilibrium
scatter(5.4310, 2*(muY/m)*N,'red','filled')
xlabel('per-capita density-dependent juvenile mortality rate','interpreter','latex');
ylabel('juvenile pop.','interpreter','latex');
set(gca, 'fontsize',15)

%% 1.6. Plots for EA system for both FSL and BSL; both numerical and 
% analytical results are calculated, then a final plot produced based on 
% these results. This is for the system with a general genetic load. 

% simulation parameters
params          = struct();
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
tVec            = 0:365; 

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
N               = 2*10^4;

params.lethality_type   = "EA";
params.lethality_case   = "BSL";

opts = odeset('RelTol',1e-8,'AbsTol',1e-9);

gVec = 0.0:0.05:1.0; 
% initial conditions (pre-release equil.)
WT_I0 = [(muY/m)*N, (muY/m)*N, (muY/muZ)*N, N];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% BI-SEX LETHAL (EA) %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gen logistic, beta = 0.5
beta        = 0.5; 
alpha       = 0.0321302;
fInv        = @(x) (x/alpha).^(1/beta);

tmp         = (1-gVec)*(m*lambda/(2*muY))-muX-m;
genLog1_an  = (m/muY)*(1/2).*fInv(tmp);

% gen logistic, beta = 1.0
beta        = 1.0; 
alpha       = 1.90085*10^(-4);
fInv        = @(x) (x/alpha).^(1/beta);

genLog2_an  = (m/muY)*(1/2).*fInv(tmp);

% gen logistic, beta = 1.5
beta        = 1.5; 
alpha       = 1.12456*10^(-6);
fInv        = @(x) (x/alpha).^(1/beta);

genLog3_an  = (m/muY)*(1/2).*fInv(tmp);

% logarithmic
beta        = 0.9; 
alpha       = 0.014543805956894;
fInv        = @(x) (1/alpha)*(exp(x) -1).^(1/beta);

log1_an     = (m/muY)*(1/2).*fInv(tmp);

% vectors
genLog1_num = zeros(1,length(gVec));
genLog2_num = zeros(1,length(gVec));
genLog3_num = zeros(1,length(gVec));
log1_num    = zeros(1,length(gVec));

% make sure LAMBDA_g or LAMBDA is greater than zero...
LAMBDA_g    = (1-gVec)*(m*lambda/(2*muY))-muX-m;

genLog1_an(find(LAMBDA_g < 0))     = nan;
genLog2_an(find(LAMBDA_g < 0))     = nan;
genLog3_an(find(LAMBDA_g < 0))     = nan;
log1_an(find(LAMBDA_g < 0))        = nan;

for i = 1:length(gVec)
    disp(i);
    params.g        = gVec(i);

    params.beta     = 0.5;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
    % store adult female population
    genLog1_num(i)  = xout(end,3);

    params.beta     = 1.0;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
    genLog2_num(i)  = xout(end,3);    

    params.beta     = 1.5;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
    genLog3_num(i)  = xout(end,3);

    params.beta     = 0.9;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_iii(t,x,params),tVec,WT_I0,opts); 
    log1_num(i)     = xout(end,3);    

    % [female_juveniles, male_juveniles, ...
    % [female_adults, male_adults]
end

subplot(1,2,2)
plot(gVec, genLog1_an/N, '-k','linewidth',1.5)
title('(b) BSL')
hold on
plot(gVec, genLog1_num/N, '--r','linewidth',1.5)
plot(gVec, genLog2_an/N, '-k','linewidth',1.5)
plot(gVec, genLog2_num/N, '--r','linewidth',1.5)
plot(gVec, genLog3_an/N, '-k','linewidth',1.5)
plot(gVec, genLog3_num/N, '--r','linewidth',1.5)
plot(gVec, log1_an/N, '-k','linewidth',1.5)
plot(gVec, log1_num/N, '--r','linewidth',1.5)
ylim([0,1])

% save results for plot comparing FSL and BSL EA performance
genLog1_num_BSL = genLog1_num/N;
genLog2_num_BSL = genLog2_num/N;
genLog3_num_BSL = genLog3_num/N;
log1_num_BSL    = log1_num/N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% FEMALE-SEX LETHAL (EA) %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.lethality_case   = "FSL";

% gen logistic, beta = 0.5
beta        = 0.5; 
alpha       = 0.0321302;
fInv        = @(x) (x/alpha).^(1/beta);

tmp         = (1-gVec)*(m*lambda/(2*muY))-muX-m;
genLog1_an  = (m/muY)*((1-gVec)./(2-gVec)).*fInv(tmp);

% gen logistic, beta = 1.0
beta        = 1.0; 
alpha       = 1.90085*10^(-4);
fInv        = @(x) (x/alpha).^(1/beta);

genLog2_an  = (m/muY)*((1-gVec)./(2-gVec)).*fInv(tmp);

% gen logistic, beta = 1.5
beta        = 1.5; 
alpha       = 1.12456*10^(-6);
fInv        = @(x) (x/alpha).^(1/beta);

genLog3_an  = (m/muY)*((1-gVec)./(2-gVec)).*fInv(tmp);

% logarithmic
beta        = 0.9; 
alpha       = 0.014543805956894;
fInv        = @(x) (1/alpha)*(exp(x) -1).^(1/beta);

log1_an     = (m/muY)*((1-gVec)./(2-gVec)).*fInv(tmp);

% vectors
genLog1_num = zeros(1,length(gVec));
genLog2_num = zeros(1,length(gVec));
genLog3_num = zeros(1,length(gVec));
log1_num    = zeros(1,length(gVec));

% make sure LAMBDA_g or LAMBDA is greater than zero...
LAMBDA_g    = (1-gVec)*(m*lambda/(2*muY))-muX-m;

genLog1_an(find(LAMBDA_g < 0))     = nan;
genLog2_an(find(LAMBDA_g < 0))     = nan;
genLog3_an(find(LAMBDA_g < 0))     = nan;
log1_an(find(LAMBDA_g < 0))        = nan;

for i = 1:length(gVec)
    disp(i);
    params.g        = gVec(i);

    params.beta     = 0.5;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
    % store adult female population
    genLog1_num(i)  = xout(end,3);

    params.beta     = 1.0;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
    genLog2_num(i)  = xout(end,3);    

    params.beta     = 1.5;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_ii(t,x,params),tVec,WT_I0,opts); 
    genLog3_num(i)  = xout(end,3);

    params.beta     = 0.9;
    [~,xout] = ode45(@(t,x) vanilla_DDGD_iii(t,x,params),tVec,WT_I0,opts); 
    log1_num(i)     = xout(end,3);    

    % [female_juveniles, male_juveniles, ...
    % [female_adults, male_adults]
end

% for a 1x2 plot... [this has been removed since we no longer show the
% inverse functions]
subplot(1,2,1)

plot(gVec, genLog1_an/N, '-k','linewidth',1.5)
hold on
plot(gVec, genLog1_num/N, '--r','linewidth',1.5)
title('(a) FSL')
plot(gVec, genLog2_an/N, '-k','linewidth',1.5)
plot(gVec, genLog2_num/N, '--r','linewidth',1.5)
plot(gVec, genLog3_an/N, '-k','linewidth',1.5)
plot(gVec, genLog3_num/N, '--r','linewidth',1.5)
plot(gVec, log1_an/N, '-k','linewidth',1.5)
plot(gVec, log1_num/N, '--r','linewidth',1.5)
% title('FSL','interpreter','latex')
xlabel('control load, $g$','interpreter','latex')
ylabel('relative adult female population','interpreter','latex')
set(gca,'fontsize',15)
ylim([0,1])

% save results for plot comparing FSL and BSL EA performance
genLog1_num_FSL = genLog1_num/N;
genLog2_num_FSL = genLog2_num/N;
genLog3_num_FSL = genLog3_num/N;
log1_num_FSL    = log1_num/N;

% plot for comparison in early-acting case
figure

subplot(1,2,1)
plot(gVec, genLog1_num_FSL, '-k','linewidth',1.5)
hold on
title('(a)')
plot(gVec, genLog2_num_FSL, '--k','linewidth',1.5)
plot(gVec, genLog3_num_FSL, '-.k','linewidth',1.5)
plot(gVec, log1_num_FSL, ':k','linewidth',1.5)
xticks(0:0.2:1.0)
xlabel('control load, $g$','interpreter','latex')
ylabel('relative female population','interpreter','latex')
set(gca,'fontsize',15)
ylim([0,1])

subplot(1,2,2)
plot(gVec, genLog1_num_BSL, '-k','linewidth',1.5)
hold on
title('(b)')
plot(gVec, genLog2_num_BSL, '--k','linewidth',1.5)
plot(gVec, genLog3_num_BSL, '-.k','linewidth',1.5)
plot(gVec, log1_num_BSL, ':k','linewidth',1.5)
xticks(0:0.2:1.0)
xlabel('control load, $g$','interpreter','latex')
ylabel('relative female population','interpreter','latex')
legend('gen. logistic, $\beta=0.5$', 'gen. logistic, $\beta=1.0$', ...
    'gen. logistic, $\beta=1.5$', 'logarithmic', ...
    'interpreter','latex');
set(gca,'fontsize',15)
ylim([0,1])

% and now plot FSL and BSL over each other
figure
plot(gVec, genLog1_num_FSL, '-k','linewidth',1.5)
hold on
plot(gVec, genLog2_num_FSL, '--k','linewidth',1.5)
plot(gVec, genLog3_num_FSL, '-.k','linewidth',1.5)
plot(gVec, log1_num_FSL, ':k','linewidth',1.5)
plot(gVec, genLog1_num_BSL, '-r','linewidth',1.5)
plot(gVec, genLog2_num_BSL, '--r','linewidth',1.5)
plot(gVec, genLog3_num_BSL, '-.r','linewidth',1.5)
plot(gVec, log1_num_BSL, ':r','linewidth',1.5)
xticks(0:0.2:1.0)
xlabel('control load, $g$','interpreter','latex')
ylabel('relative female population','interpreter','latex')
set(gca,'fontsize',15)
ylim([0,1])

%% 2.1. Determination of release ratios for successful invasion for each 
% case. Invasion is successful if it occurs within 10 years of release.
% This is for the single-locus UD case. We assume that the toxin is always
% fully lethal (s_t = 1). Also, recall that a BSR is required for a
% one-locus UD.

% NOTE: This section can simulate both FSL and BSL; plot at end can handle
% both, assuming BSL results are plotted after in a 2x2 fig. 

% close all 

params          = struct();

% specify lethality 
params.lethality_case   = "FSL";

% set parameters
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.s_t      = 1.0;
N               = 2*(10^4); 

params.lethality_type   = "LA";


releaseType             = "BSR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

% study invasion over range of fitness costs
s_aVec                  = 0:0.01:0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
beta            = 0.5;
alpha           = 0.0321302;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

% [AA_larvae_males, Aa_larvae_males,        2
%  Ab_larvae_males, ab_larvae_males,        4
%  aa_larvae_males, bb_larvae_males,        6
%  AA_larvae_females, Aa_larvae_females,    8
%  Ab_larvae_females, ab_larvae_females,    10
%  aa_larvae_females, bb_larvae_females,    12
%  AA_adults_males, Aa_adults_males,        14
%  Ab_adults_males, ab_adults_males,        16
%  aa_adults_males, bb_adults_males,        18
%  AA_adults_females, Aa_adults_females,    20
%  Ab_adults_females, ab_adults_females,    22
%  aa_adults_females, bb_adults_females]    24

wt_I0                   = zeros(1,24);    
wt_I0(19)               = N;
wt_I0(7)                = (1/2)*fInv_LAMBDA;
wt_I0(1)                = (1/2)*fInv_LAMBDA; 
wt_I0(13)               = (m/(2*muZ))*fInv_LAMBDA;
% release transgenic males
release_I0              = wt_I0;


% matrices for storing everything
invasionThresh          = zeros(4,length(s_aVec)); 
controlEquil            = zeros(4,length(s_aVec)); 

tmax = 5*365;
tVec = 0:1:tmax; 
% study release thresholds up to 12
rFlagMax = 12;     
%%% counter = 1;

for i = 1:length(s_aVec) 
    % loading message
    fprintf("Running sim %.0f of %.0f: s = %0.2f.\n", i, length(s_aVec), ...
        s_aVec(i));

    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;
    %%% releaseMult = 1.5; 

    % set s_a
    params.s_a = s_aVec(i);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 0.5) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    params.beta             = 0.5;
    params.alpha            = 0.0321302;

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24    

    % check allelic frequency
    ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
    ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(16) = midFlag*wt_I0(13);
            elseif releaseType == "BSR"
                release_I0(16) = midFlag*wt_I0(13);
                release_I0(22) = midFlag*wt_I0(19);
            elseif releaseType == "FOR"
                release_I0(22) = midFlag*wt_I0(19);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(1,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(1,i)   = sum(xout(end,19:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(1,i)   = nan; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 1.0) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 2;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;

    % set s_a
    params.s_a = s_aVec(i);    

    params.beta             = 1.0;
    params.alpha            = 1.90085*10^(-4);

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24    

    % check allelic frequency
    ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
    ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(16) = midFlag*wt_I0(13);
            elseif releaseType == "BSR"
                release_I0(16) = midFlag*wt_I0(13);
                release_I0(22) = midFlag*wt_I0(19);
            elseif releaseType == "FOR"
                release_I0(22) = midFlag*wt_I0(19);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,19:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 1.5) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 3;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;

    % set s_a
    params.s_a = s_aVec(i);    

    params.beta             = 1.5;
    params.alpha            = 1.12456*10^(-6);

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24    

    % check allelic frequency
    ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
    ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(16) = midFlag*wt_I0(13);
            elseif releaseType == "BSR"
                release_I0(16) = midFlag*wt_I0(13);
                release_I0(22) = midFlag*wt_I0(19);
            elseif releaseType == "FOR"
                release_I0(22) = midFlag*wt_I0(19);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,19:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% GEN LOGARITHMIC %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 4;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;

    % set s_a
    params.s_a = s_aVec(i);    

    params.beta             = 0.9;
    params.alpha            = 0.014543805956894;

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts); 

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24    

    % check allelic frequency
    ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
    ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(16) = midFlag*wt_I0(13);
            elseif releaseType == "BSR"
                release_I0(16) = midFlag*wt_I0(13);
                release_I0(22) = midFlag*wt_I0(19);
            elseif releaseType == "FOR"
                release_I0(22) = midFlag*wt_I0(19);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,19:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (params.lethality_case == "FSL")
    subplot(2,2,1)
    % relative female pops at control equil (FSL)
    plot(s_aVec, controlEquil(1,:),'-k','linewidth',1.5);
    ylim([0,1]);
    xticks([0, 0.05, 0.1, 0.15, 0.2]);
    % xlim([0,0.11])
    % xlim([0,s_aVec(min(find(isnan(controlEquil(1,:))))-1)]);
    title('(a)')
    hold on
    plot(s_aVec, controlEquil(2,:),'--k','linewidth',1.5);
    plot(s_aVec, controlEquil(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, controlEquil(4,:),':k','linewidth',1.5);
    ylabel('relative female pop.','interpreter','latex');
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]);
    xlabel('fitness cost, $c_a$','interpreter','latex');
    legend('gen. logistic, $\beta=0.5$', 'gen. logistic, $\beta=1.0$', ...
        'gen. logistic, $\beta=1.5$', 'logarithmic', ...
        'interpreter','latex');
    set(gca,'fontsize',16)
    
    subplot(2,2,3)
    % invasion threshold (FSL)
    plot(s_aVec, invasionThresh(1,:),'-k','linewidth',1.5);
    ylim([0,12]);
    xticks([0, 0.05, 0.1, 0.15, 0.2]);
    % xlim([0,0.11])
    % xlim([0,s_aVec(min(find(isnan(controlEquil(1,:))))-1)]);
    title('(c)')
    hold on
    plot(s_aVec, invasionThresh(2,:),'--k','linewidth',1.5);
    plot(s_aVec, invasionThresh(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, invasionThresh(4,:),':k','linewidth',1.5);
    ylabel('invasion threshold','interpreter','latex');
    yticks([0, 2, 4, 6, 8, 10, 12]);
    xlabel('fitness cost, $c_a$','interpreter','latex');
    set(gca,'fontsize',16)    
elseif (params.lethality_case == "BSL")
    subplot(2,2,2)
    plot(s_aVec, controlEquil(1,:),'-k','linewidth',1.5);
    ylim([0,1]);
    xticks([0, 0.05, 0.1, 0.15, 0.2]);
    % xlim([0,0.11])
    % xlim([0,s_aVec(min(find(isnan(controlEquil(1,:))))-1)]);
    title('(b)')
    hold on
    plot(s_aVec, controlEquil(2,:),'--k','linewidth',1.5);
    plot(s_aVec, controlEquil(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, controlEquil(4,:),':k','linewidth',1.5);
    ylabel('relative female pop.','interpreter','latex');
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]);
    xlabel('fitness cost, $c_a$','interpreter','latex');
    set(gca,'fontsize',16)
    
    subplot(2,2,4)
    plot(s_aVec, invasionThresh(1,:),'-k','linewidth',1.5);
    ylim([0,12]);
    xticks([0, 0.05, 0.1, 0.15, 0.2]);
    % xlim([0,0.11])
    % xlim([0,s_aVec(min(find(isnan(controlEquil(1,:))))-1)]);
    title('(d)')
    hold on
    plot(s_aVec, invasionThresh(2,:),'--k','linewidth',1.5);
    plot(s_aVec, invasionThresh(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, invasionThresh(4,:),':k','linewidth',1.5);
    ylabel('invasion threshold','interpreter','latex');
    yticks([0, 2, 4, 6, 8, 10, 12]);
    xlabel('fitness cost, $c_a$','interpreter','latex');
    set(gca,'fontsize',16)    
end



%% 2.2. As above, but now for the bi-sex lethal case. Invasion is 
% successful if it occurs within 10 years of release. This is for the 
% single-locus UD case. We assume that the toxin is always fully lethal 
% (s_t = 1). Also, recall that a BSR is required for a one-locus UD. 
% 
% NOTE: ALL METRICS WILL BE THE SAME AS IN THE FSL CASE. 

% close all 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BI-SEX LETHAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set parameters
params          = struct();
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.s_t      = 0.5;
N               = 2*(10^4); 

params.lethality_case   = "BSL";
params.lethality_type   = "LA";


releaseType             = "BSR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

% study invasion over range of fitness costs
s_aVec                  = 0:0.01:0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
beta            = 0.5;
alpha           = 0.0321302;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

% [AA_larvae_males, Aa_larvae_males,        2
%  Ab_larvae_males, ab_larvae_males,        4
%  aa_larvae_males, bb_larvae_males,        6
%  AA_larvae_females, Aa_larvae_females,    8
%  Ab_larvae_females, ab_larvae_females,    10
%  aa_larvae_females, bb_larvae_females,    12
%  AA_adults_males, Aa_adults_males,        14
%  Ab_adults_males, ab_adults_males,        16
%  aa_adults_males, bb_adults_males,        18
%  AA_adults_females, Aa_adults_females,    20
%  Ab_adults_females, ab_adults_females,    22
%  aa_adults_females, bb_adults_females]    24

wt_I0                   = zeros(1,24);    
wt_I0(19)               = N;
wt_I0(7)                = (1/2)*fInv_LAMBDA;
wt_I0(1)                = (1/2)*fInv_LAMBDA; 
wt_I0(13)               = (m/(2*muZ))*fInv_LAMBDA;
% release transgenic males
release_I0              = wt_I0;


% matrices for storing everything
invasionThresh          = zeros(4,length(s_aVec)); 
controlEquil            = zeros(4,length(s_aVec)); 

tmax = 10*365;
tVec = 0:1:tmax; 
% study release thresholds up to 12
rFlagMax = 12;     
%%% counter = 1;

for i = 1:length(s_aVec) 
    % loading message
    fprintf("Running sim %.0f of %.0f: s = %0.2f.\n", i, length(s_aVec), ...
        s_aVec(i));

    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;
    %%% releaseMult = 1.5; 

    % set s_a
    params.s_a = s_aVec(i);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 0.5) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    params.beta             = 0.5;
    params.alpha            = 0.0321302;

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24    

    % check allelic frequency
    ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
    ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(16) = midFlag*wt_I0(13);
            elseif releaseType == "BSR"
                release_I0(16) = midFlag*wt_I0(13);
                release_I0(22) = midFlag*wt_I0(19);
            elseif releaseType == "FOR"
                release_I0(22) = midFlag*wt_I0(19);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(1,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(1,i)   = sum(xout(end,19:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(1,i)   = nan; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 1.0) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 2;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;

    % set s_a
    params.s_a = s_aVec(i);    

    params.beta             = 1.0;
    params.alpha            = 1.90085*10^(-4);

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24    

    % check allelic frequency
    ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
    ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(16) = midFlag*wt_I0(13);
            elseif releaseType == "BSR"
                release_I0(16) = midFlag*wt_I0(13);
                release_I0(22) = midFlag*wt_I0(19);
            elseif releaseType == "FOR"
                release_I0(22) = midFlag*wt_I0(19);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,19:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 1.5) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 3;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;

    % set s_a
    params.s_a = s_aVec(i);    

    params.beta             = 1.5;
    params.alpha            = 1.12456*10^(-6);

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24    

    % check allelic frequency
    ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
    ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(16) = midFlag*wt_I0(13);
            elseif releaseType == "BSR"
                release_I0(16) = midFlag*wt_I0(13);
                release_I0(22) = midFlag*wt_I0(19);
            elseif releaseType == "FOR"
                release_I0(22) = midFlag*wt_I0(19);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,19:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% GEN LOGARITHMIC %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 4;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;

    % set s_a
    params.s_a = s_aVec(i);    

    params.beta             = 0.9;
    params.alpha            = 0.014543805956894;

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts); 

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24    

    % check allelic frequency
    ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
    ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(16) = midFlag*wt_I0(13);
            elseif releaseType == "BSR"
                release_I0(16) = midFlag*wt_I0(13);
                release_I0(22) = midFlag*wt_I0(19);
            elseif releaseType == "FOR"
                release_I0(22) = midFlag*wt_I0(19);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,19:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2)
% relative female pops at control equil (BSL)
plot(s_aVec, controlEquil(1,:),'-k','linewidth',1.5);
title('(b)')
hold on
plot(s_aVec, controlEquil(2,:),'--k','linewidth',1.5);
plot(s_aVec, controlEquil(3,:),'-.k','linewidth',1.5);
plot(s_aVec, controlEquil(4,:),':k','linewidth',1.5);
ylabel('relative female pop. at control equil.','interpreter','latex');
xlabel('fitness cost, $s_a$','interpreter','latex');

subplot(2,2,4)
% invasion threshold (BSL)
plot(s_aVec, invasionThresh(1,:),'-k','linewidth',1.5);
title('(d)')
hold on
plot(s_aVec, invasionThresh(2,:),'--k','linewidth',1.5);
plot(s_aVec, invasionThresh(3,:),'-.k','linewidth',1.5);
plot(s_aVec, invasionThresh(4,:),':k','linewidth',1.5);
ylabel('invasion threshold','interpreter','latex');
xlabel('fitness cost, $s_a$','interpreter','latex');    

%% 3.1. For a 2L UD case, determination of release ratios for successful 
% invasion for each case. Invasion is successful if it occurs within 10 
% years of release. We assume that the toxin is always fully lethal 
% (s_t = 1). For a 2L UD, a male-only release occurs. This is for a
% female-sex lethal drive with strongly suppressed alleles. 
%
% NOTE: This section can simulate both FSL and BSL; plot at end can handle
% both, assuming BSL results are plotted after in a 2x2 fig. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% SET FSL OR BSL %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params                  = struct();
params.lethality_case   = "FSL";

% set parameters
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.s_t      = 1.0;
N               = 2*(10^4); 

params.lethality_type   = "LA";
params.supp             = "SS";


releaseType             = "BSR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

% study invasion over range of fitness costs
s_aVec                  = 0:0.01:0.2; % 0.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
beta            = 0.5;
alpha           = 0.0321302;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

% [AABB_larvae_male, AABB_larvae_female    2
%  AABb_larvae_male, AABb_larvae_female    4
%  AAbb_larvae_male, AAbb_larvae_female    6
%  AaBB_larvae_male, AaBB_larvae_female    8
%  AaBb_larvae_male, AaBb_larvae_female    10
%  Aabb_larvae_male, Aabb_larvae_female    12
%  aaBB_larvae_male, aaBB_larvae_female    14
%  aaBb_larvae_male, aaBb_larvae_female    16
%  aabb_larvae_male, aabb_larvae_female    18
%  AABB_adult_males, AABb_adult_males      20
%  AAbb_adult_males, AaBB_adult_males      22
%  AaBb_adult_males, Aabb_adult_males      24
%  aaBB_adult_males, aaBb_adult_males      26
%  aabb_adult_males, AABB_adult_females    28
%  AABb_adult_females, AAbb_adult_females  30
%  AaBB_adult_females, AaBb_adult_females  32
%  Aabb_adult_females, aaBB_adult_females  34
%  aaBb_adult_females, aabb_adult_females] 36

wt_I0                   = zeros(1,36);    
wt_I0(28)               = N;
wt_I0(1)                = (1/2)*fInv_LAMBDA;
wt_I0(2)                = (1/2)*fInv_LAMBDA; 
wt_I0(19)               = (m/(2*muZ))*fInv_LAMBDA;

% release transgenic males
release_I0              = wt_I0;

% matrices for storing everything
invasionThresh          = zeros(4,length(s_aVec)); 
controlEquil            = zeros(4,length(s_aVec)); 

tmax = 5*365;
tVec = 0:1:tmax; 
% study release thresholds up to 12
rFlagMax = 100;     
%%% counter = 1;

for i = 1:length(s_aVec) 
    % loading message
    fprintf("Running sim %.0f of %.0f: s = %0.2f.\n", i, length(s_aVec), ...
        s_aVec(i));

    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;
    %%% releaseMult = 1.5; 

    % set s_a
    params.s_a = s_aVec(i);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 0.5) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    params.beta             = 0.5;
    params.alpha            = 0.0321302;

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 

    % [AABB_larvae_male, AABB_larvae_female    2
    %  AABb_larvae_male, AABb_larvae_female    4
    %  AAbb_larvae_male, AAbb_larvae_female    6
    %  AaBB_larvae_male, AaBB_larvae_female    8
    %  AaBb_larvae_male, AaBb_larvae_female    10
    %  Aabb_larvae_male, Aabb_larvae_female    12
    %  aaBB_larvae_male, aaBB_larvae_female    14
    %  aaBb_larvae_male, aaBb_larvae_female    16
    %  aabb_larvae_male, aabb_larvae_female    18
    %  AABB_adult_males, AABb_adult_males      20
    %  AAbb_adult_males, AaBB_adult_males      22
    %  AaBb_adult_males, Aabb_adult_males      24
    %  aaBB_adult_males, aaBb_adult_males      26
    %  aabb_adult_males, AABB_adult_females    28
    %  AABb_adult_females, AAbb_adult_females  30
    %  AaBB_adult_females, AaBb_adult_females  32
    %  Aabb_adult_females, aaBB_adult_females  34
    %  aaBb_adult_females, aabb_adult_females] 36  

    % check allelic frequency
    ind1 = [3,4,7,8,20,22,29,31];
        % one copy of a or b
    ind2 = [5,6,9,10,13,14,21,23,25,30,32,34];
        % two copies of a or b
    ind3 = [11,12,15,16,24,26,33,35];
        % three copies of a or b
    ind4 = [17,18,27,36];
        % four copies of a or b        

    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
        4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
        4*sum(xout(:,ind4),2))./(4*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    %%%
    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
            4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
            4*sum(xout(:,ind4),2))./(4*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(27) = midFlag*wt_I0(19);
            elseif releaseType == "BSR"
                release_I0(27) = midFlag*wt_I0(19);
                release_I0(36) = midFlag*wt_I0(28);
            elseif releaseType == "FOR"
                release_I0(36) = midFlag*wt_I0(28);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                4*sum(xout(:,ind4),2))./(4*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                    4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                    4*sum(xout(:,ind4),2))./(4*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(1,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(1,i)   = sum(xout(end,28:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(1,i)   = nan; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 1.0) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 2;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;
    %%% releaseMult = 1.5; 

    % set s_a
    params.s_a = s_aVec(i);   

    params.beta             = 1.0;
    params.alpha            = 1.90085*10^(-4);

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 

    % [AABB_larvae_male, AABB_larvae_female    2
    %  AABb_larvae_male, AABb_larvae_female    4
    %  AAbb_larvae_male, AAbb_larvae_female    6
    %  AaBB_larvae_male, AaBB_larvae_female    8
    %  AaBb_larvae_male, AaBb_larvae_female    10
    %  Aabb_larvae_male, Aabb_larvae_female    12
    %  aaBB_larvae_male, aaBB_larvae_female    14
    %  aaBb_larvae_male, aaBb_larvae_female    16
    %  aabb_larvae_male, aabb_larvae_female    18
    %  AABB_adult_males, AABb_adult_males      20
    %  AAbb_adult_males, AaBB_adult_males      22
    %  AaBb_adult_males, Aabb_adult_males      24
    %  aaBB_adult_males, aaBb_adult_males      26
    %  aabb_adult_males, AABB_adult_females    28
    %  AABb_adult_females, AAbb_adult_females  30
    %  AaBB_adult_females, AaBb_adult_females  32
    %  Aabb_adult_females, aaBB_adult_females  34
    %  aaBb_adult_females, aabb_adult_females] 36  

    % check allelic frequency
    ind1 = [3,4,7,8,20,22,29,31];
        % one copy of a or b
    ind2 = [5,6,9,10,13,14,21,23,25,30,32,34];
        % two copies of a or b
    ind3 = [11,12,15,16,24,26,33,35];
        % three copies of a or b
    ind4 = [17,18,27,36];
        % four copies of a or b        

    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
        4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
        4*sum(xout(:,ind4),2))./(4*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    %%%
    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
            4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
            4*sum(xout(:,ind4),2))./(4*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(27) = midFlag*wt_I0(19);
            elseif releaseType == "BSR"
                release_I0(27) = midFlag*wt_I0(19);
                release_I0(36) = midFlag*wt_I0(28);
            elseif releaseType == "FOR"
                release_I0(36) = midFlag*wt_I0(28);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                4*sum(xout(:,ind4),2))./(4*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                    4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                    4*sum(xout(:,ind4),2))./(4*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,28:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 1.5) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 3;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;
    %%% releaseMult = 1.5; 

    % set s_a
    params.s_a = s_aVec(i);   

    params.beta             = 1.5;
    params.alpha            = 1.12456*10^(-6);

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 

    % [AABB_larvae_male, AABB_larvae_female    2
    %  AABb_larvae_male, AABb_larvae_female    4
    %  AAbb_larvae_male, AAbb_larvae_female    6
    %  AaBB_larvae_male, AaBB_larvae_female    8
    %  AaBb_larvae_male, AaBb_larvae_female    10
    %  Aabb_larvae_male, Aabb_larvae_female    12
    %  aaBB_larvae_male, aaBB_larvae_female    14
    %  aaBb_larvae_male, aaBb_larvae_female    16
    %  aabb_larvae_male, aabb_larvae_female    18
    %  AABB_adult_males, AABb_adult_males      20
    %  AAbb_adult_males, AaBB_adult_males      22
    %  AaBb_adult_males, Aabb_adult_males      24
    %  aaBB_adult_males, aaBb_adult_males      26
    %  aabb_adult_males, AABB_adult_females    28
    %  AABb_adult_females, AAbb_adult_females  30
    %  AaBB_adult_females, AaBb_adult_females  32
    %  Aabb_adult_females, aaBB_adult_females  34
    %  aaBb_adult_females, aabb_adult_females] 36  

    % check allelic frequency
    ind1 = [3,4,7,8,20,22,29,31];
        % one copy of a or b
    ind2 = [5,6,9,10,13,14,21,23,25,30,32,34];
        % two copies of a or b
    ind3 = [11,12,15,16,24,26,33,35];
        % three copies of a or b
    ind4 = [17,18,27,36];
        % four copies of a or b        

    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
        4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
        4*sum(xout(:,ind4),2))./(4*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    %%%
    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
            4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
            4*sum(xout(:,ind4),2))./(4*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(27) = midFlag*wt_I0(19);
            elseif releaseType == "BSR"
                release_I0(27) = midFlag*wt_I0(19);
                release_I0(36) = midFlag*wt_I0(28);
            elseif releaseType == "FOR"
                release_I0(36) = midFlag*wt_I0(28);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                4*sum(xout(:,ind4),2))./(4*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                    4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                    4*sum(xout(:,ind4),2))./(4*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,28:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% LOGARITHMIC %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 4;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;
    %%% releaseMult = 1.5; 

    % set s_a
    params.s_a = s_aVec(i);   

    params.beta             = 0.9;
    params.alpha            = 0.014543805956894;

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end

    [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts); 

    % [AABB_larvae_male, AABB_larvae_female    2
    %  AABb_larvae_male, AABb_larvae_female    4
    %  AAbb_larvae_male, AAbb_larvae_female    6
    %  AaBB_larvae_male, AaBB_larvae_female    8
    %  AaBb_larvae_male, AaBb_larvae_female    10
    %  Aabb_larvae_male, Aabb_larvae_female    12
    %  aaBB_larvae_male, aaBB_larvae_female    14
    %  aaBb_larvae_male, aaBb_larvae_female    16
    %  aabb_larvae_male, aabb_larvae_female    18
    %  AABB_adult_males, AABb_adult_males      20
    %  AAbb_adult_males, AaBB_adult_males      22
    %  AaBb_adult_males, Aabb_adult_males      24
    %  aaBB_adult_males, aaBb_adult_males      26
    %  aabb_adult_males, AABB_adult_females    28
    %  AABb_adult_females, AAbb_adult_females  30
    %  AaBB_adult_females, AaBb_adult_females  32
    %  Aabb_adult_females, aaBB_adult_females  34
    %  aaBb_adult_females, aabb_adult_females] 36  

    % check allelic frequency
    ind1 = [3,4,7,8,20,22,29,31];
        % one copy of a or b
    ind2 = [5,6,9,10,13,14,21,23,25,30,32,34];
        % two copies of a or b
    ind3 = [11,12,15,16,24,26,33,35];
        % three copies of a or b
    ind4 = [17,18,27,36];
        % four copies of a or b        

    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
        4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
        4*sum(xout(:,ind4),2))./(4*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    %%%
    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
            4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
            4*sum(xout(:,ind4),2))./(4*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(27) = midFlag*wt_I0(19);
            elseif releaseType == "BSR"
                release_I0(27) = midFlag*wt_I0(19);
                release_I0(36) = midFlag*wt_I0(28);
            elseif releaseType == "FOR"
                release_I0(36) = midFlag*wt_I0(28);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                4*sum(xout(:,ind4),2))./(4*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                    4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                    4*sum(xout(:,ind4),2))./(4*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,28:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (params.lethality_case == "FSL")
    % FSL case
    subplot(2,2,1)
    % relative female pops at control equil (FSL)
    plot(s_aVec, controlEquil(1,:),'-k','linewidth',1.5);
    %%% ylim([0,1])
    %%% xlim([0,0.1])
    % xlim([0,s_aVec(min(find(isnan(controlEquil(1,:))))-1)])
    title('(a)')
    hold on
    plot(s_aVec, controlEquil(2,:),'--k','linewidth',1.5);
    plot(s_aVec, controlEquil(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, controlEquil(4,:),':k','linewidth',1.5);
    ylabel('relative female pop.','interpreter','latex');
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]);    
    xticks([0, 0.02, 0.04, 0.06, 0.08, 0.1]);
    xlabel('fitness cost, $c_a$','interpreter','latex');
    set(gca,'fontsize',16)
    
    subplot(2,2,3)
    % invasion threshold (FSL)
    plot(s_aVec, invasionThresh(1,:),'-k','linewidth',1.5);
    % xlim([0,s_aVec(min(find(isnan(controlEquil(1,:))))-1)]);
    %%% xlim([0,0.1])
    % ylim([0,12]);
    title('(c)')
    hold on
    plot(s_aVec, invasionThresh(2,:),'--k','linewidth',1.5);
    plot(s_aVec, invasionThresh(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, invasionThresh(4,:),':k','linewidth',1.5);
    ylabel('invasion threshold','interpreter','latex');
    % yticks([0, 2, 4, 6, 8, 10, 12]); 
    yticks([0, 5, 10, 15, 20, 25, 30]); 
    xlabel('fitness cost, $c_a$','interpreter','latex'); 
    set(gca,'fontsize',16)
elseif (params.lethality_case == "BSL")
    figure
    % BSL
    subplot(2,2,2)
    plot(s_aVec, controlEquil(1,:),'-k','linewidth',1.5);
    %%% ylim([0,1])
    % xlim([0,s_aVec(min(find(isnan(controlEquil(1,:))))-1)])
    %%% xlim([0,0.1])
    title('(b)')
    hold on
    plot(s_aVec, controlEquil(2,:),'--k','linewidth',1.5);
    plot(s_aVec, controlEquil(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, controlEquil(4,:),':k','linewidth',1.5);
    ylabel('relative female pop.','interpreter','latex');
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); 
    xlabel('fitness cost, $c_a$','interpreter','latex');
    set(gca,'fontsize',16)
    
    subplot(2,2,4)
    % invasion threshold (BSL)
    plot(s_aVec, invasionThresh(1,:),'-k','linewidth',1.5);
    % xlim([0,s_aVec(min(find(isnan(controlEquil(1,:))))-1)]);
    %%% xlim([0,0.1])
    %%% ylim([0,12]);
    title('(d)')
    hold on
    plot(s_aVec, invasionThresh(2,:),'--k','linewidth',1.5);
    plot(s_aVec, invasionThresh(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, invasionThresh(4,:),':k','linewidth',1.5);
    yticks([0, 2, 4, 6, 8, 10, 12]); 
    % yticks([0, 5, 10, 15, 20, 25, 30]); 
    ylabel('invasion threshold','interpreter','latex');
    xlabel('fitness cost, $c_a$','interpreter','latex'); 
    set(gca,'fontsize',16)
else
    error("Unrecognized lethality case; did you mean BSL or FSL?")
end
%% 3.2. Comparison of BSL and FSL system performance for the 2L UD 
% 
% This script assumes that the above section was run and the matrices
% controlEquil and invasionThresh have been saved under appropriate names,
% i.e. "controlEquil_FSL" and "invasionThresh_FSL", etc. This script is NOT
% standalone and requires some variables to exist in memory, such as
% s_aVec.

% comparing FSL performance to BSL--black is BSL, red is FSL
subplot(1,2,1)
plot(s_aVec, controlEquil_BSL(1,:),'-k','linewidth',1.5);
title('(a) control equil comp')
hold on
ylim([0,1])
% xlim([0,s_aVec(min(find(isnan(controlEquil_BSL(1,:))))-1)])
plot(s_aVec, controlEquil_BSL(2,:),'--k','linewidth',1.5);
plot(s_aVec, controlEquil_BSL(3,:),'-.k','linewidth',1.5);
plot(s_aVec, controlEquil_BSL(4,:),':k','linewidth',1.5);
% now plot FSL
plot(s_aVec, controlEquil_FSL(1,:),'-r','linewidth',1.5);
plot(s_aVec, controlEquil_FSL(2,:),'--r','linewidth',1.5);
plot(s_aVec, controlEquil_FSL(3,:),'-.r','linewidth',1.5);
plot(s_aVec, controlEquil_FSL(4,:),':r','linewidth',1.5);
ylabel('relative female pop.','interpreter','latex');
yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); 
xlabel('fitness cost, $c_a$','interpreter','latex');
set(gca,'fontsize',16)

% invasion thresh
subplot(1,2,2)
plot(s_aVec, invasionThresh_BSL(1,:),'-k','linewidth',1.5);
title('(b) invasion thresh comp')
hold on
ylim([0,12])
% xlim([0,s_aVec(min(find(isnan(invasionThresh_BSL(1,:))))-1)])
plot(s_aVec, invasionThresh_BSL(2,:),'--k','linewidth',1.5);
plot(s_aVec, invasionThresh_BSL(3,:),'-.k','linewidth',1.5);
plot(s_aVec, invasionThresh_BSL(4,:),':k','linewidth',1.5);
% now plot FSL
plot(s_aVec, invasionThresh_FSL(1,:),'-r','linewidth',1.5);
plot(s_aVec, invasionThresh_FSL(2,:),'--r','linewidth',1.5);
plot(s_aVec, invasionThresh_FSL(3,:),'-.r','linewidth',1.5);
plot(s_aVec, invasionThresh_FSL(4,:),':r','linewidth',1.5);
yticks([0, 2, 4, 6, 8, 10, 12]); 
ylabel('invasion threshold','interpreter','latex');
xlabel('fitness cost, $c_a$','interpreter','latex'); 
set(gca,'fontsize',16)

%% 3.3. Introgression is much more successful in the BSL vs. the FSL case. 
% Hypothesize that in the BSL scenario, males must have complementary
% transgenes, so that allele frequency is greater in BSL case. In the FSL
% case, wild-type alleles can "lurk" in males a lot longer. 

% how do allele frequencies compare between two simulations? 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% BSL %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set parameters
params          = struct();
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.s_t      = 1.0;
N               = 2*(10^4); 

params.lethality_case   = "BSL";
params.lethality_type   = "LA";
params.supp             = "SS";


releaseType             = "MOR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

% [AABB_larvae_male, AABB_larvae_female    2
%  AABb_larvae_male, AABb_larvae_female    4
%  AAbb_larvae_male, AAbb_larvae_female    6
%  AaBB_larvae_male, AaBB_larvae_female    8
%  AaBb_larvae_male, AaBb_larvae_female    10
%  Aabb_larvae_male, Aabb_larvae_female    12
%  aaBB_larvae_male, aaBB_larvae_female    14
%  aaBb_larvae_male, aaBb_larvae_female    16
%  aabb_larvae_male, aabb_larvae_female    18
%  AABB_adult_males, AABb_adult_males      20
%  AAbb_adult_males, AaBB_adult_males      22
%  AaBb_adult_males, Aabb_adult_males      24
%  aaBB_adult_males, aaBb_adult_males      26
%  aabb_adult_males, AABB_adult_females    28
%  AABb_adult_females, AAbb_adult_females  30
%  AaBB_adult_females, AaBb_adult_females  32
%  Aabb_adult_females, aaBB_adult_females  34
%  aaBb_adult_females, aabb_adult_females] 36

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
beta            = 0.5;
alpha           = 0.0321302;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

wt_I0                   = zeros(1,36);    
wt_I0(28)               = N;
wt_I0(1)                = (1/2)*fInv_LAMBDA;
wt_I0(2)                = (1/2)*fInv_LAMBDA; 
wt_I0(19)               = (m/(2*muZ))*fInv_LAMBDA;

% release transgenic males
release_I0              = wt_I0;

tmax = 2*365;
tVec = 0:1:tmax; 
% study release thresholds up to 12
rFlag = 8;

% set s_a
params.s_a = 0.05;

params.beta             = 0.9;
params.alpha            = 0.014543805956894;

% adjusted release based on releaseType
if releaseType == "MOR"
    release_I0(27) = rFlag*wt_I0(19);
elseif releaseType == "BSR"
    release_I0(27) = rFlag*wt_I0(19);
    release_I0(36) = rFlag*wt_I0(28);
elseif releaseType == "FOR"
    release_I0(36) = rFlag*wt_I0(28);
end

[~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts);

% check allele frequency
ind1 = [3,4,7,8,20,22,29,31];
    % one copy of a or b
ind2 = [5,6,9,10,13,14,21,23,25,30,32,34];
    % two copies of a or b
ind3 = [11,12,15,16,24,26,33,35];
    % three copies of a or b
ind4 = [17,18,27,36];
    % four copies of a or b      

allele_freq_BSL = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
    4*sum(xout(:,ind4),2))./(4*sum(xout,2));

plot(allele_freq_BSL,'-k','linewidth',1.5)
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% FSL %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.lethality_case   = "FSL";

[~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts);    

allele_freq_FSL = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
    4*sum(xout(:,ind4),2))./(4*sum(xout,2));

plot(allele_freq_FSL,'--r','linewidth',1.5)
ylim([0,1])
xlim([0,tmax])
xlabel('time','interpreter','latex');
ylabel('allele frequency','interpreter','latex'); 
set(gca,'fontsize',16)

%% 3.4. Comparing BSL vs FSL for two locus case when release ratios are the 
% same... does BSL perform better at suppression after 4 years, as others 
% have found?

s_aVec = 0:0.01:0.1;
suppBSL = zeros(1,length(s_aVec));
suppFSL = zeros(1,length(s_aVec));

% set parameters
params          = struct();
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.s_t      = 1.0;
N               = 2*(10^4); 

params.lethality_case   = "BSL";
params.lethality_type   = "LA";
params.supp             = "SS";


releaseType             = "MOR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

% [AABB_larvae_male, AABB_larvae_female    2
%  AABb_larvae_male, AABb_larvae_female    4
%  AAbb_larvae_male, AAbb_larvae_female    6
%  AaBB_larvae_male, AaBB_larvae_female    8
%  AaBb_larvae_male, AaBb_larvae_female    10
%  Aabb_larvae_male, Aabb_larvae_female    12
%  aaBB_larvae_male, aaBB_larvae_female    14
%  aaBb_larvae_male, aaBb_larvae_female    16
%  aabb_larvae_male, aabb_larvae_female    18
%  AABB_adult_males, AABb_adult_males      20
%  AAbb_adult_males, AaBB_adult_males      22
%  AaBb_adult_males, Aabb_adult_males      24
%  aaBB_adult_males, aaBb_adult_males      26
%  aabb_adult_males, AABB_adult_females    28
%  AABb_adult_females, AAbb_adult_females  30
%  AaBB_adult_females, AaBb_adult_females  32
%  Aabb_adult_females, aaBB_adult_females  34
%  aaBb_adult_females, aabb_adult_females] 36

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
params.beta             = 0.9;
params.alpha            = 0.014543805956894;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

wt_I0                   = zeros(1,36);    
wt_I0(28)               = N;
wt_I0(1)                = (1/2)*fInv_LAMBDA;
wt_I0(2)                = (1/2)*fInv_LAMBDA; 
wt_I0(19)               = (m/(2*muZ))*fInv_LAMBDA;

% release transgenic males
release_I0              = wt_I0;

tmax = 5*365;
tVec = 0:1:tmax; 

for i = 1:length(s_aVec)
    % set s_a
    params.s_a = s_aVec(i);

    % BSL
    params.lethality_case   = "BSL";
    % study release thresholds up to 12
    rFlag = invasionThresh_BSL(4,i);        
    release_I0(27) = rFlag*wt_I0(19);
    
    [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts);

    % store female pop (BSL)
    suppBSL(i) = sum(xout(end,28:end))/(2*10^4);

    % FSL
    params.lethality_case   = "FSL";

    rFlag = invasionThresh_FSL(4,i);
    release_I0(27) = rFlag*wt_I0(19);

    [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts);    
    % store female pop (FSL)
    suppFSL(i) = sum(xout(end,28:end))/(2*10^4);

end

figure
plot(s_aVec,suppBSL,'-k','linewidth',1.5)
hold on
plot(s_aVec,suppFSL,'--r','linewidth',1.5)
set(gca,'fontsize',16)



%% 4.1. Section recreating the figures in Khamis to make sure the functions 
% are behaving correctly. (This section recreates Fig 2.)

%%%%%%%%%%%%%%%%%
%%%% Fig 2a %%%%%
%%%%%%%%%%%%%%%%%

params                  = struct();
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 16;
params.muX      = 0.03;
params.muZ      = log(10/9);
params.muY      = log(10/9);
params.m        = 0.1;
params.s_t      = 0.99;
params.s_a      = 0.0;
N               = 2*(10^3); 

params.lethality_type   = "EA";
params.supp             = "SS";

releaseType             = "BSR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

beta            = 0.9;
alpha           = (m/(2*2*1000*muZ))*((exp(lambda*m/(2*muZ)-m-muX) - 1)^(1/beta));
fInv            = @(x) (1/alpha)*((exp(x)-1).^(1/beta));
fInv_LAMBDA     = fInv(LAMBDA); 

params.beta     = beta;
params.alpha    = alpha;                

wt_I0                   = zeros(1,24);    
wt_I0(19)               = N;
wt_I0(7)                = (1/2)*fInv_LAMBDA;
wt_I0(1)                = (1/2)*fInv_LAMBDA; 
wt_I0(13)               = (m/(2*muZ))*fInv_LAMBDA;

% check that the system is at equilibrium prior to any releases...
% SUCCESS
% tVec = 0:365; 
% [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,wt_I0,opts); 
% plot(xout)

releaseVec      = 0.63:0.001:0.69;
outputMat       = zeros(5,length(releaseVec));
dayVec          = [90,180,360,720,1440];

for i = 1:length(dayVec)
    disp(i); 
    tVec = 0:dayVec(i);
    for j = 1:length(releaseVec)
        release_I0          = wt_I0;
        release_ratio       = releaseVec(j)/(1-releaseVec(j));

        if releaseType == "MOR"
            release_I0(16) = release_ratio*wt_I0(13);
        elseif releaseType == "BSR"
            release_I0(16) = release_ratio*wt_I0(13);
            release_I0(22) = release_ratio*wt_I0(19);
        elseif releaseType == "FOR"
            release_I0(22) = release_ratio*wt_I0(19);
        end

        [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts); 
        % calculate drive success metric
        outputMat(i,j) = (xout(end,19) - xout(end,22))/sum(xout(end,19:end));

        % [AA_larvae_males, Aa_larvae_males,        2
        %  Ab_larvae_males, ab_larvae_males,        4
        %  aa_larvae_males, bb_larvae_males,        6
        %  AA_larvae_females, Aa_larvae_females,    8
        %  Ab_larvae_females, ab_larvae_females,    10
        %  aa_larvae_females, bb_larvae_females,    12
        %  AA_adults_males, Aa_adults_males,        14
        %  Ab_adults_males, ab_adults_males,        16
        %  aa_adults_males, bb_adults_males,        18
        %  AA_adults_females, Aa_adults_females,    20
        %  Ab_adults_females, ab_adults_females,    22
        %  aa_adults_females, bb_adults_females]    24

    end
end

subplot(1,3,1)
plot(releaseVec, outputMat,'linewidth',1.5)
xlim([0.63,0.69])
yticks([-1.0, -0.5, 0, 0.5, 1.0])
grid on
title('(a)')
set(gca, 'fontsize',16)

%%%%%%%%%%%%%%%%%
%%%% Fig 2b %%%%%
%%%%%%%%%%%%%%%%%

params                  = struct();
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 16;
params.muX      = 0.03;
params.muZ      = log(10/9);
params.muY      = log(10/9);
params.m        = 0.1;
params.s_t      = 0.99;
params.s_a      = 0.0;
N               = 2*(10^3); 

params.lethality_type   = "EA";
params.supp             = "WS";

releaseType             = "BSR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

beta            = 0.9;
alpha           = (m/(2*2*1000*muZ))*((exp(lambda*m/(2*muZ)-m-muX) - 1)^(1/beta));
fInv            = @(x) (1/alpha)*((exp(x)-1).^(1/beta));
fInv_LAMBDA     = fInv(LAMBDA); 

params.beta     = beta;
params.alpha    = alpha;                

wt_I0                   = zeros(1,36);    
wt_I0(28)               = N;
wt_I0(1)                = (1/2)*fInv_LAMBDA;
wt_I0(2)                = (1/2)*fInv_LAMBDA; 
wt_I0(19)               = (m/(2*muZ))*fInv_LAMBDA;

% release transgenic males
release_I0              = wt_I0;

% check that the system is at equilibrium prior to any releases...
% SUCCESS
% tVec = 0:365; 
% [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,wt_I0,opts); 
% plot(xout)

releaseVec      = 0.46:0.001:0.52;
outputMat       = zeros(5,length(releaseVec));
dayVec          = [90,180,360,720,1440];

for i = 1:length(dayVec)
    disp(i); 
    tVec = 0:dayVec(i);
    for j = 1:length(releaseVec)
        release_I0          = wt_I0;
        release_ratio       = releaseVec(j)/(1-releaseVec(j));

        % adjusted release based on releaseType
        if releaseType == "MOR"
            release_I0(27) = release_ratio*wt_I0(19);
        elseif releaseType == "BSR"
            release_I0(27) = release_ratio*wt_I0(19);
            release_I0(36) = release_ratio*wt_I0(28);
        elseif releaseType == "FOR"
            release_I0(36) = release_ratio*wt_I0(28);
        end

        [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts); 
        % calculate drive success metric
        outputMat(i,j) = (xout(end,28) - xout(end,36))/sum(xout(end,28:end));

        % [AABB_larvae_male, AABB_larvae_female    2
        %  AABb_larvae_male, AABb_larvae_female    4
        %  AAbb_larvae_male, AAbb_larvae_female    6
        %  AaBB_larvae_male, AaBB_larvae_female    8
        %  AaBb_larvae_male, AaBb_larvae_female    10
        %  Aabb_larvae_male, Aabb_larvae_female    12
        %  aaBB_larvae_male, aaBB_larvae_female    14
        %  aaBb_larvae_male, aaBb_larvae_female    16
        %  aabb_larvae_male, aabb_larvae_female    18
        %  AABB_adult_males, AABb_adult_males      20
        %  AAbb_adult_males, AaBB_adult_males      22
        %  AaBb_adult_males, Aabb_adult_males      24
        %  aaBB_adult_males, aaBb_adult_males      26
        %  aabb_adult_males, AABB_adult_females    28
        %  AABb_adult_females, AAbb_adult_females  30
        %  AaBB_adult_females, AaBb_adult_females  32
        %  Aabb_adult_females, aaBB_adult_females  34
        %  aaBb_adult_females, aabb_adult_females] 36

    end
end

subplot(1,3,2)
plot(releaseVec, outputMat,'linewidth',1.5)
xlim([0.46,0.52])
yticks([-1.0, -0.5, 0, 0.5, 1.0])
grid on
title('(b)')
set(gca, 'fontsize',16)

%%%%%%%%%%%%%%%%%
%%%% Fig 2c %%%%%
%%%%%%%%%%%%%%%%%

params                  = struct();
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 16;
params.muX      = 0.03;
params.muZ      = log(10/9);
params.muY      = log(10/9);
params.m        = 0.1;
params.s_t      = 0.99;
params.s_a      = 0.0;
N               = 2*(10^3); 

params.lethality_type   = "EA";
params.supp             = "SS";

releaseType             = "BSR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

beta            = 0.9;
alpha           = (m/(2*2*1000*muZ))*((exp(lambda*m/(2*muZ)-m-muX) - 1)^(1/beta));
fInv            = @(x) (1/alpha)*((exp(x)-1).^(1/beta));
fInv_LAMBDA     = fInv(LAMBDA); 

params.beta     = beta;
params.alpha    = alpha;                

wt_I0                   = zeros(1,36);    
wt_I0(28)               = N;
wt_I0(1)                = (1/2)*fInv_LAMBDA;
wt_I0(2)                = (1/2)*fInv_LAMBDA; 
wt_I0(19)               = (m/(2*muZ))*fInv_LAMBDA;

% release transgenic males
release_I0              = wt_I0;

% check that the system is at equilibrium prior to any releases...
% SUCCESS
% tVec = 0:365; 
% [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,wt_I0,opts); 
% plot(xout)

releaseVec      = 0.23:0.001:0.29;
outputMat       = zeros(5,length(releaseVec));
dayVec          = [90,180,360,720,1440];

for i = 1:length(dayVec)
    disp(i); 
    tVec = 0:dayVec(i);
    for j = 1:length(releaseVec)
        release_I0          = wt_I0;
        release_ratio       = releaseVec(j)/(1-releaseVec(j));

        % adjusted release based on releaseType
        if releaseType == "MOR"
            release_I0(27) = release_ratio*wt_I0(19);
        elseif releaseType == "BSR"
            release_I0(27) = release_ratio*wt_I0(19);
            release_I0(36) = release_ratio*wt_I0(28);
        elseif releaseType == "FOR"
            release_I0(36) = release_ratio*wt_I0(28);
        end

        [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts); 
        % calculate drive success metric
        outputMat(i,j) = (xout(end,28) - xout(end,36))/sum(xout(end,28:end));

        % [AABB_larvae_male, AABB_larvae_female    2
        %  AABb_larvae_male, AABb_larvae_female    4
        %  AAbb_larvae_male, AAbb_larvae_female    6
        %  AaBB_larvae_male, AaBB_larvae_female    8
        %  AaBb_larvae_male, AaBb_larvae_female    10
        %  Aabb_larvae_male, Aabb_larvae_female    12
        %  aaBB_larvae_male, aaBB_larvae_female    14
        %  aaBb_larvae_male, aaBb_larvae_female    16
        %  aabb_larvae_male, aabb_larvae_female    18
        %  AABB_adult_males, AABb_adult_males      20
        %  AAbb_adult_males, AaBB_adult_males      22
        %  AaBb_adult_males, Aabb_adult_males      24
        %  aaBB_adult_males, aaBb_adult_males      26
        %  aabb_adult_males, AABB_adult_females    28
        %  AABb_adult_females, AAbb_adult_females  30
        %  AaBB_adult_females, AaBb_adult_females  32
        %  Aabb_adult_females, aaBB_adult_females  34
        %  aaBb_adult_females, aabb_adult_females] 36

    end
end

subplot(1,3,3)
plot(releaseVec, outputMat,'linewidth',1.5)
xlim([0.23,0.29])
yticks([-1.0, -0.5, 0, 0.5, 1.0])
grid on
legend('$t_N = 90$','$t_N = 180$','$t_N = 360$','$t_N = 720$','$t_N = 1440$',...
    'location','best','interpreter','latex')
title('(c)')
set(gca, 'fontsize',16)

%% 4.2.1. Recreating Fig. 3 in Khamis, or just studying the frequency of transgene 
% homozygous females relative to the wild-type female equilibrium as
% initial release frequency is varied. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% First part studies a bi-sex release for each UD case study. %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,1,1)

releaseVec      = 0.0:0.001:1.0;         
outputMat       = zeros(3,length(releaseVec));
tVec            = 0:365;        

%%%%%%%%%%%%%%%%%%%%%
%%%% HOMOLOGOUS %%%%%
%%%%%%%%%%%%%%%%%%%%%

matInd                  = 1;
params                  = struct();
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 16;
params.muX      = 0.03;
params.muZ      = log(10/9);
params.muY      = log(10/9);
params.m        = 0.1;
params.s_t      = 0.9;
params.s_a      = 0.05;
N               = 2*(10^3); 

params.lethality_type   = "EA";
params.supp             = "SS";

releaseType             = "BSR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

beta            = 0.9;
alpha           = (m/(2*2*1000*muZ))*((exp(lambda*m/(2*muZ)-m-muX) - 1)^(1/beta));
fInv            = @(x) (1/alpha)*((exp(x)-1).^(1/beta));
fInv_LAMBDA     = fInv(LAMBDA); 

params.beta     = beta;
params.alpha    = alpha;                

wt_I0                   = zeros(1,24);    
wt_I0(19)               = N;
wt_I0(7)                = (1/2)*fInv_LAMBDA;
wt_I0(1)                = (1/2)*fInv_LAMBDA; 
wt_I0(13)               = (m/(2*muZ))*fInv_LAMBDA;

for j = 1:length(releaseVec)
    disp(j)
    release_I0          = wt_I0;
    release_ratio       = releaseVec(j)/(1-releaseVec(j));

    if releaseType == "MOR"
        release_I0(16) = release_ratio*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = release_ratio*wt_I0(13);
        release_I0(22) = release_ratio*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = release_ratio*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts); 
    % calculate drive success metric
    outputMat(matInd,j) = sum(xout(end,19:end))/2000;

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24

end

%%%%%%%%%%%%%%%%%%%%%
% WS NON-HOMOLOGOUS %
%%%%%%%%%%%%%%%%%%%%%

matInd                  = 2;
params                  = struct();
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 16;
params.muX      = 0.03;
params.muZ      = log(10/9);
params.muY      = log(10/9);
params.m        = 0.1;
params.s_t      = 0.9;
params.s_a      = 0.05;
N               = 2*(10^3); 

params.lethality_type   = "EA";
params.supp             = "WS";

releaseType             = "BSR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

beta            = 0.9;
alpha           = (m/(2*2*1000*muZ))*((exp(lambda*m/(2*muZ)-m-muX) - 1)^(1/beta));
fInv            = @(x) (1/alpha)*((exp(x)-1).^(1/beta));
fInv_LAMBDA     = fInv(LAMBDA); 

params.beta     = beta;
params.alpha    = alpha;                

wt_I0                   = zeros(1,36);    
wt_I0(28)               = N;
wt_I0(1)                = (1/2)*fInv_LAMBDA;
wt_I0(2)                = (1/2)*fInv_LAMBDA; 
wt_I0(19)               = (m/(2*muZ))*fInv_LAMBDA;

for j = 1:length(releaseVec)
    release_I0          = wt_I0;
    release_ratio       = releaseVec(j)/(1-releaseVec(j));

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(27) = release_ratio*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = release_ratio*wt_I0(19);
        release_I0(36) = release_ratio*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = release_ratio*wt_I0(28);
    end

    [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts); 
    % freq. of (female) drive homozygotes as fraction of wt equil pop
    outputMat(matInd,j) = sum(xout(end,28:end))/2000;
end

%%%%%%%%%%%%%%%%%%%%%
% SS NON-HOMOLOGOUS %
%%%%%%%%%%%%%%%%%%%%%

matInd                  = 3;
params                  = struct();
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 16;
params.muX      = 0.03;
params.muZ      = log(10/9);
params.muY      = log(10/9);
params.m        = 0.1;
params.s_t      = 0.9;
params.s_a      = 0.05;
N               = 2*(10^3); 

params.lethality_type   = "EA";
params.supp             = "SS";

releaseType             = "BSR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

beta            = 0.9;
alpha           = (m/(2*2*1000*muZ))*((exp(lambda*m/(2*muZ)-m-muX) - 1)^(1/beta));
fInv            = @(x) (1/alpha)*((exp(x)-1).^(1/beta));
fInv_LAMBDA     = fInv(LAMBDA); 

params.beta     = beta;
params.alpha    = alpha;                

wt_I0                   = zeros(1,36);    
wt_I0(28)               = N;
wt_I0(1)                = (1/2)*fInv_LAMBDA;
wt_I0(2)                = (1/2)*fInv_LAMBDA; 
wt_I0(19)               = (m/(2*muZ))*fInv_LAMBDA;

for j = 1:length(releaseVec)
    release_I0          = wt_I0;
    release_ratio       = releaseVec(j)/(1-releaseVec(j));

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(27) = release_ratio*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = release_ratio*wt_I0(19);
        release_I0(36) = release_ratio*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = release_ratio*wt_I0(28);
    end

    [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts); 
    % freq. of (female) drive homozygotes as fraction of wt equil pop
    outputMat(matInd,j) = sum(xout(end,28:end))/2000;
end

plot(releaseVec, outputMat,'linewidth',1.5)
ylim([0,1])
yticks([0,0.2,0.4,0.6,0.8,1.0])
xticks([0,0.2,0.4,0.6,0.8,1.0])
grid on
legend('H','WSNH','SSNH','location','best','interpreter','latex')
title('(a)')
set(gca, 'fontsize',16)



%% 4.2.2. As above but with MOR instead of BSR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% First part studies a bi-sex release for each UD case study. %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,1,2)

releaseVec      = 0.0:0.001:1.0;         
outputMat       = zeros(3,length(releaseVec));
tVec            = 0:365;        

%%%%%%%%%%%%%%%%%%%%%
%%%% HOMOLOGOUS %%%%%
%%%%%%%%%%%%%%%%%%%%%

matInd                  = 1;
params                  = struct();
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 16;
params.muX      = 0.03;
params.muZ      = log(10/9);
params.muY      = log(10/9);
params.m        = 0.1;
params.s_t      = 0.9;
params.s_a      = 0.05;
N               = 2*(10^3); 

params.lethality_type   = "EA";
params.supp             = "SS";

releaseType             = "MOR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

beta            = 0.9;
alpha           = (m/(2*2*1000*muZ))*((exp(lambda*m/(2*muZ)-m-muX) - 1)^(1/beta));
fInv            = @(x) (1/alpha)*((exp(x)-1).^(1/beta));
fInv_LAMBDA     = fInv(LAMBDA); 

params.beta     = beta;
params.alpha    = alpha;                

wt_I0                   = zeros(1,24);    
wt_I0(19)               = N;
wt_I0(7)                = (1/2)*fInv_LAMBDA;
wt_I0(1)                = (1/2)*fInv_LAMBDA; 
wt_I0(13)               = (m/(2*muZ))*fInv_LAMBDA;

for j = 1:length(releaseVec)
    release_I0          = wt_I0;
    release_ratio       = releaseVec(j)/(1-releaseVec(j));

    if releaseType == "MOR"
        release_I0(16) = release_ratio*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = release_ratio*wt_I0(13);
        release_I0(22) = release_ratio*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = release_ratio*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts); 
    % calculate drive success metric
    outputMat(matInd,j) = sum(xout(end,19:end))/2000;

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24

end

%%%%%%%%%%%%%%%%%%%%%
% WS NON-HOMOLOGOUS %
%%%%%%%%%%%%%%%%%%%%%

matInd                  = 2;
params                  = struct();
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 16;
params.muX      = 0.03;
params.muZ      = log(10/9);
params.muY      = log(10/9);
params.m        = 0.1;
params.s_t      = 0.9;
params.s_a      = 0.05;
N               = 2*(10^3); 

params.lethality_type   = "EA";
params.supp             = "WS";

releaseType             = "MOR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

beta            = 0.9;
alpha           = (m/(2*2*1000*muZ))*((exp(lambda*m/(2*muZ)-m-muX) - 1)^(1/beta));
fInv            = @(x) (1/alpha)*((exp(x)-1).^(1/beta));
fInv_LAMBDA     = fInv(LAMBDA); 

params.beta     = beta;
params.alpha    = alpha;                

wt_I0                   = zeros(1,36);    
wt_I0(28)               = N;
wt_I0(1)                = (1/2)*fInv_LAMBDA;
wt_I0(2)                = (1/2)*fInv_LAMBDA; 
wt_I0(19)               = (m/(2*muZ))*fInv_LAMBDA;

for j = 1:length(releaseVec)
    release_I0          = wt_I0;
    release_ratio       = releaseVec(j)/(1-releaseVec(j));

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(27) = release_ratio*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = release_ratio*wt_I0(19);
        release_I0(36) = release_ratio*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = release_ratio*wt_I0(28);
    end

    [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts); 
    % freq. of (female) drive homozygotes as fraction of wt equil pop
    outputMat(matInd,j) = sum(xout(end,28:end))/2000;
end

%%%%%%%%%%%%%%%%%%%%%
% SS NON-HOMOLOGOUS %
%%%%%%%%%%%%%%%%%%%%%

matInd                  = 3;
params                  = struct();
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 16;
params.muX      = 0.03;
params.muZ      = log(10/9);
params.muY      = log(10/9);
params.m        = 0.1;
params.s_t      = 0.9;
params.s_a      = 0.05;
N               = 2*(10^3); 

params.lethality_type   = "EA";
params.supp             = "SS";

releaseType             = "MOR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

beta            = 0.9;
alpha           = (m/(2*2*1000*muZ))*((exp(lambda*m/(2*muZ)-m-muX) - 1)^(1/beta));
fInv            = @(x) (1/alpha)*((exp(x)-1).^(1/beta));
fInv_LAMBDA     = fInv(LAMBDA); 

params.beta     = beta;
params.alpha    = alpha;                

wt_I0                   = zeros(1,36);    
wt_I0(28)               = N;
wt_I0(1)                = (1/2)*fInv_LAMBDA;
wt_I0(2)                = (1/2)*fInv_LAMBDA; 
wt_I0(19)               = (m/(2*muZ))*fInv_LAMBDA;

for j = 1:length(releaseVec)
    release_I0          = wt_I0;
    release_ratio       = releaseVec(j)/(1-releaseVec(j));

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(27) = release_ratio*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = release_ratio*wt_I0(19);
        release_I0(36) = release_ratio*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = release_ratio*wt_I0(28);
    end

    [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts); 
    % freq. of (female) drive homozygotes as fraction of wt equil pop
    outputMat(matInd,j) = sum(xout(end,28:end))/2000;
end

plot(releaseVec, outputMat,'linewidth',1.5)
ylim([0,1])
yticks([0,0.2,0.4,0.6,0.8,1.0])
xticks([0,0.2,0.4,0.6,0.8,1.0])
grid on
legend('H','WSNH','SSNH','location','best','interpreter','latex')
title('(b)')
set(gca, 'fontsize',16)

%% 5.1. Simulating a homing gene drive for a MOR for both BSL and FSL 
% cases. Since homing drives are not threshold drives, we only study level
% of suppression achieved and keep the release ratio fixed at 1%. 

s_aVec          = 0.0:0.01:1.0;
outputMat       = zeros(3,length(s_aVec));
tVec            = 0:(10*365);        

release_ratio   = 0.01;    

% set parameters
params          = struct();
params.lambda   = 8;
params.CONV_EFF = 1; 
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.h        = 0.5;
N               = 2*(10^4); 

params.lethality_type   = "LA";
params.lethality_case   = "FSL";

opts                    = odeset('RelTol',1e-10,'AbsTol',1e-12,'NonNegative',1:12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
beta            = 0.5;
alpha           = 0.0321302;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

% [AA_male_juveniles, Aa_male_juveniles,      2
%  aa_male_juveniles, AA_female_juveniles,    4
%  Aa_female_juveniles, aa_female_juveniles,  6
%  AA_adult_males, Aa_adult_males,            8
%  aa_adult_males, AA_adult_females,          10
%  Aa_adult_females, aa_adult_females]        12

wt_I0                   = zeros(1,12);    
wt_I0(10)               = N;
wt_I0(4)                = (1/2)*fInv_LAMBDA;
wt_I0(1)                = (1/2)*fInv_LAMBDA; 
wt_I0(7)                = (m/(2*muZ))*fInv_LAMBDA;
% release transgenic males
release_I0              = wt_I0;
release_I0(8)           = ceil(release_ratio*wt_I0(7));

for i = 1:length(s_aVec)
    disp(i)
    params.s                = s_aVec(i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 0.5) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    matInd                  = 1;
    
    params.beta             = 0.5;
    params.alpha            = 0.0321302;
    
    % see how the homing drive performs
    [~,xout] = ode15s(@(t,x) DDGD_ii_HOM(t,x,params),tVec,release_I0,opts);
    % store relative adult female pop.
    outputMat(matInd, i)    = sum(xout(end,10:12))/N; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 1.0) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    matInd                  = 2;
    
    params.beta             = 1.0;
    params.alpha            = 1.90085*10^(-4);
    
    % see how the homing drive performs
    [~,xout] = ode15s(@(t,x) DDGD_ii_HOM(t,x,params),tVec,release_I0,opts);
    % store relative adult female pop.
    outputMat(matInd, i)    = sum(xout(end,10:12))/N;  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 1.5) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    matInd                  = 3;
    
    params.beta             = 1.5;
    params.alpha            = 1.12456*10^(-6);
    
    % see how the homing drive performs
    [~,xout] = ode15s(@(t,x) DDGD_ii_HOM(t,x,params),tVec,release_I0,opts);
    % store relative adult female pop.
    outputMat(matInd, i)    = sum(xout(end,10:12))/N;  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% LOGARITHMIC %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    matInd                  = 4;
    
    params.beta             = 0.9;
    params.alpha            = 0.014543805956894;    
    
    % see how the homing drive performs
    [~,xout] = ode15s(@(t,x) DDGD_iii_HOM(t,x,params),tVec,release_I0,opts);
    % store relative adult female pop.
    outputMat(matInd, i)    = sum(xout(end,10:12))/N;      

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (params.lethality_case == "FSL")
    % FSL case
    subplot(1,2,1)
    % relative female pops at control equil (FSL)
    plot(s_aVec, outputMat(1,:),'-k','linewidth',1.5);
    ylim([0,1])
    xlim([0,1])
    title('(a)')
    hold on
    plot(s_aVec, outputMat(2,:),'--k','linewidth',1.5);
    plot(s_aVec, outputMat(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, outputMat(4,:),':k','linewidth',1.5);
    ylabel('relative female pop.','interpreter','latex');
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]);     
    xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]);     
    xlabel('fitness cost, $c_a$','interpreter','latex');
    set(gca,'fontsize',16)
elseif (params.lethality_case == "BSL")
    % BSL case
    subplot(1,2,2)
    % relative female pops at control equil (FSL)
    plot(s_aVec, outputMat(1,:),'-k','linewidth',1.5);
    ylim([0,1])
    xlim([0,1])
    title('(b)')
    hold on
    plot(s_aVec, outputMat(2,:),'--k','linewidth',1.5);
    plot(s_aVec, outputMat(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, outputMat(4,:),':k','linewidth',1.5);
    ylabel('relative female pop.','interpreter','latex');
    yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); 
    xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); 
    xlabel('fitness cost, $c_a$','interpreter','latex');
    set(gca,'fontsize',16)
end

% compare FSL to BSL results in the same plot; assumes that both
% outputMat_FSL and outputMat_BSL exist
%%% performance is identical when I last checked (09/27)

% plot(s_aVec, outputMat_BSL(1,:),'-k','linewidth',1.5);
% ylim([0,1])
% xlim([0,1])
% title('FSL (red) and BSL (black) comp')
% hold on
% plot(s_aVec, outputMat_BSL(2,:),'--k','linewidth',1.5);
% plot(s_aVec, outputMat_BSL(3,:),'-.k','linewidth',1.5);
% plot(s_aVec, outputMat_BSL(4,:),':k','linewidth',1.5);
% plot(s_aVec, outputMat_FSL(1,:),'-r','linewidth',1.5);
% plot(s_aVec, outputMat_FSL(2,:),'--r','linewidth',1.5);
% plot(s_aVec, outputMat_FSL(3,:),'-.r','linewidth',1.5);
% plot(s_aVec, outputMat_FSL(4,:),':r','linewidth',1.5);
% ylabel('relative female pop.','interpreter','latex');
% yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); 
% xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); 
% xlabel('fitness cost, $c_a$','interpreter','latex');
% set(gca,'fontsize',16)

%% 6.1. Exploring how DD strength affects invasion threshold for 2LUD. 
% Current hypothesis is that non-viable genotypes are formed on release, so
% that a BSR will remove difference in invasion thresholds. 
%
% Both FSL and BSL scenarios are studied.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% SET FSL OR BSL %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params                  = struct();
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.1; % 0.28;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NOTE MORTALITY RATES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.muY      = 0.1;
params.m        = 0.14;
params.s_t      = 1.0;
N               = 2*(10^4); 

params.lethality_type   = "LA";
params.supp             = "SS";


releaseType             = "MOR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

% study invasion over range of fitness costs
s_aVec                  = 0:0.01:0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
beta            = 0.5;
alpha           = 0.0321302;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

% [AABB_larvae_male, AABB_larvae_female    2
%  AABb_larvae_male, AABb_larvae_female    4
%  AAbb_larvae_male, AAbb_larvae_female    6
%  AaBB_larvae_male, AaBB_larvae_female    8
%  AaBb_larvae_male, AaBb_larvae_female    10
%  Aabb_larvae_male, Aabb_larvae_female    12
%  aaBB_larvae_male, aaBB_larvae_female    14
%  aaBb_larvae_male, aaBb_larvae_female    16
%  aabb_larvae_male, aabb_larvae_female    18
%  AABB_adult_males, AABb_adult_males      20
%  AAbb_adult_males, AaBB_adult_males      22
%  AaBb_adult_males, Aabb_adult_males      24
%  aaBB_adult_males, aaBb_adult_males      26
%  aabb_adult_males, AABB_adult_females    28
%  AABb_adult_females, AAbb_adult_females  30
%  AaBB_adult_females, AaBb_adult_females  32
%  Aabb_adult_females, aaBB_adult_females  34
%  aaBb_adult_females, aabb_adult_females] 36

wt_I0                   = zeros(1,36);    
wt_I0(28)               = N;
wt_I0(1)                = (1/2)*fInv_LAMBDA;
wt_I0(2)                = (1/2)*fInv_LAMBDA; 
wt_I0(19)               = (m/(2*muZ))*fInv_LAMBDA;

% release transgenic males
release_I0              = wt_I0;

% matrices for storing everything
invasionThresh          = zeros(4,length(s_aVec)); 
controlEquil            = zeros(4,length(s_aVec)); 

tmax = 5*365;
tVec = 0:1:tmax; 
% study release thresholds up to 12
rFlagMax = 12;     
%%% counter = 1;

for i = 1:length(s_aVec) 
    % loading message
    fprintf("Running sim %.0f of %.0f: s = %0.2f.\n", i, length(s_aVec), ...
        s_aVec(i));

    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;
    %%% releaseMult = 1.5; 

    % set s_a
    params.s_a = s_aVec(i);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 0.5) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    params.beta             = 0.5;
    params.alpha            = 0.0321302;

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 

    % [AABB_larvae_male, AABB_larvae_female    2
    %  AABb_larvae_male, AABb_larvae_female    4
    %  AAbb_larvae_male, AAbb_larvae_female    6
    %  AaBB_larvae_male, AaBB_larvae_female    8
    %  AaBb_larvae_male, AaBb_larvae_female    10
    %  Aabb_larvae_male, Aabb_larvae_female    12
    %  aaBB_larvae_male, aaBB_larvae_female    14
    %  aaBb_larvae_male, aaBb_larvae_female    16
    %  aabb_larvae_male, aabb_larvae_female    18
    %  AABB_adult_males, AABb_adult_males      20
    %  AAbb_adult_males, AaBB_adult_males      22
    %  AaBb_adult_males, Aabb_adult_males      24
    %  aaBB_adult_males, aaBb_adult_males      26
    %  aabb_adult_males, AABB_adult_females    28
    %  AABb_adult_females, AAbb_adult_females  30
    %  AaBB_adult_females, AaBb_adult_females  32
    %  Aabb_adult_females, aaBB_adult_females  34
    %  aaBb_adult_females, aabb_adult_females] 36  

    % check allelic frequency
    ind1 = [3,4,7,8,20,22,29,31];
        % one copy of a or b
    ind2 = [5,6,9,10,13,14,21,23,25,30,32,34];
        % two copies of a or b
    ind3 = [11,12,15,16,24,26,33,35];
        % three copies of a or b
    ind4 = [17,18,27,36];
        % four copies of a or b        

    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
        4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
        4*sum(xout(:,ind4),2))./(4*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    %%%
    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
            4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
            4*sum(xout(:,ind4),2))./(4*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(27) = midFlag*wt_I0(19);
            elseif releaseType == "BSR"
                release_I0(27) = midFlag*wt_I0(19);
                release_I0(36) = midFlag*wt_I0(28);
            elseif releaseType == "FOR"
                release_I0(36) = midFlag*wt_I0(28);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                4*sum(xout(:,ind4),2))./(4*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                    4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                    4*sum(xout(:,ind4),2))./(4*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(1,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(1,i)   = sum(xout(end,28:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(1,i)   = nan; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 1.0) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 2;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;
    %%% releaseMult = 1.5; 

    % set s_a
    params.s_a = s_aVec(i);   

    params.beta             = 1.0;
    params.alpha            = 1.90085*10^(-4);

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 

    % [AABB_larvae_male, AABB_larvae_female    2
    %  AABb_larvae_male, AABb_larvae_female    4
    %  AAbb_larvae_male, AAbb_larvae_female    6
    %  AaBB_larvae_male, AaBB_larvae_female    8
    %  AaBb_larvae_male, AaBb_larvae_female    10
    %  Aabb_larvae_male, Aabb_larvae_female    12
    %  aaBB_larvae_male, aaBB_larvae_female    14
    %  aaBb_larvae_male, aaBb_larvae_female    16
    %  aabb_larvae_male, aabb_larvae_female    18
    %  AABB_adult_males, AABb_adult_males      20
    %  AAbb_adult_males, AaBB_adult_males      22
    %  AaBb_adult_males, Aabb_adult_males      24
    %  aaBB_adult_males, aaBb_adult_males      26
    %  aabb_adult_males, AABB_adult_females    28
    %  AABb_adult_females, AAbb_adult_females  30
    %  AaBB_adult_females, AaBb_adult_females  32
    %  Aabb_adult_females, aaBB_adult_females  34
    %  aaBb_adult_females, aabb_adult_females] 36  

    % check allelic frequency
    ind1 = [3,4,7,8,20,22,29,31];
        % one copy of a or b
    ind2 = [5,6,9,10,13,14,21,23,25,30,32,34];
        % two copies of a or b
    ind3 = [11,12,15,16,24,26,33,35];
        % three copies of a or b
    ind4 = [17,18,27,36];
        % four copies of a or b        

    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
        4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
        4*sum(xout(:,ind4),2))./(4*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    %%%
    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
            4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
            4*sum(xout(:,ind4),2))./(4*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(27) = midFlag*wt_I0(19);
            elseif releaseType == "BSR"
                release_I0(27) = midFlag*wt_I0(19);
                release_I0(36) = midFlag*wt_I0(28);
            elseif releaseType == "FOR"
                release_I0(36) = midFlag*wt_I0(28);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                4*sum(xout(:,ind4),2))./(4*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                    4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                    4*sum(xout(:,ind4),2))./(4*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,28:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 1.5) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 3;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;
    %%% releaseMult = 1.5; 

    % set s_a
    params.s_a = s_aVec(i);   

    params.beta             = 1.5;
    params.alpha            = 1.12456*10^(-6);

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 

    % [AABB_larvae_male, AABB_larvae_female    2
    %  AABb_larvae_male, AABb_larvae_female    4
    %  AAbb_larvae_male, AAbb_larvae_female    6
    %  AaBB_larvae_male, AaBB_larvae_female    8
    %  AaBb_larvae_male, AaBb_larvae_female    10
    %  Aabb_larvae_male, Aabb_larvae_female    12
    %  aaBB_larvae_male, aaBB_larvae_female    14
    %  aaBb_larvae_male, aaBb_larvae_female    16
    %  aabb_larvae_male, aabb_larvae_female    18
    %  AABB_adult_males, AABb_adult_males      20
    %  AAbb_adult_males, AaBB_adult_males      22
    %  AaBb_adult_males, Aabb_adult_males      24
    %  aaBB_adult_males, aaBb_adult_males      26
    %  aabb_adult_males, AABB_adult_females    28
    %  AABb_adult_females, AAbb_adult_females  30
    %  AaBB_adult_females, AaBb_adult_females  32
    %  Aabb_adult_females, aaBB_adult_females  34
    %  aaBb_adult_females, aabb_adult_females] 36  

    % check allelic frequency
    ind1 = [3,4,7,8,20,22,29,31];
        % one copy of a or b
    ind2 = [5,6,9,10,13,14,21,23,25,30,32,34];
        % two copies of a or b
    ind3 = [11,12,15,16,24,26,33,35];
        % three copies of a or b
    ind4 = [17,18,27,36];
        % four copies of a or b        

    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
        4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
        4*sum(xout(:,ind4),2))./(4*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    %%%
    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
            4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
            4*sum(xout(:,ind4),2))./(4*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(27) = midFlag*wt_I0(19);
            elseif releaseType == "BSR"
                release_I0(27) = midFlag*wt_I0(19);
                release_I0(36) = midFlag*wt_I0(28);
            elseif releaseType == "FOR"
                release_I0(36) = midFlag*wt_I0(28);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                4*sum(xout(:,ind4),2))./(4*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                    4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                    4*sum(xout(:,ind4),2))./(4*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_2LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,28:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% LOGARITHMIC %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 4;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;
    %%% releaseMult = 1.5; 

    % set s_a
    params.s_a = s_aVec(i);   

    params.beta             = 0.9;
    params.alpha            = 0.014543805956894;

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end

    [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts); 

    % [AABB_larvae_male, AABB_larvae_female    2
    %  AABb_larvae_male, AABb_larvae_female    4
    %  AAbb_larvae_male, AAbb_larvae_female    6
    %  AaBB_larvae_male, AaBB_larvae_female    8
    %  AaBb_larvae_male, AaBb_larvae_female    10
    %  Aabb_larvae_male, Aabb_larvae_female    12
    %  aaBB_larvae_male, aaBB_larvae_female    14
    %  aaBb_larvae_male, aaBb_larvae_female    16
    %  aabb_larvae_male, aabb_larvae_female    18
    %  AABB_adult_males, AABb_adult_males      20
    %  AAbb_adult_males, AaBB_adult_males      22
    %  AaBb_adult_males, Aabb_adult_males      24
    %  aaBB_adult_males, aaBb_adult_males      26
    %  aabb_adult_males, AABB_adult_females    28
    %  AABb_adult_females, AAbb_adult_females  30
    %  AaBB_adult_females, AaBb_adult_females  32
    %  Aabb_adult_females, aaBB_adult_females  34
    %  aaBb_adult_females, aabb_adult_females] 36  

    % check allelic frequency
    ind1 = [3,4,7,8,20,22,29,31];
        % one copy of a or b
    ind2 = [5,6,9,10,13,14,21,23,25,30,32,34];
        % two copies of a or b
    ind3 = [11,12,15,16,24,26,33,35];
        % three copies of a or b
    ind4 = [17,18,27,36];
        % four copies of a or b        

    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
        4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
        4*sum(xout(:,ind4),2))./(4*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    %%%
    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
            4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
            4*sum(xout(:,ind4),2))./(4*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(27) = midFlag*wt_I0(19);
            elseif releaseType == "BSR"
                release_I0(27) = midFlag*wt_I0(19);
                release_I0(36) = midFlag*wt_I0(28);
            elseif releaseType == "FOR"
                release_I0(36) = midFlag*wt_I0(28);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                4*sum(xout(:,ind4),2))./(4*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)) + 3*sum(xout(end,ind3)) + ...
                    4*sum(xout(end,ind4)))/(4*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2) + 3*sum(xout(:,ind3),2) + ...
                    4*sum(xout(:,ind4),2))./(4*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(27) = rFlag*wt_I0(19);
    elseif releaseType == "BSR"
        release_I0(27) = rFlag*wt_I0(19);
        release_I0(36) = rFlag*wt_I0(28);
    elseif releaseType == "FOR"
        release_I0(36) = rFlag*wt_I0(28);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_iii_2LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,28:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% only plot invasion thresholds

if (params.lethality_case == "FSL")
    % FSL case    
    subplot(1,2,1)
    % invasion threshold (FSL)
    plot(s_aVec, invasionThresh(1,:),'-k','linewidth',1.5);
    xlim([0,s_aVec(min(find(isnan(controlEquil(1,:))))-1)]);
    ylim([0,12]);
    title('(a)')
    hold on
    plot(s_aVec, invasionThresh(2,:),'--k','linewidth',1.5);
    plot(s_aVec, invasionThresh(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, invasionThresh(4,:),':k','linewidth',1.5);
    ylabel('invasion threshold','interpreter','latex');
    yticks([0, 2, 4, 6, 8, 10, 12]); 
    xlabel('fitness cost, $c_a$','interpreter','latex'); 
    set(gca,'fontsize',16)
elseif (params.lethality_case == "BSL")    
    subplot(1,2,2)
    % invasion threshold (BSL)
    plot(s_aVec, invasionThresh(1,:),'-k','linewidth',1.5);
    xlim([0,max(s_aVec)]);
    ylim([0,12]);
    title('(b)')
    hold on
    plot(s_aVec, invasionThresh(2,:),'--k','linewidth',1.5);
    plot(s_aVec, invasionThresh(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, invasionThresh(4,:),':k','linewidth',1.5);
    yticks([0, 2, 4, 6, 8, 10, 12]); 
    ylabel('invasion threshold','interpreter','latex');
    xlabel('fitness cost, $c_a$','interpreter','latex'); 
    set(gca,'fontsize',16)
else
    error("Unrecognized lethality case; did you mean BSL or FSL?")
end

%% 6.2. Exploring how DD strength affects invasion threshold for 1LUD. 
% Considers single-sex releases for c_t < 1.
%
% Both FSL and BSL scenarios are studied.

% close all 

params          = struct();

% specify lethality 
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.s_t      = 1.0;
N               = 2*(10^4); 

params.lethality_type   = "LA";


releaseType             = "BSR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

% study invasion over range of fitness costs
s_aVec                  = 0:0.01:0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
beta            = 0.5;
alpha           = 0.0321302;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

% [AA_larvae_males, Aa_larvae_males,        2
%  Ab_larvae_males, ab_larvae_males,        4
%  aa_larvae_males, bb_larvae_males,        6
%  AA_larvae_females, Aa_larvae_females,    8
%  Ab_larvae_females, ab_larvae_females,    10
%  aa_larvae_females, bb_larvae_females,    12
%  AA_adults_males, Aa_adults_males,        14
%  Ab_adults_males, ab_adults_males,        16
%  aa_adults_males, bb_adults_males,        18
%  AA_adults_females, Aa_adults_females,    20
%  Ab_adults_females, ab_adults_females,    22
%  aa_adults_females, bb_adults_females]    24

wt_I0                   = zeros(1,24);    
wt_I0(19)               = N;
wt_I0(7)                = (1/2)*fInv_LAMBDA;
wt_I0(1)                = (1/2)*fInv_LAMBDA; 
wt_I0(13)               = (m/(2*muZ))*fInv_LAMBDA;
% release transgenic males
release_I0              = wt_I0;


% matrices for storing everything
invasionThresh          = zeros(4,length(s_aVec)); 
controlEquil            = zeros(4,length(s_aVec)); 

tmax = 5*365;
tVec = 0:1:tmax; 
% study release thresholds up to 12
rFlagMax = 12;     
%%% counter = 1;

for i = 1:length(s_aVec) 
    % loading message
    fprintf("Running sim %.0f of %.0f: s = %0.2f.\n", i, length(s_aVec), ...
        s_aVec(i));

    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;
    %%% releaseMult = 1.5; 

    % set s_a
    params.s_a = s_aVec(i);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 0.5) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    params.beta             = 0.5;
    params.alpha            = 0.0321302;

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24    

    % check allelic frequency
    ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
    ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(16) = midFlag*wt_I0(13);
            elseif releaseType == "BSR"
                release_I0(16) = midFlag*wt_I0(13);
                release_I0(22) = midFlag*wt_I0(19);
            elseif releaseType == "FOR"
                release_I0(22) = midFlag*wt_I0(19);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(1,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(1,i)   = sum(xout(end,19:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(1,i)   = nan; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 1.0) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 2;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;

    % set s_a
    params.s_a = s_aVec(i);    

    params.beta             = 1.0;
    params.alpha            = 1.90085*10^(-4);

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24    

    % check allelic frequency
    ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
    ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(16) = midFlag*wt_I0(13);
            elseif releaseType == "BSR"
                release_I0(16) = midFlag*wt_I0(13);
                release_I0(22) = midFlag*wt_I0(19);
            elseif releaseType == "FOR"
                release_I0(22) = midFlag*wt_I0(19);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,19:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% GEN LOGISTIC (BETA = 1.5) %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 3;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;

    % set s_a
    params.s_a = s_aVec(i);    

    params.beta             = 1.5;
    params.alpha            = 1.12456*10^(-6);

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24    

    % check allelic frequency
    ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
    ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(16) = midFlag*wt_I0(13);
            elseif releaseType == "BSR"
                release_I0(16) = midFlag*wt_I0(13);
                release_I0(22) = midFlag*wt_I0(19);
            elseif releaseType == "FOR"
                release_I0(22) = midFlag*wt_I0(19);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_ii_1LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,19:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% GEN LOGARITHMIC %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matInd = 4;
    % left flag, right flag
    lFlag = 0.01;
    rFlag = rFlagMax;

    % set s_a
    params.s_a = s_aVec(i);    

    params.beta             = 0.9;
    params.alpha            = 0.014543805956894;

    % adjusted release based on releaseType
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end

    [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts); 

    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24    

    % check allelic frequency
    ind1 = [2,3,8,9,14,15,20,21];                   % one copy of a or b
    ind2 = [4,5,6,10,11,12,16,17,18,22,23,24];      % two copies of a or b
    % is it above frequency of 5%?
    fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
    % is it unchanging?
    freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
    freq_steady = var(freq_steady((end-99):end)) < 0.01;
    thresh_invade = 0.05; 

    increaseCount = 2;    
    while (~freq_steady)
        disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
        tVecTmp = 0:1:(increaseCount*tmax); 
        
        [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
        % is it above frequency of 5%?
        fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
        % is it unchanging?
        freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
        freq_steady = var(freq_steady((end-99):end)) < 0.01;
        
        increaseCount = increaseCount + 1;
    end

    if ((fixProp_tmp > thresh_invade) && freq_steady) 
        % convergence is successful!
        
        while ~(abs(lFlag - rFlag) < 10^(-2))
            % disp("Error: " + abs(lFlag-rFlag)); 
            
            midFlag = (lFlag + rFlag)/2;
            % adjusted release based on releaseType
            if releaseType == "MOR"
                release_I0(16) = midFlag*wt_I0(13);
            elseif releaseType == "BSR"
                release_I0(16) = midFlag*wt_I0(13);
                release_I0(22) = midFlag*wt_I0(19);
            elseif releaseType == "FOR"
                release_I0(22) = midFlag*wt_I0(19);
            end
            
            [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts); 
            
            % is it above frequency of 5%?
            fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
            % is it unchanging?
            freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
            freq_steady = var(freq_steady((end-99):end)) < 0.01;
            
            increaseCount = 2;
            while (~freq_steady)
                disp('Numerical equil. convergence unsuccessful --- increasing time.\n');
                tVecTmp = 0:1:(increaseCount*tmax); 
                
                [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVecTmp,release_I0,opts); 
                % is it above frequency of 5%?
                fixProp_tmp = (sum(xout(end,ind1)) + 2*sum(xout(end,ind2)))/(2*sum(xout(end,:)));
                % is it unchanging?
                freq_steady = (sum(xout(:,ind1),2) + 2*sum(xout(:,ind2),2))./(2*sum(xout,2));
                freq_steady = var(freq_steady((end-99):end)) < 0.01;
                
                increaseCount = increaseCount + 1;
            end
        
            % boolean variable -> 1 if "fixed," 0 otherwise
            if ((fixProp_tmp > thresh_invade) && freq_steady) 
                tmpBool = 1;
            else
                tmpBool = 0;
            end
            % no fixation, right flag is right flag
            % fixation, right flag is mid flag
            rFlag = (tmpBool)*midFlag + (1-tmpBool)*rFlag;
            % no fixation, left flag is mid flag
            % fixation, left flag is left flag
            lFlag = (tmpBool)*lFlag + (1-tmpBool)*midFlag;
        end % end of while flag loop
    else 
        % convergence unsuccessful, return NaN
        rFlag = NaN; 
    end 
    % end of if thresh_invade statement    

    % store data
    invasionThresh(matInd,i)     = rFlag;
    % run a simulation using the invasion threshold
    if releaseType == "MOR"
        release_I0(16) = rFlag*wt_I0(13);
    elseif releaseType == "BSR"
        release_I0(16) = rFlag*wt_I0(13);
        release_I0(22) = rFlag*wt_I0(19);
    elseif releaseType == "FOR"
        release_I0(22) = rFlag*wt_I0(19);
    end
        
    [~,xout] = ode45(@(t,x) DDGD_iii_1LUD(t,x,params),tVec,release_I0,opts);

    % after convergence is achieved, what is control equil. of adult female
    % pop?
    if ~isnan(rFlag)
        % invasion successful (release < 12)
        controlEquil(matInd,i)   = sum(xout(end,19:end))/(2*(10^4)); 
    else
        % invasion unsuccessful (release > 12)
        controlEquil(matInd,i)   = nan; 
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (params.lethality_case == "FSL")    
    subplot(1,2,1)
    % invasion threshold (FSL)
    plot(s_aVec, invasionThresh(1,:),'-k','linewidth',1.5);
    ylim([0,12]);
    xticks([0, 0.05, 0.1, 0.15, 0.2]);
    title('(a)')
    hold on
    plot(s_aVec, invasionThresh(2,:),'--k','linewidth',1.5);
    plot(s_aVec, invasionThresh(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, invasionThresh(4,:),':k','linewidth',1.5);
    ylabel('invasion threshold','interpreter','latex');
    yticks([0, 2, 4, 6, 8, 10, 12]);
    xlabel('fitness cost, $c_a$','interpreter','latex');
    set(gca,'fontsize',16)    
elseif (params.lethality_case == "BSL")    
    subplot(1,2,2)
    plot(s_aVec, invasionThresh(1,:),'-k','linewidth',1.5);
    ylim([0,12]);
    xticks([0, 0.05, 0.1, 0.15, 0.2]);
    % xlim([0,0.11])
    % xlim([0,s_aVec(min(find(isnan(controlEquil(1,:))))-1)]);
    title('(b)')
    hold on
    plot(s_aVec, invasionThresh(2,:),'--k','linewidth',1.5);
    plot(s_aVec, invasionThresh(3,:),'-.k','linewidth',1.5);
    plot(s_aVec, invasionThresh(4,:),':k','linewidth',1.5);
    ylabel('invasion threshold','interpreter','latex');
    yticks([0, 2, 4, 6, 8, 10, 12]);
    xlabel('fitness cost, $c_a$','interpreter','latex');
    set(gca,'fontsize',16)    
end

%% 7.1. For the 1LUD case, invasion threshold is very sensitive to the 
% lifespan of males. This script calculates allele freq. after 1 year 
% for different assumptions of male lifespan for all strengths of DD.
%
% Only a male-only release is considered, and c_t = 0.9 as in Khamis et 
% al.

params          = struct();

% specify lethality 
params.lethality_case   = "BSL";

% set parameters
params.lambda   = 8;
params.muX      = 0.029;
params.muZ      = 0.28;
params.muY      = 0.1;
params.m        = 0.14;
params.s_t      = 0.9;
N               = 2*(10^4); 

params.lethality_type   = "LA";


releaseType             = "MOR"; % adjusts release type for ALL sims in 
                                 % this section
opts                    = odeset('RelTol',1e-8,'AbsTol',1e-9);

% study invasion over range of fitness costs
s_aVec                  = 0:0.01:0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% calculate pre-release equilibrium %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda          = params.lambda;
muX             = params.muX;
muZ             = params.muZ;
muY             = params.muY;
m               = params.m;
LAMBDA          = (m*lambda/(2*muY))-muX-m;

% pre-release equil. equal between all functional forms, choose one to
% calculate this equil.
beta            = 0.5;
alpha           = 0.0321302;
fInv            = @(x) (x/alpha).^(1/beta);
fInv_LAMBDA     = fInv(LAMBDA); 

% [AA_larvae_males, Aa_larvae_males,        2
%  Ab_larvae_males, ab_larvae_males,        4
%  aa_larvae_males, bb_larvae_males,        6
%  AA_larvae_females, Aa_larvae_females,    8
%  Ab_larvae_females, ab_larvae_females,    10
%  aa_larvae_females, bb_larvae_females,    12
%  AA_adults_males, Aa_adults_males,        14
%  Ab_adults_males, ab_adults_males,        16
%  aa_adults_males, bb_adults_males,        18
%  AA_adults_females, Aa_adults_females,    20
%  Ab_adults_females, ab_adults_females,    22
%  aa_adults_females, bb_adults_females]    24

wt_I0                   = zeros(1,24);    
wt_I0(19)               = N;
wt_I0(7)                = (1/2)*fInv_LAMBDA;
wt_I0(1)                = (1/2)*fInv_LAMBDA; 
wt_I0(13)               = (m/(2*muZ))*fInv_LAMBDA;
% release transgenic males
release_I0              = wt_I0;


% matrices for storing everything
invasionThresh          = zeros(4,length(s_aVec)); 
controlEquil            = zeros(4,length(s_aVec)); 

tmax = 5*365;
tVec = 0:1:tmax; 
% study release thresholds up to 12
rFlagMax = 12;     
%%% counter = 1;


%% 7.2. Compare invasion thresholds for a bi-sex release for 1LUD and 2LUD,
% a fairer comparison than presently made. 



%% individual functions for each DD case
function [t,xout] = run_DDGD_ii(params,  tVec, beta, r_ratio, rel_case, lethal_case, supp_case)
    % lambda = params(1);
    % muX = params(2);
    muZ = params(3);
    muY = params(4);
    m = params(5);
    a = params(6);
    b = params(7);
    c = params(8);
    gamma = params(9);
    M = params(10);
    % s_a = params(11);
    % s_t = params(12);
    N = params(13); 
    
    switch beta
        case 0.5
            alpha = 0.0321302;
        case 1
            alpha = 1.90085*10^(-4);
        case 1.5
            alpha = 1.12456*10^(-6);
        otherwise
            error('Error! Beta must be 0.5, 1, or 1.5.')
    end
    
    if (rel_case == "BSR") % bi-sex release
        IC = [14285.71429, (muY/m)*N, ... % 2
            0, 0, ...             % 4  
            0, 0, ...             % 6  
            0, 0, ...             % 8
            0, 0, ...             % 10
            0, 0, ...             % 12
            0, 0, ...             % 14
            0, 0, ...             % 16
            0, 0, ...             % 18
            7142.86, 0, ...       % 20
            0, 0,  ...            % 22  
            0, 0, ...             % 24  
            0, 0, ...             % 26
            (r_ratio*N)/2, N, ... % 28
            0, 0, ...             % 30
            0, 0, ...             % 32
            0, 0, ...             % 34
            0, (r_ratio*N)/2, ... % 36
            18535.13565, 8879.397446]; % 38
    else % male only release
        IC = [14285.71429, (muY/m)*N, ... % 2
            0, 0, ...             % 4  
            0, 0, ...             % 6  
            0, 0, ...             % 8
            0, 0, ...             % 10
            0, 0, ...             % 12
            0, 0, ...             % 14
            0, 0, ...             % 16
            0, 0, ...             % 18
            7142.86, 0, ...       % 20
            0, 0,  ...            % 22  
            0, 0, ...             % 24  
            0, 0, ...             % 26
            r_ratio*N, N, ...     % 28
            0, 0, ...             % 30
            0, 0, ...             % 32
            0, 0, ...             % 34
            0, 0, ...             % 36
            18535.13565, 8879.397446]; % 38 
    end

    % ((a^2)*b*c*N-M*gamma*muY)/((a^2)*b*c+a*b*muY),((a^2)*b*c*M*N-(M^2)*gamma*muY)/((a^2)*b*c*N+a*c*M*gamma)]; 

    [t,xout] = ode45(@(t,x) DDGD_ii(t,x,params,alpha,beta,lethal_case,supp_case),tVec,IC); 

end

function [t,xout] = run_DDGD_iii(params, tVec, r_ratio, rel_case, lethal_case, supp_case)
    % lambda = params(1);
    % muX = params(2);
    muZ = params(3);
    muY = params(4);
    m = params(5);
    a = params(6);
    b = params(7);
    c = params(8);
    gamma = params(9);
    M = params(10);
    % s_a = params(11);
    % s_t = params(12);
    N = params(13); 

    alpha = 0.014543805956894;  
    beta = 0.9;

    if (rel_case == "BSR") % bi-sex release
        IC = [14285.71429, (muY/m)*N, ... % 2
            0, 0, ...             % 4  
            0, 0, ...             % 6  
            0, 0, ...             % 8
            0, 0, ...             % 10
            0, 0, ...             % 12
            0, 0, ...             % 14
            0, 0, ...             % 16
            0, 0, ...             % 18
            7142.86, 0, ...       % 20
            0, 0,  ...            % 22  
            0, 0, ...             % 24  
            0, 0, ...             % 26
            (r_ratio*N)/2, N, ... % 28
            0, 0, ...             % 30
            0, 0, ...             % 32
            0, 0, ...             % 34
            0, (r_ratio*N)/2, ... % 36
            18535.13565, 8879.397446]; % 38
    else % male only release
        IC = [14285.71429, (muY/m)*N, ... % 2
            0, 0, ...             % 4  
            0, 0, ...             % 6  
            0, 0, ...             % 8
            0, 0, ...             % 10
            0, 0, ...             % 12
            0, 0, ...             % 14
            0, 0, ...             % 16
            0, 0, ...             % 18
            7142.86, 0, ...       % 20
            0, 0,  ...            % 22  
            0, 0, ...             % 24  
            0, 0, ...             % 26
            r_ratio*N, N, ...     % 28
            0, 0, ...             % 30
            0, 0, ...             % 32
            0, 0, ...             % 34
            0, 0, ...             % 36
            18535.13565, 8879.397446]; % 38       
    end

    % ((a^2)*b*c*N-M*gamma*muY)/((a^2)*b*c+a*b*muY),((a^2)*b*c*M*N-(M^2)*gamma*muY)/((a^2)*b*c*N+a*c*M*gamma)]; 
    opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
    [t,xout] = ode45(@(t,x) DDGD_iii(t,x,params,alpha,beta,lethal_case,supp_case),tVec,IC,opts); 

end

% single locus UD code
function dx = DDGD_ii_1LUD(t,x,params)
    % Case study for one-locus engineered underdominance (male-sex
    % lethality is ignored). Lower-case letts are transgenic.
    % ---------------------------------------------------------------------
    
    % function simulating the generalized logistic case---note that the
    % function is programmed for beta+1 and not beta

    % point check---pops are positive, up to numerical error
    %%% x(x<0) = 0;
    % % % if any(x < 0)
    % % %     if (abs(min(x)) > 10^(-6))
    % % %         disp(x(x<0))
    % % %         error("Error: Populations are negative!")
    % % %     else
    % % %         % numerical error
    % % %         x(x < 0) = 0; 
    % % %     end
    % % % end


    lambda      = params.lambda;
    muX         = params.muX;
    muZ         = params.muZ;
    muY         = params.muY;
    m           = params.m;
    s_a         = params.s_a;
    s_t         = params.s_t;

    beta        = params.beta;
    alpha       = params.alpha;
    
    lethality_case = params.lethality_case;
    lethality_type = params.lethality_type;
    
    % compute relative fitnesses for ambient carriage and transgene
    c_a = 1 - s_a;
    c_t = 1 - s_t;
    
    % memory allocation
    dx=zeros(1,24);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change beta so that things make sense
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta = beta - 1; 

    x(x < 1e-10) = 0; 
    
    Nl = sum(x(1:12));      % mosquito larvae pop
    N  = sum(x(13:end));    % adult mosquito pop
    Nf = sum(x(19:end));    % adult female mosquito pop
    Nm = sum(x(13:18));
    
    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24

    % adult counts
    % NAA = x(13) + x(19);
    % NAa = x(14) + x(20);
    % NAb = x(15) + x(21);
    % Nab = x(16) + x(22);
    % Naa = x(17) + x(23);
    % Nbb = x(18) + x(24);
    % % gamete probabilities
    % fA = (2*NAA + NAa + NAb)/(2*N);
    % fa = (2*Naa + NAa + Nab)/(2*N);
    % fb = (2*Nbb + NAb + Nab)/(2*N); 

    % calculate probabilities of offspring having each genotype
    pVec = zeros(1,6);
    % pVec(1) = fA^2;
    % pVec(2) = 2*fA*fa;
    % pVec(3) = 2*fA*fb;
    % pVec(4) = 2*fa*fb;
    % pVec(5) = fa^2;
    % pVec(6) = fb^2; 

    % if pop is extinct, set pVec to zero
    if ~((Nm < 1e-10 || Nf < 1e-10))
        pVec(1) = (1/(Nm*Nf))*(x(19)*x(13) + (1/2)*(x(19)*x(14) + x(13)*x(20)) + (1/2)*(x(19)*x(15) + x(13)*x(21)) + (1/4)*x(20)*x(14) + (1/4)*x(21)*x(15) + (1/4)*(x(20)*x(15) + x(14)*x(21)));
        pVec(2) = (1/(Nm*Nf))*((1/2)*(x(19)*x(14) + x(13)*x(20)) + (1/2)*(x(13)*x(22) + x(19)*x(16)) + (x(13)*x(23) + x(19)*x(17)) + (1/2)*x(20)*x(14) + (1/4)*(x(20)*x(15) + x(14)*x(21)) + (1/4)*(x(14)*x(22) + x(20)*x(16)) + (1/2)*(x(14)*x(23) + x(20)*x(17)) + (1/2)*(x(21)*x(17) + x(15)*x(23)) + (1/4)*(x(15)*x(22) + x(21)*x(16)));
        pVec(3) = (1/(Nm*Nf))*((1/2)*(x(19)*x(15) + x(13)*x(21)) + (1/2)*(x(13)*x(22) + x(19)*x(16)) + (x(19)*x(18) + x(13)*x(24)) + (1/2)*(x(20)*x(18) + x(14)*x(24)) + (1/4)*(x(20)*x(15) + x(14)*x(21)) + (1/4)*(x(14)*x(22) + x(20)*x(16)) + (1/2)*x(21)*x(15) + (1/2)*(x(21)*x(18) + x(15)*x(24)) + (1/4)*(x(15)*x(22) + x(21)*x(16)));
        pVec(4) = (1/(Nm*Nf))*((1/4)*(x(20)*x(15) + x(14)*x(21)) + (1/4)*(x(14)*x(22) + x(20)*x(16)) + (1/2)*(x(20)*x(18) + x(14)*x(24)) + (1/4)*(x(15)*x(22) + x(21)*x(16)) + (1/2)*(x(21)*x(17) + x(15)*x(23)) + (x(23)*x(18) + x(17)*x(24)) + (1/2)*(x(17)*x(22) + x(23)*x(16)) + (1/2)*(x(18)*x(22) + x(24)*x(16)) + (1/2)*x(22)*x(16));
        pVec(5) = (1/(Nm*Nf))*((1/4)*x(20)*x(14) + (1/4)*(x(14)*x(22) + x(20)*x(16)) + (1/2)*(x(20)*x(17) + x(14)*x(23)) + (1/4)*x(22)*x(16) + (1/2)*(x(17)*x(22) + x(23)*x(16)) + x(23)*x(17));
        pVec(6) = (1/(Nm*Nf))*((1/4)*x(21)*x(15) + (1/4)*(x(15)*x(22) + x(21)*x(16)) + (1/2)*(x(21)*x(18) + x(15)*x(24)) + (1/4)*x(22)*x(16) + (1/2)*(x(18)*x(22) + x(24)*x(16)) + x(24)*x(18));
    end

    % correct for numerical error
    pVec(pVec < 0) = 0;

    % sum(pVec)
    % disp(pVec)
    % point check---probabilities sum to 1
    if (sum(abs(pVec)) - 1) > 10^(-9)
        error("Error: Offspring genotype probabilities do not sum to 1!")
    end
    % point check---probabilities are positive, up to numerical error
    if any(pVec < 0)
        if (abs(min(pVec)) > 10^(-12))
            error("Error: Offspring genotype probabilities are negative!")
        else
            % numerical error
            pVec(pVec < 0) = 0; 
        end
    end
    
    %%%%%%%%%%%%% coupled genetic/population dynamics model %%%%%%%%%%%%% 
    if (lethality_type == "EA")
        if (lethality_case == "BSL") % bi-sex lethality
            % larvae for each genotype
            % males
            dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);                % AA
            dx(2) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(2);        % Aa
            dx(3) = (lambda/2)*(Nf)*pVec(3)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(3);        % Ab
            dx(4) = (lambda/2)*(Nf)*pVec(4)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(4);        % ab
            dx(5) = (lambda/2)*(Nf)*pVec(5)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(5);    % aa
            dx(6) = (lambda/2)*(Nf)*pVec(6)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(6);    % bb
            % females
            dx(7) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(7);                % AA
            dx(8) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(8);        % Aa
            dx(9) = (lambda/2)*(Nf)*pVec(3)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(9);        % Ab
            dx(10) = (lambda/2)*(Nf)*pVec(4)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(10);      % ab
            dx(11) = (lambda/2)*(Nf)*pVec(5)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(11);  % aa
            dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(12);  % bb
        elseif (lethality_case == "FSL") % female-sex lethality
            % larvae for each genotype
            % males
            dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);                % AA
            dx(2) = (lambda/2)*(Nf)*pVec(2)*c_a - (muX+m+alpha*Nl^(beta+1))*x(2);            % Aa
            dx(3) = (lambda/2)*(Nf)*pVec(3)*c_a - (muX+m+alpha*Nl^(beta+1))*x(3);            % Ab
            dx(4) = (lambda/2)*(Nf)*pVec(4)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(4);        % ab
            dx(5) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(5);        % aa
            dx(6) = (lambda/2)*(Nf)*pVec(6)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(6);        % bb
            % females
            dx(7) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(7);                % AA
            dx(8) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(8);        % Aa
            dx(9) = (lambda/2)*(Nf)*pVec(3)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(9);        % Ab
            dx(10) = (lambda/2)*(Nf)*pVec(4)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(10);      % ab
            dx(11) = (lambda/2)*(Nf)*pVec(5)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(11);  % aa
            dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(12);  % bb            
        else
            error("Lethality not recognized. Note: Male-sex lethality not included.")
        end

        % adults for each genotype
        % males
        dx(13) = m*x(1) - muZ*x(13);    % AA male
        dx(14) = m*x(2) - muZ*x(14);    % Aa male
        dx(15) = m*x(3) - muZ*x(15);    % Ab male
        dx(16) = m*x(4) - muZ*x(16);    % ab male
        dx(17) = m*x(5) - muZ*x(17);    % aa male 
        dx(18) = m*x(6) - muZ*x(18);    % bb male
        % females
        dx(19) = m*x(7) - muY*x(19);     % AA female
        dx(20) = m*x(8) - muY*x(20);     % Aa female
        dx(21) = m*x(9) - muY*x(21);     % Ab female
        dx(22) = m*x(10) - muY*x(22);    % ab female
        dx(23) = m*x(11) - muY*x(23);    % aa female 
        dx(24) = m*x(12) - muY*x(24);    % bb female
    else % fitness costs are late-acting
        if (lethality_case == "BSL") % bi-sex lethality
            % adults for each genotype
            % males
            dx(13) = m*x(1) - muZ*x(13);                % AA male
            dx(14) = m*x(2)*c_a*c_t - muZ*x(14);        % Aa male
            dx(15) = m*x(3)*c_a*c_t - muZ*x(15);        % Ab male
            dx(16) = m*x(4)*(c_a)^2 - muZ*x(16);        % ab male
            dx(17) = m*x(5)*(c_a*c_t)^2 - muZ*x(17);    % aa male 
            dx(18) = m*x(6)*(c_a*c_t)^2 - muZ*x(18);    % bb male
            % females
            dx(19) = m*x(7) - muY*x(19);                % AA female
            dx(20) = m*x(8)*c_a*c_t - muY*x(20);        % Aa female
            dx(21) = m*x(9)*c_a*c_t - muY*x(21);        % Ab female
            dx(22) = m*x(10)*(c_a)^2 - muY*x(22);       % ab female
            dx(23) = m*x(11)*(c_a*c_t)^2 - muY*x(23);   % aa female 
            dx(24) = m*x(12)*(c_a*c_t)^2 - muY*x(24);   % bb female
        elseif (lethality_case == "FSL") % female-sex lethality
            % adults for each genotype
            % males
            dx(13) = m*x(1) - muZ*x(13);                % AA male
            dx(14) = m*x(2)*c_a - muZ*x(14);            % Aa male
            dx(15) = m*x(3)*c_a - muZ*x(15);            % Ab male
            dx(16) = m*x(4)*(c_a)^2 - muZ*x(16);        % ab male
            dx(17) = m*x(5)*(c_a)^2 - muZ*x(17);        % aa male 
            dx(18) = m*x(6)*(c_a)^2 - muZ*x(18);        % bb male
            % females
            dx(19) = m*x(7) - muY*x(19);                % AA female
            dx(20) = m*x(8)*c_a*c_t - muY*x(20);        % Aa female
            dx(21) = m*x(9)*c_a*c_t - muY*x(21);        % Ab female
            dx(22) = m*x(10)*(c_a)^2 - muY*x(22);       % ab female
            dx(23) = m*x(11)*(c_a*c_t)^2 - muY*x(23);   % aa female 
            dx(24) = m*x(12)*(c_a*c_t)^2 - muY*x(24);   % bb female  
        else
            error("Lethality not recognized. Note: Male-sex lethality not included.")
        end

        % larvae for each genotype
        % males
        dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);       % AA
        dx(2) = (lambda/2)*(Nf)*pVec(2) - (muX+m+alpha*Nl^(beta+1))*x(2);       % Aa
        dx(3) = (lambda/2)*(Nf)*pVec(3) - (muX+m+alpha*Nl^(beta+1))*x(3);       % Ab
        dx(4) = (lambda/2)*(Nf)*pVec(4) - (muX+m+alpha*Nl^(beta+1))*x(4);       % ab
        dx(5) = (lambda/2)*(Nf)*pVec(5) - (muX+m+alpha*Nl^(beta+1))*x(5);       % aa
        dx(6) = (lambda/2)*(Nf)*pVec(6) - (muX+m+alpha*Nl^(beta+1))*x(6);       % bb
        % females
        dx(7) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(7);       % AA
        dx(8) = (lambda/2)*(Nf)*pVec(2) - (muX+m+alpha*Nl^(beta+1))*x(8);       % Aa
        dx(9) = (lambda/2)*(Nf)*pVec(3) - (muX+m+alpha*Nl^(beta+1))*x(9);       % Ab
        dx(10) = (lambda/2)*(Nf)*pVec(4) - (muX+m+alpha*Nl^(beta+1))*x(10);     % ab
        dx(11) = (lambda/2)*(Nf)*pVec(5) - (muX+m+alpha*Nl^(beta+1))*x(11);     % aa
        dx(12) = (lambda/2)*(Nf)*pVec(6) - (muX+m+alpha*Nl^(beta+1))*x(12);     % bb
    end % end of if lethality_type statement
    dx = dx';
end

function dx = DDGD_iii_1LUD(t,x,params)
    % Case study for one-locus engineered underdominance (male-sex
    % lethality is ignored). Lower-case letts are transgenic.
    % ---------------------------------------------------------------------
    
    % function simulating the generalized logistic case---note that the
    % function is programmed for beta+1 and not beta

    % point check---pops are positive, up to numerical error
    %%% x(x<0) = 0;
    % % % if any(x < 0)
    % % %     if (abs(min(x)) > 10^(-6))
    % % %         disp(x(x<0))
    % % %         error("Error: Populations are negative!")
    % % %     else
    % % %         % numerical error
    % % %         x(x < 0) = 0; 
    % % %     end
    % % % end


    lambda      = params.lambda;
    muX         = params.muX;
    muZ         = params.muZ;
    muY         = params.muY;
    m           = params.m;
    s_a         = params.s_a;
    s_t         = params.s_t;

    beta        = params.beta;
    alpha       = params.alpha;
    
    lethality_case = params.lethality_case;
    lethality_type = params.lethality_type;
    
    % compute relative fitnesses for ambient carriage and transgene
    c_a = 1 - s_a;
    c_t = 1 - s_t;
    
    % memory allocation
    dx=zeros(1,24);
    
    Nl = sum(x(1:12));      % mosquito larvae pop
    N  = sum(x(13:end));    % adult mosquito pop
    Nf = sum(x(19:end));    % adult female mosquito pop
    Nm = sum(x(13:18));
    
    % [AA_larvae_males, Aa_larvae_males,        2
    %  Ab_larvae_males, ab_larvae_males,        4
    %  aa_larvae_males, bb_larvae_males,        6
    %  AA_larvae_females, Aa_larvae_females,    8
    %  Ab_larvae_females, ab_larvae_females,    10
    %  aa_larvae_females, bb_larvae_females,    12
    %  AA_adults_males, Aa_adults_males,        14
    %  Ab_adults_males, ab_adults_males,        16
    %  aa_adults_males, bb_adults_males,        18
    %  AA_adults_females, Aa_adults_females,    20
    %  Ab_adults_females, ab_adults_females,    22
    %  aa_adults_females, bb_adults_females]    24

    % adult counts
    % NAA = x(13) + x(19);
    % NAa = x(14) + x(20);
    % NAb = x(15) + x(21);
    % Nab = x(16) + x(22);
    % Naa = x(17) + x(23);
    % Nbb = x(18) + x(24);
    % % gamete probabilities
    % fA = (2*NAA + NAa + NAb)/(2*N);
    % fa = (2*Naa + NAa + Nab)/(2*N);
    % fb = (2*Nbb + NAb + Nab)/(2*N); 

    % calculate probabilities of offspring having each genotype
    pVec = zeros(1,6);
    % pVec(1) = fA^2;
    % pVec(2) = 2*fA*fa;
    % pVec(3) = 2*fA*fb;
    % pVec(4) = 2*fa*fb;
    % pVec(5) = fa^2;
    % pVec(6) = fb^2; 

    pVec(1) = (1/(Nm*Nf))*(x(19)*x(13) + (1/2)*(x(19)*x(14) + x(13)*x(20)) + (1/2)*(x(19)*x(15) + x(13)*x(21)) + (1/4)*x(20)*x(14) + (1/4)*x(21)*x(15) + (1/4)*(x(20)*x(15) + x(14)*x(21)));
    pVec(2) = (1/(Nm*Nf))*((1/2)*(x(19)*x(14) + x(13)*x(20)) + (1/2)*(x(13)*x(22) + x(19)*x(16)) + (x(13)*x(23) + x(19)*x(17)) + (1/2)*x(20)*x(14) + (1/4)*(x(20)*x(15) + x(14)*x(21)) + (1/4)*(x(14)*x(22) + x(20)*x(16)) + (1/2)*(x(14)*x(23) + x(20)*x(17)) + (1/2)*(x(21)*x(17) + x(15)*x(23)) + (1/4)*(x(15)*x(22) + x(21)*x(16)));
    pVec(3) = (1/(Nm*Nf))*((1/2)*(x(19)*x(15) + x(13)*x(21)) + (1/2)*(x(13)*x(22) + x(19)*x(16)) + (x(19)*x(18) + x(13)*x(24)) + (1/2)*(x(20)*x(18) + x(14)*x(24)) + (1/4)*(x(20)*x(15) + x(14)*x(21)) + (1/4)*(x(14)*x(22) + x(20)*x(16)) + (1/2)*x(21)*x(15) + (1/2)*(x(21)*x(18) + x(15)*x(24)) + (1/4)*(x(15)*x(22) + x(21)*x(16)));
    pVec(4) = (1/(Nm*Nf))*((1/4)*(x(20)*x(15) + x(14)*x(21)) + (1/4)*(x(14)*x(22) + x(20)*x(16)) + (1/2)*(x(20)*x(18) + x(14)*x(24)) + (1/4)*(x(15)*x(22) + x(21)*x(16)) + (1/2)*(x(21)*x(17) + x(15)*x(23)) + (x(23)*x(18) + x(17)*x(24)) + (1/2)*(x(17)*x(22) + x(23)*x(16)) + (1/2)*(x(18)*x(22) + x(24)*x(16)) + (1/2)*x(22)*x(16));
    pVec(5) = (1/(Nm*Nf))*((1/4)*x(20)*x(14) + (1/4)*(x(14)*x(22) + x(20)*x(16)) + (1/2)*(x(20)*x(17) + x(14)*x(23)) + (1/4)*x(22)*x(16) + (1/2)*(x(17)*x(22) + x(23)*x(16)) + x(23)*x(17));
    pVec(6) = (1/(Nm*Nf))*((1/4)*x(21)*x(15) + (1/4)*(x(15)*x(22) + x(21)*x(16)) + (1/2)*(x(21)*x(18) + x(15)*x(24)) + (1/4)*x(22)*x(16) + (1/2)*(x(18)*x(22) + x(24)*x(16)) + x(24)*x(18));
    
    % point check---probabilities sum to 1
    if (sum(abs(pVec)) - 1) > 10^(-9)
        error("Error: Offspring genotype probabilities do not sum to 1!")
    end
    % point check---probabilities are positive, up to numerical error
    if any(pVec < 0)
        if (abs(min(pVec)) > 10^(-12))
            error("Error: Offspring genotype probabilities are negative!")
        else
            % numerical error
            pVec(pVec < 0) = 0; 
        end
    end
    
    %%%%%%%%%%%%% coupled genetic/population dynamics model %%%%%%%%%%%%% 
    if (lethality_type == "EA")
        if (lethality_case == "BSL") % bi-sex lethality
            % larvae for each genotype
            % males
            dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);                % AA
            dx(2) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(2);        % Aa
            dx(3) = (lambda/2)*(Nf)*pVec(3)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(3);        % Ab
            dx(4) = (lambda/2)*(Nf)*pVec(4)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(4);        % ab
            dx(5) = (lambda/2)*(Nf)*pVec(5)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(5);    % aa
            dx(6) = (lambda/2)*(Nf)*pVec(6)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(6);    % bb
            % females
            dx(7) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(7);                % AA
            dx(8) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(8);        % Aa
            dx(9) = (lambda/2)*(Nf)*pVec(3)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(9);        % Ab
            dx(10) = (lambda/2)*(Nf)*pVec(4)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(10);      % ab
            dx(11) = (lambda/2)*(Nf)*pVec(5)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(11);  % aa
            dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(12);  % bb
        elseif (lethality_case == "FSL") % female-sex lethality
            % larvae for each genotype
            % males
            dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);                % AA
            dx(2) = (lambda/2)*(Nf)*pVec(2)*c_a - (muX+m+log(1+(alpha*Nl)^beta))*x(2);            % Aa
            dx(3) = (lambda/2)*(Nf)*pVec(3)*c_a - (muX+m+log(1+(alpha*Nl)^beta))*x(3);            % Ab
            dx(4) = (lambda/2)*(Nf)*pVec(4)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(4);        % ab
            dx(5) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(5);        % aa
            dx(6) = (lambda/2)*(Nf)*pVec(6)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(6);        % bb
            % females
            dx(7) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(7);                % AA
            dx(8) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(8);        % Aa
            dx(9) = (lambda/2)*(Nf)*pVec(3)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(9);        % Ab
            dx(10) = (lambda/2)*(Nf)*pVec(4)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(10);      % ab
            dx(11) = (lambda/2)*(Nf)*pVec(5)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(11);  % aa
            dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(12);  % bb            
        else
            error("Lethality not recognized. Note: Male-sex lethality not included.")
        end

        % adults for each genotype
        % males
        dx(13) = m*x(1) - muZ*x(13);    % AA male
        dx(14) = m*x(2) - muZ*x(14);    % Aa male
        dx(15) = m*x(3) - muZ*x(15);    % Ab male
        dx(16) = m*x(4) - muZ*x(16);    % ab male
        dx(17) = m*x(5) - muZ*x(17);    % aa male 
        dx(18) = m*x(6) - muZ*x(18);    % bb male
        % females
        dx(19) = m*x(7) - muY*x(19);     % AA female
        dx(20) = m*x(8) - muY*x(20);     % Aa female
        dx(21) = m*x(9) - muY*x(21);     % Ab female
        dx(22) = m*x(10) - muY*x(22);    % ab female
        dx(23) = m*x(11) - muY*x(23);    % aa female 
        dx(24) = m*x(12) - muY*x(24);    % bb female
    else % fitness costs are late-acting
        if (lethality_case == "BSL") % bi-sex lethality
            % adults for each genotype
            % males
            dx(13) = m*x(1) - muZ*x(13);                % AA male
            dx(14) = m*x(2)*c_a*c_t - muZ*x(14);        % Aa male
            dx(15) = m*x(3)*c_a*c_t - muZ*x(15);        % Ab male
            dx(16) = m*x(4)*(c_a)^2 - muZ*x(16);        % ab male
            dx(17) = m*x(5)*(c_a*c_t)^2 - muZ*x(17);    % aa male 
            dx(18) = m*x(6)*(c_a*c_t)^2 - muZ*x(18);    % bb male
            % females
            dx(19) = m*x(7) - muY*x(19);                % AA female
            dx(20) = m*x(8)*c_a*c_t - muY*x(20);        % Aa female
            dx(21) = m*x(9)*c_a*c_t - muY*x(21);        % Ab female
            dx(22) = m*x(10)*(c_a)^2 - muY*x(22);       % ab female
            dx(23) = m*x(11)*(c_a*c_t)^2 - muY*x(23);   % aa female 
            dx(24) = m*x(12)*(c_a*c_t)^2 - muY*x(24);   % bb female
        elseif (lethality_case == "FSL") % female-sex lethality
            % adults for each genotype
            % males
            dx(13) = m*x(1) - muZ*x(13);                % AA male
            dx(14) = m*x(2)*c_a - muZ*x(14);            % Aa male
            dx(15) = m*x(3)*c_a - muZ*x(15);            % Ab male
            dx(16) = m*x(4)*(c_a)^2 - muZ*x(16);        % ab male
            dx(17) = m*x(5)*(c_a)^2 - muZ*x(17);        % aa male 
            dx(18) = m*x(6)*(c_a)^2 - muZ*x(18);        % bb male
            % females
            dx(19) = m*x(7) - muY*x(19);                % AA female
            dx(20) = m*x(8)*c_a*c_t - muY*x(20);        % Aa female
            dx(21) = m*x(9)*c_a*c_t - muY*x(21);        % Ab female
            dx(22) = m*x(10)*(c_a)^2 - muY*x(22);       % ab female
            dx(23) = m*x(11)*(c_a*c_t)^2 - muY*x(23);   % aa female 
            dx(24) = m*x(12)*(c_a*c_t)^2 - muY*x(24);   % bb female  
        else
            error("Lethality not recognized. Note: Male-sex lethality not included.")
        end

        % larvae for each genotype
        % males
        dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);       % AA
        dx(2) = (lambda/2)*(Nf)*pVec(2) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);       % Aa
        dx(3) = (lambda/2)*(Nf)*pVec(3) - (muX+m+log(1+(alpha*Nl)^beta))*x(3);       % Ab
        dx(4) = (lambda/2)*(Nf)*pVec(4) - (muX+m+log(1+(alpha*Nl)^beta))*x(4);       % ab
        dx(5) = (lambda/2)*(Nf)*pVec(5) - (muX+m+log(1+(alpha*Nl)^beta))*x(5);       % aa
        dx(6) = (lambda/2)*(Nf)*pVec(6) - (muX+m+log(1+(alpha*Nl)^beta))*x(6);       % bb
        % females
        dx(7) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(7);       % AA
        dx(8) = (lambda/2)*(Nf)*pVec(2) - (muX+m+log(1+(alpha*Nl)^beta))*x(8);       % Aa
        dx(9) = (lambda/2)*(Nf)*pVec(3) - (muX+m+log(1+(alpha*Nl)^beta))*x(9);       % Ab
        dx(10) = (lambda/2)*(Nf)*pVec(4) - (muX+m+log(1+(alpha*Nl)^beta))*x(10);     % ab
        dx(11) = (lambda/2)*(Nf)*pVec(5) - (muX+m+log(1+(alpha*Nl)^beta))*x(11);     % aa
        dx(12) = (lambda/2)*(Nf)*pVec(6) - (muX+m+log(1+(alpha*Nl)^beta))*x(12);     % bb
    end % end of if lethality_type statement
    dx = dx';
end

% two locus UD code
function dx = DDGD_ii_2LUD(t,x,params)
    % Case study for two-locus engineered underdominance (male-sex
    % lethality is ignored)
    % ---------------------------------------------------------------------
    
    % function simulating the generalized logistic case---note that the
    % function is programmed for beta+1 and not beta


    lambda      = params.lambda;
    muX         = params.muX;
    muZ         = params.muZ;
    muY         = params.muY;
    m           = params.m;
    s_a         = params.s_a;
    s_t         = params.s_t;

    beta        = params.beta;
    alpha       = params.alpha;
    
    lethality_case = params.lethality_case;
    lethality_type = params.lethality_type;
    supp           = params.supp;
    
    % compute relative fitnesses for ambient carriage and transgene
    c_a = 1 - s_a;
    c_t = 1 - s_t;
    
    % memory allocation
    dx=zeros(1,36);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change beta so that things make sense
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta = beta - 1; 
    
    Nl = sum(x(1:18)); % mosquito larvae pop
    Nm = sum(x(19:27)); % mosquito male pop
    Nf = sum(x(28:36)); % mosquito female pop

    % [AABB_larvae_male, AABB_larvae_female    2
    %  AABb_larvae_male, AABb_larvae_female    4
    %  AAbb_larvae_male, AAbb_larvae_female    6
    %  AaBB_larvae_male, AaBB_larvae_female    8
    %  AaBb_larvae_male, AaBb_larvae_female    10
    %  Aabb_larvae_male, Aabb_larvae_female    12
    %  aaBB_larvae_male, aaBB_larvae_female    14
    %  aaBb_larvae_male, aaBb_larvae_female    16
    %  aabb_larvae_male, aabb_larvae_female    18
    %  AABB_adult_males, AABb_adult_males      20
    %  AAbb_adult_males, AaBB_adult_males      22
    %  AaBb_adult_males, Aabb_adult_males      24
    %  aaBB_adult_males, aaBb_adult_males      26
    %  aabb_adult_males, AABB_adult_females    28
    %  AABb_adult_females, AAbb_adult_females  30
    %  AaBB_adult_females, AaBb_adult_females  32
    %  Aabb_adult_females, aaBB_adult_females  34
    %  aaBb_adult_females, aabb_adult_females] 36
    
    % calculate probabilities of offspring having each genotype
    pVec = zeros(1,9);
    pVec(1) = (1/(Nm*Nf))*(x(28)*x(19) + (1/2)*(x(28)*x(22) + x(19)*x(31)) + (1/2)*(x(28)*x(20) + x(19)*x(29)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + (1/4)*x(31)*x(22) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/4)*x(29)*x(20) + (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/16)*(x(32)*x(23)));
    pVec(4) = (1/(Nm*Nf))*((1/2)*(x(28)*x(22) + x(19)*x(31)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + (x(28)*x(25) + x(19)*x(34)) + (1/2)*(x(28)*x(26) + x(19)*x(35)) + (1/2)*x(31)*x(22) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/4)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(25) + x(22)*x(34)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(25) + x(20)*x(34)) + (1/4)*(x(29)*x(26) + x(20)*x(35)) + (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/8)*(x(32)*x(26) + x(23)*x(35)));
    pVec(2) = (1/(Nm*Nf))*((1/2)*(x(28)*x(20) + x(19)*x(29)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + (x(28)*x(21) + x(19)*x(30)) + (1/2)*(x(28)*x(24) + x(19)*x(33)) + (1/2)*x(29)*x(20) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(21) + x(22)*x(30)) + (1/4)*(x(31)*x(24) + x(22)*x(33)) + (1/4)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(21) + x(20)*x(30)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(24) + x(23)*x(33)));
    pVec(5) = (1/(Nm*Nf))*((1/4)*(x(28)*x(23) + x(19)*x(32)) + (1/2)*(x(28)*x(26) + x(19)*x(35))  + (1/2)*(x(28)*x(24) + x(19)*x(33)) + (x(28)*x(27) + x(19)*x(36)) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/4)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(21) + x(22)*x(30)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/2)*(x(31)*x(24) + x(22)*x(33)) + (1/2)*(x(31)*x(27) + x(22)*x(36)) + (1/4)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(25) + x(20)*x(34)) + (1/2)*(x(29)*x(26) + x(20)*x(35)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/2)*(x(29)*x(27) + x(20)*x(36)) + (1/4)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/4)*(x(32)*x(26) + x(23)*x(35)) + (1/4)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + (x(32)*x(21) + x(25)*x(30)) + (1/2)*(x(34)*x(24) + x(25)*x(33)) + (1/2)*(x(30)*x(26) + x(21)*x(35)) + (1/4)*(x(35)*x(24) + x(26)*x(33)));
    pVec(7) = (1/(Nm*Nf))*((1/4)*x(31)*x(22) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(25) + x(22)*x(34)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/16)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + x(34)*x(25) + (1/2)*(x(34)*x(26) + x(25)*x(35)) + (1/4)*x(35)*x(26));
    pVec(3) = (1/(Nm*Nf))*((1/4)*x(29)*x(20) + (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(21) + x(20)*x(30)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/16)*x(32)*x(23) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + x(30)*x(21) + (1/2)*(x(30)*x(24) + x(21)*x(33)) + (1/4)*x(33)*x(24));
    pVec(8) = (1/(Nm*Nf))*((1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/4)*(x(31)*x(24) + x(22)*x(33)) + (1/2)*(x(31)*x(27) + x(22)*x(36)) + (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/4)*(x(32)*x(26) + x(23)*x(35)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + (1/2)*(x(34)*x(26) + x(25)*x(35))  + (1/2)*(x(34)*x(24) + x(25)*x(33))  + (x(34)*x(27) + x(25)*x(36)) + (1/2)*x(35)*x(26) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + (1/2)*(x(35)*x(27) + x(26)*x(36)));
    pVec(6) = (1/(Nm*Nf))*((1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/4)*(x(29)*x(26) + x(20)*x(35)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/2)*(x(29)*x(27) + x(20)*x(36)) + (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + (1/4)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + (1/2)*(x(30)*x(26) + x(21)*x(35)) + (1/2)*(x(30)*x(24) + x(21)*x(33)) + (x(30)*x(27) + x(21)*x(36)) + (1/2)*x(33)*x(24) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + (1/2)*(x(33)*x(27) + x(24)*x(36)));
    pVec(9) = (1/(Nm*Nf))*((1/16)*x(32)*x(23) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + (1/4)*x(35)*x(26) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + (1/2)*(x(35)*x(27) + x(26)*x(36)) + (1/4)*x(33)*x(24) + (1/2)*(x(33)*x(27) + x(24)*x(36)) + x(36)*x(27));

    % pVec(1) = (1/(Nm*Nf))*(x(28)*x(19) + (1/2)*(x(28)*x(22) + x(19)*x(31)) + (1/2)*(x(28)*x(20) + x(19)*x(29)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + (1/4)*x(31)*x(22) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/4)*x(29)*x(20) + (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/16)*(x(32)*x(23)));
    % pVec(4) = (1/(Nm*Nf))*((1/2)*(x(28)*x(22) + x(19)*x(31)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + (x(28)*x(25) + x(19)*x(34)) + (1/2)*(x(28)*x(26) + x(19)*x(35)) + (1/2)*x(31)*x(22) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/4)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(25) + x(22)*x(34)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(25) + x(20)*x(34)) + (1/4)*(x(29)*x(26) + x(20)*x(35)) + (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/8)*(x(32)*x(26) + x(23)*x(35)));
    % pVec(2) = (1/(Nm*Nf))*((1/2)*(x(28)*x(20) + x(19)*x(29)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + (x(28)*x(21) + x(19)*x(30)) + (1/2)*(x(28)*x(24) + x(19)*x(33)) + (1/2)*x(29)*x(20) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(21) + x(22)*x(30)) + (1/4)*(x(31)*x(24) + x(22)*x(33)) + (1/4)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(21) + x(20)*x(30)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(24) + x(23)*x(33)));
    % pVec(5) = (1/(Nm*Nf))*((1/4)*(x(28)*x(23) + x(19)*x(32)) + (1/2)*(x(28)*x(26) + x(19)*x(35))  + (1/2)*(x(28)*x(24) + x(19)*x(33)) + (x(28)*x(27) + x(19)*x(36)) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/4)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(21) + x(22)*x(30)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/2)*(x(31)*x(24) + x(22)*x(33)) + (1/2)*(x(31)*x(27) + x(22)*x(36)) + (1/4)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(25) + x(20)*x(34)) + (1/2)*(x(29)*x(26) + x(20)*x(35)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/2)*(x(29)*x(27) + x(20)*x(36)) + (1/4)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/4)*(x(32)*x(26) + x(23)*x(35)) + (1/4)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + (x(32)*x(21) + x(25)*x(30)) + (1/2)*(x(34)*x(24) + x(25)*x(33)) + (1/2)*(x(30)*x(26) + x(21)*x(35)) + (1/4)*(x(35)*x(24) + x(26)*x(33)));
    % pVec(7) = (1/(Nm*Nf))*((1/4)*x(31)*x(22) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(25) + x(22)*x(34)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/16)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + x(34)*x(25) + (1/2)*(x(34)*x(26) + x(25)*x(35)) + (1/4)*x(35)*x(26));
    % pVec(3) = (1/(Nm*Nf))*((1/4)*x(29)*x(20) + (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(21) + x(20)*x(30)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/16)*x(32)*x(23) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + x(30)*x(21) + (1/2)*(x(30)*x(24) + x(21)*x(33)) + (1/4)*x(33)*x(24));
    % pVec(8) = (1/(Nm*Nf))*((1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/4)*(x(31)*x(24) + x(22)*x(33)) + (1/2)*(x(31)*x(27) + x(22)*x(36)) + (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/4)*(x(32)*x(26) + x(23)*x(35)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + (1/8)*(x(32)*x(27) + x(23)*x(36)) + (1/2)*(x(34)*x(26) + x(25)*x(35))  + (1/2)*(x(34)*x(24) + x(25)*x(33))  + (x(34)*x(27) + x(25)*x(36)) + (1/2)*x(35)*x(26) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + (1/2)*(x(35)*x(27) + x(26)*x(36)));
    % pVec(6) = (1/(Nm*Nf))*((1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/4)*(x(29)*x(26) + x(20)*x(35)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/2)*(x(29)*x(27) + x(20)*x(36)) + (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + (1/4)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + (1/2)*(x(30)*x(26) + x(21)*x(35)) + (1/2)*(x(30)*x(24) + x(21)*x(33)) + (x(30)*x(27) + x(21)*x(36)) + (1/2)*x(33)*x(24) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + (1/2)*(x(33)*x(27) + x(24)*x(36)));
    % pVec(9) = (1/(Nm*Nf))*((1/16)*x(32)*x(23) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + (1/4)*x(35)*x(26) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + (1/2)*(x(35)*x(27) + x(26)*x(36)) + (1/4)*x(33)*x(24) + (1/2)*(x(33)*x(27) + x(24)*x(36)) + x(36)*x(27));
    % % normalize to account for roundoff
    % pVec = pVec/sum(pVec); 


    
    % pVec(1) = (1/(Nm*Nf))*(x(19)*x(28)+(1/2)*(x(28)*x(22)+x(19)*x(31)) + ...
    %     (1/2)*(x(28)*x(20)+x(19)*x(29)) + 0.25*(x(28)*x(23)+x(19)*x(32)) + ...
    %     0.25*x(22)*x(31) + 0.25*(x(20)*x(31) + x(22)*x(29)) + ...
    %     (1/8)*(x(23)*x(31) + x(32)*x(22)) + (1/4)*x(20)*x(29) + ...
    %     (1/8)*(x(23)*x(29) + x(32)*x(20)) + (1/16)*x(23)*x(32));
    % pVec(2) = (1/(Nm*Nf))*((1/2)*(x(28)*x(20) + x(29)*x(19)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + ...
    %     (x(28)*x(21) + x(19)*x(30)) + (1/2)*(x(28)*x(24) + x(19)*x(33)) + ...
    %     (1/2)*x(20)*x(29) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + ...
    %     (1/2)*(x(31)*x(21) + x(22)*x(30)) + (1/4)*(x(31)*x(24) + x(22)*x(32)) + ...
    %     (1/4)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(21) + x(20)*x(30)) + ...
    %     (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/8)*x(23)*x(32) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + ...
    %     (1/8)*(x(32)*x(24)+x(23)*x(32))); 
    % pVec(3) = (1/(Nm*Nf))*((1/4)*x(29)*x(20) + (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(21)+x(20)*x(30)) + ...
    %     (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/16)*x(32)*x(23) + ...
    %     (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + x(30)*x(21) + ...
    %     (1/2)*(x(30)*x(24) + x(21)*x(33)) + (1/4)*x(33)*x(24));
    % pVec(4) = (1/(Nm*Nf))*((1/2)*(x(28)*x(22) + x(19)*x(31)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + ...
    %     (x(28)*x(25) + x(19)*x(34)) + (1/2)*(x(28)*x(26) + x(19)*x(35)) + ...
    %     (1/2)*x(31)*x(22) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/4)*(x(31)*x(23) + x(22)*x(32)) + ...
    %     (1/2)*(x(31)*x(25) + x(22)*x(34)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + ...
    %     (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(25) + x(20)*x(34)) + ...
    %     (1/4)*(x(29)*x(26) + x(20)*x(35)) + (1/8)*x(23)*x(32) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + ...
    %     (1/8)*(x(32)*x(26) + x(23)*x(35)));
    % pVec(5) = (1/(Nm*Nf))*((1/4)*(x(28)*x(23) + x(19)*x(32)) + (1/2)*(x(28)*x(26) + x(19)*x(35)) + ...
    %     (1/2)*(x(28)*x(24) + x(19)*x(33)) + (x(28)*x(27) + x(19)*x(36)) + ...
    %     (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/4)*(x(31)*x(23) + x(22)*x(32)) + ...
    %     (1/2)*(x(31)*x(21) + x(22)*x(30)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + ...
    %     (1/2)*(x(31)*x(24) + x(22)*x(33)) + (1/2)*(x(31)*x(27) + x(22)*x(36)) + ...
    %     (1/4)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(25) + x(20)*x(34)) + ...
    %     (1/2)*(x(29)*x(26) + x(20)*x(35)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + ...
    %     (1/2)*(x(29)*x(27) + x(20)*x(36)) + (1/4)*x(23)*x(32) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + ...
    %     (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/4)*(x(32)*x(26) + x(23)*x(35)) + ...
    %     (1/4)*(x(32)*x(24) + x(23)*x(32)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + ...
    %     (x(34)*x(21) + x(25)*x(30)) + (1/2)*(x(34)*x(24) + x(25)*x(33)) + ...
    %     (1/2)*(x(30)*x(26) + x(21)*x(35)) + (1/4)*(x(35)*x(24) + x(26)*x(33)));
    % pVec(6) = (1/(Nm*Nf))*((1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/4)*(x(29)*x(26) + x(20)*x(35)) + ...
    %     (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/2)*(x(29)*x(27) + x(20)*x(36)) + ...
    %     (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + ...
    %     (1/4)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + ...
    %     (1/2)*(x(30)*x(26) + x(21)*x(35)) + (1/2)*(x(30)*x(24) + x(21)*x(33)) + (x(30)*x(27) + x(21)*x(36)) + ...
    %     (1/2)*x(33)*x(24) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + ...
    %     (1/2)*(x(33)*x(27) + x(24)*x(36)));
    % pVec(7) = (1/(Nm*Nf)) *((1/4)*x(22)*x(31) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(25) + x(22)*x(34)) + ...
    %     (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/16)*x(23)*x(32) + ...
    %     (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + x(34)*x(25) + ...
    %     (1/2)*(x(34)*x(26) + x(25)*x(35)) + (1/4)*x(35)*x(26));
    % pVec(8) = (1/(Nm*Nf))*((1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + ...
    %     (1/4)*(x(31)*x(24) + x(22)*x(33)) + (1/2)*(x(31)*x(27) + x(22)*x(36)) + ...
    %     (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/4)*(x(32)*x(26) + x(23)*x(35)) + ...
    %     (1/8)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + ...
    %     (1/2)*(x(34)*x(26) + x(25)*x(35)) + (1/2)*(x(34)*x(24) + x(25)*x(33)) + (x(34)*x(27) + x(25)*x(36)) + ...
    %     (1/2)*x(35)*x(26) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + ...
    %     (1/2)*(x(35)*x(27) + x(26)*x(36)));
    % pVec(9) = (1/(Nm*Nf))*((1/16)*x(32)*x(23) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + ...
    %     (1/4)*(x(32)*x(27) + x(23)*x(36)) + (1/4)*x(35)*x(26) + ...
    %     (1/4)*(x(35)*x(24) + x(26)*x(33)) + (1/2)*(x(35)*x(27) + x(26)*x(36)) + (1/4)*x(33)*x(24) + ...
    %     (1/2)*(x(33)*x(27) + x(24)*x(36)) + x(36)*x(27));
    
    %%% disp(abs(sum(pVec)))
    % point check---probabilities sum to 1
    % if abs(sum(pVec) - 1) > 10^(-9)
    %     error("Error: Offspring genotype probabilities do not sum to 1!")
    % end
    % point check---probabilities are positive
    if any(pVec < 0)
        error("Error: Offspring genotype probabilities are negative!")
    end
    
    %%%%%%%%%%%%% coupled genetic/population dynamics model %%%%%%%%%%%%% 
    if (lethality_type == "EA")
        if (lethality_case == "BSL") % bi-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);               % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(3);       % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(5);   % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(7);       % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(9);       % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3)*c_t - (muX+m+alpha*Nl^(beta+1))*x(11); % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(13); % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3)*c_t - (muX+m+alpha*Nl^(beta+1))*x(15); % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(17);     % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(2);               % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(4);       % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(6);   % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(8);       % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(10);     % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3)*c_t - (muX+m+alpha*Nl^(beta+1))*x(12); % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(14); % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3)*c_t - (muX+m+alpha*Nl^(beta+1))*x(16); % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(18);     % aabb
            else % strongly suppressed 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);               % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(3);       % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(5);   % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(7);       % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(9);       % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(11);     % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(13); % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(15);     % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(17);     % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(2);               % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(4);       % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(6);   % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(8);       % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(10);     % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(12);     % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(14); % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(16);     % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(18);     % aabb
            end % end of if supp statement
        else % female-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);               % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a - (muX+m+alpha*Nl^(beta+1))*x(3);           % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a)^2 - (muX+m+alpha*Nl^(beta+1))*x(5);       % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a - (muX+m+alpha*Nl^(beta+1))*x(7);           % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(9);       % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(11);     % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a)^2 - (muX+m+alpha*Nl^(beta+1))*x(13);     % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(15);     % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(17);     % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(2);               % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(4);       % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(6);   % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(8);       % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(10);     % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3)*c_t - (muX+m+alpha*Nl^(beta+1))*x(12); % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(14); % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3)*c_t - (muX+m+alpha*Nl^(beta+1))*x(16); % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(18);     % aabb
            else % strongly suppressed 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);               % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a - (muX+m+alpha*Nl^(beta+1))*x(3);           % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a)^2 - (muX+m+alpha*Nl^(beta+1))*x(5);       % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a - (muX+m+alpha*Nl^(beta+1))*x(7);           % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(9);       % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(11);     % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a)^2 - (muX+m+alpha*Nl^(beta+1))*x(13);     % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(15);     % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(17);     % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(2);               % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(4);       % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(6);   % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(8);       % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(10);     % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(12);     % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(14); % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(16);     % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(18);     % aabb
            end % end of if supp statement
        end
        
        % adults for each genotype
        % males
        dx(19) = m*x(1) - muZ*x(19);    % AABB male
        dx(20) = m*x(3) - muZ*x(20);    % AABb male
        dx(21) = m*x(5) - muZ*x(21);    % AAbb male
        dx(22) = m*x(7) - muZ*x(22);    % AaBB male
        dx(23) = m*x(9) - muZ*x(23);    % AaBb male 
        dx(24) = m*x(11) - muZ*x(24);   % Aabb male
        dx(25) = m*x(13) - muZ*x(25);   % aaBB male
        dx(26) = m*x(15) - muZ*x(26);   % aaBb male
        dx(27) = m*x(17) - muZ*x(27);   % aabb male
        % females
        dx(28) = m*x(2) - muY*x(28);    % AABB female
        dx(29) = m*x(4) - muY*x(29);    % AABb female
        dx(30) = m*x(6) - muY*x(30);    % AAbb female
        dx(31) = m*x(8) - muY*x(31);    % AaBB female
        dx(32) = m*x(10) - muY*x(32);   % AaBb female 
        dx(33) = m*x(12) - muY*x(33);   % Aabb female
        dx(34) = m*x(14) - muY*x(34);   % aaBB female
        dx(35) = m*x(16) - muY*x(35);   % aaBb female
        dx(36) = m*x(18) - muY*x(36);   % aabb female
    else % fitness costs are late-acting
        if (lethality_case == "BSL") % bi-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a*c_t - muZ*x(20);        % AABb male
                dx(21) = m*x(5)*(c_a*c_t)^2 - muZ*x(21);    % AAbb male
                dx(22) = m*x(7)*c_a*c_t - muZ*x(22);        % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3)*c_t - muZ*x(24);   % Aabb male
                dx(25) = m*x(13)*(c_a*c_t)^2 - muZ*x(25);   % aaBB male
                dx(26) = m*x(15)*(c_a^3)*c_t - muZ*x(26);   % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3)*c_t - muY*x(33);   % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3)*c_t - muY*x(35);   % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            else % strongly suppressed 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a*c_t - muZ*x(20);        % AABb male
                dx(21) = m*x(5)*(c_a*c_t)^2 - muZ*x(21);    % AAbb male
                dx(22) = m*x(7)*c_a*c_t - muZ*x(22);        % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3) - muZ*x(24);       % Aabb male
                dx(25) = m*x(13)*(c_a*c_t)^2 - muZ*x(25);   % aaBB male
                dx(26) = m*x(15)*(c_a^3) - muZ*x(26);       % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3) - muY*x(33);       % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3) - muY*x(35);       % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            end % end of if supp statement
        else % female-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a - muZ*x(20);            % AABb male
                dx(21) = m*x(5)*(c_a)^2 - muZ*x(21);        % AAbb male
                dx(22) = m*x(7)*c_a - muZ*x(22);            % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3) - muZ*x(24);       % Aabb male
                dx(25) = m*x(13)*(c_a)^2 - muZ*x(25);       % aaBB male
                dx(26) = m*x(15)*(c_a^3) - muZ*x(26);       % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3)*c_t - muY*x(33);   % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3)*c_t - muY*x(35);   % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            else % strongly suppressed 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a - muZ*x(20);            % AABb male
                dx(21) = m*x(5)*(c_a)^2 - muZ*x(21);        % AAbb male
                dx(22) = m*x(7)*c_a - muZ*x(22);            % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3) - muZ*x(24);       % Aabb male
                dx(25) = m*x(13)*(c_a)^2 - muZ*x(25);       % aaBB male
                dx(26) = m*x(15)*(c_a^3) - muZ*x(26);       % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3) - muY*x(33);       % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3) - muY*x(35);       % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            end % end of if supp statement
        end

        % larvae for each genotype
        % males
        dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);       % AABB 
        dx(3) = (lambda/2)*(Nf)*pVec(2) - (muX+m+alpha*Nl^(beta+1))*x(3);       % AABb
        dx(5) = (lambda/2)*(Nf)*pVec(3) - (muX+m+alpha*Nl^(beta+1))*x(5);       % AAbb
        dx(7) = (lambda/2)*(Nf)*pVec(4) - (muX+m+alpha*Nl^(beta+1))*x(7);       % AaBB
        dx(9) = (lambda/2)*(Nf)*pVec(5) - (muX+m+alpha*Nl^(beta+1))*x(9);       % AaBb
        dx(11) = (lambda/2)*(Nf)*pVec(6) - (muX+m+alpha*Nl^(beta+1))*x(11);     % Aabb
        dx(13) = (lambda/2)*(Nf)*pVec(7) - (muX+m+alpha*Nl^(beta+1))*x(13);     % aaBB
        dx(15) = (lambda/2)*(Nf)*pVec(8) - (muX+m+alpha*Nl^(beta+1))*x(15);     % aaBb
        dx(17) = (lambda/2)*(Nf)*pVec(9) - (muX+m+alpha*Nl^(beta+1))*x(17);     % aabb
        % females
        dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(2);       % AABB 
        dx(4) = (lambda/2)*(Nf)*pVec(2) - (muX+m+alpha*Nl^(beta+1))*x(4);       % AABb
        dx(6) = (lambda/2)*(Nf)*pVec(3) - (muX+m+alpha*Nl^(beta+1))*x(6);       % AAbb
        dx(8) = (lambda/2)*(Nf)*pVec(4) - (muX+m+alpha*Nl^(beta+1))*x(8);       % AaBB
        dx(10) = (lambda/2)*(Nf)*pVec(5) - (muX+m+alpha*Nl^(beta+1))*x(10);     % AaBb
        dx(12) = (lambda/2)*(Nf)*pVec(6) - (muX+m+alpha*Nl^(beta+1))*x(12);     % Aabb
        dx(14) = (lambda/2)*(Nf)*pVec(7) - (muX+m+alpha*Nl^(beta+1))*x(14);     % aaBB
        dx(16) = (lambda/2)*(Nf)*pVec(8) - (muX+m+alpha*Nl^(beta+1))*x(16);     % aaBb
        dx(18) = (lambda/2)*(Nf)*pVec(9) - (muX+m+alpha*Nl^(beta+1))*x(18);     % aabb
    end % end of if lethality_type statement
    
    dx = dx';

end

function dx = DDGD_iii_2LUD(t,x,params)
    % Case study for two-locus engineered underdominance (male-sex
    % lethality is ignored)
    % ---------------------------------------------------------------------

    lambda      = params.lambda;
    muX         = params.muX;
    muZ         = params.muZ;
    muY         = params.muY;
    m           = params.m;
    s_a         = params.s_a;
    s_t         = params.s_t;

    beta        = params.beta;
    alpha       = params.alpha;
    
    lethality_case = params.lethality_case;
    lethality_type = params.lethality_type;
    supp           = params.supp;
    
    % compute relative fitnesses for ambient carriage and transgene
    c_a = 1 - s_a;
    c_t = 1 - s_t;
    
    % memory allocation
    dx=zeros(1,36);
    
    Nl = sum(x(1:18)); % mosquito larvae pop
    Nm = sum(x(19:27)); % mosquito male pop
    Nf = sum(x(28:36)); % mosquito female pop

    % [AABB_larvae_male, AABB_larvae_female    2
    %  AABb_larvae_male, AABb_larvae_female    4
    %  AAbb_larvae_male, AAbb_larvae_female    6
    %  AaBB_larvae_male, AaBB_larvae_female    8
    %  AaBb_larvae_male, AaBb_larvae_female    10
    %  Aabb_larvae_male, Aabb_larvae_female    12
    %  aaBB_larvae_male, aaBB_larvae_female    14
    %  aaBb_larvae_male, aaBb_larvae_female    16
    %  aabb_larvae_male, aabb_larvae_female    18
    %  AABB_adult_males, AABb_adult_males      20
    %  AAbb_adult_males, AaBB_adult_males      22
    %  AaBb_adult_males, Aabb_adult_males      24
    %  aaBB_adult_males, aaBb_adult_males      26
    %  aabb_adult_males, AABB_adult_females    28
    %  AABb_adult_females, AAbb_adult_females  30
    %  AaBB_adult_females, AaBb_adult_females  32
    %  Aabb_adult_females, aaBB_adult_females  34
    %  aaBb_adult_females, aabb_adult_females] 36    
    
    % calculate probabilities of offspring having each genotype
    pVec = zeros(1,9);
    pVec(1) = (1/(Nm*Nf))*(x(28)*x(19) + (1/2)*(x(28)*x(22) + x(19)*x(31)) + (1/2)*(x(28)*x(20) + x(19)*x(29)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + (1/4)*x(31)*x(22) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/4)*x(29)*x(20) + (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/16)*(x(32)*x(23)));
    pVec(4) = (1/(Nm*Nf))*((1/2)*(x(28)*x(22) + x(19)*x(31)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + (x(28)*x(25) + x(19)*x(34)) + (1/2)*(x(28)*x(26) + x(19)*x(35)) + (1/2)*x(31)*x(22) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/4)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(25) + x(22)*x(34)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(25) + x(20)*x(34)) + (1/4)*(x(29)*x(26) + x(20)*x(35)) + (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/8)*(x(32)*x(26) + x(23)*x(35)));
    pVec(2) = (1/(Nm*Nf))*((1/2)*(x(28)*x(20) + x(19)*x(29)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + (x(28)*x(21) + x(19)*x(30)) + (1/2)*(x(28)*x(24) + x(19)*x(33)) + (1/2)*x(29)*x(20) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(21) + x(22)*x(30)) + (1/4)*(x(31)*x(24) + x(22)*x(33)) + (1/4)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(21) + x(20)*x(30)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(24) + x(23)*x(33)));
    pVec(5) = (1/(Nm*Nf))*((1/4)*(x(28)*x(23) + x(19)*x(32)) + (1/2)*(x(28)*x(26) + x(19)*x(35))  + (1/2)*(x(28)*x(24) + x(19)*x(33)) + (x(28)*x(27) + x(19)*x(36)) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/4)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(21) + x(22)*x(30)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/2)*(x(31)*x(24) + x(22)*x(33)) + (1/2)*(x(31)*x(27) + x(22)*x(36)) + (1/4)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(25) + x(20)*x(34)) + (1/2)*(x(29)*x(26) + x(20)*x(35)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/2)*(x(29)*x(27) + x(20)*x(36)) + (1/4)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/4)*(x(32)*x(26) + x(23)*x(35)) + (1/4)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + (x(32)*x(21) + x(25)*x(30)) + (1/2)*(x(34)*x(24) + x(25)*x(33)) + (1/2)*(x(30)*x(26) + x(21)*x(35)) + (1/4)*(x(35)*x(24) + x(26)*x(33)));
    pVec(7) = (1/(Nm*Nf))*((1/4)*x(31)*x(22) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(25) + x(22)*x(34)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/16)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + x(34)*x(25) + (1/2)*(x(34)*x(26) + x(25)*x(35)) + (1/4)*x(35)*x(26));
    pVec(3) = (1/(Nm*Nf))*((1/4)*x(29)*x(20) + (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(21) + x(20)*x(30)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/16)*x(32)*x(23) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + x(30)*x(21) + (1/2)*(x(30)*x(24) + x(21)*x(33)) + (1/4)*x(33)*x(24));
    pVec(8) = (1/(Nm*Nf))*((1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/4)*(x(31)*x(24) + x(22)*x(33)) + (1/2)*(x(31)*x(27) + x(22)*x(36)) + (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/4)*(x(32)*x(26) + x(23)*x(35)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + (1/2)*(x(34)*x(26) + x(25)*x(35))  + (1/2)*(x(34)*x(24) + x(25)*x(33))  + (x(34)*x(27) + x(25)*x(36)) + (1/2)*x(35)*x(26) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + (1/2)*(x(35)*x(27) + x(26)*x(36)));
    pVec(6) = (1/(Nm*Nf))*((1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/4)*(x(29)*x(26) + x(20)*x(35)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/2)*(x(29)*x(27) + x(20)*x(36)) + (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + (1/4)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + (1/2)*(x(30)*x(26) + x(21)*x(35)) + (1/2)*(x(30)*x(24) + x(21)*x(33)) + (x(30)*x(27) + x(21)*x(36)) + (1/2)*x(33)*x(24) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + (1/2)*(x(33)*x(27) + x(24)*x(36)));
    pVec(9) = (1/(Nm*Nf))*((1/16)*x(32)*x(23) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + (1/4)*x(35)*x(26) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + (1/2)*(x(35)*x(27) + x(26)*x(36)) + (1/4)*x(33)*x(24) + (1/2)*(x(33)*x(27) + x(24)*x(36)) + x(36)*x(27));
    
    % disp(abs(sum(pVec)))
    % point check---probabilities sum to 1
    % if abs(sum(pVec) - 1) > 10^(-9)
    %     error("Error: Offspring genotype probabilities do not sum to 1!")
    % end
    % point check---probabilities are positive
    if any(pVec < 0)
        error("Error: Offspring genotype probabilities are negative!")
    end
    
    %%%%%%%%%%%%% coupled genetic/population dynamics model %%%%%%%%%%%%% 
    if (lethality_type == "EA")
        if (lethality_case == "BSL") % bi-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);               % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(3);       % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(5);   % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(7);       % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(9);       % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3)*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(11); % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(13); % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3)*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(15); % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(17);     % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);               % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(4);       % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(6);   % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(8);       % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(10);     % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3)*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(12); % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(14); % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3)*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(16); % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(18);     % aabb
            else % strongly suppressed 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);               % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(3);       % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(5);   % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(7);       % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(9);       % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(11);     % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(13); % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(15);     % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(17);     % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);               % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(4);       % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(6);   % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(8);       % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(10);     % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(12);     % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(14); % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(16);     % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(18);     % aabb
            end % end of if supp statement
        else % female-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);               % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a - (muX+m+log(1+(alpha*Nl)^beta))*x(3);           % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(5);       % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a - (muX+m+log(1+(alpha*Nl)^beta))*x(7);           % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(9);       % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(11);     % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(13);     % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(15);     % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(17);     % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);               % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(4);       % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(6);   % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(8);       % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(10);     % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3)*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(12); % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(14); % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3)*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(16); % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(18);     % aabb
            else % strongly suppressed 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);               % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a - (muX+m+log(1+(alpha*Nl)^beta))*x(3);           % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(5);       % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a - (muX+m+log(1+(alpha*Nl)^beta))*x(7);           % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(9);       % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(11);     % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(13);     % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(15);     % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(17);     % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);               % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(4);       % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(6);   % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(8);       % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(10);     % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(12);     % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(14); % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(16);     % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(18);     % aabb
            end % end of if supp statement
        end
        
        % adults for each genotype
        % males
        dx(19) = m*x(1) - muZ*x(19);    % AABB male
        dx(20) = m*x(3) - muZ*x(20);    % AABb male
        dx(21) = m*x(5) - muZ*x(21);    % AAbb male
        dx(22) = m*x(7) - muZ*x(22);    % AaBB male
        dx(23) = m*x(9) - muZ*x(23);    % AaBb male 
        dx(24) = m*x(11) - muZ*x(24);   % Aabb male
        dx(25) = m*x(13) - muZ*x(25);   % aaBB male
        dx(26) = m*x(15) - muZ*x(26);   % aaBb male
        dx(27) = m*x(17) - muZ*x(27);   % aabb male
        % females
        dx(28) = m*x(2) - muY*x(28);    % AABB female
        dx(29) = m*x(4) - muY*x(29);    % AABb female
        dx(30) = m*x(6) - muY*x(30);    % AAbb female
        dx(31) = m*x(8) - muY*x(31);    % AaBB female
        dx(32) = m*x(10) - muY*x(32);   % AaBb female 
        dx(33) = m*x(12) - muY*x(33);   % Aabb female
        dx(34) = m*x(14) - muY*x(34);   % aaBB female
        dx(35) = m*x(16) - muY*x(35);   % aaBb female
        dx(36) = m*x(18) - muY*x(36);   % aabb female
    else % fitness costs are late-acting
        if (lethality_case == "BSL") % bi-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a*c_t - muZ*x(20);        % AABb male
                dx(21) = m*x(5)*(c_a*c_t)^2 - muZ*x(21);    % AAbb male
                dx(22) = m*x(7)*c_a*c_t - muZ*x(22);        % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3)*c_t - muZ*x(24);   % Aabb male
                dx(25) = m*x(13)*(c_a*c_t)^2 - muZ*x(25);   % aaBB male
                dx(26) = m*x(15)*(c_a^3)*c_t - muZ*x(26);   % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3)*c_t - muY*x(33);   % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3)*c_t - muY*x(35);   % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            else % strongly suppressed 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a*c_t - muZ*x(20);        % AABb male
                dx(21) = m*x(5)*(c_a*c_t)^2 - muZ*x(21);    % AAbb male
                dx(22) = m*x(7)*c_a*c_t - muZ*x(22);        % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3) - muZ*x(24);       % Aabb male
                dx(25) = m*x(13)*(c_a*c_t)^2 - muZ*x(25);   % aaBB male
                dx(26) = m*x(15)*(c_a^3) - muZ*x(26);       % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3) - muY*x(33);       % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3) - muY*x(35);       % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            end % end of if supp statement
        else % female-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a - muZ*x(20);            % AABb male
                dx(21) = m*x(5)*(c_a)^2 - muZ*x(21);        % AAbb male
                dx(22) = m*x(7)*c_a - muZ*x(22);            % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3) - muZ*x(24);       % Aabb male
                dx(25) = m*x(13)*(c_a)^2 - muZ*x(25);       % aaBB male
                dx(26) = m*x(15)*(c_a^3) - muZ*x(26);       % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3)*c_t - muY*x(33);   % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3)*c_t - muY*x(35);   % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            else % strongly suppressed 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a - muZ*x(20);            % AABb male
                dx(21) = m*x(5)*(c_a)^2 - muZ*x(21);        % AAbb male
                dx(22) = m*x(7)*c_a - muZ*x(22);            % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3) - muZ*x(24);       % Aabb male
                dx(25) = m*x(13)*(c_a)^2 - muZ*x(25);       % aaBB male
                dx(26) = m*x(15)*(c_a^3) - muZ*x(26);       % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3) - muY*x(33);       % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3) - muY*x(35);       % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            end % end of if supp statement
        end

        % larvae for each genotype
        % males
        dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);       % AABB 
        dx(3) = (lambda/2)*(Nf)*pVec(2) - (muX+m+log(1+(alpha*Nl)^beta))*x(3);       % AABb
        dx(5) = (lambda/2)*(Nf)*pVec(3) - (muX+m+log(1+(alpha*Nl)^beta))*x(5);       % AAbb
        dx(7) = (lambda/2)*(Nf)*pVec(4) - (muX+m+log(1+(alpha*Nl)^beta))*x(7);       % AaBB
        dx(9) = (lambda/2)*(Nf)*pVec(5) - (muX+m+log(1+(alpha*Nl)^beta))*x(9);       % AaBb
        dx(11) = (lambda/2)*(Nf)*pVec(6) - (muX+m+log(1+(alpha*Nl)^beta))*x(11);     % Aabb
        dx(13) = (lambda/2)*(Nf)*pVec(7) - (muX+m+log(1+(alpha*Nl)^beta))*x(13);     % aaBB
        dx(15) = (lambda/2)*(Nf)*pVec(8) - (muX+m+log(1+(alpha*Nl)^beta))*x(15);     % aaBb
        dx(17) = (lambda/2)*(Nf)*pVec(9) - (muX+m+log(1+(alpha*Nl)^beta))*x(17);     % aabb
        % females
        dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);       % AABB 
        dx(4) = (lambda/2)*(Nf)*pVec(2) - (muX+m+log(1+(alpha*Nl)^beta))*x(4);       % AABb
        dx(6) = (lambda/2)*(Nf)*pVec(3) - (muX+m+log(1+(alpha*Nl)^beta))*x(6);       % AAbb
        dx(8) = (lambda/2)*(Nf)*pVec(4) - (muX+m+log(1+(alpha*Nl)^beta))*x(8);       % AaBB
        dx(10) = (lambda/2)*(Nf)*pVec(5) - (muX+m+log(1+(alpha*Nl)^beta))*x(10);     % AaBb
        dx(12) = (lambda/2)*(Nf)*pVec(6) - (muX+m+log(1+(alpha*Nl)^beta))*x(12);     % Aabb
        dx(14) = (lambda/2)*(Nf)*pVec(7) - (muX+m+log(1+(alpha*Nl)^beta))*x(14);     % aaBB
        dx(16) = (lambda/2)*(Nf)*pVec(8) - (muX+m+log(1+(alpha*Nl)^beta))*x(16);     % aaBb
        dx(18) = (lambda/2)*(Nf)*pVec(9) - (muX+m+log(1+(alpha*Nl)^beta))*x(18);     % aabb
    end % end of if lethality_type statement
    
    dx = dx';

end

% homing code
function dx = DDGD_ii_HOM(t,x,params)
    % Case study for homing drive (male-sex lethality is ignored)
    % ---------------------------------------------------------------------
    
    % function simulating the generalized logistic case---note that the
    % function is programmed for beta+1 and not beta

    lambda      = params.lambda;
    muX         = params.muX;
    muZ         = params.muZ;
    muY         = params.muY;
    m           = params.m;
    CONV_EFF    = params.CONV_EFF; 
    h           = params.h;     % dominance 
    s           = params.s;     % fitness cost

    beta        = params.beta;
    alpha       = params.alpha;
    
    lethality_case = params.lethality_case;
    lethality_type = params.lethality_type;

    % correct for numerical mumbo jumbo when pops are small
    % x(x < 0) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change beta so that things make sense
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta = beta - 1;     

    % memory allocation
    dx=zeros(1,12);
 
    % pop counts
    Nl = sum(x(1:6)); 
    Nm = sum(x(7:9));
    Nf = sum(x(10:12));
    
    % genotype probabilities, [AA, Aa, aa]
    % if either sex is extinct, set pVec to zero
    pVec = zeros(1,3);    
    if (Nm > 0) && (Nf > 0)
        pVec(1) = (1/(Nm*Nf))*(x(7)*x(10) + (1/2)*(1-CONV_EFF)*(x(7)*x(11) + x(8)*x(10)) + (1/4)*((1-CONV_EFF)^2)*x(8)*x(11));
        pVec(2) = (1/(Nm*Nf))*((x(7)*x(12) + x(9)*x(10)) + (1/2)*(1+CONV_EFF)*(x(7)*x(11) + x(8)*x(10)) + (1/2)*(1-CONV_EFF^2)*x(8)*x(11) + (1/2)*(1-CONV_EFF)*(x(8)*x(12) + x(9)*x(11)));
        pVec(3) = (1/(Nm*Nf))*(x(9)*x(12) + (1/2)*(1+CONV_EFF)*(x(8)*x(12) + x(9)*x(11)) + (1/4)*((1+CONV_EFF)^2)*x(8)*x(11));
    end
    
    % vector of relative fitnesses
    phiVec = [1, 1-h*s, 1-s];
    
    % [AA_male_juveniles, Aa_male_juveniles,      2
    %  aa_male_juveniles, AA_female_juveniles,    4
    %  Aa_female_juveniles, aa_female_juveniles,  6
    %  AA_adult_males, Aa_adult_males,            8
    %  aa_adult_males, AA_adult_females,          10
    %  Aa_adult_females, aa_adult_females]        12
    
    %%%%%%%%%%%%% coupled genetic/population dynamics model %%%%%%%%%%%%% 
    if (lethality_type == "EA") % early-acting

        if (lethality_case == "BSL")
            % larvae for each genotype
            % males
            dx(1) = (lambda/2)*(Nf)*pVec(1)*phiVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);    % AA
            dx(2) = (lambda/2)*(Nf)*pVec(2)*phiVec(2) - (muX+m+alpha*Nl^(beta+1))*x(2);    % Aa
            dx(3) = (lambda/2)*(Nf)*pVec(3)*phiVec(3) - (muX+m+alpha*Nl^(beta+1))*x(3);    % aa
            % females
            dx(4) = (lambda/2)*(Nf)*pVec(1)*phiVec(1) - (muX+m+alpha*Nl^(beta+1))*x(4);    % AA
            dx(5) = (lambda/2)*(Nf)*pVec(2)*phiVec(2) - (muX+m+alpha*Nl^(beta+1))*x(5);    % Aa
            dx(6) = (lambda/2)*(Nf)*pVec(3)*phiVec(3) - (muX+m+alpha*Nl^(beta+1))*x(6);    % aa
        elseif (lethality_case == "FSL")
            % males
            dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);              % AA
            dx(2) = (lambda/2)*(Nf)*pVec(2) - (muX+m+alpha*Nl^(beta+1))*x(2);              % Aa
            dx(3) = (lambda/2)*(Nf)*pVec(3) - (muX+m+alpha*Nl^(beta+1))*x(3);              % aa
            % females
            dx(4) = (lambda/2)*(Nf)*pVec(1)*phiVec(1) - (muX+m+alpha*Nl^(beta+1))*x(4);    % AA
            dx(5) = (lambda/2)*(Nf)*pVec(2)*phiVec(2) - (muX+m+alpha*Nl^(beta+1))*x(5);    % Aa
            dx(6) = (lambda/2)*(Nf)*pVec(3)*phiVec(3) - (muX+m+alpha*Nl^(beta+1))*x(6);    % aa            
        else
            error("Lethality not recognized. Note: Male-sex lethality not included.");
        end
    
        % adults for each genotype
        % males
        dx(7) = m*x(1) - muZ*x(7);      % AA male
        dx(8) = m*x(2) - muZ*x(8);      % Aa male
        dx(9) = m*x(3) - muZ*x(9);      % aa male
        % females
        dx(10) = m*x(4) - muY*x(10);    % AA female
        dx(11) = m*x(5) - muY*x(11);    % Aa female
        dx(12) = m*x(6) - muY*x(12);    % aa female

    elseif (lethality_type == "LA") % late-acting

        if (lethality_case == "BSL")
            % adults for each genotype
            % males
            dx(7) = m*x(1)*phiVec(1) - muZ*x(7);      % AA male
            dx(8) = m*x(2)*phiVec(2) - muZ*x(8);      % Aa male
            dx(9) = m*x(3)*phiVec(3) - muZ*x(9);      % aa male
            % females
            dx(10) = m*x(4)*phiVec(1) - muY*x(10);    % AA female
            dx(11) = m*x(5)*phiVec(2) - muY*x(11);    % Aa female
            dx(12) = m*x(6)*phiVec(3) - muY*x(12);    % aa female
        elseif (lethality_case == "FSL")
            % adults for each genotype
            % males
            dx(7) = m*x(1) - muZ*x(7);                % AA male
            dx(8) = m*x(2) - muZ*x(8);                % Aa male
            dx(9) = m*x(3) - muZ*x(9);                % aa male
            % females
            dx(10) = m*x(4)*phiVec(1) - muY*x(10);    % AA female
            dx(11) = m*x(5)*phiVec(2) - muY*x(11);    % Aa female
            dx(12) = m*x(6)*phiVec(3) - muY*x(12);    % aa female
        else
            error("Lethality not recognized. Note: Male-sex lethality not included.");
        end
    
        % larvae for each genotype
        % males
        dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);    % AA
        dx(2) = (lambda/2)*(Nf)*pVec(2) - (muX+m+alpha*Nl^(beta+1))*x(2);    % Aa
        dx(3) = (lambda/2)*(Nf)*pVec(3) - (muX+m+alpha*Nl^(beta+1))*x(3);    % aa
        % females
        dx(4) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(4);    % AA
        dx(5) = (lambda/2)*(Nf)*pVec(2) - (muX+m+alpha*Nl^(beta+1))*x(5);    % Aa
        dx(6) = (lambda/2)*(Nf)*pVec(3) - (muX+m+alpha*Nl^(beta+1))*x(6);    % aa        
    end % end of if lethality_type statement
    
    dx = dx';

end

function dx = DDGD_iii_HOM(t,x,params)
    % Case study for homing drive (male-sex lethality is ignored)
    % ---------------------------------------------------------------------
    
    % function simulating the logarithmic case

    lambda      = params.lambda;
    muX         = params.muX;
    muZ         = params.muZ;
    muY         = params.muY;
    m           = params.m;
    CONV_EFF    = params.CONV_EFF; 
    h           = params.h;     % dominance 
    s           = params.s;     % fitness cost

    beta        = params.beta;
    alpha       = params.alpha;
    
    lethality_case = params.lethality_case;
    lethality_type = params.lethality_type;

    % correct for numerical mumbo jumbo when pops are small
    x(x < 0) = 0;    

    % memory allocation
    dx=zeros(1,12);
 
    % pop counts
    Nl = sum(x(1:6)); 
    Nm = sum(x(7:9));
    Nf = sum(x(10:12));
    
    % genotype probabilities, [AA, Aa, aa]
    % if either sex is extinct, set pVec to zero    
    pVec = zeros(1,3);
    if (Nm > 0) && (Nf > 0)
        pVec(1) = (1/(Nm*Nf))*(x(7)*x(10) + (1/2)*(1-CONV_EFF)*(x(7)*x(11) + x(8)*x(10)) + (1/4)*((1-CONV_EFF)^2)*x(8)*x(11));
        pVec(2) = (1/(Nm*Nf))*((x(7)*x(12) + x(9)*x(10)) + (1/2)*(1+CONV_EFF)*(x(7)*x(11) + x(8)*x(10)) + (1/2)*(1-CONV_EFF^2)*x(8)*x(11) + (1/2)*(1-CONV_EFF)*(x(8)*x(12) + x(9)*x(11)));
        pVec(3) = (1/(Nm*Nf))*(x(9)*x(12) + (1/2)*(1+CONV_EFF)*(x(8)*x(12) + x(9)*x(11)) + (1/4)*((1+CONV_EFF)^2)*x(8)*x(11));
    end
    
    % vector of relative fitnesses
    phiVec = [1, 1-h*s, 1-s];
    
    % [AA_male_juveniles, Aa_male_juveniles,      2
    %  aa_male_juveniles, AA_female_juveniles,    4
    %  Aa_female_juveniles, aa_female_juveniles,  6
    %  AA_adult_males, Aa_adult_males,            8
    %  aa_adult_males, AA_adult_females,          10
    %  Aa_adult_females, aa_adult_females]        12
    
    %%%%%%%%%%%%% coupled genetic/population dynamics model %%%%%%%%%%%%% 
    if (lethality_type == "EA") % early-acting

        if (lethality_case == "BSL")
            % larvae for each genotype
            % males
            dx(1) = (lambda/2)*(Nf)*pVec(1)*phiVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);    % AA
            dx(2) = (lambda/2)*(Nf)*pVec(2)*phiVec(2) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);    % Aa
            dx(3) = (lambda/2)*(Nf)*pVec(3)*phiVec(3) - (muX+m+log(1+(alpha*Nl)^beta))*x(3);    % aa
            % females
            dx(4) = (lambda/2)*(Nf)*pVec(1)*phiVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(4);    % AA
            dx(5) = (lambda/2)*(Nf)*pVec(2)*phiVec(2) - (muX+m+log(1+(alpha*Nl)^beta))*x(5);    % Aa
            dx(6) = (lambda/2)*(Nf)*pVec(3)*phiVec(3) - (muX+m+log(1+(alpha*Nl)^beta))*x(6);    % aa
        elseif (lethality_case == "FSL")
            % males
            dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);              % AA
            dx(2) = (lambda/2)*(Nf)*pVec(2) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);              % Aa
            dx(3) = (lambda/2)*(Nf)*pVec(3) - (muX+m+log(1+(alpha*Nl)^beta))*x(3);              % aa
            % females
            dx(4) = (lambda/2)*(Nf)*pVec(1)*phiVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(4);    % AA
            dx(5) = (lambda/2)*(Nf)*pVec(2)*phiVec(2) - (muX+m+log(1+(alpha*Nl)^beta))*x(5);    % Aa
            dx(6) = (lambda/2)*(Nf)*pVec(3)*phiVec(3) - (muX+m+log(1+(alpha*Nl)^beta))*x(6);    % aa            
        else
            error("Lethality not recognized. Note: Male-sex lethality not included.");
        end
    
        % adults for each genotype
        % males
        dx(7) = m*x(1) - muZ*x(7);      % AA male
        dx(8) = m*x(2) - muZ*x(8);      % Aa male
        dx(9) = m*x(3) - muZ*x(9);      % aa male
        % females
        dx(10) = m*x(4) - muY*x(10);    % AA female
        dx(11) = m*x(5) - muY*x(11);    % Aa female
        dx(12) = m*x(6) - muY*x(12);    % aa female

    elseif (lethality_type == "LA") % late-acting

        if (lethality_case == "BSL")
            % adults for each genotype
            % males
            dx(7) = m*x(1)*phiVec(1) - muZ*x(7);      % AA male
            dx(8) = m*x(2)*phiVec(2) - muZ*x(8);      % Aa male
            dx(9) = m*x(3)*phiVec(3) - muZ*x(9);      % aa male
            % females
            dx(10) = m*x(4)*phiVec(1) - muY*x(10);    % AA female
            dx(11) = m*x(5)*phiVec(2) - muY*x(11);    % Aa female
            dx(12) = m*x(6)*phiVec(3) - muY*x(12);    % aa female
        elseif (lethality_case == "FSL")
            % adults for each genotype
            % males
            dx(7) = m*x(1) - muZ*x(7);                % AA male
            dx(8) = m*x(2) - muZ*x(8);                % Aa male
            dx(9) = m*x(3) - muZ*x(9);                % aa male
            % females
            dx(10) = m*x(4)*phiVec(1) - muY*x(10);    % AA female
            dx(11) = m*x(5)*phiVec(2) - muY*x(11);    % Aa female
            dx(12) = m*x(6)*phiVec(3) - muY*x(12);    % aa female
        else
            error("Lethality not recognized. Note: Male-sex lethality not included.");
        end
    
        % larvae for each genotype
        % males
        dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);    % AA
        dx(2) = (lambda/2)*(Nf)*pVec(2) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);    % Aa
        dx(3) = (lambda/2)*(Nf)*pVec(3) - (muX+m+log(1+(alpha*Nl)^beta))*x(3);    % aa
        % females
        dx(4) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(4);    % AA
        dx(5) = (lambda/2)*(Nf)*pVec(2) - (muX+m+log(1+(alpha*Nl)^beta))*x(5);    % Aa
        dx(6) = (lambda/2)*(Nf)*pVec(3) - (muX+m+log(1+(alpha*Nl)^beta))*x(6);    % aa        
    end % end of if lethality_type statement
    
    dx = dx';

end

% bellows DD function
function dx = bellows_5(t,x,a,b)
    % define mu expression
    dx = -x*log(1+(a*x)^b);
end

% gene drive, version 2 density-dependence 
function dx = DDGD_ii(t,x,params,alpha,beta,lethality_case,supp,lethality_type)
    % Case study for two-locus engineered underdominance (male-sex
    % lethality is ignored)
    % ---------------------------------------------------------------------
    
    % function simulating the generalized logistic case---note that the
    % function is programmed for beta+1 and not beta

    lambda      = params.lambda;
    muX         = params.muX;
    muZ         = params.muZ;
    muY         = params.muY;
    m           = params.m;
    a           = params.a;
    b           = params.b;
    c           = params.c;
    gamma       = params.gamma;
    M           = params.M;
    s_a         = params.s_a;
    s_t         = params.s_t;
    % N = params(13); 
    
    % compute relative fitnesses for ambient carriage and transgene
    c_a = 1 - s_a;
    c_t = 1 - s_t;
    
    % memory allocation
    dx=zeros(1,38);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change beta so that things make sense
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta = beta - 1; 
    
    Nl = sum(x(1:18)); % mosquito larvae pop
    Nm = sum(x(19:27)); % mosquito male pop
    Nf = sum(x(28:36)); % mosquito female pop
    
    % calculate probabilities of offspring having each genotype
    pVec = zeros(1,9);
    pVec(1) = (1/(Nm*Nf))*(x(19)*x(28)+(1/2)*(x(28)*x(22)+x(19)*x(31)) + ...
        (1/2)*(x(28)*x(20)+x(19)*x(29)) + 0.25*(x(28)*x(23)+x(19)*x(32)) + ...
        0.25*x(22)*x(31) + 0.25*(x(20)*x(31) + x(22)*x(29)) + ...
        (1/8)*(x(23)*x(31) + x(32)*x(22)) + (1/4)*x(20)*x(29) + ...
        (1/8)*(x(23)*x(29) + x(32)*x(20)) + (1/16)*x(23)*x(32));
    pVec(2) = (1/(Nm*Nf))*((1/2)*(x(28)*x(20) + x(29)*x(19)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + ...
        (x(28)*x(21) + x(19)*x(30)) + (1/2)*(x(28)*x(24) + x(19)*x(33)) + ...
        (1/2)*x(20)*x(29) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + ...
        (1/2)*(x(31)*x(21) + x(22)*x(30)) + (1/4)*(x(31)*x(24) + x(22)*x(32)) + ...
        (1/4)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(21) + x(20)*x(30)) + ...
        (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/8)*x(23)*x(32) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + ...
        (1/8)*(x(32)*x(24)+x(23)*x(32))); 
    pVec(3) = (1/(Nm*Nf))*((1/4)*x(29)*x(20) + (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(21)+x(20)*x(30)) + ...
        (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/16)*x(32)*x(23) + ...
        (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + x(30)*x(21) + ...
        (1/2)*(x(30)*x(24) + x(21)*x(33)) + (1/4)*x(33)*x(24));
    pVec(4) = (1/(Nm*Nf))*((1/2)*(x(28)*x(22) + x(19)*x(31)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + ...
        (x(28)*x(25) + x(19)*x(34)) + (1/2)*(x(28)*x(26) + x(19)*x(35)) + ...
        (1/2)*x(31)*x(22) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/4)*(x(31)*x(23) + x(22)*x(32)) + ...
        (1/2)*(x(31)*x(25) + x(22)*x(34)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + ...
        (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(25) + x(20)*x(34)) + ...
        (1/4)*(x(29)*x(26) + x(20)*x(35)) + (1/8)*x(23)*x(32) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + ...
        (1/8)*(x(32)*x(26) + x(23)*x(35)));
    pVec(5) = (1/(Nm*Nf))*((1/4)*(x(28)*x(23) + x(19)*x(32)) + (1/2)*(x(28)*x(26) + x(19)*x(35)) + ...
        (1/2)*(x(28)*x(24) + x(19)*x(33)) + (x(28)*x(27) + x(19)*x(36)) + ...
        (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/4)*(x(31)*x(23) + x(22)*x(32)) + ...
        (1/2)*(x(31)*x(21) + x(22)*x(30)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + ...
        (1/2)*(x(31)*x(24) + x(22)*x(33)) + (1/2)*(x(31)*x(27) + x(22)*x(36)) + ...
        (1/4)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(25) + x(20)*x(34)) + ...
        (1/2)*(x(29)*x(26) + x(20)*x(35)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + ...
        (1/2)*(x(29)*x(27) + x(20)*x(36)) + (1/4)*x(23)*x(32) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + ...
        (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/4)*(x(32)*x(26) + x(23)*x(35)) + ...
        (1/4)*(x(32)*x(24) + x(23)*x(32)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + ...
        (x(34)*x(21) + x(25)*x(30)) + (1/2)*(x(34)*x(24) + x(25)*x(33)) + ...
        (1/2)*(x(30)*x(26) + x(21)*x(35)) + (1/4)*(x(35)*x(24) + x(26)*x(33)));
    pVec(6) = (1/(Nm*Nf))*((1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/4)*(x(29)*x(26) + x(20)*x(35)) + ...
        (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/2)*(x(29)*x(27) + x(20)*x(36)) + ...
        (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + ...
        (1/4)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + ...
        (1/2)*(x(30)*x(26) + x(21)*x(35)) + (1/2)*(x(30)*x(24) + x(21)*x(33)) + (x(30)*x(27) + x(21)*x(36)) + ...
        (1/2)*x(33)*x(24) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + ...
        (1/2)*(x(33)*x(27) + x(24)*x(36)));
    pVec(7) = (1/(Nm*Nf)) *((1/4)*x(22)*x(31) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(25) + x(22)*x(34)) + ...
        (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/16)*x(23)*x(32) + ...
        (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + x(34)*x(25) + ...
        (1/2)*(x(34)*x(26) + x(25)*x(35)) + (1/4)*x(35)*x(26));
    pVec(8) = (1/(Nm*Nf))*((1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + ...
        (1/4)*(x(31)*x(24) + x(22)*x(33)) + (1/2)*(x(31)*x(27) + x(22)*x(36)) + ...
        (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/4)*(x(32)*x(26) + x(23)*x(35)) + ...
        (1/8)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + ...
        (1/2)*(x(34)*x(26) + x(25)*x(35)) + (1/2)*(x(34)*x(24) + x(25)*x(33)) + (x(34)*x(27) + x(25)*x(36)) + ...
        (1/2)*x(35)*x(26) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + ...
        (1/2)*(x(35)*x(27) + x(26)*x(36)));
    pVec(9) = (1/(Nm*Nf))*((1/16)*x(32)*x(23) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + ...
        (1/4)*(x(32)*x(27) + x(23)*x(36)) + (1/4)*x(35)*x(26) + ...
        (1/4)*(x(35)*x(24) + x(26)*x(33)) + (1/2)*(x(35)*x(27) + x(26)*x(36)) + (1/4)*x(33)*x(24) + ...
        (1/2)*(x(33)*x(27) + x(24)*x(36)) + x(36)*x(27));
    
    % point check---probabilities sum to 1
    if abs(sum(pVec) - 1) > 10^(-9)
        error("Error: Offspring genotype probabilities do not sum to 1!")
    end
    % point check---probabilities are positive
    if any(pVec < 0)
        error("Error: Offspring genotype probabilities are negative!")
    end
    
    %%%%%%%%%%%%% coupled genetic/population dynamics model %%%%%%%%%%%%% 
    if (lethality_type == "EA")
        if (lethality_case == "BSL") % bi-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);               % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(3);       % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(5);   % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(7);       % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(9);       % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3)*c_t - (muX+m+alpha*Nl^(beta+1))*x(11); % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(13); % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3)*c_t - (muX+m+alpha*Nl^(beta+1))*x(15); % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(17);     % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(2);               % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(4);       % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(6);   % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(8);       % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(10);     % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3)*c_t - (muX+m+alpha*Nl^(beta+1))*x(12); % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(14); % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3)*c_t - (muX+m+alpha*Nl^(beta+1))*x(16); % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(18);     % aabb
            else % strongly suppressed 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);               % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(3);       % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(5);   % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(7);       % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(9);       % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(11);     % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(13); % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(15);     % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(17);     % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(2);               % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(4);       % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(6);   % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(8);       % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(10);     % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(12);     % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(14); % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(16);     % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(18);     % aabb
            end % end of if supp statement
        else % female-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);               % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a - (muX+m+alpha*Nl^(beta+1))*x(3);           % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a)^2 - (muX+m+alpha*Nl^(beta+1))*x(5);       % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a - (muX+m+alpha*Nl^(beta+1))*x(7);           % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(9);       % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(11);     % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a)^2 - (muX+m+alpha*Nl^(beta+1))*x(13);     % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(15);     % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(17);     % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(2);               % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(4);       % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(6);   % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(8);       % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(10);     % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3)*c_t - (muX+m+alpha*Nl^(beta+1))*x(12); % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(14); % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3)*c_t - (muX+m+alpha*Nl^(beta+1))*x(16); % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(18);     % aabb
            else % strongly suppressed 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);               % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a - (muX+m+alpha*Nl^(beta+1))*x(3);           % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a)^2 - (muX+m+alpha*Nl^(beta+1))*x(5);       % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a - (muX+m+alpha*Nl^(beta+1))*x(7);           % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(9);       % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(11);     % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a)^2 - (muX+m+alpha*Nl^(beta+1))*x(13);     % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(15);     % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(17);     % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(2);               % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(4);       % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(6);   % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+alpha*Nl^(beta+1))*x(8);       % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+alpha*Nl^(beta+1))*x(10);     % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(12);     % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+alpha*Nl^(beta+1))*x(14); % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+alpha*Nl^(beta+1))*x(16);     % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+alpha*Nl^(beta+1))*x(18);     % aabb
            end % end of if supp statement
        end
        
        % adults for each genotype
        % males
        dx(19) = m*x(1) - muZ*x(19);    % AABB male
        dx(20) = m*x(3) - muZ*x(20);    % AABb male
        dx(21) = m*x(5) - muZ*x(21);    % AAbb male
        dx(22) = m*x(7) - muZ*x(22);    % AaBB male
        dx(23) = m*x(9) - muZ*x(23);    % AaBb male 
        dx(24) = m*x(11) - muZ*x(24);   % Aabb male
        dx(25) = m*x(13) - muZ*x(25);   % aaBB male
        dx(26) = m*x(15) - muZ*x(26);   % aaBb male
        dx(27) = m*x(17) - muZ*x(27);   % aabb male
        % females
        dx(28) = m*x(2) - muY*x(28);    % AABB female
        dx(29) = m*x(4) - muY*x(29);    % AABb female
        dx(30) = m*x(6) - muY*x(30);    % AAbb female
        dx(31) = m*x(8) - muY*x(31);    % AaBB female
        dx(32) = m*x(10) - muY*x(32);   % AaBb female 
        dx(33) = m*x(12) - muY*x(33);   % Aabb female
        dx(34) = m*x(14) - muY*x(34);   % aaBB female
        dx(35) = m*x(16) - muY*x(35);   % aaBb female
        dx(36) = m*x(18) - muY*x(36);   % aabb female
    else % fitness costs are late-acting
        if (lethality_case == "BSL") % bi-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a*c_t - muZ*x(20);        % AABb male
                dx(21) = m*x(5)*(c_a*c_t)^2 - muZ*x(21);    % AAbb male
                dx(22) = m*x(7)*c_a*c_t - muZ*x(22);        % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3)*c_t - muZ*x(24);   % Aabb male
                dx(25) = m*x(13)*(c_a*c_t)^2 - muZ*x(25);   % aaBB male
                dx(26) = m*x(15)*(c_a^3)*c_t - muZ*x(26);   % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3)*c_t - muY*x(33);   % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3)*c_t - muY*x(35);   % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            else % strongly suppressed 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a*c_t - muZ*x(20);        % AABb male
                dx(21) = m*x(5)*(c_a*c_t)^2 - muZ*x(21);    % AAbb male
                dx(22) = m*x(7)*c_a*c_t - muZ*x(22);        % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3) - muZ*x(24);       % Aabb male
                dx(25) = m*x(13)*(c_a*c_t)^2 - muZ*x(25);   % aaBB male
                dx(26) = m*x(15)*(c_a^3) - muZ*x(26);       % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3) - muY*x(33);       % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3) - muY*x(35);       % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            end % end of if supp statement
        else % female-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a - muZ*x(20);            % AABb male
                dx(21) = m*x(5)*(c_a)^2 - muZ*x(21);        % AAbb male
                dx(22) = m*x(7)*c_a - muZ*x(22);            % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3) - muZ*x(24);       % Aabb male
                dx(25) = m*x(13)*(c_a)^2 - muZ*x(25);       % aaBB male
                dx(26) = m*x(15)*(c_a^3) - muZ*x(26);       % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3)*c_t - muY*x(33);   % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3)*c_t - muY*x(35);   % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            else % strongly suppressed 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a - muZ*x(20);            % AABb male
                dx(21) = m*x(5)*(c_a)^2 - muZ*x(21);        % AAbb male
                dx(22) = m*x(7)*c_a - muZ*x(22);            % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3) - muZ*x(24);       % Aabb male
                dx(25) = m*x(13)*(c_a)^2 - muZ*x(25);       % aaBB male
                dx(26) = m*x(15)*(c_a^3) - muZ*x(26);       % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3) - muY*x(33);       % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3) - muY*x(35);       % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            end % end of if supp statement
        end

        % larvae for each genotype
        % males
        dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(1);       % AABB 
        dx(3) = (lambda/2)*(Nf)*pVec(2) - (muX+m+alpha*Nl^(beta+1))*x(3);       % AABb
        dx(5) = (lambda/2)*(Nf)*pVec(3) - (muX+m+alpha*Nl^(beta+1))*x(5);       % AAbb
        dx(7) = (lambda/2)*(Nf)*pVec(4) - (muX+m+alpha*Nl^(beta+1))*x(7);       % AaBB
        dx(9) = (lambda/2)*(Nf)*pVec(5) - (muX+m+alpha*Nl^(beta+1))*x(9);       % AaBb
        dx(11) = (lambda/2)*(Nf)*pVec(6) - (muX+m+alpha*Nl^(beta+1))*x(11);     % Aabb
        dx(13) = (lambda/2)*(Nf)*pVec(7) - (muX+m+alpha*Nl^(beta+1))*x(13);     % aaBB
        dx(15) = (lambda/2)*(Nf)*pVec(8) - (muX+m+alpha*Nl^(beta+1))*x(15);     % aaBb
        dx(17) = (lambda/2)*(Nf)*pVec(9) - (muX+m+alpha*Nl^(beta+1))*x(17);     % aabb
        % females
        dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+alpha*Nl^(beta+1))*x(2);       % AABB 
        dx(4) = (lambda/2)*(Nf)*pVec(2) - (muX+m+alpha*Nl^(beta+1))*x(4);       % AABb
        dx(6) = (lambda/2)*(Nf)*pVec(3) - (muX+m+alpha*Nl^(beta+1))*x(6);       % AAbb
        dx(8) = (lambda/2)*(Nf)*pVec(4) - (muX+m+alpha*Nl^(beta+1))*x(8);       % AaBB
        dx(10) = (lambda/2)*(Nf)*pVec(5) - (muX+m+alpha*Nl^(beta+1))*x(10);     % AaBb
        dx(12) = (lambda/2)*(Nf)*pVec(6) - (muX+m+alpha*Nl^(beta+1))*x(12);     % Aabb
        dx(14) = (lambda/2)*(Nf)*pVec(7) - (muX+m+alpha*Nl^(beta+1))*x(14);     % aaBB
        dx(16) = (lambda/2)*(Nf)*pVec(8) - (muX+m+alpha*Nl^(beta+1))*x(16);     % aaBb
        dx(18) = (lambda/2)*(Nf)*pVec(9) - (muX+m+alpha*Nl^(beta+1))*x(18);     % aabb
    end % end of if lethality_type statement


    
    %%%%%%%%%%%%% disease model %%%%%%%%%%%%% 
    dx(37) = a*c*(x(38)/M)*(Nf - x(37)) - muY*x(37);    % infected adult mosquitoes
    dx(38) = a*b*(x(37)/M)*(M - x(38)) - gamma*x(38);   % infected humans 
    
    dx = dx';

end


% gene drive, version 3 density-dependence 
function dx = DDGD_iii(t,x,params,alpha,beta,lethality_case,supp,lethality_type)
    % Case study for two-locus engineered underdominance (male-sex
    % lethality is ignored)
    % ---------------------------------------------------------------------

    % [AABB_larvae_male, AABB_larvae_female 2
    %  AABb_larvae_male, AABb_larvae_female 4
    %  AAbb_larvae_male, AAbb_larvae_female 6
    %  AaBB_larvae_male, AaBB_larvae_female 8
    %  AaBb_larvae_male, AaBb_larvae_female 10
    %  Aabb_larvae_male, Aabb_larvae_female 12
    %  aaBB_larvae_male, aaBB_larvae_female 14
    %  aaBb_larvae_male, aaBb_larvae_female 16
    %  aabb_larvae_male, aabb_larvae_female 18
    %  AABB_adult_males, AABb_adult_males 20
    %  AAbb_adult_males, AaBB_adult_males 22
    %  AaBb_adult_males, Aabb_adult_males 24
    %  aaBB_adult_males, aaBb_adult_males 26
    %  aabb_adult_males, AABB_adult_females 28
    %  AABb_adult_females, AAbb_adult_females 30
    %  AaBB_adult_females, AaBb_adult_females 32
    %  Aabb_adult_females, aaBB_adult_females 34
    %  aaBb_adult_females, aabb_adult_females 36
    %  infected_mosqs, infected_humans] 38
    
    lambda      = params.lambda;
    muX         = params.muX;
    muZ         = params.muZ;
    muY         = params.muY;
    m           = params.m;
    a           = params.a;
    b           = params.b;
    c           = params.c;
    gamma       = params.gamma;
    M           = params.M;
    s_a         = params.s_a;
    s_t         = params.s_t;
    % N = params(13); 
    
    % impose hard cut off at zero
    % x(x < 0) = 0; 
    
    % compute relative fitnesses for ambient carriage and transgene
    c_a = 1 - s_a;
    c_t = 1 - s_t;
    
    % memory allocation
    dx=zeros(1,38);
    
    Nl = sum(x(1:18)); % mosquito larvae pop
    Nm = sum(x(19:27)); % mosquito male pop
    Nf = sum(x(28:36)); % mosquito female pop
    
    % calculate probabilities of offspring having each genotype
    pVec = zeros(1,9);
    pVec(1) = (1/(Nm*Nf))*(x(19)*x(28)+(1/2)*(x(28)*x(22)+x(19)*x(31)) + ...
        (1/2)*(x(28)*x(20)+x(19)*x(29)) + 0.25*(x(28)*x(23)+x(19)*x(32)) + ...
        0.25*x(22)*x(31) + 0.25*(x(20)*x(31) + x(22)*x(29)) + ...
        (1/8)*(x(23)*x(31) + x(32)*x(22)) + (1/4)*x(20)*x(29) + ...
        (1/8)*(x(23)*x(29) + x(32)*x(20)) + (1/16)*x(23)*x(32));
    pVec(2) = (1/(Nm*Nf))*((1/2)*(x(28)*x(20) + x(29)*x(19)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + ...
        (x(28)*x(21) + x(19)*x(30)) + (1/2)*(x(28)*x(24) + x(19)*x(33)) + ...
        (1/2)*x(20)*x(29) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + ...
        (1/2)*(x(31)*x(21) + x(22)*x(30)) + (1/4)*(x(31)*x(24) + x(22)*x(32)) + ...
        (1/4)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(21) + x(20)*x(30)) + ...
        (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/8)*x(23)*x(32) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + ...
        (1/8)*(x(32)*x(24)+x(23)*x(32))); 
    pVec(3) = (1/(Nm*Nf))*((1/4)*x(29)*x(20) + (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(21)+x(20)*x(30)) + ...
        (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/16)*x(32)*x(23) + ...
        (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + x(30)*x(21) + ...
        (1/2)*(x(30)*x(24) + x(21)*x(33)) + (1/4)*x(33)*x(24));
    pVec(4) = (1/(Nm*Nf))*((1/2)*(x(28)*x(22) + x(19)*x(31)) + (1/4)*(x(28)*x(23) + x(19)*x(32)) + ...
        (x(28)*x(25) + x(19)*x(34)) + (1/2)*(x(28)*x(26) + x(19)*x(35)) + ...
        (1/2)*x(31)*x(22) + (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/4)*(x(31)*x(23) + x(22)*x(32)) + ...
        (1/2)*(x(31)*x(25) + x(22)*x(34)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + ...
        (1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(25) + x(20)*x(34)) + ...
        (1/4)*(x(29)*x(26) + x(20)*x(35)) + (1/8)*x(23)*x(32) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + ...
        (1/8)*(x(32)*x(26) + x(23)*x(35)));
    pVec(5) = (1/(Nm*Nf))*((1/4)*(x(28)*x(23) + x(19)*x(32)) + (1/2)*(x(28)*x(26) + x(19)*x(35)) + ...
        (1/2)*(x(28)*x(24) + x(19)*x(33)) + (x(28)*x(27) + x(19)*x(36)) + ...
        (1/4)*(x(31)*x(20) + x(22)*x(29)) + (1/4)*(x(31)*x(23) + x(22)*x(32)) + ...
        (1/2)*(x(31)*x(21) + x(22)*x(30)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + ...
        (1/2)*(x(31)*x(24) + x(22)*x(33)) + (1/2)*(x(31)*x(27) + x(22)*x(36)) + ...
        (1/4)*(x(29)*x(23) + x(20)*x(32)) + (1/2)*(x(29)*x(25) + x(20)*x(34)) + ...
        (1/2)*(x(29)*x(26) + x(20)*x(35)) + (1/4)*(x(29)*x(24) + x(20)*x(33)) + ...
        (1/2)*(x(29)*x(27) + x(20)*x(36)) + (1/4)*x(23)*x(32) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + ...
        (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/4)*(x(32)*x(26) + x(23)*x(35)) + ...
        (1/4)*(x(32)*x(24) + x(23)*x(32)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + ...
        (x(34)*x(21) + x(25)*x(30)) + (1/2)*(x(34)*x(24) + x(25)*x(33)) + ...
        (1/2)*(x(30)*x(26) + x(21)*x(35)) + (1/4)*(x(35)*x(24) + x(26)*x(33)));
    pVec(6) = (1/(Nm*Nf))*((1/8)*(x(29)*x(23) + x(20)*x(32)) + (1/4)*(x(29)*x(26) + x(20)*x(35)) + ...
        (1/4)*(x(29)*x(24) + x(20)*x(33)) + (1/2)*(x(29)*x(27) + x(20)*x(36)) + ...
        (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(21) + x(23)*x(30)) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + ...
        (1/4)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + ...
        (1/2)*(x(30)*x(26) + x(21)*x(35)) + (1/2)*(x(30)*x(24) + x(21)*x(33)) + (x(30)*x(27) + x(21)*x(36)) + ...
        (1/2)*x(33)*x(24) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + ...
        (1/2)*(x(33)*x(27) + x(24)*x(36)));
    pVec(7) = (1/(Nm*Nf)) *((1/4)*x(22)*x(31) + (1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/2)*(x(31)*x(25) + x(22)*x(34)) + ...
        (1/4)*(x(31)*x(26) + x(22)*x(35)) + (1/16)*x(23)*x(32) + ...
        (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + x(34)*x(25) + ...
        (1/2)*(x(34)*x(26) + x(25)*x(35)) + (1/4)*x(35)*x(26));
    pVec(8) = (1/(Nm*Nf))*((1/8)*(x(31)*x(23) + x(22)*x(32)) + (1/4)*(x(31)*x(26) + x(22)*x(35)) + ...
        (1/4)*(x(31)*x(24) + x(22)*x(33)) + (1/2)*(x(31)*x(27) + x(22)*x(36)) + ...
        (1/8)*x(32)*x(23) + (1/4)*(x(32)*x(25) + x(23)*x(34)) + (1/4)*(x(32)*x(26) + x(23)*x(35)) + ...
        (1/8)*(x(32)*x(24) + x(23)*x(33)) + (1/4)*(x(32)*x(27) + x(23)*x(36)) + ...
        (1/2)*(x(34)*x(26) + x(25)*x(35)) + (1/2)*(x(34)*x(24) + x(25)*x(33)) + (x(34)*x(27) + x(25)*x(36)) + ...
        (1/2)*x(35)*x(26) + (1/4)*(x(35)*x(24) + x(26)*x(33)) + ...
        (1/2)*(x(35)*x(27) + x(26)*x(36)));
    pVec(9) = (1/(Nm*Nf))*((1/16)*x(32)*x(23) + (1/8)*(x(32)*x(26) + x(23)*x(35)) + (1/8)*(x(32)*x(24) + x(23)*x(33)) + ...
        (1/4)*(x(32)*x(27) + x(23)*x(36)) + (1/4)*x(35)*x(26) + ...
        (1/4)*(x(35)*x(24) + x(26)*x(33)) + (1/2)*(x(35)*x(27) + x(26)*x(36)) + (1/4)*x(33)*x(24) + ...
        (1/2)*(x(33)*x(27) + x(24)*x(36)) + x(36)*x(27));
    
    % point check
    if abs(sum(pVec) - 1) > 10^(-9)
        error("Error: Offspring genotype probabilities do not sum to 1!")
    end
    % point check---probabilities are positive
    if any(pVec < 0)
        error("Error: Offspring genotype probabilities are negative!")
    end    
    
    % [AABB_larvae_male, AABB_larvae_female 2
    %  AABb_larvae_male, AABb_larvae_female 4
    %  AAbb_larvae_male, AAbb_larvae_female 6
    %  AaBB_larvae_male, AaBB_larvae_female 8
    %  AaBb_larvae_male, AaBb_larvae_female 10
    %  Aabb_larvae_male, Aabb_larvae_female 12
    %  aaBB_larvae_male, aaBB_larvae_female 14
    %  aaBb_larvae_male, aaBb_larvae_female 16
    %  aabb_larvae_male, aabb_larvae_female 18
    %  AABB_adult_males, AABb_adult_males 20
    %  AAbb_adult_males, AaBB_adult_males 22
    %  AaBb_adult_males, Aabb_adult_males 24
    %  aaBB_adult_males, aaBb_adult_males 26
    %  aabb_adult_males, AABB_adult_females 28
    %  AABb_adult_females, AAbb_adult_females 30
    %  AaBB_adult_females, AaBb_adult_females 32
    %  Aabb_adult_females, aaBB_adult_females 34
    %  aaBb_adult_females, aabb_adult_females 36
    %  infected_mosqs, infected_humans] 38
    %%%%%%%%%%%%% coupled genetic/population dynamics model %%%%%%%%%%%%% 
    if (lethality_type == "EA")
        if (lethality_case == "BSL") % bi-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1); % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(3);          % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(5);      % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(7);          % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(9);          % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3)*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(11);    % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(13);    % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3)*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(15);    % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(17);        % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);                  % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(4);          % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(6);      % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(8);          % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(10);        % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3)*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(12);    % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(14);    % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3)*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(16);    % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(18);        % aabb
            else % strongly suppressed 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);                  % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(3);          % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(5);      % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(7);          % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(9);          % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(11);        % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(13);    % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(15);        % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(17);        % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);                  % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(4);          % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(6);      % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(8);          % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(10);        % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(12);        % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(14);    % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(16);        % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(18);        % aabb
            end % end of if supp statement
        else % female-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);                  % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a - (muX+m+log(1+(alpha*Nl)^beta))*x(3);              % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(5);          % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a - (muX+m+log(1+(alpha*Nl)^beta))*x(7);              % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(9);          % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(11);        % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(13);        % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(15);        % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(17);        % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);                  % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(4);          % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(6);      % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(8);          % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(10);        % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3)*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(12);    % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(14);    % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3)*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(16);    % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(18);        % aabb
            else % strongly suppressed 
                % larvae for each genotype
                % males
                dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);                  % AABB 
                dx(3) = (lambda/2)*(Nf)*pVec(2)*c_a - (muX+m+log(1+(alpha*Nl)^beta))*x(3);              % AABb
                dx(5) = (lambda/2)*(Nf)*pVec(3)*(c_a)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(5);          % AAbb
                dx(7) = (lambda/2)*(Nf)*pVec(4)*c_a - (muX+m+log(1+(alpha*Nl)^beta))*x(7);              % AaBB
                dx(9) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(9);          % AaBb
                dx(11) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(11);        % Aabb
                dx(13) = (lambda/2)*(Nf)*pVec(7)*(c_a)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(13);        % aaBB
                dx(15) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(15);        % aaBb
                dx(17) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(17);        % aabb
                % females
                dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);                  % AABB 
                dx(4) = (lambda/2)*(Nf)*pVec(2)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(4);          % AABb
                dx(6) = (lambda/2)*(Nf)*pVec(3)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(6);      % AAbb
                dx(8) = (lambda/2)*(Nf)*pVec(4)*c_a*c_t - (muX+m+log(1+(alpha*Nl)^beta))*x(8);          % AaBB
                dx(10) = (lambda/2)*(Nf)*pVec(5)*(c_a^2) - (muX+m+log(1+(alpha*Nl)^beta))*x(10);        % AaBb
                dx(12) = (lambda/2)*(Nf)*pVec(6)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(12);        % Aabb
                dx(14) = (lambda/2)*(Nf)*pVec(7)*(c_a*c_t)^2 - (muX+m+log(1+(alpha*Nl)^beta))*x(14);    % aaBB
                dx(16) = (lambda/2)*(Nf)*pVec(8)*(c_a^3) - (muX+m+log(1+(alpha*Nl)^beta))*x(16);        % aaBb
                dx(18) = (lambda/2)*(Nf)*pVec(9)*(c_a^4) - (muX+m+log(1+(alpha*Nl)^beta))*x(18);        % aabb
            end % end of if supp statement
        end
        
        % adults for each genotype
        % males
        dx(19) = m*x(1) - muZ*x(19);    % AABB male
        dx(20) = m*x(3) - muZ*x(20);    % AABb male
        dx(21) = m*x(5) - muZ*x(21);    % AAbb male
        dx(22) = m*x(7) - muZ*x(22);    % AaBB male
        dx(23) = m*x(9) - muZ*x(23);    % AaBb male 
        dx(24) = m*x(11) - muZ*x(24);   % Aabb male
        dx(25) = m*x(13) - muZ*x(25);   % aaBB male
        dx(26) = m*x(15) - muZ*x(26);   % aaBb male
        dx(27) = m*x(17) - muZ*x(27);   % aabb male
        % females
        dx(28) = m*x(2) - muY*x(28);    % AABB female
        dx(29) = m*x(4) - muY*x(29);    % AABb female
        dx(30) = m*x(6) - muY*x(30);    % AAbb female
        dx(31) = m*x(8) - muY*x(31);    % AaBB female
        dx(32) = m*x(10) - muY*x(32);   % AaBb female 
        dx(33) = m*x(12) - muY*x(33);   % Aabb female
        dx(34) = m*x(14) - muY*x(34);   % aaBB female
        dx(35) = m*x(16) - muY*x(35);   % aaBb female
        dx(36) = m*x(18) - muY*x(36);   % aabb female
    else % late-acting fitness costs
        if (lethality_case == "BSL") % bi-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a*c_t - muZ*x(20);        % AABb male
                dx(21) = m*x(5)*(c_a*c_t)^2 - muZ*x(21);    % AAbb male
                dx(22) = m*x(7)*c_a*c_t - muZ*x(22);        % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3)*c_t - muZ*x(24);   % Aabb male
                dx(25) = m*x(13)*(c_a*c_t)^2 - muZ*x(25);   % aaBB male
                dx(26) = m*x(15)*(c_a^3)*c_t - muZ*x(26);   % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3)*c_t - muY*x(33);   % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3)*c_t - muY*x(35);   % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            else % strongly suppressed 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a*c_t - muZ*x(20);        % AABb male
                dx(21) = m*x(5)*(c_a*c_t)^2 - muZ*x(21);    % AAbb male
                dx(22) = m*x(7)*c_a*c_t - muZ*x(22);        % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3) - muZ*x(24);       % Aabb male
                dx(25) = m*x(13)*(c_a*c_t)^2 - muZ*x(25);   % aaBB male
                dx(26) = m*x(15)*(c_a^3) - muZ*x(26);       % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3) - muY*x(33);       % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3) - muY*x(35);       % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            end % end of if supp statement
        else % female-sex lethality
            % are lethal alleles weakly or strongly suppressed?
            if (supp == "WS") 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a - muZ*x(20);            % AABb male
                dx(21) = m*x(5)*(c_a)^2 - muZ*x(21);        % AAbb male
                dx(22) = m*x(7)*c_a - muZ*x(22);            % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3) - muZ*x(24);       % Aabb male
                dx(25) = m*x(13)*(c_a)^2 - muZ*x(25);       % aaBB male
                dx(26) = m*x(15)*(c_a^3) - muZ*x(26);       % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3)*c_t - muY*x(33);   % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3)*c_t - muY*x(35);   % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            else % strongly suppressed 
                % adults for each genotype
                % males
                dx(19) = m*x(1) - muZ*x(19);                % AABB male
                dx(20) = m*x(3)*c_a - muZ*x(20);            % AABb male
                dx(21) = m*x(5)*(c_a)^2 - muZ*x(21);        % AAbb male
                dx(22) = m*x(7)*c_a - muZ*x(22);            % AaBB male
                dx(23) = m*x(9)*(c_a^2) - muZ*x(23);        % AaBb male 
                dx(24) = m*x(11)*(c_a^3) - muZ*x(24);       % Aabb male
                dx(25) = m*x(13)*(c_a)^2 - muZ*x(25);       % aaBB male
                dx(26) = m*x(15)*(c_a^3) - muZ*x(26);       % aaBb male
                dx(27) = m*x(17)*(c_a^4) - muZ*x(27);       % aabb male
                % females
                dx(28) = m*x(2) - muY*x(28);                % AABB female
                dx(29) = m*x(4)*c_a*c_t - muY*x(29);        % AABb female
                dx(30) = m*x(6)*(c_a*c_t)^2 - muY*x(30);    % AAbb female
                dx(31) = m*x(8)*c_a*c_t - muY*x(31);        % AaBB female
                dx(32) = m*x(10)*(c_a^2) - muY*x(32);       % AaBb female 
                dx(33) = m*x(12)*(c_a^3) - muY*x(33);       % Aabb female
                dx(34) = m*x(14)*(c_a*c_t)^2 - muY*x(34);   % aaBB female
                dx(35) = m*x(16)*(c_a^3) - muY*x(35);       % aaBb female
                dx(36) = m*x(18)*(c_a^4) - muY*x(36);       % aabb female
            end % end of if supp statement
        end
        % larvae for each genotype
        % males
        dx(1) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(1);      % AABB 
        dx(3) = (lambda/2)*(Nf)*pVec(2) - (muX+m+log(1+(alpha*Nl)^beta))*x(3);      % AABb
        dx(5) = (lambda/2)*(Nf)*pVec(3) - (muX+m+log(1+(alpha*Nl)^beta))*x(5);      % AAbb
        dx(7) = (lambda/2)*(Nf)*pVec(4) - (muX+m+log(1+(alpha*Nl)^beta))*x(7);      % AaBB
        dx(9) = (lambda/2)*(Nf)*pVec(5) - (muX+m+log(1+(alpha*Nl)^beta))*x(9);      % AaBb
        dx(11) = (lambda/2)*(Nf)*pVec(6) - (muX+m+log(1+(alpha*Nl)^beta))*x(11);    % Aabb
        dx(13) = (lambda/2)*(Nf)*pVec(7) - (muX+m+log(1+(alpha*Nl)^beta))*x(13);    % aaBB
        dx(15) = (lambda/2)*(Nf)*pVec(8) - (muX+m+log(1+(alpha*Nl)^beta))*x(15);    % aaBb
        dx(17) = (lambda/2)*(Nf)*pVec(9) - (muX+m+log(1+(alpha*Nl)^beta))*x(17);    % aabb
        % females
        dx(2) = (lambda/2)*(Nf)*pVec(1) - (muX+m+log(1+(alpha*Nl)^beta))*x(2);      % AABB 
        dx(4) = (lambda/2)*(Nf)*pVec(2) - (muX+m+log(1+(alpha*Nl)^beta))*x(4);      % AABb
        dx(6) = (lambda/2)*(Nf)*pVec(3) - (muX+m+log(1+(alpha*Nl)^beta))*x(6);      % AAbb
        dx(8) = (lambda/2)*(Nf)*pVec(4) - (muX+m+log(1+(alpha*Nl)^beta))*x(8);      % AaBB
        dx(10) = (lambda/2)*(Nf)*pVec(5) - (muX+m+log(1+(alpha*Nl)^beta))*x(10);    % AaBb
        dx(12) = (lambda/2)*(Nf)*pVec(6) - (muX+m+log(1+(alpha*Nl)^beta))*x(12);    % Aabb
        dx(14) = (lambda/2)*(Nf)*pVec(7) - (muX+m+log(1+(alpha*Nl)^beta))*x(14);    % aaBB
        dx(16) = (lambda/2)*(Nf)*pVec(8) - (muX+m+log(1+(alpha*Nl)^beta))*x(16);    % aaBb
        dx(18) = (lambda/2)*(Nf)*pVec(9) - (muX+m+log(1+(alpha*Nl)^beta))*x(18);    % aabb        
    end % end of if lethality_type statement
    
    %%%%%%%%%%%%% disease model %%%%%%%%%%%%% 
    dx(37) = a*c*(x(38)/M)*(Nf - x(37)) - muY*x(37); % infected adult mosquitoes
    dx(38) = a*b*(x(37)/M)*(M - x(38)) - gamma*x(38); % infected humans 
    
    % end
    
    dx = dx';

end

% general genetic load functions
function dx = vanilla_DDGD_ii(t,x,params)
    % General genetic load ODE system
    % ---------------------------------------------------------------------
    % [female_juveniles, male_juveniles, ...
    % [female_adults, male_adults]
    
    % function simulating the generalized logistic case---note that the
    % function is programmed for beta+1 and not beta
    lambda          = params.lambda;
    muX             = params.muX;
    muZ             = params.muZ;
    muY             = params.muY;
    m               = params.m;
    g               = params.g;
    % other params
    beta            = params.beta;
    lethality_case  = params.lethality_case;
    lethality_type  = params.lethality_type;

    switch beta
    case 0.5
        alpha = 0.0321302;
    case 1
        alpha = 1.90085*10^(-4);
    case 1.5
        alpha = 1.12456*10^(-6);
    otherwise
        error('Error! Beta must be 0.5, 1, or 1.5.')
    end
    
    % memory allocation
    dx=zeros(1,4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change beta so that things make sense
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta = beta - 1; 
    
    J = sum(x(1:2)); % mosquito larvae pop
    
    %%%%%%%%%%%%% coupled genetic/population dynamics model %%%%%%%%%%%%% 
    if (lethality_type == "EA")
        if (lethality_case == "BSL") % bi-sex lethality
            dx(1) = (lambda/2)*x(3)*(1-g) - (muX+m+alpha*J^(beta+1))*x(1);
            dx(2) = (lambda/2)*x(3)*(1-g) - (muX+m+alpha*J^(beta+1))*x(2);
        elseif (lethality_case == "FSL") % female-sex lethality
            dx(1) = (lambda/2)*x(3)*(1-g) - (muX+m+alpha*J^(beta+1))*x(1);
            dx(2) = (lambda/2)*x(3) - (muX+m+alpha*J^(beta+1))*x(2);
        elseif (lethality_case == "MSL") % male-sex lethality
            dx(1) = (lambda/2)*x(3) - (muX+m+alpha*J^(beta+1))*x(1);
            dx(2) = (lambda/2)*x(3)*(1-g) - (muX+m+alpha*J^(beta+1))*x(2);
        end

        % adults
        dx(3) = m*x(1) - muY*x(3); % female adults
        dx(4) = m*x(2) - muZ*x(4); % male adults

    else % fitness costs are late-acting
        if (lethality_case == "BSL") % bi-sex lethality
            dx(3) = m*x(1)*(1-g) - muY*x(3);
            dx(4) = m*x(2)*(1-g) - muZ*x(4);
        elseif (lethality_case == "FSL") % female-sex lethality
            dx(3) = m*x(1)*(1-g) - muY*x(3);
            dx(4) = m*x(2) - muZ*x(4);
        elseif (lethality_case == "MSL") % male-sex lethality
            dx(3) = m*x(1) - muY*x(3);
            dx(4) = m*x(2)*(1-g) - muZ*x(4);
        end

        % juveniles
        dx(1) = (lambda/2)*x(3) - (muX+m+alpha*J^(beta+1))*x(1);
        dx(2) = (lambda/2)*x(3) - (muX+m+alpha*J^(beta+1))*x(2);

    end % end of if lethality_type statement
    
    dx = dx';

end

function dx = vanilla_DDGD_iii(t,x,params)
    % General genetic load ODE system
    % ---------------------------------------------------------------------
    % [female_juveniles, male_juveniles, ...
    % [female_adults, male_adults]
    
    % function simulating the generalized logistic case---note that the
    % function is programmed for beta+1 and not beta
    lambda          = params.lambda;
    muX             = params.muX;
    muZ             = params.muZ;
    muY             = params.muY;
    m               = params.m;
    g               = params.g;
    % other params
    % alpha           = params.alpha;
    alpha           = 0.014543805956894;
    beta            = params.beta;
    lethality_case  = params.lethality_case;
    lethality_type  = params.lethality_type;
    
    % memory allocation
    dx=zeros(1,4);
    
    J = sum(x(1:2)); % mosquito larvae pop
    
    %%%%%%%%%%%%% coupled genetic/population dynamics model %%%%%%%%%%%%% 
    if (lethality_type == "EA")
        if (lethality_case == "BSL") % bi-sex lethality
            dx(1) = (lambda/2)*x(3)*(1-g) - (muX+m+log(1+(alpha*J)^beta))*x(1);
            dx(2) = (lambda/2)*x(3)*(1-g) - (muX+m+log(1+(alpha*J)^beta))*x(2);
        elseif (lethality_case == "FSL") % female-sex lethality
            dx(1) = (lambda/2)*x(3)*(1-g) - (muX+m+log(1+(alpha*J)^beta))*x(1);
            dx(2) = (lambda/2)*x(3) - (muX+m+log(1+(alpha*J)^beta))*x(2);
        elseif (lethality_case == "MSL") % male-sex lethality
            dx(1) = (lambda/2)*x(3) - (muX+m+log(1+(alpha*J)^beta))*x(1);
            dx(2) = (lambda/2)*x(3)*(1-g) - (muX+m+log(1+(alpha*J)^beta))*x(2);
        end

        % adults
        dx(3) = m*x(1) - muY*x(3); % female adults
        dx(4) = m*x(2) - muZ*x(4); % male adults

    else % fitness costs are late-acting
        if (lethality_case == "BSL") % bi-sex lethality
            dx(3) = m*x(1)*(1-g) - muY*x(3);
            dx(4) = m*x(2)*(1-g) - muZ*x(4);
        elseif (lethality_case == "FSL") % female-sex lethality
            dx(3) = m*x(1)*(1-g) - muY*x(3);
            dx(4) = m*x(2) - muZ*x(4);
        elseif (lethality_case == "MSL") % male-sex lethality
            dx(3) = m*x(1) - muY*x(3);
            dx(4) = m*x(2)*(1-g) - muZ*x(4);
        end

        % juveniles
        dx(1) = (lambda/2)*x(3) - (muX+m+log(1+(alpha*J)^beta))*x(1);
        dx(2) = (lambda/2)*x(3) - (muX+m+log(1+(alpha*J)^beta))*x(2);

    end % end of if lethality_type statement
    
    dx = dx';
end
