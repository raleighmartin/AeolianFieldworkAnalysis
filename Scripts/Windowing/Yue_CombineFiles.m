%combine windows for Yue
clearvars;
theta_bias_Oceano = 9.2428; %value of bias in theta for Oceano for adjusting values (based on mean of all times with fQ = 1)

%load first set of windows
load('/Users/raleighmartin/GitHub/AeolianFieldworkAnalysis/AnalysisData/Windowing/DataWindowCalcs_Oceano_Yue_1.mat');
StartTimes = StartTimes_all{1};
EndTimes = EndTimes_all{1};
Q = Q_all{1};
sigma_Q = sigma_Q_all{1};
tauRe = tauRe_all{1};
sigma_tauRe = sigma_tauRe_all{1};
ustRe = ustRe_all{1};
sigma_ustRe = sigma_ustRe_all{1};
zq = zq_all{1};
sigma_zq = sigma_zq_all{1};
theta = theta_all{1};
zL = zL_all{1};

%load and add second set of windows
load('/Users/raleighmartin/GitHub/AeolianFieldworkAnalysis/AnalysisData/Windowing/DataWindowCalcs_Oceano_Yue_2.mat');
StartTimes = [StartTimes; StartTimes_all{1}];
EndTimes = [EndTimes; EndTimes_all{1}];
Q = [Q; Q_all{1}];
sigma_Q = [sigma_Q; sigma_Q_all{1}];
tauRe = [tauRe; tauRe_all{1}];
sigma_tauRe = [sigma_tauRe; sigma_tauRe_all{1}];
ustRe = [ustRe; ustRe_all{1}];
sigma_ustRe = [sigma_ustRe; sigma_ustRe_all{1}];
zq = [zq; zq_all{1}];
sigma_zq = [sigma_zq; sigma_zq_all{1}];
theta = [theta; theta_all{1}];
zL = [zL; zL_all{1}];

%load and add third set of windows
load('/Users/raleighmartin/GitHub/AeolianFieldworkAnalysis/AnalysisData/Windowing/DataWindowCalcs_Oceano_Yue_3.mat');
StartTimes = [StartTimes; StartTimes_all{1}];
EndTimes = [EndTimes; EndTimes_all{1}];
Q = [Q; Q_all{1}];
sigma_Q = [sigma_Q; sigma_Q_all{1}];
tauRe = [tauRe; tauRe_all{1}];
sigma_tauRe = [sigma_tauRe; sigma_tauRe_all{1}];
ustRe = [ustRe; ustRe_all{1}];
sigma_ustRe = [sigma_ustRe; sigma_ustRe_all{1}];
zq = [zq; zq_all{1}];
sigma_zq = [sigma_zq; sigma_zq_all{1}];
theta = [theta; theta_all{1}];
zL = [zL; zL_all{1}];

%load and add fourth set of windows
load('/Users/raleighmartin/GitHub/AeolianFieldworkAnalysis/AnalysisData/Windowing/DataWindowCalcs_Oceano_Yue_4.mat');
StartTimes = [StartTimes; StartTimes_all{1}];
EndTimes = [EndTimes; EndTimes_all{1}];
Q = [Q; Q_all{1}];
sigma_Q = [sigma_Q; sigma_Q_all{1}];
tauRe = [tauRe; tauRe_all{1}];
sigma_tauRe = [sigma_tauRe; sigma_tauRe_all{1}];
ustRe = [ustRe; ustRe_all{1}];
sigma_ustRe = [sigma_ustRe; sigma_ustRe_all{1}];
zq = [zq; zq_all{1}];
sigma_zq = [sigma_zq; sigma_zq_all{1}];
theta = [theta; theta_all{1}];
zL = [zL; zL_all{1}];

%get dates
y = year(StartTimes);
m = month(StartTimes);
d = day(StartTimes);
Dates = datetime(y,m,d);

%adjust theta
theta = theta - theta_bias_Oceano;

%save files
save('DataWindowCalcs_Oceano_Yue.mat','Dates','StartTimes','EndTimes',...
    'Q','sigma_Q','tauRe','sigma_tauRe','ustRe','sigma_ustRe',...
    'zq','sigma_zq','theta','zL');