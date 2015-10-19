function [R lags T_L] = autocorrelation(X,dt,maxlag)

%set maxlag if not specified by user
if nargin == 2
  maxlag = 500;
end

%initialize output variables
R = zeros(maxlag+1,1);
lags = (0:maxlag)*dt;

%compute autocorrelation for each lag
Xbar = mean(X);
Xvar = var(X);
for i = 0:maxlag
    X0 = X(1:(end-i))-Xbar;
    X1 = X((1+i):end)-Xbar;
    R(i+1) = mean(X0.*X1)/Xvar;
end


%Lagrangian time scale computation

% %Method 1: perform integration - trapezoidal method
% i_max = min(find(R<0)); %index of first negative entry
% delta_t = diff(lags(1:i_max));
% R_delta_t = mean([R(1:i_max-1)';R(2:i_max)']);
% T_L = sum(delta_t.*R_delta_t); %Lagrangian timescale

% %Method 2: estimate based on first crossing of R=0.1
% i_crossing = min(find(R<0.1));
% T_L = lags(i_crossing)/log(5);

% %Method 3: estimate based on first crossing of R=0.1, or first positive
% %excursion, whichever comes first
% i_crossing = min(union(find(R<0.1),find(diff(R)>0)));
% R_crossing = R(i_crossing);
% T_L = lags(i_crossing)/log(1/R_crossing);

%Method 4: estimate based on first 0.5 seconds
T_max = 0.5; %max t, seconds
N = floor(T_max/dt);
T_L = mean(-lags(2:N+1)'./log(R(2:N+1)));

% %plot autocorrelation
% figure;
% plot(lags,R,'+')
% hold on;
% plot(lags,exp(-lags/T_L));
% xlabel('lags (s)');
% ylabel('autocorrelation');
% set(gca,'Yscale','log');
% legend('data','fit');