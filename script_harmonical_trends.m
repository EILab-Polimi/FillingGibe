% script for fitting harmonical trend to cumulated precipitation (Fig 3a)
clc
clear

%% load and initialize data
addpath('data')
load('monthly_precOTB_Jan1999_Nov2018.mat') 
% contains:
% p: monthly average precipitation in the OTB, and 
% datem: montly calendar dates of the evaluation horizon
% datemExt: monthly calendar dates of the evaluation + extrapolation
% horizon as the harmonical trend is extrapolated 6 years in the future. 


Ncal = length(datem);    % length of available timeseries
Next = length(datemExt); % harmonical trend is evaluated for the available time series and extrapolated 6 years in the future as per vector datemExt
E = 3;                   % number of independent climate signals


%% estimate harmonical trends
four = cell(E,1);            % preallocate matrix hosting fourier estimation result
harm = cell(E,1);            % preallocate vectors of harminical trends
x = [1:Ncal]';               % harmonical trend is calibrated for the length of available timeseries
x2 = [1:Next]';              % harmonical trend is evaluated for the elongated time horizon
prec_residual = cell(E+1,1); % preallocate precipitation residual matrix
prec_residual{1} = p;        % initialize first value to precipitation 
h_sum = 0;                   % sum of harmonical trends

% recursively fit the precipitation residual to a fourier harmonic
for i = 1:E
    f = fit(x,prec_residual{i},'fourier1');                 
    four{i} = f;
    harm{i} = f.a0 + f.a1*cos(x2*f.w) + f.b1*sin(x2*f.w);
    h =  f.a0 + f.a1*cos(x*f.w) + f.b1*sin(x*f.w) ;
    prec_residual{i+1} = prec_residual{i} - h;
    h_sum = h_sum + harm{i};
end

% plot precipitation and sum of harmonics trend
figure; 
bar(datem, p, 'Facecolor', [.7 .7 .7]); hold on;
plot(datemExt, h_sum, 'r', 'LineWidth', 2)

legend('1 year cumulated precpitation','H = sum of harmonics') 
ylim([400 1400])
ylabel('Precipitation [mm]')
set(gca,'FontSize', 16)
xticks([datemExt(1:24:end)])
grid on
ax.XMinorGrid = 'on';






