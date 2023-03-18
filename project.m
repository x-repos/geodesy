clc; clear; close;

% read data and pre-process
gpsdata = importdata("P780.tenv3");
t = gpsdata.data(:,1);          % time
east = gpsdata.data(:,7)*1000;
north = gpsdata.data(:,9)*1000;
up = gpsdata.data(:,11)*1000;
eaststd = gpsdata.data(:,13)*1000;
northstd = gpsdata.data(:,14)*1000;
upstd = gpsdata.data(:,15)*1000;



east = east - median(east);     % displacement
north = north - median(north);  % displacement
up = up - median(up);           % displacement
qeast = diag(eaststd.^2);   % covariance matrix
qnorth = diag(northstd.^2); % covariance matrix
qup = diag(upstd.^2);       % covariance matrix
eqtime = [2014.034, 2019.7290, 2020.0137, 2020.0164, 2020.0274]; % earthquake time
% process
[eEst, eResidual, epar, eStdEstPar] = lsquare(east, qeast, t, eqtime);
[nEst, nResidual, npar, nStdEstPar] = lsquare(north, qnorth, t, eqtime);
[uEst, uResidual, upar, uStdEstPar] = lsquare(up, qup, t, eqtime);
fprintf("parameters, rms of pars for East Displacement\n")
for i = 1:length(epar)
  fprintf("%15.5f\t%15.5f\n", epar(i), eStdEstPar(i));
end
fprintf("parameters, rms of pars for North Displacement\n")
for i = 1:length(npar)
  fprintf("%15.5f\t%15.5f\n", npar(i), nStdEstPar(i));
end
fprintf("parameters, rms of pars for Up Displacement\n")
for i = 1:length(upar)
  fprintf("%15.5f\t%15.5f\n", upar(i), uStdEstPar(i));
end


% plots fitting curve

textposition = [-60, -90, -35];

velocity(1) = num2str(epar(2)) + "+/-" + num2str(eStdEstPar(2));
velocity(2) = num2str(npar(2)) + "+/-" + num2str(nStdEstPar(2));
velocity(3) = num2str(upar(2)) + "+/-" + num2str(uStdEstPar(2));
ax(3) = 0;
for i=1:3
  ax(i) = subplot(3,1,i);
  xline(ax(i),eqtime,"k--",LineWidth=1)
  hold(ax(i), "on");
  pbaspect(ax(i),[10, 4, 1])
  box(ax(i),"on")
  xlim(ax(i), [min(t), max(t)])
  text(2016, textposition(i), "V = " + velocity(i))
end
ylabel(ax(1),"North (mm)")
ylabel(ax(2),"East (mm)")
ylabel(ax(3),"Up (mm)")
xlabel(ax(3),"Time (year)")
plot(ax(1), t, east, "bo","MarkerSize",2,"MarkerFaceColor","b")
plot(ax(2), t, north, "bo","MarkerSize",2,"MarkerFaceColor","b")
plot(ax(3), t, up, "bo","MarkerSize",2,"MarkerFaceColor","b")
plot(ax(1), t, eEst,LineWidth=2,Color="r");
plot(ax(2), t, nEst,LineWidth=2,Color="r");
plot(ax(3), t, uEst,LineWidth=2,Color="r");



% fitting function
function res = fitfunc(t, eqtime)
  n = length(t);
  o = 2*pi; % omega
  res = [ones(n,1), t, sin(o*t), cos(o*t), sin(2*o*t), cos(2*o*t)];
  for i = 1:length(eqtime)
    res = [res, stepFunc(t, eqtime(i))]; %#ok
  end
end

% least square scheme
function [res, residual, coef, stdEstPar] = lsquare(motion, Q, t, eqtime)
  A = fitfunc(t, eqtime);
  N = A'*Q^-1*A;
  coef = N^-1*A'*Q^-1*motion;
  A = fitfunc(t, eqtime); % earthquake
  res = A*coef;
  residual = motion - res;
  posteriorvar = (residual'*Q^-1*residual)./(length(t)-length(coef));
  corvar = posteriorvar*N^-1;
  stdEstPar = diag(corvar);
  stdEstPar = sqrt(stdEstPar);
  
end

% heaviside function
function res = stepFunc(t, tpoint)
  res = zeros(length(t), 1.);
  res = res(t<tpoint);
  res = [res; 0.5];
  n = length(res);
  res(n+1: length(t)) = 1.;
end
