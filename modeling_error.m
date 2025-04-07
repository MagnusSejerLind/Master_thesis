function [k,m,snr] = modeling_error(k,m)

rng("default")


error_mag = 0.05;
alpha_k = 1-error_mag + ((1+error_mag) - (1-error_mag)).*rand(size(k));  % From uniform distribution
alpha_m = 1-error_mag + ((1+error_mag) - (1-error_mag)).*rand(size(m));
k = k.*alpha_k;
m = m.*alpha_m;


snr = 25;   % signal to noise raio [db], % good value: 30
% snr = 'none';  % -- ~no noise

if error_mag == 0; disp('No modeling error implemented'); end
if snr == 'none'; disp('No measuring noise implemented'); end