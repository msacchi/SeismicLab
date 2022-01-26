%SPIKING_DECON_DEMO: Apply spiking deconvolution to seismic traces.
%
%
%            -----------------------------------
%             A demo for Spiking Deconvolution
%
%                 M.D.Sacchi, SeismicLab
%            -----------------------------------
%
%
% Spiking deconvolution by-pass estimating the wavelet by
% assuming that the autocorrelation of the trace is
% an estimator of the autocorrelation of the wavelet. This
% is true if the reflectivity is white. In addition,
% the wavelet is assumed to be minimum phase. With these 
% assumptions spiking deconvolution should be able
% to remove the wavelet from the seismic trace.
%
% This demo uses the following SeismicLab functions:
% 
%               laplace_mixtue.m
%               spiking.m
%               taper.m
%               smooth_spectrum.m
%               plot_wb.m

  clear;
  close all;


% Read a min phase wavelet;

 [w,H] = readsegy('min_phase_wavelet.su');

% Reflectivity series

  nr = 400;
  r = laplace_mixture(nr, [0.001,0.1,0.7]);

% Trace

  s = conv(r,w);
  s = s(1:nr);
  s = taper(s,10,10);
  s = s +   0*max(s)*randn(size(s))/12;

  prewhitening = 0.1;                 % 0.1 percent
  lf = 50;                            % 50 points filter

 [f,o] = spiking(s, lf, prewhitening);


% Plot as wiggles but repeating the i/o 8 times

  figure(1);
  dt = 4./1000;
  time = (0:1:nr-1)*dt;

  subplot(211); plot_wb(time,[s,o],0);
  title('Seismogram (top: After decon, bottom: Input trace)');
  xlabel('Time (s)');

% Compute power spectrum to visualize gain in frequency
% content after spiking decon

  [Ps,f] = smooth_spectrum(s,dt,30,'li');
  [Po,f] = smooth_spectrum(o,dt,30,'li');

  subplot(212); plot(f,Po,f,Ps);
  title('Power spectum');
  xlabel('Frequency (Hz)');
  ylabel('Normalized PSD');
  legend('After', 'Before');
  grid
