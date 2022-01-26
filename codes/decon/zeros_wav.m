function [z] = zeros_wav(w);
%ZEROS_WAV: computes the zeros of a wavelet.
%
%  [z] = zeros(w);
%
%  IN   w: a wavelet;
%
%  OUT  z: complex zeros
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


[N,M] = size(w);
if N == 1; wr = fliplr(w); z = roots(wr); end;
if M == 1; wr = flipud(w); z = roots(wr); end;

return;
 

