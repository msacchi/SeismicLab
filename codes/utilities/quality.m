function Q = quality(d,d0);
%QUALITY: Compute quality of reconstruction in db
%
%  [Q] = quality(d,d0);
%
%  IN   d:  estimated signal
%       d0: true signal
%
%  OUT  Q:  reconstruction quality in db
%
%  Copyright (C) 2015, Signal Analysis and Imaging Group
%  For more information: http://seismic-lab.physics.ualberta.ca/
%  Author: M.D.Sacchi
%
%  SeismicLab is licensed under the MIT License.


e = (d-d0).^2;
Q = -10*log10(sum(e(:))/sum(d0(:).^2));

return;
