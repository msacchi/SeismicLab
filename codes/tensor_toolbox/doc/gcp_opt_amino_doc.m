%% GCP-OPT Examples with Amino Acids Dataset
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="cp.html">CP Decompositions</a> 
% &#62;&#62; <a href="gcp_opt_doc.html">GCP-OPT</a>
% &#62;&#62; <a href="gcp_opt_amino_doc.html">GCP-OPT and Amino Acids Dataset</a>
% </p>
% </html>
%

%% Setup
% We use the well known amino acids dataset for some tests. This data has
% some negative values, but the factorization itself should be nonnegative.

% Load the data
load(fullfile(getfield(what('tensor_toolbox'),'path'),'doc','aminoacids.mat'))

clear M fit

vizopts = {'PlotCommands',{@bar,@(x,y) plot(x,y,'r'),@(x,y) plot(x,y,'g')},...
    'BottomSpace',0.1, 'HorzSpace', 0.04, 'Normalize', @(x) normalize(x,'sort',2)};

%% CP-ALS
% Just a reminder of what CP-ALS does.

cnt = 1;

tic, M{cnt} = cp_als(X,3,'printitn',10); toc

fit(cnt) = 1 - norm(full(M{cnt})-X)/norm(X);
fprintf('Fit: %g\n', fit(cnt));

viz(M{cnt},'Figure',cnt,vizopts{:});

%% GCP with Gaussian
% We can instead call the GCP with the Gaussian function. 

cnt = 2;
M{cnt} = gcp_opt(X,3,'type','Gaussian','printitn',10);

fit(cnt) = 1 - norm(full(M{cnt})-X)/norm(X);
fprintf('Fit: %g\n', fit(cnt));

viz(M{cnt},'Figure',cnt,vizopts{:});

%% GCP with Gaussian and Missing Data
% What if some data is missing? 
cnt = 3;

% Proportion of missing data
p = 0.35; 

% Create a mask with the missing entries set to 0 and everything else 1
W = tensor(double(rand(size(X))>p)); 

% Fit the model, using the 'mask' option
M{cnt} = gcp_opt(X.*W,3,'type','Gaussian','mask',W,'printitn',10);

fit(cnt) = 1 - norm(full(M{cnt})-X)/norm(X);
fprintf('Fit: %g\n', fit(cnt));

viz(M{cnt},'Figure',cnt,vizopts{:});

%% GCP with ADAM
% We can also use stochastic gradient, though it's pretty slow for such a
% small tensor.
cnt = 4;

% Specify 'opt' = 'adam'
M{cnt} = gcp_opt(X,3,'type','Gaussian','opt','adam','printitn',1,'fsamp',5000,'gsamp',500);

fit(cnt) = 1 - norm(full(M{cnt})-X)/norm(X);
fprintf('Fit: %g\n', fit(cnt));

viz(M{cnt},'Figure',cnt,vizopts{:});

%% GCP with Gamma (terrible!)
% We can try Gamma, but it's not really the right distribution and produces
% a terrible result.
cnt = 5;

Y = tensor(X(:) .* (X(:) > 0), size(X));
M{cnt} = gcp_opt(Y,3,'type','Gamma','printitn',25);

fit(cnt) = 1 - norm(full(M{cnt})-X)/norm(X);
fprintf('Fit: %g\n', fit(cnt));

viz(M{cnt},'Figure',cnt,vizopts{:});

%% GCP with Huber + Lower Bound
% Huber works well. By default, Huber has no lower bound. To add one, we
% have to pass in the func/grad/lower information explicitly. We can use
% |gcp_fg_setup| to get the func/grad parameters.
cnt = 6;

% Call helper function tt_gcp_fg_setup to get the function and gradient handles
[fh,gh] = tt_gcp_fg_setup('Huber (0.25)');

% Pass the func/grad/lower explicitly.
M{cnt} = gcp_opt(X,3,'func',fh,'grad',gh,'lower',0,'printitn',25);

fit(cnt) = 1 - norm(full(M{cnt})-X)/norm(X);
fprintf('Fit: %g\n', fit(cnt));

viz(M{cnt},'Figure',cnt,vizopts{:});

%% GCP with Beta
% This is also pretty bad, which gives an idea of the struggle of choosing
% the wrong distribution. It can work a little bit, but it's clearly the
% wrong objective.
cnt = 7;

M{cnt} = gcp_opt(X,3,'type','beta (0.75)','printitn',25);

fit(cnt) = 1 - norm(full(M{cnt})-X)/norm(X);
fprintf('Fit: %g\n', fit(cnt));
viz(M{cnt},'Figure',cnt,vizopts{:});
