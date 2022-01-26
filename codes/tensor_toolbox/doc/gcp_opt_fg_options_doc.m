%% Function Types for GCP
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="cp.html">CP Decompositions</a> 
% &#62;&#62; <a href="gcp_opt_doc.html">GCP-OPT</a>
% &#62;&#62; <a href="gcp_opt_gf_options_doc.html">Function Types for GCP</a>
% </p>
% </html>
%
% The GCP capability of the Tensor Toolbox allows the user to specify a fit
% function. There are a number of ''standard'' choices that we provide via
% the helper function |tt_gcp_fg_setup| function. These choices are
% presented in detail below. Motivations and details for these choices can
% be found in: 
%
% * D. Hong, T. G. Kolda, J. A. Duersch, Generalized Canonical
%   Polyadic Tensor Decomposition, SIAM Review, 62:133-163, 2020,
%   <https://doi.org/10.1137/18M1203626>
%
% These choices can be passed directly to gcp_opt via the 'type' option. To
% test the options, call the hidden function:
%
% |[f,g,lowerbnd] = tt_gcp_fg_setup(type)|
%
% We discuss the choices for the type below.
%% Gaussian (real-valued data)
% This is indicated by specifying the type as either |'normal'| or
% |'gaussian'|. This choice correspond to standard CP, which is implemented
% in  |cp_als| and |cp_opt|. It is useful for continuous real-valued data
% tensors. This choice specifies 
% 
% $$f(x,m) = (x-m)^2, \quad g(x,m) = 2(m-x), \quad \ell=-\infty$$
% 

[f,g,lowerbnd] = tt_gcp_fg_setup('normal')

%% Poisson (count data)
% This is indicated by specifying the type as either |'count'| or
% |'poisson'|. This choice is useful for count data tensors, i.e.,
% tensors that have only entries in {0,1,2,...}. This choice corresponds to
% Poisson CP, which is implemente din |cp_apr|. This choice specifies  
% 
% $$f(x,m) = m - x \log(m + 10^{-10}), 
% \quad g(x,m) = 1 - \frac{x}{m+10^{-10}}, 
% \quad \ell=0$$
% 
% The quantity $10^{-10}$ is a fudge factor to avoid divide-by-zero errors.

[f,g,lowerbnd] = tt_gcp_fg_setup('count')

%% Poisson with Log Link (count data)
% This is indicated by specifying the type as |'poisson-log'|. This choice
% is useful for count data tensors, i.e., tensors that have only entries in
% {0,1,2,...}.  This choice specifies    
% 
% $$f(x,m) = e^m - x m, 
% \quad g(x,m) = e^m - x, 
% \quad \ell=-\infty$$
% 

[f,g,lowerbnd] = tt_gcp_fg_setup('poisson-log')

%% Bernoulli with Odds Link (binary data)
% This is indicated by specifying the type as either |'binary'| or
% |'bernoulli-odds'|. This choice is useful for binary data tensors, i.e.,
% tensors that have only 0 or 1 entries. This choice specifies  
% 
% $$f(x,m) = \log(m+1) - x \log(m + 10^{-10}), 
% \quad g(x,m) = \frac{1}{m+1} - \frac{x}{m+10^{-10}}, 
% \quad \ell=0$$
% 
% The quantity $10^{-10}$ is a fudge factor to avoid divide-by-zero errors.

[f,g,lowerbnd] = tt_gcp_fg_setup('binary')

%% Bernoulli with Logit Link (binary data)
% This is indicated by specifying the type as |'bernoulli-logit'|. This
% choice is useful for binary data tensors, i.e., tensors that have only 0
% or 1 entries. This choice specifies   
% 
% $$f(x,m) = \log(e^m+1) - x m, 
% \quad g(x,m) = \frac{e^m}{e^m+1} - x, 
% \quad \ell=-\infty$$
% 

[f,g,lowerbnd] = tt_gcp_fg_setup('bernoulli-logit')

%% Rayleigh (real-valued data)
% This is indicated by specifying the type |'rayleigh'|. This choice is
% useful for nonnegative real-values data tensors, i.e., 
% tensors that have only nonnegative. This choice specifies  
% 
% $$f(x,m) = 2 \log(m+10^{-10}) - \frac{\pi}{4} \frac{x}{(m + 10^{-10})^2}, 
% \quad g(x,m) = \frac{1}{m+10^{-10}} - \frac{\pi}{2} \frac{x}{(m + 10^{-10})^3}, 
% \quad \ell=0$$
% 
% The quantity $10^{-10}$ is a fudge factor to avoid divide-by-zero errors.

[f,g,lowerbnd] = tt_gcp_fg_setup('rayleigh')

%% Gamma (nonnegative real-valued data)
% This is indicated by specifying the type |'gamma'|. This choice is
% useful for nonnegative real-values data tensors, i.e., 
% tensors that have only nonnegative. This choice specifies  
% 
% $$f(x,m) = \frac{x}{m+10^{-10}} + \log(m + 10^{-10}), 
% \quad g(x,m) = \frac{-x}{(m+10^{-10})^2} - \frac{1}{m + 10^{-10}}, 
% \quad \ell=0$$
% 
% The quantity $10^{-10}$ is a fudge factor to avoid divide-by-zero errors.

[f,g,lowerbnd] = tt_gcp_fg_setup('gamma')

%% Huber (nonnegative real-valued data)
% This is indicated by specifying the type |'huber (DELTA)'|, where |DELTA|
% is $\Delta$ in the equations below. This choice is useful for
% nonnegative real-values data tensors, i.e., tensors that 
% have only nonnegative. This choice specifies   
% 
% $$f(x,m) = \left\{ \begin{array}{ll}(x-m)^2 & \mbox{if } |x-m| \leq \Delta, \\ 
% 2\Delta|x-m|-\Delta^2 & \mbox{otherwise}\end{array}\right., 
% \quad
% g(x,m) = \left\{ \begin{array}{ll}-2(x-m) & \mbox{if } |x-m| \leq \Delta, \\ 
% 2\Delta\mbox{sgn}(x-m) & \mbox{otherwise}\end{array}\right., 
% \quad
% \ell = 0
% $$

[f,g,lowerbnd] = tt_gcp_fg_setup('huber (0.25)')

%% Negative Binomial (count data)
% This is indicated by specifying the type |'negative-binomial (r)'|, where |r|
% is $r$ in the equations below. This choice is useful for
% count data tensors. This choice specifies  
% 
% $$f(x,m) = (r+x) \log(1+m) - x \log(m+10^{-10}), 
% \quad
% g(x,m) = \frac{(r+x)}{1+m} - \frac{x}{m+10^{-10}}, 
% \quad
% \ell = 0
% $$

[f,g,lowerbnd] = tt_gcp_fg_setup('negative-binomial (4)')

%% Beta (nonnegative real-valued data)
% This is indicated by specifying the type |'beta (BETA)'|, where |BETA|
% is $\beta$ in the equations below. This choice is useful for
% nonnegative data tensors. Choices of $\beta=0$ or $\beta=1$ are not
% allowed because these correspond to 'gamma' or 'rayleigh'. 
% This choice specifies  
% 
% $$f(x,m) = \frac{ (m+10^{-10})^\beta }{\beta} - \frac{x(m+10^{-10})^{(\beta-1)} }{\beta-1}, 
% \quad
% g(x,m) = (m+10^{-10})^{(\beta-1)} - x(m+10^{-10})^{(\beta-2)}, 
% \quad
% \ell = 0
% $$

[f,g,lowerbnd] = tt_gcp_fg_setup('beta (0.3)')
