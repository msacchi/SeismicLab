%% Optimization Methods for Tensor Toolbox
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="opt_options_doc.html">Optimization Methods</a> 
% </p>
% </html>
%
% Most MATLAB optimization methods have different interfaces or none at all.
% In Tensor Toolbox, we are adopting wrappers for moderate consistency, and
% these will be called from within other Tensor Toolbox functions. Here we
% outline the choices, their installation instructions, and then the
% user-tunable parameters.
%
% A few general notes:
%
% * These methods all require explicit initial guesses as well as handles for
%   the function and gradient calculations. These are handled by the calling routines.
% * The parameters marked with asterisks should generally not be modified
%   by the user because they will be set by the calling routine.
%
% For more information on the details of these methods, see
% <tt_opt_doc.html Developer Information for Optimization Methods in Tensor Toolbox>.

%% |lbfgsb|: Limited-Memory Quasi-Newton with Bound Constraints 
% In most methods, setting |'opt'| to |'lbfgsb'| will enable this method,
% and then any of the settings in the table below can be modified by
% passing these as additional options to the method.
%
% <html>
% <table>
% <tr><td>Name</td><td>Description</td><td>Default</td></tr>
% <tr><td><tt>lower</tt></td><td>Lower bounds, can be vector or scalar</td><td><tt>-Inf</tt></td></tr>
% <tr><td><tt>upper</tt></td><td>Upper bounds, can be vector or scalar</td><td><tt>+Inf</tt></td></tr>
% <tr><td><tt>maxiters</tt></td><td>Max outer iterations</td><td>1000</td></tr>
% <tr><td><tt>printitn</tt></td><td>Printing frequency by iteration (0=no output)</td><td>1</td></tr>
% <tr><td><tt>m</tt></td><td>Limited-memory parameter</td><td>5</td></tr>
% <tr><td><tt>subiters</tt></td><td>Controls maximum calls to function-gradient evalations</td><td>10</td></tr>
% <tr><td><tt>ftol</tt></td><td>Stopping condition based on relative function change</td><td>1e-10</td></tr>
% <tr><td><tt>gtol</tt></td><td>Stopping condition based on gradient norm</td><td>1e-5</td></tr>
% <tr><td><tt>mdesc</tt> (*)</td><td>Method description printed out before run</td><td><tt>'L-BFGS-B Optimization'</tt></td></tr>
% <tr><td><tt>xdesc</tt> (*)</td><td>Variable description</td><td>auto-generated</td></tr>
% </table>
% </html>
% 
% *Installation Instructions.*
% Download and install
% <https://github.com/stephenbeckr/L-BFGS-B-C *L-BFGS-B* by Stephen Becker>.   
% Please see that web page for full details on references, installation, etc. 
% Here we provide cursory instructions for installation:
% 
% # Download the zip file <https://github.com/stephenbeckr/L-BFGS-B-C/archive/master.zip https://github.com/stephenbeckr/L-BFGS-B-C/archive/master.zip>
% # Unzip and goto the |Matlab/| subdirectoy with MATLAB
% # Type |compile_mex|
% # _Add this directory to your saved path!_
%
% *Detailed notes.*
% The wrapper for this method is |tt_lbfgsb| in the Tensor Toolbox. 
% Notes regarding mappings to the parameters of Becker's L-BFGS-B code:
% 
% * |maxIts| maps to |maxiters| and the default is increased from 100 to 1000
% * |printEvery| maps to |printitn|
% * |maxTotalIts| is set to |maxiters*subiters| and this effectively changes the default from 5000 to 10000
% * |factr| is set to |ftol| / eps and this effectively changes the default from 1e7 to 4.5e5
% * |pgtol| is set to |gtol| 

%% |lbfgs|: Limited-Memory Quasi-Newton 
% In most methods, setting |'opt'| to |'lbfgs'| will enable this method,
% and then any of the settings in the table below can be modified by
% passing these as additional options to the method.
%
% <html>
% <table>
% <tr><td>Name</td><td>Description</td><td>Default</td></tr>
% <tr><td><tt>maxiters</tt></td><td>Max outer iterations</td><td>1000</td></tr>
% <tr><td><tt>printitn</tt></td><td>Printing frequency by iteration (0=no output)</td><td>1</td></tr>
% <tr><td><tt>m</tt></td><td>Limited-memory parameter</td><td>5</td></tr>
% <tr><td><tt>subiters</tt></td><td>Controls maximum calls to function-gradient evalations</td><td>10</td></tr>
% <tr><td><tt>ftol</tt></td><td>Stopping condition based on relative function change</td><td>1e-10</td></tr>
% <tr><td><tt>gtol</tt></td><td>Stopping condition based on gradient norm</td><td>1e-5</td></tr>
% <tr><td><tt>mdesc</tt> (*)</td><td>Method description printed out before run</td><td><tt>'Poblano L-BFGS Optimization'</tt></td></tr>
% <tr><td><tt>xdesc</tt> (*)</td><td>Variable description</td><td>auto-generated</td></tr>
% </table>
% </html>
% 
% *Installation Instructions.*
% Download and install
% <https://github.com/sandialabs/poblano_toolbox/releases/tag/v1.2 *Poblano Toolbox*, v1.2>.   
% Please see that web page for full details on references, installation, etc. 
% Here we provide cursory instructions for installation:
% 
% # Download the zip file <https://github.com/sandialabs/poblano_toolbox/archive/v1.2.zip https://github.com/sandialabs/poblano_toolbox/archive/v1.2.zip>
% # Unzip and goto the |poblano_toolbox-1.2/| subdirectoy within MATLAB
% # Type |install_poblano| to save this directory to your path
%
% *Detailed notes.*
% The wrapper for this method is |tt_lbfgs| in the Tensor Toolbox. 
% Notes regarding mappings to the parameters of Poblano's L-BFGS code:
% 
% * |MaxIters| maps to |maxiters|
% * |Display| maps to |printitn|
% * |MaxFuncEvals| is set to |maxiters*subiters| 
% * |RelFuncTol| is set to |ftol| 
% * |StopTol| is set to |gtol| 

%% |fminunc|: Optimizaton Toolbox Quasi-Newton Method
% In most methods, setting |'opt'| to |'fminunc'| will enable this method,
% and then any of the settings in the table below can be modified by
% passing these as additional options to the method. This requires the
% MATLAB Optimization Toolbox.
%
% <html>
% <table>
% <tr><td>Name</td><td>Description</td><td>Default</td></tr>
% <tr><td><tt>maxiters</tt></td><td>Max outer iterations</td><td>1000</td></tr>
% <tr><td><tt>printitn</tt></td><td>Display (0=no output)</td><td>1</td></tr>
% <tr><td><tt>subiters</tt></td><td>Controls maximum calls to function-gradient evalations</td><td>10</td></tr>
% <tr><td><tt>gtol</tt></td><td>Stopping condition based on gradient norm</td><td>1e-5</td></tr>
% <tr><td><tt>mdesc</tt> (*)</td><td>Method description printed out before run</td><td><tt>'Quasi-Newton Optimization (via Optimization Toolbox)'</tt></td></tr>
% <tr><td><tt>xdesc</tt> (*)</td><td>Variable description</td><td>auto-generated</td></tr>
% </table>
% </html>


%% |adam|: Stochastic Gradient Descent with Momentum
% In most methods, setting |'opt'| to |'adam'| will enable this method,
% and then any of the settings in the table below can be modified by
% passing these as additional options to the method.
%
% This is our own implementation of Adam. A _failed epoch_ is one where the
% function value does not decrease. After a failed epoch, the method either
% reduces the learning rate (by |decay|) or exits (once the number of
% failed epochs exceeds |maxfails|). 
%
%
% <html>
% <table>
% <tr><td>Name</td><td>Description</td><td>Default</td></tr>
% <tr><td><tt>lower</tt></td><td>Lower bounds, can be vector or scalar</td><td><tt>-Inf</tt></td></tr>
% <tr><td><tt>subiters</tt></td><td>Number of iterations per epoch</td><td>100</td></tr>
% <tr><td><tt>maxiters</tt></td><td>Maximum number of epochs</td><td>100</td></tr>
% <tr><td><tt>rate</tt></td><td>Initial learning rate</td><td>1e-2</td></tr>
% <tr><td><tt>maxfails</tt></td><td>Maximum number of failed epochs</td><td>1</td></tr>
% <tr><td><tt>decay</tt></td><td>Decay of learning rate after failed epoch</td><td>0.1</td></tr>
% <tr><td><tt>backup</tt></td><td>Revert to end of previous epoch after failure</td><td>true</td></tr>
% <tr><td><tt>ftol</tt></td><td>uit if function value goes below this value</td><td><tt>-Inf</tt></td></tr>
% <tr><td><tt>beta1</tt></td><td>Adam parameter</td><td>0.9</td></tr>
% <tr><td><tt>beta2</tt></td><td>Adam parameter</td><td>0.999</td></tr>
% <tr><td><tt>epsilon</tt></td><td>Adam parameter</td><td>1e-8</td></tr>
% <tr><td><tt>printitn</tt></td><td>Printing frequency by epoch (0=no output)</td><td>1</td></tr>
% <tr><td><tt>state</tt> (*)</td><td>State of random number generator</td><td>current state</td></tr>
% <tr><td><tt>mdesc</tt> (*)</td><td>Method description printed out before run</td><td><tt>'Adam Stochastic Optimization'</tt></td></tr>
% <tr><td><tt>xdesc</tt> (*)</td><td>Variable description</td><td>auto-generated</td></tr>
% <tr><td><tt>fdesc</tt> (*)</td><td>Description of (approximate) function computation</td><td>none</td></tr>
% <tr><td><tt>gdesc</tt> (*)</td><td>Description of stochastic gradient computation</td><td>none</td></tr>
% <tr><td><tt>fexact</tt> (*)</td><td>Boolean if function is computed exactly</td><td>true</td></tr>
% </table>
% </html>
% 



