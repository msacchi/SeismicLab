%% Alternating least squares for CANDECOMP/PARAFAC (CP) Decomposition
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="cp.html">CP Decompositions</a> 
% &#62;&#62; <a href="cp_als_doc.html">CP-ALS</a>
% </p>
% </html>
%
% The function |cp_als| computes an estimate of the best rank-R
% CP model of a tensor X using the well-known alternating least-squares
% algorithm (see, e.g., Kolda and Bader, SIAM Review, 2009, for more
% information).  The input X can be almost any type of tensor inclusing a 
% |tensor|, |sptensor|, |ktensor|, or |ttensor|. The output CP model is a
% |ktensor|. 

%% Load some data 
% We use the well-known _amino acids data set_ from Andersson and Bro.
% It contains fluorescence measurements of 5 samples containing 3 amino
% acids: Tryptophan, Tyrosine, and Phenylalanine.Each amino acid
% corresponds to a rank-one component. The tensor is of size 5 x 51 x 201
% from  5 samples, 51 excitations, and 201 emissions. 
% Further details can be found here: 
% <http://www.models.life.ku.dk/Amino_Acid_fluo>.
% Please cite the following paper for this data: 
% Rasmus Bro, PARAFAC: Tutorial and applications, Chemometrics and 
% Intelligent Laboratory Systems, 1997, 38, 149-171.  
% This dataset can be found in the |doc| directory.
load aminoacids
%% Basic call to the method, specifying the data tensor and its rank
% This uses a _random_ initial guess. At each iteration, it reports the 'fit'
% which is defined as |1-(norm(X-M)/norm(X))| and is loosely the proportion
% of the data described by the CP model, i.e., a fit of 1 is perfect.
rng('default') %<- Setting random seed for reproducibility of this script
M1 = cp_als(X,3); %<- Call the method
%%
% We typically can achieve a final fit of f = 0.97. The method stops when
% the change in the fit becomes less than the specified
% tolerance, which defaults to 1-e4. 

%% Visualize the results
% Use the |ktensor/viz| function to visualize the results.
vizopts = {'PlotCommands',{'bar','line','line'},...
    'ModeTitles',{'Concentration','Emission','Excitation'},...
    'BottomSpace',0.10,'HorzSpace',0.04,'Normalize',0};
info1 = viz(M1,'Figure',1,vizopts{:});

%% Run again with a different initial guess, output the initial guess.
% This time we have two outputs. The first output is the solution as a
% ktensor. The second output is a cell array containing the initial guess.
% Since the first mode is not needed, it is omitted from the cell array.
[M2bad,U2] = cp_als(X,3);

%% Increase the maximium number of iterations
% Note that the previous run kicked out at only 50 iterations before
% reaching the specified convegence tolerance. Let's increate the maximum
% number of iterations and try again, using the same initial guess.
M2 = cp_als(X,3,'maxiters',100,'init',U2);

%%
% This solution looks more or less the same as the previous one.
info2 = viz(M2,'Figure',2,vizopts{:});

%% Compare the two solutions 
% Use the |ktensor/score| function to compare the two solutions. A score of
% 1 indicates a perfect match. These are not exactly the same, but they are
% pretty close.
score(M1,M2)

%% Rerun with same initial guess
% Using the same initial guess (and all other parameters) gives the exact
% same solution. 
M2alt = cp_als(X,3,'maxiters',100,'init',U2);
score(M2, M2alt) %<- Score of 1 indicates the same solution

%% Changing the output frequency
% Using the |'printitn'| option to change the output frequency.
M2alt2 = cp_als(X,3,'maxiters',100,'init',U2,'printitn',10); 

%% Suppress all output
% Set |'printitn'| to zero to suppress all output.
M2alt3 = cp_als(X,3,'maxiters',100,'init',U2,'printitn',0); % <-No output

%% Use HOSVD initial guess
% Use the |'nvecs'| option to use the leading mode-n singular vectors as
% the initial guess.
M3 = cp_als(X,3,'init','nvecs','printitn',10);

%%
% Compare to the first solution using score, and see they are nearly the
% same because the score is close to 1.
score(M1,M3)

%% Change the order of the dimensions in CP
[M4,~,info] = cp_als(X,3,'dimorder',[2 3 1],'init','nvecs','printitn',10);
score(M1,M4)

%%
% In the last example, we also collected the third output argument which
% has some extra information in it. The field |info.iters| has the total
% number of iterations. The field |info.params| has the information used to
% run the method. Unless the initialization method is 'random', passing the
% parameters back to the method will yield the exact same results.
M4alt = cp_als(X,3,info.params);
score(M4,M4alt)

%% Change the tolerance
% It's also possible to loosen or tighten the tolerance on the change in
% the fit. You may need to increase the number of iterations for it to
% converge.
M5 = cp_als(X,3,'init','nvecs','tol',1e-6,'maxiters',1000,'printitn',10);

%% Control sign ambiguity of factor matrices
% The default behavior of |cp_als| is to make a call to |fixsigns| to fix
% the sign ambiguity of the factor matrices. You can turn off this behavior
% by passing the |'fixsigns'| parameter value of |false| when calling |cp_als|.
X = ktensor([1;1], {[1, 1; 1, -10],[1, 1; 1, -10]});
M = cp_als(X, 2, 'printitn', 0, 'init', X.U) % <-default behavior, fixsigns called
M = cp_als(X, 2, 'printitn', 0, 'init', X.U, 'fixsigns', false) % <-fixsigns not called

%% Recommendations
% * Run multiple times with different guesses and select the solution with
% the best fit. 
% * Try different ranks and choose the solution that is the best descriptor
% for your data based on the combination of the fit and the interpretaton
% of the factors, e.g., by visualizing the results.
