# Tests for Tensor Toolbox for MATLAB

Name each test file as `Test_Something.m`.

## Running the Tests
``` matlab
r=runtests; % Runs all tests
table(r) % View results
r=runtests(Test_Somthing); % Runs just the tests in that file
r=runtests('Test_NewTTM','ProcedureName','Compare','Verbosity',4); % Particular tests
```


## Creating New Tests
Copy one of the existing tests as a guide and modify to do tests relvant for the 
m-files being created. See the MATLAB documentation on `Class-Based Unit Tests` 
for more information.
