# Tensor Toolbox Contribution Guide

## Getting Started

We do not accept Developer requests, but we do accept Reporter requests.
Before you can make a contribution, you need to create a **fork** of this repository.
Then you will be able to make a merge request following the instructions below.

## Checklist

- [ ] **Issue** Please link to any relevant issues that this merge request fixes.

- [ ] **Fork** Create a branch or fork of the code and make your changes. It's helpful if you create a branch on your fork.

- [ ] **Help Comments** Create or update comments for the m-files, following the style of the existing files. Be sure to explain all code options.

- [ ] **HTML Documentation** For any major new functionality, please follow the following steps.
  - [ ] Add HTML documentation in the `doc\html` directory with the name `XXX_doc.html`
  - [ ] Publish the documentation into `doc\html` via `cd doc; publish('XXX_doc.m','stylesheet','ttb.xsl')`
  - [ ] Add a pointer to this documentation file in `doc\html\helptoc.xml`
  - [ ] Add pointers in any related higher-level files, e.g., a new method for CP should be referenced in the `cp.html` file
  - [ ] Add link to HTML documentation from help comments in function
  - [ ] Update search database by running: builddocsearchdb('[full path to tensor_toolbox/doc/html directory]')
  
- [ ] **Tests** Create or update tests in the `tests` directory, especially for bug fixes or strongly encouraged for new code.

- [ ] **Contents** If new functions were added to a class, go to the `maintenance` directory and run `update_classlist('Class',XXX)` to add the new functions to the class XXX help information. If new functions were added at 
top level, go to `maintenance` and run `update_topcontents` to update the Contents.m file at the top level.

- [ ] **Release Notes** 
Update `README.txt` (under "Changes from [MOST RECENT VERSION]") with any significant bug fixes or additions.

- [ ] **Contributors List**
Update `CONTRIBUTORS.md` with your name and a brief description of the contribution.

- [ ] **Pass All Tests**
Confirm that all tests (including existing tests) pass in `tests` directory.

- [ ] **Merge Request** At any point, create a work-in-progress merge request, referencing the issue number and with this checklist and WIP in the header. To do this within the GITLAB website...
  * Start in your private branch
  * Go to Respoitory->Branches, select "Merge request"
  * On the "New Merge Request" screen, select "Change branches"
  * For the Target branch, enter `tensors/tensor_toolbox` and `master`
  * Give a description of the merge request, referencing any issues or other information.
  * Include _this_ checklist as the _first comment_ on the merge request.
  


