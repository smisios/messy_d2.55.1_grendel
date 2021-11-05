# How to contribute to the MESSy development in the GitLab repository

## Basics

- If you are new to git, read https://git-scm.com/book/en/v2 first.
- Do not mix git with GitLab! GitLab is more, see https://gitlab.com/

## Branches

We have four protected branches:

 - *master*: contains only tagged official releases which are tagged accordingly
 - *rc*: (release candidate) contains a close-to-release code with feature freeze and open only for bug-fixes
 - *hotfix*: bug-fixes and machine environment adaptions
 - *devel*: main developing branch towards next release

For the further development and sharing the work additional branches 
(branching off from devel) must be opened, i.e.

 - XYZ-*: branches off from devel and deals with solving issue
   number XYZ

These branches will be closed after successful merge into devel.

## Workflow

### First time

 1. Setup your local git environment:

    % git config --global user.name "John Doe"

    % git config --global user.email john.doe@institution.org

 2. Starting from scratch, the central repository needs to be cloned once
    on your local machine:

    % git clone git@gitlab.dkrz.de:MESSy/MESSy.git

    This will generate the subdirectory MESSy.

 3. Check-out the branch 'devel':

    % git checkout devel

### Daily workflow

 1. Change to your local clone of the repository

    % cd MESSy

 2. Update changes from the central repository:

    % git pull

 3. Either checkout an existing issue branch and work with it (see below),
    or create a new issue on the GitLab web server and open
    a branch connected to it:  
    Click on the arrow next to "Create Merge Request", name new branch
    according to the issue (XYZ) and make sure it branches off
    from devel.

 4. Checkout new issue branch in your local repository:

    % git pull  
    % git checkout XYZ  
    % git pull

 5. Update your issue branch regularly with other developments
    from the devel branch (after you pulled also devel):

    % git merge devel

 6. Resolve potential conflicts and commit them locally:

    % git add [...]

    % git commit 

 7. ... work on code and test it ...

    % git status 

    % git add [...]

    % git commit [...]

 8. Push the changes in your branch to the central repository:

    % git push

 9. Once the issue is solved and all changes of the issue branch
    have been pushed to central repository, create a merge request into the
    devel branch on the GitLab web server and assign the reviewers.
    The reviewers will check your code modifications and potentially
    comment or revise it and finally accept them to be merged.

    For this, click on "Merge Requests", select "New Merge Request", select
    source branch (XYZ) and destination branch (devel), click on
    "Compare branches and continue", fill in form (maybe remove "WIP:"),
    name co-referees (by @-mentions in the comment), assign it the
    main referee and submit it.

    Important:
    Each merge rquest triggers a GitLab CI pipeline job to automatically
    build the distribution on a mistral node. You can see this
    directly on the merge-request page or under CI/CD. The build test takes
    almost one hour. If it fails, you can download and inspect the so-called
    job artifacts (basically a zip file containing the log-output of the
    build). Due to storage restrictions on the node, the artifacts
    will be deleted automatically after 12 hours.
    **Note: Each time you push a change to an open merge request, a new
    pipeline job will be triggered.**

### Commit messages

Make sure your commit messages contain:

 1. A first line, starting with the issue number (plus colon) folled by a meaningful description of what has been done.

 2. An empty second line.

 3. Additional explanatory text with as much details as desired from the
    third line onwards.

    * Note that there is no need to list modified files, because
      git traces these anyway.

    * The text should also contain the issue number in the form
       issue #XYZ.
