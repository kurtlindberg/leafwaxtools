.. _contributing_to_leafwaxtools:

#############################
Contributing to leafwaxtools
#############################


GitHub
======
All contributions, bug reports, bug fixes, documentation improvements, enhacenments,
and ideas are welcome, and take place through `GitHub <https://github.com/kurtlindberg/leafwaxtools/issues>`_

There are several levels of contributions to an open development software package like leafwaxtools, incuding:

#.  Reporting bugs
#.  Updating the documentation
#.  Updating existing functionalities
#.  Contributing new functionalities

All of that takes place through GitHub `issues <https://docs.github.com/en/issues/tracking-your-work-with-issues/quickstart>`_.

When you start working on an issue, it's a good idea to assign the issue to yourself, again to limit duplication. If you can't think of an issue on your own, we have you covered:  check the list of unassigned issues and assign yourself one you like.
If, for whatever reason, you are not able to continue working with the issue, please tyo to unassign it, so other people know it's available again. You can chekc the list of assigned issues, since people may not be working on them anymore. If you want to work on one that is assigned, feel free to kindly ask the current assignee (on GitHub) if you can take it (please allow at least a week of inactivity before considering work in the issue discontinued).

Bug reports and enhancement requests
====================================

Bug reports are an important part of improving any software. Having a complete bug report will allow others to reproduce the bug and provide insight into fixing. See this `stackoverflow article <https://stackoverflow.com/help/mcve>`_ and `this blog post <https://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports>`_ for tips on writing a good bug report.
Trying the bug-producing code out on the master branch is often a worthwhile exercise to confirm the bug still exists. It is also worth searching existing bug reports and pull requests to see if the issue has already been reported and/or fixed.
Bug reports must:

1. Include a minimal working example (a short, self-contained Python snippet reproducing the problem). You can format the code nicely by using GitHub-flavored Markdown::

   test_df = leafwax.WaxData(test_data)
   test_total_conc = test_df.acl()

2. Include the full version string of leafwaxtools, which you can obtain through::

   leafwax.__version__

3. Explain why the current behavior is wrong/not desired and what you expect or would like to see instead.


Working with the leafwaxtools codebase
=======================================
Version control, Git, and GitHub
""""""""""""""""""""""""""""""""

To the new programmer, working with Git is one of the more daunting aspects of contributing to open ousrce projects.
It Cna very quickly become overwhelming, but sticking to the guidelines below will help keep the process staightforward and mostly trouble free. As always, if you are having difficulties, plesae feel free to ask for help.
The code is hosted on `GitHub <https://github.com/kurtlindberg/leafwaxtools>`_. To contribute, you will need to `sing up for a (free) GitHub account <https://github.com/signup/free>`_. `Git <https://git*scm.com/>`_ is the industry standard for version control to allow many people to work together on the project, keep track of issues, manage the project, and much more.

Some great resources for learning Git:
  * the `GitHub help pages <https://help.github.com/>`_
  * the `Numpy documentation <https://numpy.org/doc/stable/dev/index.html>`_
  * Matthew Brett's `Pydagogue <https://matthew-brett.github.io/pydagogue/>`_

GitHub has `instructions <https://help.github.com/set-up-git-redirect>`_ for installing git, settingup your SSH key, and configuring git. All these steps need to be completed before you can work seamlessly between your local repository and GitHub.

Forking
"""""""
You will need your own fork to work on the code. Go to the leafwaxtools repository and hit the Fork button. You will then want to clone your fork (i.e. download all the code to your local machine so you can edit it locally).
At the command line. this would look something like::

    git clone https://github.com/your-user-name/leafwaxtools.git leafwaxtools-yourname
    cd leafwaxtools-yourname
    git remote add upstream https://github.com/kurtlindberg/leafwaxtools.git

This creates the directory `leafwaxtools-yourname` and connects your repository to the upstream (main project) leafwaxtools repository. However, most Git first-timers may find it easier to do so thorugh the GitHub web interface or desktop app (where there is a proverbial "button for that").

Creating a development environment
""""""""""""""""""""""""""""""""""
We recommend developing in the same conda environment in which you installed leafwaxtools.

Creating a branch
"""""""""""""""""
You want your master branch to reflect only production-ready code, so create a feature branch for making your changes. For example::

    git branch shiny-new-feature
    git checkout shiny-new-feature

The above can be simplified to::

    git checkout -b shiny-new-feature

This changes your working directory to the `shiny-new-feature` branch. Keep any changes in this branch specific to one bug or feature so it is clear what the branch rbrings to leafwaxtools. You can have many `shiny-new-features` and switch between them using the `git checkout` command.
When creating this branch, make sure your master branch is up to date with the latest upstream master version. To update your local master branch, you can do::

    git checkout main
    git pull upstream master --ff-only

When you want to update the feature branch with change in master after you created the branch, check the section on updating a pull request.

leafwaxtools Protocol
""""""""""""""""""""""

Contributing new functionalities
********************************

1.  Open an issue on GitHub (See above)
2.  Implement outside of leafwaxtools

    Before incorporating any code into leafwaxtools, make sure you have a solution that works outside leafwaxtools. Demonstrate this in a notebook, which can be hosted on GitHub, as well, so it is easy for the maintainers to check out. The notebook should be organized as follows:

    - dependencies (package names and versions)
    - body of the function
    - example usage
3.  Integrate the new functionality

    Now you may implement the new functionality inside leafwaxtools. In so doing. make sure you:

    * Re-use as many of leafwaxtools' existing utilites as you can, introducing new package dependencies only as necessary.
    * Create a docstrong fro your new function, describing arguments and returned variables, and showing an example of use. (Use an existing docstring for inspiration).
    * If possible, also include a unit test for `continuous integration <https://youtu.be/_WvjhrZR01U>`_ (leafwaxtools uses `pytest`). Feel free to ask for help from the package developers.

4.  Expose the new functionality in the leafwaxtools user API (files located in the `core` folder)


Updating existing functionalities
**********************************

1.  Open an issue on GitHub (same advice as above)
2.  Implement outside of leafwaxtools, including a benchmark of how the existing function performs vs the proposed upgrade (e.g. with `timeit`). Take into consideration memory requirements and describe on what architechture/OS you ran the test.
3.  Integrate the new functionality within leafwaxtools (same advice as above)
4. Update the unit test(s) to make sure they still pass muster. Depending on the complexity of the feature, there may be more than one test to update.

Testing
"""""""

Testing is hugely important, as you don't want your "upgrades" to break the whole package by introducing errors. Thankfully, there is a proverbial app for that: *unit testing*. Write a test of your code using the naming rules:

1.  class: `Test{filename}{Class}{method}` with appropriate camel case convention
2.  function: `test_{method}_t{test_id}`

(see e.g. test_api_WaxData.py for example)

Your test should be as minimal as possible; it is aimed to see if the function you wrote/updated works as advertised given a reasonably comprehensive list of possible arguments. leafwaxtools' tests rely on data already included in the data directory, and we strongly recommend that you do the same; inly introduce a new dataset if the existing ones are insufficient to properly test your code. In general, the simpler the test, the better, as it will run in less time.

To run the test(s):

1.  Make sure the `pytest package <https://docs.pytest.org>`_ is installed on your system; run `pip install pytest` if not.
2.  In your terminal, switch to the "tests" subdirectory of your leafwaxtools forked repository. If you wish to test a specific class/method inside a specified file, run `pytest {file_path}\::{TestClass}\::{test_method}`
3.  To run *all* tests in the specified file, run `pytest {file_path}`
4.  To perform all test in all testing files inside the specified directory, execute `pytest {directory_path}`

To order above is somewhat loose, but goes from least complex (time-consuming) to more complex.


Stylistic considerations
""""""""""""""""""""""""
Guide van Rossum's great insight is that code is read far more often than it is written, so it is important for the code to be of a somewhat uniform style, so that people can read and understand it with relative ease. leafwaxtools strives to use fairly consistent notation, including:

  * Capital letters for matrices, lowercase for vectors
  * Independent variable is called ys, the dependent variable xs.
  * Function names use CamelCase convention

Conventions
"""""""""""

- leafwaxtools uses Numpy doc for documentation

Contributing your changes to leafwaxtools
==========================================

Comitting your code
"""""""""""""""""""
Once you've made changes, you can see them by typing::

    git status

If you created a new file, it is not being tracked by git. Add it by typing::

    git add path/to/file-to-be-added.py

Typing `git status` again should give something like::

    On branch shiny-new-feature
    modified:   /relative/path/to/file-you-added.py

Finally, commit your changes to your local repository with an explanatory message. The message need not be encyclopedic, but it should say what you did, what GitHub issue it refers to, and what part of the code it is expected to affect.
The preferred style is:

  * A subject line with < 80 characters.
  * One blank line.
  * Optionally, a commit message body.

Now you can commit your changes in your local repository::

    git commit -m 'type your message here'

Pushing your changes
""""""""""""""""""""

When you want your changes to appear publicly on your GitHub page, push your forked feature branch's commits::

    git push origin shiny-new-feature

Here `origin` is the default name given to your remote repository on GitHub. You can see the remote repositories::

    git remote -v

If you added the upstream repository as described above, you will see something like::

    origin  git@github.com:yourname/leafwaxtools.git (fetch)
    origin  git@github.com:yourname/leafwaxtools.git (push)
    upstream    git://github.com/kurtlindberg/leafwaxtools.git (fetch)
    upstream    git://github.com/kurtlindberg/leafwaxtools.git (push)

Now your code is on GitHub, but it is not yet a part of the leafwaxtools project. For that to happen, a pull request needs to be submitted on GitHub.

Filing a Pull Request
"""""""""""""""""""""
When you're ready to ask for a code review, file a pull request. But before you do, please double-check that you have followed all the guidelines outlined in this document regarding code style, tests, performace tests, and documentation. You should also double check your branch changes against the branch it was based on:

  * Navigate to your repository on GitHub
  * Click on Branches
  * Click on the Compare button for your feature branch
  * Selct the base and compare branches, if necessary. This will be *Development* and *shiny-new-feature*, respectively.

If everything looks good, you are ready to make a pull request. A pull request is how code fomr a local repository becomes available to the GitHub community and can be reviewed by a project's owners/developers and eventually merged in the master version. This pull request and its associated changes will eventually be committed to the master branch and available in the next release. To submit a pull request:

  * Navigate to your repository on GitHub
  * Click on the Pull Request button
  * You can then click on Commits and Files Changed to make sure everything looks okay one last time
  * Write a description of your change in the Preview Discussion tab
  * Click Send Pull Request

This request then goes to the repository maintainers, and they will review the code.

Updating your pull request
""""""""""""""""""""""""""

Based on the review you get on your pull request, you will probably need to make some change to the code. In that case, you can make them in your branch, add a new commit to that branch, push it to GitHub, and the pull request will be automatically updated. Pushing them to GitHub again is done by::

    git push origin shiny-new-feature

This will automatically update your pull request with the latest code and restart the Continuous Integration tests (which is why is is important to provide a test for your code).
Another reason you might need to update your pull request is to solve conflicts with changes that have been merged into the main branch since you opened you pull request.
To do this, you need to `merge upstream main` in your branch::

    git checkout shiny-new-feature
    git fetch upstream
    git merge upstream/master

If there are no conflicts (or the could be fixed automatically), a file with a default commit message will open, and you can simply save and quit this file.
If there are merge conflicts, you need to solve those conflicts. See `this example <https://help.github.com/articles/resolving-a-merge-conflict-using-the-command-line/>`_ for an explanation on how to do this. Once the conflicts are merged and the files where the conflicts were solved are added, you can run git commit to save those fixes.
If you have uncommitted changes at the moment you want to update the branch with master, you will need to stash them prior to updating (see the stash docs). This will effectively store your changes and they can be reapplied after updating.
After the feature branch has been updated locally, you can now update your pull request by pushing to the branch on GitHub::

    git push origin shiny-new-feature

Delete your merged branch (optional)
""""""""""""""""""""""""""""""""""""

Once your feature branch is accepted into upstream, you'll probably want to get rid of the branch. First, merge upstream master into your branch so git knows it is safe to delete your branch::

    git fetch upstream
    git checkout main
    git merge upstream/main

Then you can do::

    git branch -d shiny-new-feature

Make sure you use a lowercase `-d`, or else git won't warn you if your feature branch has not actually been merged.
The branch will still exist on GitHub, so to delete it there do::

    git push origin --delete shiny-new-feature

Tips for a successful pull request
""""""""""""""""""""""""""""""""""
If you have made it to the "Review your code" phase, one of the core contributors will take a look. Please note, however, that reponse time will be variable (e.g. don't try the week before AGU).
To improve the change of your pull request being reviewed, you should:

  * Reference an open issue for non-trivial changes to clarify the PR's purpose
  * Ensure you have apporpriate tests. These should be the first part of any PR
  * Keep your pull requests as simple as possible. Largers PRs take longer to review
  * If you need to add on to what you submitted, keep updating your original pull request, either by request or every few days

Documentation
=============

About the leafwaxtools documentation
"""""""""""""""""""""""""""""""""""""
leafwaxtools' documentation is built automatically from the function and class docstrings, via `Read The Docs <https://readthedocs.org>`_. It is, therefore, especially important for your code to include a docstring, and to modify the docstrings of the functions/classes you modify to make sure the documentation is current.

Updating a leafwaxtools docstring
""""""""""""""""""""""""""""""""""
You may use existing docstrings as examples. A good docstring explains:

  * What the function/class is about
  * What it does, with what properties/inputs/outputs
  * How to use it, via a minimal working exmaple

For the latter, make sure the example is prefaced by:

    .. jupyter-execute::

and properly indented (look at other docstrings for inspiration).

How to build the leafwaxtools documentation
""""""""""""""""""""""""""""""""""""""""""""

Navigate to the leafwaxtools/docs directory and type `make html`. This may require installing other packages (sphinx, numpydoc, nbsphinx, sphinx_search, jupyter-sphinx, sphinx-copybutton, sphinx_rtd_theme).


You are done! Thanks for reading.
