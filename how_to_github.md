# Basic git commands 

```
$ git clone [put repo url here]

$ git status # Now you should see that you are on the master branch.

$ git checkout -b [any branch name] # -b means to create a new repo 

$ git status # Now you should see that you have been switched to the created branch. 

Then you can make any modification of the files now. 

$ git status # After modifying the files, you will see that the github system will notice some file has been changed

$ git add [file name] # This is to put the given file "on stage" for commits. 

$ git add --all # This is to stage all the untracked files. 

$ git reset -- [file] # "unstage" the file while keeping the modifications since the last commit. 

$ git commit -m "Any message you want to add for this commit" # Commit the changes. 

$ git commit -a -m "abc" # commit all the unstaged modifications directly with comments "abc". 

$ git revert --no-commit [commit hash]..HEAD
$ git commit        # revert to the status corresponding to the [commit hash] with history intact. 

$ git push origin [current branch name] # Push the new branch to the remote repo for review and merge. 

$ git diff [file] # check what modifications you made in [file] since the last commit. 

$ git checkout . # undo all the uncomitted changes. 

$ git submodule add [module git url] # add a submodule to the current module

$ git submodule init

$ git submodule update --init --recursive # Besides recursive update, it will also recursively initialize all the uninitialized modules. 
```
If you want to push local branch A to remote branch B, make sure you first have local branch B, and then push. 


