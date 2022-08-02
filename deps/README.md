Dependencies
============

Dependencies are maintained as [git subrepos](https://github.com/ingydotnet/git-subrepo)

```bash
c-dilks/largex-eic-asym
eic/delphes_EIC
delphes/delphes
```

Subdirectories are local copies of repositories, called "subrepos". They are automatically
included when you clone `sidis-eic`.

Maintainers of `sidis-eic` will be repsonsible for keeping these subrepos up-to-date.

Do not push changes to these subrepos, unless you are using `git subrepo`; see the next
section for details.


Notes For Maintainers
=====================

Each subrepo directory contains the file `.gitrepo`, which contains details
about the subrepo configuration. Use `git subrepo` commands when making changes
to any subrepos; run them from the top-level directory (`../`).

### Subrepo Updates
Pull subrepo updates, and push them to `sidis-eic` remote:
```bash
git checkout -b <new-branch-name>  # make a new sidis-eic branch

# pull subrepo updates with either of:
git subrepo pull deps/<subrepo>  # pull updates from a particular <subrepo>
git subrepo pull --all           # pull updates from all subrepos

git push -u origin <new-branch-name>  # push the new sidis-eic branch (which already has the subrepo pull commit(s))
```
Then open a pull request for `<new-branch-name> -> main`

### Subrepo Contributions
If you have `push` access to a `subrepo` remote, you can make your changes here. If you do not
have `push` access, you will need to fork and make your changes separately (or reconfigure subrepo remotes here).

Assuming you have `push` access to a subrepo `<subrepo>`:
```bash

# *make changes* #
git commit ...
git push ...  # push those changes to sidis-eic
git subrepo push deps/<subrepo>   # push your changes to the subrepo remote
```

