# Git model

The git branching strategy for developing and releasing the RASCAS code is the following. There are two main and protected branches. The **master** branch is the latest stable version. This is the one that you get when cloning the RASCAS repository. The **develop** branch is the branch where new developments are merged. When the new features in the **develop** branch have been tested, debugged, and used by as many users as possible, we merge the **develop** branch into the **master** branch and release a new stable version. Between two releases in the **master** branch, the only activity is bug fixes and documentation updates.

---

**A short cheat sheet to get the code**


Download the **master** branch:

```
git clone https://git-cral.univ-lyon1.fr/rascas/rascas.git
```

Download the **master** branch with ssh (with a gitlab account and SSH key): 

```
git clone git@git-cral.univ-lyon1.fr:rascas/rascas.git
```

Download the **develop** branch:

```
git clone --branch develop https://git-cral.univ-lyon1.fr/rascas/rascas.git
```

Go to the **develop** branch from the **master** branch: 

```
git checkout -b develop --track origin/develop
```

