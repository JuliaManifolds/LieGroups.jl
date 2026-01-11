# A short Search for software packages

This is a short collection of packages I could find that define Lie groups in either of two criteria:

* The are a package that defines Lie groups explicitly and that is the main purpose of the package.
* The are a package that defines Lie groups as part of a larger purpose (e.g., differential geometry, algebraic groups, etc).

Among the packages, those that can be cited are preferred.

## DiffMan

A package developed about 1999-2001 to solve differential equations on manifolds. A main ingerdient there are Lie groups, which the package defines and implements.
The package is object oriented and implemented in Matlab.
[EngoeMarthinsenMunthe-Kaas:2001]

## Manifolds.jl

A direct predecessor of `LieGroups.jl` was the `GroupManifold` type in `Manifolds.jl` up to version 0.10. This type extended a manifold with a group operation, but stayed closer to the manifold perspective.

## Sophus

Is a package in C++ package, that provides headers for Lie groups.
There does not seem to be a (web-based/rendered/available) documentation. The Readme states, that the package is feature complete.

URL: https://github.com/strasdat/Sophus

Available Lie groups seem to be \(\mathrm{SO}(2), \mathrm{SO}(3), \mathrm{SE}(2), \mathrm{SE}(3)\).

## jaxlie

Inspired by `Sophus`, the package `jaxlie` implements Lie groups in Python using JAX for automatic differentiation and GPU acceleration,
though it provides only very specialized Lie groups \(\mathrm{SO}(n), \mathrm{SE}(n)\) for \(n=2,3\).
[YiLeeKlossMartinMartinBohg:2021]

## Lie++

A C++ Library again with a focus on robotics.
It seems to provide \(\mathrm{SO}(3)\) and \(\mathrm{SE}(3)\), but also the 2-product or n-product of \(\mathrm{SE}(3)\) an the Gallilean group

URL: https://github.com/aau-cns/Lie-plusplus
Citations: [fornasier2023msceqf], [fornasier2023equivariant]

## manif

A Python package that provides headers for Lie groups,
again \(\mathrm{SO}(2), \mathrm{SO}(3), \mathrm{SE}(2), \mathrm{SE}(3)\), the Gallilean group.
Since it is header only and the ¡ÈQuick start¡É links do not work, it is a bit unclear, how easy these are to use.

URL: https://github.com/artivis/manif
Ciations: They have a page, thatone could check.


## Further packages

### LiE

This seems to be or have been a Library in C.
There is a manual available at
http://www-math.univ-poitiers.fr/~maavl/LiEman/manual.pdf

but their project page linked in the document seems to be down.

### lie_learn

At first glance seems to provide one manifold, the sphere, and two Lie groups \(\mathrm{SO}(3)\) and \(\mathrm{SE}(3)\), but does not seem to provide a documentation not a way to cite the package. THey claim to be inspired by `LiE`.

URL: https://github.com/AMLab-Amsterdam/lie_learn/tree/master/lie_learn

### smooth

A C++ package for Lie theory with a focus on robtics, besides the usual metnioned Lie groups from before it also seems to provide \(\mathrm{SE}(3)\) with multiuple translations and products of Lie groups.

It is again a bit unclear how to cite this package

URL: https://github.com/pettni/smooth

### LieLab

From the code files it seems they also have \(\mathrm{SU}(n)\) and a ¡ÈSP¡É Lie group, that is a bit unclear, what it is, and general liear in both real and complex.

Citation: As a citation one should just cite the URL
URL: https://github.com/sandialabs/Lielab

### rotations

While it seems to be limited to rotation groups \(\mathrm{SO}(n)\), \(n=2,3\), this is added here, since it seems the only one in R that I could find.

Citation: Maybe the DOI from CRAN? https://cran.r-project.org/web/packages/rotations/index.html
URL: https://github.com/stanfill/rotationsC