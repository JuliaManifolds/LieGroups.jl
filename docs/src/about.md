# About LieGroups.jl

the package `LieGroups.jl` in its current form started in 2024 as a “spin-off” from [`Manifolds.jl](@extref Manifolds :std:doc:`index`), where a type `GroupManifolds` was a predecessor of the
main [`LieGroup`](@ref) type of this package. That approach started around 2021 by [Seth Axen](https://sethaxen.com).
At about the same, [Yueh-Hua Tu](https://github.com/yuehhua) started a package `LieGroups.jl`,
which was continued than here with a full rewrite in order to base all Lie groups on a
manifold from [`Manifolds.jl](@extref Manifolds :std:doc:`index`).

The current main developers are [Ronny Bergmann](https://ronnybergmann.net) and [Mateusz Baran](https://github.com/mateuszbaran).

## Contributors

Thanks to the following contributors to `LieGroups.jl`:

* [Paula Conrad](https://github.com/plc99) and [Leonard Schreiter](https://github.com/LeoSchreiter) implemented the [Circle group](@ref CircleGroup) as part of their student assistant project.
* [Olivier Verdier](https://www.olivierverdier.com) helped in the design and some mathematical explanations

as well as [contributors](https://github.com/JuliaManifolds/LieGroups.jl/graphs/contributors) providing small extensions, finding small bugs and mistakes and fixing them by opening [PR](https://github.com/JuliaManifolds/LieGroups.jl/pulls)s. Thanks to all of you.

If you want to contribute a manifold or algorithm or have any questions, visit
the [GitHub repository](https://github.com/JuliaManifolds/LieGroups.jl/)
to clone/fork the repository or open an issue.

## Work using LieGroups.jl

If `LieGroups.jl` is useful within another package or in your project, we would like to list that here.
Please [open an issue](https://github.com/JuliaManifolds/LieGroups.jl/issues/new).
It would be great to collect anything and anyone using `LieGroups.jl` in this list.