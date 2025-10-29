# Contributing to `LieGroups.jl`

First, thanks for taking the time to contribute.
We appreciate and welcome any contribution.

The following is a set of guidelines to [`LieGroups.jl`](https://juliamanifolds.github.io/LieGroups.jl/).

#### Table of contents

- [Contributing to `LieGroups.jl`](#Contributing-to-liegroupsjl)
      - [Table of Contents](#Table-of-Contents)
  - [How to just ask a question](#How-to-ask-a-question)
  - [How to file an issue](#How-to-file-an-issue)
  - [How to contribute](#How-to-contribute)
  - [Code style](#Code-style)

## How to just ask a question

You can most easily reach the developers in the Julia Slack channel [#manifolds](https://julialang.slack.com/archives/CP4QF0K5Z).
You can apply for the Julia Slack workspace [here](https://julialang.org/slack/) if you haven't joined yet.
You can also ask your question on [discourse.julialang.org](https://discourse.julialang.org).

## How to file an issue

If you found a bug or want to propose a feature, please open an issue in within the [GitHub repository](https://github.com/JuliaManifolds/LieGroups.jl/issues).

## How to contribute

Currently most details are still work-in-progress.
Feel free to contribute ideas, features you would like to see, Lie groups you want to have or would like to contribute, or any other idea for `LieGroups.jl`. For these, use either the [discussions](https://github.com/JuliaManifolds/LieGroups.jl/discussions) or [issues](https://github.com/JuliaManifolds/LieGroups.jl/issues) in the [GitHub repository](https://github.com/JuliaManifolds/LieGroups.jl)

## Code style

Please follow the [documentation guidelines](https://docs.julialang.org/en/v1/manual/documentation/) from the Julia documentation and use [Runic.jl](https://github.com/fredrikekre/Runic.jl) for code formatting.

Please consider a few internal conventions:

- Include the mathematical formulae for any implemented function if a closed form exists.
- Define a Lie group, a Lie group action, or a Lie algebra in its own file. Include all related functions in the same file
- an alphabetical order of functions in every file is preferable.
- The preceding implies that the mutating variant of a function follows the non-mutating variant.
- Document both the allocating and the mutating variants of a function. To avoid duplication, attach one doc string defined before both functions and attach it to both.
- There should be no dangling `=` signs.
- Add a newline between things of different types (struct/method/const).
- Add a newline between methods for different functions (including allocating/mutating variants).
- Prefer to have no newline between methods for the same function; when reasonable, merge the documentation strings.
- All `import`/`using`/`include` should be in the main module file.
- Avoid using `import` and use the explicit full name, like `Base.exp`, when implementing functions, that extend functions of other packages.
- if possible provide both mathematical formulae and literature references using [DocumenterCitations.jl](https://juliadocs.org/DocumenterCitations.jl/stable/) and BibTeX where possible
- Always document all input variables and keyword arguments