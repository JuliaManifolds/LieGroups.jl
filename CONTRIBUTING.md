# Contributing to `LieGroups.jl`

First, thanks for taking the time to contribute.
Any contribution is appreciated and welcome.

The following is a set of guidelines to [`LieGroups.jl`](https://juliamanifolds.github.io/LieGroups.jl/).

#### Table of contents

- [Contributing to `LieGroups.jl`](#Contributing-to-manoptjl)
      - [Table of Contents](#Table-of-Contents)
  - [I just have a question](#I-just-have-a-question)
  - [How can I file an issue?](#How-can-I-file-an-issue)
  - [How can I contribute?](#How-can-I-contribute)
  - [Code style](#Code-style)

## I just have a question

The developer can most easily be reached in the Julia Slack channel [#manifolds](https://julialang.slack.com/archives/CP4QF0K5Z).
You can apply for the Julia Slack workspace [here](https://julialang.org/slack/) if you haven't joined yet.
You can also ask your question on [discourse.julialang.org](https://discourse.julialang.org).

## How can I file an issue?

If you found a bug or want to propose a feature, please open an issue in within the [GitHub repository](https://github.com/JuliaManifolds/LieGroups.jl/issues).

## How can I contribute?

Currently we are still consolidating most details and defaults, but feel free to contribute to these discussions within this repository.

## Code style

Try to follow the [documentation guidelines](https://docs.julialang.org/en/v1/manual/documentation/) from the Julia documentation as well as [Blue Style](https://github.com/invenia/BlueStyle).
Run [`JuliaFormatter.jl`](https://github.com/domluna/JuliaFormatter.jl) on the repository in the way set in the `.JuliaFormatter.toml` file, which enforces a number of conventions consistent with the Blue Style. Furthermore [vale](https://vale.sh) is run on both Markdown and code files, affecting documentation and source code comments

Please follow a few internal conventions:

- It is preferred that the `AbstractManoptProblem`'s struct contains information about the general structure of the problem.
- Any implemented function should be accompanied by its mathematical formulae if a closed form exists.
- `AbstractManoptProblem` and helping functions are stored within the `plan/` folder and sorted by properties of the problem and/or solver at hand.
- the solver state is usually stored with the solver itself
- Within the source code of one algorithm, following the state, the high level interface should be next, then the initialization, then the step.
- Otherwise an alphabetical order of functions is preferable.
- The preceding implies that the mutating variant of a function follows the non-mutating variant.
- There should be no dangling `=` signs.
- Always add a newline between things of different types (struct/method/const).
- Always add a newline between methods for different functions (including mutating/nonmutating variants).
- Prefer to have no newline between methods for the same function; when reasonable, merge the documentation strings.
- All `import`/`using`/`include` should be in the main module file.

Concerning documentation

- if possible provide both mathematical formulae and literature references using [DocumenterCitations.jl](https://juliadocs.org/DocumenterCitations.jl/stable/) and BibTeX where possible
- Always document all input variables and keyword arguments

If you implement an algorithm with a certain numerical example in mind, it would be great, if this could be added to the [ManoptExamples.jl](https://github.com/JuliaManifolds/ManoptExamples.jl) package as well.