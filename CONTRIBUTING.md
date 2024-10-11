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
Run [`JuliaFormatter.jl`](https://github.com/domluna/JuliaFormatter.jl) on the repository in the way set in the `.JuliaFormatter.toml` file, which enforces a number of conventions consistent with the Blue Style.

Please follow a few internal conventions:

- Any implemented function should be accompanied by its mathematical formulae if a closed form exists.
- A Lie group, a Lie group action or a Lie algebra should be implemented in it own file
- an alphabetical order of functions in every file is preferable.
- The preceding implies that the mutating variant of a function follows the non-mutating variant.
- usually, both the allocating and the mutating variants of a function should be documented. For avoiding duplication, one doc string attached to both is preferrable
- There should be no dangling `=` signs.
- Always add a newline between things of different types (struct/method/const).
- Always add a newline between methods for different functions (including allocating/mutating variants).
- Prefer to have no newline between methods for the same function; when reasonable, merge the documentation strings.
- All `import`/`using`/`include` should be in the main module file.
- `import` is discouraged and the explicit full name, like `Base.exp` should be used when implementing functions

Concerning documentation

- if possible provide both mathematical formulae and literature references using [DocumenterCitations.jl](https://juliadocs.org/DocumenterCitations.jl/stable/) and BibTeX where possible
- Always document all input variables and keyword arguments