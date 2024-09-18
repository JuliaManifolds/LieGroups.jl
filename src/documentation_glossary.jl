#
# Lie groups Glossary
# ===
#
# This file collects
# * LaTeX snippets
# * math formulae
# * Variable names
# * links
# * notes
#
# to keep naming, notation, and formatting in a unifed way

# In general every dictionary here can be either :Symbol-> String or :Symbol -> Dictionary enrties

_LIEGROUPS_DOC_TYPE = Dict{Symbol,Union{String,Dict,Function}}

_manopt_glossary = _LIEGROUPS_DOC_TYPE()

# easier access functions
"""
    glossary(s::Symbol, args...; kwargs...)
    glossary(g::Dict, s::Symbol, args...; kwargs...)

Access an entry in the glossary at `Symbol` s
if that entrs is
* a string, this is returned
* a function, it is called with `args...` and `kwargs...` passed
* a dictionary, then the arguments and keyword arguments are passed to this dictionary, assuming `args[1]` is a symbol
"""
#do not document for now, until we have an internals section
glossary(s::Symbol, args...; kwargs...) = glossary(_manopt_glossary, s, args...; kwargs...)
function glossary(g::_LIEGROUPS_DOC_TYPE, s::Symbol, args...; kwargs...)
    return glossary(g[s], args...; kwargs...)
end
glossary(s::String, args...; kwargs...) = s
glossary(f::Function, args...; kwargs...) = f(args...; kwargs...)

define!(s::Symbol, args...) = define!(_manopt_glossary, s, args...)
function define!(g::_LIEGROUPS_DOC_TYPE, s::Symbol, e::Union{String,Function})
    g[s] = e
    return g
end
function define!(g::_LIEGROUPS_DOC_TYPE, s1::Symbol, s2::Symbol, args...)
    !(haskey(g, s1)) && (g[s1] = _LIEGROUPS_DOC_TYPE())
    define!(g[s1], s2, args...)
    return g
end

# ---
# LaTeX
# Define LaTeX shortcuts
_tex(args...; kwargs...) = glossary(:LaTeX, args...; kwargs...)
define!(:LaTeX, :bigl, raw"\bigl")
define!(:LaTeX, :bigr, raw"\bigr")
define!(:LaTeX, :Big, raw"\Big")
define!(:LaTeX, :Bigl, raw"\Bigl")
define!(:LaTeX, :Bigr, raw"\Bigr")
define!(:LaTeX, :Cal, (letter) -> raw"\mathcal " * "$letter")
define!(:LaTeX, :Frak, (letter) -> raw"\mathfrak " * "$letter")
define!(:LaTeX, :rm, (letter) -> raw"\mathrm " * "$letter")
#
# ---
# Mathematics and semantic symbols
# :symbol the symbol,
# :description the description
_math(args...; kwargs...) = glossary(:Math, args...; kwargs...)
define!(:Math, :GroupOp, :symbol, "∘")
define!(:Math, :GroupOp, :descrption, "the Lie Group operation")
define!(:Math, :op, _math(:GroupOp, :symbol))
define!(:Math, :e, _tex(:rm, "e"))
define!(:Math, :LieAlgebra, :symbol, (; g="g") -> _tex(:Frak, g))
define!(:Math, :LieAlgebra, :descrption, "the ie Algebra")
define!(:Math, :𝔤, (; G="G") -> _math(:LieAlgebra, :symbol; g="g"))
define!(:Math, :LieGroup, :symbol, (; G="G") -> _tex(:Cal, G))
define!(:Math, :LieGroup, :descrption, "the Lie Group")
define!(:Math, :G, (; G="G") -> _math(:LieGroup, :symbol; G=G))
define!(:Math, :Manifold, :symbol, (; M="M") -> _tex(:Cal, M))
define!(:Math, :Manifold, :descrption, "the Riemannian manifold")
define!(:Math, :M, (; M="M") -> _math(:Manifold, :symbol; M=M))

#
# ---
# Links
# Collect short forms for links, especially Interdocs ones.
_link(args...; kwargs...) = glossary(:Link, args...; kwargs...)

define!(
    :Link,
    :AbstractManifold,
    "[`AbstractManifold`](@extref `ManifoldsBase.AbstractManifold`)",
)
define!(:Link, :TangentSpace, "[`TangentSpace`](@extref `ManifoldsBase.TangentSpace`)")
