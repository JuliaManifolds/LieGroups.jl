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
# to keep naming, notation, and formatting in a unified way

# In general every dictionary here can be either `:Symbol-> String` or `:Symbol -> Dictionary enrties`

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
#Do not attach the doc string here for now, since there is no internals section in the documentation yet
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
define!(:LaTeX, :big, raw"\big")
define!(:LaTeX, :bigl, raw"\bigl")
define!(:LaTeX, :bigr, raw"\bigr")
define!(:LaTeX, :Big, raw"\Big")
define!(:LaTeX, :Bigl, raw"\Bigl")
define!(:LaTeX, :Bigr, raw"\Bigr")
define!(:LaTeX, :def, raw"\coloneqq")
define!(:LaTeX, :Cal, (letter) -> raw"\mathcal " * "$letter")
define!(:LaTeX, :exp, raw"\exp")
define!(:LaTeX, :frac, (a, b) -> raw"\frac" * "{$a}{$b}")
define!(:LaTeX, :Frak, (letter) -> raw"\mathfrak " * "$letter")
define!(:LaTeX, :l, "") # lazy fallback for sets
define!(:LaTeX, :r, "") # lazy fallback for sets
define!(:LaTeX, :log, raw"\log")
define!(:LaTeX, :qquad, raw"\qquad")
define!(:LaTeX, :quad, raw"\quad")
define!(:LaTeX, :rm, (letter) -> raw"\mathrm" * "{$letter}")
define!(
    :LaTeX,
    :Set,
    (content, size="") ->
        _tex(Symbol("$(size)l")) *
        raw"\{ " *
        "$(content)" *
        _tex(Symbol("$(size)r")) *
        raw" \}",
)
define!(
    :LaTeX,
    :SetDef,
    (elem, cond, size="") ->
        _tex(:Set, elem * raw"\ " * _tex(Symbol("$(size)")) * raw"|\ " * "$(cond)", size),
)
define!(:LaTeX, :sum, raw"\sum")
define!(:LaTeX, :transp, raw"\mathrm{T}")
define!(:LaTeX, :text, (text) -> raw"\text" * "{$text}")

#
# ---
# Mathematics and semantic symbols
# :symbol the symbol,
# :description the description
_math(args...; kwargs...) = glossary(:Math, args...; kwargs...)
define!(:Math, :Adjoint, :symbol, raw"\mathrm{Ad}")
define!(:Math, :Adjoint, :descrption, "the adjoint operation")
define!(:Math, :Ad, _math(:Adjoint, :symbol))
define!(:Math, :GroupAction, :symbol, "⋅")
define!(:Math, :GroupAction, :descrption, "a Lie Group Action")
define!(:Math, :act, _math(:GroupAction, :symbol))
define!(:Math, :GroupOp, :symbol, "∘")
define!(:Math, :GroupOp, :descrption, "the Lie Group operation")
define!(:Math, :∘, _math(:GroupOp, :symbol))
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
# Collect short forms for links, especially `DocumenterInterlinks` ones.
_link(args...; kwargs...) = glossary(:Link, args...; kwargs...)

define!(
    :Link,
    :AbstractManifold,
    "[`AbstractManifold`](@extref `ManifoldsBase.AbstractManifold`)",
)
define!(
    :Link,
    :isapprox,
    "[`isapprox`](@extref `Base.isapprox-Tuple{AbstractManifold, Any, Any, Any}`)",
)
define!(:Link, :TangentSpace, "[`TangentSpace`](@extref `ManifoldsBase.TangentSpace`)")
define!(
    :Link,
    :zero_vector,
    "[`zero_vector`](@extref `ManifoldsBase.zero_vector-Tuple{AbstractManifold, Any}`)",
)

#
# ---
# Links
# Collect certain formulae, short texts or even admonitions
_note(args...; kwargs...) = glossary(:Note, args...; kwargs...)

define!(
    :Note,
    :LeftInverseActionIsRight,
    """
    ```math
    τ_g(τ_h(p))
    = σ^{$(_tex(:rm,"L"))}_{g^{-1}}(σ^{$(_tex(:rm,"L"))}_{h^{-1}}(p))
    = σ^{$(_tex(:rm,"L"))}_{g^{-1}h^{-1}}(p)
    = σ^{$(_tex(:rm,"L"))}_{(hg)^{-1}}(p)
    τ_{hg}(p).
    ```
    """,
)

define!(
    :Note,
    :RightInverseActionIsLeft,
    """
    ```math
    τ_g(τ_h(p))
    = σ^{$(_tex(:rm,"R"))}_{g^{-1}}(σ^{$(_tex(:rm,"R"))}_{h^{-1}}(p))
    = σ^{$(_tex(:rm,"R"))}_{h^{-1}g^{-1}}(p)
    = σ^{$(_tex(:rm,"R"))}_{(gh)^{-1}}(p)
    τ_{gh}(p).
    ```
    """,
)
