# The interface for Lie Groups and Lie algebras

```@docs
LieGroup
LieAlgebra
```

## Functions on Lie Groups

```@docs
adjoint
adjoint!
base_manifold
conjugate
conjugate!
diff_conjugate
diff_conjugate!
exp
exp!
identity_element
identity_element!
is_identity
isapprox
is_point(::LieGroup, ::Any)
is_vector(::LieGroup, ::Any)
is_vector(::LieGroup{𝔽,O}, ::Identity{O}, ::Any) where {𝔽,O}
log
log!
```


## Functions on Lie algebras

```@docs
is_point(::LieAlgebra, ::Any)
lie_bracket
lie_bracket!
```

## Literature

```@bibliography
Pages = ["group.md"]
Canonical=false
```