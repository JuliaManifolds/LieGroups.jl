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
compose
compose!
diff_left_compose
diff_left_compose!
diff_right_compose
diff_right_compose!
inv_left_compose
inv_left_compose!
inv_right_compose
inv_right_compose!
conjugate
conjugate!
diff_conjugate
diff_conjugate!
exp
exp!
identity_element
identity_element!
is_identity
inv(::LieGroup, ::Any)
inv!
diff_inv
diff_inv!
log
log!
```


## Functions on Lie Algebras

```@docs
lie_bracket
lie_bracket!
```

## Literature

```@bibliography
Pages = ["group.md"]
Canonical=false
```