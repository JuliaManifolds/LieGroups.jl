
# Lie group operations

```@docs
AbstractGroupOperation
Identity
compose
compose!
diff_inv
diff_inv!
diff_left_compose
diff_left_compose!
diff_right_compose
diff_right_compose!
inv_left_compose
inv_left_compose!
inv_right_compose
inv_right_compose!
inv(::LieGroup, ::Any)
inv!
```

For some generic group operations, an easy generic implementation can be provided. This section lists these specific cases and what is implemented for these.

## Additive Group operation

```@docs
AdditionGroupOperation
```
