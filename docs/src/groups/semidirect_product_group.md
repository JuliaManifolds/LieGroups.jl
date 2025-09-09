# The semidirect product Lie group

The semidirect product has a few choices regarding left and right:

* for the order of the product: [`LeftSemidirectProductLieGroup`](@ref) ``\mathcal G ⋉ \mathcal H`` vs.  [`RightSemidirectProductLieGroup`](@ref) ``\mathcal H ⋊ \mathcal G``
* for the [`GroupAction`](@ref) ``α``: [`AbstractLeftGroupActionType`](@ref) ``σ`` vs. [`AbstractRightGroupActionType`](@ref) ``τ``
* for the [`GroupAction`](@ref) ``α_g`` w.r.t. a fixed ``g ∈ \mathcal G``, within a group operation ``h_1⋄h_2``: [`ActionActsOnLeft`](@ref) ``α_g(h_1)⋄h_2`` vs. [`ActionActsOnRight`](@ref) ``h_1⋄α_g(h_2)``

These choices lead to different formulae, usually even all eight cases are different. We still try to document them.

```@autodocs
Modules = [LieGroups]
Pages = ["groups/semidirect_product_group.jl"]
Order = [:type, :function]
```