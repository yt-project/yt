Todo
----

Known Issues:

* ~~FRB Off-axis projections are broken I think. Currently should raise not-implemented error.~~
* Parallelism
  * Need to write parallel z-buffer reduce.
  * Need to verify brick ordering
* Alpha blending level for opaque sources such as grid lines/domains/etc may
  not currently be ideal. Difficult to get it right when the transparent VRs
  have wildly different levels. One approach would be to normalize the transfer
  function such that the integral of the TF multiplied by the depth of the 
  rendering is equal to 1. With grey opacity on, all of these things get a bit
  easier, in my opinion

Documentation:

* ~~Scene~~
* ~~Camera~~
* Lens
* Narrative
  * Have started, but more work to do. Replaced at least the tutorial
    rendering, which saves a number of lines!
* Cookbooks
  * All relevant cookbooks have been updated
* Parallelism
* OpaqueSource
* RenderSource
* Narrative Developer Documentation
