# ray_tracing
Ray tracing graphics rendering project in Rust.

TODO:
Try rendering pixel tiles (e.g. 4x4) instead of lines for better temporal data locality
Mesh fix topology should handle non-connected surfaces somehow
Retriangulation/uptimisation (of bad-topology meshes) (Hoppe et al. w. subdivision surfaces has good references)
Add vertex normals and uv-mapping to triangle mesh
Add internal crease and corner logic to subdivision surface alg.