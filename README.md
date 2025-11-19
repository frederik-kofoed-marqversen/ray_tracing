# Ray Tracing in Rust

This repository is a personal project focused on developing a ray tracing graphics renderer using the Rust programming language. Initiated out of personal interest, the project has been developed intermittently as a hobby. Currently, the renderer operates solely on the CPU, but I do plan on exploring GPU acceleration in the future.

## Rendering an Image

Currently, to render an image, simply run the `main.rs` file and direct standard output to a `.ppm` file. From within the binary (not module) folder, run the following command in the terminal:

```bash
cargo run --release > image.ppm
```

## TODO

- **Image Post Processing**: Adding low-pass filters. Maybe use Mitchell–Netravali? Maybe use 2D gaussian splats method?

- **Cleanup**
    - **fn hit_bool**: Combine hit and hit_bool into a single function to reduce code duplication
    - **lib.rs**: Reconsider the organisation of lib.rs, traits.rs, utils.rs, bsdf.rs...
    - **Surface trait**: Move sampling methods from Light trait to Surface trait.
    - **Light trait**: Remove the need of this trait, and move functionality to object/material, and simply let engine know which objects should be considered as lights by supplying a vec of references. This way lights simply become a part of the objects vector.

- **Upgrade Engine**
    - **Stratified sampling**: As a start, use a low-discrepancy sequence for jittering on each pixel.
    - **Pixel Tiles**: Rendering of pixel tiles (e.g., 4x4) instead of lines could improve temporal data locality.
    - **BVH**: Once BVH struct is upgraded, engine should build BVH of all objects and lights in the scene before rendering.
    - **Light propagation**: Include history of refraction indices somehow plus include exponential absorbtion of light through dielectric.

- **Additional Scene Features**
    - **Additional Material Properties**: Such as [Frenel reflectance](https://pbr-book.org/3ed-2018/Reflection_Models/Fresnel_Incidence_Effects) for non-dielectrics and also some [Microfacet models](https://pbr-book.org/3ed-2018/Reflection_Models/Microfacet_Models).

- **Additional Procedural Features**
    - **Mesh Topology Fixes**: Fix the mesh topology fixing algorithm to also allow non-connected surfaces.
    - **Vertex Normals and UV Mapping**: Add support for vertex normals and UV mapping to triangle meshes.
    - **Mesh Retriangulation/Optimization**: Implement retriangulation and optimization techniques for meshes with poor topology (Hoppe et al. w. subdivision surfaces has good references).
    - **Subdivision Surface Enhancements**: Incorporate internal crease and corner logic into the subdivision surface algorithm.

- **Accelerators**
    - **GPU**: Self explanatory
    - **General BVH**: Make BVH take a Vec of some general objects instead of only Triangle.

- **GUI and Real Time Rendering**: Using e.g. [egui](https://docs.rs/egui) or [imgui](https://docs.rs/imgui)
