# Ray Tracing in Rust

This repository is a personal project focused on developing a ray tracing graphics renderer using the Rust programming language. Initiated out of personal interest, the project has been developed intermittently as a hobby. Currently, the renderer operates solely on the CPU, but I do plan on exploring GPU acceleration in the future.

## Rendering an Image

Currently, to render an image, simply run the `main.rs` file and direct standard output to a `.ppm` file. From within the binary (not module) folder, run the following command in the terminal:

```bash
cargo run --release > image.ppm
```

## TODO

- **Image post processing**: Adding low-pass filters. Maybe use Mitchell–Netravali?

- **Additional Material Properties**: Such as [Frenel reflectance](https://pbr-book.org/3ed-2018/Reflection_Models/Fresnel_Incidence_Effects) for non-dielectrics and also some [Microfacet models](https://pbr-book.org/3ed-2018/Reflection_Models/Microfacet_Models).

- **GUI and Real Time Rendering**: Using e.g. [egui](https://docs.rs/egui) or [imgui](https://docs.rs/imgui)

- **Pixel Tiles**: Implement rendering of pixel tiles (e.g., 4x4) instead of lines to improve temporal data locality.

- **Mesh Retriangulation/Optimization**: Implement retriangulation and optimization techniques for meshes with poor topology (Hoppe et al. w. subdivision surfaces has good references).

- **Vertex Normals and UV Mapping**: Add support for vertex normals and UV mapping to triangle meshes.

- **Mesh Topology Fixes**: Fix the mesh topology fixing algorithm to also allow non-connected surfaces.

- **Subdivision Surface Enhancements**: Incorporate internal crease and corner logic into the subdivision surface algorithm.

