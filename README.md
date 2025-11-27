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
    - **More unit tests**: Add more unit tests for e.g. bsdf.rs, lights.rs, 
    - **Surface trait**: Move sampling methods from Light trait to Surface trait.
    - **Light trait**: Remove the need of this trait, and move functionality to object/material, and simply let engine know which objects should be considered as lights by supplying a vec of references. This way lights simply become a part of the objects vector.

- **Upgrade Engine**
    - **Eta correction**: Use eta corrected beta value for Russian rulette.
    - **Area lights**: Add more types of area lights e.g. TriangleMesh, BVH
    - **Surface roughness**: Add roughness modsel like microfacet model
    - **Depth of field**: Include depth of field to camera struct
    - **Parametrised media**: Build medium class that takes parametrised properties and then builds some effecient structure of constant majorant segments from that.
    
- **Diffraction**
    - **RGB to spectrum**: Implement RGB to spectrum and spectrum to RGB methods.
    - **Wavelength sampling**: Sample wavelengths from a camara activation function for each ray instansiation.
    - **Wavelength collections**: Make engine take a collection of wavelengths and integrate those simultaneously using MIS. Singular wavelength integration can be used as an initial (ineffecient) implementation.

- **Additional Procedural Features**
    - **Mesh Topology Fixes**: Fix the mesh topology fixing algorithm to also allow non-connected surfaces.
    - **Vertex Normals and UV Mapping**: Add support for vertex normals and UV mapping to triangle meshes.
    - **Mesh Retriangulation/Optimization**: Implement retriangulation and optimization techniques for meshes with poor topology (Hoppe et al. w. subdivision surfaces has good references).
    - **Subdivision Surface Enhancements**: Incorporate internal crease and corner logic into the subdivision surface algorithm.

- **Accelerators**
    - **GPU**: Self explanatory
    - **Mesh to BVH**: After building BVH, optimise the vertex array of the mesh to match BVH for better data locality
    - **BVH**: Engine should build BVH of all objects and lights and materials in the scene to use while rendering.

- **GUI and Real Time Rendering**: Using e.g. [egui](https://docs.rs/egui) or [imgui](https://docs.rs/imgui)
