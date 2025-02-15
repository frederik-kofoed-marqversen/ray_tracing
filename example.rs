use std::collections::HashMap;

use super::triangle_mesh::TriangleMesh;
use super::{Point3D, Vec3D};

struct SDVertex {
    point: Vec3D,
    reference_face: usize,
}

#[derive(Hash, PartialEq, Eq, Clone)]
struct SDEdge {
    v1: usize,
    v2: usize,
}

impl SDEdge {
    fn new(v1: usize, v2: usize) -> Self {
        let mut vertices = [v1, v2];
        vertices.sort();
        return Self {
            v1: vertices[0],
            v2: vertices[1],
        };
    }
}

struct SDFace {
    vertices: [usize; 3],
    neighbours: [Option<usize>; 3],
}

impl SDFace {
    fn vertex_index(&self, vertex: usize) -> usize {
        if let Some(index) = self.vertices.iter().position(|other| other == &vertex) {
            return index;
        } else {
            // We should never ask for a vertex that is not part of a face
            panic!("Logic error!")
        }
    }
}

struct SDSurface {
    vertices: Vec<SDVertex>,
    faces: Vec<SDFace>,
}

fn next(i: usize) -> usize {
    (i + 1) % 3
}

fn prev(i: usize) -> usize {
    (i + 2) % 3
}

impl SDSurface {
    fn next_vertex(&self, vertex: usize, face: usize) -> usize {
        let face = &self.faces[face];
        face.vertices[next(face.vertex_index(vertex))]
    }

    fn prev_vertex(&self, vertex: usize, face: usize) -> usize {
        let face = &self.faces[face];
        face.vertices[prev(face.vertex_index(vertex))]
    }

    // Given a vertex `i` of a face the next face is defined as that which shares the edge
    // (`i`, `next(i)`) which is simply the face `neighbours[i]`.
    fn next_face(&self, vertex: usize, face: usize) -> Option<usize> {
        let face = &self.faces[face];
        face.neighbours[face.vertex_index(vertex)]
    }

    fn prev_face(&self, vertex: usize, face: usize) -> Option<usize> {
        let face = &self.faces[face];
        face.neighbours[prev(face.vertex_index(vertex))]
    }

    /// Computes the neighbouring vertices for a given vertex along with a bool
    /// which is true if the given vertex is a boundary vertex and false if not.
    fn vertex_neighbours(&self, vertex: usize) -> (Vec<usize>, bool) {
        let start_face = self.vertices[vertex].reference_face;
        let mut neighbours = Vec::new();

        // Loop through incident faces until return to start or hit boundary.
        // For each face, store the next vertex along that face.
        let mut face = start_face;
        neighbours.push(self.next_vertex(vertex, face));
        dbg!(vertex);
        while let Some(next_face) = self.next_face(vertex, face) {
            if next_face == start_face {
                // Has returned to starting point, and vertex is thus not a boundary point
                return (neighbours, false);
            } else {
                face = next_face;
                neighbours.push(self.next_vertex(vertex, next_face));
            }
        }
        
        // A boundary must have been hit.
        // Reverse the order of neighbours so that boundary vertices will be at either end
        neighbours.reverse();
        // Start again but go in opposite direction.
        face = start_face;
        neighbours.push(self.prev_vertex(vertex, face));
        while let Some(prev_face) = self.prev_face(vertex, face) {
            face = prev_face;
            neighbours.push(self.prev_vertex(vertex, prev_face));
        }

        return (neighbours, true);
    }
}