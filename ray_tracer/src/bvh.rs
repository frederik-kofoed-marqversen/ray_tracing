use super::primitives::AxisAlignedBoundingBox as AABB;
use super::traits::{Bounded, Surface};
use super::triangle_mesh::Triangle;
use super::{Ray, Vec3D};

impl Triangle {
    #[inline]
    fn centroid_axis(&self, axis: usize) -> f32 {
        self.vertices().iter().map(|p| p[axis]).sum::<f32>() / 3.0
    }

    #[inline]
    fn centroid(&self) -> Vec3D {
        self.vertices().iter().sum::<Vec3D>() / 3.0
    }
}

/// Nodes have a bounding box, a counter for the number of primitives it contains,
/// and an index. If `num_prim == 0` then the Node is not a leaf node and `index` is
/// the index of the first child node. The second child node always has index equal
/// to `index + 1`. If on the other hand `num_prim > 0`, then the node contains
/// primitives and is a leaf node. In that case `index` is the index of the first
/// primitive.
/// We use u32 instead of usize to ensure Node has size 32 bytes, which should help
/// with caching.
#[derive(Debug)]
pub struct Node {
    bounding_box: AABB,
    index: u32,
    num_prim: u32,
}

impl Node {
    fn is_leaf(&self) -> bool {
        self.num_prim > 0
    }
}

/// Bounding Volume Hierarchy contains a reference to a Vec of primitives (for now
/// only triangles), and a Vec of nodes. The `indices` is a list of indices into
/// `triangles` such that if a node points to some `index`, then the corresponding
/// triangle is `triangles[indices[index]]`.
pub struct BoundingVolumeHierarchy {
    triangles: Vec<Triangle>,
    indices: Vec<usize>,
    nodes: Vec<Node>,
}

impl BoundingVolumeHierarchy {
    pub fn build(triangles: Vec<Triangle>) -> Self {
        // Initialise BVH.
        // The number of nodes is upper bounded by 2n-1 since the maximal number of
        // leafs is n, and for each layer of m nodes, there are a maximum of 2m
        // children. Thus, the maximum number of nodes is 1 + 2 +... + n = 2n-1.
        let n = triangles.len();
        assert!(2 * n - 1 < u32::MAX as usize);
        let mut bvh = Self {
            triangles,
            indices: (0..n).collect(),
            nodes: Vec::with_capacity(2 * n - 1),
        };

        // Add root node
        let root_node = bvh.new_node(0, n as u32);
        bvh.nodes.push(root_node);
        // Recursively subdivide the root node
        bvh.subdivide_recursive(0);
        bvh.nodes.shrink_to_fit();

        return bvh;
    }

    #[inline]
    fn get_triangles(&self, first_prim: u32, num_prim: u32) -> impl Iterator<Item = &Triangle> {
        (first_prim..first_prim + num_prim).map(|i| &self.triangles[self.indices[i as usize]])
    }

    #[inline]
    fn new_node(&self, first_prim: u32, num_prim: u32) -> Node {
        let bounding_box = self
            .get_triangles(first_prim, num_prim)
            .map(|tri| tri.bounding_box())
            .fold(AABB::empty(), |a, b| AABB::combine(&a, &b));

        return Node {
            bounding_box,
            index: first_prim,
            num_prim,
        };
    }

    fn subdivide_recursive(&mut self, node_index: usize) {
        // Compute optimal partition plane
        let (axis, split_point, cost) = self.binned_sah_partition_plane(node_index);

        // Get node data
        let node = &self.nodes[node_index];
        let first_prim = node.index;
        let num_prim = node.num_prim;

        // Only do partitioning if it improves the BVH
        let current_cost = num_prim as f32 * node.bounding_box.half_area();
        if cost >= current_cost {
            return;
        }

        // Perform partitioning
        let partition_point = self.partition_sort(first_prim, num_prim, axis, split_point);

        // Catch failed partitioning
        let left_count = partition_point - first_prim;
        let right_count = num_prim - left_count;
        if left_count == 0 || right_count == 0 {
            // One of the partitions are empty and
            // so the node cannot be subdivided.
            return;
        };

        // Update current node since it is now being subdivided
        let left_child_index = self.nodes.len();
        let node = &mut self.nodes[node_index];
        node.index = left_child_index as u32;
        node.num_prim = 0;

        // Spawn children
        let left_node = self.new_node(first_prim, left_count);
        let right_node = self.new_node(partition_point, right_count);
        self.nodes.push(left_node);
        self.nodes.push(right_node);

        // Recursively subdivide child nodes
        self.subdivide_recursive(left_child_index);
        self.subdivide_recursive(left_child_index + 1);
    }

    /// Compute the optimal partition plane using binned Surface Area Heuristic
    fn binned_sah_partition_plane(&mut self, node_index: usize) -> (usize, f32, f32) {
        // Hardcoded number of bins
        const BINS: usize = 100;

        let node = &self.nodes[node_index];
        let first_prim = node.index;
        let num_prim = node.num_prim;

        // Construct centroid bounding box
        let centroid_aabb = self
            .get_triangles(first_prim, num_prim)
            .fold(AABB::empty(), |res, tri| res.grow(tri.centroid()));

        // Determine optimal partition plane for each axis using binned Surface Area Heuristic
        let mut axis = 0;
        let mut split_point = 0.0;
        let mut cost = f32::INFINITY;
        for plane_axis in 0..3 {
            let max = centroid_aabb.upper[plane_axis];
            let min = centroid_aabb.lower[plane_axis];
            if (max - min).abs() < f32::EPSILON {
                continue;
            }
            let scale = BINS as f32 / (max - min);

            // Build bins: (AABB, tri_count)
            let mut bins = vec![(AABB::empty(), 0); BINS];
            for triangle in self.get_triangles(first_prim, num_prim) {
                let bin_index = ((triangle.centroid_axis(plane_axis) - min) * scale) as usize;
                let bin_index = bin_index.min(BINS - 1);
                let bin = &mut bins[bin_index];
                bin.0 = AABB::combine(&bin.0, &triangle.bounding_box());
                bin.1 += 1;
            }

            // Initialise data Vecs
            let mut left_area = vec![0.0; BINS - 1];
            let mut right_area = vec![0.0; BINS - 1];
            let mut left_count = vec![0; BINS - 1];
            let mut right_count = vec![0; BINS - 1];

            // Compute bin data
            let mut left_box = AABB::empty();
            let mut right_box = AABB::empty();
            let mut left_sum = 0;
            let mut right_sum = 0;
            for i in 0..BINS - 1 {
                left_sum += bins[i].1;
                left_count[i] = left_sum;
                left_box = AABB::combine(&left_box, &bins[i].0);
                left_area[i] = left_box.half_area();

                right_sum += bins[BINS - 1 - i].1;
                right_count[BINS - 2 - i] = right_sum;
                right_box = AABB::combine(&right_box, &bins[BINS - 1 - i].0);
                right_area[BINS - 2 - i] = right_box.half_area();
            }

            // Determine best partition plane
            let scale = (max - min) / BINS as f32;
            for i in 0..BINS - 1 {
                let plane_cost =
                    left_count[i] as f32 * left_area[i] + right_count[i] as f32 * right_area[i];
                if plane_cost < cost {
                    axis = plane_axis;
                    split_point = min + scale * (i + 1) as f32;
                    cost = plane_cost;
                }
            }
        }

        return (axis, split_point, cost);
    }

    fn partition_sort(
        &mut self,
        first_prim: u32,
        num_prim: u32,
        axis: usize,
        split_point: f32,
    ) -> u32 {
        // Move triangles contained in this node to match the partitioning
        let mut i = first_prim;
        let mut j = first_prim + num_prim - 1;
        while i <= j {
            if self.triangles[self.indices[i as usize]].centroid_axis(axis) < split_point {
                i += 1;
            } else {
                self.indices.swap(i as usize, j as usize);
                j -= 1;
            }
        }
        // Return the partition point (beginning of second partition)
        return i;
    }
}

impl Bounded for BoundingVolumeHierarchy {
    fn bounding_box(&self) -> AABB {
        self.nodes[0].bounding_box.clone()
    }
}

impl Surface for BoundingVolumeHierarchy {
    fn hit(&self, ray: &Ray, t_min: f32, mut t_max: f32) -> Option<(f32, Vec3D)> {
        let mut hit = None;

        let mut queue: Vec<&Node> = vec![&self.nodes[0]];
        while let Some(node) = queue.pop() {
            if node.is_leaf() {
                // Intersect triangles contained in node.
                for triangle in self.get_triangles(node.index, node.num_prim) {
                    if let Some(new_hit) = triangle.hit(ray, t_min, t_max) {
                        t_max = new_hit.0;
                        hit = Some(new_hit);
                    }
                }
                continue;
            }

            // Intersect ray with both children
            let mut child1 = &self.nodes[node.index as usize];
            let mut child2 = &self.nodes[node.index as usize + 1];
            let hit1 = child1.bounding_box.hit(ray, t_min, t_max);
            let hit2 = child2.bounding_box.hit(ray, t_min, t_max);

            // Check if either child is hit and add them to the queue.
            match (hit1, hit2) {
                (None, None) => continue,
                (Some(_), None) => {
                    queue.push(child1);
                    continue;
                }
                (None, Some(_)) => {
                    queue.push(child2);
                    continue;
                }
                (Some(t1), Some(t2)) => {
                    // Both bounding boxes were hit => Add both to the queue, the closest
                    // is added last so that it is checked first.
                    if t1 > t2 {
                        std::mem::swap(&mut child1, &mut child2);
                    }
                    queue.push(child2);
                    queue.push(child1);
                }
            };
        }

        return hit;
    }
}
