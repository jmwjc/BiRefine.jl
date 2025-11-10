// square.geo - Gmsh geometry for a square mesh

// Define square corner points (side length = 1)
lc = 0.1;
Point(1) = {0, 0, 0, lc};
Point(2) = {4, 0, 0, lc};
Point(3) = {4, 4, 0, lc};
Point(4) = {0, 4, 0, lc};

// Define lines between points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Create curve loop and surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Curve {1, 3} = 11;
Transfinite Curve {2, 4} = 11;
Transfinite Surface {1} Right;

// Define physical groups for the surface
Physical Surface("Ω") = {1};
Physical Curve("Γ¹") = {1};
Physical Curve("Γ²") = {2};
Physical Curve("Γ³") = {3};
Physical Curve("Γ⁴") = {4};

// Mesh settings (optional)
// SetOrder = 2;
// Mesh.Algorithm = 1; // Use frontal algorithm

// To generate the mesh, run: gmsh square.geo -2