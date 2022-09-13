# g-poly
Polytopal complexes and their g- (or h-) polynomials

# Overview
This code was written as a part of a project with Marc Besson, Samuel Jeralds, and Joshua Kiers investigating the connection between Lie theory and certain toric varieties.The key features are:
1) simple implementation of regular polytopal complexes,
2) cones and barycentric subdivisions,
3) algorithm for generating type A permutohedra,
4) recursive algorithm for computing h-polynomials.

# Building a complex
The building blocks of polytopal complexes are convex polytopes. The interior of a convex polytope is called a "cell". Cells are constructed recursively as Cell objects. The constructor for Cell takes two arguments, a set of boundary facets (codimension 1 faces), and a dimension. For instance, a vertex a is constructed by a = Cell(set(), 0). Alternatively, one can simply call the vertex() method: a = vertex(). Cells are compared at the object level, so to create four distinct vertices, call vertex() four times.

a = vertex()
b = vertex()
c = vertex()
d = vertex()

Once the vertices are created, edges can be created by specifying their two boundary vertices. The following would create edges between a and b, b and c, c and d, a and d, and a and c.

ab = Cell({a, b}, 1)
bc = Cell({b, c}, 1)
cd = Cell({c, d}, 1)
ad = Cell({a, d}, 1)
ac = Cell({a, c}, 1)

Once the edges have been created, 2-cells can be created by specifying the codimension 1 cells on the boundary and so on for higher dimensional cells. The following code would create two triangular cells.

abc = Cell({ab, bc, ac}, 2)
acd = Cell({ac, cd, ad}, 2)

The Cell object retains information about its boundary. One can use Cell.get_faces(d) to get a set of all d-dimensional cells in the cell closure. Cell.get_all_faces() returns a set of all cells in the cell closure.

ab.getfaces(0) = {a, b}
abc.getfaces(1) = {ab, bc, ac}
acd.get_all_faces() = {acd, ac, cd, ad, a, c, d}

Once the desired cells are built, one can assemble them into a complex using the Complex object. Its constructor takes a set of maximal cells. The complex

C = Complex({abc, acd})

is a square cut in half along a diagonal.


As an alternative to manually constructing cells and complexes this way, the complex_from_hulls() method builds a complex from a combinatorial encoding.

C = complex_from_hulls([[{'a'}, {'b'}, {'c'}], [{'a', 'b'}, {'b', 'c'}, {'c', 'd'}, {'a', 'd'}, {'a', 'c'}], [{'a', 'b', 'c'}, {'a', 'c', 'd'}]])

results in the same complex as the manual approach above. We represent cells as convex hulls of vertices. A convex hull is given by a set of strings, where each string represents a vertex. These sets are grouped into lists according to dimension of the hull. The input is the list of all of these lists in ascending order of dimension.

# Calculating invariants and other operations on complexes

Calling Complex.f_vec(), Complex.h_poly(), and Complex.g_poly() returns a sympy expression for the f-vector, h-polynomial, and g_polynomial respectively. One can also take the cone over vertex v by calling cone(C, v) and the barycentric subdivision of C by calling bary_sub(C). These will be Complex objects as well. The product() method takes two complexes and returns their Cartesian product. Finally, permutohedron(n) returns the type A_n permutohedron as a Complex object.
