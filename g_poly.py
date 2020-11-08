# g-Polynomial
from sympy import symbols, expand, poly, init_printing
from math import floor
from IPython.display import display
from itertools import permutations, combinations

#general machinery for polytopal complexes

#parameter for h,g-polynomials
t = symbols('t')

#Polytopal complex with specified set of maximal cells
class Complex:
    def __init__(self, max_cells):
        self.max_cells = max_cells
        self.dim = max([cell.dim for cell in max_cells])
        self.f = tuple()
        self.g = -1
        self.h = -1
        
        self.cells = [set() for i in range(self.dim + 1)]
        for cell in max_cells:
            for i in range(cell.dim + 1):
                self.cells[i] |= cell.get_faces(i)
        return
    
    #Set of cells of dimension dim
    def get_cells(self, dim):
        return self.cells[dim]
    
    #Set of all cells
    def get_all_cells(self):
        all_cells = set()
        for s in self.cells:
            all_cells |= s
        return all_cells
    
    #Set of maximal cells
    def get_max_cells(self):
        return self.max_cells
    
    #List of the names of all cells organized according to dimension
    def display_cells(self):
        for i in range(self.dim + 1):
            print("***** " + str(i) + "-cells *****")
            for cell in self.get_cells(i):
                print(cell.name)
        return
    
    #Returns True if all maximal cells are simplices
    def is_simplicial(self):
        for cell in self.max_cells:
            if not cell.is_simplex():
                return False
        return True
    
    #Removed non-maximal cells from max_cells. Must be called if non-maximal
    #cells were passed to constructor
    def clean(self):
        red_cells = set()
        for mc in self.max_cells:
            for c in self.get_all_cells():
                if (mc.get_all_faces() < c.get_all_faces()):
                    red_cells.add(mc)
        self.max_cells -= red_cells
        return self
    
    #Returns the f-vector of the complex as a tuple of ints   
    def f_vec(self):
        if self.f == tuple():
            self.f = tuple([len(self.get_cells(i)) for i in range(self.dim + 1)])
        return self.f
    
    #Returns the h-polynomial as a sympy expression in t
    def h_poly(self):
        if self.h == -1:
            self.h = (t - 1)**self.dim
            facet_comp = Complex(self.get_cells(self.dim - 1))
            if facet_comp.is_simplicial():
                for i in range(self.dim):
                    self.h += self.f_vec()[i]*(t - 1)**(self.dim - i - 1)
            else:
                for cell in facet_comp.get_all_cells():
                    self.h += Complex({cell}).g_poly() * (t - 1)**(self.dim 
                                    - 1 - cell.dim)
        return self.h
    
    #Returns the g-polynomial as a sympy expression in t
    def g_poly(self):
        coeffs_h = poly(self.h_poly(), t).all_coeffs()
        if self.g == -1:
            self.g = coeffs_h[0]
            for i in range(floor(self.dim / 2)):
                self.g += (coeffs_h[i + 1] - coeffs_h[i]) * t**(i + 1)
        return self.g
    
#Polytope with specified dimension and set of facets. The set of facets of a
#vertex is empty. Higher dimensional polytopes are built inductively. The name
#attribute is for labeling only. Cells are compared at the object level, so
#distinct cells may have the same name.
class Cell:
    def __init__(self, dim, facets, name = ""):
        self.dim = dim
        self.name = name
        
        self.faces = [{self}]
        for i in range(dim):
            self.faces.insert(0, set())
            
        for facet in facets:
            for i in range(dim):
                self.faces[i] |= facet.get_faces(i)
        return

    #Set of faces of dimension dim                    
    def get_faces(self, dim):
        return self.faces[dim]
    
    #Set of all faces
    def get_all_faces(self):
        all_faces = set()
        for s in self.faces:
            all_faces |= s
        return all_faces
    
    #Set of codimension 1 faces
    def get_facets(self):
        if self.dim == 0:
            return set()
        return self.faces[self.dim - 1]
    
    #Returns true if polytope is a simplex
    def is_simplex(self):
        if len(self.faces[0]) == self.dim + 1:
            return True
        else:
            return False
    
#Constructs a vertex with specified name
def vertex(name = ""):
    return Cell(0, set(), name)

#Forms the cone of complex comp with cone point v. If no v is passed, a new
#vertex is created to serve as the cone point. The optional cones argument is
#a dictionary of cones of complexes, {(cell, vertex):cone of cell over v}).
#This is useful for merging multiple cones together.
def cone(comp, v = None, cones = {}):
    new_cells = set()
    if v == None:
        v = vertex()
    for vert in comp.get_cells(0):
        if (vert, v) not in cones:
            cones[(vert, v)] = Cell(1, {vert, v}, v.name + 'C' + vert.name)
        new_cells.add(cones[(vert, v)])
    for i in range(1, comp.dim + 1):
        for cell in comp.get_cells(i):
            if (cell, v) not in cones:
                cones[(cell, v)] = Cell(i + 1, {cones[(facet, v)] for facet 
                                         in cell.get_facets()} | {cell}, 
                                        v.name + 'C' + cell.name)
            new_cells.add(cones[(cell, v)])
    return Complex(new_cells).clean()
    
def bary_sub(comp, level = 0, B = {}, C = {}):
    if comp.dim <= level:
        for max_cell in comp.get_max_cells():
            if max_cell not in B:
                B[max_cell] = vertex('B' + max_cell.name)
        return comp
    else:
        new_cells = set()
        for max_cell in comp.get_max_cells():
            if max_cell not in B:
                B[max_cell] = vertex('B' + max_cell.name)
            facets = Complex(max_cell.get_facets())
            new_cells |= cone(bary_sub(facets, level, B, C), 
                              B[max_cell], C).get_max_cells()
        return Complex(new_cells)

#Barycentric subdivision of the boundary of a complex (up to the level-
#skeleton. When the complex consists of a single polytope P, refine(P) is the 
#polytope corresponding to the resolution of the toric variety associated to P.
def refine(comp, level = 0):
    facet_comp = Complex(comp.get_cells(comp.dim - 1))
    subdiv = bary_sub(facet_comp, level)
    open_cell = Cell(comp.dim, subdiv.get_max_cells(), 
                     comp.get_max_cells().copy().pop().name)
    return Complex({open_cell})

#Returns a tuple of step-by-step refinements of comp. The first element is the
#full refinement, and the last is comp. (Naive implementation, slow)
def refine_in_steps(comp):
    return tuple([refine(comp, level) for level in range(comp.dim)])

#Builds a complex from the convex hulls of vertices. Input is a list,
#[0-cells, 1-cells,..., d-cells], where k-cells is a list of sets of vertices.
#A set of vertices corresponds to the convex hull of those vertices. A vertex 
#can be defined as a set containing a single object of any type.
def complex_from_hulls(skeleta):
    cells = [set() for i in range(len(skeleta))]
    cells_dict = {}
    for s in skeleta[0]:
        cells_dict[frozenset(s)] = vertex(str(s))
        cells[0].add(cells_dict.get(frozenset(s)))
        
    for i in range(1, len(skeleta)):
        for hull in skeleta[i]:
            facets = set()
            for f in skeleta[i - 1]:
                if f <= hull:
                    facets.add(cells_dict.get(frozenset(f)))
            cells_dict[frozenset(hull)] = Cell(i, facets, str(hull))
            cells[i].add(cells_dict.get(frozenset(hull)))
    return Complex(cells[len(skeleta) - 1])

#Cartesian product of two complexes
def product(comp1, comp2):
    D = {}
    d = comp1.dim + comp2.dim
    pairs = [set() for i in range(d + 1)]
    
    for i in range(d + 1):
        for j in range(min(i + 1, comp1.dim + 1)):
            if i - j <= comp2.dim:
                pairs[i] |= set([(cell1, cell2) 
                                for cell1 in comp1.get_cells(j)
                                for cell2 in comp2.get_cells(i - j)])
    for pair in pairs[0]:
        D[pair] = vertex(pair[0].name + ';' + pair[1].name)
    for i in range(1, d+1):
        for cell1, cell2 in pairs[i]:
            facets = set()
            for facet in cell1.get_facets():
                facets.add(D[(facet, cell2)])
            for facet in cell2.get_facets():
                facets.add(D[(cell1, facet)])
            D[(cell1, cell2)] = Cell(i, facets, 
                                      cell1.name + ';' + cell2.name)
    return Complex(set(D.values())).clean()

#Unit n-cube as a complex
def cube(n):
    if n == 0:
        return Complex({vertex("0")})
    else:
        L = complex_from_hulls([[{0}, {1}], [{0, 1}]])
    if n == 1:
        return L
    return product(cube(n - 1), L)


# type A permutohedron

#Returns a set of weak orderings of 1,...,n with classes many equivalence
#classes. Weak orerings are encoded as tuples of sets of ints. The sets are the
#equivalence classes of the reflexive closure of the relation given by
#incomparability.
def weak_orderings(n, classes):
    perms = set(permutations(range(1, n + 1)))
    if classes == 1:
        return {(frozenset(range(1, n + 1)),)}
    
    bars = set(combinations(range(1, n), classes - 1))
    ords = set()
    for p in perms:
        for b in bars:
            ords.add((frozenset(p[:b[0]]),)
                + tuple([frozenset(p[b[i - 1]:b[i]]) 
                            for i in range(1, classes - 1)])
                + (frozenset(p[b[classes - 2]:]),))
    return ords

#Returns true if string of ints s is ordered with respect to weak ordering w    
def is_ordered(s, w):
    for i in range(len(w) - 1):
        for j in range(i + 1, len(w)):
            for c1 in w[i]:
                for c2 in w[j]:
                    if s[c1 - 1] > s[c2 - 1]:
                        return False
    return True

#String representation of weak ordering w
def str_weak_ordering(w):
    s = ""
    for c in w[0]:
        s += str(c)
    if len(w) == 1:
        return s
    s += "|"
    for i in range(1, len(w) - 1):
        for c in w[i]:
            s += str(c)
        s += "|"
    for c in w[-1]:
            s += str(c)
    return s

#A_n permutohedron as a complex
def permutohedron(n):
    cell = {}
    facets = {}
    for w in weak_orderings(n + 1, n + 1):
        cell[w] = vertex(str_weak_ordering(w))
    for d in range(1, n + 1):
        for w in weak_orderings(n + 1, n + 2 - d):
            for i in range(len(w) - 1):
                join = list(range(len(w) - 1))
                if i > 0:
                    join[:i] = w[:i]
                join[i] = w[i] | w[i + 1]
                join[i + 1:] = w[i + 2:]
                
                join = tuple(join)
                if join in facets:
                    facets[join].add(cell[w])
                else:
                    facets[join] = {cell[w]}
                    
        for w in weak_orderings(n + 1, n + 1 - d):
            cell[w] = Cell(d, facets[w], str_weak_ordering(w))
                
    return Complex({cell[(frozenset(range(1, n + 2)),)]})

#---------------------code ends here----------------------------


init_printing()

# 1-simplex
L = complex_from_hulls([[{'a'}, {'b'}], [{'a', 'b'}]])

# 2-simplex
T2 = complex_from_hulls([[{'a'}, {'b'}, {'c'}], 
                            [{'a', 'b'}, {'b', 'c'}, {'a', 'c'}], 
                            [{'a', 'b', 'c'}]])

# square
C2 = complex_from_hulls([[{'a'}, {'b'}, {'c'}, {'d'}], [{'a', 'b'}, {'b', 'c'}, 
                                                        {'c', 'd'}, {'a', 'd'}],
                            [{'a', 'b', 'c', 'd'}]])

# 3-simplex
T3 = cone(T2)    

# 3-cube
C3 = complex_from_hulls([[{'a1'}, {'b1'}, {'c1'}, {'d1'}, 
                             {'a2'}, {'b2'}, {'c2'}, {'d2'}], 
                            [{'a1', 'b1',}, {'b1', 'c1'}, {'c1', 'd1'},
                             {'d1', 'a1'}, {'a2', 'b2',}, {'b2', 'c2'}, 
                             {'c2', 'd2'}, {'d2', 'a2'}, {'a1', 'a2'}, 
                             {'b1', 'b2'}, {'c1', 'c2'}, {'d1', 'd2'}], 
                            [{'a1', 'b1', 'c1', 'd1'}, {'a2', 'b2', 'c2', 'd2'}, 
                             {'a1', 'b1', 'a2', 'b2'}, {'d1', 'c1', 'd2', 'c2'}, 
                             {'c1', 'b1', 'c2', 'b2'}, {'a1', 'd1', 'a2', 'd2'}], 
                            [{'a1', 'b1', 'c1', 'd1', 'a2', 'b2', 'c2', 'd2'}]])

# 4-cube
C4 = product(C3, L)

# 5-cube
C5 = product(C4, L)

# A2 permutohedron
A2 = complex_from_hulls([[{1}, {2}, {3}, {4}, {5}, {6}], 
                         [{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {1, 6}], 
                         [{1, 2, 3, 4, 5, 6}]])

# square with diagonal
S = complex_from_hulls([[{1}, {2}, {3}, {4}], [{1,2}, {2,3}, {3,4}, {1,4},
                        {1,3}], [{1, 2, 3}, {1, 4, 3}]])

    
    
# h-polynomial tests
n = 4
C = permutohedron(n)

#step-by-step subdivision of A_n permutohedron

#print("\nf-vectors and h-polynomials of refinements of A_" + str(n) 
#        + " permutohedron")
#print("------------------------------")
#for i in range(n):
#    R = refine(C, n - 1 - i)
#    print(R.f_vec())
#    display(expand(R.h_poly()))
#    print("------------------------------")
    
#h-polynomials of facets of A_n permutohedron
    
#print("\nh-polynomials of facets of A_" + str(n) + " permutohedron")
#print("------------------------------")
#for i in range(floor((n - 1) / 2) + 1):
#    l = product(permutohedron(i), permutohedron(n - 1 - i))
#    m = [Complex({facet}).f_vec() 
#            for facet in C.get_cells(n - 1)].count(l.f_vec())
#    print("A_" + str(i) + " x " + "A_" + str(n - 1 - i) + " (" + str(m) 
#            + " many)")
#    display(expand(l.h_poly()))
#    print("------------------------------")

R2C = refine(C, 2)
print(expand(R2C.h_poly())


    
#A4 = 0
#A31 = 0
#A22 = 0
#for cell in C.get_cells(4):
#    if len(cell.get_faces(0)) == permutohedron(4).f_vec()[0]:
#        A4 += 1
#    elif len(cell.get_faces(0)) == product(permutohedron(3), 
#                                     permutohedron(1)).f_vec()[0]:
#        A31 += 1
#    elif len(cell.get_faces(0)) == product(permutohedron(2), 
#                                 permutohedron(2)).f_vec()[0]:
#        A22 += 1
#    else:
#        print("oops")
#print(str(A4) + " many A4")
#print(str(A31) + " many A3xA1")
#print(str(A22) + " many A2xA2")



#print("A2 permutohedron and its refinement")
#A2 = permutohedron(2)
#print("\n f = " + str(A2.f_vec()))
#display(expand(A2.h_poly()))
#
#RA2 = refine(A2)
#print("\n f = " + str(RA2.f_vec()))
#display(expand(RA2.h_poly()))
#
#print("\n\n")
#    
#print("A3 permutohedron and its refinement")
#A3 = permutohedron(3)
#print("\n f = " + str(A3.f_vec()))
#display(expand(A3.h_poly()))
#
#RA3 = refine(A3)
#print("\n f = " + str(RA3.f_vec()))
#display(expand(RA3.h_poly()))
#
#print("\n\n")
#
#print("A4 permutohedron and its refinement")
#A4 = permutohedron(4)
#print("\n f = " + str(A4.f_vec()))
#display(expand(A4.h_poly()))
#
#RA4 = refine(A4)
#print("\n f = " + str(RA4.f_vec()))
#display(expand(RA4.h_poly()))
#
#print("\n\n")
#
#print("A5 permutohedron and its refinement")
#A5 = permutohedron(5)
#print("\n f = " + str(A5.f_vec()))
#display(expand(A5.h_poly()))
#
#RA5 = refine(A5)
#print("\n f = " + str(RA5.f_vec()))
#display(expand(RA5.h_poly()))
#
#print("\n\n")