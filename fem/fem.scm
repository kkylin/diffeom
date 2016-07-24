#|

This file is part of DIFFEOM, a system for solving
differential equations on manifolds.

Copyright (C) 2016 by Kevin K Lin
<kkylin@alum.mit.edu>

This program is free software; you can redistribute
it and/or modify it under the terms of the GNU
General Public License as published by the Free
Software Foundation; either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General
Public License along with this program; if not, write
to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301 USA.

|#

;;; This program solves linear partial differential equations (with variable
;;; coefficients) using the finite element method (FEM).  It takes as input the
;;; source function (the RHS of the equation) and a list of nodes to be
;;; used. The boundary nodes should be initialized to the desired values.


;;; Internal dependencies among data structures:

;;; * In the final output matrix, each row corresponds to a node.  Hence, we
;;;   use nodes to organize the computation of coefficients in each row, and
;;;   the actual integrals are computed by the elements.  As a result, the
;;;   procedures ASSEMBLY and MAKE-NODE are fairly general, and one should
;;;   seldom need to change them.

;;; * The differential operator is encapsulated in the constructor for
;;;   elements, since the operator only affects the computation of the
;;;   integrals.  Many properties of the system (such as the differential
;;;   operator, element shapes/sizes, dimension of the domain) are parametrized
;;;   through this.

;;; Mathematical assumptions:

;;; * Linearity of the differential operator appears to be necessary for these
;;;   methods.  Otherwise nodal assembly wouldn't work, and we wouldn't be able
;;;   to construct linear equations out of the matrix of inner products.

;;; * We cannot handle systems of equations yet. 

;;; * MAKE-NODE and ELEMENT-MAKER appear to be the only procedures that
;;;   restrict the dimensions to which this program applies.  The former does
;;;   so by requiring X and Y as arguments, while the latter needs to handle
;;;   higher-dimensional (dim > 1) faces of simplices.

;;; * None of this code actually does anything, of course, since it's
;;;   completely abstract.  To generalize this to higher dimensions, the
;;;   appropriate integrators and constructors would still need to be supplied,
;;;   which can be a non-trivial task.  How much of what's in 2d-domains.scm
;;;   generalized to higher dimensions (assuming regions simple enough to be
;;;   triangulated efficiently)?

;;; Our basic reference is:
;;; Vichnevetsky, Robert. _Computer Methods for Partial Differential Equations,
;;; Volume 1: Elliptic Equatins and the Finite-Element Method_. Prentice-Hall:
;;; Englewood Cliffs, New Jersey, 1981.

;;; Note that elements are implicitly accessed through nodes, so that we should
;;; never have to directly refer to elements.  Also, the program only computes
;;; the system of linear algebraic equations; it does not attempt to solve
;;; them.

;;; A possible future direction is to extend this program to systems of
;;; equations.  How does one handle non-linear equations in general?

(declare (usual-integrations))


;;; Use FEM to produce a matrix:

(define (fem source nodes potential)
  (initialize-values nodes potential)
  (assemble-equations source nodes))


;;; Set boundary values:

(define (initialize-values nodes f)
  (let ((size (vector-length nodes)))
    (do ((i 0 (+ i 1)))
        ((>= i size))
      (let ((node (vector-ref nodes i)))
        (node:set-value! node (f node))))))


;;; The Rayleigh-Ritz method, as described in Vichnevetsky.  Actually, since
;;; MAKE-ELEMENT already incurs most of the cost of discretization up front, we
;;; only need to assemble the equations.

;;; NODE:ASSEMBLE returns a SPARSE-MATRIX data structure, as described in
;;; sparse.scm.  It can be used directly as input to SOR, or converted into a
;;; matrix and solved by LU decomposition.

(define (assemble-equations source nodes)

  ;; SOURCE is a function from R^2 to R, and NODES is expected to be a vector.

  (let* ((ncount (vector-length nodes))
         (bcount 0)
         (index-map (make-vector ncount)))

    ;; First, assign each node an index and count the number of boundary nodes.

    (do ((i 0 (+ i 1)))
        ((>= i ncount))
      (node:set-id! (vector-ref nodes i) i)
      (if (node:boundary? (vector-ref nodes i))
          (set! bcount (+ bcount 1))))

    ;; Next, create a mapping from node indices into matrix row number. (The
    ;; matrix has one row per interior node.)

    ;; This enforces the constraint that the equations satisfy the boundary
    ;; conditions.  Note that we can enforce the constraint before *or* after
    ;; minimizing the action.  If we do it before, everything is fine.  If we
    ;; do it in the other order, then we have to justify dropping the
    ;; constraint equations.  (Why does this sound vaguely familiar?  Does it
    ;; have anything at all to do with nonholonomic constraints?)

    ;; In any case, we can drop the equations associated with boundary nodes by
    ;; enforcing the constraint *before* we differentiate.  Perhaps this cannot
    ;; be done with the wave equation?  Would the Neumann condition along the
    ;; initial line look like some kind of nonholonomic constraint when we
    ;; perform the constrained minimization of the action?  But ODEs do not
    ;; have this problem...

    (let loop ((i 0) (row 0))
      (if (< i ncount)
	  (if (node:boundary? (vector-ref nodes i))
	      (begin
		(vector-set! index-map i #f)
		(loop (+ i 1) row))
	      (begin
		(vector-set! index-map i row)
		(loop (+ i 1) (+ row 1))))))

    ;; Loop over the nodes to create row entries:

    (let* ((icount (- ncount bcount))
           (big-matrix (make-sparse-matrix icount (1+ icount))))

      (do ((i 0 (+ i 1)))
          ((>= i ncount))

        (if (not (node:boundary? (vector-ref nodes i)))
            (let ((row (vector-ref index-map i)))

              ;; Compute the source term for this row:

	      (sparse-matrix-set! big-matrix row icount
				  (node:compute-source (vector-ref nodes i)
						       source))

              ;; Combine boundary values:

              (for-each
               (lambda (pair)
                 (let ((id (car pair))
                       (val (cadr pair)))
                   (if (node:boundary? (vector-ref nodes id))
                       (sparse-matrix-set!
			big-matrix row icount
			(- (sparse-matrix-ref big-matrix row icount)
			   (* val (node:get-value (vector-ref nodes id)))))
                       (sparse-matrix-set! big-matrix row
					   (vector-ref index-map id) val))))
               (node:assemble (vector-ref nodes i))))))

      big-matrix)))


;;; These procedures localize some of the assembly process in nodes; they are
;;; defined separately from nodes themselves to isolate the definitions
;;; specific to this particular application, thus increasing the generality of
;;; the definitions in this file.

(define (node:assemble node)
  (let ((l (append-map
	    (lambda (element index)
	      (element:compute-integrals element index))
	    (node:get-elements node)
	    (node:get-local-ids node))))
    (merge-terms l + (lambda (x y) (< (car x) (car y))))))

(define (node:compute-source node source)
  (apply + (map (lambda (element index)
		  (element:compute-source element source index))
		(node:get-elements node)
		(node:get-local-ids node))))


;;; This is a useful helper procedure: Given a list L of the form L = ((index1
;;; val1) (index2 val2) ...) and a procedure COMBINE, use COMBINE to
;;; concatenate the values of elements of L with the same index.

(define (merge-terms l combine <)

  ;; Sort first, then accumulate.  This is O(n log n), which is after than the
  ;; obvious O(n^2) algorithm.

  (if (null? l)
      '()
      (let* ((l (sort l <))
	     (indices (map car l))
	     (values (map cadr l)))
	(let loop ((indices (cdr indices))
		   (vals (cdr values))
		   (result '())
		   (id (car indices))
		   (accum (car values)))
	  (if (null? indices)
	      (cons (list id accum) result)
	      (if (eq? (car indices) id)
		  (loop (cdr indices) (cdr vals) result id
			(combine accum (car vals)))
		  (loop (cdr indices) (cdr vals)
			(cons (list id accum) result)
			(car indices) (car vals))))))))



;;; Let's now define elements.  To be completely general (in terms of
;;; dimensions of applicability), we should allow the construction of nodes on
;;; higher-dimensional (> 1) faces.

;;; Note that this implicitly assumes that elements are the convex hull of
;;; their vertices.

;;; The constructor for element-constructors:

(define (element-maker make-operator
		       make-integrator
		       make-basis-function)

  ;; MAKE-INTEGRATOR should take as argument a list of nodes, and returns a
  ;; procedure that takes a variable number of functions (at least 1) and
  ;; integrates their product over the domain specified implicitly as the
  ;; convex hull of the vertex nodes.

  ;; MAKE-BASIS-FUNCTION should take as argument a list of nodes and the index
  ;; of the node that is to be the center of the basis function, and return
  ;; some structure representing basis functions.

  ;; We place no restrictions on the representation of functions over elements,
  ;; so long as the particular instances of MAKE-BASIS-FUNCTION and
  ;; MAKE-INTEGRATOR agree a-priori on the representation.

  ;; MAKE-OPERATOR should take a list of nodes and return LEFT-OP, RIGHT-OP,
  ;; and COMBINE procedure, satisfying (INTEGRATE (COMBINE (LEFT-OP F)
  ;; (RIGHT-OP G))) = (INTEGRATE F (OP G)), i.e. implement integration by parts
  ;; so that basis functions can be less smooth.

  ;; The list of nodes facilitates the interpolation of variable coefficients
  ;; in the operator.  This may not be a good interface, as it makes artificial
  ;; assumptions on the contract between basis functions and operators (as is
  ;; the explicit use of LEFT-OP and RIGHT-OP).

  (define (make-element vertex-nodes other-nodes)

    ;; The first part stores the coefficients, the second part the source
    ;; terms.  What about coefficients?  Maybe we should incorporate the
    ;; source term into the differential operator.

    (let* ((nodes (append vertex-nodes other-nodes))
	   (number-of-nodes (length nodes))
	   (n-choose-2 (choose (+ number-of-nodes 2) 2))
	   (element
	    (vector (make-vector n-choose-2 0)
		    (make-vector n-choose-2 0)
		    vertex-nodes
		    other-nodes
		    (make-vector number-of-nodes #f)))
	   (op (make-operator nodes)))

      ;; Add the element to the nodes:

      (let loop ((nodes nodes) (i 0))
	(if (not (null? nodes))
	    (begin
	      (node:add-element (car nodes) element i)
	      (loop (cdr nodes) (+ i 1)))))

      ;; Initiailize elements (and hiding the hair)...

      (let ((integrate (make-integrator vertex-nodes))
	    (local-form (operator:get-local-form op)))

	(do ((i 0 (+ i 1)))
	    ((>= i number-of-nodes))
	  (element:set-basis-function!
	   element i (make-basis-function nodes i)))

	(do ((i 0 (+ i 1)))
	    ((>= i number-of-nodes))

	  (let ((f (element:get-basis-function element i)))

	    (do ((j i (+ j 1)))
		((>= j number-of-nodes))
	      (let ((g (element:get-basis-function element j)))
		(element:set-coeff! element i j
				    (integrate (local-form f g)))
		(element:set-source! element i j (integrate f g)))))))

      element))
  make-element)


;;; Methods for accessing the data structure:

(define (element:get-coeff element i j)
  (vector-ref (vector-ref element 0) (symmetric->vector-index i j)))

(define (element:set-coeff! element i j val)
  (vector-set! (vector-ref element 0) (symmetric->vector-index i j) val))

(define (element:get-source element i j)
  (vector-ref (vector-ref element 1) (symmetric->vector-index i j)))

(define (element:set-source! element i j val)
  (vector-set! (vector-ref element 1) (symmetric->vector-index i j) val))

(define (element:get-vertex-nodes element)
  (vector-ref element 2))

(define (element:get-non-vertex-nodes element)
  (vector-ref element 3))

(define (element:get-nodes element)
  (append (element:get-vertex-nodes element)
	  (element:get-non-vertex-nodes element)))

(define (element:node-count element)
  (length (element:get-nodes element)))

(define (element:set-basis-function! element i basis-function)
  (vector-set! (vector-ref element 4) i basis-function))

(define (element:get-basis-function element i)
  (vector-ref (vector-ref element 4) i))


;;; Computing the data needed in assembly:

(define (element:compute-source element source i)
  (let loop ((nodes (element:get-nodes element)) (sum 0.) (j 0))
    (if (null? nodes)
	sum
	(loop (cdr nodes)
	      (* (source (car nodes)) (element:get-source element i j))
	      (+ j 1)))))

(define (element:compute-integrals element i)
  (let loop ((nodes (element:get-nodes element)) (l '()) (j 0))
    (if (null? nodes)
	l
	(loop (cdr nodes)
	      (cons (list (node:get-id (car nodes))
			  (element:get-coeff element i j))
		    l)
	      (+ j 1)))))


;;; It's useful to have constant source functions:

(define (make-constant-function const)
  (lambda (node) const))

(define 0-function (make-constant-function 0))
