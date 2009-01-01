;;; This file defines procedures that need to be easily modifiable:


;;; Methods for generating equations from triangulated domains:

(define combine-equations-without-overlap
  (pde:equation-maker merge-equations))

(define combine-equations-with-overlap
  (pde:equation-maker
   (append-constraint-equations make-ordered-boundary-constraints)))

(define combine-equations-with-overlap1
  (pde:equation-maker
   (append-constraint-equations make-all-constraints)))

(define combine-equations-with-overlap2
  (pde:equation-maker
   (append-constraint-equations make-all-ordered-constraints)))

(define combine-equations-using-CMPGRD
  (pde:equation-maker cmpgrd:combine-equations))


;;; Methods for triangulating domains:

(define pde:make-domain-without-overlaps
  (pde:domain-maker generate-node-lists exact-overlap))

(define pde:make-domain-with-small-overlaps
  (pde:domain-maker copy-between-node-lists exact-overlap))

(define pde:make-domain-with-overlaps
  (pde:domain-maker make-nodes-for-each-chart reduce-overlap))

(define pde:make-domain-with-larger-overlaps
  (pde:domain-maker make-nodes-for-each-chart extended-overlap))

(define pde:make-simple-domain
  (pde:domain-maker make-nodes-for-each-chart do-nothing-to-complex))
