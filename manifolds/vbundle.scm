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

;;; Definitely need cotangent bundles!

;;; OLD COMMENTS THAT ARE STILL RELEVANT:

;;; We Shouldn't be too particular about vector representations, and should
;;; leave several different options.
;;;
;;; The differential operator representation is a good one, though, since via
;;; compositions we can generate most other differential operators of interest.


;;; OLD COMMENTS THAT APPEAR TO BE IRRELEVANT:

;;; Some of the issues involved in making tangent bundles, etc.:
;;;
;;; Something that looks like vector bundles will allow us to implement vector
;;; fields and differential forms.  Maybe differential operators, too.
;;;
;;; Should we make this more like a manifold?
;;;
;;; Maybe the right thing is to make local trivializations a special kind of
;;; chart (just as product charts are).  We need standard ways to attach and
;;; access optional structures on charts (and manifolds), such as product
;;; structures, metric structures, symplectic structures, etc.
;;;
;;; Is the vector space structure necessary on the coordinate-free level?
;;;
;;; And when do we ever need to treat the bundle as a manifold?
;;;
;;; Make tangent and cotangent bundles?  Metric tensors, symplectic forms,
;;; differential operators are all related to these bundles...
;;;
;;; GJS used nscmutils to directly differentiate the functions.  Maybe we
;;; should consider using that at some point, too.  Tangent vectors appear to
;;; be represented by first-order differential operators on functions.
;;;
;;; BUT, we don't know that charts always map to an imbedding of the manifold
;;; in a euclidean space!  What if it uses some other representation on the
;;; other end?  We can numerically differentiate transition maps, but in
;;; general not charts.  So use some other representation of tangent vectors
;;; and vector fields?
;;;
;;; Chart and vector may be the best represenatation; that appears to be what
;;; GJS uses, too.

(declare (usual-integrations))
(load "tangent")
(load "imbedding")
(load "cotangent")

;;; Abstract vector bundles:

(define (make-vector-bundle M E proj fiber)

  ;; PROJ takes a pair (x,v) and returns x, FIBER yields the operations for the
  ;; vector space structure on the fiber above x.

  (let ((dim (- (manifold:dimension E) (manifold:dimension M))))
    (manifold:install-extra E 'vector-bundle (vector (delay M)
						     proj
						     fiber
						     dim)))
  E)

(define (vbundle:get-manifold E)
  (let ((structs (manifold:get-extra E 'vector-bundle)))
    (if structs
	(force (vector-ref structs 0))
	#f)))

(define (vbundle:get-fiber-map E)
  (let ((structs (manifold:get-extra E 'vector-bundle)))
    (if structs
	(vector-ref structs 2)
	#f)))

(define (vbundle:get-projection E)
  (let ((structs (manifold:get-extra E 'vector-bundle)))
    (if structs
	(vector-ref structs 1)
	#f)))

(define (vbundle:dimension E)
  (let ((structs (manifold:get-extra E 'vector-bundle)))
    (if structs
	(vector-ref structs 3)
	#f)))


;;; Abstract vector spaces (fibers of the bundle) above each point of the
;;; manifold:

(define (make-fiber + - * member?)
  (vector + - * member?))

(define (fiber:get+ fiber)
  (vector-ref fiber 0))

(define (fiber:get- fiber)
  (vector-ref fiber 1))

(define (fiber:get* fiber)
  (vector-ref fiber 2))

(define (fiber:get-membership-test fiber)
  (vector-ref fiber 3))

(define (fiber:+ fiber v w)
  ((fiber:get+ fiber) v w))

(define (fiber:- fiber v w)
  ((fiber:get- fiber) v w))

(define (fiber:* fiber a v)
  ((fiber:get* fiber) a v))

(define (fiber:member? v fiber)
  ((fiber:get-membership-test fiber) v))
