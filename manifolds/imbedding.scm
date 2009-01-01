;;; This file defines some procedures that are useful for working with
;;; imbedding representations of manifolds:

(declare (usual-integrations))


;;; For tangent vectors:

(define (imbedding->tangent M p v)
  (let ((U (manifold:find-best-chart M p)))
    (make-tangent U p (push-forward-in-coords (chart:get-coord-map U) p v))))

(define (tangent->imbedding v)
  (let* ((U (tangent:get-chart v))
	 (p (tangent:get-anchor v)))
    (list p
	  (push-forward-in-coords (chart:get-inverse-map U)
				  (chart:point->coords p U)
				  (tangent:get-coords v)))))

(define tangent->imbedded-velocity
  (compose cadr tangent->imbedding))


;;; For cotangent vectors:

(define (imbedding->cotangent M p v)
  (let ((U (manifold:find-best-chart M p)))
    (make-cotangent U p (pull-back-in-coords (chart:get-inverse-map U)
					     (chart:point->coords p U)
					     v))))

(define (cotangent->imbedding v)
  (let* ((U (cotangent:get-chart v))
	 (p (cotangent:get-anchor v)))
    (list p

	  ;; Need to project the pulled-back functional onto the imbedded
	  ;; surface because it's represented in the standard basis of the
	  ;; ambient space, not the basis of the tangent space to M imbedded
	  ;; inside R^n.

	  (project-onto-basis
	   (make-imbedded-basis U p)
	   (pull-back-in-coords (chart:get-coord-map U)
				p
				(cotangent:get-coords v))))))

(define (make-imbedded-basis chart x)
  (let ((dim (chart:dimension chart)))
    (let loop ((i 0) (vlist '()))
      (if (< i dim)
	  (loop (+ i 1) (cons (tangent->imbedded-velocity
			       (make-tangent chart x (vector:basis dim i 1)))
			      vlist))
	  vlist))))
