;;; This file defines some examples of things we can do on manifolds.  First,
;;; load some files that need to be compiled:


;;; Let's make a torus!

(define circle (make-sphere 1))
(define torus (product-manifold circle circle))


;;; Now we need the tangent bundle of the circle, or this wouldn't make very
;;; much sense...

(define TS^1 (make-tangent-bundle circle))

;;; The tangent bundle of the circle is trivial, because the circle is a Lie
;;; group!

(define cylinder TS^1)

;;; Here's something to integrate:

(define circle-path
  (v.field->flow circle
		 (v.field->local-field-maker circle-field)
		 (make-rk4-integrator (* 2 pi 1e-3))
		 local-distortion))

(define (real-circ t)
  (vector (cos t) (sin t)))

;;; Try this:
;(vector:distance (cadar (circle-path (vector 1 0) (* 2 pi))) (vector 1 0))


;;; Make a nonlinear pendulum:

(define pend-field (make-pendulum 1 1 1))

(define pend-path
  (v.field->flow cylinder
		 (v.field->local-field-maker pend-field)
		 (make-rk4-integrator 1e-3)
		 local-distortion))

(define pend-energy (make-pendulum-energy-function 1 1 1))

(define pend-init (imbedding->tangent circle (vector 0 1) (vector 0 0)))

;;; Try this:
;(define pend-results (pend-path pend-init 1))

;;; We should check things like the conservation of energy.  Also, do it with
;;; other methods (such as the imbedding) so that we have something to compare.

;;; The spherical pendulum isn't much harder to make:

(define S^2 (make-sphere 2))
(define TS^2 (make-tangent-bundle S^2))

(define make-sphere-tangent
  (let* ((chart (make-tangent-chart (make-spherical-chart 2 '(2 0 1) 0)))
	 (coords->point (chart:get-inverse-map chart)))

    (lambda (latitude longitude dlat dlong)
      (let ((p (coords->point (vector latitude longitude dlat dlong))))
	;(write-line `(initial position: ,(tangent:get-anchor p)))
	;(write-line `(initial velocity: ,(tangent->imbedded-velocity p)))
	p))))

(define spherical-init (make-sphere-tangent (* 3/4 pi) 0 0 1))

(define spherical-field (make-spherical-pendulum 1 1 1))

(define spherical-path
  (v.field->flow TS^2
		 (v.field->local-field-maker spherical-field)
		 (make-rk4-integrator 1e-3)
		 local-distortion))


;;; An example of something defined using Lagrangian methods:

(define free-L-field (lagrangian->v.field (make-free-particle-lagrangian 1)))


;;; An example of something defined using Hamiltonian methods:

(define free-H-field (hamiltonian->v.field (make-free-particle-hamiltonian 1)))


;;; Let's try the spherical pendulum again:

(define spherical-inclusion
  (smooth-map:diff (make-simple-map S^2 R^3 identity)))

(define spherical-lagrangian
  (smooth-map:compose falling-lagrangian spherical-inclusion))

(define spherical-L-field

  ;; Note that this works only in numerical mode -- Due to some structural
  ;; problems, ScmUtils won't do the right thing on smooth functions on tangent
  ;; bundles that are *not* derived from maps on the base space.

  (lagrangian->v.field spherical-lagrangian))

(define spherical-L-init
  (imbedding->tangent S^2 (vector 1 0 0) (vector 0 1 0)))

(define spherical-L-path
  (v.field->flow TS^2
		 (v.field->local-field-maker spherical-L-field)
		 (make-rk4-integrator 1e-3)
		 local-distortion))


;;; And again...

(define T*S^2 (make-cotangent-bundle S^2))

(define spherical-inclusion*
  (let* ((chart (car (manifold:get-finite-atlas R^3)))
	 (f (lambda (v)
	      (apply make-cotangent
		     (cons chart (cotangent->imbedding v))))))
    (make-simple-map T*S^2 T*R^3 f)))

(define spherical-hamiltonian
  (smooth-map:compose falling-hamiltonian spherical-inclusion*))

(define spherical-H-field
  (if #t
      (make-spherical-H-pendulum 1 1 1)
      (hamiltonian->v.field spherical-hamiltonian)))

(define spherical-H-init
  (imbedding->cotangent S^2 (vector 1 0 0) (vector 0 1 .5)))

(define spherical-H-path
  (v.field->flow T*S^2
		 (v.field->local-field-maker spherical-H-field)
		 (make-rk4-integrator 1e-3)
		 (check-vector-conservation-law
		  (smooth-map:get-point-function spherical-hamiltonian)
		  spherical-H-init)))

(define (spherical-H-angular-momentum cv)
  (let ((p (cadr (cotangent->imbedding cv)))
	(q (cotangent:get-anchor cv)))
    (let ((px (vector-ref p 0))
	  (py (vector-ref p 1))
	  (x (vector-ref q 0))
	  (y (vector-ref q 1)))
      (- (* x py) (* y px)))))

(define (print-spherical-H-state pair port)
  (let ((p (cadr (cotangent->imbedding (cadr pair))))
	(q (cotangent:get-anchor (cadr pair))))
    (let ((t (car pair))
	  (px (vector-ref p 0))
	  (py (vector-ref p 1))
	  (pz (vector-ref p 2))
	  (x (vector-ref q 0))
	  (y (vector-ref q 1))
	  (z (vector-ref q 2)))
      (display t port)
      (for-each (lambda (val)
		  (display " " port)
		  (display val port))
		(list x y z px py pz))
      (newline port))))


;;; Here's a test of rigid bodies and Euler angles.

(define euler-angles (make-euler-angles 0 1 (/ -pi 2) (/ pi 2)))

(define rigid-body-energy
  (make-free-rigid-body-lagrangian 1. (sqrt 2) 2.))

(define rigid-body-momentum
  (make-free-rigid-body-angular-momentum 1. (sqrt 2) 2.))

(define rigid-body-init
  (make-tangent euler-angles
		(chart:coords->point (vector 0 1 0) euler-angles)
		(vector .1 .1 .1)))

(define singular-init
  (chart:coords->point (vector 0. 1. 0. -.01 -.1 -.01)
		       (make-tangent-chart euler-angles)))

(define bad-init
  (chart:coords->point
   #(-1.309711193385365 .12149475001297212 1.1518832285401293
    -.19763392062291368 .09649536172931211 .18861249508985967)
   (list-ref (manifold:get-finite-atlas TSO3) 0)))

(define make-rigid-body-field
  (free-rigid-body-field-maker 1 (sqrt 2) 2))

(define energy+momentum
  (let ((E (smooth-map:get-point-function rigid-body-energy))
	(L (smooth-map:get-point-function rigid-body-momentum)))
    (lambda (p)
      (vector-append (E p) (L p)))))

(define (correct->traditional-order v)
  (let ((psi (vector-ref v 0))
	(theta (vector-ref v 1))
	(phi (vector-ref v 2))
	(psidot (vector-ref v 3))
	(thetadot (vector-ref v 4))
	(phidot (vector-ref v 5)))
    (vector 0. theta phi psi thetadot phidot psidot)))

(define scmutils-energy+momentum
  (let ((energy (make-rigid-body-energy 1. (sqrt 2) 2.))
	(momentum (make-rigid-body-momentum 1. (sqrt 2) 2.)))
    (lambda (state)
      (let ((L (momentum state))
	    (E (energy state)))
	(vector E (vector-first L) (vector-second L) (vector-third L))))))

(define rigid-body-path
  (v.field->flow
   TSO3
   make-rigid-body-field
   (if *using-scmutils?*
       (begin
	 (set! *ode-integration-method* 'bulirsch-stoer)
	 (make-scmutils-integrator .01 1e-12))
       (make-rk4-integrator .01))
   (check-vector-conservation-law
    (if *using-scmutils?*
	scmutils-energy+momentum
	energy+momentum)
    singular-init)))
