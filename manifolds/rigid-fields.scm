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

;;; These vector fields are machine-generated:

(define (rigid-field-0 a b c x v)
  (let ((psi (vector-ref x 0))
	(theta (vector-ref x 1))
	(phi (vector-ref x 2))
	(psidot (vector-ref v 0))
	(thetadot (vector-ref v 1))
	(phidot (vector-ref v 2)))
    (vector psidot thetadot phidot
	    (/ (+ (* (cos theta) (+ (* (cos theta) (+ (* (cos psi) (+ (* (cos
	     psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* a b c phidot thetadot)
	     (* (expt b 2) c phidot thetadot))) (* (sin theta) (sin psi) (+ (*
	     -1 (expt a 2) c (expt phidot 2)) (* (expt b 2) c (expt phidot
	     2)))))) (* (expt (sin psi) 2) (+ (* a (+ (* a c phidot thetadot)
	     (* 2 b c phidot thetadot))) (* (expt b 2) c phidot thetadot))) (*
	     -1 b (expt c 2) phidot thetadot))) (* (sin theta) (sin psi) (+ (*
	     (expt (sin psi) 2) (+ (* -1 (expt a 2) c (expt phidot 2)) (* (expt
	     b 2) c (expt phidot 2)))) (* a (expt c 2) (expt phidot 2)) (* -1 b
	     (expt c 2) (expt phidot 2)))))) (* (expt (sin psi) 2) (+ (* (expt
	     (sin psi) 2) a (+ (* a c phidot thetadot) (* b c phidot
	     thetadot))) (* -1 a (expt c 2) phidot thetadot))))) (* (cos psi)
	     (+ (* (cos psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* -1 a b c
	     psidot thetadot) (* (expt b 2) c psidot thetadot))) (* (sin theta)
	     (sin psi) (+ (* -1 (expt a 2) c phidot psidot) (* (expt b 2) c
	     phidot psidot))))) (* (expt (sin psi) 2) (+ (* a (+ (* a c psidot
	     thetadot) (* -2 b c psidot thetadot))) (* (expt b 2) c psidot
	     thetadot))) (* -1 b (expt c 2) psidot thetadot))) (* (sin theta)
	     (sin psi) (+ (* (expt (sin psi) 2) (+ (* -1 (expt a 2) c phidot
	     psidot) (* (expt b 2) c phidot psidot))) (* a (expt c 2) phidot
	     psidot) (* -1 b (expt c 2) phidot psidot))))) (* (expt (sin psi)
	     2) (+ (* (expt (sin psi) 2) a (+ (* a c psidot thetadot) (* -1 b c
	     psidot thetadot))) (* -1 a (expt c 2) psidot thetadot))))) (* (cos
	     psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* (cos
	     psi) (+ (* (cos psi) (expt (sin theta) 2) a (+ (* -1 a b phidot
	     thetadot) (* (expt b 2) phidot thetadot))) (* (sin theta) (+ (*
	     (expt (sin theta) 2) (sin psi) a (+ (* -1 a b (expt phidot 2)) (*
	     (expt b 2) (expt phidot 2)))) (* (sin psi) a (+ (* a b (expt
	     thetadot 2)) (* -1 (expt b 2) (expt thetadot 2)))))))) (* (expt
	     (sin theta) 2) (+ (* (expt (sin psi) 2) a (+ (* -1 a b phidot
	     thetadot) (* (expt b 2) phidot thetadot))) (* a b c phidot
	     thetadot))))) (* (sin theta) (+ (* (expt (sin theta) 2) (expt (sin
	     psi) 3) a (+ (* -2 a b (expt phidot 2)) (* 2 (expt b 2) (expt
	     phidot 2)))) (* (expt (sin psi) 3) a (+ (* 2 a b (expt thetadot
	     2)) (* -2 (expt b 2) (expt thetadot 2)))))))) (* (expt (sin theta)
	     2) (expt (sin psi) 2) (+ (* (expt (sin psi) 2) a (+ (* a b phidot
	     thetadot) (* -1 (expt b 2) phidot thetadot))) (* 2 a b c phidot
	     thetadot))))) (* (sin theta) (+ (* (expt (sin theta) 2) (expt (sin
	     psi) 5) a (+ (* -1 a b (expt phidot 2)) (* (expt b 2) (expt phidot
	     2)))) (* (expt (sin psi) 5) a (+ (* a b (expt thetadot 2)) (* -1
	     (expt b 2) (expt thetadot 2)))))))) (* (expt (sin theta) 2) (expt
	     (sin psi) 4) (+ (* (expt (sin psi) 2) a (+ (* a b phidot thetadot)
	     (* -1 (expt b 2) phidot thetadot))) (* a b c phidot thetadot))))
	     (+ (* (expt (cos psi) 2) (+ (* (expt (cos psi) 2) (sin theta) a b
	     c) (* 2 (sin theta) (expt (sin psi) 2) a b c))) (* (sin theta)
	     (expt (sin psi) 4) a b c))) (/ (+ (* (cos theta) (+ (* (sin theta)
	     (+ (* (expt (cos psi) 2) a (+ (* a (expt phidot 2)) (* -1 c (expt
	     phidot 2)))) (* (expt (sin psi) 2) b (+ (* b (expt phidot 2)) (*
	     -1 c (expt phidot 2)))))) (* (cos psi) (sin psi) (+ (* a (+ (* -1
	     a phidot thetadot) (* c phidot thetadot))) (* b (+ (* b phidot
	     thetadot) (* -1 c phidot thetadot))))))) (* (sin theta) (+ (*
	     (expt (cos psi) 2) (+ (* -1 (expt (cos psi) 2) a b phidot psidot)
	     (* -2 (expt (sin psi) 2) a b phidot psidot) (* a (+ (* a phidot
	     psidot) (* -1 c phidot psidot))))) (* (expt (sin psi) 2) (+ (* -1
	     (expt (sin psi) 2) a b phidot psidot) (* b (+ (* b phidot psidot)
	     (* -1 c phidot psidot))))))) (* (cos psi) (sin psi) (+ (* a (+ (*
	     -1 a psidot thetadot) (* c psidot thetadot))) (* b (+ (* b psidot
	     thetadot) (* -1 c psidot thetadot)))))) (+ (* (expt (cos psi) 2)
	     (+ (* (expt (cos psi) 2) a b) (* 2 (expt (sin psi) 2) a b))) (*
	     (expt (sin psi) 4) a b))) (/ (+ (* (cos theta) (+ (* (cos psi) (+
	     (* (cos psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* -1 a b phidot
	     thetadot) (* -1 (expt b 2) phidot thetadot))) (* (sin theta) (sin
	     psi) (+ (* (expt a 2) (expt phidot 2)) (* -1 (expt b 2) (expt
	     phidot 2)))))) (* (expt (sin psi) 2) (+ (* a (+ (* -1 a phidot
	     thetadot) (* -2 b phidot thetadot))) (* -1 (expt b 2) phidot
	     thetadot))) (* b c phidot thetadot))) (* (sin theta) (sin psi) (+
	     (* (expt (sin psi) 2) (+ (* (expt a 2) (expt phidot 2)) (* -1
	     (expt b 2) (expt phidot 2)))) (* -1 a c (expt phidot 2)) (* b c
	     (expt phidot 2)))))) (* (expt (sin psi) 2) (+ (* (expt (sin psi)
	     2) a (+ (* -1 a phidot thetadot) (* -1 b phidot thetadot))) (* a c
	     phidot thetadot))))) (* (cos psi) (+ (* (cos psi) (+ (* (cos psi)
	     (+ (* (cos psi) (+ (* a b psidot thetadot) (* -1 (expt b 2) psidot
	     thetadot))) (* (sin theta) (sin psi) (+ (* (expt a 2) phidot
	     psidot) (* -1 (expt b 2) phidot psidot))))) (* (expt (sin psi) 2)
	     (+ (* a (+ (* -1 a psidot thetadot) (* 2 b psidot thetadot))) (*
	     -1 (expt b 2) psidot thetadot))) (* b c psidot thetadot))) (* (sin
	     theta) (sin psi) (+ (* (expt (sin psi) 2) (+ (* (expt a 2) phidot
	     psidot) (* -1 (expt b 2) phidot psidot))) (* -1 a c phidot psidot)
	     (* b c phidot psidot))))) (* (expt (sin psi) 2) (+ (* (expt (sin
	     psi) 2) a (+ (* -1 a psidot thetadot) (* b psidot thetadot))) (* a
	     c psidot thetadot)))) (+ (* (expt (cos psi) 2) (+ (* (expt (cos
	     psi) 2) (sin theta) a b) (* 2 (sin theta) (expt (sin psi) 2) a
	     b))) (* (sin theta) (expt (sin psi) 4) a b))))))

(define (rigid-field-1 a b c x v)
  (let ((psi (vector-ref x 0))
	(theta (vector-ref x 1))
	(phi (vector-ref x 2))
	(psidot (vector-ref v 0))
	(thetadot (vector-ref v 1))
	(phidot (vector-ref v 2)))
    (vector psidot thetadot phidot
	    (/ (+ (* (cos theta) (+ (* (cos theta) (+ (* (cos psi) (+ (* (cos
	     psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* a b c phidot thetadot)
	     (* (expt b 2) c phidot thetadot))) (* (sin theta) (sin psi) (+ (*
	     -1 (expt a 2) c (expt phidot 2)) (* (expt b 2) c (expt phidot
	     2)))))) (* (expt (sin psi) 2) (+ (* a (+ (* a c phidot thetadot)
	     (* 2 b c phidot thetadot))) (* (expt b 2) c phidot thetadot))) (*
	     -1 b (expt c 2) phidot thetadot))) (* (sin theta) (sin psi) (+ (*
	     (expt (sin psi) 2) (+ (* -1 (expt a 2) c (expt phidot 2)) (* (expt
	     b 2) c (expt phidot 2)))) (* a (expt c 2) (expt phidot 2)) (* -1 b
	     (expt c 2) (expt phidot 2)))))) (* (expt (sin psi) 2) (+ (* (expt
	     (sin psi) 2) a (+ (* a c phidot thetadot) (* b c phidot
	     thetadot))) (* -1 a (expt c 2) phidot thetadot))))) (* (cos psi)
	     (+ (* (cos psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* -1 a b c
	     psidot thetadot) (* (expt b 2) c psidot thetadot))) (* (sin theta)
	     (sin psi) (+ (* -1 (expt a 2) c phidot psidot) (* (expt b 2) c
	     phidot psidot))))) (* (expt (sin psi) 2) (+ (* a (+ (* a c psidot
	     thetadot) (* -2 b c psidot thetadot))) (* (expt b 2) c psidot
	     thetadot))) (* -1 b (expt c 2) psidot thetadot))) (* (sin theta)
	     (sin psi) (+ (* (expt (sin psi) 2) (+ (* -1 (expt a 2) c phidot
	     psidot) (* (expt b 2) c phidot psidot))) (* a (expt c 2) phidot
	     psidot) (* -1 b (expt c 2) phidot psidot))))) (* (expt (sin psi)
	     2) (+ (* (expt (sin psi) 2) a (+ (* a c psidot thetadot) (* -1 b c
	     psidot thetadot))) (* -1 a (expt c 2) psidot thetadot))))) (* (cos
	     psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* (cos
	     psi) (+ (* (cos psi) (expt (sin theta) 2) a (+ (* -1 a b phidot
	     thetadot) (* (expt b 2) phidot thetadot))) (* (sin theta) (+ (*
	     (expt (sin theta) 2) (sin psi) a (+ (* -1 a b (expt phidot 2)) (*
	     (expt b 2) (expt phidot 2)))) (* (sin psi) a (+ (* a b (expt
	     thetadot 2)) (* -1 (expt b 2) (expt thetadot 2)))))))) (* (expt
	     (sin theta) 2) (+ (* (expt (sin psi) 2) a (+ (* -1 a b phidot
	     thetadot) (* (expt b 2) phidot thetadot))) (* a b c phidot
	     thetadot))))) (* (sin theta) (+ (* (expt (sin theta) 2) (expt (sin
	     psi) 3) a (+ (* -2 a b (expt phidot 2)) (* 2 (expt b 2) (expt
	     phidot 2)))) (* (expt (sin psi) 3) a (+ (* 2 a b (expt thetadot
	     2)) (* -2 (expt b 2) (expt thetadot 2)))))))) (* (expt (sin theta)
	     2) (expt (sin psi) 2) (+ (* (expt (sin psi) 2) a (+ (* a b phidot
	     thetadot) (* -1 (expt b 2) phidot thetadot))) (* 2 a b c phidot
	     thetadot))))) (* (sin theta) (+ (* (expt (sin theta) 2) (expt (sin
	     psi) 5) a (+ (* -1 a b (expt phidot 2)) (* (expt b 2) (expt phidot
	     2)))) (* (expt (sin psi) 5) a (+ (* a b (expt thetadot 2)) (* -1
	     (expt b 2) (expt thetadot 2)))))))) (* (expt (sin theta) 2) (expt
	     (sin psi) 4) (+ (* (expt (sin psi) 2) a (+ (* a b phidot thetadot)
	     (* -1 (expt b 2) phidot thetadot))) (* a b c phidot thetadot))))
	     (+ (* (expt (cos psi) 2) (+ (* (expt (cos psi) 2) (sin theta) a b
	     c) (* 2 (sin theta) (expt (sin psi) 2) a b c))) (* (sin theta)
	     (expt (sin psi) 4) a b c))) (/ (+ (* (cos theta) (+ (* (sin theta)
	     (+ (* (expt (cos psi) 2) a (+ (* a (expt phidot 2)) (* -1 c (expt
	     phidot 2)))) (* (expt (sin psi) 2) b (+ (* b (expt phidot 2)) (*
	     -1 c (expt phidot 2)))))) (* (cos psi) (sin psi) (+ (* a (+ (* -1
	     a phidot thetadot) (* c phidot thetadot))) (* b (+ (* b phidot
	     thetadot) (* -1 c phidot thetadot))))))) (* (sin theta) (+ (*
	     (expt (cos psi) 2) (+ (* -1 (expt (cos psi) 2) a b phidot psidot)
	     (* -2 (expt (sin psi) 2) a b phidot psidot) (* a (+ (* a phidot
	     psidot) (* -1 c phidot psidot))))) (* (expt (sin psi) 2) (+ (* -1
	     (expt (sin psi) 2) a b phidot psidot) (* b (+ (* b phidot psidot)
	     (* -1 c phidot psidot))))))) (* (cos psi) (sin psi) (+ (* a (+ (*
	     -1 a psidot thetadot) (* c psidot thetadot))) (* b (+ (* b psidot
	     thetadot) (* -1 c psidot thetadot)))))) (+ (* (expt (cos psi) 2)
	     (+ (* (expt (cos psi) 2) a b) (* 2 (expt (sin psi) 2) a b))) (*
	     (expt (sin psi) 4) a b))) (/ (+ (* (cos theta) (+ (* (cos psi) (+
	     (* (cos psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* -1 a b phidot
	     thetadot) (* -1 (expt b 2) phidot thetadot))) (* (sin theta) (sin
	     psi) (+ (* (expt a 2) (expt phidot 2)) (* -1 (expt b 2) (expt
	     phidot 2)))))) (* (expt (sin psi) 2) (+ (* a (+ (* -1 a phidot
	     thetadot) (* -2 b phidot thetadot))) (* -1 (expt b 2) phidot
	     thetadot))) (* b c phidot thetadot))) (* (sin theta) (sin psi) (+
	     (* (expt (sin psi) 2) (+ (* (expt a 2) (expt phidot 2)) (* -1
	     (expt b 2) (expt phidot 2)))) (* -1 a c (expt phidot 2)) (* b c
	     (expt phidot 2)))))) (* (expt (sin psi) 2) (+ (* (expt (sin psi)
	     2) a (+ (* -1 a phidot thetadot) (* -1 b phidot thetadot))) (* a c
	     phidot thetadot))))) (* (cos psi) (+ (* (cos psi) (+ (* (cos psi)
	     (+ (* (cos psi) (+ (* a b psidot thetadot) (* -1 (expt b 2) psidot
	     thetadot))) (* (sin theta) (sin psi) (+ (* (expt a 2) phidot
	     psidot) (* -1 (expt b 2) phidot psidot))))) (* (expt (sin psi) 2)
	     (+ (* a (+ (* -1 a psidot thetadot) (* 2 b psidot thetadot))) (*
	     -1 (expt b 2) psidot thetadot))) (* b c psidot thetadot))) (* (sin
	     theta) (sin psi) (+ (* (expt (sin psi) 2) (+ (* (expt a 2) phidot
	     psidot) (* -1 (expt b 2) phidot psidot))) (* -1 a c phidot psidot)
	     (* b c phidot psidot))))) (* (expt (sin psi) 2) (+ (* (expt (sin
	     psi) 2) a (+ (* -1 a psidot thetadot) (* b psidot thetadot))) (* a
	     c psidot thetadot)))) (+ (* (expt (cos psi) 2) (+ (* (expt (cos
	     psi) 2) (sin theta) a b) (* 2 (sin theta) (expt (sin psi) 2) a
	     b))) (* (sin theta) (expt (sin psi) 4) a b))))))

(define (rigid-field-2 a b c x v)
  (let ((psi (vector-ref x 0))
	(theta (vector-ref x 1))
	(phi (vector-ref x 2))
	(psidot (vector-ref v 0))
	(thetadot (vector-ref v 1))
	(phidot (vector-ref v 2)))

    (vector psidot thetadot phidot
	    (/ (+ (* (cos theta) (+ (* (cos theta) (+ (* (cos psi) (+ (* (cos
	     psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* -1 a b c phidot
	     thetadot) (* -1 (expt b 2) c phidot thetadot))) (* (sin theta)
	     (sin psi) (+ (* -1 (expt a 2) c (expt phidot 2)) (* (expt b 2) c
	     (expt phidot 2)))))) (* (expt (sin psi) 2) (+ (* a (+ (* -1 a c
	     phidot thetadot) (* -2 b c phidot thetadot))) (* -1 (expt b 2) c
	     phidot thetadot))) (* b (expt c 2) phidot thetadot))) (* (sin
	     theta) (sin psi) (+ (* (expt (sin psi) 2) (+ (* -1 (expt a 2) c
	     (expt phidot 2)) (* (expt b 2) c (expt phidot 2)))) (* a (expt c
	     2) (expt phidot 2)) (* -1 b (expt c 2) (expt phidot 2)))))) (*
	     (expt (sin psi) 2) (+ (* (expt (sin psi) 2) a (+ (* -1 a c phidot
	     thetadot) (* -1 b c phidot thetadot))) (* a (expt c 2) phidot
	     thetadot))))) (* (cos psi) (+ (* (cos psi) (+ (* (cos psi) (+ (*
	     (cos psi) (+ (* -1 a b c psidot thetadot) (* (expt b 2) c psidot
	     thetadot))) (* (sin theta) (sin psi) (+ (* (expt a 2) c phidot
	     psidot) (* -1 (expt b 2) c phidot psidot))))) (* (expt (sin psi)
	     2) (+ (* a (+ (* a c psidot thetadot) (* -2 b c psidot thetadot)))
	     (* (expt b 2) c psidot thetadot))) (* -1 b (expt c 2) psidot
	     thetadot))) (* (sin theta) (sin psi) (+ (* (expt (sin psi) 2) (+
	     (* (expt a 2) c phidot psidot) (* -1 (expt b 2) c phidot psidot)))
	     (* -1 a (expt c 2) phidot psidot) (* b (expt c 2) phidot
	     psidot))))) (* (expt (sin psi) 2) (+ (* (expt (sin psi) 2) a (+ (*
	     a c psidot thetadot) (* -1 b c psidot thetadot))) (* -1 a (expt c
	     2) psidot thetadot))))) (* (cos psi) (+ (* (cos psi) (+ (* (cos
	     psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* (cos psi) (expt (sin
	     theta) 2) a (+ (* a b phidot thetadot) (* -1 (expt b 2) phidot
	     thetadot))) (* (sin theta) (+ (* (expt (sin theta) 2) (sin psi) a
	     (+ (* -1 a b (expt phidot 2)) (* (expt b 2) (expt phidot 2)))) (*
	     (sin psi) a (+ (* a b (expt thetadot 2)) (* -1 (expt b 2) (expt
	     thetadot 2)))))))) (* (expt (sin theta) 2) (+ (* (expt (sin psi)
	     2) a (+ (* a b phidot thetadot) (* -1 (expt b 2) phidot
	     thetadot))) (* -1 a b c phidot thetadot))))) (* (sin theta) (+ (*
	     (expt (sin theta) 2) (expt (sin psi) 3) a (+ (* -2 a b (expt
	     phidot 2)) (* 2 (expt b 2) (expt phidot 2)))) (* (expt (sin psi)
	     3) a (+ (* 2 a b (expt thetadot 2)) (* -2 (expt b 2) (expt
	     thetadot 2)))))))) (* (expt (sin theta) 2) (expt (sin psi) 2) (+
	     (* (expt (sin psi) 2) a (+ (* -1 a b phidot thetadot) (* (expt b
	     2) phidot thetadot))) (* -2 a b c phidot thetadot))))) (* (sin
	     theta) (+ (* (expt (sin theta) 2) (expt (sin psi) 5) a (+ (* -1 a
	     b (expt phidot 2)) (* (expt b 2) (expt phidot 2)))) (* (expt (sin
	     psi) 5) a (+ (* a b (expt thetadot 2)) (* -1 (expt b 2) (expt
	     thetadot 2)))))))) (* (expt (sin theta) 2) (expt (sin psi) 4) (+
	     (* (expt (sin psi) 2) a (+ (* -1 a b phidot thetadot) (* (expt b
	     2) phidot thetadot))) (* -1 a b c phidot thetadot)))) (+ (* (expt
	     (cos psi) 2) (+ (* (expt (cos psi) 2) (sin theta) a b c) (* 2 (sin
	     theta) (expt (sin psi) 2) a b c))) (* (sin theta) (expt (sin psi)
	     4) a b c))) (/ (+ (* (cos theta) (+ (* (sin theta) (+ (* (expt
	     (cos psi) 2) a (+ (* a (expt phidot 2)) (* -1 c (expt phidot 2))))
	     (* (expt (sin psi) 2) b (+ (* b (expt phidot 2)) (* -1 c (expt
	     phidot 2)))))) (* (cos psi) (sin psi) (+ (* a (+ (* a phidot
	     thetadot) (* -1 c phidot thetadot))) (* b (+ (* -1 b phidot
	     thetadot) (* c phidot thetadot))))))) (* (sin theta) (+ (* (expt
	     (cos psi) 2) (+ (* (expt (cos psi) 2) a b phidot psidot) (* 2
	     (expt (sin psi) 2) a b phidot psidot) (* a (+ (* -1 a phidot
	     psidot) (* c phidot psidot))))) (* (expt (sin psi) 2) (+ (* (expt
	     (sin psi) 2) a b phidot psidot) (* b (+ (* -1 b phidot psidot) (*
	     c phidot psidot))))))) (* (cos psi) (sin psi) (+ (* a (+ (* -1 a
	     psidot thetadot) (* c psidot thetadot))) (* b (+ (* b psidot
	     thetadot) (* -1 c psidot thetadot)))))) (+ (* (expt (cos psi) 2)
	     (+ (* (expt (cos psi) 2) a b) (* 2 (expt (sin psi) 2) a b))) (*
	     (expt (sin psi) 4) a b))) (/ (+ (* (cos theta) (+ (* (cos psi) (+
	     (* (cos psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* -1 a b phidot
	     thetadot) (* -1 (expt b 2) phidot thetadot))) (* (sin theta) (sin
	     psi) (+ (* -1 (expt a 2) (expt phidot 2)) (* (expt b 2) (expt
	     phidot 2)))))) (* (expt (sin psi) 2) (+ (* a (+ (* -1 a phidot
	     thetadot) (* -2 b phidot thetadot))) (* -1 (expt b 2) phidot
	     thetadot))) (* b c phidot thetadot))) (* (sin theta) (sin psi) (+
	     (* (expt (sin psi) 2) (+ (* -1 (expt a 2) (expt phidot 2)) (*
	     (expt b 2) (expt phidot 2)))) (* a c (expt phidot 2)) (* -1 b c
	     (expt phidot 2)))))) (* (expt (sin psi) 2) (+ (* (expt (sin psi)
	     2) a (+ (* -1 a phidot thetadot) (* -1 b phidot thetadot))) (* a c
	     phidot thetadot))))) (* (cos psi) (+ (* (cos psi) (+ (* (cos psi)
	     (+ (* (cos psi) (+ (* -1 a b psidot thetadot) (* (expt b 2) psidot
	     thetadot))) (* (sin theta) (sin psi) (+ (* (expt a 2) phidot
	     psidot) (* -1 (expt b 2) phidot psidot))))) (* (expt (sin psi) 2)
	     (+ (* a (+ (* a psidot thetadot) (* -2 b psidot thetadot))) (*
	     (expt b 2) psidot thetadot))) (* -1 b c psidot thetadot))) (* (sin
	     theta) (sin psi) (+ (* (expt (sin psi) 2) (+ (* (expt a 2) phidot
	     psidot) (* -1 (expt b 2) phidot psidot))) (* -1 a c phidot psidot)
	     (* b c phidot psidot))))) (* (expt (sin psi) 2) (+ (* (expt (sin
	     psi) 2) a (+ (* a psidot thetadot) (* -1 b psidot thetadot))) (*
	     -1 a c psidot thetadot)))) (+ (* (expt (cos psi) 2) (+ (* (expt
	     (cos psi) 2) (sin theta) a b) (* 2 (sin theta) (expt (sin psi) 2)
	     a b))) (* (sin theta) (expt (sin psi) 4) a b))))))

(define (rigid-field-3 a b c x v)
  (let ((psi (vector-ref x 0))
	(theta (vector-ref x 1))
	(phi (vector-ref x 2))
	(psidot (vector-ref v 0))
	(thetadot (vector-ref v 1))
	(phidot (vector-ref v 2)))
    (vector psidot thetadot phidot
	    (/ (+ (* (cos theta) (+ (* (cos theta) (+ (* (cos psi) (+ (* (cos
	     psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* -1 a b c phidot
	     thetadot) (* -1 (expt b 2) c phidot thetadot))) (* (sin theta)
	     (sin psi) (+ (* -1 (expt a 2) c (expt phidot 2)) (* (expt b 2) c
	     (expt phidot 2)))))) (* (expt (sin psi) 2) (+ (* a (+ (* -1 a c
	     phidot thetadot) (* -2 b c phidot thetadot))) (* -1 (expt b 2) c
	     phidot thetadot))) (* b (expt c 2) phidot thetadot))) (* (sin
	     theta) (sin psi) (+ (* (expt (sin psi) 2) (+ (* -1 (expt a 2) c
	     (expt phidot 2)) (* (expt b 2) c (expt phidot 2)))) (* a (expt c
	     2) (expt phidot 2)) (* -1 b (expt c 2) (expt phidot 2)))))) (*
	     (expt (sin psi) 2) (+ (* (expt (sin psi) 2) a (+ (* -1 a c phidot
	     thetadot) (* -1 b c phidot thetadot))) (* a (expt c 2) phidot
	     thetadot))))) (* (cos psi) (+ (* (cos psi) (+ (* (cos psi) (+ (*
	     (cos psi) (+ (* -1 a b c psidot thetadot) (* (expt b 2) c psidot
	     thetadot))) (* (sin theta) (sin psi) (+ (* (expt a 2) c phidot
	     psidot) (* -1 (expt b 2) c phidot psidot))))) (* (expt (sin psi)
	     2) (+ (* a (+ (* a c psidot thetadot) (* -2 b c psidot thetadot)))
	     (* (expt b 2) c psidot thetadot))) (* -1 b (expt c 2) psidot
	     thetadot))) (* (sin theta) (sin psi) (+ (* (expt (sin psi) 2) (+
	     (* (expt a 2) c phidot psidot) (* -1 (expt b 2) c phidot psidot)))
	     (* -1 a (expt c 2) phidot psidot) (* b (expt c 2) phidot
	     psidot))))) (* (expt (sin psi) 2) (+ (* (expt (sin psi) 2) a (+ (*
	     a c psidot thetadot) (* -1 b c psidot thetadot))) (* -1 a (expt c
	     2) psidot thetadot))))) (* (cos psi) (+ (* (cos psi) (+ (* (cos
	     psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* (cos psi) (expt (sin
	     theta) 2) a (+ (* a b phidot thetadot) (* -1 (expt b 2) phidot
	     thetadot))) (* (sin theta) (+ (* (expt (sin theta) 2) (sin psi) a
	     (+ (* -1 a b (expt phidot 2)) (* (expt b 2) (expt phidot 2)))) (*
	     (sin psi) a (+ (* a b (expt thetadot 2)) (* -1 (expt b 2) (expt
	     thetadot 2)))))))) (* (expt (sin theta) 2) (+ (* (expt (sin psi)
	     2) a (+ (* a b phidot thetadot) (* -1 (expt b 2) phidot
	     thetadot))) (* -1 a b c phidot thetadot))))) (* (sin theta) (+ (*
	     (expt (sin theta) 2) (expt (sin psi) 3) a (+ (* -2 a b (expt
	     phidot 2)) (* 2 (expt b 2) (expt phidot 2)))) (* (expt (sin psi)
	     3) a (+ (* 2 a b (expt thetadot 2)) (* -2 (expt b 2) (expt
	     thetadot 2)))))))) (* (expt (sin theta) 2) (expt (sin psi) 2) (+
	     (* (expt (sin psi) 2) a (+ (* -1 a b phidot thetadot) (* (expt b
	     2) phidot thetadot))) (* -2 a b c phidot thetadot))))) (* (sin
	     theta) (+ (* (expt (sin theta) 2) (expt (sin psi) 5) a (+ (* -1 a
	     b (expt phidot 2)) (* (expt b 2) (expt phidot 2)))) (* (expt (sin
	     psi) 5) a (+ (* a b (expt thetadot 2)) (* -1 (expt b 2) (expt
	     thetadot 2)))))))) (* (expt (sin theta) 2) (expt (sin psi) 4) (+
	     (* (expt (sin psi) 2) a (+ (* -1 a b phidot thetadot) (* (expt b
	     2) phidot thetadot))) (* -1 a b c phidot thetadot)))) (+ (* (expt
	     (cos psi) 2) (+ (* (expt (cos psi) 2) (sin theta) a b c) (* 2 (sin
	     theta) (expt (sin psi) 2) a b c))) (* (sin theta) (expt (sin psi)
	     4) a b c))) (/ (+ (* (cos theta) (+ (* (sin theta) (+ (* (expt
	     (cos psi) 2) a (+ (* a (expt phidot 2)) (* -1 c (expt phidot 2))))
	     (* (expt (sin psi) 2) b (+ (* b (expt phidot 2)) (* -1 c (expt
	     phidot 2)))))) (* (cos psi) (sin psi) (+ (* a (+ (* a phidot
	     thetadot) (* -1 c phidot thetadot))) (* b (+ (* -1 b phidot
	     thetadot) (* c phidot thetadot))))))) (* (sin theta) (+ (* (expt
	     (cos psi) 2) (+ (* (expt (cos psi) 2) a b phidot psidot) (* 2
	     (expt (sin psi) 2) a b phidot psidot) (* a (+ (* -1 a phidot
	     psidot) (* c phidot psidot))))) (* (expt (sin psi) 2) (+ (* (expt
	     (sin psi) 2) a b phidot psidot) (* b (+ (* -1 b phidot psidot) (*
	     c phidot psidot))))))) (* (cos psi) (sin psi) (+ (* a (+ (* -1 a
	     psidot thetadot) (* c psidot thetadot))) (* b (+ (* b psidot
	     thetadot) (* -1 c psidot thetadot)))))) (+ (* (expt (cos psi) 2)
	     (+ (* (expt (cos psi) 2) a b) (* 2 (expt (sin psi) 2) a b))) (*
	     (expt (sin psi) 4) a b))) (/ (+ (* (cos theta) (+ (* (cos psi) (+
	     (* (cos psi) (+ (* (cos psi) (+ (* (cos psi) (+ (* -1 a b phidot
	     thetadot) (* -1 (expt b 2) phidot thetadot))) (* (sin theta) (sin
	     psi) (+ (* -1 (expt a 2) (expt phidot 2)) (* (expt b 2) (expt
	     phidot 2)))))) (* (expt (sin psi) 2) (+ (* a (+ (* -1 a phidot
	     thetadot) (* -2 b phidot thetadot))) (* -1 (expt b 2) phidot
	     thetadot))) (* b c phidot thetadot))) (* (sin theta) (sin psi) (+
	     (* (expt (sin psi) 2) (+ (* -1 (expt a 2) (expt phidot 2)) (*
	     (expt b 2) (expt phidot 2)))) (* a c (expt phidot 2)) (* -1 b c
	     (expt phidot 2)))))) (* (expt (sin psi) 2) (+ (* (expt (sin psi)
	     2) a (+ (* -1 a phidot thetadot) (* -1 b phidot thetadot))) (* a c
	     phidot thetadot))))) (* (cos psi) (+ (* (cos psi) (+ (* (cos psi)
	     (+ (* (cos psi) (+ (* -1 a b psidot thetadot) (* (expt b 2) psidot
	     thetadot))) (* (sin theta) (sin psi) (+ (* (expt a 2) phidot
	     psidot) (* -1 (expt b 2) phidot psidot))))) (* (expt (sin psi) 2)
	     (+ (* a (+ (* a psidot thetadot) (* -2 b psidot thetadot))) (*
	     (expt b 2) psidot thetadot))) (* -1 b c psidot thetadot))) (* (sin
	     theta) (sin psi) (+ (* (expt (sin psi) 2) (+ (* (expt a 2) phidot
	     psidot) (* -1 (expt b 2) phidot psidot))) (* -1 a c phidot psidot)
	     (* b c phidot psidot))))) (* (expt (sin psi) 2) (+ (* (expt (sin
	     psi) 2) a (+ (* a psidot thetadot) (* -1 b psidot thetadot))) (*
	     -1 a c psidot thetadot)))) (+ (* (expt (cos psi) 2) (+ (* (expt
	     (cos psi) 2) (sin theta) a b) (* 2 (sin theta) (expt (sin psi) 2)
	     a b))) (* (sin theta) (expt (sin psi) 4) a b))))))
