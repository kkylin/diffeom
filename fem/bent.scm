#|

This file is part of DIFFEOM ("DIFFerential Equations
On Manifolds"), a system for solving differential
equations on manifolds.

Copyright (C) 2016 by Kevin K Lin
<klin@math.arizona.edu>

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

;;; Drawing the bent characteristics:


;;; First, open graphics device and set scale:

(define dev (make-graphics-device 'x))
(graphics-set-coordinate-limits dev -.1 -.1 1.1 1.1)

(if #f
    (begin
      (graphics-operation dev 'set-background-color "black")
      (graphics-operation dev 'set-foreground-color "red")
      (graphics-operation dev 'set-mouse-color "white")))


;;; Next, define a procedure to integrate and draw the characteristics (using a
;;; silly foward-Euler integrator).

(define (draw-characteristic slope x0 t0 dt)
  (graphics-move-cursor dev x0 t0)

  (let loop ((i 1) (x x0))
    (let ((t (+ (* i dt) t0)))
      (if (and (<= 0 t) (<= t 1)
	       (<= 0 x) (<= x 1))
	  (let* ((dt/dx (slope t))
		 (new-x (+ x (/ dt dt/dx))))
	    (cond ((> new-x 1) (graphics-drag-cursor dev 1 t))
		  ((< new-x 0) (graphics-drag-cursor dev 0 t))
		  (else (graphics-drag-cursor dev new-x t)))
	    (loop (+ i 1) new-x)))))
  'done)


;;; The functions we want to use:

(define (cut-off t)
  (cond ((<= t 0) 1.)
	((>= t 1) 0.)
	(else (+ (* 2 (expt t 3)) (* -3 (expt t 2)) 1))))

(define (f t)
  (cut-off (- (* 2 t) 1)))

(define (g t)
  (- (f t)))


;;; Something to generate sample points:

(define (samples min max count)
  (let ((dx (/ (- max min) (- count 1))))
    (let loop ((result '()) (i (- count 1)))
      (if (< i 0)
	  result
	  (loop (cons (+ (* i dx) min) result) (- i 1))))))

(define (constant-list val count)
  (vector->list (make-vector count val)))


;;; Do it!

(graphics-clear dev)

(graphics-draw-line dev 0. 0. 0. 1.)
(graphics-draw-line dev 0. 0. 1. 0.)
(graphics-draw-line dev 0. 0. 1. 0.)
(graphics-draw-line dev 1. 0. 1. 1.)
(graphics-draw-line dev 0. 1. 1. 1.)

(for-each
 (lambda (x t)
   (draw-characteristic f x t .01))
 (samples 0. 1. 11)
 (constant-list 0. 11))

(for-each
 (lambda (x t)
   (draw-characteristic g x t .01))
 (samples 0. 1. 11)
 (constant-list 0. 11))

(for-each
 (lambda (x t)
   (draw-characteristic f x t .01))
 (constant-list 0. 11)
 (samples 0. .9 11))

(for-each
 (lambda (x t)
   (draw-characteristic g x t .01))
 (constant-list 1. 11)
 (samples 0. .9 11))

;(graphics-close dev)
