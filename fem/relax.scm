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

;;; Matrix inversion using successive overrelaxation methods (SOR):
;;;
;;; We use classical SOR methods to solve linear systems of equations.
;;; The representation for matrices is defined in sparse.scm.

(declare (usual-integrations))
;(load "sparse")


;;; Note that this procedure first modifies the matrix by dividing through by
;;; the diagonal...

(define (sor sm n . aux)
  (let* ((nrows (sparse-matrix-row-count sm))
         (ncols-1 (- (sparse-matrix-column-count sm) 1))
         (rhs (make-vector nrows))
         (state (make-vector nrows 0))
         (factor 1)
         (residual 0))

    ;; Parse auxiliary arguments, if any:
    ;; (The first should be SOR factor, the second an alternate state.)

    (if (not (null? aux))
        (if (not (null? (cdr aux)))
            (set! state (cadr aux))
            (set! factor (car aux))))

    ;; Normalize the matrix:

    (do ((i 0 (+ i 1)))
        ((>= i nrows))
      (let ((diag (sparse-matrix-ref sm i i)))
        (for-each
         (lambda (pair)
           (sparse-matrix-set! sm i (car pair) (/ (cadr pair) diag)))
         (sparse-matrix-get-row sm i))))

    ;; Collect the right hand side of the equation Ax = b:

    (do ((i 0 (+ i 1)))
        ((>= i nrows))
      (vector-set! rhs i (sparse-matrix-ref sm i ncols-1)))

    ;; Now perform SOR (note that we have normalized the matrix so that the
    ;; diagonal terms are all unity):

    (do ((n n (- n 1)))
        ((<= n 0))

      (set! residual 0.)

      (do ((i 0 (+ i 1)))
          ((>= i nrows))

        ;; Compute the row sum:

        (let ((sum 0))
          (for-each
           (lambda (pair)
             (let ((index (car pair))
                   (val (cadr pair)))
               (if (< index ncols-1)
                   (set! sum (+ sum (* (vector-ref state index) val))))))
           (sparse-matrix-get-row sm i))

          ;; Step forward:

          (let ((step (- (vector-ref rhs i) sum)))
            (if (> (abs step) residual)
                (set! residual (abs step)))

            (vector-set! state i (+ (vector-ref state i) (* factor step)))))))

    (write-line `(residual: ,residual))
    state))
