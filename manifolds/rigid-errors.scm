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

(define (check-relative-error list)
  (let ((ref (car list)))
    (write-line `(reference value = ,ref))
    (map (lambda (val) (relative-error val ref)) list)))

(define (check-absolute-error list)
  (let ((ref (car list)))
    (write-line `(reference value = ,ref))
    (map (lambda (val) (abs (- val ref))) list)))

#|
(define reg-list
  (show-time
   (lambda ()
     (read-regular-file "rigid-reg.data"))))

(define reg-e-list
  (show-time
   (lambda ()
     (map (t-rigid-body 1. (sqrt 2.) 2.) regular))))

(define reg-l-list
  (show-time
   (lambda ()
     (map (state->L-space 1. (sqrt 2) 2.) regular))))

(define reg-l1-list (map vector-first reg-l-list))
(define reg-l2-list (map vector-second reg-l-list))
(define reg-l3-list (map vector-third reg-l-list))

(define reg-e-errors (check-absolute-error reg-e-list))
(define reg-l1-errors (check-absolute-error reg-l1-list))
(define reg-l2-errors (check-absolute-error reg-l2-list))
(define reg-l3-errors (check-absolute-error reg-l3-list))

(define man-list
  (show-time
   (lambda ()
     (read-manifold-file "rigid-man.data"))))

(define man-e-list
  (show-time
   (lambda ()
     (map (t-rigid-body 1. (sqrt 2) 2.) manifold))))

(define man-l-list
  (show-time
   (lambda ()
     (map (state->L-space 1. (sqrt 2) 2.) man-list))))

(define man-l1-list (map vector-first man-l-list))
(define man-l2-list (map vector-second man-l-list))
(define man-l3-list (map vector-third man-l-list))

(define man-e-errors (check-absolute-error man-e-list))
(define man-l1-errors (check-absolute-error man-l1-list))
(define man-l2-errors (check-absolute-error man-l2-list))
(define man-l3-errors (check-absolute-error man-l3-list))
|#
