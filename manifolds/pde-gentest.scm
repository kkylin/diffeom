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

;;; Useful for generating the test cases in pde-thesis.scm:

(define (generate-test-case ilist rectangular spherical count sor-fact string)
  (for-each
   (lambda (i)
     (newline)
     (display
      (string-append
       "(test-" (number->string i) " "
       "'(rectangular " (number->string (cadr rectangular)) " "
       (number->string (caddr rectangular)) ") "
       "'(spherical " (number->string (cadr spherical)) " "
       (number->string (caddr spherical)) ") "
       (number->string count) " "
       (number->string sor-fact) " "
       "\"Data/thesis/test" (number->string i) string "\")")))
   ilist)
  (newline))


;;; Generate the test cases:

(define (generate-test-set indices)
  (newline)
  (display "*** Test cases:")
  (newline)
  (for-each
   (lambda (args)
     (apply generate-test-case (cons indices args)))
   (list (list '(rectangular 10 5) '(spherical 5 10) 10000 1.9 "a")
	 (list '(rectangular 14 7) '(spherical 7 14) 11000 1.9 "b")
	 (list '(rectangular 18 9) '(spherical 9 18) 12000 1.9 "c")
	 (list '(rectangular 22 11) '(spherical 11 22) 13000 1.9 "d")
	 (list '(rectangular 26 13) '(spherical 13 26) 14000 1.9 "e")
	 (list '(rectangular 30 15) '(spherical 15 30) 15000 1.9 "f")
	 (list '(rectangular 34 17) '(spherical 17 34) 16000 1.9 "g")
	 (list '(rectangular 38 19) '(spherical 19 38) 17000 1.9 "h")
	 (list '(rectangular 42 21) '(spherical 21 42) 18000 1.9 "i")
	 (list '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "j")
	 (list '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "k"))))
