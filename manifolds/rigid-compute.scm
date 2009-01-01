(load "rigid")


;;; Compute some expressions:

(for-each
 (lambda (make-sysder filename)
   (let ((port (open-output-file filename)))
     (pp (traditional->correct-order
	  (vector-tail
	   (show-time
	    (lambda ()
	      (*sysder-simplify*
	       ((make-sysder 'a 'b 'c) rigid-qqdot))))
	   1))
	 port)
     (close-output-port port)))
 (list rigid-sysder-0 rigid-sysder-1 rigid-sysder-2 rigid-sysder-3)
 (list "rigid-field-0" "rigid-field-1" "rigid-field-2" "rigid-field-3"))
