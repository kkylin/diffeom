(load "load-ode")

(show-time
 (lambda ()

   (let ((segment-size 1.)
	 (count 101)
	 (filename "rigid-man.data")
	 (charts (manifold:get-finite-atlas TSO3)))

     (let loop ((i 0) (x singular-init) (t 0.))
       (if (< i count)
	   (begin

	     (write-line `(step ,i t = ,t))

	     (let ((results
		    (show-time
		     (lambda ()
		       (rigid-body-path x (+ t segment-size) t)))))

	       (let ((port (open-output-file filename #t)))

		 (for-each

		  (lambda (l)
		    (let ((t (car l))
			  (p (cadr l)))

		      (display t port)

		      (let loop ((i 0) (charts charts))
			(if (null? charts)
			    (display "No chart!" port)
			    (let ((chart (car charts)))
			      (if (chart:member? p chart)
				  (begin
				    (display " " port)
				    (display i port)
				    (display " " port)
				    (display (chart:point->coords p chart)
					     port))
				  (loop (+ i 1) (cdr charts))))))

		      (newline port)))

		  (let ((l (sort results (lambda (x y) (< (car x) (car y))))))
		    (if (> t 0)
			(cdr l)
			l)))

		 (close-output-port port))

	       (loop (+ i 1) (cadar results) (caar results)))))))))
