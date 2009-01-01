;;; Figure out if ScmUtils is loaded (by checking if a few key procedures are
;;; defined):

(define *using-scmutils?*
  (let ((procs '(derivative vector:scalar*vector)))
    (and mit-scheme?
	 (not (memq #f (map
			(lambda (proc)
			  (environment-bound? (the-environment) proc))
			procs))))))

(newline)
(display
 (if *using-scmutils?*
     "*** It looks like we're running ScmUtils..."
     "*** ScmUtils doesn't seem to be running.  Going numerical..."))
(newline)

(let ((preload '("misc"
		 "lshared"))

      (core '("charts"
	      "manifold"
	      "vbundle"
	      "product"
	      "boundary"
	      "ranges"
	      "smooth"
	      "spaces"))

      (numdiff '("misc-math"
		 "fem/matlib"
		 "linear"
		 "richardson"))

      (scmutils '("stubs")))

  (for-each load preload)
  (for-each load (if *using-scmutils?* scmutils numdiff))
  (for-each load core))
