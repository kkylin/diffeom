;;; Load the PDE stuff separately (so that the parts of the program *not* under
;;; development can still be used without these definitions).

(let ((pde-core
       '("pde-aux"
	 "pde-charts"
	 "pde-cmpgrd"
	 "pde-debug"
	 "pde-elements"
	 "pde-main"
	 "pde-mergers"
	 "pde-nodes"
	 "pde-ops"
	 "pde-tools"
	 "pde-config"

	 "basis-imb"
	 "basis-poly"
	 "basis-real"

	 "pde-examples"))

      (fem-stuff
       '("basis"
	 "2d-poly-basis"
	 "2d-real-basis"
	 "2d-trapezoid"
	 "delaunay"
	 "delaux"
	 "dyntable"
	 "edge"
	 "fem"
	 "relax"
	 "sparse"
	 "util-too"))

      (fem-dir "fem/"))

  (load "load-main")

  (if *using-scmutils?*
      (begin
	(newline)
	(display "*** Warning: PDE code does not work well with ScmUtils!")
	(newline)))

  (for-each (lambda (file)
	      (load (string-append fem-dir file)))
	    fem-stuff)

  (for-each load pde-core))
