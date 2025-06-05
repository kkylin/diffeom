This is DIFFEOM ("DIFFerential Equations On Manifolds"), a
Scheme library for numerically solving differential
equations on manifolds.  It is the source code for MIT AI
Tech Memo 1631 (see below), which basically documents my
MEng thesis, written at Project MAC under the direction of
Gerald Jay Sussman.  The goal was to provide a software
framework for solving differential equations (ordinary and
partial) on manifolds, with a focus on those equations
arising from classical mechanics.

Main features:

1. Implementation of smooth manifolds as computational
   objects.  Basic manifolds are defined via charts;
   procedures are provided to compose simpler manifolds to
   form more complex ones.  For example, given a manifold,
   one can construct its tangent and cotangent bundles.
   Basic manifolds provided include circles, spheres, and
   the rotation group SO(3) (useful for modeling rigid body
   motion).

2. Numerical ODE solver on manifolds.

3. Numerical PDE solver on manifolds, using finite element
   methods.  These are applied to the Laplace and wave
   equations, the latter solved on spacetime meshes.  I also
   implemented a 2d Delaunay triangulation generator.

A few notes on various quirks:

1. MIT Scheme.  The code was written for MIT Scheme with the
   ScmUtils library.  I haven't ran this in many years, and
   it may require some tweaking to get working again.

2. Other Schemes.  At the time, I had a Mac SE/30 at home
   that only ran Gambit.  So I wrote stubs for the few
   ScmUtils functions I really relied on, e.g., DIFF.  But
   the replacements are not as reliable.

3. Duplications.  At the time, I really believed in writing
   my own code rather than using libraries and other
   standard tools.  For example, I had my own little script
   to compile all the code, rather than just write a
   Makefile.

4. Bugs.  There are many bugs, I'm sure.  One that I
   remember is the caching of pre-computed results wasn't
   done in a very systematic or reliable ways.  Users
   beware.

The original abstract is below.

I welcome questions, comments, what you use it for.

Kevin K Lin
lin1@arizona.edu
July 24, 2016
updated June 4, 2025

----------------------------------------------------------------------------

AIM-1631

Coordinate-Independent Computations on Differential
Equations

Author[s]: Kevin K. Lin

Date: March 1998

PS Download: ftp://publications.ai.mit.edu/ai-publications/1500-1999/AIM-1631.ps

PDF Download: ftp://publications.ai.mit.edu/ai-publications/pdf/AIM-1631.pdf

Abstract: This project investigates the computational
representation of differentiable manifolds, with the primary
goal of solving partial differential equations using
multiple coordinate systems on general n- dimensional
spaces. In the process, this abstraction is used to perform
accurate integrations of ordinary differential equations
using multiple coordinate systems. In the case of linear
partial differential equations, however, unexpected
difficulties arise even with the simplest equations.
