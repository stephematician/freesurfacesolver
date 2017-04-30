Free surface continuation
=========================

_Stephen Wade_

This is the code for the thesis
[Very steep forced solitary waves in two-dimensional free surface flow](thesis_link).

The original plan was to have most of the computational work done in C++, with
the residual for the system of equations specified as a Lua script. This ended
up being a somewhat poor choice, as the residual function can occupy a lot of
computational time.

The code was use to produce all the figures in the thesis, and the figures for
the publication [On the free-surface flow of steep forced solitary waves](jfm_link).

Another publication based on the results in Chapter 4 of the thesis has been
accepted, with link/details pending.

Anyway, don't judge me too harshly, it was a big learning process. I have
updated the code to be slightly more memory friendly, and attempted to
consolidate a lot of the Lua scripts.

The MATLAB scripts are, in my opinion, a rushed together hack. I would not be
surprised if they break easily, and it is best to ignore them as they are only
helpers that I used to make some plots.

The basic directory structure

  - `include` - header files for C libraries required by Lua modules for free-
      surface problem.
  - `input` - scripts (Lua) batch solve free surface problem
  - `lua` - Lua modules for free surface problem
  - `matlab` - scripts (MATLAB) used to generate .tikz output
  - `src` - source files for C libraries required by Lua modules for free-
      surface problem.

## Instructions

Build requirements:

  - [GNU Scientific Library](gsl_link);
  - [Lua 5.2](lua_link); and
  - C++0x and GNU tool chain

### Build Lua libraries

The Lua libraries can be built using the GNU make tool chain

  - on Ubuntu 14.04 LTS (tested)
  ```
  make linux
  ```
  - Max OS-X with MacPorts
  ```
  make macosx
  ```

The macosx build has not been tested since 2011, and is definitely broken. For
example, build commands are not updated to reflect C++0x requirements now in
place.

### Generate output

All the main scripts are run from the command line, e.g.

```
lua input/htopography_F110_L100.lua
```

This will generate output in the `output` directory. Most scripts will try to
create one or more branches of free surfaces in the amplitude of forcing-free
surface height plane, for either fixed values of the Froude number, or fixed
values of the speed at the crest.

The various `matlab` scripts supplied can be used to generate .tikz output as
per the publications and thesis. These scripts are _not_ Octave compatible.

[thesis_link]: http://hdl.handle.net/2440/97905 
[jfm_link]: https://doi.org/10.1017/jfm.2013.590
[gsl_link]: https://www.gnu.org/software/gsl/
[lua_link]: https://www.lua.org/manual/5.2/
