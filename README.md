Free surface continuation
=========================

_Stephen Wade_

Directory structure

  - `include` - header files for C libraries required by Lua modules for free-
      surface problem.
  - `input` - scripts (Lua) batch solve free surface problem
  - `lua` - Lua modules for free surface problem
  - `matlab` - scripts (MATLAB) used to generate .tikz output
  - `src` - source files for C libraries required by Lua modules for free-
      surface problem.

## Instructions

Build requirements:
  - GNU Scientific Library and
  - Lua 5.2

The C libraries can be built using the GNU make tool chain

  - on Ubuntu 14.04 LTS (tested)
  ```
  make linux
  ```
  - Max OS-X with MacPorts
  ```
  make macosx
  ```

Then all the scripts can be run e.g.

```
lua input/htopography_F110_L100.lua
```

This will generate output in the `output` directory

Various `matlab` scripts can then be used to generate .tikz output as per the
publications and thesis.


