--[[
-- Surface object
--
-- A surface object is an interface to the continuator module with additional
-- functions for progressing along a continuation curve in a given 'domain'
-- (e.g. the amplitude of forcing/free-surface height domain), or for
-- approaching a desired limit in the continuation variable.
--
-- Author        : Stephen Wade
-- Date created  : ?
-- Last modified : 13/05/2017 
--]]


local surface = {}


surface.new = require 'surface.new'


surface.display = function(prefix, s, elements)

    if not prefix then
        prefix = ''
    end

    _unknowns = s.expand()

    for k = 1, #elements do
        io.write(prefix,
                 '  - ', elements[k], '; ', _unknowns[elements[k]], '\n')
    end

end


return surface

