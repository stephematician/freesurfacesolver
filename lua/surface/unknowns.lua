local pack
local expand

--[[
-- This function selects the 'free' parts of a vector and arranges them in
-- sequence, placing the continuation variable last in
-- the final vector
--]]
pack = function(current, is_free, is_continued)

    local packed = {};
    local index_map = {}

    local continued_value = nil;
    local continued_found = false;

    i = 1;
    for k, k_value in pairs(current) do

        if type(k_value) == "number" then

            if is_free[k] then
                if not is_continued[k] then
                    packed[i] = k_value
                    index_map[k] = i
                    i = i + 1;
                else -- continued
                    if continued_found then
                        error(
                            [=[Too many variables selected for continuation.]=]
                        )
                    end
                    continued_found = true
                    continued_value = k_value
                    index_map[k] = 0
                end
            else -- not free
                if is_continued[k] then
                    error([=[Non-free variable selected for continuation.]=])
                end
            end

        elseif type(k_value) == "table" then

            index_map[k] = {}

            for kk, kk_value in ipairs(k_value) do

                if type(kk_value) == "number" then
                    if is_free[k][kk] then
                        if not is_continued[k][kk] then
                            packed[i] = kk_value
                            index_map[k][kk] = i
                            i = i + 1
                        else -- continued
                            if continued_found then
                                error(
                                    [=[Too many variables selected for ]=] ..
                                    [=[continuation.]=]
                                )
                            end
                            continued_found = true
                            continued_value = kk_value
                            index_map[k][kk] = 0
                        end
                    else -- not free
                        if (is_continued[k])[kk] then
                            error(
                                [=[Non-free variable selected for ]=] ..
                                [=[continuation.]=]
                            )
                        end
                    end
                else -- type value not == number
                    error(
                        [=[Expected the unpacked surface description to ]=] ..
                        [=[be a table containing numbers and tables of ]=] ..
                        [=[numbers.]=]
                    )
                end

            end

        else -- type of cvec not == number or table 
            error(
                [=[Expected the unpacked surface description to ]=] ..
                [=[be a table.]=]
            )
        end
    
    end
  
    packed[i] = continued_value

    for k, k_value in pairs(index_map) do
        if type(k_value) == "table" then
            for kk, kk_value in ipairs(k_value) do
                if kk_value == 0 then
                    index_map[k][kk] = i
                end
            end
        else
            if k_value == 0 then
                index_map[k] = i
            end
        end
    end

    return packed, index_map 
  
end

--[[
-- Creates a surface description from the packed 'unknowns' (free variables) and
-- inserts initial values for non-free variables.
-- Does not perform any checks on the arguments, requires a valid packed,
-- is_free and index_map from a surface obejct.
--]]
expand = function(packed, initial, is_free, index_map)

    local current = {}

    for k, k_value in pairs(initial) do

        if type(k_value) == "number" then
    
            if is_free[k] then
                current[k] = packed[index_map[k]]
            else 
                current[k] = k_value;
            end

        else
 
            current[k] = {}

            for kk, kk_value in ipairs(k_value) do

                if is_free[k][kk] then
                    current[k][kk] = packed[index_map[k][kk]]
                else -- not free
                    current[k][kk] = kk_value;
                end

            end

        end

    end

    return current

end

return {
    pack = pack,
    expand = expand
}
