local deepcopy
local create_file
local append_file
local get_function_name

-- Taken from lua-users wiki; uses a lookup-table/closure to ensure that tables
-- that appear multiple times are only duplicated once
deepcopy = function(object)
    local lookup_table = {}

    local function _copy(_object)
        if type(_object) ~= "table" then
            return _object
        elseif lookup_table[_object] then
            return lookup_table[_object]
        end

        local lt = {}
        lookup_table[_object] = lt

        for k, v in pairs(_object) do
            lt[_copy(k)] = _copy(v)
        end

        return setmetatable(lt, getmetatable(_object))
    end

    return _copy(object)

end

get_function_name = function(file)

    local ext_j = file:reverse():find(".", 1, true)
    if ext_j then
        ext_j = 1 + file:len() - ext_j
    else
        ext_j = file:len()
    end

    local path_j = file:reverse():find("/", 1, true)
    if path_j then
        path_j = 1 + file:len() - path_j
    else
        path_j = 0
    end

    return file:sub(path_j+1, ext_j-1)

end

create_file = function(file)

    local fn_name = get_function_name(file)

    local of = io.open(file, "w")
    of:write(string.format("%% %08d", 0) .. "\n\n")
    of:write("function [ o ] = " .. fn_name .. "()\n\n")
    of:write("\nend\n")
    of:close()

end

append_file = function(file,
                       output)

    -- Write down the new number of surfaces contained in file
    local of = io.open(file, "r+")

    local fn_name = get_function_name(file)
    
    of:seek("set", string.len("%% "))
    local n = of:read("*n") + 1
    
    of:seek("set", 0)
    of:write(string.format("%% %08d", n) .. "\n\n")
    of:write("function [ o ] = " .. fn_name .. "()\n\n")
    of:close()

    -- Write the new surface data into file
    of = io.open(file, "r+")
    of:seek("end", -string.len("end\n"))

    for k, v in pairs(output) do
        if type(v) == "number" then
            of:write("  o." .. k .. "{" .. n .. "} = " .. v .. ";\n")
        elseif type(v) == "table" then
            of:write("  o." .. k .. "{" .. n .. "} = [")
            for kk, vv in ipairs(v) do
                of:write(vv .. " ")
            end
        of:write("];\n")
        else
            error([=[Output computed of surface should be vectors and ]=] ..
                  [=[scalars only.]=])
        end
    end

    of:write("\nend\n")

    of:close()

end

return {
    deepcopy = deepcopy,
    create_file = create_file,
    append_file = append_file
}
