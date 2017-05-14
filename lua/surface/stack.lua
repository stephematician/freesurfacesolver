--[[
-- Basic (default) copy-semantics stack which silently drops overflow elements
--]]

local stack = {}

local new

new = function(n)

    --[[--------------------------------------------------------------
    --                                         Forward declarations --
    ----------------------------------------------------------------]]
    local get_count, push, peek, pop, flush, copy


    --[[--------------------------------------------------------------
    --                                  Private data initialisation --
    ----------------------------------------------------------------]]
    local self = {
        size = n,
        top = 0,
        count = 0,
        elements = {}
    }


    --[[--------------------------------------------------------------
    --                                                     Closures --
    ----------------------------------------------------------------]]
    local get_count = function()
        return self.count
    end


    local push = function(v)
        -- silently drop the oldest element on the stack
        if self.count < self.size then
            self.count = self.count + 1
        end
        self.top = (self.top + 1) % self.size
        self.elements[self.top] = v
    end


    local peek = function(k)
        -- need to determine if index is valid
        if (k > self.count) or (-k > self.count) then
            return nil
        end

        if k > 0 then
            return self.elements[(self.top-(k-1)) % self.size]
        elseif k < 0 then
            return self.elements[(self.top-(self.count+k)) % self.size]
        else
            return nil
        end

    end


    local pop = function()
        if self.count == 0 then
            return nil
        else
            self.top = (self.top - 1) % self.size
            self.count = self.count - 1
            return self.elements[(self.top + 1) % self.size]
        end
    end


    local flush = function()
        for k = 1, #self.elements do self.elements[k] = nil end
        self.count = 0
    end


    local copy = function()
        _stack = new(self.size)
        for k = 1, self.get_count() do
            _stack.push(peek(-k))
        end
        return _stack
    end


    return {
        get_count = get_count,
        push = push,
        peek = peek,
        pop = pop,
        flush = flush,
        copy = copy
   }

end


return {
    new = new
}

