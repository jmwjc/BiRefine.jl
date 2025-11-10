using BiRefine
using Test

@testset "BiRefine.jl" begin
    function f(x,t)
        if x < t - 1
            return 2/π
        elseif x > t
            return 0.0
        else
            return (1-cos(π*(c*t - x)))/π
        end
    end
    birefine("msh/square.msh", f)
end
