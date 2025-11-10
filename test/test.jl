using BiRefine

function f(x,t)
    if x < t - 1
        return 2/π
    elseif x > t
        return 0.0
    else
        return (1-cos(π*(t - x)))/π
    end
end

birefine("msh/square.msh", f, order=10, tol=1e-13, maxiter=4)