module BiRefine

import Gmsh: gmsh
export birefine

function birefine(filename::String, f::Function;order::Int=2,tol::Float64=1e-6,maxiter::Int=5)
    gmsh.initialize()
    gmsh.open(filename)

    removedElementTags = Set{UInt}([1])
    iter = 0
    while !isempty(removedElementTags) && iter < maxiter
        iter += 1
        for (elementTypes, elementTags, nodeTags, tag) in getElements2D()
            if elementTypes[1] == 2
                empty!(removedElementTags)
                elementTags = elementTags[1]
                nodeTags = reshape(nodeTags[1], (3, :))
                localCoord, weights = gmsh.model.mesh.getIntegrationPoints(2, "Gauss"*string(order))
                localCoord = reshape(localCoord, (3, :))
                for (elementTag, nodeTag) in zip(elementTags, eachcol(nodeTags))
                    if elementTag ∉ removedElementTags && issplit(nodeTag, localCoord, weights, f, tol)
                        𝐱₁, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[1])
                        𝐱₂, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[2])
                        𝐱₃, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[3])
                        longestTag, ~ = getLongestEdge(𝐱₁, 𝐱₂, 𝐱₃)
                        𝐱₁, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[longestTag[1]])
                        𝐱₂, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[longestTag[2]])
                        𝐱₄ = 0.5.*(𝐱₁ .+ 𝐱₂)
                        nt4 = gmsh.model.mesh.getMaxNodeTag() + 1
                        gmsh.model.mesh.addNodes(2, tag, [nt4], 𝐱₄)
                        splitElement2D(elementTag, nodeTag, longestTag, nt4, removedElementTags, tag)
                    end
                end
            end
        end
    end

    filename_components = split(filename, ".")
    gmsh.write(filename_components[1]*"_refined."*filename_components[2])
    gmsh.finalize()
end

function getElements1D()
    elements1D = []
    dim_tags = gmsh.model.getEntities()
    for (dim, tag) in dim_tags
        if dim == 1  # Assuming we are refining 2D elements
            elements = gmsh.model.mesh.getElements(dim, tag)
            push!(elements1D, (elements..., tag))
        end
    end
    return elements1D
end

function getElements2D()
    elements2D = []
    dim_tags = gmsh.model.getEntities()
    for (dim, tag) in dim_tags
        if dim == 2  # Assuming we are refining 1D elements
            elements = gmsh.model.mesh.getElements(dim, tag)
            push!(elements2D, (elements..., tag))
        end
    end
    return elements2D
end

function getLongestEdge(𝐱₁::Vector{Float64}, 𝐱₂::Vector{Float64}, 𝐱₃::Vector{Float64})
    L₁ = sqrt((𝐱₃[1]-𝐱₂[1])^2 + (𝐱₃[2]-𝐱₂[2])^2 + (𝐱₃[3]-𝐱₂[3])^2)
    L₂ = sqrt((𝐱₁[1]-𝐱₃[1])^2 + (𝐱₁[2]-𝐱₃[2])^2 + (𝐱₁[3]-𝐱₃[3])^2)
    L₃ = sqrt((𝐱₂[1]-𝐱₁[1])^2 + (𝐱₂[2]-𝐱₁[2])^2 + (𝐱₂[3]-𝐱₁[3])^2)
    if L₁ >= L₂ && L₁ >= L₃
        return (2,3), L₁
    elseif L₂ >= L₁ && L₂ >= L₃
        return (3,1), L₂
    else
        return (1,2), L₃
    end
end

function getPairSection(sectionTag::Set{UInt}, addNodeTag::UInt)
    for (elementTypes1D, elementTags1D, nodeTags1D, tag1D) in getElements1D()
        if elementTypes1D[1] == 1
            elementTags1D = elementTags1D[1]
            nodeTags1D = reshape(nodeTags1D[1], (2, :))
            for (elementTag1D, nodeTag1D) in zip(elementTags1D, eachcol(nodeTags1D))
                if sectionTag ⊆ nodeTag1D
                    gmsh.model.mesh.removeElements(1, tag1D,[elementTag1D])
                    nt1 = nodeTag1D[1]
                    nt2 = nodeTag1D[2]
                    nt3 = addNodeTag
                    gmsh.model.mesh.addElements(1, tag1D, [1], [[gmsh.model.mesh.getMaxElementTag() + 1, gmsh.model.mesh.getMaxElementTag() + 2]], [[nt1, nt3, nt3, nt2]])
                    return false
                end
            end
        end
    end
    return true
end

function getPairTriangle(elementTag_::UInt, nodeTag_::Tuple{UInt, UInt})
    for (elementTypes, elementTags, nodeTags, tag) in getElements2D()
        if elementTypes[1] == 2
            elementTags = elementTags[1]
            nodeTags = reshape(nodeTags[1], (3, :))
            for (elementTag, nodeTag) in zip(elementTags, eachcol(nodeTags))
                if nodeTag_ ⊆ nodeTag && elementTag ≠ elementTag_
                    sectionTag = (indexin(nodeTag_[2], nodeTag)[1], indexin(nodeTag_[1], nodeTag)[1])
                    return elementTag, nodeTag, sectionTag
                end
            end
        end
        # nodeTag_ = Int.(nodeTag_)
        # error("Cannot find the pair triangle element for elementTag $elementTag_, nodeTag $nodeTag_")
    end
end

function issplit(nodeTag::SubArray{UInt}, localCoord::Matrix{Float64}, weights::Vector{Float64}, f::Function, tol::Float64)
    𝐱₁, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[1])
    𝐱₂, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[2])
    𝐱₃, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[3])
    val = 0.0
    𝐴 = 0.5*(𝐱₁[1]*𝐱₂[2]+𝐱₂[1]*𝐱₃[2]+𝐱₃[1]*𝐱₁[2]-𝐱₁[2]*𝐱₂[1]-𝐱₂[2]*𝐱₃[1]-𝐱₃[2]*𝐱₁[1])
    for ((ξ, η, ~), w) in zip(eachcol(localCoord), weights)
        N₁ = 1.0 - ξ - η
        N₂ = ξ
        N₃ = η
        x = N₁*𝐱₁[1] + N₂*𝐱₂[1] + N₃*𝐱₃[1]
        y = N₁*𝐱₁[2] + N₂*𝐱₂[2] + N₃*𝐱₃[2]
        val += f(x,y) * w * 𝐴
    end

    longestTag, ~ = getLongestEdge(𝐱₁, 𝐱₂, 𝐱₃)
    peakTag = setdiff([1,2,3],longestTag)[1]
    𝐱₁, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[peakTag])
    𝐱₂, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[longestTag[1]])
    𝐱₃, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[longestTag[2]])
    𝐱₄ = 0.5.*(𝐱₂ .+ 𝐱₃)
    val_ = 0.0
    𝐴 /= 2
    for ((ξ, η, ~), w) in zip(eachcol(localCoord), weights)
        N₁ = 1.0 - ξ - η
        N₂ = ξ
        N₃ = η
        x = N₁*𝐱₁[1] + N₂*𝐱₂[1] + N₃*𝐱₄[1]
        y = N₁*𝐱₁[2] + N₂*𝐱₂[2] + N₃*𝐱₄[2]
        val_ += f(x,y) * w * 𝐴
        x = N₁*𝐱₁[1] + N₂*𝐱₄[1] + N₃*𝐱₃[1]
        y = N₁*𝐱₁[2] + N₂*𝐱₄[2] + N₃*𝐱₃[2]
        val_ += f(x,y) * w * 𝐴
    end
    return abs(val-val_)/abs(val) > tol ? true : false
end

function splitElement2D(elementTag::UInt, nodeTag::SubArray{UInt}, sectionTag::Tuple{Int,Int}, addNodeTag::UInt, removedElementTags::Set{UInt}, tag::Int32)
    push!(removedElementTags, elementTag)
    gmsh.model.mesh.removeElements(2, tag,[elementTag])

    nt1 = nodeTag[setdiff([1,2,3],sectionTag)...]
    nt2 = nodeTag[sectionTag[1]]
    nt3 = nodeTag[sectionTag[2]]
    nt4 = addNodeTag
    𝐱₁, ~, ~, ~ = gmsh.model.mesh.getNode(nt1)
    𝐱₂, ~, ~, ~ = gmsh.model.mesh.getNode(nt2)
    𝐱₃, ~, ~, ~ = gmsh.model.mesh.getNode(nt3)
    𝐱₄, ~, ~, ~ = gmsh.model.mesh.getNode(nt4)
    longestTag, ~ = getLongestEdge(𝐱₁, 𝐱₂, 𝐱₃)

    if longestTag ≠ (2,3)
        nt5 = gmsh.model.mesh.getMaxNodeTag() + 1
        𝐱₄ = 0.5.*(𝐱₂ .+ 𝐱₃)
        if longestTag == (1,2)
            𝐱₅ = 0.5.*(𝐱₁ .+ 𝐱₂)
            gmsh.model.mesh.addNodes(2, tag, [nt5], 𝐱₅)
            gmsh.model.mesh.addElements(2, tag, [2], [[(gmsh.model.mesh.getMaxElementTag() + 1:gmsh.model.mesh.getMaxElementTag() + 3)...]], [[nt1,nt5,nt3,nt2,nt4,nt5,nt3,nt5,nt4]])
            if getPairSection(Set{UInt}([nt1,nt2]), nt5)
                if getPairTriangle(elementTag, (nt1,nt2)) ≠ nothing
                    elementTag_, nodeTag_, sectionTag_ = getPairTriangle(elementTag, (nt1,nt2))
                    if elementTag_ ∉ removedElementTags
                        splitElement2D(elementTag_, nodeTag_, sectionTag_, nt5, removedElementTags, tag)
                    end
                end
            end
        else
            𝐱₅ = 0.5.*(𝐱₁ .+ 𝐱₃)
            gmsh.model.mesh.addNodes(2, tag, [nt5], 𝐱₅)
            gmsh.model.mesh.addElements(2, tag, [2], [[(gmsh.model.mesh.getMaxElementTag() + 1:gmsh.model.mesh.getMaxElementTag() + 3)...]], [[nt1,nt2,nt5,nt2,nt4,nt5,nt3,nt5,nt4]])
            if getPairSection(Set{UInt}([nt3,nt1]), nt5)
                if getPairTriangle(elementTag, (nt3,nt1)) ≠ nothing
                    elementTag_, nodeTag_, sectionTag_ = getPairTriangle(elementTag, (nt3,nt1))
                    if elementTag_ ∉ removedElementTags
                        splitElement2D(elementTag_, nodeTag_, sectionTag_, nt5, removedElementTags, tag)
                    end
                end
            end
        end
    else
        gmsh.model.mesh.addElements(2, tag, [2], [[(gmsh.model.mesh.getMaxElementTag() + 1:gmsh.model.mesh.getMaxElementTag() + 2)...]], [[nt1,nt2,nt4,nt1,nt4,nt3]])
    end
    if getPairSection(Set{UInt}([nt2,nt3]), nt4)
        if getPairTriangle(elementTag, (nt2,nt3)) ≠ nothing
        elementTag_, nodeTag_, sectionTag_ = getPairTriangle(elementTag, (nt2,nt3))
        if elementTag_ ∉ removedElementTags
            splitElement2D(elementTag_, nodeTag_, sectionTag_, nt4, removedElementTags, tag)
        end
        end
    end
end

end
