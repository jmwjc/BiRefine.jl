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
        elements1D = []
        elements2D = []
        dim_tags = gmsh.model.getEntities()
        for (dim, tag) in dim_tags
            if dim == 1  # Assuming we are refining 2D elements
                elements = gmsh.model.mesh.getElements(dim, tag)
                push!(elements1D, (elements..., tag))
            elseif dim == 2  # Assuming we are refining 1D elements
                elements = gmsh.model.mesh.getElements(dim, tag)
                push!(elements2D, (elements..., tag))
            end
        end

        for (elementTypes, elementTags, nodeTags, tag) in elements2D
            if elementTypes[1] == 2
                empty!(removedElementTags)
                elementTags = elementTags[1]
                nodeTags = reshape(nodeTags[1], (3, :))
                localCoord, weights = gmsh.model.mesh.getIntegrationPoints(2, "Gauss"*string(order))
                localCoord = reshape(localCoord, (3, :))
                for (elementTag, nodeTag) in zip(elementTags, eachcol(nodeTags))
                    if elementTag ∉ removedElementTags
                        isPairOnBoundary = false
                        𝐱₁, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[1])
                        𝐱₂, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[2])
                        𝐱₃, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[3])
                        longestTag, L = getLongestEdge(𝐱₁, 𝐱₂, 𝐱₃)
                        val = calculate_integration([(𝐱₁, 𝐱₂, 𝐱₃)], localCoord, weights, f)

                        for (elementTypes1D, elementTags1D, nodeTags1D, tag1D) in elements1D
                            if elementTypes1D[1] == 1
                                elementTags1D = elementTags1D[1]
                                nodeTags1D = reshape(nodeTags1D[1], (2, :))
                                for (elementTag1D, nodeTag1D) in zip(elementTags1D, eachcol(nodeTags1D))
                                    if nodeTag[collect(longestTag)] ⊆ nodeTag1D
                                        isPairOnBoundary = true
                                        𝒙 = (𝐱₁,𝐱₂,𝐱₃)
                                        𝐱̃₁ = 𝒙[longestTag[1]]
                                        𝐱̃₂ = 𝒙[longestTag[2]]
                                        𝐱̃₃ = 𝒙[setdiff([1,2,3],longestTag)...]
                                        𝐱̃₄ = 0.5.*[
                                            𝐱̃₁[1]+𝐱̃₂[1],
                                            𝐱̃₁[2]+𝐱̃₂[2],
                                            𝐱̃₁[3]+𝐱̃₂[3]
                                        ]
                                        val_ = calculate_integration([(𝐱̃₁, 𝐱̃₄, 𝐱̃₃),(𝐱̃₄,𝐱̃₂, 𝐱̃₃)], localCoord, weights, f)
                                        if abs(val-val_)/abs(val) > tol
                                            push!(removedElementTags, elementTag)
                                            gmsh.model.mesh.removeElements(2, tag,[elementTag])
                                            gmsh.model.mesh.removeElements(1, tag1D,[elementTag1D])
                                            nt1, nt2 = nodeTag[collect(longestTag)]
                                            nt3 = nodeTag[setdiff([1,2,3],longestTag)...]
                                            nt4 = gmsh.model.mesh.getMaxNodeTag() + 1
                                            gmsh.model.mesh.addNodes(1, tag1D, [nt4], 𝐱̃₄)
                                            gmsh.model.mesh.addNodes(2, tag, [nt4], 𝐱̃₄)
                                            gmsh.model.mesh.addElements(1, tag1D, [1], [[gmsh.model.mesh.getMaxElementTag() + 1, gmsh.model.mesh.getMaxElementTag() + 2]], [[nt1, nt4, nt4, nt2]])
                                            gmsh.model.mesh.addElements(2, tag, [2], [[(gmsh.model.mesh.getMaxElementTag() + 1:gmsh.model.mesh.getMaxElementTag() + 2)...]], [[nt1,nt4,nt3,nt4,nt2,nt3]])
                                        end
                                        break
                                    end
                                end
                            end
                        end

                        if !isPairOnBoundary
                            index = setdiff(1:length(elementTags),indexin(removedElementTags, elementTags))
                            pairTriangle = getPairTriangle(elementTags[index], nodeTags[:,index], elementTag, nodeTag[collect(longestTag)])
                            if pairTriangle ≠ nothing
                                elementTag_, nodeTag_, longestTag_, (𝐱̄₁, 𝐱̄₂, 𝐱̄₃), L_ = pairTriangle

                                if abs(L-L_)/L < 1e-2
                                    refinedElements, addNodes, index = refine_mode_1((𝐱₁, 𝐱₂, 𝐱₃), nodeTag, longestTag, (𝐱̄₁, 𝐱̄₂, 𝐱̄₃), nodeTag_, longestTag_)
                                else
                                    refinedElements, addNodes, index = refine_mode_2((𝐱₁, 𝐱₂, 𝐱₃), nodeTag, longestTag, (𝐱̄₁, 𝐱̄₂, 𝐱̄₃), nodeTag_, longestTag_)
                                end
                                val_ = calculate_integration(refinedElements[1:2], localCoord, weights, f)
                                if abs(val-val_)/abs(val) > tol
                                    gmsh.model.mesh.removeElements(2, tag,[elementTag, elementTag_])
                                    push!(removedElementTags, (elementTag, elementTag_ )...)
                                    if length(addNodes) == 3
                                        nt1,nt2,nt3 = nodeTag[index[1:3]]
                                        nt4 = nodeTag_[index[4]]
                                        nt5 = gmsh.model.mesh.getMaxNodeTag() + 1
                                        gmsh.model.mesh.addNodes(2, tag, [nt5], addNodes)
                                        gmsh.model.mesh.addElements(2, tag, [2], [[(gmsh.model.mesh.getMaxElementTag() + 1:gmsh.model.mesh.getMaxElementTag() + 4)...]], [[nt1,nt2,nt5,nt1,nt5,nt3,nt4,nt3,nt5,nt4,nt5,nt2]])
                                    else
                                        nt1,nt2,nt3 = nodeTag[index[1:3]]
                                        nt4 = nodeTag_[index[4]]
                                        nt5 = gmsh.model.mesh.getMaxNodeTag() + 1
                                        nt6 = gmsh.model.mesh.getMaxNodeTag() + 2
                                        indexin_nt4 = indexin([nt4], nodeTag_[collect(longestTag_)])[1]
                                        gmsh.model.mesh.addNodes(2, tag, [nt5, nt6], addNodes)
                                        if indexin_nt4 == 1
                                            gmsh.model.mesh.addElements(2, tag, [2], [[(gmsh.model.mesh.getMaxElementTag() + 1:gmsh.model.mesh.getMaxElementTag() + 5)...]], [[nt1,nt2,nt5,nt1,nt5,nt3,nt2,nt6,nt5,nt3,nt5,nt6,nt2,nt4,nt6]])
                                        elseif indexin_nt4 == 2
                                            gmsh.model.mesh.addElements(2, tag, [2], [[(gmsh.model.mesh.getMaxElementTag() + 1:gmsh.model.mesh.getMaxElementTag() + 5)...]], [[nt1,nt2,nt5,nt1,nt5,nt3,nt2,nt6,nt5,nt3,nt5,nt6,nt4,nt3,nt6]])
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    filename_components = split(filename, ".")
    gmsh.write(filename_components[1]*"_refined."*filename_components[2])
    gmsh.finalize()
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

function getPairSection()
    
end

function getPairTriangle(elementTags::Vector{UInt}, nodeTags::Matrix{UInt}, elementTag_::UInt, nodeTag_::Vector{UInt})
    for (elementTag, nodeTag) in zip(elementTags, eachcol(nodeTags))
        if nodeTag_ ⊆ nodeTag && elementTag ≠ elementTag_
            𝐱₁, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[1])
            𝐱₂, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[2])
            𝐱₃, ~, ~, ~  = gmsh.model.mesh.getNode(nodeTag[3])
            longestTag, L = getLongestEdge(𝐱₁, 𝐱₂, 𝐱₃)
            return elementTag, nodeTag, longestTag, (𝐱₁, 𝐱₂, 𝐱₃), L
        end
    end

end

function calculate_integration(elements::Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}, localCoord::Matrix{Float64}, weights::Vector{Float64}, f::Function)
    val = 0.0
    for (𝐱₁, 𝐱₂, 𝐱₃) in elements
        𝐴 = 0.5*(𝐱₁[1]*𝐱₂[2]+𝐱₂[1]*𝐱₃[2]+𝐱₃[1]*𝐱₁[2]-𝐱₁[2]*𝐱₂[1]-𝐱₂[2]*𝐱₃[1]-𝐱₃[2]*𝐱₁[1])
        for ((ξ, η, ~), w) in zip(eachcol(localCoord), weights)
            N₁ = 1.0 - ξ - η
            N₂ = ξ
            N₃ = η
            x = N₁*𝐱₁[1] + N₂*𝐱₂[1] + N₃*𝐱₃[1]
            y = N₁*𝐱₁[2] + N₂*𝐱₂[2] + N₃*𝐱₃[2]
            val += f(x,y) * w * 𝐴
        end
    end
    return val
end

function refine_mode_1(element::NTuple{3, Vector{Float64}}, nodeTag::SubArray{UInt}, longestTag::Tuple{Int,Int}, element_::NTuple{3, Vector{Float64}}, nodeTag_::SubArray{UInt}, longestTag_::Tuple{Int,Int})
    if nodeTag[collect(longestTag)] ≠ nodeTag_[[longestTag_[2],longestTag_[1]]]
        longestNodeTag1 = nodeTag[collect(longestTag)]
        longestNodeTag2 = nodeTag_[collect(longestTag_)]
        error("Longest edge nodes do not match: $longestNodeTag1 vs $longestNodeTag2")
    end
    𝐱₁ = element[setdiff([1,2,3],longestTag)...]
    𝐱₂ = element[longestTag[1]]
    𝐱₃ = element[longestTag[2]]
    𝐱₄ = element_[setdiff([1,2,3],longestTag_)...]
    𝐱₅ = [0.5*(𝐱₂[1] + 𝐱₃[1]), 0.5*(𝐱₂[2] + 𝐱₃[2]), 0.5*(𝐱₂[3] + 𝐱₃[3])]
    return [(𝐱₁,𝐱₂,𝐱₅), (𝐱₁,𝐱₅,𝐱₃), (𝐱₄,𝐱₃,𝐱₅), (𝐱₄,𝐱₅,𝐱₂)], 𝐱₅, [setdiff([1,2,3],longestTag)...,longestTag[1],longestTag[2],setdiff([1,2,3],longestTag_)...]
end

function refine_mode_2(element::NTuple{3, Vector{Float64}}, nodeTag::SubArray{UInt}, longestTag::Tuple{Int,Int}, element_::NTuple{3, Vector{Float64}}, nodeTag_::SubArray{UInt}, longestTag_::Tuple{Int,Int})
    nt4 = setdiff(nodeTag_, nodeTag[collect(longestTag)])[1]
    indexin_nt4 = indexin([nt4], nodeTag_[collect(longestTag_)])[1]
    𝐱₁ = element[setdiff([1,2,3],longestTag)...]
    𝐱₂ = element[longestTag[1]]
    𝐱₃ = element[longestTag[2]]
    𝐱₄ = element_[indexin([nt4],collect(nodeTag_))[1]]
    𝐱₅ = [0.5*(𝐱₂[1] + 𝐱₃[1]), 0.5*(𝐱₂[2] + 𝐱₃[2]), 0.5*(𝐱₂[3] + 𝐱₃[3])]
    if indexin_nt4 == 1
        𝐱₆ = [0.5*(𝐱₃[1] + 𝐱₄[1]), 0.5*(𝐱₃[2] + 𝐱₄[2]), 0.5*(𝐱₃[3] + 𝐱₄[3])]
        return [(𝐱₁,𝐱₂,𝐱₅), (𝐱₁,𝐱₅,𝐱₃), (𝐱₂,𝐱₆,𝐱₅), (𝐱₃,𝐱₅,𝐱₆), (𝐱₂,𝐱₄,𝐱₆)], [𝐱₅...,𝐱₆...], [setdiff([1,2,3],longestTag)...,longestTag[1],longestTag[2],indexin([nt4],collect(nodeTag_))[1]]
    elseif indexin_nt4 == 2
        𝐱₆ = [0.5*(𝐱₂[1] + 𝐱₄[1]), 0.5*(𝐱₂[2] + 𝐱₄[2]), 0.5*(𝐱₂[3] + 𝐱₄[3])]
        return [(𝐱₁,𝐱₂,𝐱₅), (𝐱₁,𝐱₅,𝐱₃), (𝐱₂,𝐱₆,𝐱₅), (𝐱₃,𝐱₅,𝐱₆), (𝐱₄,𝐱₃,𝐱₆)], [𝐱₅...,𝐱₆...], [setdiff([1,2,3],longestTag)...,longestTag[1],longestTag[2],indexin([nt4],collect(nodeTag_))[1]]
    end
end

end
