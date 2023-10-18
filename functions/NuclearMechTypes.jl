module NuclearMechTypes

using WriteVTK
using Meshes
using SparseArrays
using Dates

export envelopeType, adherensType, inputParametersType, scaledParametersType, chromatinType, pipetteType, micromanipulationType, exportSettingsType, simulationSettingsType, replicationCompartmentType, afmType
Base.@kwdef mutable struct envelopeForcesType
    volume::Vector{Vec{3,Float64}} = []
    area::Vector{Vec{3,Float64}} = []
    bending::Vector{Vec{3,Float64}} = []
    elastic::Vector{Vec{3,Float64}} = []
    envelopeRepulsion::Vector{Vec{3,Float64}} = []
    chromationRepulsion::Vector{Vec{3,Float64}} = []
    replRepulsion::Vector{Vec{3,Float64}} = []
    ladEnveForces::Vector{Vec{3,Float64}} = []
    planeRepulsion::Vector{Vec{3,Float64}} = []
    pipetteRepulsion::Vector{Vec{3,Float64}} = []
    afmRepulsion::Vector{Vec{3,Float64}} = []
    aspiration::Vector{Vec{3,Float64}} = []
    micromanipulation::Vector{Vec{3,Float64}} = []
    total::Vector{Vec{3,Float64}} = []
end
Base.@kwdef mutable struct envelopeType
    vert::Vector{Vec{3,Float64}} = []
    neighbors::Vector{Vector{Int64}} = []
    tri::Vector{Vector{Int64}} = []
    edges::Vector{Vector{Int64}} = []
    edges2::Vector{Vector{Int64}} = []
    edgeVectors::Vector{Vec{3,Float64}} = []
    edgeVectorNorms::Vector{Float64} = []
    edgeUnitVectors::Vector{Vec{3,Float64}} = []
    mirrorEdges::Vector{Int64} = []
    firstEdges::Vector{Int64} = []
    vertexEdges::Vector{Vector{Int64}} = []
    vertexTri::Vector{Vector{Int64}} = []
    edgesTri::Vector{Vector{Int64}} = []
    edges3Vertex::Vector{Vector{Int64}} = []
    voronoiAreas::Vector{Float64} = []
    triangleAreas::Vector{Float64} = []
    curvatures::Vector{Float64} = []
    vertexNormalUnitVectors::Vector{Vec{3,Float64}} = []
    edgeNormalUnitVectors::Vector{Vec{3,Float64}} = []
    triangleNormalUnitVectors::Vector{Vec{3,Float64}} = []
    areaUnitVectors::Vector{Vector{Vec{3,Float64}}} = []
    normalVolume::Float64 = 0
    normalArea::Float64 = 0
    normalAngle::Float64 = 0
    normalLengths::Vector{Float64} = []
    normalTriangleAreas::Vector{Float64} = []
    edgeThirdVertices::Vector{Vector{Int64}} = []
    forces = envelopeForcesType()
    triEdge1::Vector{Int64} = []
    triEdge2::Vector{Int64} = []
    lads::Vector{Vector{Int64}} = []
    envelopeMultipliers::Vector{Float64} = []
end

Base.@kwdef mutable struct inputParametersType 
    freeNucleusRadius::Float64 = 0;
    laminaStiffness::Float64 = 0;
    laminaFriction::Float64 = 0;
    laminaBendingStiffness::Float64 = 0;
    osmoticPressure::Float64 = 0;
    areaCompressionModulus::Float64 = 0;
    poissonsRatio::Float64 = 0;
    bulkModulus::Float64 = 0;
    chromatinBendingModulus::Float64 = 0;
    chromatinStiffness::Float64 = 0;
    ladStrength::Float64 = 0;
    crosslinkStrength::Float64 = 0;
    chromatinNormalAngle::Float64 = 0
    viscosity::Float64 = 0;
    repulsionConstant::Float64 = 0;
    repulsionDistance::Float64 = 0;
    chroVertexDistance::Float64 = 0;
    chromatinNumber::UInt64 = 0;
    chromatinLength::UInt64 = 0;
    scalingLength::Float64 = 0;
    scalingTime::Float64 = 0;
    nSubdivisions::Float64 = 0;
    dt::Float64 = 0;
    boltzmannConst::Float64 = 1.380649e-23
    temperature::Float64 = 0;
    crosslinkingBindingProbability::Float64 = 0
    crosslinkingUnbindingProbability::Float64 = 0
    pullingForce::Float64 = 0
    iLUCutoff::Float64 = 0
    exportStep::Int64 = 0
    pipetteRadius::Float64 = 0
    aspirationPressure::Float64 = 0
    maxMovement::Float64 = 0
    minMovement::Float64 = 0
    planeForce::Float64 = 0
    replSizeMultiplier::Float64 = 0
    maxCrosslinkDistance::Float64 = 0
    chromatinFriction::Float64 = 0
    ladFriction::Float64 = 0
    replFriction::Float64 = 0
    crosslinkFriction::Float64 = 0
    maxChromosomeOverlap::Float64 = 0
    maxLadDistanceMultiplier::Float64 = 0
    minLadVertexDistance::Int64 = 0
    minLadNumber::Int64 = 0
    maxLadNumber::Int64 = 0
    staticPullingForceMultiplier::Float64 = 0
    ladLength::Float64 = 0
    crosslinkLength::Float64 = 0
    cytoskeletonPlaneRadius::Float64 = 0
    planeRepulsionConstant::Float64 = 0
    replStiffness::Float64 = 0
    replBendingStiffness::Float64 = 0
    outsideRepulsionMultiplier::Float64 = 0
    replPressure::Float64 = 0
    laminaVariabilityMultiplier::Float64 = 0
    replBulkModulus::Float64 = 0
    cantileverSpeed::Float64 = 0
    cantileverFriction::Float64 = 0
    cantileverSpring::Float64 = 0
    cantileverMaxForce::Float64 = 0
    nucleusHeight::Float64 = 0

end

Base.@kwdef mutable struct scaledParametersType
    bulkModulus::Float64 = 0;
    laminaStiffness::Float64 = 0;
    osmoticPressure::Float64 = 0;
    areaCompressionStiffness::Float64 = 0;
    laminaBendingStiffness::Float64 = 0;
    bendingStiffness::Float64 = 0;
    repulsionConstant::Float64 = 0;
    repulsionDistance::Float64 = 0;
    laminaFriction::Float64 = 0;
    chromatinBendingModulus::Float64 = 0;
    chromatinStiffness::Float64 = 0;
    ladStrength::Float64 = 0;
    crosslinkStrength::Float64 = 0;
    chromatinNormalAngle::Float64 = 0
    freeNucleusRadius::Float64 = 0;
    chroVertexDistance::Float64 = 0;
    chromatinNumber::UInt64 = 0;
    chromatinLength::UInt64 = 0;
    scalingLength::Float64 = 0
    scalingTime::Float64 = 0
    viscosity::Float64 = 0;
    nSubdivisions::Float64 = 0;
    dt::Float64 = 0;
    boltzmannConst::Float64 = 0;
    temperature::Float64 = 0;
    crosslinkingBindingProbability::Float64 = 0
    crosslinkingUnbindingProbability::Float64 = 0
    meanLaminaLength::Float64 = 0;
    pullingForce::Float64 = 0
    iLUCutoff::Float64 = 0
    exportStep::Int64 = 0
    pipetteRadius::Float64 = 0
    aspirationPressure::Float64 = 0
    maxMovement::Float64 = 0
    minMovement::Float64 = 0
    planeForce::Float64 = 0
    replSizeMultiplier::Float64 = 0
    maxCrosslinkDistance::Float64 = 0
    chromatinFriction::Float64 = 0
    ladFriction::Float64 = 0
    crosslinkFriction::Float64 = 0
    replFriction::Float64 = 0
    maxChromosomeOverlap::Float64 = 0
    maxLadDistanceMultiplier::Float64 = 0
    minLadVertexDistance::Int64 = 0
    minLadNumber::Int64 = 0
    maxLadNumber::Int64 = 0
    staticPullingForceMultiplier::Float64 = 0
    ladLength::Float64 = 0
    crosslinkLength::Float64 = 0
    cytoskeletonPlaneRadius::Float64 = 0
    planeRepulsionConstant::Float64 = 0
    replStiffness::Float64 = 0
    replBendingStiffness::Float64 = 0
    outsideRepulsionMultiplier::Float64 = 0
    replPressure::Float64 = 0
    laminaVariabilityMultiplier::Float64 = 0
    replBulkModulus::Float64 = 0
    cantileverSpeed::Float64 = 0
    cantileverFriction::Float64 = 0
    cantileverSpring::Float64 = 0
    cantileverMaxForce::Float64 = 0
    nucleusHeight::Float64 = 0
end

Base.@kwdef mutable struct chromatinForceType
    linear::Vector{Vec{3,Float64}} = []
    bending::Vector{Vec{3,Float64}} = []
    crosslink::Vector{Vec{3,Float64}} = []
    chroRepulsion::Vector{Vec{3,Float64}} = []
    enveRepulsion::Vector{Vec{3,Float64}} = []
    replRepulsion::Vector{Vec{3,Float64}} = []
    ladChroForces::Vector{Vec{3,Float64}} = []
    strandLinear::Vector{Any} = []
    strandBending::Vector{Any} = []
    strandCrosslink::Vector{Any} = []
    strandChroRepulsion::Vector{Any} = []
    strandEnveRepulsion::Vector{Any} = []
    strandReplRepulsion::Vector{Any} = []
    strandLadChroForces::Vector{Any} = []
    total::Vector{Vec{3,Float64}} = []
end
Base.@kwdef mutable struct chromatinType
    vert::Vector{Vec{3,Float64}} = []
    strandIdx::Vector{Vector{Int64}} = []
    strandVert::Vector{Any} = []
    forces = chromatinForceType()
    vectors::Vector{Vector{Vec{3,Float64}}} = []
    vectorNorms::Vector{Vector{Float64}} = []
    lads::Vector{Vector{Int64}} = []
    crosslinks::Vector{Vector{Int64}} = []
    crosslinked::Vector{Int64} = []
    neighbors::Vector{Vector{Int64}} = []
    heterochro::Vector{Bool} = []
    strandHeterochro::Vector{Any} = []
end
Base.@kwdef mutable struct pipetteType
    vert::Vector{Vec{3,Float64}} = []
    tri::Vector{Vector{Int64}} = []
    vertexTri::Array{Vector{Int64}} = Array{Int64}[]
    neighbors::Vector{Vector{Int64}} = []
    edges::Vector{Vector{Int64}} = []
    vertexNormalUnitVectors::Vector{Vec{3,Float64}} = []
    edgeNormalUnitVectors::Vector{Vec{3,Float64}} = []
    triangleNormalUnitVectors::Vector{Vec{3,Float64}} = []
    edgeVectors::Vector{Vec{3,Float64}} = []
    firstEdges::Vector{Int64} = []
    triEdge1::Vector{Int64} = []
    triEdge2::Vector{Int64} = []
    mirrorEdges::Vector{Int64} = []
    edgesTri::Vector{Vector{Int64}} = []
    pipetteTree::Any = []
end

Base.@kwdef mutable struct micromanipulationType
    rightmostVertex::Int64 = 0
    leftmostVertex::Int64 = 0
    leftmostVertexPosition::Vec{3,Float64} = Vec(0.,0.,0.)
    pipettePosition::Vec{3,Float64} = Vec(0.,0.,0.)
    pipetteMovements::Array{Float64} = []
end

Base.@kwdef mutable struct afmType
    beadPosition::Vec{3,Float64} = Vec(0.,0.,0.)
    topPosition::Vec{3,Float64} = Vec(0.,0.,0.)
    touchingIdx::Vector{Bool} = []
    normDistance::Float64 = 0;
    directionDown::Bool = true
    forceOnBead::Float64 = 0;
end

Base.@kwdef mutable struct exportSettingsType
    exportData::Bool = true
    enveCells::Vector{MeshCell{VTKCellType, Vector{Int64}}} = []
    chroCells::Vector{MeshCell{PolyData.Lines, Vector{Int64}}} = []
    ladCells::Vector{MeshCell{PolyData.Lines, Vector{Int64}}} = []
    ladIdx::Vector{Int} = []
    ladEnveVertices::Vector{Int} = []
    ladChroVertices::Vector{Int} = []
    folderName::String = ""
    step::Int64 = 1
end

Base.@kwdef mutable struct adherensType
    adherent::Bool = false
    topPlane::Float64 = 0
    bottomPlane::Float64 = 0
    touchingTop::Vector{Bool} = []
    cellForcesOnPlane::Float64 = 0
    static::Bool = false
end

Base.@kwdef mutable struct simulationSettingsType
    frictionMatrix::Any = []
    envelopeTree::Any = []
    chromatinTree::Any = []
    simType::String = ""
    prog::Any = []
    initType::String = ""
    timeStepProgress::Rational = 0;
    timeStepMultiplier::Rational = 1;
    timeStepTiming::DateTime = now();
    noEnveSolve::Bool = false;
    noChromatin::Bool = true
    iLU::Any = []
    adh = adherensType();
    laminaRemodel::String = ""
    newVolumeSimulation::Bool = false
end
Base.@kwdef mutable struct virusforcesType
    volume::Vector{Vec{3,Float64}} = []
    area::Vector{Vec{3,Float64}} = []
    bending::Vector{Vec{3,Float64}} = []
    elastic::Vector{Vec{3,Float64}} = []
    envelopeRepulsion::Vector{Vec{3,Float64}} = []
    chromationRepulsion::Vector{Vec{3,Float64}} = []
    total::Vector{Vec{3,Float64}} = []
end

Base.@kwdef mutable struct replicationCompartmentType
    infected::Bool = false
    vert::Vector{Vec{3,Float64}} = []
    neighbors::Vector{Vector{Int64}} = []
    tri::Vector{Vector{Int64}} = []
    edges::Vector{Vector{Int64}} = []
    edges2::Vector{Vector{Int64}} = []
    edgeVectors::Vector{Vec{3,Float64}} = []
    edgeVectorNorms::Vector{Float64} = []
    edgeUnitVectors::Vector{Vec{3,Float64}} = []
    mirrorEdges::Vector{Int64} = []
    firstEdges::Vector{Int64} = []
    vertexEdges::Vector{Vector{Int64}} = []
    vertexTri::Vector{Vector{Int64}} = []
    edgesTri::Vector{Vector{Int64}} = []
    edges3Vertex::Vector{Vector{Int64}} = []
    voronoiAreas::Vector{Float64} = []
    triangleAreas::Vector{Float64} = []
    curvatures::Vector{Float64} = []
    vertexNormalUnitVectors::Vector{Vec{3,Float64}} = []
    edgeNormalUnitVectors::Vector{Vec{3,Float64}} = []
    triangleNormalUnitVectors::Vector{Vec{3,Float64}} = []
    areaUnitVectors::Vector{Vector{Vec{3,Float64}}} = []
    normalVolume::Float64 = 0
    normalArea::Float64 = 0
    normalAngle::Float64 = 0
    normalLengths::Vector{Float64} = []
    normalTriangleAreas::Vector{Float64} = []
    initVolume::Float64 = 0
    edgeThirdVertices::Vector{Vector{Int64}} = []
    forces = virusforcesType()
    triEdge1::Vector{Int64} = []
    triEdge2::Vector{Int64} = []
    frictionMatrix::Any = []
    iLU::Any = []
    baseArea::Float64 = 0;
    tree::Any = []
    
end

end