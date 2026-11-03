"""
    Handle

    Handle(path::AbstractString, [indexpath::AbstractString])

Load a shapefile into GeoInterface compatible objects. This can be plotted
with Plots.jl `plot`.

The Vector of shape object can be accessed with `shapes(handle)`.

`Handle` may have a known bounding box, which can be retrieved with `GeoInterface.bbox`.
"""
mutable struct Handle{T<:Union{<:AbstractShape,Missing}}
    header::Header
    shapes::Vector{T}
    crs::Union{Nothing, GFT.ESRIWellKnownText{GFT.CRS}}
end
function Handle(path::AbstractString)
    shx = splitext(path)[1] * ".shx"
    if isfile(shx)
        Handle(path, shx)
    else
        Handle(path, nothing)
    end
end
function Handle(path::AbstractString, index)
    open(path) do io
        read(io, Handle, index; path = path)
    end
end
function Handle(path::AbstractString, indexpath::AbstractString)
    index = open(indexpath) do io
        read(io, IndexHandle)
    end
    Handle(path, index)
end

shapes(h::Handle) = h.shapes

# GeoInterface
GI.isgeometry(::Type{<:Handle}) = true
GI.geomtrait(::Handle) = GI.GeometryCollectionTrait()
GI.ncoord(::GI.GeometryCollectionTrait, h::Handle) = GI.ncoord(first(h.shapes))
GI.ngeom(::GI.GeometryCollectionTrait, h::Handle) = length(h.shapes)
GI.getgeom(::GI.GeometryCollectionTrait, h::Handle, i) = h.shapes[i]
GI.getgeom(::GI.GeometryCollectionTrait, h::Handle) = h.shapes
GI.crs(h::Handle) = h.crs

Base.length(shp::Handle) = length(shapes(shp))

function Base.read(io::IO, ::Type{Handle}, index = nothing; path = nothing)
    header = read(io, Header)
    read_shapes(::Type{T}) where {T} = _read_handle_inner(io::IO, T, header, shapes, index; path)
    # One call per shape type, instead of looking the type up in SHAPETYPE, so the
    # compiler knows every call. This is needed to compile with `juliac --trim`.
    code = header.shapecode
    return code == 0 ? read_shapes(Missing) :
        code == 1 ? read_shapes(Point) :
        code == 3 ? read_shapes(Polyline) :
        code == 5 ? read_shapes(Polygon) :
        code == 8 ? read_shapes(MultiPoint) :
        code == 11 ? read_shapes(PointZ) :
        code == 13 ? read_shapes(PolylineZ) :
        code == 15 ? read_shapes(PolygonZ) :
        code == 18 ? read_shapes(MultiPointZ) :
        code == 21 ? read_shapes(PointM) :
        code == 23 ? read_shapes(PolylineM) :
        code == 25 ? read_shapes(PolygonM) :
        code == 28 ? read_shapes(MultiPointM) :
        code == 31 ? read_shapes(MultiPatch) :
        throw(KeyError(code))
end

function _read_handle_inner(io::IO, ::Type{T}, header, shapes, index=nothing; path) where T
    shapes = Vector{Union{T,Missing}}(undef, 0)
    num = Int32(0)
    while (!eof(io))
        seeknext(io, num, index)
        num = bswap(read(io, Int32))
        rlength = bswap(read(io, Int32))
        shapecode = read(io, Int32)
        if shapecode === Int32(0)
            push!(shapes, missing)
        else
            push!(shapes, read(io, T))
        end
    end
    crs = nothing
    if path !== nothing
        prjpath = _shape_paths(path).prj
        if isfile(prjpath)
            try
                crs = GFT.ESRIWellKnownText(GFT.CRS(), read(open(prjpath), String))
            catch
                @warn "Projection file $prjpath appears to be corrupted. `nothing` used for `crs`"
            end
        end
    end
    return Handle(header, shapes, crs)
end

function seeknext(io, num, index::IndexHandle)
    seek(io, index.indices[num + 1].offset * 2)
end
seeknext(io, num, ::Nothing) = nothing
