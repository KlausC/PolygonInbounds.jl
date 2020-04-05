
"""
    mean_overlap(a, y0, x0, x1, y1)

Mean overlap of an interval of length `a` and interval `[x0,x1]`,
when moving inside `[y0,y1]`.

Assumptions `a <= y1 - y0` and ` y0 <= x0 <= x1 <= y1`.

Return `result` with `min(x1 - x0, a) <= result <= a`.
"""
function mean_overlap(a::T, y0::T, x0::T, x1::T, y1::T) where T <: Real
    x0 = max(x0, y0)
    x1 = min(x1, y1)
    y0 <= x0 <= x1 <= y1 || throw(ArgumentError("violating ordering of boundaries"))
    a = min(a, y1 - y0)
    y10 = y1 - y0 - a
    if y10 <= (eps())*a
        return min(x1 - x0, a)
    end
    m11 = max(x0 - y0 - a, 0)
    m12 = max(min(x0+a, x1) - max(y0+a, x0), 0)
    d12 = (min(x0+a, x1) + max(y0+a, x0)) / 2 - x0
    m13 = max(min(x0+a, y1) - max(y0+a, x1), 0)
    m22 = max(x1 - x0 - a, 0)
    m23 = max(min(x1+a, y1) - max(x0+a, x1), 0)
    d23 = -(min(x1+a, y1) + max(x0+a, x1)) /2 + x1 + a
    m33 = max(y1 - x1 - a, 0)

    #println(m11 + m12 + m13 + m22 + m23 + m33 - y1 + y0 + a)
    #display([m11 m12 m13 m22 m23 m33;0 d12 x1 - x0 a d23 0])
    (m12*d12 + m13*(x1-x0) + m22*a + m23*d23) / y10
end

function sumeff(points, vmin, vmax, poly, pmin, pmax, xy, tol)
    pmi = pmin[xy] - tol
    pma = pmax[xy] + tol
    vmi = vmin[xy] - tol
    vma = vmax[xy] + tol
    sum(mean_overlap.(all_boxlengths(poly, xy) .+ 2tol, pmi, vmi, vma, pma))
end

function countmid(points, pmin, pmax, xy, tol)
    npoints = length(points)
    v = (vertex(points, k, xy) for k = 1:npoints)
    pmi = pmin[xy] - tol
    pma = pmax[xy] + tol
    low, mid, hi = sum([x < pmi, pmi <= x <= pma, x > pma] for x in v)
    mid
end

function estimate_effort(points, vmin, vmax, poly, pmin, pmax, tol) 
    sumx = sumeff(points, vmin, vmax, poly, pmin, pmax, 1, tol)
    sumy = sumeff(points, vmin, vmax, poly, pmin, pmax, 2, tol)
    sumx0 = sumeff(points, vmin, vmax, poly, pmin, pmax, 1, 0.0)
    sumy0 = sumeff(points, vmin, vmax, poly, pmin, pmax, 2, 0.0)
    
    xmid = countmid(points, pmin, pmax, 1, tol)
    ymid = countmid(points, pmin, pmax, 2, tol)
    xmid0 = countmid(points, pmin, pmax, 1, 0.0)
    ymid0 = countmid(points, pmin, pmax, 2, 0.0)
    
    xmid * sumx, ymid * sumy, xmid0 * sumx0, ymid0 * sumy0
end

