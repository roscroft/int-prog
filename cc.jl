using JuMP, CDDLib, Polyhedra
# Generate the Big-M formulation.
k = 4
d = [0,1,2,3,4]
m = [-1,2,-2,1]
c = [1,-2,6,-3]

function f(pt)
    if 0<=pt && pt<=1
        return -1*pt + 1
    elseif 1<pt && pt<=2
        return 2*pt - 2
    elseif 2<pt && pt<=3
        return -2*pt + 6
    elseif 3<pt && pt<=4
        return pt - 3
    end
end

function make_model()
    model = Model()
    @variable(model, x)
    @variable(model, z)
    @variable(model, y[1:k], Bin)
    @variable(model, 0<=la[1:k+1]<=1)
    @constraint(model, sum(la[i]*d[i] for i=1:k+1) == x)
    @constraint(model, sum(la[i]*f(d[i]) for i=1:k+1) ==z)
    @constraint(model, sum(la[i] for i=1:k+1) == 1)
    @constraint(model, la[1] <= y[1])
    @constraint(model, [i=2:k], la[i] <= y[i-1] + y[i])
    @constraint(model, la[k+1] <= y[k])
    @constraint(model, sum(y[i] for i=1:k) == 1)

    poly = polyhedron(model, CDDLibrary(:exact))
    vertices = SimpleVRepresentation(poly)
    removevredundancy!(poly)
    vertices = SimpleVRepresentation(poly)
    return vertices
end

function eik(i,k)
    e = zeros(Int64, 1, k)
    for j = 1:k
        if i == j
            e[j] = 1
        else
            e[j] = 0
        end
    end
    return e
end

function conv_hull()
    pts = Array{Int64,2}[]
    for i in 1:k
        y = [0 0 0 0]
        y[i] = 1
        first = [d[i] f(d[i]) y eik(i,k+1)]
        pts = [pts; first]
        second = [d[i+1] f(d[i+1]) y eik(i+1,k+1)]
        pts = [pts; second]
    end
    pts = convert(Array{Int64,2},pts)
    vertices = SimpleVRepresentation(pts)
    poly = polyhedron(vertices, CDDLibrary(:exact))
    return poly
end

function make_fixed_model()
    model = Model()
    @variable(model, x)
    @variable(model, z)
    @variable(model, y[1:k], Bin)
    @variable(model, 0<=la[1:k+1]<=1)
    @constraint(model, sum(la[i]*d[i] for i=1:k+1) == x)
    @constraint(model, sum(la[i]*f(d[i]) for i=1:k+1) ==z)
    @constraint(model, sum(la[i] for i=1:k+1) == 1)
    @constraint(model, la[1] <= y[1])
    @constraint(model, [i=2:k], la[i] <= y[i-1] + y[i])
    @constraint(model, la[k+1] <= y[k])
    @constraint(model, sum(y[i] for i=1:k) == 1)
    
    new_poly = conv_hull()
    ineq = SimpleHRepresentation(new_poly)
    ineqA = ineq.A
    ineqb = ineq.b
    lcmA = lcm(denominator.(ineqA))
    A = convert.(Int64, lcmA*ineqA)
    b = convert.(Int64, lcmA*ineqb)
    for i=1:size(A,1)
        @constraint(model, A[i,1]*x + A[i,2]*z + sum(A[i,j+2]*y[j] for j=1:k) + sum(A[i,j+6]*la[j] for j=1:k) <= b[i])
    end

    for eq in eqs(new_poly)
        eqa = eq.a
        eqB = eq.β
        lcma = lcm(denominator.(eqa))
        LHS = convert.(Int64, lcma*eqa)
        RHS = convert.(Int64, lcma*eqB)
        @constraint(model, LHS[1]*x + LHS[2]*z + sum(LHS[j+2]*y[j] for j=1:k) + sum(LHS[j+6]*la[j] for j=1:k) == RHS)
    end
    poly = polyhedron(model, CDDLibrary(:exact))
    vertices = SimpleVRepresentation(poly)
    removevredundancy!(poly)
    vertices = SimpleVRepresentation(poly)
    return vertices
end
