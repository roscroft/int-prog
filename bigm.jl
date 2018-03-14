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

function bigM(i)
    lowM = m[i]*d[1] + c[i] - f(d[1])
    highM = f(d[1]) - m[i]*d[1] - c[i]
    for j in 2:k+1
        new_lowM = m[i]*d[j] + c[i] - f(d[j])
        new_highM = f(d[j]) - m[i]*d[j] - c[i]
        if new_lowM > lowM
            lowM = new_lowM
        end
        if new_highM > highM
            highM = new_highM
        end
    end
    return (lowM, highM)
end

function make_model()
    model = Model()
    @variable(model, x)
    @variable(model, z)
    @variable(model, y[1:k], Bin)
    for i in 1:k
        (lowM, highM) = bigM(i)
        @constraint(model, d[i] - (d[i] - d[1])*(1-y[i]) <= x)
        @constraint(model, d[i+1] + (d[k+1]-d[i+1])*(1-y[i]) >= x)
        @constraint(model, m[i]*x + c[i] - lowM*(1-y[i]) <= z)
        @constraint(model, m[i]*x + c[i] + highM*(1-y[i]) >= z)
    end
    @constraint(model, sum(y[i] for i=1:k) == 1)

    poly = polyhedron(model, CDDLibrary(:exact))
    vertices = SimpleVRepresentation(poly)
    removevredundancy!(poly)
    vertices = SimpleVRepresentation(poly)
    return vertices
end

function conv_hull()
    pts = Array{Int64,2}[]
    for i in 1:k
        y = [0 0 0 0]
        y[i] = 1
        first = [d[i] f(d[i]) y]
        pts = [pts; first]
        second = [d[i+1] f(d[i+1]) y]
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
    for i in 1:k
        (lowM, highM) = bigM(i)
        @constraint(model, d[i] - (d[i] - d[1])*(1-y[i]) <= x)
        @constraint(model, d[i+1] + (d[k+1]-d[i+1])*(1-y[i]) >= x)
        @constraint(model, m[i]*x + c[i] - lowM*(1-y[i]) <= z)
        @constraint(model, m[i]*x + c[i] + highM*(1-y[i]) >= z)
    end
    @constraint(model, sum(y[i] for i=1:k) == 1)
    
    new_poly = conv_hull()
    ineq = SimpleHRepresentation(new_poly)
    A = convert.(Int64,ineq.A)
    b = convert.(Int64,ineq.b)
    for i=1:size(A,1)
        @constraint(model, A[i,1]*x + A[i,2]*z + sum(A[i,j+2]*y[j] for j=1:k) <= b[i])
    end

    for eq in eqs(new_poly)
        LHS = convert.(Int64,eq.a)
        RHS = convert.(Int64,eq.β)
        @constraint(model, LHS[1]*x + LHS[2]*z + sum(LHS[j+2]*y[j] for j=1:k) == RHS)
    end
    #= println(model) =#
    poly = polyhedron(model, CDDLibrary(:exact))
    vertices = SimpleVRepresentation(poly)
    removevredundancy!(poly)
    vertices = SimpleVRepresentation(poly)
    return vertices
end
