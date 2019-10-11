module OpSel

using Printf
using DataFrames
using JuMP
using Ipopt
using DelimitedFiles
using SparseArrays
using LinearAlgebra

# function testEqual(sp_csv, N, theta)
function testEqual(Z, N = 50)
    sp_csv = loadDataBySize(Z)
    (x_result, info_result) = steepestAscent(sp_csv, N)
    return (x_result, info_result, sp_csv)
end


thetaTable = spzeros(300000, 30000)
thetaTable[ 50,   200] = 3.3429999999999994e-02/2
thetaTable[ 50,  1050] = 6.2748061171875091e-02/2
thetaTable[ 50,  2045] = 7.1083439999999956e-02/2
thetaTable[ 50,  5050] = 1.0806060937499347e-01/2
thetaTable[ 50,  5255] = 2.4439005000000000e-02/2
thetaTable[ 50, 10100] = 7.0094325234373603e-02/2
thetaTable[ 50, 15100] = 6.7641368203122818e-02/2
thetaTable[ 50, 15222] = 3.8808090000000316e-02/2
thetaTable[100,   200] = 2.5817499999999993e-02/2
thetaTable[100,  1050] = 5.3925029121093773e-02/2
thetaTable[100,  2045] = 6.2823779999999996e-02/2
thetaTable[100,  5050] = 9.9463449550774755e-02/2
thetaTable[100,  5255] = 1.5510273749999984e-02/2
thetaTable[100, 10100] = 6.0998935781248584e-02/2
thetaTable[100, 15100] = 5.8678192968747774e-02/2
thetaTable[100, 15222] = 3.0043923750000281e-02/2

function steepestAscent(sp_csv, N=50)
    Z = size(sp_csv,1)
    try
        theta = thetaTable[N, Z]
    catch error1
        error(@sprintf("theta for N = %d, Z = %d is not available in opsel.thetaTable\n", N, Z))
    end
        
    info_result = Dict()
    info_result["build_time"] =  @elapsed begin
        Z = size(sp_csv,1)
        Ainv = makeAinv(sp_csv)
        g = sp_csv[:,4]
        l = zeros(Z)
        u = min.(sp_csv[:,5], 1/N)
        # A = inv(full(Ainv))
        A = inv(Matrix(Ainv))
        A = (A+A')/2
        if theta < 0
            println("*** Negative theta is given. Setting up new theta")
            theta = newTheta(Z, A, u, N)
        end
        
        println("Z = $Z, N = $N, theta = $theta")


        # JuMP
        model1 = Model(with_optimizer(Ipopt.Optimizer, acceptable_tol=1e-5))

        @variable(model1, y[1:Z])
        @constraint(model1, Ainv*y .<=u)
        @constraint(model1, l.<=Ainv*y)
        @constraint(model1, sum(Ainv*y) == 1.0)
        (I1,J1,V1) = findnz(Ainv)
        @constraint(model1, sum(V1[k]*y[I1[k]]*y[J1[k]] for k=1:length(V1)) <=2*theta)

        @objective(model1, Max, (Ainv*g)'*y)

    end # build _time 
        
        
    info_result["solver_time"] = @elapsed optimize!(model1)
    
    y_value = value.(y)
    x_value = Ainv*y_value
    x_value = slimData(x_value,1.0e-5)
    
    @printf("=== SOCP Result Summary ===\n")
    @printf("JuMP status = %s\n", termination_status(model1))
    info_result["JuMPstatus"] = @sprintf("%s", termination_status(model1))
    info_result["socp_time"] = info_result["build_time"] + info_result["solver_time"]
    @printf("sum(x_socp) = %lf, gx = %lf, xAx = %lf, 2theta = %lf\n",
            sum(x_value), dot(g, x_value),
            dot(y_value, Ainv*y_value), 2*theta)


    info_result["steep_time"] = @elapsed begin
        # compute Lambda0 by Meuwissen's lagrange multiplier
        lambda0 = computeLambda0(Z, Ainv, g, theta);
        lambda0 = 2*lambda0; # mulitply with 2 is import for Z = 1050
        (x_sa, iter_sa, beforecobj_socp, aftercobj_socp) = localSwap(x_value, sp_csv, A, Ainv, N, theta, lambda0)
        info_result["steep_iter"] = iter_sa
    end
    
    @printf("sum(x_sa) = %lf, gx = %lf, xAx = %lf, 2theta = %lf\n",
            sum(x_sa), dot(g,x_sa),
            dot(x_sa, Ainv\x_sa), 2*theta)

    gx_value = dot(g, x_value)
    gx_sa    = dot(g, x_sa)

    @printf("SOCP = %f, STEEP = %f, gap = %f%%\n",
        gx_value, gx_sa, abs(gx_value - gx_sa)/max(gx_value, gx_sa)*100)
    
    info_result["total_time"] = info_result["socp_time"] + info_result["steep_time"]

    @printf("time(s): build = %.3f, solver = %.3f, steep = %.3f, total = %.3f\n",
            info_result["build_time"],
            info_result["solver_time"],
            info_result["steep_time"],
            info_result["total_time"])
    return (x_sa, info_result)
end

function testUnequal(Z)
    sp_csv = loadDataBySize(Z)
    if Z == 200
        println("Z = 200 is not solvable for N_s = 14 and N = 2800")
        return (0, 0, 0)
    end    
    (x_result, info_result) = compactSOCP(sp_csv, 14, 2800)
    return (x_result, info_result, sp_csv)
end


function compactSOCP(sp_csv, N_s = 14, N = 2800)

    # N_s is the effective number
    # In the data, we fixed N_s = 14

    info_result = Dict()
    info_result["build_time"] =  @elapsed begin
        Z = size(sp_csv,1)
        theta = 0.5/N_s
        println("Z = $Z, N_s = $N_s, N = $N, theta = $theta")

        p = sp_csv[:,2]
        q = sp_csv[:,3]
        g = sp_csv[:,4]
        l = zeros(Z)
        u = sp_csv[:,5]/N
        inbreeding = sp_csv[:,6]
        Ainv = makeAinv(sp_csv)

        # JuMP
        model1 = Model(with_optimizer(Ipopt.Optimizer, acceptable_tol=1e-5))

        @variable(model1, y[1:Z])
        @constraint(model1, Ainv*y .<=u)
        @constraint(model1, l.<=Ainv*y)
        @constraint(model1, sum(Ainv*y) == 1.0)
        (I1,J1,V1) = findnz(Ainv)
        @constraint(model1, sum(V1[k]*y[I1[k]]*y[J1[k]] for k=1:length(V1)) <=2*theta)
        @objective(model1, Max, dot(Ainv*g, y))
        # @objective(model1, Max, dot(Ainv*g, y))
    end # build _time 
        
        
    info_result["solver_time"] = @elapsed optimize!(model1)
    info_result["total_time"] = info_result["build_time"] + info_result["solver_time"]
    
    y_value = value.(y)
    x_value = slimData(Ainv*y_value,0)
    
    @printf("=== Result Summary ===\n")
    @printf("JuMP status = %s\n", termination_status(model1))
    info_result["JuMPstatus"] = @sprintf("%s", termination_status(model1))
        
    
    println("Z = $Z, N_s = $N_s, N = $N")
    @printf("gx = %lf, xAx = %lf, 2theta = %lf\n", dot(g, x_value),
            dot(y_value, Ainv*y_value), 2*theta)
    @printf("time(s): build = %.3f, solver = %.3f, total = %.3f\n",
            info_result["build_time"],
            info_result["solver_time"],
            info_result["total_time"])
    return (x_value, info_result)
end


function loadDataBySize(Z)
    filename = @sprintf("%s/../data/sorted%06d.csv",@__DIR__, Z)
    sp_csv = loadFile(filename)
    return sp_csv
end

function loadFile(filename)
    sp_csv = readdlm(filename,',')
    sp_csv = convert(Array{Float64,2}, sp_csv[2:end,:])
    # ignore the column names in the first line
    return sp_csv
end


function slimData(A, epsilon)
    if length(size(A)) == 1
        return A.*(abs.(A).>=epsilon)
    end
    (m,n) = size(A)
    (I1,J1,V1) = findnz(A)
    V2 = V1.*(abs(V1).>=epsilon)
    A2 = sparse(I1,J1,V2,m,n)
    return A2
end

# function [Ainv] = makeAinv2( sp_csv )
#   compute Ainverse by Henderson's formula

function makeAinv(sp_csv)
    Z = size(sp_csv,1)
    inbreeding = sp_csv[:,6]

    eleNo = 0
    for i = 1:Z
        p = convert(Int, sp_csv[i,2])
        q = convert(Int, sp_csv[i,3])
        if p != 0 && q != 0
            eleNo = eleNo + 9
            continue
        end
        if (p != 0 && q == 0) || (p == 0 && q != 0 )
            eleNo = eleNo + 4
            continue
        end
        if p == 0 && q == 0
            eleNo = eleNo + 1
        end
    end

    I = zeros(Int, eleNo)
    J = zeros(Int, eleNo)
    V = zeros(eleNo)

    eleNo = 0
    for i = 1:Z
        p = convert(Int, sp_csv[i,2])
        q = convert(Int, sp_csv[i,3])
        if p != 0 && q != 0
            b = 4.0/( (1-inbreeding[p]) + (1-inbreeding[q]) )
            I[eleNo + 1] = i; J[eleNo + 1] = i; V[eleNo + 1] =  1.0*b;
            I[eleNo + 2] = p; J[eleNo + 2] = i; V[eleNo + 2] = -0.5*b;
            I[eleNo + 3] = q; J[eleNo + 3] = i; V[eleNo + 3] = -0.5*b;
            I[eleNo + 4] = i; J[eleNo + 4] = p; V[eleNo + 4] = -0.5*b;
            I[eleNo + 5] = i; J[eleNo + 5] = q; V[eleNo + 5] = -0.5*b;
            I[eleNo + 6] = p; J[eleNo + 6] = p; V[eleNo + 6] = +0.25*b;
            I[eleNo + 7] = p; J[eleNo + 7] = q; V[eleNo + 7] = +0.25*b;
            I[eleNo + 8] = q; J[eleNo + 8] = p; V[eleNo + 8] = +0.25*b;
            I[eleNo + 9] = q; J[eleNo + 9] = q; V[eleNo + 9] = +0.25*b;
            eleNo = eleNo + 9
            continue
        end
        if (p != 0 && q == 0) || (p == 0 && q != 0 )
            p = max(p,q)
            b = 4.0/( 1.0*(1-inbreeding[p]) + 2.0*(1-0) )
            I[eleNo + 1] = i; J[eleNo + 1] = i; V[eleNo + 1] =  1.0*b;
            I[eleNo + 2] = p; J[eleNo + 2] = i; V[eleNo + 2] = -0.5*b;
            I[eleNo + 3] = i; J[eleNo + 3] = p; V[eleNo + 3] = -0.5*b;
            I[eleNo + 4] = p; J[eleNo + 4] = p; V[eleNo + 4] = +0.25*b;
            eleNo = eleNo + 4
            continue
        end
        if p == 0 && q == 0
            b = 4.0/( 2.0*(1-0) + 2.0*(1-0) )
            I[eleNo + 1] = i; J[eleNo + 1] = i; V[eleNo + 1] =  1.0*b;
            eleNo = eleNo + 1
        end
    end

    Ainv = sparse(I,J,V,Z,Z)
    # Ainv = Ainv + Ainv' - diag(diag(Ainv));
    return Ainv
end # end of makeAinv

function makeB(sp_csv, f)

    Z = size(sp_csv,1)
    p = convert(Array{Int,1}, sp_csv[:,2])
    q = convert(Array{Int,1}, sp_csv[:,3])

    indI = zeros(Int,3*Z)
    indJ = zeros(Int,3*Z)
    indV = zeros(3*Z)
    ind0 = 1

    for i=1:Z
        if  p[i] == 0 && q[i] == 0
            b = 4.0/(2*(1-0) + 2*(1-0))
            indI[ind0] = i
            indJ[ind0] = i
            indV[ind0] = 1*sqrt(b)
            ind0 = ind0+1
        elseif q[i] == 0
            b = 4.0/(1.0*(1-f[p[i]]) + 2.0*(1-0))
            indI[ind0] = i
            indJ[ind0] = i
            indV[ind0] = 1*sqrt(b)
            ind0 = ind0+1
            indI[ind0] = p[i]
            indJ[ind0] = i
            indV[ind0] = -0.5*sqrt(b)
            ind0 = ind0+1
        else
            b = 4.0/((1-f[p[i]])+(1-f[q[i]]));
            indI[ind0] = i
            indJ[ind0] = i
            indV[ind0] = 1.0*sqrt(b)
            ind0 = ind0+1
            indI[ind0] = p[i]
            indJ[ind0] = i
            indV[ind0] = -0.5*sqrt(b)
            ind0 = ind0+1
            indI[ind0] = q[i]
            indJ[ind0] = i
            indV[ind0] = -0.5*sqrt(b)
            ind0 = ind0+1
        end
    end

    B = sparse(indJ[1:ind0-1],indI[1:ind0-1],indV[1:ind0-1],Z,Z)
    return B
end # end of function makeB


function localSwap(x_value, sp_csv, A, Ainv, N, theta_relax, lambda0)

    Z = size(sp_csv,1)
    l = zeros(Z,1)
    g = sp_csv[:,4]
    u = min.(sp_csv[:,5],1/N)
    u_orig = u

    sort_index = sortperm(x_value, rev = true)
    sort_x = x_value[sort_index]
    xf = zeros(Z,1)
    F = convert(Array{Int,1},sort_index[1:N])
    xf[F] .= 1/N
    V = setdiff(findall(x -> x>=1/N, u_orig), F);

    take_max = 0; # usually zero
    # take_max = -inf; # usually zero
    c_obj = dot(g,xf) - lambda0*max(take_max, dot(xf,A*xf)/2 - theta_relax)
    before_cobj = c_obj
    @printf("======= Before swap c_obj = %f\n", c_obj);

    iter2 = 0
    xAx_prev = dot(xf, A*xf)/2
    diagA = diag(A)
    while true
        changed = 0
        iter2 = iter2+1
        A_FF = sum(A[F,F])
        gF = sum(g[F])
        best_f = 0
        best_v = 0
        for index_f = 1:length(F)
            f = F[index_f]
            xAx_vec =  (A_FF + diagA[f] - 2*sum(A[F,f])) .+
              (diagA[V] -2*A[V,f]+2*sum(A[V,F],dims=2))
            c_obj3 = 1/N*((gF - g[f]) .+ g[V]) -
                lambda0 * max.(take_max, 1/N^2*xAx_vec/2 .- theta_relax)
            (c_obj2, max_v) = findmax(c_obj3)
            if c_obj2 > c_obj
                changed = 1
                best_f = f
                best_v = V[max_v]
                c_obj = c_obj2
            end
        end
        if changed == 0
            @printf("No swap c_obj = %f, gx = %f, xAx = %f\n",
                    c_obj, dot(g,xf), dot(xf,A*xf));
            break
        end
        xf[best_f] = 0
        xf[best_v] = 1/N
        F = union(setdiff(F,best_f),best_v)
        V = union(setdiff(V,best_v),best_f)
        
        @printf("iter2 = %d, c_obj = %f, gx = %f, xAx = %f, best_f = %d, best_v = %d\n",
                iter2, c_obj, dot(g,xf), dot(xf, Ainv\xf), best_f, best_v)
    end
    after_cobj = c_obj
    @printf("Finished ... iter2 = %d, before_cobj = %f, after_cobj = %f\n",
            iter2, before_cobj, after_cobj)

    return  (xf, iter2, before_cobj, after_cobj)

end # end of localSswap


function computeLambda0(Z, Ainv, g, theta)
#  compute lambda0 based on Meuwissen's lagrange multiplier

    Q = ones(Z,1)
    r = ones(1,1)
    theta_min = dot(r,(Q'*Ainv*Q)\r)/2
    if theta_min > theta
        @printf("theta (%f) is too small. It must be >= %f\n", theta, theta_min)
        error("theta error")
    end
    Ainv_g = Ainv*g
    B = inv(Q'*Ainv*Q)
    l_bunshi = dot(g,Ainv_g) - dot(Q'*Ainv_g, B*(Q'*Ainv_g))
    l_bunbo = 8*theta - 4*dot(r, B*r)
    lambda0 = 2*sqrt(l_bunshi/l_bunbo)*2 # mulitply with 2 is import for Z = 1050

    return lambda0
end # end of computeLambda0

end # end of module 
