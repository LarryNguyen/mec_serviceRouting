using JuMP, GLPKMathProgInterface
using CSV
using JSON

##########################################
### Function Utility
##########################################
function getIndexFromBsId(bsId)
    for i in 1:sizeN
        if N[i]==bsId
            return i
        end
    end
    return 0
end
function getBsIndexFromUser(u)
    for i in 1:sizeN
        if E_u[i,u]==1
            return i
        end
    end
    return 0
end
##########################################
### Dataset
##########################################
N = 7
U = 10
# Base station
#bs_file = CSV.read("dataset/BS_US_PN_Location.csv")
#println(bs_file)
#N = bs_file[1]
sizeN = length(N)

## User id U
#user_file = CSV.read("dataset/Twitter_User_US_PN_Location.csv")
#println(user_file)
#U = user_file[1]
sizeU = length(U)

## List of Edge connected with user E_u[nodeId, userId]
user_bs_file=readdlm("dataset/mapping_user_bs.txt")
#println(length(user_bs_file))
E_u=zeros(sizeN, sizeU)
r_u=zeros(sizeU)

#e.g. user 1 -> edge 1; user 2,3 --> edge 2
#E_u[1,1] = 1
#E_u[2,2] = 1
#E_u[2,3] = 1
for i in 1:sizeU
    tmp_bs = user_bs_file[i,2]
    r_u[i] = user_bs_file[i,4]
    E_u[getIndexFromBsId(tmp_bs),i]=1
end
println(E_u)

## Load Service data
service_file=readdlm("dataset/service.txt")
rowS = service_file[1,:]
sizeS = length(service_file)/length(rowS)
sizeS = convert(Int, sizeS)

## Service resource usage r_cpu[service] , r_mem[service] (cpu, mem)
r_cpu=zeros(sizeS)
r_mem=zeros(sizeS)

## Cloud Node Resource Capacity R_cpu[node], R_mem[node]
R_cpu=zeros(sizeN)
R_mem=zeros(sizeN)

## QoS threshold QoS[service]
QoS=zeros(sizeS)

## Image size img[service]
img=zeros(sizeS)

#Input/Output size
input=zeros(sizeS)
output=zeros(sizeS)
## synched data size db[service]
db=zeros(sizeS)
## processing time for 1 request
mu=zeros(sizeS)
## alpha, beta
alpha=zeros(sizeS)
beta=zeros(sizeS)

S = zeros(sizeS)
for i in 1:sizeS
    S[i] = service_file[i,1]
    img[i] = service_file[i,2]
    QoS[i] = service_file[i,3]
    db[i] = service_file[i,4]
    input[i] = service_file[i,5]
    output[i] = service_file[i,6]
    r_cpu[i] = service_file[i,7]
    r_mem[i] = service_file[i,8]
    mu[i] = service_file[i,9]
    alpha = service_file[i,10]
    beta = service_file[i,11]
end

## Service request binary Q[userId,service]
user_service_file=readdlm("dataset/user_service.txt")
q=zeros(sizeU, sizeS)
#q[1,1]=1 #user 1 request service 1
#q[2,1]=1 #user 2 request service 1
#q[3,1]=1 #user 3 request service 1
for i in 1:sizeU
    tmp_u_s = convert(Int, user_service_file[i,2])
    q[i,tmp_u_s]=1
end
#println(q)

## Network delay D[nodeId, nodeId]
delay_file=readdlm("dataset/network_delay.txt")
D=zeros(sizeN, sizeN)

for i in 1:sizeN
    for j in 1:sizeN
        D[i,j] = delay_file[i,j]
    end
end
#println(D)

println("=============================")
println("Number of service: ", sizeS)
println("Number of user: ", sizeU)
println("Number of BS: ", sizeN)
println("=============================")

##########################################
function centralized_solver()
    #model = Model(solver=SCSSolver())
    model = Model(solver=GLPKSolverMIP())

    #variables
    @variable(model, x[1:sizeS, 1:sizeN, 1:sizeU], Bin)  # Define a variable for routing decision
    @variable(model, y[1:sizeS, 1:sizeN], Bin)  # Define a variable for service deployemnt

    #objective
    @objective(model, Min, sum(y[s,n] for s=1:sizeS for n=1:sizeN))

    #constraints a
    for u in 1:sizeU, s in 1:sizeS
        @constraint(model, q[u,s]*((input[s]+output[s])/r_u[u] +
                                sum(x[s,n,u]*(input[s]+output[s])*D[getBsIndexFromUser(u),n] for n=1:sizeN) +
                                sum(x[s,n,u]*y[s,n]*(mu[s] + beta[s] + alpha[s]*sum(x[s,n,v] for v=1:sizeU)) for n=1:sizeN)+
                                sum(x[s,n,u]*(sum(sum(x[s,m,v]*db[s]*D[getBsIndexFromUser(v),m] for m=1:sizeN) for v=1:sizeG)) for n=1:sizeN)
                                ) <= QoS[s])
    end
    #constraints b
    for n in 1:sizeN
        @constraint(model, sum(y[s,n]*r_cpu[s] for s=1:sizeS) <= 10)
        @constraint(model, sum(y[s,n]*r_mem[s] for s=1:sizeS) <= 10)
    end
    #constraints c
    for u in 1:sizeU, s in 1:sizeS
        @constraint(model, sum(x[s,n,u] for n=1:sizeN) == q[u,s])
    end
    #constraint d
    for u in 1:sizeU, s in 1:sizeS, n in 1:sizeN
        @constraint(model, x[s,n,u] <= y[s,n])
    end

    solve(model)
    #println((getvalue(x)))
    println("====== Service Deployment Result =======")
    #println((getvalue(y)))
    return getvalue(x),getvalue(y)
end

function totalCost(y)
    cost = 0
    for i in 1:sizeS, j in 1:sizeN
        cost = cost + y[i,j]
    end
    return cost
end

#x,y=centralized_solver()
println(y)
cost = totalCost(y)
println("Total Cost = ", cost)
for n in 1:sizeN, s in 1:sizeS
    println(y[s,n]*r_cpu[s], y[s,n]*r_mem[s])
end
