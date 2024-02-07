import time
import numpy as np
import gurobipy as grb
#from gurobipy import GRB
GRB = grb.GRB

"""
Sets up and solves the Kantorovich problem using the Gurobi optimizer.
"""

def sparse_constraints_dual(M, N):
    constraints = []

    for i in range(M):
        for j in range(N):
            constraints.append([[str(i), str(M + j)], [1.0, 1.0]])
    return constraints

def sparse_constraints_dual_rhs(cost, M, N):
    rhs = []

    for i in range(M):
        for j in range(N):
            rhs.append(cost[i, j])
    return rhs

def sparse_constraints_first(M, N):  # this one produces constraints regarding the first marginal
    constraints = []

    for j in range(M):
        constraints.append([["{0}".format(j * N + k) for k in range(N)], [1] * N])

    return constraints

def sparse_constraints_second(M, N):  # this one produces constraints regarding the second marginal
    constraints = []

    for i in range(N):
        constraints.append([["{0}".format(k * N + i) for k in range(M)], [1] * M])

    return constraints

## New code

def solve_Kantorovich(cost, mu, nu, epr=0, JustObjValue=False):
    # sanity checking... make sure all districts are assigned a polling center
    # and that all polling centers are assigned a district
    districtList = []
    centerList = []
    totalAssigned = 0;
    
    # an effort to translate this to the academic formulation of the problem
    model = grb.Model()
    warehouses = range(len(mu))
    customers = range(len(nu))
    print('warehouses: ')
    print(warehouses)
    print('customers: ')
    print(customers)
    to_ship = get_ship_vars(model, warehouses, customers)
    model.update()
    capacity_constrs = get_capacity_constrs(model, to_ship, mu, nu)
    demand_constrs = get_demand_constrs(model, to_ship, mu, nu)
    model.setObjective(grb.quicksum(cost[warehouse][customer]*to_ship[warehouse, customer]
                                for warehouse in warehouses
                                for customer in customers))
    model.optimize()
    if model.Status == GRB.OPTIMAL:
        for (warehouse, customer), var in sorted(to_ship.items()):
            if var.X > 1e-4:
                totalAssigned = totalAssigned + var.X
                if warehouse not in districtList:
                    districtList.append(warehouse)
                if customer not in centerList:
                    centerList.append(customer)
                print ("Send", var.X, "households from district", warehouse, "to polling center", customer)
    else:
        model.dispose()
        raise Exception("Optimization problem with Gurobi.")
    
    print("Total households assigned: ")
    print(totalAssigned)
    print("Unassigned districts: ")
    for i in range(len(mu)):
        if i not in districtList:
            print("Tract/group/block",i,"is not assigned to a target.")
    print("Polling centers assigned: ")
    centerList.sort()
    print(centerList)
    
    return model.X
                                
def get_demand_constrs(model, to_ship, capacities, demands):
    warehouses = range(len(capacities))
    customers = range(len(demands))
    return [model.addConstr(grb.quicksum(to_ship[warehouse, customer] for warehouse in warehouses) == demand,
                            name='demand.' + str(customer))
           for customer, demand in enumerate(demands)]
    
def get_capacity_constrs(model, to_ship, capacities, demands):
    warehouses = range(len(capacities))
    customers = range(len(demands))
    return [model.addConstr(grb.quicksum(to_ship[warehouse, customer] for customer in customers) <= capacity,
                            name='capacity.' + str(warehouse))
            for warehouse, capacity in enumerate(capacities)]
    
def get_ship_vars(model, warehouses, customers):
    return {(warehouse, customer): model.addVar(name='ship.w' + str(warehouse) + '.c' + str(customer))
                                   for customer in customers
                                   for warehouse in warehouses}

## End new code

def solve_Kantorovich_old(cost, mu, nu, epr=0, JustObjValue=False):
    M = len(mu)
    N = len(nu)
    # epr = percentage of error allowed when matching district populations
    nu_top = nu * (1 + epr * (0.01))
    nu_bottom = nu * (1 - epr * (0.01))

    print("")
    print("(solve_Kantorovich)")
    print("")
    
    ## Matthew start new code here
    print("Cost matrix shape: ")
    print(cost.shape)
    c = cost.flatten()
    A1 = np.kron(np.ones((1, M)), np.eye(N))
    A2 = np.kron(np.eye(M), np.ones((1, N)))
    A = np.vstack((A1, A2))
    
    p = np.ones(M)
    q = np.ones(N)
    
    d = np.concatenate((p, q))

    model = gp.Model()

    x = model.addMVar(shape=N*M, vtype=GRB.CONTINUOUS, name="x")
    model.setObjective(c @ x, GRB.MINIMIZE)

    model.addConstr(A @ x == d, name="constraint")
    model.setParam("OutputFlag", 0)
    gurobi_start = time.time()
    model.optimize()
    gurobi_time = time.time() - gurobi_start
    if model.status == GRB.OPTIMAL:
        #print(model.getAttr("x"))
        pi = np.array(model.getAttr("x"))
        print(pi.shape)
        #print(pi)
        print(N)
        pi = np.reshape(pi, (N,M))
        print(pi.shape)
        u = pi[:N]
        print("u")
        print(u)
        v = pi[N:N+M]
        print("v")
        print(v)
        val = model.objVal
        print("val")
        print(val)
    else:
        model.dispose()
        raise Exception("Optimization problem with Gurobi.")
    ## Matthew end new code here
    
    return 0;

    #### Initialize the Gurobi model
    OT_prob = gp.Model()

    #### Add variables to model and objective function
    empty_list = list(range(M * N))
    #Var_names = [for i in range(M * N)]
    Var_names = empty_list
    OT_vars = OT_prob.addVars(Var_names)
    
    #### Set up constraints
    OT_prob.addConstrs(gp.quicksum(OT_vars[i*j] for j in range(1,N)) == mu[i] for i in range(1,M))
    OT_prob.addConstrs(gp.quicksum(OT_vars[i*j] for i in range(1,M)) <= nu_top[j] for j in range(1,N))
    OT_prob.addConstrs(gp.quicksum(OT_vars[i*j] for i in range(1,M)) >= nu_bottom[j] for j in range(1,N))

    #### Set objective function
    OT_prob.setObjective(gp.quicksum(cost[i, j] * OT_vars[i*j] for i in range(M) for j in range(N)), GRB.MINIMIZE)

    #### Find the minimizer
    #OT_prob.setParam('OutputFlag', 0)  # Suppress Gurobi output
    OT_prob.optimize()

    if OT_prob.status == GRB.OPTIMAL:
        #### Reversing the normalization:
        #sol = np.array(OT_prob.getAttr("x"))
        obj = OT_prob.getObjective()
        sol = obj.getValue()
        #sol = np.array([OT_vars[i*j].x for i in range(M) for j in range(N)])

        print('(Kantorovich) Total Cost = ' + str(OT_prob.objVal))

        if JustObjValue == True:
            return sol
            #return OT_prob.objVal
        else:
            return sol
    else:
        OT_prob.dispose()
        raise Exception("Optimization problem with Gurobi.")


def solve_Dual(cost, mu, nu):
    M = len(mu)
    N = len(nu)

    print("")
    print("(solve_Dual)")
    print("")

    #### Initialize the Gurobi model
    OT_prob = gp.Model()

    #### Add variables to model and objective functional
    Var_names = ["{0}".format(i) for i in range(M + N)]
    OT_vars = OT_prob.addVars(Var_names)

    #### Set up constraints
    OT_prob.addConstrs(gp.quicksum(OT_vars[i] for i in range(M * N) if (i // N) == j) <= cost[j] for j in range(M))
    OT_prob.addConstrs(gp.quicksum(OT_vars[i] for i in range(M * N) if (i % N) == j) <= cost[M + j] for j in range(N))

    #### Set objective function
    OT_prob.setObjective(gp.quicksum(mu[i] * OT_vars[i] for i in range(M)) + gp.quicksum(nu[j] * OT_vars[M + j] for j in range(N)), GRB.MAXIMIZE)

    #### Find the maximizer
    OT_prob.setParam('OutputFlag', 0)  # Suppress Gurobi output
    OT_prob.optimize()

    print('(Dual) Total Cost = ' + str(OT_prob.objVal))

    return np.array([OT_vars[i].x for i in range(M + N)])

def main():
    pass

if __name__ == '__main__':
    main()
