import gurobipy as gp, numpy as np, math
import sys, os, argparse

def bi_cluster(M,alpha,beta):

    n,m = M.shape
    n_cells, n_sites = math.floor(alpha*n), math.floor(beta*m)

    model = gp.Model()
    model.Params.Threads = args.threads
    model.Params.TimeLimit = args.run_time

    # initialize variables
    A = np.empty((n,m), dtype=object)
    C = np.empty(n, dtype=object)
    S = np.empty(m, dtype=object)
    for i in range(n):
        for j in range(m):
            A[i][j] = model.addVar(vtype=gp.GRB.BINARY)
    for i in range(n):
        C[i] = model.addVar(vtype=gp.GRB.BINARY)
    for j in range(m):
        S[j] = model.addVar(vtype=gp.GRB.BINARY)

    # set constraints
    for i in range(n):
        for j in range(m):
            model.addConstr(A[i][j] <= C[i])
            model.addConstr(A[i][j] <= S[j])
            model.addConstr(C[i]+S[j]-1 <= A[i][j])
    model.addConstr(gp.quicksum(C[i] for i in range(n)) == n_cells)
    model.addConstr(gp.quicksum(S[j] for j in range(m)) == n_sites)

    model.setObjective(gp.quicksum(A[i][j]*M[i][j] for i in range(n) for j in range(m)),
                       gp.GRB.MAXIMIZE)

    model.optimize()

    return np.array([i for i in range(n) if C[i].X > 0]), \
           np.array([j for j in range(m) if S[j].X > 0])


if __name__=="__main__":
    #parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('-a', '--alpha', type=float, required=True)
    parser.add_argument('-b', '--beta', type=float, required=True)
    parser.add_argument('-c', '--threads', type=int, required=True)
    parser.add_argument('-t', '--run_time', type=int, required=True)
    args = parser.parse_args(sys.argv[1:])

    # biclustering

    assert args.alpha <= 1 and args.beta <=1

    m = np.load(args.input, allow_pickle=True)['m']

    # run bicluster
    cells_idx, sites_idx = bi_cluster(m,args.alpha,args.beta)
    np.savez(args.output, rows=cells_idx, cols=sites_idx)

