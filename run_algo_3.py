from scipy.optimize import minimize, milp, Bounds, LinearConstraint
from opti_classes import *
from statistics import mean
import numpy as np
import pandas as pd
import os

#------------------------------------------------------------------ Input data --------------------------------------------------------------------#

directory = r"\\smb.isibhv.dmawi.de\projects\p_Bionik_memb\Projekte\04 Laufende Projekte\Fiona\AP 220 - Elise Development\Toolpath_Generation_Algorithm"
inputDataFileName = "graph_data_3.xlsx"
inputDataFilePath = os.path.join(directory, inputDataFileName) 
sheet_names = ["sheet1", "sheet2", "sheet3"]    # Available sheets 
Fw = 0.1          # Fiber width
eta = 1           # Set exponential value of history term
num_layers = 10   # Set number of layers
use_constraint = False
write_output = False

# Set up dict for history layers the layer ids for required layers
history_layers = {}
layers = dict((x, {"id" : x}) for x in range(num_layers))

#--------------------------------------------------------------- Important functions --------------------------------------------------------------#

# Cost function for milp without history
def costFunction(x, A, Y1):
    error = Y1 - np.matmul(A, x)
    cost = np.linalg.norm(error, 2)
    return cost

# Cost function for milp with history
def costFunction_hist(A_curr, Y_curr, hist_layers, n, eta):
    # Extract the history and calculate summation term
    summation_hist = np.zeros(len(Y_curr))
    for i in range(n):
        A_hist = hist_layers[i]["A"]
        x_hist = hist_layers[i]["x"]
        summation_hist = summation_hist + np.matmul(A_hist, x_hist)
    # Calculate the coefficients of cost function
    history_term = np.power((n+1)*Y_curr - summation_hist, eta)
    cost = -np.transpose(A_curr).dot(history_term)
    return cost

# Function to calculate optimal fiber orientation
def solve_layer(current_layer, history_layers, Fw, eta, sheet_names):
    # Extract current layer's sheet
    n = current_layer["id"]
    
    sheet_costs = {}
    for i in range(len(sheet_names)):
        sheet = createNewSheet(inputDataFilePath, sheet_names[i])
    
        # Get connection loop relations; A - connection-loop matrix, Y1 - connection target
        A, Y1 = sheet.get_connection_loop_relations(Fw)
        
        # Get edge loop relation; B : edge-loop matrix, Y2 : edge target
        B, Y2 = sheet.get_edge_loop_relations(Fw)
        
        # Get edge loop relation; V : vertex-loop matrix, Y3 : vertex target
        V, Y3 = sheet.get_vertex_loop_relations(Fw)
        
        # Define optmisation bounds
        loops_list = [sheet.loops[x] for x in sorted(sheet.loops.keys())]
        lowerBounds = np.zeros(len(loops_list))
        #upperBounds = np.array([(x.target/Fw) for x in loops_list])
        upperBounds = np.array([0.5*(x.target/Fw) for x in loops_list])
        boundData=Bounds(lowerBounds,upperBounds)
        
        # Define linear constraints for optimisation
        if not use_constraint:
            # Optimisation linear constraints without vertex constraints
            matrixLinearConstraints = B
            lowerBoundLinearConstraints = np.zeros(len(Y2))
            upperBoundLinearConstraints = Y2
            linearConstraints = LinearConstraint(matrixLinearConstraints, lowerBoundLinearConstraints, upperBoundLinearConstraints)
        else:
            # Optimisation linear constraints with vertex constraints
            matrixLinearConstraints = np.vstack((B, V))
            lowerBoundLinearConstraints = np.zeros(len(Y2) + len(Y3))
            upperBoundLinearConstraints = np.hstack((Y2, Y3))
            linearConstraints = LinearConstraint(matrixLinearConstraints, lowerBoundLinearConstraints, upperBoundLinearConstraints)
        
        # Minimize cost function using milp
        cost_hist = costFunction_hist(A, Y1, history_layers, n, eta)
        integrality = np.ones_like(cost_hist)
        result = milp(cost_hist, constraints=linearConstraints, bounds=boundData, integrality=integrality)
        x = np.array(result.x)
        
        # Calculate cost for current layer
        cost = - cost_hist.dot(x)
        
        # Calculate value of Ax for current layer
        Ax = A.dot(x)
        
        # Fill in the calculated 
        sheet_costs[i] = {"sheet_name" : sheet_names[i], "cost" : cost ,
                          "A" : A , "B" : B, "Y1" : Y1 ,"x" : x, "Ax" : Ax}
                          
    # Save the values in a dictionary
    cost_key = max(sheet_costs, key = lambda k : sheet_costs[k]["cost"])
    history_layers[n] = {"sheet_name" : sheet_costs[cost_key]["sheet_name"],
                         "A" : sheet_costs[cost_key]["A"], "B" : sheet_costs[cost_key]["B"],
                         "Y1" : sheet_costs[cost_key]["Y1"], "x" : sheet_costs[cost_key]["x"],
                         "cost" : sheet_costs[cost_key]["cost"], "Ax" : sheet_costs[cost_key]["Ax"]}
                         
    # Layer outputs
    x = sheet_costs[cost_key]["x"]
    sum_connections = np.matmul(sheet_costs[cost_key]["A"], sheet_costs[cost_key]["x"])
    cost = sheet_costs[cost_key]["cost"]
    
    return x, sum_connections, cost

#---------------------------------------------------------------- Run algorithm ------------------------------------------------------------------#

# Solve all layers
cum_summation = 0
cumulative = np.zeros((33, len(layers)))
print("Cost and assigned sheet names for all layers : ")
for i in range(len(layers.keys())):
    x, cum_sum, cost = solve_layer(layers[i], history_layers, Fw, eta, sheet_names)
    print("cost of layer " + str(i) + ":", f'{cost : .2f}', "sheet name : ", history_layers[i]["sheet_name"])
    print("Y :", history_layers[i]["Y1"])
    print("Ax :", history_layers[i]["Ax"])
    print("Y-Ax :", history_layers[i]["Y1"] - history_layers[i]["Ax"])
    print("\n")
    
    cum_summation = cum_summation + cum_sum
    cumulative[:, i] = cum_summation
    if (i == len(layers.keys())-1):
        print("Cumulative strength of all connections :\n", cum_summation)
        print("Max  : ", round(max(cum_summation)), "\nMean : ", round(mean(cum_summation)))
   
#----------------------------------------------------------------- Write output ------------------------------------------------------------------#     

# Print and save the solution in an xl sheet
solutionFileName = "solution_" + inputDataFileName
solutionFilePath = os.path.join(directory, solutionFileName)

if write_output:
    writer = pd.ExcelWriter(solutionFilePath)
    for i in range(len(layers.keys())):
        # Write an xl sheet of layer solution
        xlSheetName = "layer_" + str(i)
        solution_df =  pd.DataFrame()
        solution_df[history_layers[i]["sheet_name"]] = np.asarray(x, dtype = 'int')
        solution_df.to_excel(writer, sheet_name = xlSheetName, index = False)
    writer.close()
        