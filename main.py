import json
import csv
from collections import defaultdict

import numpy as np
import pandas
import geopandas
import folium
import matplotlib.pyplot as plt

#from graph_tools import district_subgraph
#from graph_tools import is_district_connected
#from graph_tools import connected_districts

from cko_solver import CKOscheme
from OT_solver import solve_Kantorovich

"""

We are trying to take a shapefile that contains population and geometric data (a graph) to 

1) generate a rectangular distance matrix from the entire graph to a subset of k elements (k<100 nodes)
2) use that rectangular data to set up an optimal transport problem from N elements to k elements (N is the size of the graph)
3) save the resulting map as a column in the shapefile 

"""

class OptimizationError(Exception):
    pass

#The next two functions take an optimal transport plan and if possible returns a corresponding optimal transport map
def plan_to_map(matrix,geoidtable):
    districts = dict()
    maxes = np.argmax(matrix, axis=1)
    #districts.update({'GEOID':'DISTRICT'})
    for i, j in enumerate(maxes):
        #districts.update({geoidtable[i][7:]:j}) #Only for data sets that have a GEOID suffix
        districts.update({geoidtable[i]:j})
    return districts    

def plan_to_map_array(size,matrix,geoidtable):
    districts = np.empty(size)
    maxes = np.argmax(matrix, axis=1)
    #districts.update({'GEOID':'DISTRICT'})
    for i, j in enumerate(maxes):
        #districts.update({geoidtable[i][7:]:j}) #Only for data sets that have a GEOID suffix
        districts[i] = j
    return districts    

def compute_districts(cost_matrix, population,number_of_districts,geoidtable):
    total_population = np.sum(population, dtype=float)
    N = len(population)
    
    c = cost_matrix.astype(float)
    print("Cost matrix shape: ")
    print(c.shape)
    
    mu = population
    #mu = population / total_population # Matthew added the divison for normalization
    nu = np.ones(number_of_districts) * total_population/number_of_districts
    #nu = np.ones(number_of_districts) / number_of_districts # Matthew divide for normalize
    print('µ (mu):')
    print(mu.shape)
    print(mu)
    print('𝜈 (nu):')
    print(nu.shape)
    print(nu)
    solution = solve_Kantorovich(c, mu, nu)
    print("Solution length: ")
    print(len(solution))

    ot_matrix = np.reshape(solution, (c.shape)) # This is the optimal transport plan pi (a matrix)
    
    print('')
    print('Number of source points assigned to multiple targets')
    print(count_splitting(ot_matrix))
    print('')
    print('Number of tracts/groups/blocks with non-zero population')
    print(count_support(mu))    
    print('')
        
    #print(ot_matrix)    
        
    ####Transport plan to transport map
    #districts = plan_to_map(ot_matrix,geoidtable)
    districts = plan_to_map_array(N,ot_matrix,geoidtable)
        
    return districts
  
#This function counts the number of points in the source assigned to multiple targets  
def count_splitting(coupling_matrix):  
  
  M = len(coupling_matrix[0])
  N = len(coupling_matrix[:,0]) #this is the large one
 
  counter = int(0)
  
  for i in range(N):
      if len(np.nonzero(coupling_matrix[i])[0])>1:
          counter = counter+1
  return counter

#This function counts the number of non-zero rows in a distribution
def count_support(distribution):  
  
  N = len(distribution)
  
  counter = int(0)
  
  for i in range(N):
      if distribution[i] > 0:
          counter = counter+1
  return counter
  
def geopdvisual(dataframe,plan_id,outfilename,color):
    dataframe.plot(column = plan_id, cmap = color)
    plt.savefig(outfilename, dpi = 1000)
    
def foliumvisual():    
    pass
    
    ####The first line is commented out if an appropriate json file already exists
    #data_provider.shpconverter(filepath,'testjson.json')
    
    #map_center = np.array([pandas.to_numeric(df['INTPTLAT']).mean(),pandas.to_numeric(df['INTPTLON']).mean()])

    #m = folium.Map(location = [map_center[0],map_center[1]], zoom_start = 8, tiles='Mapbox bright' )        
    
    #m.choropleth(??)
    #m.save(outfile = "foliumtest.html")
  
        
def DetermineDistricts(N,c,alpha):

    nodes = np.empty((N,len(alpha)))
    district = np.empty(N)
    for i in range(N):
        Aff = []
        for j in range(len(alpha)):
            Aff.append(c[i,j]-alpha[j])
        nodes[i] = Aff
        
    for i in range(N):
        counter = 0
        
        Aff = nodes[i]
        for j in range(len(alpha)):
            if Aff[j] == min(Aff):
                counter = counter + 1
        if counter == 1:
                district[i] = np.argmin(np.asarray(Aff))
        if counter >1:
                district[i] = len(alpha)

    return district 
  
def demo_optimal_transport(num_districts=36, state_geoid=42):

    filepath = "./data/TX_vtds.shp"    
    #filepath = "./data/2012_42_test.shp"
    metric_path = "./data/TX_vtds_%s.json" %num_districts
    #metric_path = "./data/VA_%s.json" %num_districts    
    #metric_path = "./data/VA_qed_k%s_1.json"  % num_districts
    #metric_path = "./data/graph_distance_2012_42_test_k16.json"    
    plan_id = 'gedk%s_1' % num_districts                
                    
    num_districts= 36
                    
    #### Sets up data frame from shapefile
    df = geopandas.read_file(filepath) 
    
    #### Temporary population gatherer
    #tracts_graph = data_provider.tracts_graph(state_geoid)
    #data = tracts_graph.data #data is a dictionary     
    #population = np.array([feature["properties"]["2013_population_estimate"] for feature in data["features"]], dtype=float)   
        
    #block_id_table = df['GEOID'].values            
    
    #### Reads population and GEOID table from the dataframe
    block_id_table = df['FIPS'].values        
    population = df['TOTPOP'].values.astype(float)
    #population = np.ones(len(block_id_table))
    Pop_district = np.sum(population) / num_districts
    g = np.ones(num_districts) * Pop_district
    delta = Pop_district * (0.01)    
    
    #import code
    #code.interact(local=locals())
       
    #### Loads distance matrix     
    with open(metric_path) as f:
      data_cost_matrix = np.array(json.load(f))
                 
    #  #### Call compute_districts to solve OT problem
    # num_districts_k = num_districts
    # alpha = np.zeros(num_districts_k)
    # alpha[0] = 50000
    # print('Length of population: ')
    # print(len(population))
    # print('CKOscheme... alpha.shape: ')
    # print(alpha.shape)
    # print('Population shape: ')
    # print(population.shape)
    # print('matrix.shape: ')
    # print(data_cost_matrix.shape)
    # alpha_sol = CKOscheme(alpha,population,g,data_cost_matrix,delta)
 
    # Matthew added actual call since it wasn't present:
    districts = compute_districts(data_cost_matrix, population, num_districts, state_geoid)
    
    # districts = DetermineDistricts(len(population), data_cost_matrix, alpha_sol)
        
    # districts = districts
    #### Adds the districts to the data frame
    df[plan_id] = pandas.Series(districts)
        
    #### Writes changes to the data frame file 
    #df.to_file(filepath)
    
    #sub_df = district_subgraph(df,plan_id,district_number)    
    #list = connected_districts(df,plan_id,num_districts)

    ####Now we call a function that generates and saves a map
    ####it uses geopandas, geojson, and pyplot
    geopdvisual(df,plan_id,'VA_qed_k%s_1_unif.png' % num_districts,'tab20')
    plt.show()
    #with open('test.csv', 'w') as f:
        #fieldnames = ['GEOID', 'DISTRICT']
        #writer = csv.DictWriter(f, fieldnames=fieldnames)
        #writer.writeheader()
        #data = [dict(zip(fieldnames, [k, v])) for k, v in districts.items()]
        #writer.writerows(data)

    import code
    code.interact(local=locals())


    #return districts
    
def main():
    
    #### Starts interactive python terminal
    #import code
    #code.interact(local=locals())
    
    demo_optimal_transport()

if __name__ == '__main__':
    main()
