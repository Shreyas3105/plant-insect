import networkx as nx
import numpy as np
import random
import pandas as pd
import itertools
import argparse
import os



df1 = pd.read_csv('Interaction_Files/Butterfly_Plant/mutualistic.csv')
df2 = pd.read_csv('Interaction_Files/Butterfly_Plant/antagonistic.csv')

count = pd.read_csv('Interaction_Files/Butterfly_Plant/Plant_Count.csv')

def filter_df_by_count_na(df1, count):
    # Identify values in "taxon" column with null "Freq" in the count DataFrame
    na_freq_taxa = count.loc[count['Count'].isna(), 'Plant']

    # Remove rows in df1 where "Target_Name" matches the identified values
    df1_filtered = df1[~df1['Target_Name'].isin(na_freq_taxa)]

    return df1_filtered

# Database truncation
df1 = filter_df_by_count_na(df1, count)
df2 = filter_df_by_count_na(df2, count)

# Differentiating the two networks
df1["Interaction_Strength"] = 1
df2["Interaction_Strength"] = -1

# Combining the networks
df_combined = pd.concat([df1, df2], ignore_index=True)

# Network creation
G = nx.from_pandas_edgelist(df_combined, 'Source_Name', 'Target_Name', create_using=nx.DiGraph())

# Metrics
def rich_club_coefficient(graph):
    degrees = dict(graph.degree())
    rich_club_coefficients = []

    for node in graph.nodes():
        neighbors = list(graph.neighbors(node))
        node_degree_sum = sum(degrees[neighbor] for neighbor in neighbors)

        if node_degree_sum > 0:
            rich_club_coefficients.append(len(neighbors) / node_degree_sum)

    average_rich_club = sum(rich_club_coefficients) / len(rich_club_coefficients) if rich_club_coefficients else 0
    return average_rich_club

def calculate_NODF(graph):

    # Convert partitions to lists
    partition1 = [node for node in graph.nodes if graph.out_degree(node) > 0]
    partition2 = [node for node in graph.nodes if graph.in_degree(node) > 0]

    # Calculate NODF
    nodf = 0
    for node in partition1:
        neighbors = set(graph.neighbors(node))
        overlap_i = len(neighbors.intersection(partition2))
        dissimilarity_i = len(neighbors.difference(partition2))
        if len(partition2) == 0:
            return 0
        nodf_i = overlap_i / len(partition2)
        nodf += (nodf_i - dissimilarity_i / len(partition2))
        
    if len(partition1) != 0:
        nodf /= len(partition1)
        return nodf
    return 0

def modularity_metric(graph):
    if len(graph) <= 2:
        return 1
    # Detect communities using greedy algorithm
    community = nx.community.greedy_modularity_communities(graph)

    # Calculate modularity
    modularity = nx.community.modularity(graph, community)
    return modularity

def graph_size(graph):
    return len(graph)

def total_avg_degree(graph):
    if len(graph) <= 1:
        return 0
    return np.mean(list(dict(graph.degree()).values()))

def pagerank(graph):
    if len(graph) <= 1:
        return 0
    return np.mean(list(dict(nx.pagerank(graph)).values()))

def calculate_connectance(graph):
    if len(graph.nodes) != 0:
        return len(graph.edges()) / (len(graph.nodes()) ** 2)
    return 0

def calculate_number_of_weakly_connected_components(graph):
    return len(list(nx.weakly_connected_components(graph)))

all_metrics_functions = [
    ("Graph Size", graph_size),
    ("Average Degree", total_avg_degree),
    ("Average PageRank", pagerank),
    ("Modularity", modularity_metric),
    ("Nestedness", calculate_NODF),
    ("Rich Club Coefficient", rich_club_coefficient),
    ("Connectance", calculate_connectance),
    ("Number of Weakly Connected Components", calculate_number_of_weakly_connected_components),
]


# Define the list of metrics you want to calculate in your run
selected_metrics = [
    "Graph Size",
    "Nestedness",
    "Connectance",
    "Average Degree",
    "Average PageRank",
    "Rich Club Coefficient",
    "Number of Weakly Connected Components",
    "Modularity"
]

# Filter the metrics_functions list based on your selection
metrics_functions = [(name, func) for name, func in all_metrics_functions if name in selected_metrics]

# Create lists to store results for each run
all_metrics_data = {metric_name: [] for metric_name, _ in metrics_functions}  # Store metrics data for each run

def perturbation_run_single(run_index, graph, selection_strategy, removal_strategy):
    # Create a copy of the original graph for this run
    perturbed_network = graph.copy()
    random.seed(run_index)
    # Create lists to store results for this run
    metrics_data = {metric_name: [] for metric_name, _ in metrics_functions}
    primary_extinctions = 0
    # Iterate through the removal of pollinators
    while perturbed_network.number_of_nodes() > 0: 
        # Define pollinator_nodes based on out-degrees
        pollinator_nodes = [node for node in perturbed_network.nodes if perturbed_network.out_degree(node) > 0]
        # Define plant_nodes based on in-degrees
        plant_nodes = [node for node in perturbed_network.nodes if perturbed_network.in_degree(node) > 0]

        if selection_strategy == "Random":
            if removal_strategy == "Pollinator":
                if not pollinator_nodes:
                    break  # No more pollinators left
                # Randomly select one pollinator to remove
                node_to_remove = random.choice(pollinator_nodes)

            if removal_strategy == "Plant":
                if not plant_nodes:
                    break  # No more plants left
                # Randomly select one pollinator to remove
                node_to_remove = random.choice(plant_nodes)

            if removal_strategy == "Random":
                if not perturbed_network.nodes:
                    break  # No more plants or pollinators left
                # Randomly select one pollinator to remove
                node_to_remove = random.choice(list(perturbed_network.nodes))

        if selection_strategy == "Generalist":
            if removal_strategy == "Pollinator":
                if not pollinator_nodes:
                    break  # No more pollinators left
                # Find pollinators with the highest out-degree
                max_out_degree = max(perturbed_network.out_degree(pollinator_nodes), key=lambda x: x[1])[1]
                highest_degree_pollinators = [node for node, out_degree in perturbed_network.out_degree(pollinator_nodes) if out_degree == max_out_degree]
                # Randomly select one from the highest degree pollinators
                node_to_remove = random.choice(highest_degree_pollinators)

            if removal_strategy == "Plant":
                if not plant_nodes:
                    break  # No more plants left
                # Find plants with the highest in-degree
                max_in_degree = max(perturbed_network.in_degree(plant_nodes), key=lambda x: x[1])[1]
                highest_degree_plants = [node for node, in_degree in perturbed_network.in_degree(plant_nodes) if in_degree == max_in_degree]
                # Randomly select one from the highest degree plants
                node_to_remove = random.choice(highest_degree_plants)

            if removal_strategy == "Random":
                if not plant_nodes and not pollinator_nodes:
                    break  # No more plants or pollinators left
                # Combine all nodes (plants and pollinators)
                all_nodes = list(perturbed_network.nodes)
                # Find nodes with the highest degree (both in-degree and out-degree)
                max_degree = max(perturbed_network.degree(all_nodes), key=lambda x: x[1])[1]
                highest_degree_nodes = [node for node, degree in perturbed_network.degree(all_nodes) if degree == max_degree]
                # Randomly select one from the highest degree nodes
                node_to_remove = random.choice(highest_degree_nodes)
                
        if selection_strategy == "Specialist":
            if removal_strategy == "Pollinator":
                if not pollinator_nodes:
                    break  # No more pollinators left
                # Find pollinators with the highest out-degree
                min_out_degree = min(perturbed_network.out_degree(pollinator_nodes), key=lambda x: x[1])[1]
                lowest_degree_pollinators = [node for node, out_degree in perturbed_network.out_degree(pollinator_nodes) if out_degree == min_out_degree]
                # Randomly select one from the highest degree pollinators
                node_to_remove = random.choice(lowest_degree_pollinators)

            if removal_strategy == "Plant":
                if not plant_nodes:
                    break  # No more plants left
                # Find plants with the highest in-degree
                min_in_degree = min(perturbed_network.in_degree(plant_nodes), key=lambda x: x[1])[1]
                lowest_degree_plants = [node for node, in_degree in perturbed_network.in_degree(plant_nodes) if in_degree == min_in_degree]
                # Randomly select one from the highest degree plants
                node_to_remove = random.choice(lowest_degree_plants)

            if removal_strategy == "Random":
                if not plant_nodes and not pollinator_nodes:
                    break  # No more plants or pollinators left
                # Combine all nodes (plants and pollinators)
                all_nodes = list(perturbed_network.nodes)
                # Find nodes with the highest degree (both in-degree and out-degree)
                min_degree = min(perturbed_network.degree(all_nodes), key=lambda x: x[1])[1]
                lowest_degree_nodes = [node for node, degree in perturbed_network.degree(all_nodes) if degree == min_degree]
                # Randomly select one from the highest degree nodes
                node_to_remove = random.choice(lowest_degree_nodes)

        if selection_strategy == "Count":
            if not plant_nodes:
                break
            # Find remaining plant nodes in the network
            remaining_plant_nodes = set(plant_nodes)
            # Filter count dataframe to include only remaining plant nodes
            count_subset = count[count['Plant'].isin(remaining_plant_nodes)]
            count_subset_sorted = count_subset.sort_values(by='Count')
            # Find all plants with the minimum count
            min_count_plants = list(count_subset_sorted[count_subset_sorted['Count'] == count_subset_sorted['Count'].min()]['Plant'])
            # Randomly select one plant if there are multiple with the minimum count
            min_count_node = random.choice(min_count_plants)
            node_to_remove = min_count_node

        # Remove the selected pollinator
        perturbed_network.remove_node(node_to_remove)
        # Find isolated nodes (nodes with no connections)
        isolates = list(nx.isolates(perturbed_network))
        # Remove isolated nodes
        perturbed_network.remove_nodes_from(isolates)

        # Calculate modularity after every 10 steps
        if primary_extinctions % 10 == 0:
            print(f'{len(plant_nodes)} plants left. {len(pollinator_nodes)} pollinators left')
            for metric_name, metric_function in metrics_functions:
                metric_value = metric_function(perturbed_network)
                metrics_data[metric_name].append(metric_value)

        primary_extinctions += 1
    

    return metrics_data

def run_simulation(parameters, task_id):
    print(f'Starting simulation with parameters: {parameters}')
    
    # df = parameters["CHOSEN_DF"]
    selection = parameters["selection_strategy"]
    removal = parameters["removal_strategy"]

    simulation = perturbation_run_single(run_index = task_id,
                                        graph = G,
                                        selection_strategy = selection,
                                        removal_strategy = removal,
                                        )
    df = pd.DataFrame(simulation)

    save_path = f'Results/Combined_Butterfly/{selection}_{removal}/'
    filename = f"Perturbation_{selection}_{removal}_{str(task_id).zfill(4)}_butterflyinsect.csv"
    
    # Check if the directory exists, and create it if not
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # Combine the save path and filename
    file_path = os.path.join(save_path, filename)
    
    df.to_csv(file_path, index=False)

def product_dict(**kwargs):
    keys = kwargs.keys()
    for instance in itertools.product(*kwargs.values()):
        yield dict(zip(keys, instance))


script_params = [
    {
        'selection_strategy': ['Random', 'Generalist', 'Specialist', 'Count'],
        'removal_strategy': ['Plant'],
        'ID': list(range(1000))  
    }
]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--task_id', type=int)
    args = parser.parse_args()
    task_id = args.task_id
    parameter_combinations = list(
        itertools.chain.from_iterable([
            list(product_dict(**script_params[i]))
            for i in range(len(script_params))
        ]))
    print(len(parameter_combinations))
    run_simulation(parameter_combinations[task_id], task_id)