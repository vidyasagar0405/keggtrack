from Bio.KEGG.KGML import KGML_parser
import networkx as nx
import csv
from os import path

# Define the source and target metabolites
source_metabolite = "C17213"
target_metabolite = "C08412"

# File paths
kgml_file = "./kgml/rsz/rsz00966.kgml"
file_basename = path.basename(kgml_file)
pathway_info = f"{file_basename.split('.')[0]}_pathway_info.tsv"
tracked_pathway_info = f"{file_basename.split('.')[0]}_tracked_pathway_info.tsv"


def load_kgml(kgml_file):
    """Load and parse the KGML file."""
    return KGML_parser.read(open(kgml_file, "r"))


def build_graph(pathway):
    """Build a directed graph from KGML data."""
    G = nx.DiGraph()
    reaction_to_genes = {}

    for reaction in pathway.reactions:
        reaction_name = reaction.name
        reaction_genes = set()

        substrates = [s.name for s in reaction.substrates]
        products = [p.name for p in reaction.products]

        # Add nodes and edges
        for sub in substrates:
            G.add_node(sub)
        for prod in products:
            G.add_node(prod)
        for sub in substrates:
            for prod in products:
                G.add_edge(sub, prod, reaction=reaction_name)

        # Find genes linked to reaction
        for entry in pathway.entries.values():
            if entry.type == "gene" and str(reaction.id) in entry.name:
                reaction_genes.add(entry.name)

        reaction_to_genes[reaction_name] = reaction_genes

    return G, reaction_to_genes


def write_pathway_info(graph, reactions, output_file):
    """Write reaction, substrate, product, and genes to TSV."""
    with open(output_file, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["id", "Reaction_id", "Substrate", "Product", "Genes_involved"])

        for reaction, genes in reactions.items():
            substrates = []
            products = []

            for sub, prod, data in graph.edges(data=True):
                if data["reaction"] == reaction:
                    substrates.append(sub)
                    products.append(prod)

            writer.writerow([reaction, ",".join(set(substrates)), ",".join(set(products)), ",".join(genes) if genes else "No genes found"])


def find_metabolite_path(graph, source, destination):
    """Find the shortest path between two metabolites."""
    if not source.startswith("cpd:"):
        source = f"cpd:{source}"
    if not destination.startswith("cpd:"):
        destination = f"cpd:{destination}"

    if source not in graph or destination not in graph:
        print(f"⚠️ Error: {source} or {destination} is missing from the graph!")
        return None

    try:
        return nx.shortest_path(graph, source=source, target=destination)
    except nx.NetworkXNoPath:
        print(f"🚫 No pathway found between {source} and {destination}.")
        return None


def get_genes_in_path(graph, path, reaction_to_genes):
    """Extract genes involved in a given path."""
    genes = set()
    for i in range(len(path) - 1):
        edge_data = graph.get_edge_data(path[i], path[i + 1])
        if edge_data:
            reaction_name = edge_data["reaction"]
            genes.update(reaction_to_genes.get(reaction_name, []))
    return genes


def write_traced_path(graph, path, output_file, reaction_to_genes):
    """Write traced pathway information to TSV."""
    with open(output_file, "w", newline="") as pathfile:
        path_writer = csv.writer(pathfile, delimiter="\t")
        path_writer.writerow(["Step", "Metabolite", "Next_Metabolite", "Reaction_id", "Genes_Involved"])

        for i in range(len(path) - 1):
            sub = path[i]
            prod = path[i + 1]
            edge_data = graph.get_edge_data(sub, prod)
            if edge_data:
                reaction_name = edge_data["reaction"]
                genes = reaction_to_genes.get(reaction_name, "No genes found")
                path_writer.writerow([i + 1, sub, prod, reaction_name, ",".join(genes)])


def main():

    pathway = load_kgml(kgml_file)
    G, reaction_to_genes = build_graph(pathway)

    # Write full pathway info
    write_pathway_info(G, reaction_to_genes, pathway_info)

    # Find and write the traced path
    path = find_metabolite_path(G, source_metabolite, target_metabolite)
    if path:
        print(f"\nMetabolite Pathway from {source_metabolite} to {target_metabolite}:")
        print(" → ".join(path))
        genes_involved = get_genes_in_path(G, path, reaction_to_genes)
        print(f"\nGenes involved in the pathway: {', '.join(genes_involved) if genes_involved else 'No genes found'}")

        write_traced_path(G, path, tracked_pathway_info, reaction_to_genes)
    else:
        print(f"No pathway found from {source_metabolite} to {target_metabolite}.")

if __name__ == "__main__":
    main()
