'''
Class to store a symmetric distance matrix

@author: Luis Carlos Ospina Restrepo
'''
import copy
from Evolution import Evolution
from Evolution import ToolsToWorkWithSequences

from Sequence import Sequence


class DistanceMatrix(object):
    
    '''
    Requires a list with the name of the species
    '''

    def __init__(self, name_of_species):
        '''
        Constructor
        '''
        self.name_of_species = name_of_species
        # create the list of lists
        self.m = []
        # create for each row the column
        for r in range(len(name_of_species)):
            #
            co = [0]*len(name_of_species)
            #
            self.m.append(co)
    
    
    '''
    Add a value in position i,j (same for j,i)
    '''
    def add_value_i_j(self,i,j,d):
        self.m[i][j] = d
        self.m[j][i] = d
        
    
    '''
    Get the value at position i,j
    '''
    def get_value_i_j(self,i,j):
        return self.m[i][j]
    
    
    '''
    Set a new name at species i
    '''
    def change_name_species_i(self,i,new_name):
        self.name_of_species[i] = new_name

    
    '''
    Get the number of rows
    '''
    def n_rows(self):
        return len(self.m)
    
    
    '''
    Remove the row and column at position i
    '''
    def remove_species_i(self,i):
        #
        self.name_of_species.pop(i)
        #
        self.m.pop(i)
        #
        for j in range(len(self.m)):
            #
            self.m[j].pop(i)
                
    
    '''
    Get the name of the species
    '''            
    def get_name_of_species(self):
        return self.name_of_species
    
    
    '''
    Generate a copy of this distance matrix
    '''
    def copy(self):
        dp = DistanceMatrix(copy.deepcopy(self.name_of_species))
        for i in range(len(self.name_of_species)-1):
            for j in range(i+1, len(self.name_of_species)):
                dp.add_value_i_j(i, j, self.get_value_i_j(i,j))
        return dp
    

class DistanceBasedAlgorithms():

    def __init__(self, d_matrix):
        """
        Initialize the DistanceBasedAlgorithms class with a distance matrix.

        :param d_matrix: A distance matrix representing pairwise distances between species.
        """
        self.d_matrix = d_matrix

    def helper_mymin(self, dm):
        """
        Finds the pair of species with the smallest distance in the matrix.

        :param dm: The distance matrix.
        :return: The indices (min_i, min_j) of the closest pair and their distance (min_d).
        """
        min_d = dm.get_value_i_j(0,1) # Initialize with the first pair's distance
        min_i = 0
        min_j = 1
        for i in range(dm.n_rows()-1): # Loop through each row
            for j in range(i+1, dm.n_rows()): # Loop through each column, starting from the next row
                if dm.get_value_i_j(i,j) < min_d:
                    # Check if the current pair has a smaller distance
                    min_d = dm.get_value_i_j(i,j) # Update smallest distance
                    min_i = i # Update index of first species
                    min_j = j # Update index of second species
        return min_i, min_j, min_d # Return the closest pair and their distance
    
    def UPGMA(self):
        """
        Performs the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) algorithm 
        to construct a phylogenetic tree using hierarchical clustering.

        :return: The final tree in Newick format.
        """
        dm = self.d_matrix.copy() # Copy the distance matrix to avoid modifying the original

        while dm.n_rows() > 1:
            # Find the closest pair of species
            min_i, min_j, min_d = DistanceBasedAlgorithms.helper_mymin(self, dm)

            species_names = dm.get_name_of_species() # Get species names from the matrix

            # Create a new node combining the two closest species, including branch lengths
            new_name = '('+species_names[min_i] +':'+ str(min_d/2) + ','+species_names[min_j]  + ':'+ str(min_d/2)+')'
            # Update the name of the merged cluster

            dm.change_name_species_i(min_i, new_name)
            # Update the distances for the newly formed cluster
            for x in range(dm.n_rows()):
                if x != min_i and x != min_j:
                    d_i = dm.get_value_i_j(min_i, x) # Distance from species min_i to x
                    d_j = dm.get_value_i_j(min_j, x) # Distance from species min_j to x
                    new_distance = (d_i + d_j)/2 # Average distance for UPGMA
                    dm.add_value_i_j(min_i,x, new_distance) # Update new distance
            dm.remove_species_i(min_j)
            # Remove the merged species from the matrix
        return dm.get_name_of_species()[0] # Return the final tree in Newick format


    
    def NJ(self):
        """
        Performs the Neighbor-Joining (NJ) algorithm to construct a phylogenetic tree 
        from a distance matrix.

        :return: The final tree in Newick format.
        """
        dm = self.d_matrix.copy() # Copy the distance matrix
        while dm.n_rows() > 2:
            # Find the closest pair of species using the helper function
            min_i, min_j, min_d = self.helper_mymin(dm)
            species_names = dm.get_name_of_species()

            # Create a new node combining the closest species (without explicit branch lengths yet)
            new_name = f"({species_names[min_i]},{species_names[min_j]})"
            dm.change_name_species_i(min_i, new_name)

            # Update distances for the new formed cluster using the NJ formula
            for x in range(dm.n_rows()):
                if x != min_i and x != min_j:
                    d_i = dm.get_value_i_j(min_i, x) # Distance from species min_i to x
                    d_j = dm.get_value_i_j(min_j, x) # Distance from species min_j to x
                    new_distance = (d_i + d_j - min_d) / 2 # Neighbor-Joining formula
                    dm.add_value_i_j(min_i, x, new_distance) # Update new distance
            # Remove the merged species from the matrix
            dm.remove_species_i(min_j)

        # Get the final two remaining species and construct the final tree with branch lengths
        species_names = dm.get_name_of_species()
        final_tree = f"({species_names[0]}:{dm.get_value_i_j(0,1)},{species_names[1]}:{dm.get_value_i_j(0,1)})"
        return final_tree # Return the final tree in Newick format
    
def main():
    sequence_ancestral = Sequence("Ancestral", "ACTGACTGACTGACTGACTG") 
    transition_probability = {
        "A": {"G": 0.03, "C": 0.03, "T": 0.03, "A": 0.91},
        "C": {"G": 0.03, "C": 0.91, "T": 0.03, "A": 0.03},
        "G": {"G": 0.91, "C": 0.03, "T": 0.03, "A": 0.03},
        "T": {"G": 0.03, "C": 0.03, "T": 0.91, "A": 0.03}
    }
    evolution = Evolution(sequence_ancestral, transition_probability)

    # Generate four sequences
    evolution.split_species_in_two("Ancestral", "SpeciesA")
    evolution.split_species_in_two("Ancestral", "SpeciesC")
    evolution.evolve(100)
    evolution.split_species_in_two("SpeciesA", "SpeciesB")
    evolution.evolve(50)
    evolution.split_species_in_two("SpeciesC", "SpeciesD")
    evolution.evolve(150)

    # Compute Pairwise Distances
    species_names = evolution.get_list_of_species_name()
    distance_matrix = DistanceMatrix(species_names) # Create a distance matrix using the species names
    for i in range(len(species_names)): # Loop through each pair of species
        for j in range(i + 1, len(species_names)): # Start from the next species to avoid repeating pairs
            # Calculate the nucleotide distance between the two species
            dist = ToolsToWorkWithSequences.observed_pairwise_nucleotide_distance(
                evolution.get_sequence_species(species_names[i]),
                evolution.get_sequence_species(species_names[j])
            )
            # Add the calculated distance to the distance matrix at the corresponding position
            distance_matrix.add_value_i_j(i, j, dist)

    # Run Algorithms
    distance_algorithms = DistanceBasedAlgorithms(distance_matrix)
    print(distance_algorithms.UPGMA())
    print(distance_algorithms.NJ())

    return 

if __name__ == '__main__':
    main()
    
