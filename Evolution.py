'''

@author: Luis Carlos Ospina Restrepo
'''
from Sequence import Sequence
from Alias_Vose import RandomMultinomial

class Evolution(object):
    '''
    This defines a class that simulates DNA sequence evolution based on a given transition matrix.    
    '''
    
    def __init__(self, ancestral_sequence, transition_matrix):


        '''
        This is the constructor, it initializes the Evolution class.
        
        Parameters:
        ancestral_sequence (Sequence): The initial sequence from which evolution starts.
        transition_matrix (dict): A dictionary that defines mutation probabilities between nucleotides.
        '''
        if not isinstance(ancestral_sequence, Sequence): # Here we are checking if the ancestra_sequence is a Sequence object
            raise TypeError("ancestral_sequence must be a Sequence object!") # If not then an error is raised
        
        # Now we are creating attributes 
        self.evolving_sequences = {} # This one is a dictionary that stores species and their corresponding sequences
        self.evolving_sequences[ancestral_sequence.get_name()] = ancestral_sequence # Adding ancestral sequence
        self.transition_matrix = transition_matrix # This is an attribute to store the given mutation probability matrix
        self.random_transition = {}# This is an attribute (dict) to hold random mutation generators for each nucleotide type

        for key in self.transition_matrix: # This converts the mutation probability distributions into random samplers
            self.random_transition[key] = RandomMultinomial(list(transition_matrix[key].values()))
        
    
    def get_sequence_species(self,name):
        '''
        This is used to retrieve the sequence object of a given species.
        
        Parameters:
        name (str): Name of the species.
        
        Returns:
        Sequence: The sequence object of the given species.
        '''
        return self.evolving_sequences[name]
    
    
    def get_list_of_species_name(self):
        '''
        This retrieve a list of all species names in the evolving population.
        
        Returns:
        list: List of species names.
        '''
        return list(self.evolving_sequences.keys())
    
    
    def split_species_in_two(self,name_of_species_that_splits, new_name_of_species):
        '''
        This creates a new species from an existing species by duplicating its sequence.
        
        Parameters:
        name_of_species_that_splits (str): Name of the original species.
        new_name_of_species (str): Name for the new species.
        '''
        new_sequence = self.evolving_sequences[name_of_species_that_splits].copy(new_name_of_species)
        self.evolving_sequences[new_sequence.get_name()] = new_sequence
    
    
    def evolve(self, generations):
        '''
        This evolves all species for a given number of generations.
        
        Parameters:
        generations (int): Number of generations to evolve.
        '''
        for gen in range(generations): # Here we iterate through each generation

            for species in self.evolving_sequences:# In this line we iterate through each species

                sequence_species = self.evolving_sequences[species]# We get the sequence object
                for n in range(sequence_species.sequence_length()):# It iterates through nucleotide positions

                    nucleotide_at_position_n = sequence_species.nucleotide_at_position(n)# This gets the nucleotide at position n
                    propose_change = self.random_transition[nucleotide_at_position_n].sample() # Now this line sample a proposed change based on mutation probabilities
                    nucleotide_propose_change = list(self.transition_matrix[nucleotide_at_position_n].keys())[propose_change] # Finally it determines the new nucleotide after mutation

                    print(nucleotide_at_position_n,nucleotide_propose_change)# It just prints what we got
                    sequence_species.mutate_nucleotide_at_position(n,nucleotide_propose_change)# This apply the mutation to the sequence

            
class ToolsToWorkWithSequences():
    '''
    This class provides different methods, such as calculating nucleotide statistics.
    '''
    
    def nucleotide_statistics(sequence):
        '''
        This method calculates the percentage of each nucleotide (A, C, T, G) in a given sequence.
        
        Parameters:
        sequence (Sequence): A Sequence object whose nucleotide statistics are to be calculated.
        
        Returns:
        dict: A dictionary with the percentage of each nucleotide.
        '''
        # Initialize counters for each nucleotide
        count_A = sequence.string_of_nucleotides.count('A')
        count_C = sequence.string_of_nucleotides.count('C')
        count_T = sequence.string_of_nucleotides.count('T')
        count_G = sequence.string_of_nucleotides.count('G')
        
        # Calculate the total length of the sequence
        total_nucleotides = len(sequence.string_of_nucleotides)
        
        # Calculate the percentage of each nucleotide
        stats = {
            'A': (count_A / total_nucleotides) * 100,
            'C': (count_C / total_nucleotides) * 100,
            'T': (count_T / total_nucleotides) * 100,
            'G': (count_G / total_nucleotides) * 100
        }
        
        return stats
    def observed_pairwise_nucleotide_distance(sequence1, sequence2):
        '''
        This method calculates the observed pairwise nucleotide distance between two sequences.
        It counts the number of positions where the nucleotides are different in the two sequences.
        
        Parameters:
        sequence1 (Sequence): The first sequence to compare.
        sequence2 (Sequence): The second sequence to compare.
        
        Returns:
        int: The number of positions where the nucleotides are different in the two sequences.
        '''
        # Ensure both sequences have the same length
        if len(sequence1.string_of_nucleotides) != len(sequence2.string_of_nucleotides):
            raise ValueError("Sequences must have the same length to compare.")
        
        # Initialize the counter for differences
        distance = 0
        
        # Iterate over the positions in the sequences and compare the nucleotides
        for i in range(len(sequence1.string_of_nucleotides)):
            if sequence1.string_of_nucleotides[i] != sequence2.string_of_nucleotides[i]:
                distance += 1
        
        return distance


def main():
    # Firstly we initialize ancestral sequence
    sequence_ancestral = Sequence("Ancestral", "ACTGACTGACTGACTGACTGACTGACTGACTGACTG")
    transition_probability = {
        "A": {"G": 0.03, "C": 0.03, "T": 0.03, "A": 0.91},
        "C": {"G": 0.03, "C": 0.91, "T": 0.03, "A": 0.03},
        "G": {"G": 0.91, "C": 0.03, "T": 0.03, "A": 0.03},
        "T": {"G": 0.03, "C": 0.03, "T": 0.91, "A": 0.03}
    }
    evolution = Evolution(sequence_ancestral, transition_probability)
    
    # Now we evolve for 10 generations
    evolution.evolve(10)
    
    # After evolving we split the species into SpeciesA and SpeciesB
    evolution.split_species_in_two("Ancestral", "SpeciesA")
    evolution.split_species_in_two("Ancestral", "SpeciesB")

    # We keep evolving
    evolution.evolve(590)
    
    # Now we create species C and D from ancestral
    evolution.split_species_in_two("Ancestral", "SpeciesC")
    evolution.split_species_in_two("Ancestral", "SpeciesD")


    # We keep evolving
    evolution.evolve(8400)

    # Now we just print the result
    for species in evolution.get_list_of_species_name():
        sequence = evolution.get_sequence_species(species)
        print(f"Species: {species}")
        print(sequence)
        print("Nucleotide Statistics:", ToolsToWorkWithSequences.nucleotide_statistics(sequence))  # Call the method using the instance
        print("-" * 40)
    # Comparing two sequences (for example, 'SpeciesA' and 'SpeciesB') and calculating their observed pairwise nucleotide distance
    species_A = evolution.get_sequence_species("SpeciesA")
    species_B = evolution.get_sequence_species("SpeciesB")
    
    distance = ToolsToWorkWithSequences.observed_pairwise_nucleotide_distance(species_A, species_B)
    print(f"Observed pairwise nucleotide distance between SpeciesA and SpeciesB: {distance}")
    


if __name__ == "__main__":
    main()
