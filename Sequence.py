'''

@author: Luis Carlos Ospina Restrepo
'''


class Sequence(object):
    '''
    Defining a class that represents a sequence of nucleotides, with methods to mutate the sequence and track mutations.
    '''


    def __init__(self, name, string_of_nucleotides,mutations=None):
        '''
        question: Why are we using list(str)?
        Answer: Because str are inmutable, and we will probably change a nucleotide. Therefore we change it to a list.
        question: What is this mutations=None doing?
        It is just a default value if it is not provided as input
        Constructor:

        Initializing the Sequence object with a name, nucleotide string or list, and an optional mutations list.
        
        Arguments:
        name-> the name of the sequence (e.g., species or sample name)
        string_of_nucleotides-> a string or list representing the sequence of nucleotides
        mutations-> a list of mutation counts (optional, defaults to a list of zeros)
        '''
        self.name = name
        
        # If the input for string_of_nucleotides is a string, convert it into a list of characters (nucleotides)
        if isinstance(string_of_nucleotides, str):
            self.string_of_nucleotides = list(string_of_nucleotides)
        elif isinstance(string_of_nucleotides,list):#if it is already a list, directly assign it
            self.string_of_nucleotides = string_of_nucleotides
        else:# If the input is neither a string nor a list, raise an error
            raise TypeError("Use a string to initialize the sequence or a list of nucleotides")
            
            
        if mutations== None:# If no mutations are provided, initialize the mutations list as a list of zeros (one for each nucleotide)
            self.mutations = [0]*len(string_of_nucleotides)
        else:
            self.mutations = mutations# Otherwise, use the provided mutations list
    
    # Method to get the nucleotide at a specific position (index) in the sequence
    def nucleotide_at_position(self,position):
        '''
        Returns the nucleotide at a specified position in the sequence.
        '''
        return self.string_of_nucleotides[position]
 
     # Method to return the length of the sequence (number of nucleotides)
    def sequence_length(self):
        '''
        Returns the length of the nucleotide sequence.
        '''
        return len(self.mutations)
    
    
    def get_name(self):
        '''
        Returns the name of the sequence.
        '''
        return self.name
   
    
    def copy(self, new_name):
        '''
        Returns a new Sequence object which is a copy of the current sequence, but with a new name.
        '''
        return Sequence(new_name,self.string_of_nucleotides.copy(),self.mutations.copy())
    
    
    def __str__(self):
        '''
        Returns a string representing the sequence object.        
        '''
        return self.name + ' : ' + str(self.string_of_nucleotides)
    
    
    def mutate_nucleotide_at_position(self,position,nucleotide):
        '''
        Mutates the nucleotide at a specified position by replacing it with a new nucleotide.
        If the nucleotide is changed, the mutation count at that position is incremented.
        '''
        if(nucleotide!=self.string_of_nucleotides[position]):
            self.mutations[position] = self.mutations[position] + 1
        
        self.string_of_nucleotides[position] = nucleotide
        
    def get_number_of_mutations_at_each_position(self):
        '''
        Returns the mutation counts for each nucleotide in the sequence.
        '''
        return self.mutations


    
def main():
    '''
    The main function where two Sequence objects are created, mutated, and displayed.
    '''
    sequence_a = Sequence("First_Species", list("ACTGACTG")) # Creating the first sequence object with a name and a nucleotide sequence
    sequence_b = sequence_a.copy("Second_Species") # It just creating a copy
    
    sequence_a.mutate_nucleotide_at_position(1, "C")# Mutating a nucleotide in the first sequence at position 1

    sequence_b.mutate_nucleotide_at_position(2, "A")# Mutating a nucleotide in the second sequence at position 2
    
    # It just printing the results
    print(sequence_a)
    print(sequence_b)
    print("")
    
    print(sequence_a.get_number_of_mutations_at_each_position())
    print(sequence_b.get_number_of_mutations_at_each_position())
    
    
if __name__ == "__main__":
    main () 
