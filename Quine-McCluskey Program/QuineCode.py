# Quine McCluskey Code Implementation
# Fernandez, Joseph Bryan Lloyd C.

import itertools

# Function to generate and verify the Minterms
def DetermineMinterms (MaxVariables):
    # The while loop will not end if the number of minterms exceeds the max number of variables 
    while True:
        n = int(input('How many minterms would you like? '))
        if n > MaxVariables:
            print('Too many minterms')
        else: break
    # Maps the input of the user and creates an array of minterms    
    minterms = list(map(int, input('\nEnter the Numbers: ').strip().split()))

    if not len(minterms) <= n:  # Returns nothing if the number of minterms exceeds the maximum number
        print('Invalid number of minterms')
        return  
    else:
        MaxMinterm = max(minterms)
        if MaxMinterm > MaxVariables:
            print('Minterm/s exceeded the limit of 0 to {}'.format(MaxVariables))
            minterms = [] # CLears the list of minterms  

    print('\nThe Minterms are: {}'.format(minterms))
    return minterms            

# Function to convert decimal to binart and then group in order of number of 1s
def DecimaltoBCD(listofDecimal, NumberOfBits):
    listofBCD = [[] for x in range(NumberOfBits + 1)]  # Makes a list where the index signifies the number of 1s 

    for i in range(len(listofDecimal)):
        listofDecimal[i] = bin(listofDecimal[i])[2:]
        # Adds leading zeros to the binary if the nummber of bits is not equal to the number of literal variables
        if len(listofDecimal[i]) < NumberOfBits:
            for j in range(NumberOfBits - len(listofDecimal[i])):
                listofDecimal[i] = '0' + listofDecimal[i]
        
        # counts the number of 1s then appends them to their respected groups based on the number of 1s
        index = listofDecimal[i].count('1')
        listofBCD[index].append(listofDecimal[i])
    # Filters out empty groups in the list
    listofBCD = [group for group in listofBCD if group]
    return listofBCD

# Compares two binary, then true and the index if there is one difference between
def CompareBCD(binary1, binary2):
    # Initialize the difference counter and the index where the difference occurs
    difference = 0
    index = 0

    # Iterate over each bit in the binary representation
    for bit in range(len(binary1)):
        # Check if the bits at the same position are different
        if binary1[bit] != binary2[bit]:
            # If different, increment the difference counter and store the index
            difference += 1
            index = bit

    # If only one difference is found, return True and the index of the difference
    if difference == 1:
        return True, index
    else:
        # If there is no difference or more than one difference, return False and None
        return False, None

# Uses the CompareBCD function to compare binary groups, then append a '-' in the position of the difference
def CombineBinaryPairs(group, remaining):
    # Determine the length of the input group
    length = len(group) - 1
    # Initialize lists to store paired and unpaired binary strings in the next iteration
    pairedBinary = []
    nextGroup = [[] for x in range(length)]

    # Iterate through adjacent subgroups of the current group
    for i in range(length):
        for binary1 in group[i]:
            for binary2 in group[i+1]:
                # Check for a single differing bit in the binary pair
                bitDifference, position = CompareBCD(binary1, binary2)
                if bitDifference == True:
                    # Add the paired binary strings to the list
                    pairedBinary.append(binary1)
                    pairedBinary.append(binary2)

                    # Mark the differing bit with a dash ("-") to form a new binary string
                    newBinary = list(binary1)
                    newBinary[position] = '-'
                    newBinary = "".join(newBinary)

                    # Add the new binary string to the next group
                    nextGroup[i].append(newBinary)
    
    # Iterate through the current group to find unpaired binary strings
    for i in group:
        for j in i:
            if j not in pairedBinary:
                # Add unpaired binary strings to the remaining list
                remaining.append(j)
    
    # Return the next group and the remaining binary strings
    return nextGroup, remaining

# Removes duplicates
def RemoveRedundant(listofGroup):
    optimizedList = []
    for group in listofGroup:
        optimized = []
        for i in group:
            if i not in optimized:
                optimized.append(i)
        optimizedList.append(optimized)
    
    return optimizedList

# Check the current binary groups if there are still elements available for pairing
def CheckforPairings(group):
    if len(group) == 0 : 
        return True
    else:
        pairCount = 0
        for pair in group:
            if pair:
                pairCount += 1
        
        if pairCount == 0 : return True

    return False

# Converts a binary stream to its corresponding literal
def BinarytoLiterals(binary):
    # Initialize an empty string to store the resulting Boolean expression
    literal = ''
    # Start with the first letter 'A'
    letter = 'A'
    # Iterate through each bit in the binary representation
    for i in range(len(binary)):
        # If the bit is '1', add the current letter to the expression and increment the letter
        if binary[i] == '1':
            literal = literal + letter
            letter = chr(ord(letter) + 1)
        # If the bit is '0', add the negation of the current letter to the expression and increment the letter
        elif binary[i] == '0':
            literal = literal + letter + '\''
            letter = chr(ord(letter) + 1)
        # If the bit is '-', increment the letter without adding anything to the expression
        elif binary[i] == '-':
            letter = chr(ord(letter) + 1)
    
    return literal

# Matchers the implicants to its corresponding minterms
def checkBinary(binary, decimal):
    for i in range(len(binary)):
        if binary[i] != '-':
            if binary[i] != decimal[i]:
                return False
    return True

# Marks the minterms that is covered by the prime implicants
def Tabulation(unchecked, decimal):
    # Initialize a 2D table with dimensions determined by the lengths of unchecked and decimal lists
    table = [[' ' for column in range(len(decimal))] for row in range(len(unchecked))]
    # Iterate through each minterm and each prime implicant
    for i in range(len(decimal)):
        for j in range(len(unchecked)):
            # Check if the binary representation of the prime implicant covers the current minterm
            if checkBinary(unchecked[j], decimal[i]):
                # Mark the corresponding cell with 'X'
                table[j][i] = 'X'
    return table

# Creates a chart that shows the prime implicants and minterms
def FinalTabulation(unchecked, minterms, numofVariables, table):
    Letter = 'A'
    Literal = []
    minterms = list(sorted(minterms))
    for i in range(numofVariables):
        Literal.append(Letter)
        Letter = chr(ord(Letter) + 1)

    variables = ' '.join(Literal)
    minterms = list(map(str, minterms))

    # Find the maximum width for each column
    max_width_literals = max(len(literal) for literal in Literal)
    max_width_minterms = max(len(str(m)) for m in minterms)
    max_width_checkers = max(max(len(str(checker)) for checker in group) for group in table)

    # Set a fixed width for each column
    width_literals = max(max_width_literals, 2)  # Minimum width of 2
    width_minterms = max(max_width_minterms, 3)  # Minimum width of 2
    width_checkers = max(max_width_checkers, 3)  # Minimum width of 2

    # Print header
    header_literals = f'{variables.center(width_literals)} '
    header_minterms = ' '.join(str(term).center(width_minterms) for term in minterms)
    print(f' {header_literals}| {header_minterms}')

    # Print data
    for implicant, checker in zip(unchecked, table):
        implicants_str = ' '.join(implicant).center(width_literals)
        checker_str = ' '.join(str(value).center(width_checkers) for value in checker)
        print(f' {implicants_str} | {checker_str}')

# Saves all essential prime implicants in a list
def EssentialPrimeImplicants(tabulation):
    EssentialPrime = []

    # Iterate through each column of the tabulation
    for column in range(len(tabulation[0])):
        count = 0  # Count of 'X' entries in the current column
        pos = 0    # Row position of the 'X' entry in the current column
        # Iterate through each row of the tabulation
        for row in range(len(tabulation)):
            if tabulation[row][column] == 'X':
                count += 1
                pos = 0  # Reset the row position if 'X' is encountered
            # If there is exactly one 'X' entry in the column, append the row position to EssentialPrime
            if count == 1:
                EssentialPrime.append(pos)
    return EssentialPrime  


# Simulates the distribution of prime implicants and generates all possible combinations.
def MultiplicationOfList(list1, list2):
    # Handle the case where either list is empty
    result = []
    if len(list1) == 0 and len(list2) == 0:
        return result
    elif len(list1) == 0:
        return list2
    elif len(list2) == 0:
        return list1
    else:
        # Iterate through each element in list1 and list2
        for i in list1:
            for j in list2:
                # If elements are equal, add to the result list
                if i == j:
                    result.append(i)
                else:
                    # If elements are not equal, take their union and add to the result list
                    result.append(list(set(i+j)))
        # Sort the result list to ensure consistent ordering
        result.sort()
        # Remove duplicates from the result list using itertools.groupby
        return list(result for result, _ in itertools.groupby(result))

# Checks and Explores combinations of prime implicants
def PetrickMethod(tabulation):
    # Initialize a list to store prime implicants
    Petrick = []

    # Iterate through each column in the tabulation
    for column in range(len(tabulation[0])):
        prime = []
        # Iterate through each row in the column
        for row in range(len(tabulation)):
            # Check for 'X' marks in the tabulation, indicating a prime implicant
            if tabulation[row][column] == 'X':
                # Append the row (minterm) to the list of prime implicants for this column
                prime.append([row])
        # Append the list of prime implicants for this column to the Petrick list
        Petrick.append(prime)
    
    # Perform Petrick's Method multiplication of prime implicant lists
    for l in range(len(Petrick) - 1):
        Petrick[l + 1] = MultiplicationOfList(Petrick[l], Petrick[l+1])
    
    # Sort the final Petrick list based on the length of the implicant combinations
    Petrick = sorted(Petrick[len(Petrick) - 1], key=len)
    
    # Initialize a result list to store the minimal implicant combinations
    result = []
    # Set the minimum length to the length of the first implicant combination
    minimum = len(Petrick[0])
    # Iterate through implicant combinations in the Petrick list
    for terms in Petrick:
        # Check if the length of the current combination matches the minimum length
        if len(terms) == minimum:
            # Append the combination to the result list
            result.append(terms)
        else:
            # Break the loop if a longer combination is encountered
            break
    return result

# Counts the number of literals (non-dash characters) 
def NumberOfLiterals(term):
    count = 0
    for i in range(len(term)):
        if term[i] != '-':
            count += 1
    return count

# Checks if all minterms are covered and apply Petrick's method if not
def TraceMinimumCost(tabulation, uncheckedBinary):
    P_final = []
    essentialPrime = EssentialPrimeImplicants(tabulation)   
    OptimizedEssentialPrime = []

    for prime in essentialPrime:
        if prime not in OptimizedEssentialPrime:
            OptimizedEssentialPrime.append(prime)

    # Remove covered minterms by essential prime implicants
    if len(OptimizedEssentialPrime) > 0:
        for i in range(len(OptimizedEssentialPrime)):
            for column in range(len(tabulation[0])):
                if tabulation[OptimizedEssentialPrime[i]][column] == 'X':
                    for row in range(len(tabulation)):
                        tabulation[row][column] = 0
    # Check if all minterms are covered
    isZero = True
    for i in tabulation:
        for j in i:
            if j != 0:
                isZero = False
                break
            else:
                isZero = True
                P_final = [OptimizedEssentialPrime]
    # If not all minterms are covered, apply Petrick's Method
    if isZero == False:
        primes = PetrickMethod(tabulation)
        minCost = []
        
        # Calculate the cost of each prime implicant
        for prime in primes:
            count = 0
            for i in range(len(uncheckedBinary)):
                for j in prime:
                    if j == 1:
                        count += NumberOfLiterals(uncheckedBinary[i])
            minCost.append(count)
        # Select prime implicants with the minimum cost
        for i in range(len(minCost)):
            if minCost[i] == min(minCost):
                P_final.append(primes[i])
        # Add essential prime implicants to the final result
        for i in P_final:
            for j in OptimizedEssentialPrime:
                if j not in i:
                    i.append(j) 
    return P_final

def QuineMcClusky():
#-------------------Start-------------------------#
    LiteralVariables = 0
# Prompts the user to input an integer value
# The while loop will not end until the try-block is executed
    while LiteralVariables == 0:
        try:
            LiteralVariables = int(input('Please indicate the number of variables: '))
            break
        except ValueError:
            print('Invalid input. Please enter the number of variables')

    NumberofMinterms = pow(2, LiteralVariables) # Number of variables is 2^n
    DecimalMinterms = DetermineMinterms(NumberofMinterms)
    minterms = set(DecimalMinterms)

#-------------------Convert to BCD Format-------------------------#
    BinaryGroup = DecimaltoBCD(DecimalMinterms, LiteralVariables)
    print('\nThe BCD when group in terms of number of 1s are: {}'.format(BinaryGroup))

    uncheckedBinaryGroup = []

#-------------------Start of Comparing BCD-------------------------#
    while CheckforPairings(BinaryGroup) == False:
        nextBinaryGroup, uncheckedBinaryGroup = CombineBinaryPairs(BinaryGroup, uncheckedBinaryGroup)
        BinaryGroup = RemoveRedundant(nextBinaryGroup)
        if len(nextBinaryGroup) != 0:
            print('\nThe new group of binary pairs after comparing are: {}'.format(nextBinaryGroup))

    print('\nThe Prime implicants binary are: {}'.format(uncheckedBinaryGroup))

#-------------------Tabulate the Results-------------------------#
    print('\nPrime Implicants Tabulation')
    tabulation = Tabulation(uncheckedBinaryGroup, DecimalMinterms)
    FinalTabulation(uncheckedBinaryGroup, minterms, LiteralVariables, tabulation)

#-------------------Identify Minimum Expression-------------------------#
    EssentialPrimes = TraceMinimumCost(tabulation, uncheckedBinaryGroup)
    EssentialPrimes = RemoveRedundant(EssentialPrimes)
    ('\nThe Essential Prime Implicants are: {}'.format(EssentialPrimes))

    for primes in EssentialPrimes:
        answer = ''
        end = False
        for i in range(len(uncheckedBinaryGroup)):
            for j in primes:
                if j == i:
                    answer += BinarytoLiterals(uncheckedBinaryGroup[i]) + ' + '

    print('\nThe Boolean Expression for the given minterms is:\n{}'.format(answer[:-3]))

def Start():
    while True:
        User = input("\nStart the Program? (yes/no): ").lower()
        if User == 'yes':
            print("Great! Starting ...")
            QuineMcClusky()  # Call the main function to run the program again
        elif User == 'no':
            print("Thank you. Exiting...")
            break  # Exit the loop and end the program
        else:
            print("Invalid input. Please enter 'yes' or 'no'.")

Start()

