import sys

CIRCUIT = '.circuit'
END = '.end'

if len(sys.argv)!=2:
    print('Error! Correct format: python  <input file> <netlist file>')
    exit()


def analyze(lines): # A function to analyze and extract information from words(tokens)
    analysis = dict() # Initialize an empty dictionary

    names = {'R':'Resistor', 'L':'Inductor', 'C':'Capacitor', 'V':'Independent Voltage source', 'I':'Independent current source', 'E':'VCVS', 'G':'VCCS', 'H':'CCVS', 'F':'CCCS'}

    for i in range(1,len(lines)+1):
       lines[i-1] = lines[i-1].split('#')[0].strip() # Discard the comment and remove any whitespace
       line = lines[i-1].split(' ') # Convert the string to list
    
       analysis['line'+ str(i)] = {} # Initialize a nest in the dictionary
       analysis['line'+ str(i)]['startnode'] = line[1] # Extract information from the list
       analysis['line'+ str(i)]['endnode'] = line[2]
       analysis['line'+ str(i)]['value'] = line[-1]
    
       for name in names.keys():
          if name == line[0][0]:
            analysis['line'+ str(i)]['typename'] = names[name] # Identify the type of element
    
       if len(line)==5:
          analysis['line'+str(i)]['Voltage_id'] = line[3]
        
       if len(line)==6:
          analysis['line'+str(i)]['dependent_startnode'] = line[3]
          analysis['line'+str(i)]['dependent_endnode'] = line[4]
           
    return analysis # return the nested dictionary




try: 
    with open(sys.argv[1], 'r') as f: # open the file in read mode
        lines = f.readlines() # Store each line of the file in a list
        start = -1; end = -2
        for line in lines:              # extracting circuit definition start and end lines
            if CIRCUIT == line[:len(CIRCUIT)]:
                start = lines.index(line) # Store the start of circuit block
            elif END == line[:len(END)]:
                end = lines.index(line) # Store the end of circuit block
                break
    
    
        analysis = analyze(lines[start+1:end]) # Analyze the tokens and save them in a nested dictionary 'analysis'
            
        if (start + 1 < end): # Checking validity of circuit block

            for line in reversed(lines[start+1:end]):
                if line.find('#') != -1:       # Check for comments and remove if present
                    line = line.split('#')[0]

                if (line != '' and line != '\n'):  # Check for the emptiness of line, if empty, skip the line
                    line = line.split()
                    line = " ".join(line[::-1]) # Join list elements in reverse into a string with space in between

                    print(line) #prints the output

        else:
            print("Error: Invalid circuit definition")
            exit()


except IOError:
    print('Invalid file! Enter a valid file name') # Error for invalid or non-existent file
    exit()



