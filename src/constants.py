# This class holds physical constants and atomic weights.
import sys
import re

###############

# Physical Constants

h = 6.625E-27
#CM/2.998E10/,EM/1.440E13/,HBC/1.4387/

###############

# Atomic Weight Information

class Element(object):
    def __init__(self, full_name, atomic_number, symbol, default_mass):
        # the name of this element, like "hydrogen"
        full_name = str(full_name)
        self.full_name = full_name
        if re.match("[^a-z]", full_name):
            print "Unexpected non-lowercase character in element name: %s" % full_name
            print "Quitting."
            sys.exit(1)

        # the symbol of this element, like "H"
        symbol = str(symbol)
        self.symbol = symbol
        if re.match("[^a-zA-Z]", symbol):
            print "Unexpected non-letter character in element symbol: %s" % symbol
            print "Quitting."
            sys.exit(1)
        if len(symbol) < 1 or len(symbol) > 2:
            print "Unexpected length of element symbol (must be 1 or 2): %s" % symbol
            print "Quitting."
            sys.exit(1)

        # the atomic number of this element, like 1
        atomic_number = int(atomic_number)
        self.atomic_number = atomic_number
        if atomic_number < 1 or atomic_number > 200:
            print "Unexpected atomic number: %d" % atomic_number
            print "Quitting."
            sys.exit(1)
        
        # the average weight for this element, like 1.00783
        default_mass = float(default_mass)
        self.default_mass = default_mass
        if default_mass < 0.0 or default_mass > 500.0:
            print "Unexpected default mass: %d" % default_mass
            print "Quitting."
            sys.exit(1)
 
        # pairs of tuples strings (like "2H") to masses (like 2.0141)
        self.replacements = []

    def __str__(self):
        string = "%s (%s, Z=%d, default mass = %.4f" % (self.full_name.capitalize(), self.symbol, self.atomic_number, self.default_mass)
        if len(self.replacements) == 0:
            string += ", no isotopic replacements possible)\n"
        else:
            string += ")\n"
            for s,m in self.replacements:
                string += "    %2s : %.4f\n" % (s,m)
        return string[:-1]

    def add_replacement(self, symbol, mass):
        symbol = str(symbol)
        if re.match("[^a-zA-Z0-9]", symbol):
            print "Unexpected non-letter character in isotopic replacement symbol: %s" % symbol
            print "Quitting."
            sys.exit(1)
        if len(symbol) < 1 or len(symbol) > 4:
            print "Unexpected length of element symbol in replacement (must be 1-4 inclusive, found %d): %s" % (len(symbol), symbol)
            print "Quitting."
            sys.exit(1)
        for s,m in self.replacements:
            if s == symbol:
                print "Must use a unique symbol for every isotopic replacement: %s" % s
                sys.exit(1)
        mass = float(mass)
        if mass < 0.0 or mass > 500.0:
            print "Unexpected isotopic replacement mass: %f" % mass
            sys.exit(1)
        self.replacements.append((symbol,mass))

# read in atomic weight data
elements = []
for line in open("weights.dat", "r"):
    # ignore comments and blank lines
    line = line.strip()
    if line[0] == "#" or len(line) == 0:
        continue
    line = line.split("#",1)[0]

    # parse
    fields = line.encode("ascii","ignore").split(",")
    if len(fields) < 4:
        print "Error: not enough data on this line of weights.dat:"
        print line
        print "\nQuitting."
        sys.exit(1)
    element = Element(*fields[0:4])
    if (len(fields)-4) % 2 != 0:
        print "Unexpected number of isotopic replacement fields on this line of weights.dat."
        print "The number of fields after the first four must be a multiple of 2 (found %d)." % (len(fields)-4)
        print line
        print "\nQuitting."
        sys.exit(1)
    if (len(fields) > 4):
        for i in range(4, len(fields), 2):
            element.add_replacement(fields[i], fields[i+1])
    elements.append(element)
    #print element
print "Read atomic weight data for %d elements." % len(elements)

