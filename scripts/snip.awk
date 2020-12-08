# Generates .snip files from g09 .out files.
#
# usage: awk -f snip.awk *.out
#
# Takes all *.out files and turns them into *.snip files.
# (Will overwrite any existing .snips)
#
# These files can be used with PyQuiver, but not with the
# old version of Quiver.
BEGIN {
    if ( ARGC < 2 ) {
        print "Error: please specify some filenames."
        print ""
        print "snip.awk: generates .snip files from Gaussian .out files"
        print "These files will be compatible with PyQuiver."
        print ""
        print "usage: awk -f snip.awk *.out"
        exit
    }
}

FNR == 1 {
    fileCount++
    filenames[fileCount]=FILENAME
}

/Standard orientation/,/Rotational constants/ {
    if ( NF == 6 && match($0,"[a-zA-Z]") == 0 )
        {
            symbol[fileCount,$1]=$2
            x[fileCount,$1]=$4
            y[fileCount,$1]=$5
            z[fileCount,$1]=$6
            atoms[fileCount]=$1
        }
}

/SCF Done/ {
    energy[fileCount]=$5
}

/NAtoms=/ {
    split($0, natomsplit, "NAtoms=")
    split(natomsplit[2], natomsplitsplit, " ")
    numberOfAtoms[fileCount]=natomsplitsplit[1]
}

/Zero-point correction/ {
    energy_field_count = 0
}

/Zero-point correction/,/Sum of electronic and thermal Free Energies/ {
    split($0, line, "=")
    #print trim(line[2])
    energy_field[fileCount, energy_field_count] = trim(line[2])
    energy_field_count ++
}

/l9999.exe/ {
    archive[fileCount] = ""
}

/l9999.exe/,/\\\@/ {
    #archive[fileCount] = archive[fileCount] trim($0)
    s = $0
    gsub(/%/, "%%", s)
    #n = split(s, terminations, "l9999")
    #archiveString = "l9999" terminations[n]
    archive[fileCount] = archive[fileCount] s
}

END {
    if ( exitCode != 0 )
        exit

    split("Zero-point correction=@ Thermal correction to Energy=@ Thermal correction to Enthalpy=@ Thermal correction to Gibbs Free Energy=@ Sum of electronic and zero-point Energies=@ Sum of electronic and thermal Energies=@ Sum of electronic and thermal Enthalpies=@ Sum of electronic and thermal Free Energies=", energy_lines, "@")
    for (i=1; i <= fileCount; i++)
        {
            filename = filenames[i]
            gsub(".out",".snip",filename)
            print filename

	    print "" > filename

	    printf "Electronic energy: %15.8f\n", energy[i] >> filename

	    for (j=0; j<=8; j++)
		printf "%s %s\n", energy_lines[j+1], energy_field[i,j] >> filename

            print "Standard orientation" >> filename
            for (j=1; j <= atoms[i]; j++)
                printf "%d %2d 0 %15.8f %15.8f %15.8f\n", j, symbol[i,j], x[i,j], y[i,j], z[i,j]  >> filename
	    print "Rotational constants (GHZ)" >> filename
	    printf "NAtoms= %d \n", numberOfAtoms[i] >> filename
	    printf "SCF Done Field3 Field4 %15.8f \n", energy[i] >> filename
	    printf archive[i] >> filename
            printf "\n\n" >> filename
        }
}

# functions to trim whitespace
function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }
function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }
function trim(s) { return rtrim(ltrim(s)); }
