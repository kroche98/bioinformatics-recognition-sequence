import sys
import re

enz_to_recseq_re = {}

def parse_rebase():
    '''Parse the file 'link_allenz.txt' into a dictionary of compiled regular expressions
    for recognition sequences keyed by enzyme names.
    The file must be in the same directory as the script.
    '''

    global enz_to_recseq_re
    
    # regexp to capture enzyme name and recognition sequence from REBASE format
    pattern = r'<1>(.*)\n(?:.*\n){3}<5>(.*)\n'
    regexp = re.compile(pattern)
    
    with open("link_allenz.txt", 'r') as file:
        hits = regexp.findall(file.read())
    
    # use the hits to create a dict of compiled regexps for recognition sequences
    # keyed by enzyme name, with enzymes with unknown recognition sequences removed
    enz_to_recseq_re = {enz[0]:re.compile(make_re(enz[1]), flags=re.I) for enz in hits if enz[1] != '?'}

def find_matches(enz, strand):
    '''Take an enzyme name and a DNA strand and returns a list of
    all positions (starting at 1) where the enzyme's recognition sequence occurs.
    Raises a ValueError exception if an invalid enzyme is given.
    '''
    
    if enz not in enz_to_recseq_re.keys():
        raise ValueError("Enzyme without known recognition sequence")
    else:
        recseq = enz_to_recseq_re[enz]
        
    positions = []
    next_match = recseq.search(strand)
    
    while next_match != None:
        i = next_match.start()
        positions.append(i+1) # since we want to start numbering at 1, not 0
        next_match = recseq.search(strand, i+1)

    return positions

def strand_to_binding_enzymes(strand):
    '''Take a strand of DNA and return a list of all enzymes that bind
    to at least one recognition site, sorted alphabetically
    '''
    binding_enzymes = [enz for enz in enz_to_recseq_re.keys() if find_matches(enz, strand)]
    binding_enzymes.sort()
    return binding_enzymes

def make_re(rec_seq):
    '''Take a recognition sequence in REBASE form and return a string with
    an equivalent regexp
    '''
    base_to_re = {'A':'A', 'T':'T', 'G':'G', 'C':'C', 'R':'[AG]', 'Y':'[CT]', 'M':'[AC]', 'K':'[GT]', 'S':'[GC]', 'W':'[AT]', 'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]', 'N':'[ACGT]', ',':'|'}
    return ''.join(base_to_re.get(char, '') for char in rec_seq)

def main(strand):
    parse_rebase()
    print('\n'.join(strand_to_binding_enzymes(strand)))

if __name__=='__main__':
    if len(sys.argv) == 1:
        main(input("Enter a DNA sequence: "))
    else:
        main(sys.argv[1])
