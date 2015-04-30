'''
function to return six frame translation-to-protein of a DNA sequence

Bio.SeqUtils.six_frame_translations seems to be broken in 1.59
and removed from 
'''

from Bio.Seq import reverse_complement, translate

def six_frames(seq, genetic_code=1):
    '''
    input a DNA sequence
    return the six possible protein translations
    
    note: frames 1,2,3 are with respect to the sequence start
    while frames -1,-2,-3 are with respect to the sequence end
    therefore 1 and -1 are not reverse complements unless the sequence
    length is a multiple of 3
    '''
    
    rev = reverse_complement(seq)
    
    frames = {}
    
    for i in [1,2,3]:
        frames[i] = translate(seq[i-1:],genetic_code)
        frames[-i] = translate(rev[i-1:],genetic_code)
        
    return frames
