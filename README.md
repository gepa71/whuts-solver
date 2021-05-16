# whuts-solver
Code to find solutions for tiling the space by any of the 3d unfoldings of a 4d hypercube. See Matt Parker's video for details: https://www.youtube.com/watch?v=Yq3P-LhlcQo

### Usage Example
`echo '[[0,0,0],[0,0,-1],[0,-1,-1],[0,-1,-2],[-1,-1,-2],[1,0,0],[1,0,1],[1,1,1]]' | python3 whuts.py`

The input is in the format of the field `coords` in https://api.whuts.org/unfoldings

To get the tiling for unfolding `<N>` having that file can be done by:

`cat unfoldings | jq '.unfoldings[] | select(.id==<N>) .coords' | python3 whuts.py`

### Results

I run the code for all possible octocubes (octominoes in 3d?). The octocubes are generated by the file `generate_octocubes.py`. The program could find solutions for all but 3 octocubes! From those 3 I could find a tiling for 2 of those using a similar code in C that runs faster. The solutions are in the directory `octocube_solutions`. These cover of course also all hypercube unfoldings, which is a subset of all octocubes.

This is the octocube for which the code could not find a tiling (yet?):

```
    ###
    #.#
    ###
```

Due to the very restrictive nature of this one (two of these have to be pairwise connected, so it reduces to finding a tiling of space using the 16-cube containing two of these instances connected together), I tried writing more optimized C code, but still didn't get any results. I will update here if something is found.
