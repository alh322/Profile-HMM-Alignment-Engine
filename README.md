# Profile HMM Alignment Engine
![image](https://github.com/alh322/Profile-HMM-Alignment-Engine/assets/117471005/44ccac32-4fa8-4527-b35e-33e96cac0b64)




## Description

This script is designed to align a set of sequences to a profile hidden Markov model (HMM) constructed from a given multiple sequence alignment (MSA) and a specified pseudocount. It takes an MSA file, a pseudocount value, and a file containing sequences as input. The output includes the details of the profile HMM and the aligned sequences.

## Requirements

Python environment.
Input files: MSA file and sequences file in a fasta format.

## Usage

Run the script from the command line with the following syntax:

```
python main.py <msa_path> <pseudocount> <sequences_path>
```

## Arguments

`msa_path`: Path to the MSA file.
`pseudocount`: The pseudocount value to be added to each count in the profile HMM.
`sequences_path`: Path to the file containing sequences to be aligned.

## Example

```
python main.py test_msa.fasta 1 test_seq.fasta
```

## Output

The script outputs:

Details of the Profile HMM:

- Pseudocount used.
- Number of states in the HMM.
- Transition matrix of the HMM.
- Emissions matrix for match states in the HMM.

Sequence Alignment:

- List of unaligned sequences.
- List of sequences aligned to the profile HMM.

## Note

Ensure that the input files are in a format compatible with the read_alignment function used in the script. The format is inferred from the file extension (we only support fasta files for now but this is easily extensible). The ProfileHMM class is used to create and use the profile HMM, and it requires the MSA and pseudocount as inputs.
