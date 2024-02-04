
import argparse

from read_alignment import read_alignment
import profile_hmm as phmm

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Takes in an MSA, a pseudocount, and a list of sequences and returns a list of the sequences aligned to a profile HMM built on the MSA with the pseudocount added to each count.')
    parser.add_argument('msa', type=str,
                        help='Path to the MSA')
    parser.add_argument('pseudocount', type=str,
                        help='Pseudocount to be added to each count')
    parser.add_argument('sequences', type=str,
                        help='Sequences to be aligned to the profile HMM')
    args = parser.parse_args()
    msa = read_alignment(
        args.msa, args.msa[args.msa.find('.') + 1:])
    sequences = read_alignment(
        args.sequences, args.sequences[args.sequences.find('.') + 1:])

    profile_hmm = phmm.ProfileHMM(msa, float(args.pseudocount))

    print(
        f"""
=-=-=-=-=-=-=-=-=-=-=-=-=-=
Profile HMM
=-=-=-=-=-=-=-=-=-=-=-=-=-=

---------------------------
Pseudocount: {args.pseudocount}
---------------------------
Number of states: {profile_hmm.num_match_states * 3 + 3}
---------------------------
Transitions matrix:
""")
    [print(f"{row}") for row in profile_hmm.transition_matrix]
    print("---------------------------")
    print("Emissions matrix (match states):\n")
    [print(f"{row}") for row in profile_hmm.m_emission_matrix]
    print("---------------------------")

    print(
        """
=-=-=-=-=-=-=-=-=-=-=-=-=-=
Sequence Alignment
=-=-=-=-=-=-=-=-=-=-=-=-=-=

---------------------------
Unaligned sequences:
"""
    )

    aligned_sequences = profile_hmm.align_sequences(sequences)

    [print(f"{s}") for s in sequences]
    print("---------------------------")
    print("Aligned sequences:\n")
    [print(f"{s}") for s in aligned_sequences]
    print("---------------------------")
