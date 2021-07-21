import argparse

def align_sequences():
    # - align 16S sequences to reference
    return

def profile_functions():
    # needs bug abundances per sample, K and 16S abundances per bug
    #   replace sequences with functional profiles using formula
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create functional profile\
                                     from 16S data")
    #add arguments for abundance table, representative sequences,
    # additional parameters: % identity

    align_sequences()

    for samp in samples:
        profile_functions(samp)

    #append results to new table, save



