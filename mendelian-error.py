
def rescale(li):
    s = float(sum(li))
    return [v / s for v in li]

def mendelian_error(mother, father, child):
    """
    Return the probability of a mendelian error given the log10 genotype
    likelihoods. A large value indicates a high probability of a mendelian
    error. Low values mean that the genotype-likelihoods indicate enough
    uncertainty that it could be a genotyping error.

    # everyone is het:
    >>> het = (-2.0, -0.1, -2.0)
    >>> mendelian_error(het, het, het)

    # parents are hom, child is het.
    >>> father = mother = [-0.6, -2.5, -2.5]
    >>> child = [-2.5, -0.6, -2.5]
    >>> mendelian_error(mother, father, child)
    0.939...

    # same as above, but more certainty in the called genotypes:
    >>> child[1] = 1e-6
    >>> mother[0] = father[0] = 1e-6
    >>> mendelian_error(mother, father, child)
    0.984...

    # everyone is confidently homozygous alt
    >>> child = father = mother = [-11.0, -11.0, -0.1]
    >>> mendelian_error(mother, father, child)
    5.03...e-11

    # everyone is less confidently homozygous refs:
    >>> child = father = mother = [-0.1, -2.0, -2.0]
    >>> mendelian_error(mother, father, child)
    0.047...

    mother and fater are homozygous alts
    >>> mother = father = [-3.0, -3.0, -0.1]

    # child is het
    >>> child = [-3., -0.1, -3.]
    >>> mendelian_error(mother, father, child)
    0.993...

    # but when the hom-alt call is close...
    >>> child = [-3., -0.1, -0.15]
    >>> mendelian_error(mother, father, child)
    0.52...

    # mother is hom_ref, dad is het, child is hom_alt
    >>> mother, father, child = (-0.1, -2, -2), (-2, -0.1, -2), (-2, -2, -0.1)
    >>> mendelian_error(mother, father, child)
    0.940...

    # mother is hom_ref, dad is hom_alt, child is hom_ref
    >>> mother, father, child = (-0.1, -2.5, -2.5), (-2.5, -2.5, -0.1), (-0.1, -2.5, -2.5)
    >>> mendelian_error(mother, father, child)
    0.984...

    # same, but child is hom_alt
    >>> mendelian_error(mother, father, (-5, -5, -0.01))
    0.988...

    # child should be het:
    >>> mendelian_error(mother, father, (-2, -0.1, -2))
    0.024...


    """
    M = rescale([10.**m for m in mother])
    F = rescale([10.**f for f in father])
    C = rescale([10.**c for c in child])
    assert sum(C) - 1.0 < 0.001
    assert sum(M) - 1.0 < 0.001
    assert sum(F) - 1.0 < 0.001

    # by ref, and alt, we mean hom_ref, hom_alt
    p_two_ref = M[0] * F[0]
    p_two_het = M[1] * F[1]
    p_two_alt = M[2] * F[2]

    # only 1 of the parents is ...
    p_one_ref = M[0] + F[0] - p_two_ref
    p_one_het = M[1] + F[1] - p_two_het
    p_one_alt = M[2] + F[2] - p_two_alt

    #################################
    # Ways to *not* have violations #
    #################################
    # 1. everyone is reference
    a = p_two_ref * C[0]
    # 2. everyone is hom alt
    b = p_two_alt * C[2]
    # 3. 1 het and 1 ref parent. child matches
    c = p_one_het * p_one_ref * (C[0] + C[1])
    # 4. 1 het and 1 alt parent. child matches
    d = p_one_het * p_one_alt * (C[1] + C[2])
    # 5. both parents hets. (child can be anything)
    e = p_two_het
    # 6. one hom ref, one home alt. child is het
    f = p_one_ref * p_one_alt * C[1]
    print a, b, c, d, e, f

    p_not_error = a + b + c + d + e + f
    return 1.0 - p_not_error

if __name__ == "__main__":
    import doctest
    import sys
    sys.stderr.write(str(doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS | doctest.REPORT_ONLY_FIRST_FAILURE, verbose=0)) + "\n")