mendelian-error
===============

This package attempts to assign probability to a mendelian error event in a trio.
For example, given a mother with genotype 'C/C' and a father with genotype 'C/C'
a child with genotype 'C/T' will be a "mendelian error", in this case, a candidate
*de novo* mutation.

We can filter candidates based on their genotype likelihoods. For example if the
child had a genotype likelihood of -2,-1,-20 then we are likely to consider this
a genotyping error because the homozygous reference ("C/C" with GL -2) is fairly
close the the GL for het (with GL -1). If the genotype likelihood field was
-20,0,-20, then the call is confidently het. We can use the genotype likelihoods
to assign a probability:

```Python

>>> from mendelianerror import
# everyone is homref. this should have a low probability of an error:
>>> father = mother = child = [-0.1, -8.0, -8.0]
>>> mendelian_error(mother, father, child)
7.55...e-08


# parents are hom, child is het. this is a likely mendelian error:
>>> father = mother = [-0.6, -2.5, -2.5]
>>> child = [-2.5, -0.6, -2.5]
>>> mendelian_error(mother, father, child)
0.987...

```

So the input is the 3 GL numbers for each of the father, mother, child.

CLI
===

After installation, one can use this on a multi-sample VCF file like:
```Shell
mendelianerror $input.vcf father_id mother_id child_id > $new.vcf
```

Limitations
===========

+ Only make sense for autosomal variants.
+ Only works on trios (doesn't consider extended pedigrees or siblings).
