"""Exact preservation of literal VCF alleles under plain-merge right padding."""


_LITERAL_BASES = frozenset("ACGTRYSWKMBDHVN")


def exact_padding_equivalence(start_a, ref_a, alt_a, start_b, ref_b, alt_b):
    if (not isinstance(start_a, int) or isinstance(start_a, bool)
            or not isinstance(start_b, int) or isinstance(start_b, bool)
            or start_a < 0 or start_b < 0 or start_a != start_b):
        return False
    alleles = (ref_a, alt_a, ref_b, alt_b)
    if any(not isinstance(allele, str) or not allele
           or set(allele.upper()) - _LITERAL_BASES for allele in alleles):
        return False
    ref_a, alt_a, ref_b, alt_b = (allele.upper() for allele in alleles)
    if len(ref_a) <= len(ref_b):
        return ref_b.startswith(ref_a) and alt_a + ref_b[len(ref_a):] == alt_b
    return ref_a.startswith(ref_b) and alt_b + ref_a[len(ref_b):] == alt_a
