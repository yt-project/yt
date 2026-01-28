import pytest

from yt.funcs import levenshtein_distance


@pytest.mark.parametrize(
    "a, b, expected",
    [
        ("abcdef", "abcdef", 0),
        # Deletions / additions
        ("abcdef", "abcde", 1),
        ("abcdef", "abcd", 2),
        ("abcdef", "abc", 3),
        ("abcdf", "abcdef", 1),
        ("cdef", "abcdef", 2),
        ("bde", "abcdef", 3),
        # Substitutions
        ("abcd", "abc_", 1),
        ("abcd", "ab__", 2),
        ("abcd", "a___", 3),
        ("abcd", "____", 4),
        # Deletion + Substitutions
        ("abcd", "abc_z", 2),
        ("abcd", "ab__zz", 4),
        ("abcd", "a___zzz", 6),
        ("abcd", "____zzzz", 8),
    ],
)
def test_levenshtein(a, b, expected):
    with pytest.deprecated_call():
        assert levenshtein_distance(a, b) == expected


@pytest.mark.parametrize(
    "a, b, max_dist, expected",
    [
        ("abcd", "", 0, 1),
        ("abcd", "", 3, 4),
        ("abcd", "", 10, 4),
        ("abcd", "", 1, 2),
        ("abcd", "a", 2, 3),
        ("abcd", "ad", 2, 2),
        ("abcd", "abd", 2, 1),
        ("abcd", "abcd", 2, 0),
    ],
)
def test_levenshtein_with_max_dist(a, b, max_dist, expected):
    with pytest.deprecated_call():
        assert levenshtein_distance(a, b, max_dist=max_dist) == expected
