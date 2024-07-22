import functools

def regexes_are_valid(regexes):
    all_groups = functools.reduce(lambda x, y: x | y.groupindex.keys(), regexes, set())

    sample = [group for group in all_groups if group.startswith("sample")]
    lid = [group for group in all_groups if group.startswith("lid")]

    return len(sample) > 0 and len(lid) > 0


