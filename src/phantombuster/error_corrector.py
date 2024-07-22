from math import floor, ceil
from collections import defaultdict, namedtuple, Counter
import sys
import logging


# Immutable representation of an UMI
Umi = namedtuple("Umi", ["start", "end", "tag"])
# A UMI, consisting of start- and end-position and tag
Result = namedtuple("Result", ["ecumi_counts", "ecumi_list"])

def hamming_distance(s1: str, s2: str) -> int:
    """Computed the (unscaled) hamming distance between two lists"""
    if len(s1) != len(s2):
        raise RuntimeError(
            "hamming_distance required sequences with the same length, but give %d and %d characters"
            % (len(s1), len(s2))
        )
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


    def increment(self):
        self.value.value += 1

    def get_value(self):
        return self.value.value


class UmiEC:
    """
    A mutable node in a UMI similarity graph/forest.

    Mutable representation of a node in the UMI forest. The field 'umi' contains
    an Umi, and 'parents' a set of parent nodes, which the Umi is similar too,
    *and* which have a higher total read count. If the total read counts of two
    similar Umis are equal, the less-than (<) operator of the Umi type is used to
    break the tie, and the "greater" Umi is made the parent. This ensures
    non-circularity of the graph in all cases.
    """

    def __init__(self, umi, parents):
        self.umi = umi
        self.parents = parents

    # Make hashable, use object identity as key (as is appropriate for a mutable type)
    def __hash__(self):
        return id(self)

    # Identity means object identity (see __hash__ above)
    def __eq__(self, other):
        return id() == other.id()


class UmiCount:
    """ "
    Per-UMI read count for data without stranded UMIs

    Tracks not only the total read count, but also the number of UMIs merged into
    one during error-correction (`rawumis`).
    """

    @classmethod
    def zero(cls):
        return cls(0, 0)

    @classmethod
    def rawumi(cls):
        return cls(0, 1)

    def __init__(self, count, rawumis):
        self.rawumis = rawumis
        self.count = count

    def total(self):
        return self.count

    def __add__(self, other):
        return UmiCount(self.count + other.count, self.rawumis + other.rawumis)

    def __eq__(self, other):
        return self.count == other.count

    def __ne__(self, other):
        return not (self == other)

    def __str__(self):
        return str(self.count)

    def __repr__(self):
        return str(self.count)


class UmiCountStranded:
    """ "
    Per-UMI read count for stranded UMIs, i.e. if both strands are detectable

    Tracks the read counts of both strands separately, and also the number of UMIs
    merged into one during error-correction (`rawumis`).
    """

    @classmethod
    def zero(cls):
        return cls(0, 0, 0)

    @classmethod
    def rawumi(cls):
        return cls(0, 0, 1)

    def __init__(self, plus, minus, rawumis):
        self.rawumis = rawumis
        self.plus = plus
        self.minus = minus

    def total(self):
        return self.plus + self.minus

    def __add__(self, other):
        return UmiCountStranded(
            self.plus + other.plus,
            self.minus + other.minus,
            self.rawumis + other.rawumis,
        )

    def __eq__(self, other):
        return (self.plus == other.plus) and (self.minus == other.minus)

    def __ne__(self, other):
        return not (self == other)

    def __str__(self):
        return "%d+,%d-" % (self.plus, self.minus)

class ErrorCorrector:
    """
    UMI error correction.

    Two Umis whose start and end-positions each differ by at most max_mapping_distance
    bases, and whose tags have hamming distance at most max_hamming_distance are merged
    into one.
    """

    def __init__(
        self,
        tag_length,
        max_mapping_distance_start,
        max_mapping_distance_end,
        max_hamming_distance,
        count_type=UmiCount,
        report=False,
    ):
        self.count_type = count_type
        self.tag_length = tag_length
        self.max_mapping_distance_start = max_mapping_distance_start
        self.max_mapping_distance_end = max_mapping_distance_end
        self.max_hamming_distance = max_hamming_distance
        self.tsw = int(floor(tag_length / (self.max_hamming_distance + 1)))
        self.report = report
        self.report_step_umis = 500
        self.report_step_basepairs = 50

    def process(self, umi_counts):
        if self.report:
            logging.debug(
                "%d UMIs with %s reads before error-correction\n",
                 len(umi_counts), sum(umi_counts.values(), self.count_type.zero())
            )

        # Build forest
        ecinfo = self.build_forest(umi_counts)
        # Find roots and sum counts
        result = self.sum_trees(umi_counts, ecinfo)
        ecumi_counts = result.ecumi_counts

        if self.report:
            sys.stderr.write(
                "%d UMIs with %s reads after error-correction\n"
                % (
                    len(ecumi_counts),
                    sum(ecumi_counts.values(), self.count_type.zero()),
                )
            )

        return result

    def link_if_similar(self, umi, umi2, umi_counts, ecinfo):
        # Debugging support.
        log = False
        if log:
            logging.debug("UMI1: %s x %s", umi, umi_counts[umi])
            logging.debug("UMI2: %s x %s", umi2, umi_counts[umi2])
        # Check mapping distance (start)
        if abs(umi2.start - umi.start) > self.max_mapping_distance_start:
            # Mapping positions lie too far apart to unify UMIs.
            if log:
                logging.debug("DEBUG: NOT MERGED, start distance > %d", self.max_mapping_distance_start)
            return
        elif umi2.start < umi.start:
            # Impossible, list is sorted!
            raise RuntimeError("internal error, UMI list corrupted")
        # Check mapping distance (end)
        if abs(umi2.end - umi.end) > self.max_mapping_distance_end:
            # Mapping positions lie too far apart to unify UMIs
            if log:
                logging.debug("NOT MERGED, end distance > %d", self.max_mapping_distance_end)
            return
        # Check hamming distance
        hd = False
        if self.max_hamming_distance > 0:
            hd = hamming_distance(umi.tag, umi2.tag) > self.max_hamming_distance
        else:
            hd = umi.tag != umi2.tag
        if hd:
            # Hamming distance too large to allow unification of UMIs
            if log:
                logging.debug("DEBUG: NOT MERGED, hamming distance > %d", self.max_hamming_distance)
            return
        # Determine which to merge into which
        umi_cnt = umi_counts[umi].total()
        umi2_cnt = umi_counts[umi2].total()
        if (umi_cnt, umi) < (umi2_cnt, umi2):
            p, c = umi2, umi
        elif (umi2_cnt, umi2) < (umi_cnt, umi):
            p, c = umi, umi2
        else:
            raise RuntimeError("Tertium non datur fails for %s and %s" % (umi, umi2))
        if log:
            logging.debug("WILL MERGE %s x %s --> %s x %s", c, umi_counts[c], p, umi_counts[p])
        # Merge UMIs as described above
        ecinfo[c].parents.add(ecinfo[p])

    # Splits a tag into enough disjoint parts such that whenever the hamming
    # distance of two tags is at most max_hamming_distance, then at least one
    # part is identical in both tags.
    def splittag(self, tag):
        if len(tag) != self.tag_length:
            raise ValueError("tag %s has wrong length, %d instead of %d" % (tag, len(tag), self.tag_length))
        return [
            tag[i * self.tsw : ((i + 1) * self.tsw if i < self.max_hamming_distance else self.tag_length)]
            for i in range(self.max_hamming_distance + 1)
        ]

    def build_forest(self, umi_counts):
        if self.report:
            logging.debug("Finding similar UMIs...\n")

        # Data-structure which represents the forest. ecinfo associates each Umi
        # with an UmiEC instance, which contains the Umi itself, and a set of
        # parents. Note that since parents is a set, the "forest" does not actually
        # contain trees, but DAGs.
        ecinfo = dict()
        for u in umi_counts.keys():
            ecinfo[u] = UmiEC(u, set())

        # Build Umi index. The first level splits Umis by their starting position.
        # Below that, there are max_hamming_distance + 1 associative maps (dicts),
        # M_k, one for each of the parts returned by splittag(). The map M_k maps a
        # string s to the set of Umis whose k-th tag part equals s. Therefore, given
        # a tag x, all tags y whose hamming distance is at most max_hamming_distance,
        # are contained in at least one of the sets M_k(x_k), where x_k is the k-th
        # tag part of x.
        umi_index = defaultdict(lambda: [defaultdict(set) for i in range(0, self.max_hamming_distance + 1)])
        for umi in umi_counts.keys():
            t = umi_index[umi.start]
            for i, p in enumerate(self.splittag(umi.tag)):
                t[i][p].add(umi)

        # Populate ecinfo
        last_pos = 0
        last_report_pos = 0
        # Scan Umis
        for i, umi in enumerate(sorted(umi_counts.keys(), key=lambda u: u.start)):
            if umi.start < last_pos:
                raise RuntimeError("UMIs not sorted by start position (%d after %d)!" % (umi.start, last_pos))
            if self.report and (i % self.report_step_umis == 0):
                logging.debug("  processed %d UMIs\n", i)
            if self.report and ((umi.start - last_report_pos) >= self.report_step_basepairs):
                last_report_pos = umi.start
                logging.debug("  processed %d basepairs\n", umi.start)
            last_pos = umi.start
            # Scan part of first index level corresponding to a similar start position.
            for s2 in range(umi.start, umi.start + self.max_mapping_distance_start + 1):
                # Scan tag parts of 1st Umi
                for k2, p2 in enumerate(self.splittag(umi.tag)):
                    # Scan Umis having the same k2-th tag part p2 as the 1st Umi
                    for umi2 in umi_index[s2][k2][p2]:
                        # Avoid calling link_if_similar twice for every pair, using the total
                        # order imbued on the Umis to decide which call to suppress. Note that
                        # whatever condition we use here *must* ensure that if umi has a lower
                        # start position than umi2, link_if_similar *will* be called, since we
                        # only let s2 range within [ umi.start, umi.start + max_dist ] above.
                        if umi < umi2:
                            self.link_if_similar(umi, umi2, umi_counts, ecinfo)

        # Return forest
        return ecinfo

    def sum_trees(self, umi_counts, ecinfo):
        if self.report:
            logging.debug("Merging similar UMIs...\n")

        # All counts within an UMI tree contribute to the roots count.
        # If the tree is actually a net, i.e. if there is more than one "root"
        # reachable from a certain leaf, that leaf's assignment is ambiguous,
        # and its count is not attributed to any root.
        processed = 0
        ambiguous = 0
        ambiguous_count = self.count_type.zero()
        ecumi_counts = defaultdict(lambda: self.count_type.zero())
        ecumi_list = {}
        for umi, umi_ec in ecinfo.items():
            if self.report and (processed % self.report_step_umis == 0):
                logging.debug("  processed %d UMIs\n" % processed)
            processed += 1
            # Start with current node as "active set"
            nodes = set([umi_ec])
            roots = []
            while nodes:
                # Pick node in active set,
                n = nodes.pop()
                if n.parents:
                    # and replace it by its parents
                    nodes.update(set(n.parents))
                else:
                    roots.append(n.umi)
            if len(roots) == 1:
                # Found a unique root, so add node's counts to root's counts
                if umi_counts[umi].rawumis != 1:
                    raise RuntimeError("before error correction, umis should have rawumis=1")
                ecumi_counts[roots[0]] += umi_counts[umi]
                ecumi_list.update({umi: roots[0]})
            elif len(roots) > 1:
                # Found more than one root, so throw node away
                if umi_counts[umi].rawumis != 1:
                    raise RuntimeError("before error correction, umis should have rawumis=1")
                ambiguous_count += umi_counts[umi]
            else:
                ecumi_list.update({umi: umi})
        # Check that we haven't missed any UMI
        if sum(umi_counts.values(), self.count_type.zero()) != (
            sum(ecumi_counts.values(), self.count_type.zero()) + ambiguous_count
        ):
            raise RuntimeError("error correction buggy!")
        if self.report:
            logging.debug("%d UMIs with %s reads were ambiguous\n", ambiguous_count.rawumis, ambiguous_count)

        result = Result(ecumi_counts, ecumi_list)
        return result
