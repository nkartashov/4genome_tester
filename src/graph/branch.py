__author__ = 'nikita_kartashov'


def does_intersect(branch, tree):
    return not frozenset(tree[0]).issuperset(branch[0]) and not frozenset(tree[1]).issuperset(branch[0])


def compute_tree_score_with_branches(branches, tree):
    return sum(score for branch, score in branches if not does_intersect(branch, tree))


if __name__ == '__main__':
    def intersect_test():
        branch1 = (frozenset(), frozenset(['A', 'B', 'C', 'D']))
        branch2 = (frozenset(['A', 'B']), frozenset(['C', 'D']))
        branch3 = (frozenset(['C', 'B']), frozenset(['C', 'B']))
        tree1 = (('A', 'B'), ('C', 'D'))
        tree2 = (('C', 'B'), ('A', 'D'))
        assert (not does_intersect(branch1, tree1))
        assert (not does_intersect(branch1, tree2))
        assert (not does_intersect(branch2, tree1))
        assert (does_intersect(branch2, tree2))
        assert (does_intersect(branch3, tree1))
        assert (not does_intersect(branch3, tree2))

    intersect_test()