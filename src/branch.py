__author__ = 'nikita_kartashov'


def does_intersect(branch, tree):
    return not frozenset(tree[0]).issuperset(branch[0])


def compute_tree_score_with_branches(branches, tree):
    return sum(score for branch, score in branches if not does_intersect(branch, tree))
